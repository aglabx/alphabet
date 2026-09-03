//! `annotate-reads` — annotate sequences with the alphabet, in one call.
//!
//! Takes FASTA or FASTQ (plain or gzipped) and writes, per record, the ordered
//! letter string and the coordinate of every monomer, plus the per-monomer
//! table if you want it. It is `cut-v2 --auto-orient` followed by `assign`,
//! without the intermediate files.
//!
//! **Platform.** `--platform hifi` (also the right mode for a FASTA of
//! assembled arrays) cuts each record into monomers and places each monomer in
//! the roster. That is the correct primitive when a record is long enough to
//! contain whole monomers.
//!
//! Illumina is deliberately **not** a mode here. A 101/151 bp read covers part
//! of one ~171 bp monomer, so there is nothing to cut: the primitive is direct
//! classification of the read against the letter consensuses, which is a
//! different computation with a different output shape (a read gets one letter,
//! not a string). Adding it as a flag that silently ran the wrong algorithm
//! would be worse than not having it.
//!
//! Orientation is handled: reads arrive in random orientation and the anchor
//! table is canonical-strand, so a record whose forward pass finds no cut is
//! retried reverse-complemented. Coordinates are always reported in the frame
//! of the input record; `orient` says which strand the monomers were read on.

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use crate::cmd::assign_letters::Assigner;
use crate::cmd::cut_v2::{cut_array, Params};
use crate::monomer::revcomp;

#[derive(Parser)]
#[command(
    name = "annotate-reads",
    about = "Annotate FASTA/FASTQ records with the alphabet: letter string + coordinates per record"
)]
struct Args {
    /// Input FASTA or FASTQ, plain or .gz. Format detected from the first byte.
    #[arg(long, value_name = "PATH")]
    input: String,

    /// Roster TSV: letter <TAB> consensus.
    #[arg(long, value_name = "PATH")]
    roster: String,

    /// Per-column profile counts TSV: letter <TAB> col <TAB> A <TAB> C <TAB> G <TAB> T.
    #[arg(long, value_name = "PATH")]
    profiles: String,

    /// Output: one line per record — letters and coordinates.
    #[arg(long = "out-reads", value_name = "PATH")]
    out_reads: String,

    /// Output: one line per monomer (optional).
    #[arg(long = "out-monomers", value_name = "PATH")]
    out_monomers: Option<String>,

    /// Sequencing platform. `hifi` also covers a FASTA of assembled arrays.
    #[arg(long, default_value = "hifi", value_parser = ["hifi", "fasta"])]
    platform: String,

    /// Expected monomer period in bp.
    #[arg(long, default_value_t = 171)]
    period: usize,

    /// ± tolerance per cut step.
    #[arg(long = "tol-period", default_value_t = 30)]
    tol_period: usize,

    /// Minimum anchor score to accept a cut.
    #[arg(long = "min-score", default_value_t = 2.5)]
    min_score: f64,

    /// Maximum consecutive interpolated monomers before giving up on a record.
    #[arg(long = "max-skip", default_value_t = 10)]
    max_skip: usize,

    /// k for the vote index.
    #[arg(long, default_value_t = 13)]
    k: usize,

    /// Candidates carried from the vote prefilter into the profile score.
    #[arg(long, default_value_t = 32)]
    topk: usize,

    /// Offset band for the profile score.
    #[arg(long, default_value_t = 8)]
    band: i64,

    /// log10LR at or below which a monomer is reported LOST.
    #[arg(long = "lost-thr", default_value_t = -0.430)]
    lost_thr: f64,

    /// Homopolymer-compress before matching. Measured neutral-to-negative on
    /// HiFi (it costs 31% of the distinct 13-mer keys); meant for ONT.
    #[arg(long, default_value_t = false)]
    hpc: bool,

    /// Records per parallel batch. Bounds memory on a large input.
    #[arg(long = "batch", default_value_t = 20000)]
    batch: usize,

    /// Rayon threads (0 = library default).
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

fn open_maybe_gz(path: &str) -> Box<dyn BufRead> {
    let mut f = File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic).unwrap_or_else(|e| panic!("cannot read {}: {}", path, e));
    let f = File::open(path).unwrap_or_else(|e| panic!("cannot reopen {}: {}", path, e));
    if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::with_capacity(1 << 24, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(1 << 24, f))
    }
}

/// Streaming FASTA/FASTQ reader. Yields `(id, sequence)`; the id is the header
/// up to the first whitespace.
struct SeqReader {
    rd: Box<dyn BufRead>,
    buf: String,
    pending: Option<String>,
    fastq: Option<bool>,
    path: String,
}

impl SeqReader {
    fn new(path: &str) -> SeqReader {
        SeqReader {
            rd: open_maybe_gz(path),
            buf: String::new(),
            pending: None,
            fastq: None,
            path: path.to_string(),
        }
    }

    fn next_record(&mut self) -> Option<(String, Vec<u8>)> {
        loop {
            let header = match self.pending.take() {
                Some(h) => h,
                None => {
                    self.buf.clear();
                    if self.rd.read_line(&mut self.buf).expect("read error") == 0 {
                        return None;
                    }
                    self.buf.trim_end().to_string()
                }
            };
            if header.is_empty() {
                continue;
            }
            if self.fastq.is_none() {
                self.fastq = Some(header.starts_with('@'));
                assert!(
                    header.starts_with('@') || header.starts_with('>'),
                    "{} does not start with '>' or '@': not FASTA or FASTQ",
                    self.path
                );
            }
            let id = header[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            if self.fastq == Some(true) {
                let mut seq = String::new();
                let mut plus = String::new();
                let mut qual = String::new();
                assert!(self.rd.read_line(&mut seq).expect("read error") > 0,
                    "truncated FASTQ record {} in {}", id, self.path);
                assert!(self.rd.read_line(&mut plus).expect("read error") > 0,
                    "truncated FASTQ record {} in {}", id, self.path);
                assert!(self.rd.read_line(&mut qual).expect("read error") > 0,
                    "truncated FASTQ record {} in {}", id, self.path);
                let s = seq.trim_end();
                assert_eq!(s.len(), qual.trim_end().len(),
                    "seq/qual length mismatch on {} in {}", id, self.path);
                return Some((id, s.as_bytes().to_ascii_uppercase()));
            }
            // FASTA: accumulate until the next header
            let mut seq: Vec<u8> = Vec::new();
            loop {
                self.buf.clear();
                if self.rd.read_line(&mut self.buf).expect("read error") == 0 {
                    break;
                }
                if self.buf.starts_with('>') {
                    self.pending = Some(self.buf.trim_end().to_string());
                    break;
                }
                seq.extend_from_slice(self.buf.trim_end().as_bytes());
            }
            seq.make_ascii_uppercase();
            return Some((id, seq));
        }
    }
}

struct MonoOut {
    idx: usize,
    start: usize,
    end: usize,
    orient: char,
    letter_km: String,
    votes_km: u32,
    letter_lr: String,
    log10lr: f64,
    lost: bool,
}

/// Cut a record, retrying the reverse complement, and place every monomer.
/// Coordinates come back in the frame of `seq`.
fn annotate_one(seq: &[u8], params: &Params, asg: &Assigner) -> (Vec<MonoOut>, &'static str) {
    let l = seq.len();
    if l < params.period {
        return (Vec::new(), "too_short");
    }
    let fwd = cut_array(seq, params);
    let (cuts, work, orient) = if fwd.cuts.len() >= 2 {
        (fwd.cuts, seq.to_vec(), '+')
    } else {
        let rc = revcomp(seq);
        let rev = cut_array(&rc, params);
        if rev.cuts.len() < 2 {
            return (Vec::new(), "no_cut_either_strand");
        }
        (rev.cuts, rc, '-')
    };

    let mut out = Vec::with_capacity(cuts.len() - 1);
    for i in 0..(cuts.len() - 1) {
        let (s, e) = (cuts[i], cuts[i + 1]);
        let a = asg.assign(&work[s..e]);
        let (rs, re) = if orient == '+' { (s, e) } else { (l - e, l - s) };
        out.push(MonoOut {
            idx: i,
            start: rs,
            end: re,
            orient,
            letter_km: a.km_letter,
            votes_km: a.km_votes,
            letter_lr: a.lr_letter,
            log10lr: a.log10lr,
            lost: a.lost,
        });
    }
    // left-to-right in the input frame, so coordinates read in file order
    out.sort_by_key(|m| m.start);
    (out, "ok")
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(argv);
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }
    assert!(
        args.platform == "hifi" || args.platform == "fasta",
        "unsupported platform {:?}; Illumina needs direct read classification, not cutting",
        args.platform
    );

    let params = Params {
        period: args.period,
        tol_period: args.tol_period,
        min_score: args.min_score,
        max_skip: args.max_skip,
    };
    let asg = Assigner::load(
        &args.roster, &args.profiles, args.k, args.topk, args.band, args.lost_thr, args.hpc,
    );

    let mut wr = BufWriter::with_capacity(
        1 << 22,
        File::create(&args.out_reads)
            .unwrap_or_else(|e| panic!("cannot create {}: {}", args.out_reads, e)),
    );
    writeln!(wr, "read_id\tread_len\tstatus\torient\tn_monomers\tstarts\tletters_km\tvotes_km\tletters_lr\tlog10LR\tlr_class").unwrap();
    let mut wm = args.out_monomers.as_ref().map(|p| {
        let mut w = BufWriter::with_capacity(
            1 << 22,
            File::create(p).unwrap_or_else(|e| panic!("cannot create {}: {}", p, e)),
        );
        writeln!(w, "read_id\tmono_idx\tstart\tend\tlength\torient\tletter_km\tvotes_km\tletter_lr\tlog10LR\tlr_class").unwrap();
        w
    });

    let mut rdr = SeqReader::new(&args.input);
    let (mut n, mut ok, mut nmono) = (0u64, 0u64, 0u64);
    let mut refused: std::collections::HashMap<&'static str, u64> = Default::default();
    loop {
        let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(args.batch);
        while batch.len() < args.batch {
            match rdr.next_record() {
                Some(r) => batch.push(r),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        let done: Vec<(String, usize, Vec<MonoOut>, &'static str)> = batch
            .par_iter()
            .map(|(id, seq)| {
                let (m, st) = annotate_one(seq, &params, &asg);
                (id.clone(), seq.len(), m, st)
            })
            .collect();
        for (id, len, ms, st) in done {
            n += 1;
            if st != "ok" {
                *refused.entry(st).or_insert(0) += 1;
                // a refused record keeps its line: the count in the file must
                // match the count in the input
                writeln!(wr, "{}\t{}\t{}\t.\t0\t\t\t\t\t\t", id, len, st).unwrap();
                continue;
            }
            ok += 1;
            nmono += ms.len() as u64;
            let orient = ms[0].orient;
            let j = |f: &dyn Fn(&MonoOut) -> String| ms.iter().map(f).collect::<Vec<_>>().join(",");
            writeln!(wr, "{}\t{}\tok\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                id, len, orient, ms.len(),
                j(&|m: &MonoOut| m.start.to_string()),
                j(&|m: &MonoOut| m.letter_km.clone()),
                j(&|m: &MonoOut| m.votes_km.to_string()),
                j(&|m: &MonoOut| m.letter_lr.clone()),
                j(&|m: &MonoOut| format!("{:.3}", m.log10lr)),
                j(&|m: &MonoOut| if m.lost { "LOST".into() } else { "ok".to_string() })
            ).unwrap();
            if let Some(w) = wm.as_mut() {
                for m in &ms {
                    writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
                        id, m.idx, m.start, m.end, m.end - m.start, m.orient,
                        m.letter_km, m.votes_km, m.letter_lr, m.log10lr,
                        if m.lost { "LOST" } else { "ok" }).unwrap();
                }
            }
        }
    }
    wr.flush().unwrap();
    if let Some(w) = wm.as_mut() {
        w.flush().unwrap();
    }
    eprintln!("records: {}  annotated: {}  monomers: {}", n, ok, nmono);
    for (k, v) in &refused {
        eprintln!("  refused {}: {}", k, v);
    }
}
