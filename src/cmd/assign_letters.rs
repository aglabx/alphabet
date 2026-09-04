//! `assign` — put monomers into a **fixed** letter roster.
//!
//! `cut` discovers letters while cutting; this command does the other thing:
//! takes an alphabet that already exists and asks which of its letters each
//! monomer is. That is what you need to apply a published alphabet to a new
//! individual, a new species, or a set of reads.
//!
//! Two rules are computed in one pass, because the project contains two and
//! they are not interchangeable:
//!
//! * **profile-LR** — per-column profiles theta_f[col][base] (Jeffreys 0.5),
//!   best score over an offset band, log10 likelihood ratio against an order-0
//!   background. Substitution-only: the band absorbs a frame shift, there is no
//!   internal indel term. Carries a calibrated reject class.
//! * **k-mer votes** — k-mer hits against the letter consensuses, argmax of the
//!   vote count. Frame-free, so an indel costs only the k-mers it touches. This
//!   is the rule that built the published `lattice_nodes.tsv`.
//!
//! Both share one k-mer index: it is the vote counter for the second and the
//! candidate prefilter for the first, so no monomer is ever scored against all
//! 868 letters.
//!
//! Thresholds are **reported, not applied**: votes and log10LR are columns, so
//! a minimum-vote cut or a different reject threshold stays a downstream filter
//! over the same file.

use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

const ALPHA: f64 = 0.5; // Jeffreys pseudocount
const LN10: f64 = std::f64::consts::LN_10;

#[derive(Parser)]
#[command(
    name = "assign",
    about = "Assign monomers to a fixed letter roster (profile-LR + k-mer votes)"
)]
struct Args {
    /// Roster TSV: letter <TAB> consensus sequence, one line per letter.
    #[arg(long, value_name = "PATH")]
    roster: String,

    /// Per-column profile counts TSV: letter <TAB> col <TAB> A <TAB> C <TAB> G <TAB> T.
    /// Required for profile-LR; without it only the k-mer columns are filled.
    #[arg(long, value_name = "PATH")]
    profiles: String,

    /// Monomer TSV as written by `cut-v2`: the sequence must be column 9.
    #[arg(long, value_name = "PATH")]
    monomers: String,

    /// Output TSV.
    #[arg(long, value_name = "PATH")]
    out: String,

    /// Offset band for the profile score, in bp. Absorbs a frame shift.
    #[arg(long, default_value_t = 8)]
    band: i64,

    /// k for the vote index.
    #[arg(long, default_value_t = 13)]
    k: usize,

    /// Candidates carried from the vote prefilter into the profile score.
    #[arg(long, default_value_t = 32)]
    topk: usize,

    /// log10LR at or below which a monomer is reported LOST rather than
    /// assigned. Default is the calibrated value (order-0 null, p99 of shuffled
    /// monomers).
    #[arg(long = "lost-thr", default_value_t = -0.430)]
    lost_thr: f64,

    /// Homopolymer-compress the monomer and the consensus before matching.
    ///
    /// The obvious answer to indel error, and on PacBio HiFi it does not pay:
    /// measured neutral-to-negative in every MAPQ stratum for both rules, and
    /// it costs 31% of the distinct 13-mer keys (35,543 -> 24,552) because run
    /// lengths carry letter identity. After CCS the indel penalty is not the
    /// binding constraint. Kept for ONT, where the residual indel rate is an
    /// order of magnitude higher and the trade may reverse.
    #[arg(long, default_value_t = false)]
    hpc: bool,

    /// Rayon threads (0 = library default).
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[inline]
pub fn bidx(b: u8) -> usize {
    match b {
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'A' | b'a' => 0,
        _ => 4,
    }
}

/// Collapse each run of identical bases to one base. Returns the compressed
/// sequence and, per compressed position, the index of the run's first base.
pub fn hpc_compress(seq: &[u8]) -> (Vec<u8>, Vec<usize>) {
    let mut out = Vec::with_capacity(seq.len());
    let mut src = Vec::with_capacity(seq.len());
    let mut prev = 0u8;
    for (i, &b) in seq.iter().enumerate() {
        let u = b.to_ascii_uppercase();
        if out.is_empty() || u != prev {
            out.push(u);
            src.push(i);
            prev = u;
        }
    }
    (out, src)
}

pub struct Letter {
    pub name: String,
    pub cons: Vec<u8>,
    pub logtheta: Vec<[f64; 4]>,
}

/// One monomer's placement in the roster under both rules.
pub struct Assignment {
    pub lr_letter: String,
    pub log10lr: f64,
    pub margin_lr: f64,
    pub lost: bool,
    pub km_letter: String,
    pub km_votes: u32,
    pub km_margin: i64,
}

/// The roster, its profiles and the shared k-mer index, ready to place monomers.
/// Built once, used from any number of threads.
pub struct Assigner {
    letters: Vec<Letter>,
    index: HashMap<u32, Vec<u16>>,
    logbg: [f64; 4],
    k: usize,
    topk: usize,
    band: i64,
    lost_thr: f64,
    hpc: bool,
}

impl Assigner {
    pub fn load(
        roster_path: &str,
        profiles_path: &str,
        k: usize,
        topk: usize,
        band: i64,
        lost_thr: f64,
        hpc: bool,
    ) -> Assigner {
        let prof = load_profiles(profiles_path);
        let roster = load_roster(roster_path);
        eprintln!("roster letters: {}, profiled letters: {}", roster.len(), prof.len());

        let mut comp = [0f64; 4];
        for cols in prof.values() {
            for c in cols {
                for j in 0..4 {
                    comp[j] += c[j] as f64;
                }
            }
        }
        let tot: f64 = comp.iter().sum();
        assert!(tot > 0.0, "{} carries no counts", profiles_path);
        let logbg: [f64; 4] = [
            (comp[0] / tot).ln(),
            (comp[1] / tot).ln(),
            (comp[2] / tot).ln(),
            (comp[3] / tot).ln(),
        ];
        eprintln!(
            "order-0 background A/C/G/T = {:.4} {:.4} {:.4} {:.4}",
            comp[0] / tot, comp[1] / tot, comp[2] / tot, comp[3] / tot
        );

        let mut letters: Vec<Letter> = Vec::with_capacity(roster.len());
        let mut missing = Vec::new();
        for (name, cons) in &roster {
            match prof.get(name) {
                // a roster letter with no profile means the two artifacts are not
                // the same alphabet; skipping it would silently shrink the alphabet
                None => missing.push(name.clone()),
                Some(counts) => {
                    let n = counts.len().min(cons.len());
                    assert!(n > 0, "letter {} has an empty profile", name);
                    let logtheta: Vec<[f64; 4]> = counts[..n]
                        .iter()
                        .map(|c| {
                            let t: f64 = c.iter().map(|&x| x as f64).sum::<f64>() + 4.0 * ALPHA;
                            [
                                ((c[0] as f64 + ALPHA) / t).ln(),
                                ((c[1] as f64 + ALPHA) / t).ln(),
                                ((c[2] as f64 + ALPHA) / t).ln(),
                                ((c[3] as f64 + ALPHA) / t).ln(),
                            ]
                        })
                        .collect();
                    if hpc {
                        let (ccons, src) = hpc_compress(&cons[..n]);
                        let clog: Vec<[f64; 4]> = src.iter().map(|&i| logtheta[i]).collect();
                        letters.push(Letter { name: name.clone(), cons: ccons, logtheta: clog });
                    } else {
                        letters.push(Letter {
                            name: name.clone(),
                            cons: cons[..n].to_vec(),
                            logtheta,
                        });
                    }
                }
            }
        }
        assert!(
            missing.is_empty(),
            "{} roster letters have no profile ({}...) -- roster and profiles are not the same alphabet",
            missing.len(),
            missing.iter().take(5).cloned().collect::<Vec<_>>().join(",")
        );
        eprintln!("letters usable: {}", letters.len());

        let mut index: HashMap<u32, Vec<u16>> = HashMap::new();
        for (li, l) in letters.iter().enumerate() {
            let mut seen: Vec<u32> = Vec::new();
            kmers(&l.cons, k, |c| seen.push(c));
            seen.sort_unstable();
            seen.dedup();
            for c in seen {
                index.entry(c).or_default().push(li as u16);
            }
        }
        eprintln!("{}-mer index: {} distinct keys (hpc={})", k, index.len(), hpc);

        Assigner { letters, index, logbg, k, topk, band, lost_thr, hpc }
    }

    pub fn assign(&self, raw: &[u8]) -> Assignment {
        let compressed;
        let seq: &[u8] = if self.hpc {
            compressed = hpc_compress(raw).0;
            &compressed
        } else {
            raw
        };

        let mut votes: HashMap<u16, u32> = HashMap::new();
        kmers(seq, self.k, |c| {
            if let Some(ls) = self.index.get(&c) {
                for &li in ls {
                    *votes.entry(li).or_insert(0) += 1;
                }
            }
        });
        let mut vv: Vec<(u16, u32)> = votes.into_iter().collect();
        // deterministic: votes desc, then letter index asc
        vv.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        let (km_letter, km_votes, km_margin) = match vv.first() {
            None => ("NONE".to_string(), 0u32, 0i64),
            Some(&(li, v)) => (
                self.letters[li as usize].name.clone(),
                v,
                v as i64 - vv.get(1).map(|x| x.1 as i64).unwrap_or(0),
            ),
        };

        let bgsum: f64 = seq
            .iter()
            .map(|&b| {
                let c = bidx(b);
                if c == 4 { 0.0 } else { self.logbg[c] }
            })
            .sum();
        let mut best = (f64::NEG_INFINITY, usize::MAX);
        let mut second = f64::NEG_INFINITY;
        for &(li, _) in vv.iter().take(self.topk) {
            let s = profile_score(seq, &self.letters[li as usize], &self.logbg, self.band);
            if s > best.0 {
                second = best.0;
                best = (s, li as usize);
            } else if s > second {
                second = s;
            }
        }
        let (lr_letter, log10lr, margin_lr) = if best.1 == usize::MAX {
            ("NONE".to_string(), f64::NAN, f64::NAN)
        } else {
            let l10 = (best.0 - bgsum) / LN10;
            let m = if second.is_finite() { (best.0 - second) / LN10 } else { f64::INFINITY };
            (self.letters[best.1].name.clone(), l10, m)
        };
        let lost = log10lr.is_finite() && log10lr <= self.lost_thr;

        Assignment { lr_letter, log10lr, margin_lr, lost, km_letter, km_votes, km_margin }
    }
}

pub fn load_profiles(path: &str) -> HashMap<String, Vec<[u32; 4]>> {
    let f = BufReader::with_capacity(
        1 << 22,
        File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e)),
    );
    let mut m: HashMap<String, Vec<[u32; 4]>> = HashMap::new();
    for line in f.lines() {
        let l = line.expect("read error");
        if l.starts_with("letter\t") || l.is_empty() {
            continue;
        }
        let p: Vec<&str> = l.split('\t').collect();
        assert!(p.len() >= 6, "profiles row has {} fields, need 6: {}", p.len(), l);
        let col: usize = p[1].parse().unwrap_or_else(|_| panic!("bad col in {}", l));
        let v = m.entry(p[0].to_string()).or_default();
        if v.len() <= col {
            v.resize(col + 1, [0u32; 4]);
        }
        for j in 0..4 {
            v[col][j] = p[2 + j].parse().unwrap_or_else(|_| panic!("bad count in {}", l));
        }
    }
    m
}

pub fn load_roster(path: &str) -> Vec<(String, Vec<u8>)> {
    let f = BufReader::new(File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e)));
    let mut v = Vec::new();
    for line in f.lines() {
        let l = line.expect("read error");
        if l.is_empty() || l.starts_with("letter\t") {
            continue;
        }
        let mut it = l.split('\t');
        let name = it.next().expect("no letter column").to_string();
        let cons = it
            .next()
            .unwrap_or_else(|| panic!("roster row without a consensus: {}", l));
        v.push((name, cons.as_bytes().to_vec()));
    }
    v
}

#[inline]
pub fn kmers(seq: &[u8], k: usize, mut emit: impl FnMut(u32)) {
    let mut code: u32 = 0;
    let mut valid = 0usize;
    let mask: u32 = if k * 2 >= 32 { u32::MAX } else { (1u32 << (2 * k)) - 1 };
    for &b in seq {
        let c = bidx(b);
        if c == 4 {
            valid = 0;
            code = 0;
            continue;
        }
        code = ((code << 2) | c as u32) & mask;
        valid += 1;
        if valid >= k {
            emit(code);
        }
    }
}

/// Best profile log-likelihood over the offset band. Positions falling outside
/// the consensus score at background, which cancels in the ratio.
#[inline]
pub fn profile_score(x: &[u8], lf: &Letter, logbg: &[f64; 4], band: i64) -> f64 {
    let (lx, lc) = (x.len() as i64, lf.cons.len() as i64);
    if (lx - lc).abs() > band + 6 {
        return f64::NEG_INFINITY;
    }
    let mut best = f64::NEG_INFINITY;
    let mut o = -band;
    while o <= band {
        let mut s = 0.0;
        for i in 0..lx as usize {
            let c = bidx(x[i]);
            if c == 4 {
                continue;
            }
            let j = i as i64 + o;
            s += if j >= 0 && j < lc {
                lf.logtheta[j as usize][c]
            } else {
                logbg[c]
            };
        }
        if s > best {
            best = s;
        }
        o += 1;
    }
    best
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(argv);
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    let asg = Assigner::load(
        &args.roster,
        &args.profiles,
        args.k,
        args.topk,
        args.band,
        args.lost_thr,
        args.hpc,
    );

    let mf = BufReader::with_capacity(
        1 << 24,
        File::open(&args.monomers)
            .unwrap_or_else(|e| panic!("cannot open {}: {}", args.monomers, e)),
    );
    let rows: Vec<String> = mf
        .lines()
        .map(|l| l.expect("read error"))
        .filter(|l| !l.starts_with("array_id") && !l.is_empty())
        .collect();
    eprintln!("monomers: {}", rows.len());

    let out: Vec<String> = rows
        .par_iter()
        .map(|l| {
            let p: Vec<&str> = l.split('\t').collect();
            assert!(p.len() >= 9, "monomer row has {} fields, need at least 9", p.len());
            // cut-v2 --auto-orient appends an `orient` column; carry it through so
            // the strand does not have to be re-joined from the monomer table.
            // Other producers put their own columns there (tiling tables carry
            // `source_letter`), so accept it only if it LOOKS like a strand --
            // taking it by position alone silently writes a letter id into the
            // orient field.
            let orient = match p.get(9).copied() {
                Some(v) if v == "+" || v == "-" || v == "." => v,
                _ => ".",
            };
            let a = asg.assign(p[8].as_bytes());
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}",
                p[0], p[1], p[2], p[3], p[4], orient,
                a.lr_letter, a.log10lr, a.margin_lr,
                if a.lost { "LOST" } else { "ok" },
                a.km_letter, a.km_votes, a.km_margin
            )
        })
        .collect();

    let mut w = BufWriter::with_capacity(
        1 << 22,
        File::create(&args.out).unwrap_or_else(|e| panic!("cannot create {}: {}", args.out, e)),
    );
    writeln!(
        w,
        "array_id\tmono_idx\tstart\tend\tlength\torient\tletter_lr\tlog10LR\tmargin_lr\tlr_class\tletter_km\tvotes_km\tmargin_km"
    )
    .unwrap();
    for l in &out {
        writeln!(w, "{}", l).unwrap();
    }
    w.flush().unwrap();
    eprintln!("wrote {} assignments to {}", out.len(), args.out);
}
