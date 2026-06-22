//! `scaffold-msa` — anchor-scaffolded MSA orchestrator.
//!
//! See `docs/anchor_scaffold_msa_brief.md`. Reads a monomer FASTA + a seed
//! consensus FASTA, runs the detect → chain → ond_align → tally → rebuild
//! → EM pipeline, and emits the brief's six output artifacts:
//!
//!   - `standard.jsonl` / `outliers.jsonl` / `exceptions.jsonl`
//!   - `consensus.fa`
//!   - `columns.tsv`
//!   - `summary.tsv`

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;

use crate::anchor_detect::SCAFFOLD_PANEL;
use crate::anchor_scaffold::{
    aligned_row, run_em, update_tally, Cause, ColumnTally, EmParams, MonomerOutcome, Stream,
};
use crate::io::read_fasta;

#[derive(Parser)]
#[command(
    name = "scaffold-msa",
    about = "Anchor-scaffolded MSA over canonical-rotation monomers"
)]
struct Args {
    /// FASTA of canonical-rotation monomers (one record per monomer).
    #[arg(long = "monomers", required = true, value_name = "PATH")]
    monomers: PathBuf,

    /// Seed consensus FASTA (single record); becomes the starting consensus.
    #[arg(long = "seed", required = true, value_name = "PATH")]
    seed: PathBuf,

    /// Output directory; all artifacts land here.
    #[arg(short = 'o', long = "outdir", default_value = ".")]
    outdir: PathBuf,

    /// Minimum primary anchors present for STANDARD/OUTLIER routing.
    /// Below this, monomers go to EXCEPTIONS.
    #[arg(long = "min-anchors", default_value_t = 6)]
    min_anchors: usize,

    /// Global cap on per-slot Hamming distance.
    #[arg(long = "max-hd", default_value_t = 0)]
    max_hd: usize,

    /// Maximum EM rounds.
    #[arg(long = "max-rounds", default_value_t = 3)]
    max_rounds: usize,

    /// Convergence threshold (fraction of columns that must remain unchanged
    /// to stop early). Brief default: < 0.1%.
    #[arg(long = "converge", default_value_t = 0.001)]
    converge: f64,

    /// Relative-chaining anchor detection: select the longest gap-consistent
    /// co-linear anchor chain instead of fixed absolute-position windows.
    /// Floats with indel drift; invariant to frame offset. Required for real
    /// data (see BTN-scaffold-panel-fixed-windows-dont-transfer).
    #[arg(long = "chain", default_value_t = false)]
    chained: bool,

    /// Also write a human-viewable MSA FASTA: one row per aligned monomer over
    /// the canonical columns (consensus base / Sub / `-` Del / `.` uncovered),
    /// header `>{field_id}#{mono} chrom={chr}`. Exceptions are omitted.
    #[arg(long = "emit-msa", value_name = "PATH")]
    emit_msa: Option<PathBuf>,
}

/// Extract `chrN` (or `chrN_...`) from a monomer id like
/// `chm13|chr12:34..|orig=R|..#14`.
fn chrom_of(id: &str) -> &str {
    if let Some(i) = id.find("chr") {
        let rest = &id[i..];
        let end = rest
            .find([':', '|', '#'].as_slice())
            .unwrap_or(rest.len());
        &rest[..end]
    } else {
        "?"
    }
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(&argv);
    if let Err(e) = run(args) {
        eprintln!("scaffold-msa error: {e:?}");
        std::process::exit(1);
    }
}

fn run(args: Args) -> Result<()> {
    let monomers = read_fasta(args.monomers.to_str().unwrap())
        .with_context(|| format!("reading monomers FASTA {}", args.monomers.display()))?;
    eprintln!("Loaded {} monomers from {}", monomers.len(), args.monomers.display());

    let seed_recs = read_fasta(args.seed.to_str().unwrap())
        .with_context(|| format!("reading seed FASTA {}", args.seed.display()))?;
    let (seed_name, seed_seq) = seed_recs
        .into_iter()
        .next()
        .context("seed FASTA has no records")?;
    eprintln!("Seed: {} ({} bp)", seed_name, seed_seq.len());

    std::fs::create_dir_all(&args.outdir)
        .with_context(|| format!("create outdir {}", args.outdir.display()))?;

    let params = EmParams {
        max_rounds: args.max_rounds,
        convergence_threshold: args.converge,
        min_anchors: args.min_anchors,
        max_hd: args.max_hd,
        chained: args.chained,
    };

    let t_em = std::time::Instant::now();
    let result = run_em(&seed_seq, &monomers, SCAFFOLD_PANEL, &params);
    let dt_em = t_em.elapsed();

    let (n_std, n_out, n_exc, cause_counts) = emit_jsonl(&args.outdir, &result.final_outcomes)?;
    emit_consensus(&args.outdir, &seed_name, &result.final_consensus)?;
    emit_columns(&args.outdir, &result.final_outcomes, &result.final_consensus)?;
    emit_summary(&args.outdir, n_std, n_out, n_exc, &cause_counts)?;
    if let Some(msa_path) = &args.emit_msa {
        emit_msa_file(msa_path, &result.final_outcomes, &result.final_consensus)?;
    }

    let n_mono = monomers.len();
    let rate = if dt_em.as_secs_f64() > 0.0 {
        n_mono as f64 / dt_em.as_secs_f64()
    } else {
        0.0
    };
    eprintln!();
    eprintln!(
        "EM rounds: {}, columns changed last round: {}",
        result.rounds.len(),
        result.rounds.last().map(|r| r.n_columns_changed).unwrap_or(0)
    );
    eprintln!("streams: standard={n_std} outlier={n_out} exception={n_exc}");
    eprintln!("wall: {:.2}s  ({:.0} monomers/s)", dt_em.as_secs_f64(), rate);
    eprintln!("outputs in {}", args.outdir.display());

    Ok(())
}

fn cause_label(c: Cause) -> &'static str {
    match c {
        Cause::Ok => "OK",
        Cause::TooFewAnchors => "TOO_FEW_ANCHORS",
        Cause::SegmentOverBand => "SEGMENT_OVER_BAND",
        Cause::GapUnalignable => "GAP_UNALIGNABLE",
        Cause::NoncanonicalAnchorOrder => "NONCANONICAL_ANCHOR_ORDER",
    }
}

fn emit_jsonl(
    outdir: &std::path::Path,
    outcomes: &[MonomerOutcome],
) -> Result<(usize, usize, usize, BTreeMap<&'static str, usize>)> {
    let std_path = outdir.join("standard.jsonl");
    let out_path = outdir.join("outliers.jsonl");
    let exc_path = outdir.join("exceptions.jsonl");
    let mut w_std =
        BufWriter::new(File::create(&std_path).with_context(|| format!("create {}", std_path.display()))?);
    let mut w_out =
        BufWriter::new(File::create(&out_path).with_context(|| format!("create {}", out_path.display()))?);
    let mut w_exc =
        BufWriter::new(File::create(&exc_path).with_context(|| format!("create {}", exc_path.display()))?);

    let mut n_std = 0usize;
    let mut n_out = 0usize;
    let mut n_exc = 0usize;
    let mut cause_counts: BTreeMap<&'static str, usize> = BTreeMap::new();

    for o in outcomes {
        let line = serde_json::to_string(o).context("serialize outcome")?;
        match o.stream {
            Stream::Standard => {
                writeln!(w_std, "{line}")?;
                n_std += 1;
            }
            Stream::Outlier => {
                writeln!(w_out, "{line}")?;
                n_out += 1;
            }
            Stream::Exception => {
                writeln!(w_exc, "{line}")?;
                n_exc += 1;
            }
        }
        *cause_counts.entry(cause_label(o.cause)).or_insert(0) += 1;
    }
    w_std.flush()?;
    w_out.flush()?;
    w_exc.flush()?;
    Ok((n_std, n_out, n_exc, cause_counts))
}

fn emit_msa_file(
    path: &std::path::Path,
    outcomes: &[MonomerOutcome],
    consensus: &[u8],
) -> Result<()> {
    let mut w = BufWriter::new(
        File::create(path).with_context(|| format!("create {}", path.display()))?,
    );
    let mut n = 0usize;
    for o in outcomes {
        if o.stream == Stream::Exception {
            continue;
        }
        let st = match o.stream {
            Stream::Standard => "standard",
            Stream::Outlier => "outlier",
            Stream::Exception => "exception",
        };
        writeln!(w, ">{} chrom={} stream={} n_anchors={}", o.id, chrom_of(&o.id), st, o.n_present)?;
        w.write_all(&aligned_row(o, consensus))?;
        writeln!(w)?;
        n += 1;
    }
    w.flush()?;
    eprintln!("MSA: wrote {n} aligned rows -> {}", path.display());
    Ok(())
}

fn emit_consensus(outdir: &std::path::Path, seed_name: &str, consensus: &[u8]) -> Result<()> {
    let path = outdir.join("consensus.fa");
    let mut w = BufWriter::new(File::create(&path).with_context(|| format!("create {}", path.display()))?);
    writeln!(w, ">{seed_name}_refined")?;
    w.write_all(consensus)?;
    writeln!(w)?;
    w.flush()?;
    Ok(())
}

fn shannon_entropy(counts: &[u32; 5], total: u32) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let mut e = 0.0f64;
    for &c in counts {
        if c == 0 {
            continue;
        }
        let p = c as f64 / total as f64;
        e -= p * p.log2();
    }
    e
}

fn emit_columns(
    outdir: &std::path::Path,
    outcomes: &[MonomerOutcome],
    consensus: &[u8],
) -> Result<()> {
    // Recompute the tally against the FINAL consensus so the column counts
    // describe the converged alignment rather than the previous round's.
    let mut tally = ColumnTally::new(consensus.len());
    for o in outcomes {
        update_tally(&mut tally, o, consensus);
    }

    let path = outdir.join("columns.tsv");
    let mut w = BufWriter::new(File::create(&path).with_context(|| format!("create {}", path.display()))?);
    writeln!(w, "canonical_col\tA\tC\tG\tT\tgap\tentropy\tindel_rate")?;
    for (col, c) in tally.counts.iter().enumerate() {
        let total: u32 = c.iter().sum();
        let entropy = shannon_entropy(c, total);
        let indel_rate = if total > 0 {
            c[4] as f64 / total as f64
        } else {
            0.0
        };
        writeln!(
            w,
            "{col}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            c[0], c[1], c[2], c[3], c[4], entropy, indel_rate
        )?;
    }
    w.flush()?;
    Ok(())
}

fn emit_summary(
    outdir: &std::path::Path,
    n_std: usize,
    n_out: usize,
    n_exc: usize,
    cause_counts: &BTreeMap<&'static str, usize>,
) -> Result<()> {
    let path = outdir.join("summary.tsv");
    let mut w = BufWriter::new(File::create(&path).with_context(|| format!("create {}", path.display()))?);
    writeln!(w, "stream\tcount")?;
    writeln!(w, "standard\t{n_std}")?;
    writeln!(w, "outlier\t{n_out}")?;
    writeln!(w, "exception\t{n_exc}")?;
    writeln!(w)?;
    writeln!(w, "cause\tcount")?;
    for (k, v) in cause_counts {
        writeln!(w, "{k}\t{v}")?;
    }
    w.flush()?;
    Ok(())
}
