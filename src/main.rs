use clap::{Parser, Subcommand};
use std::path::PathBuf;

use alphabet::cmd;

#[derive(Parser)]
#[command(
    name = "alphabet",
    version,
    about = "Ab initio satellite DNA monomer alphabet discovery",
    arg_required_else_help = true
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// End-to-end pipeline: discover → cut → annotate
    Run {
        /// Input .10kb.fasta of satellite arrays
        input: PathBuf,
        /// Output directory (chains.json, monomers.tsv, annotated.tsv)
        #[arg(short = 'o', long, default_value = ".")]
        outdir: PathBuf,
        /// Threads (0 = rayon default)
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },

    /// Chain-first anchor discovery (periodic motif chains)
    Discover {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Cut arrays into monomers at chain-site boundaries; classify by letter
    Cut {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Multi-anchor position-aware monomer cutter (alpha satellite v2)
    #[command(name = "cut-v2")]
    CutV2 {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Annotate monomer TSV with CENP-B box columns
    Annotate {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },

    /// Generic box pattern search (e.g. TIGD4/CENP-B)
    FindBox {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// CENP-B box spacing analysis
    Spacing {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// De novo periodic box discovery
    FindPeriodicBoxes {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },

    /// ONT reads pipeline (HPC-HMM classification)
    #[command(subcommand)]
    Reads(ReadsCmd),

    /// Experimental / research subcommands (unstable)
    #[command(subcommand, hide = true)]
    Dev(DevCmd),
}

#[derive(Subcommand)]
enum ReadsCmd {
    /// Build HPC-HMMs from an annotated TSV
    BuildHmms {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Scan raw reads for chain-site motifs
    Scan {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Classify satellite-containing reads against per-letter HMMs
    Classify {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Extract satellite-containing reads
    Extract {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Per-read HPC alphabet stats
    Alphabet {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
}

#[derive(Subcommand)]
enum DevCmd {
    /// Phase optimization experiment (negative result, see labjournal E2)
    FindPhase {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
    /// Motif transition graph (exploratory)
    MotifGraph {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },
}

fn prepend_prog(subcmd: &str, args: Vec<String>) -> Vec<String> {
    std::iter::once(format!("alphabet {}", subcmd))
        .chain(args)
        .collect()
}

fn main() {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Run { input, outdir, threads } => run_pipeline(input, outdir, threads),

        Cmd::Discover { args } => cmd::discover_chains::run_from_args(prepend_prog("discover", args)),
        Cmd::Cut { args } => cmd::motif_cut::run_from_args(prepend_prog("cut", args)),
        Cmd::CutV2 { args } => cmd::cut_v2::run_from_args(prepend_prog("cut-v2", args)),
        Cmd::Annotate { args } => cmd::annotate_cenpb::run_from_args(prepend_prog("annotate", args)),

        Cmd::FindBox { args } => cmd::find_box::run_from_args(prepend_prog("find-box", args)),
        Cmd::Spacing { args } => cmd::cenpb_spacing::run_from_args(prepend_prog("spacing", args)),
        Cmd::FindPeriodicBoxes { args } => {
            cmd::find_periodic_boxes::run_from_args(prepend_prog("find-periodic-boxes", args))
        }

        Cmd::Reads(r) => match r {
            ReadsCmd::BuildHmms { args } => {
                cmd::build_hpc_hmms::run_from_args(prepend_prog("reads build-hmms", args))
            }
            ReadsCmd::Scan { args } => cmd::scan_reads::run_from_args(prepend_prog("reads scan", args)),
            ReadsCmd::Classify { args } => {
                cmd::classify_reads::run_from_args(prepend_prog("reads classify", args))
            }
            ReadsCmd::Extract { args } => {
                cmd::reads_extract::run_from_args(prepend_prog("reads extract", args))
            }
            ReadsCmd::Alphabet { args } => {
                cmd::reads_alphabet::run_from_args(prepend_prog("reads alphabet", args))
            }
        },

        Cmd::Dev(d) => match d {
            DevCmd::FindPhase { args } => {
                cmd::find_phase::run_from_args(prepend_prog("dev find-phase", args))
            }
            DevCmd::MotifGraph { args } => {
                cmd::motif_graph::run_from_args(prepend_prog("dev motif-graph", args))
            }
        },
    }
}

fn run_pipeline(input: PathBuf, outdir: PathBuf, threads: usize) {
    std::fs::create_dir_all(&outdir).expect("create outdir");
    let chains = outdir.join("chains.json");
    let monomers = outdir.join("monomers.tsv");
    let annotated = outdir.join("annotated.tsv");
    let threads_s = threads.to_string();
    let input_s = input.to_string_lossy().into_owned();
    let chains_s = chains.to_string_lossy().into_owned();
    let monomers_s = monomers.to_string_lossy().into_owned();
    let annotated_s = annotated.to_string_lossy().into_owned();

    eprintln!("[1/3] discover → {}", chains_s);
    cmd::discover_chains::run_from_args(vec![
        "alphabet discover".into(),
        input_s.clone(),
        "-o".into(), chains_s.clone(),
        "-t".into(), threads_s.clone(),
    ]);

    eprintln!("[2/3] cut → {}", monomers_s);
    cmd::motif_cut::run_from_args(vec![
        "alphabet cut".into(),
        input_s,
        "-m".into(), chains_s,
        "-o".into(), monomers_s.clone(),
        "-t".into(), threads_s,
    ]);

    eprintln!("[3/3] annotate → {}", annotated_s);
    cmd::annotate_cenpb::run_from_args(vec![
        "alphabet annotate".into(),
        monomers_s,
        annotated_s,
    ]);
}
