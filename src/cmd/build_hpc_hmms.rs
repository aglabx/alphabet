use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use crate::monomer::hpc;

/// Step 1: Extract HPC monomers per letter from CHM13 annotated TSV.
/// Outputs one FASTA per letter (top N letters by count).
///
/// Usage: build_hpc_hmms <annotated.tsv> <outdir> [max_letters=30] [max_per_letter=200]
pub fn run_from_args(args: Vec<String>) {
    if args.len() < 3 {
        eprintln!("Usage: build_hpc_hmms <annotated.tsv> <outdir> [max_letters=30] [max_per_letter=200]");
        std::process::exit(1);
    }
    let tsv_path = &args[1];
    let outdir = &args[2];
    let max_letters: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(30);
    let max_per_letter: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(200);

    std::fs::create_dir_all(outdir).unwrap();

    // Pass 1: count letters
    let file = std::fs::File::open(tsv_path).unwrap();
    let reader = BufReader::new(file);
    let mut letter_counts: HashMap<String, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 { continue; }
        let letter = fields[4]; // column 5 = our letter
        if letter == "letter" { continue; } // header
        *letter_counts.entry(letter.to_string()).or_insert(0) += 1;
    }

    // Sort by count, take top N
    let mut sorted: Vec<_> = letter_counts.into_iter().collect();
    sorted.sort_by_key(|b| std::cmp::Reverse(b.1));
    sorted.truncate(max_letters);

    let target_letters: HashMap<String, usize> = sorted.iter().cloned().collect();
    eprintln!("Top {} letters:", sorted.len());
    for (letter, count) in &sorted {
        eprintln!("  {}: {}", letter, count);
    }

    // Pass 2: extract sequences, HPC compress, write per-letter FASTA
    let file = std::fs::File::open(tsv_path).unwrap();
    let reader = BufReader::new(file);
    let mut per_letter: HashMap<String, Vec<(String, String)>> = HashMap::new(); // letter -> [(name, hpc_seq)]

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 { continue; }
        let letter = fields[4];
        if !target_letters.contains_key(letter) { continue; }

        let seq = fields[fields.len() - 1]; // last column = sequence
        // Filter: only canonical-length monomers (150-190bp), skip dimers and fragments
        if seq.len() < 150 || seq.len() > 190 { continue; }

        let entry = per_letter.entry(letter.to_string()).or_default();
        if entry.len() >= max_per_letter { continue; }

        let hpc_seq = hpc(&seq.bytes().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>());
        let name = if fields.len() >= 3 {
            format!("{}_{}", fields[0], fields[1]) // chr_start
        } else {
            format!("{}_{}", letter, entry.len())
        };
        entry.push((name, String::from_utf8(hpc_seq).unwrap()));
    }

    // Write FASTA files
    for (letter, seqs) in &per_letter {
        let path = format!("{}/letter_{}.fasta", outdir, letter);
        let mut f = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
        for (name, seq) in seqs {
            writeln!(f, ">{}", name).unwrap();
            writeln!(f, "{}", seq).unwrap();
        }
        eprintln!("  {} → {} ({} seqs, avg {}bp HPC)",
            letter, path, seqs.len(),
            seqs.iter().map(|(_, s)| s.len()).sum::<usize>() / seqs.len().max(1));
    }

    // Write build script
    let script_path = format!("{}/build_hmms.sh", outdir);
    let mut script = std::io::BufWriter::new(std::fs::File::create(&script_path).unwrap());
    writeln!(script, "#!/bin/bash").unwrap();
    writeln!(script, "set -euo pipefail").unwrap();
    writeln!(script, "cd \"$(dirname \"$0\")\"").unwrap();
    writeln!(script).unwrap();

    for (letter, _) in &sorted {
        if !per_letter.contains_key(letter) { continue; }
        writeln!(script, "echo \"Building HMM for letter {}\"", letter).unwrap();
        writeln!(script, "muscle -align letter_{}.fasta -output letter_{}.afa 2>/dev/null", letter, letter).unwrap();
        writeln!(script, "hmmbuild --dna letter_{}.hmm letter_{}.afa >/dev/null 2>&1", letter, letter).unwrap();
        writeln!(script).unwrap();
    }

    // Concat all HMMs
    writeln!(script, "echo \"Concatenating HMMs...\"").unwrap();
    writeln!(script, "cat letter_*.hmm > all_letters.hmm").unwrap();
    writeln!(script, "echo \"Done: $(ls letter_*.hmm | wc -l) HMMs built\"").unwrap();

    eprintln!("\nBuild script: {}", script_path);
    eprintln!("Run: bash {}", script_path);
}

