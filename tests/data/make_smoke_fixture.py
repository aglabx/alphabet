#!/usr/bin/env python3
"""Deterministic ground-truth smoke fixture for EXP-anchor-scaffold-msa.

Builds a synthetic 171 bp canonical consensus with the 15-slot V5.1 anchors at
their canonical positions, then derives monomers that each exercise ONE code
path of the anchor-scaffold MSA. Every monomer's expected output stream/column
is determined by construction, so the programmer's smoke test can assert exact
invariants (not just "runs without error").

NOTE: synthetic test fixture for CODE correctness only. Filler is an ATGC cycle
(no CpG, no homopolymers) so planted SNP/indel columns are unambiguous. Not a
biological consensus. Deterministic: fixed seed; identical output every run.

Outputs (under this script's dir):
  fixtures/smoke.fasta            - 24 monomers
  fixtures/smoke_consensus.fa     - the synthetic 171 bp consensus
  fixtures/smoke_invariants.json  - expected per-monomer stream + planted edits
"""
import json
import os
import random

L = 171
ANCHORS = {  # canonical_pos -> (slot_label, kmer)  (primary slots + boundary)
    0:   ("boundary", "CTTTGTGATGT"),
    14:  ("slot0", "CATTC"),
    24:  ("slot1", "CAGAG"),
    40:  ("slot2", "CTTTT"),
    59:  ("slot3", "GAAAC"),
    86:  ("slot4", "GTGGA"),
    117: ("slot5", "TGGAA"),
    141: ("slot6", "AAACT"),
    148: ("slot7", "ACAGA"),
    155: ("slot8", "CATTC"),
    161: ("slot9", "CAGAA"),
}
N_PRIMARY = 10           # boundary excluded from the min_anchors count
MIN_ANCHORS = 6          # mirror ASM-msa-parameter-set (>=6 of 10 primary)


def build_consensus():
    seq = ["ATGC"[i % 4] for i in range(L)]   # no CpG, no homopolymer
    for pos, (_lbl, kmer) in ANCHORS.items():
        for j, b in enumerate(kmer):
            seq[pos + j] = b
    return seq


def primary_anchor_positions():
    return [p for p, (lbl, _k) in ANCHORS.items() if lbl != "boundary"]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "fixtures")
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(20260622)            # fixture-determinism only

    cons = build_consensus()
    cons_s = "".join(cons)
    assert len(cons_s) == L

    monomers = []   # (id, seq, expected_stream, note)

    # --- 6 clean copies: standard, 0 edits ---
    for i in range(1, 7):
        monomers.append((f"M{i:02d}_clean", cons_s, "standard", "exact consensus; 0 edits"))

    # --- 4 single substitutions in inter-anchor filler (unambiguous columns) ---
    sub_sites = [12, 33, 70, 130]
    for k, p in enumerate(sub_sites, start=7):
        s = list(cons); old = s[p]
        s[p] = {"A": "T", "T": "G", "G": "A", "C": "A"}[old]
        monomers.append((f"M{k:02d}_sub{p}", "".join(s), "standard",
                         f"substitution col {p}: {old}->{s[p]} (D=1)"))

    # --- 3 single-base deletions at unique (non-homopolymer) filler sites ---
    del_sites = [12, 70, 130]
    idx = 11
    for p in del_sites:
        s = cons[:p] + cons[p + 1:]
        monomers.append((f"M{idx:02d}_del{p}", "".join(s), "standard",
                         f"1bp deletion at col {p} (len {L-1}); gap unambiguous"))
        idx += 1

    # --- 2 single-base insertions in filler (inserted base != both neighbours) ---
    ins_sites = [(33, "C"), (130, "A")]   # 130: left=T(129) right=G(130) -> A unambiguous
    for p, b in ins_sites:
        s = cons[:p] + [b] + cons[p:]
        monomers.append((f"M{idx:02d}_ins{p}", "".join(s), "standard",
                         f"1bp insertion before col {p} (len {L+1})"))
        idx += 1

    # --- 3 anchor-knockouts: forces across-missing-anchor O(ND) path, stays standard ---
    # overwrite with ATGC filler (NOT poly-A) so no homopolymer is introduced
    for ko_pos in (86, 117, 141):                  # knock out one anchor each
        s = list(cons)
        for j in range(len(ANCHORS[ko_pos][1])):
            s[ko_pos + j] = "ATGC"[(ko_pos + j) % 4]   # anchor lost, filler-consistent
        monomers.append((f"M{idx:02d}_knockout{ko_pos}", "".join(s), "standard",
                         f"anchor at {ko_pos} destroyed (->filler); 9/10 primary present; "
                         f"gap aligned across missing anchor (D<=band)"))
        idx += 1

    # --- 3 outliers: one inter-anchor segment fully randomized (D>band), anchors intact ---
    for _ in range(3):
        s = list(cons)
        for p in range(91, 117):                    # segment slot4(90)..slot5(117)
            s[p] = rng.choice("ACGT")
        monomers.append((f"M{idx:02d}_outlier", "".join(s), "outlier",
                         "segment 91-116 randomized; D>band on that segment; "
                         "all 10 primary anchors present -> outliers stream"))
        idx += 1

    # --- 3 exceptions: destroy 6 of 10 primary anchors -> 4 present < min_anchors (guaranteed) ---
    prim = primary_anchor_positions()               # 10 canonical positions
    destroy_sets = [prim[:6], prim[2:8], prim[4:]]  # each leaves exactly 4 intact
    for ds in destroy_sets:
        s = list(cons)
        for pos in ds:                              # revert anchor span to ATGC filler
            for j in range(len(ANCHORS[pos][1])):
                s[pos + j] = "ATGC"[(pos + j) % 4]
        monomers.append((f"M{idx:02d}_exception", "".join(s), "exception",
                         f"{len(ds)} primary anchors reverted to filler; "
                         f"4 present < {MIN_ANCHORS} -> exceptions stream"))
        idx += 1

    # write FASTA
    with open(os.path.join(outdir, "smoke.fasta"), "w") as fh:
        for mid, s, _stream, _note in monomers:
            fh.write(f">{mid}\n{s}\n")
    with open(os.path.join(outdir, "smoke_consensus.fa"), "w") as fh:
        fh.write(f">synthetic_canonical_consensus_171bp\n{cons_s}\n")

    # invariants
    counts = {}
    for _id, _s, stream, _n in monomers:
        counts[stream] = counts.get(stream, 0) + 1
    inv = {
        "n_input": len(monomers),
        "consensus_len": L,
        "min_anchors": MIN_ANCHORS,
        "n_primary_slots": N_PRIMARY,
        "stream_counts": counts,
        "monomers": [{"id": mid, "expected_stream": st, "note": nt, "len": len(s)}
                     for (mid, s, st, nt) in monomers],
    }
    with open(os.path.join(outdir, "smoke_invariants.json"), "w") as fh:
        json.dump(inv, fh, indent=2)

    # structural self-checks (what we CAN verify without the pipeline)
    assert len(monomers) == 24
    assert counts == {"standard": 18, "outlier": 3, "exception": 3}
    assert sum(counts.values()) == 24
    print("consensus (171bp):", cons_s)
    print("n_input:", inv["n_input"], "stream_counts:", counts)
    print("wrote:", outdir)


if __name__ == "__main__":
    main()
