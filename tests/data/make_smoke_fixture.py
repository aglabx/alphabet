#!/usr/bin/env python3
"""Deterministic ground-truth smoke fixture for anchor-scaffold MSA.

Writes fixtures/{smoke.fasta, smoke_consensus.fa, smoke_invariants.json}
next to this script. Re-running yields byte-identical output (fixed seed).

See docs/anchor_scaffold_msa_brief.md for the fixture design.
"""
import json
import os
import random

L = 171
ANCHORS = {
    0: ("boundary", "CTTTGTGATGT"),
    14: ("slot0", "CATTC"),
    24: ("slot1", "CAGAG"),
    40: ("slot2", "CTTTT"),
    59: ("slot3", "GAAAC"),
    86: ("slot4", "GTGGA"),
    117: ("slot5", "TGGAA"),
    141: ("slot6", "AAACT"),
    148: ("slot7", "ACAGA"),
    155: ("slot8", "CATTC"),
    161: ("slot9", "CAGAA"),
}
N_PRIMARY, MIN_ANCHORS = 10, 6  # mirror ASM-msa-parameter-set


def build_consensus():
    seq = ["ATGC"[i % 4] for i in range(L)]  # no CpG, no homopolymer
    for pos, (_l, kmer) in ANCHORS.items():
        for j, b in enumerate(kmer):
            seq[pos + j] = b
    return seq


def primary_positions():
    return [p for p, (l, _k) in ANCHORS.items() if l != "boundary"]


def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
    os.makedirs(out, exist_ok=True)
    rng = random.Random(20260622)
    cons = build_consensus()
    cons_s = "".join(cons)
    mon = []

    # clean x6
    for i in range(1, 7):
        mon.append((f"M{i:02d}_clean", cons_s, "standard", "exact; D=0"))

    # substitution x4
    for k, p in enumerate([12, 33, 70, 130], 7):
        s = list(cons)
        old = s[p]
        s[p] = {"A": "T", "T": "G", "G": "A", "C": "A"}[old]
        mon.append((f"M{k:02d}_sub{p}", "".join(s), "standard", f"sub col {p} {old}->{s[p]}"))

    idx = 11
    # deletion x3
    for p in [12, 70, 130]:
        mon.append((f"M{idx:02d}_del{p}", "".join(cons[:p] + cons[p + 1 :]), "standard", f"del col {p}"))
        idx += 1

    # insertion x2
    for p, b in [(33, "C"), (130, "A")]:
        mon.append((f"M{idx:02d}_ins{p}", "".join(cons[:p] + [b] + cons[p:]), "standard", f"ins before {p}"))
        idx += 1

    # anchor-knockout x3 (replace with ATGC filler, not poly-A)
    for ko in (86, 117, 141):
        s = list(cons)
        for j in range(len(ANCHORS[ko][1])):
            s[ko + j] = "ATGC"[(ko + j) % 4]
        mon.append((f"M{idx:02d}_knockout{ko}", "".join(s), "standard", f"anchor {ko} lost; across-gap"))
        idx += 1

    # outlier x3
    for _ in range(3):
        s = list(cons)
        for p in range(91, 117):
            s[p] = rng.choice("ACGT")
        mon.append((f"M{idx:02d}_outlier", "".join(s), "outlier", "seg 91-116 random; D>band"))
        idx += 1

    # exception x3 (each leaves exactly 4 primary anchors intact)
    prim = primary_positions()
    for ds in [prim[:6], prim[2:8], prim[4:]]:
        s = list(cons)
        for pos in ds:
            for j in range(len(ANCHORS[pos][1])):
                s[pos + j] = "ATGC"[(pos + j) % 4]
        mon.append((f"M{idx:02d}_exception", "".join(s), "exception", f"4 present < {MIN_ANCHORS}"))
        idx += 1

    with open(os.path.join(out, "smoke.fasta"), "w") as fh:
        for mid, s, _st, _n in mon:
            fh.write(f">{mid}\n{s}\n")
    with open(os.path.join(out, "smoke_consensus.fa"), "w") as fh:
        fh.write(f">consensus_171bp\n{cons_s}\n")

    cnt = {}
    for _i, _s, st, _n in mon:
        cnt[st] = cnt.get(st, 0) + 1
    with open(os.path.join(out, "smoke_invariants.json"), "w") as fh:
        json.dump(
            {
                "n_input": len(mon),
                "consensus_len": L,
                "min_anchors": MIN_ANCHORS,
                "stream_counts": cnt,
                "monomers": [
                    {"id": m, "expected_stream": st, "note": nt, "len": len(s)}
                    for m, s, st, nt in mon
                ],
            },
            fh,
            indent=2,
        )

    assert len(mon) == 24 and cnt == {"standard": 18, "outlier": 3, "exception": 3}
    print("n_input", len(mon), "counts", cnt)


if __name__ == "__main__":
    main()
