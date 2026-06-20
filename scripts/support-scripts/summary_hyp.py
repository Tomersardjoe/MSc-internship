#!/usr/bin/env python3

from pathlib import Path
from statistics import median
import csv

SYN_A = "hydrogenase_maturation_nickel_metallochaperone_HypA"
SYN_B = "hydrogenase_nickel_incorporation_protein_HypB"


def summarize(tsv_file):

    alpha_alpha = []
    alpha_syn = []

    seen = set()

    with open(tsv_file) as f:
        for line in f:
            q, s, pid, *_ = line.rstrip().split("\t")
            pid = float(pid)

            # remove reciprocal duplicates
            pair = tuple(sorted((q, s)))
            if pair in seen:
                continue
            seen.add(pair)

            is_syn = (
                SYN_A in q or SYN_A in s or
                SYN_B in q or SYN_B in s
            )

            if is_syn:
                alpha_syn.append(pid)
            else:
                alpha_alpha.append(pid)

    result = {
        "family": Path(tsv_file).stem.replace("_pairs", ""),
    }

    if alpha_alpha:
        result.update({
            "alpha_alpha_n": len(alpha_alpha),
            "alpha_alpha_min": min(alpha_alpha),
            "alpha_alpha_median": median(alpha_alpha),
            "alpha_alpha_max": max(alpha_alpha),
        })

    if alpha_syn:
        result.update({
            "alpha_syn_n": len(alpha_syn),
            "alpha_syn_min": min(alpha_syn),
            "alpha_syn_median": median(alpha_syn),
            "alpha_syn_max": max(alpha_syn),
        })

    return result


rows = []

for tsv in sorted(Path(".").glob("hypothetical_*_pairs.tsv")):

    row = summarize(tsv)
    rows.append(row)

    print(f"\n{row['family']}")

    print(
        f"  Alpha vs Alpha: "
        f"{row['alpha_alpha_min']:.1f}-"
        f"{row['alpha_alpha_max']:.1f}% "
        f"(median {row['alpha_alpha_median']:.1f}%)"
    )

    print(
        f"  Alpha vs Synechocystis: "
        f"{row['alpha_syn_min']:.1f}-"
        f"{row['alpha_syn_max']:.1f}% "
        f"(median {row['alpha_syn_median']:.1f}%)"
    )

with open("hypothetical_identity_summary.tsv", "w", newline="") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=[
            "family",
            "alpha_alpha_n",
            "alpha_alpha_min",
            "alpha_alpha_median",
            "alpha_alpha_max",
            "alpha_syn_n",
            "alpha_syn_min",
            "alpha_syn_median",
            "alpha_syn_max",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

print("\nWrote hypothetical_identity_summary.tsv")