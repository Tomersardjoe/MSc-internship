#!/usr/bin/env python3
"""
Concatenates ribosomal protein sequences across mutiple fasta files into a single fasta file per genome.
"""

import os
import re
import csv
import argparse
from Bio import SeqIO

RIBOSOMAL_PROTEINS = [
    "L2", "L3", "L4", "L5", "L6", "L14", "L16", "L18",
    "L22", "L24", "S3", "S8", "S10", "S17", "S19"
]


def clean_field(value):
    return str(value).strip().replace(" ", "_").replace(";", "_")


def load_taxonomy(tsv_file):
    """Return dict mapping Genome ID -> taxonomy header."""
    mapping = {}

    pattern = re.compile(r"^GCF_.*\.*$")

    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            genome = row.get("Genome", "").strip()

            if not genome:
                continue

            if not pattern.match(genome):
                continue

            class_  = clean_field(row.get("Class", ""))
            order   = clean_field(row.get("Order", ""))
            family  = clean_field(row.get("Family", ""))
            genus   = clean_field(row.get("Genus", ""))
            species = clean_field(row.get("Species", ""))

            header = (
                f"Bacteria_Proteobacteria_"
                f"{class_}_{order}_{family}_{genus}_{species}"
            )

            mapping[genome] = header

    return mapping


def concatenate_proteins(input_dir, output_file, taxonomy_file):
    genome_to_header = load_taxonomy(taxonomy_file)

    protein_data = {}
    genome_ids = set()

    genome_pattern = re.compile(r"^(GCF_\d+\.\d+)")

    for rp in RIBOSOMAL_PROTEINS:

        faa_file = os.path.join(input_dir, f"{rp}.faa")

        if not os.path.exists(faa_file):
            print(f"Missing file: {faa_file}, skipping.")
            continue

        print(f"Loading {faa_file}")

        rp_dict = {}

        for record in SeqIO.parse(faa_file, "fasta"):

            match = genome_pattern.match(record.id)

            if not match:
                continue

            genome_id = match.group(1)

            if genome_id in rp_dict:
                continue

            rp_dict[genome_id] = str(record.seq)
            genome_ids.add(genome_id)

        protein_data[rp] = rp_dict

    genome_ids = sorted(genome_ids)

    print(f"Found {len(genome_ids)} genomes.")

    missing_taxonomy = set()
    written = 0

    with open(output_file, "w") as out_fasta:

        for genome in genome_ids:

            concatenated = ""
            present = 0

            for rp in RIBOSOMAL_PROTEINS:
                if rp in protein_data and genome in protein_data[rp]:
                    concatenated += protein_data[rp][genome]
                    present += 1

            if not concatenated:
                continue

            if present < len(RIBOSOMAL_PROTEINS):
                print(
                    f"Warning: {genome} contains "
                    f"{present}/{len(RIBOSOMAL_PROTEINS)} proteins."
                )

            if genome not in genome_to_header:
                missing_taxonomy.add(genome)
                continue

            header = genome_to_header[genome]

            out_fasta.write(f">{header}\n{concatenated}\n")
            written += 1

    if missing_taxonomy:
        print(
            f"\nExcluded {len(missing_taxonomy)} genomes (missing taxonomy):"
        )
        for g in sorted(missing_taxonomy):
            print(g)

    print(
        f"\nWrote {written} concatenated genomes to {output_file} "
        f"(excluded {len(missing_taxonomy)} missing taxonomy)"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate ribosomal proteins by genome."
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing ribosomal protein FASTA files"
    )

    parser.add_argument(
        "--taxonomy",
        required=True,
        help="Taxonomy TSV file"
    )

    parser.add_argument(
        "--output",
        default="concatenated_proteins.faa",
        help="Output FASTA"
    )

    args = parser.parse_args()

    concatenate_proteins(
        args.input_dir,
        args.output,
        args.taxonomy
    )


if __name__ == "__main__":
    main()