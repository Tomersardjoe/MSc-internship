#!/usr/bin/env python3
"""
Concatenate trimmed ribosomal protein alignments into a single supermatrix.
Missing proteins for a genome are replaced with gap sequences.
Order of proteins follows the Hug et al. (10.1038/nmicrobiol.2016.48) paper.
Genome IDs are replaced with taxonomy-derived headers.
"""

import os
import csv
import re
import warnings
from Bio import SeqIO

# Order of ribosomal proteins for concatenation
RIBOSOMAL_PROTEINS = [
    "L2","L3","L4","L5","L6","L14","L16","L18","L22","L24",
    "S3","S8","S10","S17","S19"
]

def clean_field(x):
    """Replace spaces with underscores and handle missing values."""
    if x is None or x.strip() == "" or x == "NA":
        return "Unknown"
    return x.strip().replace(" ", "_")

def load_taxonomy(tsv_file):
    """Return dict mapping Genome ID -> formatted taxonomy header."""
    mapping = {}

    pattern = re.compile(r"^GCF_.*\.*$")

    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader, start=1):
            genome = row.get("Genome", "").strip()

            # Warn if Genome is missing or malformed
            if not genome:
                warnings.warn(
                    f"Missing Genome value at line {i}. Row will be skipped.",
                    RuntimeWarning
                )
                continue

            if not pattern.match(genome):
                warnings.warn(
                    f"Unexpected Genome format at line {i}: '{genome}'. "
                    "Expected pattern 'GCF_*.1'.",
                    RuntimeWarning
                )
                continue

            class_  = clean_field(row.get("Class", ""))
            order   = clean_field(row.get("Order", ""))
            family  = clean_field(row.get("Family", ""))
            genus   = clean_field(row.get("Genus", ""))
            species = clean_field(row.get("Species", ""))

            header = f"Bacteria_Proteobacteria_{class_}_{order}_{family}_{genus}_{species}"
            mapping[genome] = header

    return mapping

def concatenate_proteins(trimmed_dir, output_file, taxonomy_file):
    # Load taxonomy
    genome_to_header = load_taxonomy(taxonomy_file)

    # Load all alignments into a dict of dicts: {protein: {genome_id: seq}}
    alignments = {}
    genome_ids = set()

    for rp in RIBOSOMAL_PROTEINS:
        rp_file = os.path.join(trimmed_dir, f"{rp}_trimmed.faa")
        if not os.path.exists(rp_file):
            print(f"Warning: trimmed alignment {rp_file} not found, skipping.")
            continue

        rp_dict = {}
        for record in SeqIO.parse(rp_file, "fasta"):
            genome_id = "_".join(record.id.split("_")[:2])
            seq = str(record.seq).replace("\n", "")
            rp_dict[genome_id] = seq
            genome_ids.add(genome_id)

        alignments[rp] = rp_dict

    genome_ids = sorted(genome_ids)
    print(f"Found {len(genome_ids)} genomes across all proteins.")

    # Determine alignment lengths for each protein
    rp_lengths = {
        rp: len(next(iter(seq_dict.values())))
        for rp, seq_dict in alignments.items()
    }

    # Concatenate for each genome
    with open(output_file, "w") as out_fasta:
        header_counter = {}

        for genome in genome_ids:
            concatenated_seq = ""

            for rp in RIBOSOMAL_PROTEINS:
                if rp in alignments and genome in alignments[rp]:
                    concatenated_seq += alignments[rp][genome]
                elif rp in rp_lengths:
                    concatenated_seq += "-" * rp_lengths[rp]
                else:
                    print(f"Warning: {rp} missing entirely, skipping gaps for {genome}")

            # Get taxonomy-based header
            header = genome_to_header.get(genome, genome)
            count = header_counter.get(header, 0) + 1
            header_counter[header] = count
            if count > 1:
                header = f"{header}_{count}"

            out_fasta.write(f">{header}\n{concatenated_seq}\n")

    print(f"Supermatrix written to {output_file}")

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Concatenate trimmed ribosomal protein alignments into a supermatrix"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Directory containing trimmed per-protein alignments")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA file for the concatenated supermatrix")
    parser.add_argument("-t", "--taxonomy", required=True,
                        help="Path to taxonomy TSV file")
    args = parser.parse_args()

    concatenate_proteins(args.input, args.output, args.taxonomy)

if __name__ == "__main__":
    main()