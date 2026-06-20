#!/usr/bin/env python3
"""
Filter an already-merged ribosomal supermatrix for Alphaproteobacteria
and convert outputs to single-line FASTA.
"""

import argparse
import subprocess
from Bio import SeqIO


def filter_alphaproteo(input_fasta, output_fasta):
    """
    Keep only sequences whose headers:
    - start with 'Bacteria_'
    - contain 'Alphaproteobacteria'
    """
    print(f"Filtering Alphaproteobacteria: {input_fasta} -> {output_fasta}")

    kept = 0
    total = 0

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            if record.id.startswith("Bacteria_") and "Alphaproteobacteria" in record.id:
                SeqIO.write(record, out_handle, "fasta")
                kept += 1

    print(f"Kept {kept}/{total} sequences")


def convert_to_single_line(input_fasta, output_fasta):
    """Convert FASTA to single-line sequences using seqkit."""
    print(f"Converting to single-line FASTA: {output_fasta}")

    with open(output_fasta, "w") as out_handle:
        subprocess.run(
            ["seqkit", "seq", "-w", "0", input_fasta],
            stdout=out_handle,
            check=True
        )

    print(f"Wrote {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter merged ribosomal supermatrix and convert format"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Already merged MAFFT supermatrix FASTA"
    )

    args = parser.parse_args()

    filtered = "filtered_ribo_alignment.faa"
    singleline = "filtered_ribo_alignment_singleline.faa"

    # filter
    filter_alphaproteo(args.input, filtered)

    # convert formatting
    convert_to_single_line(filtered, singleline)


if __name__ == "__main__":
    main()