#!/usr/bin/env python3
"""
Align ribosomal protein FASTA files using MAFFT, optionally trim with TrimAl.
Sequences in output FASTA will be on a single line each.
"""

import os
import subprocess

# List of ribosomal proteins
RIBOSOMAL_PROTEINS = [
    "L2","L3","L4","L5","L6","L14","L16","L18","L22","L24",
    "S3","S8","S10","S17","S19"
]

def reformat_fasta_single_line(fasta_path):
    """Rewrites a FASTA file so each sequence is on a single line."""
    lines = []
    with open(fasta_path) as f:
        seq = ""
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    lines.append(header)
                    lines.append(seq)
                header = line
                seq = ""
            else:
                seq += line
        if header:
            lines.append(header)
            lines.append(seq)
    with open(fasta_path, "w") as f:
        for i in range(0, len(lines), 2):
            f.write(lines[i] + "\n" + lines[i+1] + "\n")

def align_proteins(input_dir, output_dir, trim_threshold=None):
    os.makedirs(output_dir, exist_ok=True)

    for rp in RIBOSOMAL_PROTEINS:
        infile = os.path.join(input_dir, f"{rp}.faa")
        if not os.path.exists(infile):
            print(f"Warning: {infile} not found, skipping.")
            continue

        aligned_file = os.path.join(output_dir, f"{rp}_aligned.faa")
        trimmed_file = os.path.join(output_dir, f"{rp}_trimmed.faa")

        # Align with MAFFT
        print(f"Aligning {rp} with MAFFT...")
        mafft_out = trimmed_file if trim_threshold is not None else aligned_file
        with open(mafft_out, "w") as out_handle:
            subprocess.run(["mafft", "--auto", infile], stdout=out_handle, check=True)

        # Trim with TrimAl if requested
        if trim_threshold is not None:
            print(f"Trimming {rp} alignment with TrimAl (-gt {trim_threshold})...")
            subprocess.run(["trimal", "-in", mafft_out, "-out", trimmed_file, "-gt", str(trim_threshold)], check=True)
            reformat_fasta_single_line(trimmed_file)
            print(f"{rp} trimmed alignment saved to {trimmed_file}\n")
        else:
            reformat_fasta_single_line(aligned_file)
            print(f"{rp} aligned alignment saved to {aligned_file}\n")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Align (and optionally trim) ribosomal protein FASTA files")
    parser.add_argument("-i", "--input", required=True, help="Directory containing per-protein FASTA files (*.faa)")
    parser.add_argument("-o", "--output", required=True, help="Directory to save output alignments")
    parser.add_argument("-t", "--trim", type=float, help="Trim alignments with TrimAl using this -gt threshold (e.g., 0.05)")
    args = parser.parse_args()

    align_proteins(args.input, args.output, trim_threshold=args.trim)


if __name__ == "__main__":
    main()