#!/usr/bin/env python3
"""
Extract ribosomal protein amino acid sequences from multiple GenBank (.gbff) files
and write 16 per-protein FASTA files containing sequences from all genomes.
Takes the first sequence found for each ribosomal protein in each genome.
Prints a warning if a genome is missing any ribosomal proteins.
"""

import os
from Bio import SeqIO

# List of ribosomal proteins to extract
RIBOSOMAL_PROTEINS = [
    "L2","L3","L4","L5","L6","L14","L16","L18","L22","L24",
    "S3","S8","S10","S17","S19"
]

def extract_all_ribosomal_proteins(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Dictionary to store sequences for each ribosomal protein
    rp_sequences = {rp: [] for rp in RIBOSOMAL_PROTEINS}

    # Find all .gbff files
    gbff_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".gbff")]
    print(f"Found {len(gbff_files)} .gbff files.")

    for gbff_file in gbff_files:
        genome_id = os.path.basename(gbff_file).replace(".gbff", "")
        seen_rps = set()  # Track which ribosomal proteins have been found for this genome

        for record in SeqIO.parse(gbff_file, "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers
                gene_name = qualifiers.get("gene", [""])[0].upper()
                product_name = qualifiers.get("product", [""])[0].upper()
                translation = qualifiers.get("translation", [""])[0]

                for rp in RIBOSOMAL_PROTEINS:
                    if rp in seen_rps:
                        continue  # Already extracted this protein for this genome
                    if rp.upper() == gene_name or f"RIBOSOMAL PROTEIN {rp}" in product_name:
                        if translation:
                            header = f">{genome_id}_{rp}"
                            rp_sequences[rp].append(f"{header}\n{translation}")
                            seen_rps.add(rp)
                        break

        # Print missing ribosomal proteins for this genome
        missing_rps = set(RIBOSOMAL_PROTEINS) - seen_rps
        if missing_rps:
            print(f"Warning: genome {genome_id} missing ribosomal proteins: {', '.join(sorted(missing_rps))}")

        print(f"Processed genome {genome_id}")

    # Write per-ribosomal protein FASTA files
    for rp, sequences in rp_sequences.items():
        if sequences:
            out_path = os.path.join(output_dir, f"{rp}.faa")
            with open(out_path, "w") as out_fasta:
                out_fasta.write("\n".join(sequences) + "\n")
            print(f"Wrote {len(sequences)} sequences to {out_path}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract ribosomal protein sequences from GenBank files into per-protein FASTAs.")
    parser.add_argument("-i", "--input", required=True, help="Directory containing .gbff files")
    parser.add_argument("-o", "--output", required=True, help="Output directory for per-protein FASTA files")
    args = parser.parse_args()

    extract_all_ribosomal_proteins(args.input, args.output)

if __name__ == "__main__":
    main()