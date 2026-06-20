#!/usr/bin/env python3
import os
import subprocess
import argparse
import shutil
import gzip
import csv
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter

# Helper functions
def clean_field(x):
    if x is None or x.strip() == "" or x == "NA":
        return "Unknown"
    return x.strip().replace(" ", "_")

def load_taxonomy(tsv_file):
    mapping = {}
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genome = row["Assembly"]
            class_  = clean_field(row["Class"])
            order   = clean_field(row["Order"])
            family  = clean_field(row["Family"])
            genus   = clean_field(row["Genus"])
            species = clean_field(row["Species"])
            header = f"Bacteria_Proteobacteria_{class_}_{order}_{family}_{genus}_{species}"
            mapping[genome] = header
    return mapping

def resolve_rna_url(gcf):
    """
    Resolve correct FTP URL by first listing parent directory.
    """
    try:
        prefix, rest = gcf.split("_")
        numeric, version = rest.split(".")
        chunks = [numeric[i:i+3] for i in range(0, len(numeric), 3)]

        parent_dir = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/" + "/".join(chunks)

        # List parent directory
        response = requests.get(parent_dir)
        if response.status_code != 200:
            return None

        assembly_dir = None

        # Find matching assembly directory
        for line in response.text.splitlines():
            if gcf in line and 'href="' in line:
                start = line.find('href="') + 6
                end = line.find('"', start)
                assembly_dir = line[start:end].rstrip("/")
                break

        if not assembly_dir:
            return None

        # Build full path
        file_url = f"{parent_dir}/{assembly_dir}/{assembly_dir}_rna_from_genomic.fna.gz"

        return file_url

    except Exception as e:
        print(f"Error resolving FTP path for {gcf}: {e}")

    return None

# Parse arguments
parser = argparse.ArgumentParser(
    description="Download 16S rRNA sequences from RefSeq FTP and extract longest <=2000 bp with taxonomy headers."
)
parser.add_argument("-d", "--parent_dir", required=True)
parser.add_argument("-o", "--output_dir", required=True)
parser.add_argument("-t", "--taxonomy", required=True)
args = parser.parse_args()

parent_dir = args.parent_dir
output_dir = args.output_dir
taxonomy_file = args.taxonomy
os.makedirs(output_dir, exist_ok=True)

# Load taxonomy mapping
genome_to_header = load_taxonomy(taxonomy_file)

# Find GCF directories
subdirs = [d for d in os.listdir(parent_dir) if d.startswith("GCF_") and os.path.isdir(os.path.join(parent_dir, d))]
if not subdirs:
    print("No GCF_* subdirectories found.")
    exit(0)

current_gcfs = set(subdirs)
header_counter = {}

# Main loop
for gcf in subdirs:
    print(f"\nProcessing {gcf}...")
    subdir_path = os.path.join(parent_dir, gcf)
    final_fasta = os.path.join(output_dir, f"{gcf}_16S.fna")

    if os.path.exists(final_fasta):
        print(f"Skipping {gcf}: already processed.")
        continue

    # Resolve FTP URL
    ftp_url = resolve_rna_url(gcf)
    if not ftp_url:
        print(f"ERROR: Could not resolve FTP URL for {gcf}")
        continue

    rna_gz_path = os.path.join(subdir_path, f"{gcf}_rna_from_genomic.fna.gz")

    # Download with retry
    max_retries = 5
    for attempt in range(1, max_retries + 1):
        print(f"Downloading RNA for {gcf} (attempt {attempt})...")

        try:
            subprocess.run(["wget", "-q", "-O", rna_gz_path, ftp_url], check=True)
            break
        except subprocess.CalledProcessError:
            print(f"WARNING: Download failed (attempt {attempt})")
            if attempt == max_retries:
                print(f"ERROR: Failed to download {gcf}")
                if os.path.exists(rna_gz_path):
                    os.remove(rna_gz_path)
                continue

    if not os.path.exists(rna_gz_path):
        continue

    # Parse RNA sequences
    try:
        with gzip.open(rna_gz_path, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        print(f"ERROR reading {gcf}: {e}")
        continue

    # Filter 16S
    r16S_candidates = [
        rec for rec in records
        if "16S ribosomal RNA" in rec.description
        and len(rec.seq) <= 2000
    ]

    if not r16S_candidates:
        print(f"WARNING: No valid 16S found for {gcf}")
        os.remove(rna_gz_path)
        continue

    longest_16S = max(r16S_candidates, key=lambda r: len(r.seq))

    # Apply taxonomy header
    header_base = genome_to_header.get(gcf, gcf)
    count = header_counter.get(header_base, 0) + 1
    header_counter[header_base] = count
    header = f"{header_base}_{count}" if count > 1 else header_base

    # Adjust DNA -> RNA
    seq_str = str(longest_16S.seq).upper().replace("T", "U")
    seq_str = "".join([b for b in seq_str if b in "ACGU"])
    longest_16S.seq = Seq(seq_str)

    longest_16S.id = header
    longest_16S.name = ""
    longest_16S.description = ""

    # Write output
    with open(final_fasta, "w") as out_handle:
        writer = FastaWriter(out_handle, wrap=None)
        writer.write_record(longest_16S)

    print(f"Saved 16S ({len(longest_16S.seq)} bp) to {final_fasta}")

    # Cleanup zipped files
    os.remove(rna_gz_path)


# Remove unnecessary outputs
for file in os.listdir(output_dir):
    if file.startswith("GCF_") and file.endswith("_16S.fna"):
        gcf_name = file.replace("_16S.fna", "")
        if gcf_name not in current_gcfs:
            os.remove(os.path.join(output_dir, file))
            print(f"Removed outdated file: {file}")

print("\nAll done!")