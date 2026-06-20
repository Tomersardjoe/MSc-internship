#!/usr/bin/env python3
import os
import subprocess
import argparse
import shutil
import sys
import tempfile

# Parse arguments
parser = argparse.ArgumentParser(
    description="Download GBFF files for GCF_* accessions from NCBI."
)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "-d", "--parent_dir",
    help="Directory containing subdirectories starting with GCF_*"
)
group.add_argument(
    "-t", "--txt_file",
    help="Text file with one GCF_* accession per line (alternative to supplying -d)"
)

parser.add_argument(
    "-o", "--output_dir",
    required=True,
    help="Directory to save the final GBFF files"
)

args = parser.parse_args()

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

using_parent_dir = args.parent_dir is not None

# Build list of GCFs
subdirs = []

if using_parent_dir:
    parent_dir = args.parent_dir

    if not os.path.isdir(parent_dir):
        print(f"ERROR: {parent_dir} is not a valid directory.")
        sys.exit(1)

    subdirs = [
        d for d in os.listdir(parent_dir)
        if d.startswith("GCF_") and os.path.isdir(os.path.join(parent_dir, d))
    ]

else:
    txt_file = args.txt_file

    if not os.path.isfile(txt_file):
        print(f"ERROR: {txt_file} does not exist.")
        sys.exit(1)

    with open(txt_file) as f:
        subdirs = [
            line.strip()
            for line in f
            if line.strip().startswith("GCF_")
        ]

if not subdirs:
    print("No valid GCF_* entries found.")
    sys.exit(0)

current_gcfs = set(subdirs)

# Process each accession
for gcf in subdirs:
    print(f"\nProcessing {gcf}...")

    final_gbff = os.path.join(output_dir, f"{gcf}.gbff")

    if os.path.exists(final_gbff):
        print(f"Skipping {gcf}: already exists")
        continue

    work_dir = tempfile.mkdtemp(prefix=f"{gcf}_")
    temp_zip = os.path.join(work_dir, f"{gcf}_ncbi.zip")

    print(f"Downloading genome {gcf} via NCBI datasets...")
    cmd_download = [
        "datasets",
        "download", "genome", "accession", gcf,
        "--include", "gbff",
        "--filename", temp_zip,
        "--no-progressbar"
    ]

    try:
        subprocess.run(cmd_download, check=True)
    except subprocess.CalledProcessError:
        print(f"ERROR: Failed to download {gcf}. Skipping.")
        shutil.rmtree(work_dir, ignore_errors=True)
        continue

    print(f"Unzipping {temp_zip} ...")
    subprocess.run(["unzip", "-o", temp_zip, "-d", work_dir], check=True)

    gbff_files = []
    for root, _, files in os.walk(work_dir):
        for f in files:
            if f.endswith(".gbff"):
                gbff_files.append(os.path.join(root, f))

    if not gbff_files:
        print(f"ERROR: No GBFF found for {gcf}")
        shutil.rmtree(work_dir, ignore_errors=True)
        continue

    gbff_file = gbff_files[0]

    shutil.move(gbff_file, final_gbff)
    print(f"Saved GBFF as {final_gbff}")

    shutil.rmtree(work_dir, ignore_errors=True)

# Cleanup unnecessary outputs
if using_parent_dir:
    for file in os.listdir(output_dir):
        if file.startswith("GCF_") and file.endswith(".gbff"):
            gcf_name = file.replace(".gbff", "")
            if gcf_name not in current_gcfs:
                file_path = os.path.join(output_dir, file)
                print(f"Removing outdated GBFF file: {file_path}")
                os.remove(file_path)

print("\nAll done!")