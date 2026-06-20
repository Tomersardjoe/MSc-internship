#!/usr/bin/env python3

"""
Build a gene presence/absence table from BLAST search results.

The script:
- Loads and filters BLAST hits using E-value, percent identity, and coverage thresholds.
- Extracts gene names from subject IDs for either the 'gdmh' or 'nif' gene set.
- Maps filtered hits to genome proteomes and taxonomy metadata.
- Generates a genome-by-gene presence/absence table, including matched protein IDs and sequences.
- Produces diagnostic BLAST plots and optionally copies genomes containing gene hits.
- Outputs a TSV summary sorted by total gene presence.

Inputs:
    - BLAST hit files (*_blast_hits.tsv)
    - Genome proteomes (*.faa)
    - Taxonomy mapping TSV

Output:
    <prefix>_<gene_set>_presence_absence.tsv
"""

import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pathlib import Path
import shutil

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build gene presence/absence table from BLAST hits"
    )
    parser.add_argument(
        "-b",
        "--blast",
        required=True,
        help="Directory containing blast results or subdirectories that do (*_blast_hits.tsv)"
    )
    parser.add_argument(
        "-p",
        "--proteomes",
        required=True,
        help="Directory containing genome proteome FASTA files or subdirectories that do (*.faa)"
    )
    parser.add_argument(
        "-t",
        "--taxonomy",
        required=True,
        help="TSV file mapping genome --> taxid, class, order, family, genus, species"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="blast",
        help="Output prefix (default: blast)"
    )
    parser.add_argument(
        "--gene_set",
        choices=["nif", "gdmh"],
        default="gdmh",
        help="Gene set to analyze (default: gdmh)"
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-20,
        help="E-value threshold (default: 1e-20)"
    )
    parser.add_argument(
        "--pident",
        type=float,
        default=35,
        help="Percent identity threshold (default: 35)"
    )
    parser.add_argument(
        "--cov",
        type=float,
        default=0.7,
        help="Coverage threshold (default: 0.7)"
    )
    return parser.parse_args()


# load BLAST files
def load_blast_files(blast_dir):
    files = list(Path(blast_dir).rglob("*_blast_hits.tsv"))
    if len(files) == 0:
        raise RuntimeError("No *_blast_hits.tsv files found")
    dfs = []
    for f in files:
        if os.path.getsize(f) == 0:
            print(f"Skipping empty BLAST file: {f}")
            continue
        genome = os.path.basename(f).replace("_blast_hits.tsv","")
        df = pd.read_csv(f, sep="\t", header=None)
        df.columns = [
            "qseqid","sseqid","pident","length","qstart","qend",
            "qlen","sstart","send","slen","evalue","bitscore"
        ]
        df["genome"] = genome
        df["qcov"] = (df["qend"] - df["qstart"] + 1) / df["qlen"]
        df["qcov"] = df["qcov"].clip(upper=1.0)
        df["scov"] = (abs(df["send"] - df["sstart"]) + 1) / df["slen"]
        df["scov"] = df["scov"].clip(upper=1.0)
        dfs.append(df)
    if len(dfs) == 0:
        raise RuntimeError("No BLAST files with hits found")
    return pd.concat(dfs, ignore_index=True)


# extract gene name
def get_extract_gene_func(gene_set):
    """
    Return the appropriate gene-extraction function based on gene set.
    """
    if gene_set == "gdmh":
        def extract_gene_gdmh(sseqid):
            parts = sseqid.split("_")
            if len(parts) >= 2:
                return parts[1]
            else:
                return ""
        return extract_gene_gdmh
    elif gene_set == "nif":
        def extract_gene_nif(sseqid):
            try:
                return sseqid.split("|")[2].split("_")[0]
            except IndexError:
                return ""
        return extract_gene_nif
    else:
        raise ValueError(f"Unknown gene set: {gene_set}")

# filter hits
def filter_hits(df, evalue, pident, cov, extract_gene_func):
    df = df[
        (df["evalue"] <= evalue) &
        (df["pident"] >= pident) &
        (df["qcov"] >= cov) &
        (df["scov"] >= cov)
    ]
    df = df.sort_values("bitscore", ascending=False)
    df = df.drop_duplicates(["genome","qseqid"])
    df["gene"] = df["sseqid"].apply(extract_gene_func)
    return df


# load proteomes

def load_proteome(fasta_path):
    if not os.path.exists(fasta_path):
        return {}
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


# build presence/absence table

def build_pa_table(df_hits, proteome_dir, taxonomy):
    """
    Build presence/absence table including all genomes that have proteomes.
    """

    genes = sorted(df_hits["gene"].unique())
    rows = []

    # list all genome IDs based on subdirectories in proteome_dir
    all_genomes = [
        d for d in os.listdir(proteome_dir)
        if os.path.isdir(os.path.join(proteome_dir, d))
    ]

    for genome in all_genomes:

        row = {"Genome": genome}

        # add taxonomy info if available
        if genome in taxonomy.index:
            row.update(taxonomy.loc[genome].to_dict())
        else:
            # fill empty if no taxonomy
            row.update({c: "" for c in taxonomy.columns})

        # load proteome sequences
        proteome_path = os.path.join(proteome_dir, genome, f"{genome}.faa")
        proteome = load_proteome(proteome_path)

        # zubset hits for this genome
        if genome in df_hits["genome"].values:
            dfg = df_hits[df_hits["genome"] == genome]
        else:
            dfg = pd.DataFrame(columns=df_hits.columns)

        for gene in genes:

            sub = dfg[dfg["gene"] == gene]

            if len(sub) == 0:
                row[f"{gene}_presence"] = 0
                row[f"{gene}_protein"] = ""
                row[f"{gene}_seq"] = ""
            else:
                hit = sub.iloc[0]
                pid = hit["qseqid"]
                seq = proteome.get(pid, "")

                row[f"{gene}_presence"] = 1
                row[f"{gene}_protein"] = pid
                row[f"{gene}_seq"] = seq

        rows.append(row)

    return pd.DataFrame(rows)

    # sort genomes by total NIF presence
    nif_cols = [c for c in df_out.columns if c.endswith("_presence")]
    df_out["total_nif"] = df_out[nif_cols].sum(axis=1)
    df_out = df_out.sort_values("total_nif", ascending=False)
    df_out = df_out.drop(columns="total_nif")

    # rename 'genome' column to 'Genome'
    df_out = df_out.rename(columns={"genome": "Genome"})
    return df_out



# diagnostic plots
def plot_diagnostics(df, output_prefix, cov_thresh, pid_thresh):
    fig, axes = plt.subplots(1, 2, figsize=(12,5), constrained_layout=True)
    df["log_evalue"] = -np.log10(df["evalue"].replace(0, 1e-300))
    sc1 = axes[0].scatter(df["qcov"], df["pident"], c=df["log_evalue"], cmap="viridis", s=5, alpha=0.5)
    axes[0].axvline(cov_thresh, color="red", linestyle="--")
    axes[0].axhline(pid_thresh, color="red", linestyle="--")
    axes[0].set_xlabel("Query coverage")
    axes[0].set_ylabel("Percent identity")
    axes[0].set_title("Identity vs Coverage")
    sc2 = axes[1].scatter(df["qcov"], df["bitscore"], c=df["log_evalue"], cmap="viridis", s=5, alpha=0.5)
    axes[1].axvline(cov_thresh, color="red", linestyle="--")
    axes[1].set_xlabel("Query coverage")
    axes[1].set_ylabel("Bitscore")
    axes[1].set_title("Bitscore vs Coverage")
    cbar = fig.colorbar(sc2, ax=axes, shrink=0.8)
    cbar.set_label("-log10(E-value)")
    plt.savefig(f"{output_prefix}_diagnostics.png", dpi=300)

def copy_hit_genomes(table, proteome_dir, gene_set, out_dir):
    """
    Copy genomes that have at least one gene hit to a new directory.
    """
    os.makedirs(out_dir, exist_ok=True)
    
    # identify genomes with any presence
    presence_cols = [c for c in table.columns if c.endswith("_presence")]
    genomes_with_hits = table[table[presence_cols].sum(axis=1) > 0]["Genome"]

    for genome in genomes_with_hits:
        src = os.path.join(proteome_dir, genome)
        dst = os.path.join(out_dir, genome)
        if os.path.exists(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
            print(f"Copied genome {genome} to {out_dir}")
        else:
            print(f"Warning: genome directory {src} not found, skipping")

def main():
    args = parse_args()

    print(f"Using gene set: {args.gene_set}")
    print("Loading BLAST files...")
    df_all = load_blast_files(args.blast)

    # choose gene extraction based on selected gene set
    extract_gene = get_extract_gene_func(args.gene_set)
    
    print("Filtering hits...")
    df_filtered = filter_hits(df_all, args.evalue, args.pident, args.cov, extract_gene)

    print("Loading taxonomy mapping...")
    taxonomy_df = pd.read_csv(args.taxonomy, sep="\t").set_index("Assembly")

    print("Building presence/absence table...")
    table = build_pa_table(df_filtered, args.proteomes, taxonomy_df)
    
    if table.empty:
        print("Presence/absence table is empty. Exiting.")
        return

    # copy genomes with hits if gene_set is specified
    if args.gene_set:
        hit_genomes_dir = f"{args.gene_set}_hits_genomes"
        print(f"Copying genomes with {args.gene_set} hits to {hit_genomes_dir} ...")
        copy_hit_genomes(table, args.proteomes, args.gene_set, hit_genomes_dir)
    
        print("Generating diagnostic plots...")
        plot_diagnostics(df_all, args.output, args.cov, args.pident)

    # sort by prevalence
    presence_cols = [c for c in table.columns if c.endswith("_presence")]
    if presence_cols:
        table["total"] = table[presence_cols].sum(axis=1)
        table = table.sort_values("total", ascending=False)
        table = table.drop(columns="total")

    out = f"{args.output}_{args.gene_set}_presence_absence.tsv"
    table.to_csv(out, sep="\t", index=False)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
