#!/usr/bin/env python3
"""
Map nif and gdmH proteins onto genome annotations and analyse their genomic
context using GenBank (.gbff) files.

The script:
- Reads a genome-level protein lookup table containing nif and gdmH hits.
- Locates corresponding CDS features in genome GBFF files.
- Records genomic positions, replicons, and sequence types (chromosome/plasmid)
  for nif and gdmH proteins.
- Calculates pairwise distances between all gdmH and nif genes within each genome.
- Summarizes chromosome and plasmid counts per genome.

Outputs:
    - NIF gene locations
    - gdmH gene locations
    - gdmH–nif genomic relationship table
    - Genome chromosome/plasmid summary
"""

import pandas as pd
import csv
from Bio import SeqIO
import argparse
import os

def clean_protein_ids(prot_entry):
    if pd.isna(prot_entry):
        return []
    
    prot_ids = []
    for p in str(prot_entry).split(','):
        p = p.strip()
        if p.startswith("WP_"):
            parts = p.split('_')
            if len(parts) >= 2:
                accession = f"{parts[0]}_{parts[1]}.1"
                prot_ids.append(accession)
            else:
                prot_ids.append(p)
        else:
            prot_ids.append(p)
    return prot_ids

parser = argparse.ArgumentParser(description="Map nif and gdmh proteins and analyse genomic context.")
parser.add_argument("-i", "--input", required=True, help="Input TSV containing genome, nif protein, and gdmH protein assignments.")
parser.add_argument("-o", "--output", required=True, help="Output CSV file for nif gene locations. Additional summary and relationship files are written to the same directory.")
parser.add_argument("-g", "--gbff_dir", default=".", help="Directory containing genome GenBank (*.gbff) files named <Genome>.gbff (default: current directory).")
args = parser.parse_args()

lookup_path = args.input
output_csv = args.output
gbff_dir = args.gbff_dir

# Outputs
summary_csv = os.path.join(os.path.dirname(output_csv), "genome_seq_summary.csv")
gdmh_csv = os.path.join(os.path.dirname(output_csv), "gdmh_locations.csv")
rel_csv = os.path.join(os.path.dirname(output_csv), "gdmh_nif_relationships.csv")

# Load table
lookup = pd.read_csv(lookup_path, sep='\t', engine='python')

# Identify NIF genes
nif_columns = sorted({col.split('_')[0] for col in lookup.columns if col.upper().startswith('NIF')})
print(f"Detected NIF genes: {nif_columns}")

# Open outputs
with open(output_csv, "w", newline="") as f_genes, \
     open(summary_csv, "w", newline="") as f_summary, \
     open(gdmh_csv, "w", newline="") as f_gdmh, \
     open(rel_csv, "w", newline="") as f_rel:

    writer_genes = csv.writer(f_genes)
    writer_summary = csv.writer(f_summary)
    writer_gdmh = csv.writer(f_gdmh)
    writer_rel = csv.writer(f_rel)

    writer_genes.writerow(["Genome", "Replicon", "SeqType", "NifGene", "GeneName", "ProteinID", "Start", "End", "Strand"])
    writer_summary.writerow(["Genome", "ChromosomeCount", "PlasmidCount"])
    writer_gdmh.writerow(["Genome", "Replicon", "SeqType", "ProteinID", "GeneName", "Start", "End", "Strand"])
    writer_rel.writerow([
        "Genome",
        "GDMH_Replicon", "NIF_Replicon",
        "GDMH_Protein", "GDMH_Start", "GDMH_End",
        "NifGene", "NIF_Protein", "NIF_Start", "NIF_End",
        "Distance", "SameReplicon",
        "GDMH_SeqType", "NIF_SeqType"
    ])

    # Iterate genomes
    for genome in lookup['Genome'].unique():
        gbff_file = os.path.join(gbff_dir, f"{genome}.gbff")

        if not os.path.isfile(gbff_file):
            print(f"Warning: GBFF file not found for {genome}, skipping.")
            continue

        df_genome = lookup[lookup['Genome'] == genome]

        # NIF mapping
        protein_map = {}
        for nif_gene in nif_columns:
            col = f"{nif_gene}_protein"
            if col not in df_genome.columns:
                continue
            for entry in df_genome[col].dropna():
                for pid in clean_protein_ids(entry):
                    protein_map[pid] = nif_gene

        # GdmH mapping
        gdmh_proteins = set()
        if "ProteinIDs" in df_genome.columns:
            for entry in df_genome["ProteinIDs"].dropna():
                for prot in str(entry).split(','):
                    prot = prot.strip()
                    if prot:
                        gdmh_proteins.add(prot)

        # Storage
        nif_hits = []
        gdmh_hits = []

        chrom_count = 0
        plasmid_count = 0

        print(f"Parsing {gbff_file}...")

        for record in SeqIO.parse(gbff_file, "genbank"):
            replicon_id = record.id
            seq_type = "plasmid" if "plasmid" in record.description.lower() else "chromosome"

            if seq_type == "plasmid":
                plasmid_count += 1
            else:
                chrom_count += 1

            for feature in record.features:
                if feature.type != "CDS":
                    continue

                gbff_prot_id = feature.qualifiers.get('protein_id', [None])[0]
                if not gbff_prot_id:
                    continue
                gbff_prot_id = gbff_prot_id.strip()

                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                gene_name = feature.qualifiers.get(
                    'gene',
                    feature.qualifiers.get('locus_tag', [gbff_prot_id])
                )[0]

                # NIF
                if gbff_prot_id in protein_map:
                    nif_gene = protein_map[gbff_prot_id]

                    writer_genes.writerow([
                        genome, replicon_id, seq_type,
                        nif_gene, gene_name, gbff_prot_id,
                        start, end, strand
                    ])

                    nif_hits.append({
                        "replicon": replicon_id,
                        "gene": nif_gene,
                        "protein": gbff_prot_id,
                        "start": start,
                        "end": end,
                        "seq_type": seq_type
                    })

                # gdmH
                if gbff_prot_id in gdmh_proteins:
                    writer_gdmh.writerow([
                        genome, replicon_id, seq_type,
                        gbff_prot_id, gene_name,
                        start, end, strand
                    ])

                    gdmh_hits.append({
                        "replicon": replicon_id,
                        "protein": gbff_prot_id,
                        "start": start,
                        "end": end,
                        "seq_type": seq_type
                    })

        # Relationship analysis
        for g in gdmh_hits:
            for n in nif_hits:

                same_replicon = (g["replicon"] == n["replicon"])

                distance = min(
                    abs(g["start"] - n["end"]),
                    abs(g["end"] - n["start"])
                )

                writer_rel.writerow([
                    genome,
                    g["replicon"],
                    n["replicon"],
                    g["protein"], g["start"], g["end"],
                    n["gene"], n["protein"], n["start"], n["end"],
                    distance,
                    "yes" if same_replicon else "no",
                    g["seq_type"],
                    n["seq_type"]
                ])

        writer_summary.writerow([genome, chrom_count, plasmid_count])

print("\nDone!")
print(f"NIF genes -> {output_csv}")
print(f"gdmH genes -> {gdmh_csv}")
print(f"Relationships -> {rel_csv}")
print(f"Summary -> {summary_csv}")