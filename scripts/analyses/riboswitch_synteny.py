#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from collections import defaultdict
from statistics import median
import json
import re

RFAM_IDS = {
    "RF00442",
    "RF01763"
}


def feature_strand(feature):
    return "+" if feature.location.strand == 1 else "-"


def get_product(feature):
    return feature.qualifiers.get("product", ["unknown"])[0]


def get_locus_tag(feature):
    return feature.qualifiers.get("locus_tag", ["unknown"])[0]


def get_riboswitches(record):
    hits = []

    for feat in record.features:
        if feat.type != "regulatory":
            continue

        dbx = feat.qualifiers.get("db_xref", [])
        inf = feat.qualifiers.get("inference", [])

        text = " ".join(dbx + inf)

        if any(rfam_id in text for rfam_id in RFAM_IDS):
            hits.append(feat)

    return hits


def get_cds_features(record):
    cds = [f for f in record.features if f.type == "CDS"]
    cds.sort(key=lambda x: int(x.location.start))
    return cds


def downstream_genes(ribo, cds_features):
    """
    Return CDS features downstream in transcriptional direction.
    """

    ribo_start = int(ribo.location.start)
    ribo_end = int(ribo.location.end)

    if ribo.location.strand == 1:

        return [
            cds
            for cds in cds_features
            if int(cds.location.start) >= ribo_end
        ]

    return [
        cds
        for cds in reversed(cds_features)
        if int(cds.location.end) <= ribo_start
    ]


def intergenic_gap(gene1, gene2):
    """
    Distance between two genes in transcriptional order.
    Overlaps return negative values.
    """

    if gene1.location.strand == 1:
        return (
            int(gene2.location.start)
            - int(gene1.location.end)
        )

    return (
        int(gene1.location.start)
        - int(gene2.location.end)
    )


def infer_operon(ribo, cds_features, max_gap=200):
    """
    Infer operon downstream of riboswitch.
    """

    downstream = downstream_genes(
        ribo,
        cds_features
    )

    if not downstream:
        return []

    operon = [downstream[0]]

    for gene in downstream[1:]:

        prev = operon[-1]

        if gene.location.strand != prev.location.strand:
            break

        gap = intergenic_gap(prev, gene)

        if gap > max_gap:
            break

        operon.append(gene)

    return operon


def leader_distance(ribo, first_gene):
    """
    Distance between riboswitch and first CDS in transcriptional direction.
    """

    if ribo.location.strand == 1:
        return (
            int(first_gene.location.start)
            - int(ribo.location.end)
        )

    return (
        int(ribo.location.start)
        - int(first_gene.location.end)
    )

def feature_coords(feature):
    start = int(feature.location.start) + 1
    end = int(feature.location.end)

    if feature.location.strand == -1:
        return f"complement({start}..{end})"

    return f"{start}..{end}"
    
def get_organism(record):
    return record.annotations.get("organism", "unknown")
    
def is_agmatinase(product):
    p = product.lower()
    return ("agmatinase" in p) or ("agmatine ureohydrolase" in p)

def is_abc(product):
    return "abc" in product.lower()


def classify_operon(genes):
    products = [get_product(g) for g in genes]

    has_agmatinase = any(is_agmatinase(p) for p in products)
    has_abc = any(is_abc(p) for p in products)

    if has_agmatinase and not has_abc:
        return "AGMATINASE_ONLY"

    if has_abc and not has_agmatinase:
        return "ABC_ONLY"

    return None

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("gbff_dir")

    parser.add_argument(
        "--max-gap",
        type=int,
        default=200,
        help="Maximum intergenic gap allowed within an operon"
    )

    args = parser.parse_args()

    gbff_dir = Path(args.gbff_dir)

    patterns = Counter()

    summary_rows = []
    representatives = []
    riboswitch_fasta = []
    mismatch_operons = []
  
    TARGET_GENOMES = {
        "GCF_050408035.1.gbff",
        "GCF_900111445.1.gbff",
        "GCF_001458695.1.gbff",
        "GCF_001751245.1.gbff",
        "GCF_021088345.1.gbff",
        "GCF_019966495.1.gbff",
        "GCF_039724765.1.gbff",
        "GCF_000168975.1.gbff"
    }

    for gbff in sorted(gbff_dir.glob("*.gbff")):

        try:

            records = SeqIO.parse(gbff, "genbank")

            for record in records:

                cds_features = get_cds_features(record)

                for ribo in get_riboswitches(record):

                    genes = infer_operon(
                        ribo,
                        cds_features,
                        max_gap=args.max_gap
                    )

                    if not genes:
                        continue
                        
                    classification = classify_operon(genes)

                    if classification:
                        mismatch_operons.append({
                            "genome": gbff.name,
                            "organism": get_organism(record),
                            "record": record.id,
                            "ribo_coords": feature_coords(ribo),
                            "classification": classification,
                            "pattern": tuple(get_product(g) for g in genes),
                            "leader_distance": leader_distance(ribo, genes[0]),
                            "operon_length": len(genes),
                            "genes": [
                                {
                                    "locus_tag": get_locus_tag(g),
                                    "product": get_product(g),
                                    "start": int(g.location.start),
                                    "end": int(g.location.end),
                                    "strand": g.location.strand
                                }
                                for g in genes
                            ]
                        })

                    first_gene = genes[0]

                    products = tuple(get_product(g) for g in genes)
                    
                    patterns[products] += 1

                    dist = leader_distance(
                        ribo,
                        first_gene
                    )

                    operon_length = len(genes)

                    operon_span = abs(
                        int(genes[0].location.start)
                        - int(genes[-1].location.end)
                    )

                    ribo_seq = str(ribo.extract(record.seq))

                    riboswitch_fasta.append(
                        (f"{gbff.stem}|{record.id}", ribo_seq)
                    )

                    representatives.append(
                        {
                            "genome": gbff.name,
                            "organism": get_organism(record),
                            "record": record.id,
                            "ribo_seq": ribo_seq,
                            "ribo_coords": feature_coords(ribo),
                            "ribo_start": int(ribo.location.start),
                            "ribo_end": int(ribo.location.end),
                            "ribo_strand": ribo.location.strand,
                            "pattern": products,
                            "leader_distance": dist,
                            "operon_length": operon_length,
                            "operon_span": operon_span,
                            "genes": [
                                {
                                    "locus_tag": get_locus_tag(g),
                                    "product": get_product(g),
                                    "start": int(g.location.start),
                                    "end": int(g.location.end),
                                    "strand": g.location.strand,
                                    "sequence": str(g.extract(record.seq)),
                                    "translation": g.qualifiers.get("translation", [""])[0]
                                }
                                for g in genes
                            ]
                        }
                    )

                    summary_rows.append(
                        (
                            gbff.name,
                            feature_strand(ribo),
                            feature_coords(ribo),
                            get_locus_tag(first_gene),
                            get_product(first_gene),
                            dist,
                            operon_length,
                            operon_span,
                            " -> ".join(products)
                        )
                    )

        except Exception as e:
            print(f"ERROR\t{gbff.name}\t{e}")

    pattern_groups = defaultdict(list)
    for hit in representatives:
        pattern_groups[hit["pattern"]].append(hit)

    selected = []
    seen_organisms = set()
    
    dominant_pattern, dominant_count = patterns.most_common(1)[0]
    dominant_hits = pattern_groups.get(dominant_pattern, [])
    
    if dominant_hits:
    
        # sort by numeric GCF ID ascending
        dominant_hits_sorted = sorted(
            dominant_hits,
            key=lambda h: int(h.get("genome", "").split('_')[1].split('.')[0]) if "_" in h.get("genome", "") else 0
        )
    
        # genome ID sorted indices to pick
        pick_indices = [1, 8, 9, 14]
    
        for idx in pick_indices:
            if idx >= len(dominant_hits_sorted):
                continue
    
            h = dominant_hits_sorted[idx]
    
            org = h.get("organism", "unknown")
            if org in seen_organisms:
                continue
    
            h["pattern_count"] = dominant_count
            # selected.append(h)            # DISABLED
            # seen_organisms.add(org)

    # include target genomes
    def score(h):
        pat_hits = pattern_groups.get(h["pattern"], [])
        if pat_hits:
            m = median(x["leader_distance"] for x in pat_hits)
            return abs(h["leader_distance"] - m)
        return float("inf")
    
    for target_genome in TARGET_GENOMES:
    
        target_hits = [
            h for h in representatives
            if h["genome"] == target_genome
        ]
    
        if not target_hits:
            continue
    
        forced_rep = min(target_hits, key=score)
        forced_rep["pattern_count"] = patterns[forced_rep["pattern"]]
    
        if forced_rep not in selected:
            selected.append(forced_rep)

    print("\n=== REPRESENTATIVE LOCI ===\n")

    for rep in selected:

        print(
            f"{rep['pattern_count']}\t"
            f"{rep['genome']}\t"
            f"{rep['record']}\t"
            + " -> ".join(rep["pattern"])
        )
        
    # sort all candidates by numeric GCF ID ascending
    representatives = sorted(
        representatives,
        key=lambda h: int(h.get("genome", "").split('_')[1].split('.')[0]) if "_" in h.get("genome", "") else 0
    )
    with open("all_candidates.json", "w") as out:
        json.dump(representatives, out, indent=2)

    with open("selected_representatives.json", "w") as out:
        json.dump(selected, out, indent=2)
        
    with open("mismatch_operons.json", "w") as out:
        json.dump(mismatch_operons, out, indent=2)

    with open("riboswitches.fa", "w") as out:
        for name, seq in riboswitch_fasta:
            out.write(f">{name}\n{seq}\n")

    with open("representative_riboswitches.fa", "w") as out:
        for rep in selected:
            out.write(f">{rep['genome']}|{rep['record']}\n")
            out.write(rep["ribo_seq"] + "\n")


if __name__ == "__main__":
    main()