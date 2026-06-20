#!/usr/bin/env python3

import argparse
import csv
import re
from Bio import SeqIO


RANKS = [
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]

WEIGHTS = {
    "Domain": 1,
    "Phylum": 2,
    "Class": 4,
    "Order": 8,
    "Family": 16,
    "Genus": 32,
    "Species": 64,
}

LOW_CONFIDENCE_THRESHOLD = 175


def tokenize(text):
    """Convert text into normalized token set."""

    if not text:
        return set()

    text = str(text).lower()

    # normalize phylum synonym
    text = text.replace(
        "proteobacteria",
        "pseudomonadota"
    )

    tokens = re.split(r"[_\s]+", text)

    return {t for t in tokens if t}


def normalize_header(header):
    """Normalize alignment header."""

    header = header.lstrip(">")

    header = re.sub(
        r"_(unfinished_sequence|complete_genome|draft_genome)$",
        "",
        header,
        flags=re.IGNORECASE,
    )

    header = header.replace(
        "Proteobacteria",
        "Pseudomonadota"
    )

    return header


def load_taxonomy(taxonomy_file):
    """Load taxonomy TSV."""

    taxonomy = []

    with open(taxonomy_file) as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            taxonomy.append(row)

    return taxonomy


def score_row(query_tokens, row):
    """Score one taxonomy row against a query."""

    score = 0
    details = []

    for rank in RANKS:

        rank_tokens = tokenize(row.get(rank, ""))

        matches = query_tokens & rank_tokens

        if not matches:
            continue

        if rank == "Species":
            rank_score = len(matches) * WEIGHTS[rank]
        else:
            rank_score = WEIGHTS[rank]

        score += rank_score

        details.append(
            f"{rank}:{','.join(sorted(matches))}"
        )

    return score, details


def find_matches(header, taxonomy):
    """Return all matches sorted by score."""

    query = normalize_header(header)

    query_tokens = tokenize(query)

    hits = []

    for row in taxonomy:

        score, details = score_row(
            query_tokens,
            row
        )

        hits.append(
            (score, row, details)
        )

    hits.sort(
        key=lambda x: x[0],
        reverse=True
    )

    return hits


def format_taxonomy_header(row):
    """Construct taxonomy header from TSV row."""

    fields = [
        row.get("Domain", ""),
        row.get("Phylum", ""),
        row.get("Class", ""),
        row.get("Order", ""),
        row.get("Family", ""),
        row.get("Genus", ""),
        row.get("Species", ""),
    ]

    fields = [
        str(x).strip().replace(" ", "_")
        for x in fields
        if str(x).strip()
    ]

    return "_".join(fields)


def main():

    parser = argparse.ArgumentParser(
        description="Match alignment headers to taxonomy table."
    )

    parser.add_argument(
        "--taxonomy",
        required=True,
        help="Taxonomy TSV file"
    )

    parser.add_argument(
        "--alignment",
        required=True,
        help="Alignment FASTA file"
    )
    
    parser.add_argument(
    "--renamed-alignment",
    default="renamed_alignment.faa",
    help="Output FASTA with renamed headers"
    )
    
    parser.add_argument(
        "--unmatched",
        default="unmatched.tsv",
        help="TSV containing unmatched queries"
    )
    
    parser.add_argument(
        "--decisions",
        default="decisions.tsv",
        help="TSV recording all manual decisions"
    )

    args = parser.parse_args()

    taxonomy = load_taxonomy(args.taxonomy)

    renamed_records = []
    unmatched_rows = []
    decision_rows = []

    for record in SeqIO.parse(
        args.alignment,
        "fasta"
    ):

        hits = find_matches(
            record.id,
            taxonomy
        )

        print("\n" + "=" * 80)
        print("QUERY:")
        print(record.id)

        print("\nTop matches:\n")

        for i in range(3):

            score, row, details = hits[i]

            print(f"{i + 1}) score={score}")

            print(
                "   "
                + format_taxonomy_header(row)
            )

            print(
                "   "
                + "; ".join(details)
            )

            print()

        while True:

            choice = input(
                "Select match "
                "(1-3) or 0 for no match: "
            ).strip()

            if choice in {"0", "1", "2", "3"}:
                break

            print(
                "Invalid selection."
            )

        if choice == "0":

            unmatched_rows.append([
                record.id,
                hits[0][0],
                format_taxonomy_header(hits[0][1]),
                format_taxonomy_header(hits[1][1]),
                format_taxonomy_header(hits[2][1]),
            ])

            decision_rows.append([
                record.id,
                0,
                ""
            ])

            continue

        selected = hits[int(choice) - 1]

        new_header = format_taxonomy_header(
            selected[1]
        )

        renamed_records.append(
            (
                new_header,
                str(record.seq)
            )
        )

        decision_rows.append([
            record.id,
            choice,
            new_header
        ])

    with open(
        args.renamed_alignment,
        "w"
    ) as out_fasta:

        for header, seq in renamed_records:

            out_fasta.write(
                f">{header}\n"
            )

            out_fasta.write(
                f"{seq}\n"
            )

    with open(
        args.unmatched,
        "w",
        newline=""
    ) as out_tsv:

        writer = csv.writer(
            out_tsv,
            delimiter="\t"
        )

        writer.writerow([
            "QueryHeader",
            "BestCandidate",
            "SecondCandidate",
            "ThirdCandidate"
        ])

        writer.writerows(
            unmatched_rows
        )

    with open(
        args.decisions,
        "w",
        newline=""
    ) as out_tsv:

        writer = csv.writer(
            out_tsv,
            delimiter="\t"
        )

        writer.writerow([
            "OriginalHeader",
            "Choice",
            "SelectedHeader"
        ])

        writer.writerows(
            decision_rows
        )

    print(
        f"\nWrote renamed alignment to "
        f"{args.renamed_alignment}"
    )

    print(
        f"Wrote unmatched queries to "
        f"{args.unmatched}"
    )

    print(
        f"Wrote decision log to "
        f"{args.decisions}"
    )

if __name__ == "__main__":
    main()