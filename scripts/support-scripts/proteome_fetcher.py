#!/usr/bin/env python3

"""
Fetch genome assemblies from NCBI for organisms listed in a taxonomy TSV.

The input TSV is expected to contain one organism per row with taxonomy
columns including label, domain, phylum, class, order, family, genus,
and species, where label is an underscore concatenation of:
domain_pylum_class_order_family_genus_species.

The script resolves organism names to TaxIDs, retrieves and filters
assemblies, downloads genome datasets, generates mapping files, and
runs downstream FASTA processing.
"""

import os
import re
import json
import glob
import shutil
import readline
import subprocess

from pathlib import Path

# Helpers
def run_cmd(cmd, input_text=None, return_stderr=False):
    """
    Run shell command and return stdout.
    Optionally return stderr too.
    """

    result = subprocess.run(
        cmd,
        input=input_text,
        text=True,
        capture_output=True
    )

    if result.returncode != 0:

        print("Command failed:", " ".join(cmd))
        print(result.stderr.strip())

        if return_stderr:
            return "", result.stderr.strip()

        return ""

    if return_stderr:
        return result.stdout.strip(), result.stderr.strip()

    return result.stdout.strip()


def load_cache(path):
    cache = {}

    if os.path.exists(path):

        with open(path) as f:

            for line in f:

                parts = line.strip().split("\t")

                if len(parts) >= 2:

                    key = parts[0]
                    value = parts[1]

                    cache.setdefault(key, []).append(value)

    return cache


def append_cache(path, key, value):

    with open(path, "a") as f:
        f.write(f"{key}\t{value}\n")


def complete_path(text, state):

    text = os.path.expanduser(text)
    text = os.path.expandvars(text)

    matches = glob.glob(text + "*")

    matches = [
        m + "/" if os.path.isdir(m) else m
        for m in matches
    ]

    try:
        return matches[state]

    except IndexError:
        return None


# Taxonomy utilities
SYNONYMS = {
    "Sinorhizobium": "Ensifer",
}


def build_candidate_names(genus, species):

    genus = genus.strip()
    species = species.strip()

    # Prevent duplicated genus names
    if species.lower().startswith(genus.lower()):
        organism = species

    else:
        organism = f"{genus} {species}"

    # Replace underscores with spaces
    organism = organism.replace("_", " ")

    # Remove unwanted suffixes
    organism = re.sub(
        r"\bunfinished sequence\b",
        "",
        organism,
        flags=re.I
    )

    organism = re.sub(
        r"\bcomplete genome\b",
        "",
        organism,
        flags=re.I
    )

    # Normalise whitespace
    organism = " ".join(organism.split())

    # Apply synonym normalisation
    parts = organism.split()

    if parts and parts[0] in SYNONYMS:
        parts[0] = SYNONYMS[parts[0]]

    organism = " ".join(parts)

    # Candidate search names
    candidate_names = []

    # Full cleaned name
    candidate_names.append(organism)

    tokens = organism.split()

    # First two words only
    if len(tokens) >= 2:
        candidate_names.append(
            f"{tokens[0]} {tokens[1]}"
        )

    # Remove biovar notation
    bv_removed = re.sub(
        r"\bbv\.\s+\S+",
        "",
        organism
    )

    bv_removed = " ".join(bv_removed.split())

    if bv_removed not in candidate_names:
        candidate_names.append(bv_removed)

    # Remove duplicates
    candidate_names = list(
        dict.fromkeys(candidate_names)
    )

    return candidate_names


def fetch_taxid(name):

    cmd = [
        "datasets",
        "summary",
        "taxonomy",
        "taxon",
        name
    ]

    output, stderr = run_cmd(
        cmd,
        return_stderr=True
    )

    # Successful exact match
    if output:

        try:

            data = json.loads(output)

            reports = data.get("reports", [])

            if reports:

                taxonomy = reports[0].get(
                    "taxonomy",
                    {}
                )

                taxid = taxonomy.get("tax_id")

                if taxid:
                    return str(taxid)

        except Exception:
            pass

    # Parse suggested taxids
    suggestions = []

    capture = False

    for line in stderr.splitlines():

        line = line.strip()

        if "Try using one of the suggested taxids" in line:
            capture = True
            continue

        if capture and line:

            m = re.search(
                r"^(.*?)\s+\(.*taxid:\s*([0-9]+)",
                line
            )

            if m:

                suggestion_name = m.group(1).strip()
                suggestion_taxid = m.group(2).strip()

                suggestions.append(
                    (
                        suggestion_name,
                        suggestion_taxid
                    )
                )

    # Interactive user selection
    if suggestions:

        print("\nSuggested taxonomy matches:\n")

        for i, (s_name, s_taxid) in enumerate(
            suggestions,
            start=1
        ):

            print(
                f"{i}. "
                f"{s_name} "
                f"(TaxID: {s_taxid})"
            )

        print("0. Skip organism")

        while True:

            choice = input(
                "\nSelect taxonomy match: "
            ).strip()

            if choice == "0":
                return None

            if (
                choice.isdigit()
                and 1 <= int(choice) <= len(suggestions)
            ):

                selected = suggestions[
                    int(choice) - 1
                ]

                selected_name = selected[0]
                selected_taxid = selected[1]

                print(
                    f"Selected: "
                    f"{selected_name} "
                    f"(TaxID: {selected_taxid})"
                )

                return selected_taxid

            print("Invalid selection.")

    return None


def fetch_assemblies_by_taxid(taxid):

    cmd = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        str(taxid),
        "--assembly-source",
        "all",
        "--report",
        "genome"
    ]

    output = run_cmd(cmd)

    if not output:
        return []

    try:
        data = json.loads(output)

    except Exception:
        return []

    assemblies = []

    for r in data.get("reports", []):

        acc = (
            r.get("accession")
            or r.get("current_accession")
        )

        if acc:
            assemblies.append(acc)

    return assemblies


def fetch_assemblies_by_text(query):

    cmd = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        query,
        "--assembly-source",
        "all",
        "--report",
        "genome"
    ]

    output = run_cmd(cmd)

    if not output:
        return []

    try:
        data = json.loads(output)

    except Exception:
        return []

    assemblies = []

    for r in data.get("reports", []):

        acc = r.get("accession")

        if acc:
            assemblies.append(acc)

    return assemblies


# Setup
readline.set_completer(complete_path)

readline.parse_and_bind(
    "tab: complete"
)

input_file = input(
    "\nProvide path to ribo_taxonomy_updated.tsv: "
).strip()

if not os.path.isfile(input_file):

    print(f"File not found: {input_file}")

    exit(1)

workdir = Path("genome_fetcher")

genomedir = workdir / "genomes"

mapping_file = (
    workdir / "taxonomy_to_assembly.txt"
)

filtered_mapping_file = (
    workdir / "filtered_taxonomy_to_assembly.txt"
)

assembly_metadata = (
    workdir / "assembly_metadata.json"
)

assemblies_file = (
    workdir / "assemblies.txt"
)

organism_tax_cache_file = (
    workdir / "organism_to_taxid.tsv"
)

tax_cache_file = (
    workdir / "taxid_to_assemblies.tsv"
)

workdir.mkdir(exist_ok=True)

genomedir.mkdir(exist_ok=True)

organism_tax_cache = load_cache(
    organism_tax_cache_file
)

tax_cache = load_cache(
    tax_cache_file
)

# Reuse metadata if previously downloaded
USE_EXISTING_METADATA = False

if assembly_metadata.exists():

    choice = input(
        "\nMetadata exists. Use it? [1=yes, 2=no]: "
    ).strip() or "1"

    USE_EXISTING_METADATA = (
        choice == "1"
    )

# Map taxonomy -> assemblies
if not USE_EXISTING_METADATA:

    print("\nExpanding taxonomy to assemblies")

    assemblies_set = set()

    with open(mapping_file, "w") as mapping_out:

        with open(input_file) as f:

            header = next(f)

            for line in f:

                parts = line.rstrip(
                    "\n"
                ).split("\t")

                if len(parts) < 8:
                    continue

                (
                    label,
                    domain,
                    phylum,
                    class_,
                    order,
                    family,
                    genus,
                    species
                ) = parts[:8]

                print(
                    f"\n{'='*60}"
                )

                print(
                    f"Processing row: {label}"
                )

                candidate_names = (
                    build_candidate_names(
                        genus,
                        species
                    )
                )

                print(
                    "Candidate names:"
                )

                for c in candidate_names:
                    print(f"  {c}")

                taxid = None
                organism = None

                # Taxonomy lookup
                for candidate in candidate_names:

                    print(
                        f"\nTrying taxonomy lookup: "
                        f"{candidate}"
                    )

                    if candidate in organism_tax_cache:

                        taxid = (
                            organism_tax_cache[
                                candidate
                            ][0]
                        )

                        print(
                            f"Using cached TaxID: "
                            f"{taxid}"
                        )

                    else:

                        taxid = fetch_taxid(
                            candidate
                        )

                        if taxid:

                            append_cache(
                                organism_tax_cache_file,
                                candidate,
                                taxid
                            )

                            organism_tax_cache[
                                candidate
                            ] = [taxid]

                    if taxid:

                        organism = candidate

                        print(
                            f"Resolved TaxID: "
                            f"{taxid}"
                        )

                        break

                assemblies = []


                # TaxID -> assemblies
                if taxid:

                    if taxid in tax_cache:

                        assemblies = (
                            tax_cache[taxid]
                        )

                        print(
                            f"Using cached assemblies "
                            f"({len(assemblies)})"
                        )

                    else:

                        print(
                            f"Fetching assemblies "
                            f"for TaxID {taxid}"
                        )

                        assemblies = (
                            fetch_assemblies_by_taxid(
                                taxid
                            )
                        )

                        for acc in assemblies:

                            append_cache(
                                tax_cache_file,
                                taxid,
                                acc
                            )

                        tax_cache[
                            taxid
                        ] = assemblies

                # Fallback text search
                if not assemblies:

                    print(
                        "\nNo assemblies found "
                        "via TaxID."
                    )

                    print(
                        "Trying direct assembly search..."
                    )

                    for candidate in candidate_names:

                        assemblies = (
                            fetch_assemblies_by_text(
                                candidate
                            )
                        )

                        if assemblies:

                            organism = candidate

                            print(
                                f"Recovered "
                                f"{len(assemblies)} "
                                f"assemblies"
                            )

                            break

                # Failed completely
                if not assemblies:

                    print(
                        "No assemblies recovered."
                    )

                    continue

                # Save mappings
                for acc in assemblies:

                    assemblies_set.add(acc)

                    mapping_out.write(
                        f"{label}\t"
                        f"{organism or 'NA'}\t"
                        f"{taxid or 'NA'}\t"
                        f"{acc}\n"
                    )

    # Write assemblies list
    with open(assemblies_file, "w") as f:

        for acc in sorted(assemblies_set):
            f.write(acc + "\n")

    print(
        f"\nTotal assemblies: "
        f"{len(assemblies_set)}"
    )

    if not assemblies_set:

        print("ERROR: No assemblies found")

        exit(1)

    # If not already present, download metadata

    print(
        "\nDownloading assembly metadata..."
    )

    with open(assembly_metadata, "w") as out:

        subprocess.run([
            "datasets",
            "summary",
            "genome",
            "accession",
            "--inputfile",
            str(assemblies_file),
            "--report",
            "genome"
        ], stdout=out)

else:

    print(
        "Using cached metadata"
    )

# Filtering
print("\nFilter options:")

print("1. No filter")
print("2. Isolate")
print("3. Isolate + Complete")
print("4. Isolate + Complete + Reference")

choice = int(
    input("Choice [4]: ") or "4"
)

FILTER_MAP = {
    1: "no_filter",
    2: "isolate",
    3: "isolate_complete",
    4: "isolate_complete_reference"
}

filter_type = FILTER_MAP[choice]

filtered_file = (
    workdir /
    f"{filter_type}_filtered_assemblies.txt"
)

if choice == 1:

    filtered_file.write_text(
        assemblies_file.read_text()
    )

else:

    data = json.load(
        open(assembly_metadata)
    )

    kept = []

    for r in data.get("reports", []):

        acc = r.get("accession")

        info = r.get(
            "assembly_info",
            {}
        )

        attrs = (
            info.get(
                "biosample",
                {}
            ).get(
                "attributes",
                []
            )
        )

        def has_attr(name):

            return any(
                a.get("name") == name
                for a in attrs
            )

        def attr_contains(name, pattern):

            return any(
                a.get("name") == name and
                pattern.lower() in
                a.get(
                    "value",
                    ""
                ).lower()
                for a in attrs
            )

        # Isolate material filter
        if choice >= 2:

            org_name = (
                r.get(
                    "organism",
                    {}
                ).get(
                    "organism_name",
                    ""
                ).lower()
            )

            if "metagenome" in org_name:
                continue

            if attr_contains(
                "env_package",
                "metagenome"
            ):
                continue

            if attr_contains(
                "env_package",
                "environmental"
            ):
                continue

            if any(
                a.get(
                    "name",
                    ""
                ).startswith("env_")
                for a in attrs
            ):
                continue

            if not (
                has_attr("strain")
                or has_attr(
                    "type-material"
                )
            ):
                continue

        # Complete genomes filter
        if choice >= 3:

            if (
                info.get(
                    "assembly_level"
                )
                != "Complete Genome"
            ):
                continue

        # Reference genomes filter
        if choice == 4:

            if (
                info.get(
                    "refseq_category"
                )
                != "reference genome"
            ):
                continue

        kept.append(acc)

    with open(filtered_file, "w") as f:

        for acc in kept:
            f.write(acc + "\n")

    print(
        f"Filtered assemblies: "
        f"{len(kept)}"
    )

# Filter mappings
valid = set(
    open(filtered_file)
    .read()
    .splitlines()
)

with open(mapping_file) as inp, \
     open(filtered_mapping_file, "w") as out:

    for line in inp:

        label, organism, taxid, acc = \
            line.strip().split("\t")

        if acc in valid:
            out.write(line)

# Download genomes
zip_file = (
    genomedir /
    f"{filter_type}_genomes.zip"
)

output_dir = (
    workdir /
    f"{filter_type}_genomes"
)

redownload = True

if (
    output_dir.exists()
    and any(output_dir.iterdir())
):

    choice = input(
        f"\nFolder '{output_dir}' exists. "
        f"Redownload? [1=no, 2=yes]: "
    ).strip() or "1"

    redownload = (
        choice == "2"
    )

    if redownload:

        print(
            f"Removing {output_dir}..."
        )

        shutil.rmtree(output_dir)

        output_dir.mkdir(
            exist_ok=True
        )

    else:

        print(
            "\nSkipping download"
        )

if redownload:

    print(
        f"\nDownloading genomes "
        f"to {zip_file}..."
    )

    subprocess.run([
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        str(filtered_file),
        "--include",
        "protein",
        "--filename",
        str(zip_file)
    ])

    tmp_dir = (
        genomedir / "tmp"
    )

    tmp_dir.mkdir(
        exist_ok=True
    )

    subprocess.run([
        "unzip",
        "-o",
        str(zip_file),
        "-d",
        str(tmp_dir)
    ])

# Reorganise downloaded genomes
if redownload:

    input_base = (
        genomedir /
        "tmp/ncbi_dataset/data"
    )

    output_dir.mkdir(
        exist_ok=True
    )

    for d in input_base.glob("GC*"):

        d.rename(
            output_dir / d.name
        )

# Cleanup of unnecessary files
if genomedir.exists():

    shutil.rmtree(genomedir)

# FASTA processing
if redownload:

    subprocess.run([
        "python3",
        "scripts/fasta_processor.py",
        str(output_dir)
    ])

print("\nDone!")

print(f"Genomes: {output_dir}")

print(f"Mapping: {mapping_file}")

print(
    f"Filtered mapping: "
    f"{filtered_mapping_file}"
)

print(
    f"Assembly metadata: "
    f"{assembly_metadata}"
)