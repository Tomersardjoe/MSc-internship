import os
import re
import argparse

# Helper function
def sanitize(text):
    return text.replace("/", "_")

# Argument parsing
parser = argparse.ArgumentParser(description="Process .faa files and rename headers")
parser.add_argument(
    "target_dir",
    help="Directory containing .faa files"
)

args = parser.parse_args()
target_dir = args.target_dir

for root, dirs, files in os.walk(target_dir):
    for filename in files:
        if filename.endswith(".faa"):
            file_path = os.path.join(root, filename)
            temp_path = file_path + ".tmp"

            print(f"Processing {file_path}")

            with open(file_path, "r") as infile, open(temp_path, "w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        header = line.strip()

                        # UniProt accession matcher
                        if header.startswith(">sp|") or header.startswith(">tr|"):
                            match = re.match(r">(?:sp|tr)\|([^|]+)\|([^\s]+)", header)
                            gene_match = re.search(r"GN=([^ ]+)", header)

                            if match:
                                accession = sanitize(match.group(1))
                                entry = sanitize(match.group(2))
                                organism_code = sanitize(entry.split("_")[-1])

                                if gene_match:
                                    gene = sanitize(gene_match.group(1).capitalize())
                                else:
                                    gene = "NA"

                                new_header = f">{accession}_{gene}_{organism_code}\n"
                                outfile.write(new_header)
                            else:
                                outfile.write(line)

                        # RefSeq accession matcher
                        elif header.startswith(">WP_") or header.startswith(">NP_") or header.startswith(">XP_"):
                            parts = header[1:].split(None, 1)
                            accession = sanitize(parts[0].replace(".", "_"))
                            description = parts[1] if len(parts) > 1 else ""

                            org_match = re.search(r"\[([^\]]+)\]", description)
                            organism = org_match.group(1) if org_match else "NA"
                            organism = sanitize(organism.replace(" ", "_"))

                            new_header = f">{accession}_{organism}\n"
                            outfile.write(new_header)

                        # Other headers
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)

            # Replace original file with modified headers
            os.replace(temp_path, file_path)

            # Rename file to match parent directory
            parent_dir = os.path.basename(root)
            new_file_path = os.path.join(root, f"{parent_dir}.faa")
            
            if file_path != new_file_path:
                os.rename(file_path, new_file_path)
                print(f"Renamed {file_path} -> {new_file_path}")
