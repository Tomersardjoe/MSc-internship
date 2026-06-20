#!/usr/bin/env python3

import sys

def parse_domtbl(domtbl_path, evalue_thresh, coverage_thresh, output_path):
    """
    Parse an HMMER domtblout file and filter hits by E-value and HMM coverage.

    Parameters:
        domtbl_path (str): Path to the domains.tbl file
        evalue_thresh (float): Maximum allowed E-value for the hit
        coverage_thresh (float): Minimum fraction of HMM covered (0-1)
        output_path (str): Path to save filtered hits
    """
    filtered_hits = []
    header_lines = []
    
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                continue
    
            parts = line.strip().split()
            if len(parts) < 23:
                continue
    
            try:
                full_seq_evalue = float(parts[6])
                hmm_from = int(parts[17])
                hmm_to   = int(parts[18])
            except ValueError:
                continue
    
            hmm_length = hmm_to - hmm_from + 1
            coverage = hmm_length / 273.0  # HMM length
    
            if full_seq_evalue <= evalue_thresh and coverage >= coverage_thresh:
                filtered_hits.append(line)
    
    # Write headers first, then filtered hits
    with open(output_path, 'w') as out:
        for h in header_lines:
            out.write(h)
        for hit in filtered_hits:
            out.write(hit)

    print(f"Filtered {len(filtered_hits)} hits out of total lines in {domtbl_path}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <domains.tbl> <evalue_thresh> <coverage_thresh> <output.tbl>")
        sys.exit(1)

    domtbl_path = sys.argv[1]
    evalue_thresh = float(sys.argv[2])
    coverage_thresh = float(sys.argv[3])
    output_path = sys.argv[4]

    parse_domtbl(domtbl_path, evalue_thresh, coverage_thresh, output_path)