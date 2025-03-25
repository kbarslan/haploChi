#!/usr/bin/env python3

import sys
import os
import gzip
import glob
from Bio import SeqIO

MAX_RECORDS_PER_FILE = 200000  # Number of reads per split file
CHUNK_SIZE = 200               # Size of each sub-read

def read_file_list(file_list_path):
    """Read list of file paths from a .list file."""
    with open(file_list_path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def print_usage():
    print("""
Usage:
    python haplochi.py split <file1.fastq.gz> [file2.fastq.gz ...]
    python haplochi.py split *.fastq.gz
    python haplochi.py split <fastq.list>

Description:
    Splits FASTQ (gzipped) files into smaller chunks of subreads.
    Each output file contains up to 200,000 reads. Long reads are split
    into 200bp subreads, saved in FASTA format.
""")


def split_reads(input_file):
    """Split reads in input FASTQ file into subreads, save in .fa files."""
    record_count = 0
    file_count = 1
    base_name = ".".join(input_file.split(".")[:-1])
    output_file = None
    fh = None

    print(f"[HaploChi] Processing: {input_file}")

    with gzip.open(input_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Create a new output file after reaching MAX_RECORDS_PER_FILE
            if record_count >= MAX_RECORDS_PER_FILE:
                record_count = 0
                file_count += 1
                if fh:
                    fh.close()
                fh = None

            if record_count == 0:
                output_file = f"{base_name}.split{file_count}.fa"
                fh = open(output_file, "w")
                print(f"[HaploChi] Writing to: {output_file}")

            rid = record.id
            rseq = record.seq

            chunks = [rseq[i:i + CHUNK_SIZE] for i in range(0, len(rseq), CHUNK_SIZE)]
            for i, chunk in enumerate(chunks, 1):
                if len(chunk) == CHUNK_SIZE:
                    fh.write(f">{rid}:sr{i}\n{chunk}\n")

            record_count += 1

    if fh:
        fh.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(1)

    # Expand file inputs
    input_files = []

    if sys.argv[1].endswith(".list"):
        input_files = read_file_list(sys.argv[1])
    else:
        input_files = sys.argv[1:]

    # Glob wildcards if used (e.g. *.fastq.gz)
    expanded_files = []
    for pattern in input_files:
        expanded_files.extend(glob.glob(pattern))
    input_files = expanded_files

    if not input_files:
        print("[HaploChi] No FASTQ files found.")
        sys.exit(1)

    for f in input_files:
        split_reads(f)

    print("[HaploChi] Splitting complete.")
