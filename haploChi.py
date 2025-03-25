#!/usr/bin/env python3

"""
HaploChi – Haplotype Chimera Inspector

Modular tool to assess haplotype resolution quality via:
- read splitting
- long/short-read mapping
- GFF filtering
- chimera analysis
- per-haplotype statistics


Author: Kubra Arslan
License: MIT
"""

import os
import sys
import subprocess

# Replace this with your actual version function once available
# from utilities import get_haplochi_version

VERSION = "0.1.0"  # Temporary; link to `get_haplochi_version()` once done

CITATION = """
HaploChi: Haplotype Chimera Inspector
Kubra Arslan et al. (A&G, 2025 hopefully)
https://github.com/yourusername/haplochi
"""


def print_help():
    print(f"""
HaploChi – Haplotype Chimera Inspector
Version: {VERSION}

Usage:
  python haplochi.py <command> [options]

Commands:
  split         Split long FASTQ reads into 200bp subreads
  LRmapping     Long-read (HiFi) mapping and chimera detection
  SRmapping     Short-read mapping and GFF creation
  process       Extract per-haplotype stats and generate metafiles

Options:
  -h, --help        Show this help message
  -v, --version     Show current version
  -c, --citation    Show citation info

Examples:
  python haplochi.py split my_reads.fastq.gz
  python haplochi.py LRmapping reference.fa my_reads.fastq.gz
  python haplochi.py SRmapping reference.fa split1.fa split2.fa
  python haplochi.py process reference.fa sorg_49
""")


def main():
    if len(sys.argv) < 2:
        print_help()
        return

    cmd = sys.argv[1]

    if cmd in ("-h", "--help"):
        print_help()

    elif cmd in ("-v", "--version"):
        # print(get_haplochi_version())
        print(VERSION)

    elif cmd in ("-c", "--citation"):
        print(CITATION)

    elif cmd == "split":
        subcmd = ["python3", "pipeline_scripts/split_v1.py"] + sys.argv[2:]
        subprocess.call(subcmd)

    elif cmd == "LRmapping":
        subcmd = ["python3", "pipeline_scripts/LRself.py"] + sys.argv[2:]
        subprocess.call(subcmd)

    elif cmd == "SRmapping":
        subcmd = ["python3", "pipeline_scripts/SRself1.py"] + sys.argv[2:]
        subprocess.call(subcmd)

    elif cmd == "process":
        subcmd = ["python3", "pipeline_scripts/merged_data_man.py"] + sys.argv[2:]
        subprocess.call(subcmd)

    else:
        print(f"[HaploChi] Unrecognized command: {cmd}")
        print_help()


if __name__ == "__main__":
    main()
