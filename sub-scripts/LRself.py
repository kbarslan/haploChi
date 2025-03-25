#!/usr/bin/env python3

import subprocess
import sys
import os

def run_minimap2_indexing(reference_fa, threads, index_size):
    """Create minimap2 index for long reads."""
    minimap_path = os.environ.get("MINIMAP2", "minimap2")  # use env var or system PATH
    output_index = f"{reference_fa}.HR.LR.mmi"

    command = [
        minimap_path,
        "-x", "map-hifi",
        "-I", str(index_size),
        "-t", str(threads),
        "-d", output_index,
        reference_fa
    ]

    print(f"[HaploChi] Running minimap2 indexing:\n  {' '.join(command)}")
    subprocess.run(command, check=True)


def run_minimap2_lr_mapping(reference_fa, input_fq, threads, index_size):
    """Map long reads and detect chimeras."""
    minimap_path = os.environ.get("MINIMAP2", "minimap2")
    samtools = os.environ.get("SAMTOOLS", "samtools")

    output_index = f"{reference_fa}.HR.LR.mmi"
    raw_bam = os.path.abspath(f"{input_fq}.raw.bam")
    chimera_bam = f"{input_fq}.chimera.bam"
    chimera_ids_file = f"{input_fq}.chimera.ids"

    # Minimap2 command
    minimap2_cmd = [
        minimap_path,
        "-ax", "map-hifi",
        "-I", str(index_size),
        "-t", str(threads),
        output_index,
        input_fq,
        "--secondary=no"
    ]

    print(f"[HaploChi] Running minimap2 mapping:\n  {' '.join(minimap2_cmd)}")
    minimap_proc = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)

    # Samtools sort command
    samtools_sort_cmd = [samtools, "sort", "-o", raw_bam]
    print(f"[HaploChi] Sorting BAM:\n  {' '.join(samtools_sort_cmd)}")
    samtools_sort_proc = subprocess.Popen(samtools_sort_cmd, stdin=minimap_proc.stdout)
    minimap_proc.stdout.close()
    samtools_sort_proc.communicate()

    # Chimera detection (flag 0x800)
    samtools_view_cmd = [
        samtools,
        "view", "-f", "0x800", "-b", "-o", chimera_bam, raw_bam
    ]
    print(f"[HaploChi] Extracting chimeric reads:\n  {' '.join(samtools_view_cmd)}")
    subprocess.run(samtools_view_cmd, check=True)

    # Extract chimera read IDs
    with open(chimera_ids_file, 'w') as output_file:
        view_cmd = [samtools, "view", chimera_bam]
        print(f"[HaploChi] Writing chimera read IDs to: {chimera_ids_file}")
        with subprocess.Popen(view_cmd, stdout=subprocess.PIPE) as p:
            for line in p.stdout:
                read_id = line.decode('utf-8').split('\t')[0]
                output_file.write(read_id + '\n')
            p.stdout.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python haplochi.py LRmapping <reference.fa> <reads.fastq.gz> [t=threads] [I=index_size]")
        sys.exit(1)

    reference_fa = sys.argv[1]
    input_fq = sys.argv[2]
    threads = 8
    index_size = "8G"

    # Optional args: t=THREADS I=INDEXSIZE
    for arg in sys.argv[3:]:
        if arg.startswith("t="):
            threads = int(arg.split("=")[1])
        elif arg.startswith("I="):
            index_size = arg.split("=")[1]

    run_minimap2_indexing(reference_fa, threads, index_size)
    run_minimap2_lr_mapping(reference_fa, input_fq, threads, index_size)
