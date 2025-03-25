#!/usr/bin/env python3

import subprocess
import sys
import os

def run_minimap2_indexing(reference_fa, threads, index_size):
    """Create minimap2 index for short-read mapping."""
    minimap_path = os.environ.get("MINIMAP2", "minimap2")
    index_output = f"{reference_fa}.HR.SR.mmi"

    command = [
        minimap_path,
        "-I", str(index_size),
        "-x", "sr",
        "-t", str(threads),
        "-d", index_output,
        reference_fa
    ]

    print(f"[HaploChi] Running SR minimap2 indexing:\n  {' '.join(command)}")
    subprocess.run(command, check=True)


def run_minimap2_sr_mapping(reference_fa, input_fa, threads, index_size):
    """Map short reads and generate sorted BAM."""
    minimap_path = os.environ.get("MINIMAP2", "minimap2")
    samtools = os.environ.get("SAMTOOLS", "samtools")

    output_index = f"{reference_fa}.HR.SR.mmi"
    filtered_bam = f"{input_fa}.filtered.bam"
    sorted_bam = f"{input_fa}.sorted.bam"

    minimap2_cmd = [
        minimap_path,
        "-ax", "sr",
        "-t", str(threads),
        "-I", str(index_size),
        output_index,
        input_fa,
        "--secondary=no"
    ]

    samtools_view_cmd = [
        samtools,
        "view", "-q", "10", "-b", "-o", filtered_bam
    ]

    samtools_sort_cmd = [
        samtools,
        "sort", "-o", sorted_bam
    ]

    print(f"[HaploChi] Running SR mapping:\n  {' '.join(minimap2_cmd)}")
    minimap_proc = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)

    print(f"[HaploChi] Filtering BAM:\n  {' '.join(samtools_view_cmd)}")
    samtools_view_proc = subprocess.Popen(samtools_view_cmd, stdin=minimap_proc.stdout, stdout=subprocess.PIPE)

    print(f"[HaploChi] Sorting BAM:\n  {' '.join(samtools_sort_cmd)}")
    samtools_sort_proc = subprocess.Popen(samtools_sort_cmd, stdin=samtools_view_proc.stdout)

    # Close pipes and wait
    minimap_proc.stdout.close()
    samtools_view_proc.stdout.close()
    samtools_sort_proc.communicate()

    print(f"[HaploChi] Mapping and BAM sorting done: {sorted_bam}")
    return sorted_bam


def process_bam_to_gff(bam_file, gff_output):
    """Convert sorted BAM to simple GFF format."""
    samtools = os.environ.get("SAMTOOLS", "samtools")

    print(f"[HaploChi] Generating GFF from: {bam_file}")
    with open(gff_output, 'w') as gff_out:
        view_cmd = [samtools, "view", bam_file]
        proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            line = line.decode('utf-8')
            if line.startswith('@'):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) >= 4:
                gff_line = "\t".join([
                    fields[2], ".", "read", fields[3], fields[3], ".", ".", ".", f"ID={fields[0]}"
                ])
                gff_out.write(gff_line + "\n")

        proc.stdout.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python haplochi.py SRmapping <reference.fa> <input1.fa> [<input2.fa> ...] [t=threads] [I=index_size]")
        sys.exit(1)

    reference_fa = sys.argv[1]
    input_files = []
    threads = 8
    index_size = "8G"

    for arg in sys.argv[2:]:
        if arg.startswith("t="):
            threads = int(arg.split("=")[1])
        elif arg.startswith("I="):
            index_size = arg.split("=")[1]
        else:
            input_files.append(arg)

    run_minimap2_indexing(reference_fa, threads, index_size)

    for input_fa in input_files:
        sorted_bam = run_minimap2_sr_mapping(reference_fa, input_fa, threads, index_size)
        gff_file = f"{input_fa}.SR.gff"
        process_bam_to_gff(sorted_bam, gff_file)
