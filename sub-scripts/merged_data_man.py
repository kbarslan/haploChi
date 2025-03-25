#!/usr/bin/env python3

import os
import sys
import numpy as np

def process_gff_files(ref_fa, prefix):
    """Split SR GFF file by haplotype and filter with contig names."""
    sr_gff_file = f"{prefix}.SR.gff"
    h1_gff_file = f"{prefix}.SR.h1.gff"
    h2_gff_file = f"{prefix}.SR.h2.gff"
    h1_contigs_file = f"{prefix}.h1.contigs"
    h2_contigs_file = f"{prefix}.h2.contigs"
    h1_nc_gff_file = f"{prefix}.h1.nc.gff"
    h2_nc_gff_file = f"{prefix}.h2.nc.gff"

    print(f"[HaploChi] Splitting {sr_gff_file} into haplotypes...")

    with open(sr_gff_file, 'r') as sr_gff, \
         open(h1_gff_file, 'w') as h1_gff, \
         open(h2_gff_file, 'w') as h2_gff:
        for line in sr_gff:
            if line.startswith("h1"):
                h1_gff.write(line)
            elif line.startswith("h2"):
                h2_gff.write(line)

    print(f"[HaploChi] Extracting contigs from {ref_fa}...")

    with open(ref_fa, 'r') as fasta, \
         open(h1_contigs_file, 'w') as h1_contigs, \
         open(h2_contigs_file, 'w') as h2_contigs:
        for line in fasta:
            if line.startswith(">h1"):
                h1_contigs.write(line[1:])
            elif line.startswith(">h2"):
                h2_contigs.write(line[1:])

    def filter_gff_by_contigs(gff_file, contigs_file, output_file):
        contigs = set(line.strip() for line in open(contigs_file))
        with open(gff_file, 'r') as gff_in, open(output_file, 'w') as gff_out:
            for line in gff_in:
                if line.split()[0] in contigs:
                    gff_out.write(line)

    print("[HaploChi] Filtering GFFs by contig names...")
    filter_gff_by_contigs(h1_gff_file, h1_contigs_file, h1_nc_gff_file)
    filter_gff_by_contigs(h2_gff_file, h2_contigs_file, h2_nc_gff_file)


def get_counts(ref_fa, prefix):
    """Read updated GFFs and chimera IDs, compute counts and metadata."""
    updated_h1_gff = f"{prefix}.h1.nc.u.gff"
    updated_h2_gff = f"{prefix}.h2.nc.u.gff"
    chimera_ids_file = f"{prefix}.chimera.ids"

    print("[HaploChi] Parsing GFF and chimera files...")

    a_s = set()
    h1_c, h2_c = {}, {}
    h1_p, h2_p = {}, {}
    h1_ctg, h2_ctg = {}, {}

    # Parse H1 GFF
    with open(updated_h1_gff) as f:
        for line in f:
            arr = line.strip().split("\t")
            read_id = arr[8].split("=")[1].split(":")[0]
            a_s.add(read_id)
            h1_c.setdefault(read_id, 0)
            h1_p.setdefault(read_id, [])
            h1_ctg.setdefault(read_id, set())
            h1_c[read_id] += 1
            h1_p[read_id].append(int(arr[3]))
            h1_ctg[read_id].add(arr[0])

    # Parse H2 GFF
    with open(updated_h2_gff) as f:
        for line in f:
            arr = line.strip().split("\t")
            read_id = arr[8].split("=")[1].split(":")[0]
            a_s.add(read_id)
            h2_c.setdefault(read_id, 0)
            h2_p.setdefault(read_id, [])
            h2_ctg.setdefault(read_id, set())
            h2_c[read_id] += 1
            h2_p[read_id].append(int(arr[3]))
            h2_ctg[read_id].add(arr[0])

    # Chimera ID list
    h_s = set(line.strip() for line in open(chimera_ids_file))

    output = []
    for read_id in a_s:
        chimera_status = "CHIMERA" if read_id in h_s else "NOCHIMERA"

        if read_id in h1_c and read_id in h2_c:
            output.append(f"{read_id}\t{h1_c[read_id]}\t{h2_c[read_id]}\tBOTH\t"
                          f"{','.join(h1_ctg[read_id])}\t{int(np.median(h1_p[read_id]))}\t{max(h1_p[read_id]) - min(h1_p[read_id])}\t"
                          f"{','.join(h2_ctg[read_id])}\t{int(np.median(h2_p[read_id]))}\t{max(h2_p[read_id]) - min(h2_p[read_id])}\t"
                          f"{chimera_status}")
        elif read_id in h1_c:
            output.append(f"{read_id}\t{h1_c[read_id]}\t0\tH1\t"
                          f"{','.join(h1_ctg[read_id])}\t{int(np.median(h1_p[read_id]))}\t{max(h1_p[read_id]) - min(h1_p[read_id])}\t"
                          f"NA\t0\t0\t{chimera_status}")
        elif read_id in h2_c:
            output.append(f"{read_id}\t0\t{h2_c[read_id]}\tH2\t"
                          f"NA\t0\t0\t{','.join(h2_ctg[read_id])}\t{int(np.median(h2_p[read_id]))}\t{max(h2_p[read_id]) - min(h2_p[read_id])}\t"
                          f"{chimera_status}")
    return output


def filter_counts(counts):
    """Apply filters for high-confidence H1/H2 calls."""
    filtered_output = []
    for line in counts:
        arr = line.strip().split("\t")
        try:
            h1_count = int(arr[1])
            h2_count = int(arr[2])
            h1_dist = int(arr[6])
            h2_dist = int(arr[9])
        except (ValueError, IndexError):
            continue

        if h1_count >= 5 and h2_count >= 5 and h1_dist < 100000 and h2_dist < 100000:
            ratio = h1_count / h2_count if h2_count != 0 else 1
            log_ratio = f"{np.log10(ratio):.2f}"
            dominant = "DH1" if h1_count > h2_count else "DH2"
            arr[5] = str(round(float(arr[5])))  # median h1
            arr[8] = str(round(float(arr[8])))  # median h2
            filtered_output.append("\t".join(arr) + f"\t{log_ratio}\t{dominant}")

    return filtered_output


def main(ref_fa, prefix):
    print(f"[HaploChi] Processing data for prefix: {prefix}")
    process_gff_files(ref_fa, prefix)
    counts = get_counts(ref_fa, prefix)

    raw_counts_file = f"{prefix}.chr.counts"
    with open(raw_counts_file, 'w') as out:
        for line in counts:
            out.write(line + "\n")

    print(f"[HaploChi] Wrote raw counts: {raw_counts_file}")

    filtered = filter_counts(counts)
    filtered_file = f"{prefix}.chr.counts.min5.dist100000.h1h2"
    with open(filtered_file, 'w') as out:
        for line in filtered:
            out.write(line + "\n")

    print(f"[HaploChi] Wrote filtered counts: {filtered_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python haplochi.py process <reference.fa> <prefix>")
        sys.exit(1)

    ref_fa = sys.argv[1]
    prefix = sys.argv[2]
    main(ref_fa, prefix)
