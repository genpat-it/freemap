#!/usr/bin/env python3
"""
Calculate coverage correlation on CONCORDANT reads only.
V2: Actually filters reads and computes concordant-only coverage.
"""

import subprocess
import sys
import tempfile
import os
import numpy as np
from scipy import stats

def get_read_info(bam_file, chrom):
    """Extract read name -> (pos, cigar_len) from BAM."""
    info = {}
    cmd = f"samtools view {bam_file} {chrom}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 6:
            continue

        name = fields[0]
        flag = int(fields[1])
        pos = int(fields[3])
        cigar = fields[5]

        # Skip unmapped and secondary
        if flag & 0x4 or flag & 0x100:
            continue

        # Calculate aligned length from CIGAR
        import re
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        ref_len = sum(int(n) for n, op in ops if op in 'MDN=X')

        info[name] = (pos, ref_len, cigar)

    return info

def find_concordant_reads(info1, info2, threshold=500):
    """Find reads that map to same position (within threshold)."""
    concordant = {}
    discordant_pos = []
    discordant_chr = []

    for name in info1:
        if name in info2:
            p1, l1, c1 = info1[name]
            p2, l2, c2 = info2[name]
            if abs(p1 - p2) <= threshold:
                concordant[name] = (p1, l1, c1, p2, l2, c2)
            else:
                discordant_pos.append(name)
        else:
            discordant_chr.append(name)

    return concordant, discordant_pos, discordant_chr

def compute_coverage_subset(bam_file, chrom, read_names, tmpdir):
    """Compute coverage for a subset of reads."""
    # Write read names to file
    namefile = os.path.join(tmpdir, "reads.txt")
    with open(namefile, 'w') as f:
        for name in read_names:
            f.write(name + '\n')

    # Filter BAM and compute depth
    filtered_bam = os.path.join(tmpdir, "filtered.bam")
    cmd = f"samtools view -h {bam_file} {chrom} | grep -F -f {namefile} -e '^@' | samtools view -b -o {filtered_bam}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    # Get depth
    cmd = f"samtools depth -a {filtered_bam}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    coverage = {}
    for line in result.stdout.strip().split('\n'):
        if line:
            fields = line.split('\t')
            if len(fields) >= 3:
                pos = int(fields[1])
                depth = int(fields[2])
                coverage[pos] = depth

    return coverage

def coverage_to_windows(coverage, window_size=1000):
    """Aggregate coverage into windows."""
    if not coverage:
        return {}

    max_pos = max(coverage.keys())
    windows = {}

    for start in range(1, max_pos + 1, window_size):
        end = start + window_size
        depths = [coverage.get(p, 0) for p in range(start, end)]
        if depths and sum(depths) > 0:  # Only non-zero windows
            windows[start] = np.mean(depths)

    return windows

def main():
    if len(sys.argv) < 4:
        print("Usage: concordant_coverage_v2.py <bam1> <bam2> <chrom>")
        sys.exit(1)

    bam1 = sys.argv[1]
    bam2 = sys.argv[2]
    chrom = sys.argv[3]

    print(f"=== Concordant Coverage Analysis ===")
    print(f"Chromosome: {chrom}")
    print(f"BAM1 (freemap): {bam1}")
    print(f"BAM2 (minimap2): {bam2}")
    print()

    # Get read info
    print("Extracting read positions...")
    info1 = get_read_info(bam1, chrom)
    info2 = get_read_info(bam2, chrom)
    print(f"  BAM1: {len(info1)} reads")
    print(f"  BAM2: {len(info2)} reads")

    # Find concordant
    concordant, disc_pos, disc_chr = find_concordant_reads(info1, info2)
    total_common = len(concordant) + len(disc_pos)
    print(f"  Concordant (same chr, ±500bp): {len(concordant)} ({100*len(concordant)/len(info1):.1f}%)")
    print(f"  Discordant (position >500bp): {len(disc_pos)}")
    print(f"  Different chromosome: {len(disc_chr)}")
    print()

    with tempfile.TemporaryDirectory() as tmpdir:
        # ALL reads coverage
        print("Computing ALL reads coverage...")
        cmd1 = f"samtools depth -a -r {chrom} {bam1}"
        result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
        cov1_all = {}
        for line in result1.stdout.strip().split('\n'):
            if line:
                f = line.split('\t')
                if len(f) >= 3:
                    cov1_all[int(f[1])] = int(f[2])

        cmd2 = f"samtools depth -a -r {chrom} {bam2}"
        result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
        cov2_all = {}
        for line in result2.stdout.strip().split('\n'):
            if line:
                f = line.split('\t')
                if len(f) >= 3:
                    cov2_all[int(f[1])] = int(f[2])

        win1_all = coverage_to_windows(cov1_all)
        win2_all = coverage_to_windows(cov2_all)

        common = set(win1_all.keys()) & set(win2_all.keys())
        v1 = [win1_all[w] for w in sorted(common)]
        v2 = [win2_all[w] for w in sorted(common)]
        r_all, _ = stats.pearsonr(v1, v2)
        print(f"  Coverage correlation (all reads): r = {r_all:.3f}")
        print()

        # CONCORDANT reads only
        print("Computing CONCORDANT reads coverage...")
        concordant_names = set(concordant.keys())

        cov1_conc = compute_coverage_subset(bam1, chrom, concordant_names, tmpdir)
        cov2_conc = compute_coverage_subset(bam2, chrom, concordant_names, tmpdir)

        win1_conc = coverage_to_windows(cov1_conc)
        win2_conc = coverage_to_windows(cov2_conc)

        common_conc = set(win1_conc.keys()) & set(win2_conc.keys())
        if common_conc:
            v1c = [win1_conc[w] for w in sorted(common_conc)]
            v2c = [win2_conc[w] for w in sorted(common_conc)]
            r_conc, _ = stats.pearsonr(v1c, v2c)
            print(f"  Coverage correlation (concordant only): r = {r_conc:.3f}")
        else:
            r_conc = float('nan')
            print("  No common windows for concordant reads")

    print()
    print("=== SUMMARY ===")
    print(f"Coverage correlation (all reads, 1kb):       r = {r_all:.3f}")
    print(f"Coverage correlation (concordant only, 1kb): r = {r_conc:.3f}")
    print(f"Improvement from isolating CIGAR quality:    Δr = {r_conc - r_all:+.3f}")

if __name__ == "__main__":
    main()
