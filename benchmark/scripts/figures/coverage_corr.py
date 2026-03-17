#!/usr/bin/env python3
"""
coverage_corr.py - Compute coverage profile correlation between aligners

Stage 2 metric: Pearson correlation of per-position coverage profiles.

Usage:
    samtools depth -a aligner1.bam > aligner1.depth
    samtools depth -a aligner2.bam > aligner2.depth
    python3 coverage_corr.py aligner1.depth aligner2.depth
"""

import sys
import numpy as np

def parse_depth(filename):
    """Parse samtools depth output into position->coverage dict."""
    coverage = {}
    with open(filename) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                pos = int(parts[1])
                depth = int(parts[2])
                coverage[(chrom, pos)] = depth
    return coverage

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 coverage_corr.py depth1.txt depth2.txt")
        sys.exit(1)

    depth1 = parse_depth(sys.argv[1])
    depth2 = parse_depth(sys.argv[2])

    # Find common positions (sorted for deterministic results)
    common_pos = sorted(set(depth1) & set(depth2))

    # Diagnostic: check overlap percentage
    print(f"Positions depth1: {len(depth1):,}")
    print(f"Positions depth2: {len(depth2):,}")
    print(f"Positions compared: {len(common_pos):,}")

    if len(depth1) > 0 and len(depth2) > 0:
        overlap_pct = 100 * len(common_pos) / min(len(depth1), len(depth2))
        print(f"Overlap: {overlap_pct:.1f}% of min()")
        if overlap_pct < 50:
            print("WARNING: Low overlap - check BAM/depth files for issues", file=sys.stderr)

    if len(common_pos) < 2:
        print("Error: Need at least 2 common positions for correlation")
        sys.exit(1)

    # Extract values (float64 for numerical stability)
    vals1 = np.array([depth1[p] for p in common_pos], dtype=np.float64)
    vals2 = np.array([depth2[p] for p in common_pos], dtype=np.float64)

    print(f"Mean coverage 1: {np.mean(vals1):.2f}")
    print(f"Mean coverage 2: {np.mean(vals2):.2f}")

    # Pearson correlation using numpy (robust implementation)
    r = np.corrcoef(vals1, vals2)[0, 1]
    print(f"Pearson correlation: {r:.4f}")

if __name__ == '__main__':
    main()
