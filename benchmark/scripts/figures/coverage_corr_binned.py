#!/usr/bin/env python3
"""
coverage_corr_binned.py - Compute coverage correlation with binning

More representative of downstream usage: most analyses work on windows,
not per-base resolution. Binning reduces noise from local CIGAR differences.

Usage:
    python3 coverage_corr_binned.py depth1.txt depth2.txt [bin_size]
"""

import sys
import numpy as np
from collections import defaultdict

def parse_depth_binned(filename, bin_size=1000):
    """Parse samtools depth -a output into binned coverage.

    Note: Requires 'samtools depth -a' (all positions including zero).
    POS is 1-based, so we use (pos-1)//bin_size for correct binning.
    """
    sum_cov = defaultdict(float)
    count = defaultdict(int)

    with open(filename) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t', 3)  # split only first 3 fields
            if len(parts) >= 3:
                chrom = parts[0]
                pos = int(parts[1])
                depth = int(parts[2])
                # pos is 1-based, convert to 0-based for binning
                bin_key = (chrom, (pos - 1) // bin_size)
                sum_cov[bin_key] += depth
                count[bin_key] += 1

    # Convert to mean coverage per bin
    binned = {k: sum_cov[k] / count[k] for k in sum_cov}

    # Warn if bins seem incomplete (suggests missing -a flag)
    n_short = sum(1 for k in count if count[k] < bin_size * 0.9)
    if n_short > len(count) * 0.1:  # >10% bins are short
        print(f"Warning: {n_short}/{len(count)} bins have <90% expected positions; "
              "did you use 'samtools depth -a'?", file=sys.stderr)

    return binned

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 coverage_corr_binned.py depth1.txt depth2.txt [bin_size]")
        sys.exit(1)

    bin_size = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

    print(f"Bin size: {bin_size} bp")

    depth1 = parse_depth_binned(sys.argv[1], bin_size)
    depth2 = parse_depth_binned(sys.argv[2], bin_size)

    # Find common bins (sorted for determinism)
    common_bins = sorted(set(depth1) & set(depth2))

    print(f"Bins depth1: {len(depth1):,}")
    print(f"Bins depth2: {len(depth2):,}")
    print(f"Bins compared: {len(common_bins):,}")

    if len(depth1) > 0 and len(depth2) > 0:
        overlap_pct = 100 * len(common_bins) / min(len(depth1), len(depth2))
        print(f"Overlap: {overlap_pct:.1f}% of min()")

    if len(common_bins) < 2:
        print("Error: Need at least 2 common bins for correlation")
        sys.exit(1)

    vals1 = np.array([depth1[b] for b in common_bins], dtype=np.float64)
    vals2 = np.array([depth2[b] for b in common_bins], dtype=np.float64)

    print(f"Mean coverage 1: {np.mean(vals1):.2f}")
    print(f"Mean coverage 2: {np.mean(vals2):.2f}")

    # Pearson correlation (handle zero variance edge case)
    if np.std(vals1) == 0 or np.std(vals2) == 0:
        print("Pearson correlation (binned): N/A (zero variance)")
    else:
        r = np.corrcoef(vals1, vals2)[0, 1]
        print(f"Pearson correlation (binned): {r:.4f}")

if __name__ == '__main__':
    main()
