#!/usr/bin/env python3
"""
cigar_stats.py - Extract CIGAR statistics from SAM files

Stage 2 metrics:
- Per-read insertion/deletion event counts
- Per-read insertion/deletion base counts
- Per-read CIGAR operation count

Usage:
    python3 cigar_stats.py alignment.sam              # Summary stats
    python3 cigar_stats.py --compare sam1 sam2        # Compare two aligners
"""

import sys
import re
import argparse
import numpy as np


def parse_cigar(cigar):
    """Parse CIGAR string and return (n_ins, n_del, n_ops, ins_bases, del_bases)."""
    if cigar == '*' or not cigar:
        return 0, 0, 0, 0, 0

    n_ins = 0      # number of insertion events
    n_del = 0      # number of deletion events
    n_ops = 0      # total CIGAR operations
    ins_bases = 0  # total inserted bases
    del_bases = 0  # total deleted bases

    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(match.group(1))
        op = match.group(2)
        n_ops += 1

        if op == 'I':
            n_ins += 1
            ins_bases += length
        elif op == 'D':
            n_del += 1
            del_bases += length

    return n_ins, n_del, n_ops, ins_bases, del_bases


def parse_sam(filename, strip_pair_suffix=False, with_position=False):
    """Parse SAM file and extract per-read CIGAR stats.

    Returns:
        read_stats: dict of read_name -> (n_ins, n_del, n_ops, ins_bases, del_bases)
                    or if with_position: dict of read_name -> (chr, pos, n_ins, n_del, n_ops, ins_bases, del_bases)
        n_duplicates: count of skipped duplicate read names
    """
    read_stats = {}
    n_duplicates = 0

    with open(filename) as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue

            try:
                flag = int(fields[1])
            except ValueError:
                continue

            # Skip unmapped
            if flag & 4:
                continue

            # Skip secondary and supplementary
            if flag & 0x900:
                continue

            qname = fields[0]

            # Optional: strip /1 /2 suffix for paired reads
            if strip_pair_suffix and (qname.endswith('/1') or qname.endswith('/2')):
                qname = qname[:-2]

            cigar = fields[5]

            if qname in read_stats:
                n_duplicates += 1
                continue

            cigar_stats = parse_cigar(cigar)

            if with_position:
                chrom = fields[2]
                # Skip unmapped or invalid reference names
                if chrom == '*' or chrom == '':
                    continue
                try:
                    pos = int(fields[3])
                except ValueError:
                    pos = 0
                read_stats[qname] = (chrom, pos) + cigar_stats
            else:
                read_stats[qname] = cigar_stats

    return read_stats, n_duplicates


def safe_pearson(a, b):
    """Compute Pearson correlation, handling edge cases."""
    if len(a) < 2:
        return float('nan')
    if np.std(a) == 0 or np.std(b) == 0:
        return float('nan')
    # Use numpy corrcoef to avoid scipy dependency
    return np.corrcoef(a, b)[0, 1]


def safe_spearman(a, b):
    """Compute Spearman rank correlation without scipy dependency."""
    if len(a) < 2:
        return float('nan')
    # Rank-transform and compute Pearson on ranks
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return safe_pearson(ra, rb)


def compare_aligners_concordant(sam1, sam2, name1, name2, tolerance=500, strip_pair_suffix=False):
    """Compare CIGAR stats only for reads with concordant placement."""
    print(f"Parsing {sam1} (with positions)...", file=sys.stderr)
    stats1, _ = parse_sam(sam1, strip_pair_suffix=strip_pair_suffix, with_position=True)
    print(f"Parsing {sam2} (with positions)...", file=sys.stderr)
    stats2, _ = parse_sam(sam2, strip_pair_suffix=strip_pair_suffix, with_position=True)

    # Check for common reads before proceeding
    common_reads = set(stats1.keys()) & set(stats2.keys())
    if not common_reads:
        print("Error: No common reads found", file=sys.stderr)
        return None

    # Filter to concordant reads (same chr, pos within tolerance)
    concordant = []
    discordant_chr = 0
    discordant_pos = 0
    for r in common_reads:
        chr1, pos1 = stats1[r][0], stats1[r][1]
        chr2, pos2 = stats2[r][0], stats2[r][1]
        if chr1 != chr2:
            discordant_chr += 1
        elif abs(pos1 - pos2) > tolerance:
            discordant_pos += 1
        else:
            concordant.append(r)

    print(f"=== Concordant CIGAR comparison: {name1} vs {name2} ===")
    print(f"Common reads: {len(common_reads):,}")
    print(f"Concordant (same chr, pos ±{tolerance}): {len(concordant):,} ({100*len(concordant)/len(common_reads):.1f}%)")
    print(f"Discordant chr: {discordant_chr:,}, Discordant pos: {discordant_pos:,}")
    print()

    if len(concordant) < 10:
        print("Too few concordant reads for correlation")
        return None

    # Extract CIGAR stats for concordant reads (indices 2-6 for with_position)
    ins_bases1 = np.array([stats1[r][5] for r in concordant])
    ins_bases2 = np.array([stats2[r][5] for r in concordant])
    del_bases1 = np.array([stats1[r][6] for r in concordant])
    del_bases2 = np.array([stats2[r][6] for r in concordant])
    ops1 = np.array([stats1[r][4] for r in concordant])
    ops2 = np.array([stats2[r][4] for r in concordant])

    # Event counts for concordant reads
    ins_events1 = np.array([stats1[r][2] for r in concordant])
    ins_events2 = np.array([stats2[r][2] for r in concordant])
    del_events1 = np.array([stats1[r][3] for r in concordant])
    del_events2 = np.array([stats2[r][3] for r in concordant])

    # Correlations
    ins_pearson = safe_pearson(ins_bases1, ins_bases2)
    del_pearson = safe_pearson(del_bases1, del_bases2)
    ins_spearman = safe_spearman(ins_bases1, ins_bases2)
    del_spearman = safe_spearman(del_bases1, del_bases2)
    ins_evt_spearman = safe_spearman(ins_events1, ins_events2)
    del_evt_spearman = safe_spearman(del_events1, del_events2)

    print("Base count correlations (concordant only):")
    print(f"  Inserted bases (Pearson):  {ins_pearson:.4f}")
    print(f"  Deleted bases (Pearson):   {del_pearson:.4f}")
    print(f"  Inserted bases (Spearman): {ins_spearman:.4f}")
    print(f"  Deleted bases (Spearman):  {del_spearman:.4f}")
    print()
    print("Event count correlations (concordant only):")
    print(f"  Insertion events (Spearman): {ins_evt_spearman:.4f}")
    print(f"  Deletion events (Spearman):  {del_evt_spearman:.4f}")
    print()

    # Global totals
    total_ins1, total_ins2 = np.sum(ins_bases1), np.sum(ins_bases2)
    total_del1, total_del2 = np.sum(del_bases1), np.sum(del_bases2)
    print("Global totals (concordant reads):")
    print(f"  Total ins bases ({name1}): {total_ins1:,}")
    print(f"  Total ins bases ({name2}): {total_ins2:,} (ratio: {total_ins1/total_ins2:.2f})")
    print(f"  Total del bases ({name1}): {total_del1:,}")
    print(f"  Total del bases ({name2}): {total_del2:,} (ratio: {total_del1/total_del2:.2f})")
    print()
    print(f"Mean ops/read ({name1}): {np.mean(ops1):.1f}")
    print(f"Mean ops/read ({name2}): {np.mean(ops2):.1f}")

    return {
        'concordant': len(concordant),
        'ins_pearson': ins_pearson,
        'del_pearson': del_pearson,
        'ins_spearman': ins_spearman,
        'del_spearman': del_spearman,
    }


def compute_stats(read_stats):
    """Compute summary statistics from read stats."""
    if not read_stats:
        return {}

    n_ins_list = [s[0] for s in read_stats.values()]
    n_del_list = [s[1] for s in read_stats.values()]
    n_ops_list = [s[2] for s in read_stats.values()]
    ins_bases_list = [s[3] for s in read_stats.values()]
    del_bases_list = [s[4] for s in read_stats.values()]

    return {
        'n_reads': len(read_stats),
        'mean_ins_events': np.mean(n_ins_list),
        'mean_del_events': np.mean(n_del_list),
        'mean_ops': np.mean(n_ops_list),
        'median_ops': np.median(n_ops_list),
        'mean_ins_bases': np.mean(ins_bases_list),
        'mean_del_bases': np.mean(del_bases_list),
        'total_ins_events': sum(n_ins_list),
        'total_del_events': sum(n_del_list),
        'total_ins_bases': sum(ins_bases_list),
        'total_del_bases': sum(del_bases_list),
    }


def compare_aligners(stats1, stats2, name1, name2, dup1, dup2):
    """Compare CIGAR stats between two aligners."""
    # Find common reads
    common_reads = set(stats1.keys()) & set(stats2.keys())

    if len(common_reads) == 0:
        print("Error: No common reads found", file=sys.stderr)
        return None

    # Extract values for common reads
    ins1 = np.array([stats1[r][0] for r in common_reads])
    ins2 = np.array([stats2[r][0] for r in common_reads])
    del1 = np.array([stats1[r][1] for r in common_reads])
    del2 = np.array([stats2[r][1] for r in common_reads])
    ops1 = np.array([stats1[r][2] for r in common_reads])
    ops2 = np.array([stats2[r][2] for r in common_reads])
    ins_bases1 = np.array([stats1[r][3] for r in common_reads])
    ins_bases2 = np.array([stats2[r][3] for r in common_reads])
    del_bases1 = np.array([stats1[r][4] for r in common_reads])
    del_bases2 = np.array([stats2[r][4] for r in common_reads])

    # Compute correlations
    ins_corr = safe_pearson(ins1, ins2)
    del_corr = safe_pearson(del1, del2)
    ops_corr = safe_pearson(ops1, ops2)
    ins_bases_corr = safe_pearson(ins_bases1, ins_bases2)
    del_bases_corr = safe_pearson(del_bases1, del_bases2)

    print(f"=== CIGAR comparison: {name1} vs {name2} ===")
    print(f"Common reads: {len(common_reads):,}")
    if dup1 > 0:
        print(f"  (skipped {dup1:,} duplicates in {name1})")
    if dup2 > 0:
        print(f"  (skipped {dup2:,} duplicates in {name2})")
    print()
    print("Event count correlations:")
    print(f"  Insertion events:  {ins_corr:.4f}")
    print(f"  Deletion events:   {del_corr:.4f}")
    print()
    print("Base count correlations:")
    print(f"  Inserted bases:    {ins_bases_corr:.4f}")
    print(f"  Deleted bases:     {del_bases_corr:.4f}")
    print()
    print(f"Operations correlation: {ops_corr:.4f}")
    print()
    print(f"Mean ops/read ({name1}): {np.mean(ops1):.1f}")
    print(f"Mean ops/read ({name2}): {np.mean(ops2):.1f}")
    print()
    print(f"Mean ins events/read ({name1}): {np.mean(ins1):.2f}")
    print(f"Mean ins events/read ({name2}): {np.mean(ins2):.2f}")
    print(f"Mean del events/read ({name1}): {np.mean(del1):.2f}")
    print(f"Mean del events/read ({name2}): {np.mean(del2):.2f}")

    return {
        'common_reads': len(common_reads),
        'ins_event_corr': ins_corr,
        'del_event_corr': del_corr,
        'ins_bases_corr': ins_bases_corr,
        'del_bases_corr': del_bases_corr,
        'ops_corr': ops_corr,
        'mean_ops_1': np.mean(ops1),
        'mean_ops_2': np.mean(ops2),
    }


def main():
    parser = argparse.ArgumentParser(
        description='Extract CIGAR statistics from SAM files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s alignment.sam                    # Summary for one file
  %(prog)s --compare fm.sam mm2.sam         # Compare two aligners
  %(prog)s --compare fm.sam mm2.sam --name1 freemap --name2 minimap2
"""
    )
    parser.add_argument('sam', nargs='+', help='SAM file(s)')
    parser.add_argument('--compare', action='store_true',
                        help='Compare two SAM files (requires exactly 2 files)')
    parser.add_argument('--concordant', action='store_true',
                        help='Only compare reads with concordant placement (same chr, pos ±500bp)')
    parser.add_argument('--name1', default='aligner1', help='Name for first aligner')
    parser.add_argument('--name2', default='aligner2', help='Name for second aligner')
    parser.add_argument('--strip-pair-suffix', action='store_true',
                        help='Strip /1 /2 suffix from read names (for paired-end)')
    parser.add_argument('--json', action='store_true', help='Output JSON format')
    args = parser.parse_args()

    if args.compare:
        if len(args.sam) != 2:
            print("Error: --compare requires exactly 2 SAM files", file=sys.stderr)
            sys.exit(1)

        if args.concordant:
            result = compare_aligners_concordant(args.sam[0], args.sam[1], args.name1, args.name2,
                                                  strip_pair_suffix=args.strip_pair_suffix)
        else:
            print(f"Parsing {args.sam[0]}...", file=sys.stderr)
            stats1, dup1 = parse_sam(args.sam[0], args.strip_pair_suffix)
            print(f"Parsing {args.sam[1]}...", file=sys.stderr)
            stats2, dup2 = parse_sam(args.sam[1], args.strip_pair_suffix)

            result = compare_aligners(stats1, stats2, args.name1, args.name2, dup1, dup2)

        if args.json and result:
            import json
            print(json.dumps(result, indent=2))
    else:
        # Single file mode - output summary stats
        for sam_file in args.sam:
            print(f"Parsing {sam_file}...", file=sys.stderr)
            read_stats, n_dup = parse_sam(sam_file, args.strip_pair_suffix)
            summary = compute_stats(read_stats)

            print(f"=== {sam_file} ===")
            print(f"Reads: {summary['n_reads']:,}")
            if n_dup > 0:
                print(f"  (skipped {n_dup:,} duplicates)")
            print()
            print("Event counts (per read):")
            print(f"  Mean insertions: {summary['mean_ins_events']:.2f}")
            print(f"  Mean deletions:  {summary['mean_del_events']:.2f}")
            print()
            print("Base counts (per read):")
            print(f"  Mean inserted bases: {summary['mean_ins_bases']:.1f}")
            print(f"  Mean deleted bases:  {summary['mean_del_bases']:.1f}")
            print()
            print("CIGAR operations:")
            print(f"  Mean ops/read:   {summary['mean_ops']:.1f}")
            print(f"  Median ops/read: {summary['median_ops']:.1f}")
            print()


if __name__ == '__main__':
    main()
