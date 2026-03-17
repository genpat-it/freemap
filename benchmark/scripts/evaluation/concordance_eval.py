#!/usr/bin/env python3
"""
concordance_eval.py - Compute pairwise concordance between aligners

Two concordance metrics are supported:
  - position (default): Same chromosome AND POS within ±tolerance (default 500bp)
  - overlap: Same chromosome AND alignment overlap >= 10% of max(len1, len2)
  - both: Report both metrics side by side

The overlap metric uses CIGAR strings to compute the reference interval covered
by each alignment, then checks whether the two intervals share sufficient
overlap. This avoids penalizing reads with different soft-clip patterns that
actually cover the same reference region — a common case with ONT long reads.

This is "conditional concordance" - only reads successfully mapped by both
tools are compared. The script reports:
  - Total reads per aligner
  - Common reads (mapped by both)
  - Exclusive reads (mapped by only one)
  - Concordance rate on common reads

Only PRIMARY alignments are considered. Secondary (0x100) and supplementary
(0x800) alignments are filtered out.
"""

import sys
import os
import re
import argparse


def normalize_read_name(qname):
    """Normalize read name by removing /1 or /2 suffix."""
    if qname.endswith('/1') or qname.endswith('/2'):
        return qname[:-2]
    return qname


def cigar_ref_length(cigar):
    """Compute reference length consumed by CIGAR string.

    Operations that consume reference: M, D, N, =, X
    """
    ref_len = 0
    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(match.group(1))
        op = match.group(2)
        if op in 'MDN=X':
            ref_len += length
    return ref_len


def compute_overlap(start1, len1, start2, len2):
    """Compute overlap between two intervals [start, start+len)."""
    end1 = start1 + len1
    end2 = start2 + len2
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def parse_sam(filename):
    """Parse SAM file to extract primary alignment positions and intervals.

    Returns:
        alignments: dict of qname -> (rname, pos, ref_end, strand)
        stats: dict with parsing statistics
    """
    alignments = {}
    stats = {
        'total_lines': 0,
        'unmapped': 0,
        'secondary_supplementary': 0,
        'duplicates_skipped': 0,
        'invalid_rname': 0,
        'parse_errors': 0
    }

    with open(filename) as f:
        for line in f:
            if line.startswith('@'):
                continue

            stats['total_lines'] += 1
            fields = line.rstrip('\n').split('\t')

            if len(fields) < 6:
                stats['parse_errors'] += 1
                continue

            try:
                flag = int(fields[1])
            except ValueError:
                stats['parse_errors'] += 1
                continue

            # Skip unmapped
            if flag & 4:
                stats['unmapped'] += 1
                continue

            # Skip secondary (0x100) and supplementary (0x800)
            if flag & 0x900:
                stats['secondary_supplementary'] += 1
                continue

            qname = normalize_read_name(fields[0])
            rname = fields[2]

            # Skip invalid reference names
            if rname == '*' or rname == '':
                stats['invalid_rname'] += 1
                continue

            try:
                pos = int(fields[3])
            except ValueError:
                stats['parse_errors'] += 1
                continue

            strand = '-' if flag & 16 else '+'

            cigar = fields[5]
            ref_len = cigar_ref_length(cigar) if cigar != '*' else 0
            ref_end = pos + ref_len

            # Skip duplicates (keep first occurrence)
            if qname in alignments:
                stats['duplicates_skipped'] += 1
                continue

            alignments[qname] = (rname, pos, ref_end, strand)

    return alignments, stats


def compute_concordance(aln1, aln2, tolerance=500, metric='position'):
    """Compute concordance between two alignment sets.

    Args:
        aln1: dict of qname -> (rname, pos, ref_end, strand)
        aln2: dict of qname -> (rname, pos, ref_end, strand)
        tolerance: position tolerance in bp (for 'position' metric)
        metric: 'position', 'overlap', or 'both'

    Returns:
        dict with concordance statistics
    """
    keys1 = set(aln1.keys())
    keys2 = set(aln2.keys())

    common = keys1 & keys2
    only1 = keys1 - keys2
    only2 = keys2 - keys1

    concordant_pos = 0
    concordant_ovl = 0
    discordant_chr = 0
    discordant_strand = 0
    discordant_pos = 0
    discordant_ovl = 0

    for qname in common:
        rname1, pos1, end1, strand1 = aln1[qname]
        rname2, pos2, end2, strand2 = aln2[qname]

        if rname1 != rname2:
            discordant_chr += 1
            continue

        if strand1 != strand2:
            discordant_strand += 1
            discordant_pos += 1
            discordant_ovl += 1
            continue

        # Position-based check
        if abs(pos1 - pos2) <= tolerance:
            concordant_pos += 1
        else:
            discordant_pos += 1

        # Overlap-based check
        len1 = end1 - pos1
        len2 = end2 - pos2
        max_len = max(len1, len2)
        if max_len > 0:
            overlap = compute_overlap(pos1, len1, pos2, len2)
            threshold = 0.1 * max_len
            if overlap >= threshold:
                concordant_ovl += 1
            else:
                discordant_ovl += 1
        else:
            # Both have zero-length CIGAR — fall back to position check
            if abs(pos1 - pos2) <= tolerance:
                concordant_ovl += 1
            else:
                discordant_ovl += 1

    total_common = len(common)
    rate_pos = 100.0 * concordant_pos / total_common if total_common > 0 else 0.0
    rate_ovl = 100.0 * concordant_ovl / total_common if total_common > 0 else 0.0

    result = {
        'total1': len(keys1),
        'total2': len(keys2),
        'common': total_common,
        'only1': len(only1),
        'only2': len(only2),
        'discordant_chr': discordant_chr,
        'discordant_strand': discordant_strand,
        'tolerance': tolerance
    }

    if metric == 'position':
        result['concordant'] = concordant_pos
        result['discordant_pos'] = discordant_pos
        result['concordance_rate'] = rate_pos
    elif metric == 'overlap':
        result['concordant'] = concordant_ovl
        result['discordant_pos'] = discordant_ovl
        result['concordance_rate'] = rate_ovl
    else:  # both
        result['concordant_pos'] = concordant_pos
        result['discordant_pos'] = discordant_pos
        result['concordance_rate_pos'] = rate_pos
        result['concordant_ovl'] = concordant_ovl
        result['discordant_ovl'] = discordant_ovl
        result['concordance_rate_overlap'] = rate_ovl
        # For backward compat, set concordance_rate to position-based
        result['concordance_rate'] = rate_pos

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Compute pairwise concordance between aligners',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --sam1 freemap.sam --sam2 minimap2.sam
  %(prog)s --sam1 freemap.sam --sam2 minimap2.sam --tolerance 100
  %(prog)s --sam1 freemap.sam --sam2 minimap2.sam --metric overlap
  %(prog)s --sam1 freemap.sam --sam2 minimap2.sam --metric both
  %(prog)s --sam1 freemap.sam --sam2 minimap2.sam --name1 freemap --name2 minimap2
"""
    )
    parser.add_argument('--sam1', required=True, help='First SAM file')
    parser.add_argument('--sam2', required=True, help='Second SAM file')
    parser.add_argument('--name1', default='aligner1', help='Name for first aligner')
    parser.add_argument('--name2', default='aligner2', help='Name for second aligner')
    parser.add_argument('--tolerance', type=int, default=500,
                        help='Position tolerance in bp (default: 500)')
    parser.add_argument('--metric', choices=['position', 'overlap', 'both'],
                        default='position',
                        help='Concordance metric: position (POS ±tolerance), '
                             'overlap (>=10%% of max interval), or both')
    parser.add_argument('--quiet', action='store_true',
                        help='Only output concordance rate')
    parser.add_argument('--json', action='store_true',
                        help='Output in JSON format')
    args = parser.parse_args()

    # Parse SAM files
    aln1, stats1 = parse_sam(args.sam1)
    aln2, stats2 = parse_sam(args.sam2)

    # Compute concordance
    results = compute_concordance(aln1, aln2, args.tolerance, args.metric)

    if args.json:
        import json
        output = {
            'aligner1': args.name1,
            'aligner2': args.name2,
            'file1': args.sam1,
            'file2': args.sam2,
            'metric': args.metric,
            'stats1': stats1,
            'stats2': stats2,
            'concordance': results
        }
        print(json.dumps(output, indent=2))
    elif args.quiet:
        if args.metric == 'both':
            print(f"pos={results['concordance_rate_pos']:.2f} overlap={results['concordance_rate_overlap']:.2f}")
        else:
            print(f"{results['concordance_rate']:.2f}")
    else:
        print(f"=== Concordance: {args.name1} vs {args.name2} ===")
        print()
        print(f"{args.name1}:")
        print(f"  Mapped reads: {results['total1']:,}")
        if stats1['duplicates_skipped'] > 0:
            print(f"  Duplicates skipped: {stats1['duplicates_skipped']:,}")
        print()
        print(f"{args.name2}:")
        print(f"  Mapped reads: {results['total2']:,}")
        if stats2['duplicates_skipped'] > 0:
            print(f"  Duplicates skipped: {stats2['duplicates_skipped']:,}")
        print()
        print(f"Comparison:")
        print(f"  Common reads (both mapped): {results['common']:,}")
        print(f"  Only {args.name1}: {results['only1']:,}")
        print(f"  Only {args.name2}: {results['only2']:,}")
        print()

        if args.metric == 'both':
            print(f"Concordance (on {results['common']:,} common reads):")
            print(f"  [Position ±{args.tolerance}bp]")
            print(f"    Same chr + pos: {results['concordant_pos']:,} ({results['concordance_rate_pos']:.2f}%)")
            print(f"    Different chr:  {results['discordant_chr']:,}")
            print(f"    Different pos:  {results['discordant_pos']:,}")
            print()
            print(f"  [Overlap >=10% of max interval]")
            print(f"    Same chr + ovl: {results['concordant_ovl']:,} ({results['concordance_rate_overlap']:.2f}%)")
            print(f"    Different chr:  {results['discordant_chr']:,}")
            print(f"    No overlap:     {results['discordant_ovl']:,}")
        elif args.metric == 'overlap':
            print(f"Concordance (on {results['common']:,} common reads, overlap>=10%):")
            print(f"  Same chr + ovl: {results['concordant']:,} ({results['concordance_rate']:.2f}%)")
            print(f"  Different chr:  {results['discordant_chr']:,}")
            print(f"  No overlap:     {results['discordant_pos']:,}")
        else:
            print(f"Concordance (on {results['common']:,} common reads, pos±{args.tolerance}bp):")
            print(f"  Same chr + pos: {results['concordant']:,} ({results['concordance_rate']:.2f}%)")
            print(f"  Different chr:  {results['discordant_chr']:,}")
            print(f"  Different pos:  {results['discordant_pos']:,}")
            if results['discordant_strand'] > 0:
                print(f"    (of which strand-only: {results['discordant_strand']:,})")


if __name__ == '__main__':
    main()
