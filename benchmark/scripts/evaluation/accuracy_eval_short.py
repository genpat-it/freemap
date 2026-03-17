#!/usr/bin/env python3
"""
accuracy_eval_short.py - Evaluate short-read alignment accuracy against Mason2 ground truth

Mason2 outputs ground truth positions in the read names or in a separate SAM file.
This script compares aligned positions against truth positions.

Two evaluation metrics are supported:
  - position (default): Read is correct if within TOLERANCE bp of truth start
  - overlap: Read is correct if overlap >= 10% of max(truth_len, aln_len)

Only PRIMARY alignments are evaluated. Secondary (0x100) and supplementary
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


def parse_truth_sam(sam_file):
    """Parse Mason2 truth SAM file to extract ground truth positions.

    Returns dict: read_name -> (ref_name, pos, length)
    """
    truth = {}
    with open(sam_file) as f:
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

            qname = normalize_read_name(fields[0])
            rname = fields[2]

            if rname == '*' or rname == '':
                continue

            try:
                pos = int(fields[3])
            except ValueError:
                continue

            cigar = fields[5]
            ref_len = cigar_ref_length(cigar) if cigar != '*' else 0

            # For paired reads, use /1 or /2 suffix or just first occurrence
            if qname not in truth:
                truth[qname] = (rname, pos, ref_len)

    return truth


def parse_alignment_sam(sam_file):
    """Parse alignment SAM file to extract positions.

    Returns dict: read_name -> (ref_name, pos, mapq, ref_length)
    """
    alignments = {}
    with open(sam_file) as f:
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

            qname = normalize_read_name(fields[0])
            rname = fields[2]

            if rname == '*' or rname == '':
                continue

            try:
                pos = int(fields[3])
                mapq = int(fields[4])
            except ValueError:
                continue

            cigar = fields[5]
            ref_len = cigar_ref_length(cigar) if cigar != '*' else 0

            if qname not in alignments:
                alignments[qname] = (rname, pos, mapq, ref_len)

    return alignments


def cigar_ref_length(cigar):
    """Compute reference length consumed by CIGAR string."""
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


def evaluate(truth, alignments, tolerance=500, metric='position'):
    """Evaluate alignment accuracy."""
    total = len(truth)
    mapped = 0
    correct = 0
    wrong_chr = 0
    wrong_pos = 0
    errors_low_mapq = 0

    for read_name, (true_chr, true_pos, true_len) in truth.items():
        if read_name not in alignments:
            continue

        mapped += 1
        aln_chr, aln_pos, mapq, aln_len = alignments[read_name]

        if aln_chr != true_chr:
            wrong_chr += 1
            if mapq < 30:
                errors_low_mapq += 1
        else:
            if metric == 'overlap':
                overlap = compute_overlap(true_pos, true_len, aln_pos, aln_len)
                threshold = 0.1 * max(true_len, aln_len) if max(true_len, aln_len) > 0 else 1
                is_correct = overlap >= threshold
            else:
                is_correct = abs(aln_pos - true_pos) <= tolerance

            if is_correct:
                correct += 1
            else:
                wrong_pos += 1
                if mapq < 30:
                    errors_low_mapq += 1

    unmapped = total - mapped
    total_errors = wrong_chr + wrong_pos

    accuracy = 100.0 * correct / mapped if mapped > 0 else 0.0
    mapping_rate = 100.0 * mapped / total if total > 0 else 0.0

    if total_errors > 0:
        error_detection = 100.0 * errors_low_mapq / total_errors
        error_detection_str = f"{error_detection:.1f}%"
    else:
        error_detection = None
        error_detection_str = "N/A (no errors)"

    return {
        'total': total,
        'mapped': mapped,
        'correct': correct,
        'wrong_chr': wrong_chr,
        'wrong_pos': wrong_pos,
        'unmapped': unmapped,
        'accuracy': accuracy,
        'mapping_rate': mapping_rate,
        'error_detection': error_detection,
        'error_detection_str': error_detection_str,
        'metric': metric
    }


def main():
    parser = argparse.ArgumentParser(description='Evaluate short-read alignment accuracy')
    parser.add_argument('--truth', required=True, help='Mason2 ground truth SAM file')
    parser.add_argument('--sam', required=True, help='Alignment SAM file to evaluate')
    parser.add_argument('--tolerance', type=int, default=500, help='Position tolerance (bp)')
    parser.add_argument('--metric', choices=['position', 'overlap'], default='position',
                        help='Correctness metric')
    parser.add_argument('--json', action='store_true', help='Output JSON format')
    args = parser.parse_args()

    truth = parse_truth_sam(args.truth)
    alignments = parse_alignment_sam(args.sam)

    results = evaluate(truth, alignments, args.tolerance, args.metric)

    if args.json:
        import json
        print(json.dumps(results, indent=2))
    else:
        metric_desc = "overlap>=10%" if args.metric == 'overlap' else f"pos±{args.tolerance}bp"
        print(f"Metric:           {metric_desc}")
        print(f"Total reads:      {results['total']:,}")
        print(f"Mapped:           {results['mapped']:,}")
        print(f"Correct:          {results['correct']:,} ({results['accuracy']:.2f}%)")
        print(f"Wrong chromosome: {results['wrong_chr']:,}")
        print(f"Wrong position:   {results['wrong_pos']:,}")
        print(f"Unmapped:         {results['unmapped']:,}")
        print(f"Mapping rate:     {results['mapping_rate']:.2f}%")
        print(f"Error detection:  {results['error_detection_str']} (MAPQ<30)")


if __name__ == '__main__':
    main()
