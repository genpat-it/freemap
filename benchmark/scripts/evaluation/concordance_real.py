#!/usr/bin/env python3
"""
concordance_real.py - Concordance analysis for real data (no ground truth)

Compares mapping positions between aligners to measure agreement.

Usage:
    python3 concordance_real.py output.pdf
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

POSITION_TOLERANCE = 100  # bp tolerance for concordant mapping

_DATA_DIR = os.environ.get('DATA_DIR', str(Path(__file__).resolve().parents[3] / 'data'))

# Real dataset configurations
DATASETS = {
    'D1 HiFi': {
        'path': f'{_DATA_DIR}/D1_giab_hg002_hifi',
        'label': 'GIAB HG002 HiFi',
    },
    'D2 ONT': {
        'path': f'{_DATA_DIR}/D2_giab_hg002_ont',
        'label': 'GIAB HG002 ONT ultralong',
    },
}

ALIGNERS = ['freemap', 'minimap2']
ALIGNER_COLORS = {
    'freemap': 'coral',
    'minimap2': 'steelblue',
}


def parse_sam_mappings(sam_path, max_reads=None):
    """Parse SAM file to get mapping positions."""
    mappings = {}
    count = 0

    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            try:
                qname = fields[0]
                flag = int(fields[1])
                rname = fields[2]
                pos = int(fields[3])
                mapq = int(fields[4])
            except (ValueError, IndexError):
                continue

            # Skip unmapped and secondary/supplementary
            if flag & 4 or flag & 256 or flag & 2048:
                continue

            strand = '-' if flag & 16 else '+'

            mappings[qname] = {
                'chrom': rname,
                'pos': pos,
                'strand': strand,
                'mapq': mapq,
            }

            count += 1
            if max_reads and count >= max_reads:
                break

    return mappings


def calculate_concordance(ref_mappings, query_mappings, tolerance=POSITION_TOLERANCE):
    """Calculate concordance between two aligners."""
    concordant = 0
    discordant = 0
    ref_only = 0
    query_only = 0

    # Reads in both
    common_reads = set(ref_mappings.keys()) & set(query_mappings.keys())

    for read_name in common_reads:
        ref = ref_mappings[read_name]
        query = query_mappings[read_name]

        # Check concordance
        if ref['chrom'] == query['chrom'] and ref['strand'] == query['strand']:
            if abs(ref['pos'] - query['pos']) <= tolerance:
                concordant += 1
            else:
                discordant += 1
        else:
            discordant += 1

    ref_only = len(set(ref_mappings.keys()) - set(query_mappings.keys()))
    query_only = len(set(query_mappings.keys()) - set(ref_mappings.keys()))

    total = concordant + discordant
    concordance_rate = concordant / total * 100 if total > 0 else 0

    return {
        'concordant': concordant,
        'discordant': discordant,
        'ref_only': ref_only,
        'query_only': query_only,
        'total': total,
        'rate': concordance_rate,
    }


def analyze_discordant_mapq(ref_mappings, query_mappings, tolerance=POSITION_TOLERANCE):
    """Analyze MAPQ distribution for discordant reads."""
    ref_mapq_discord = []
    query_mapq_discord = []

    common_reads = set(ref_mappings.keys()) & set(query_mappings.keys())

    for read_name in common_reads:
        ref = ref_mappings[read_name]
        query = query_mappings[read_name]

        # Check if discordant
        is_concordant = (ref['chrom'] == query['chrom'] and
                        ref['strand'] == query['strand'] and
                        abs(ref['pos'] - query['pos']) <= tolerance)

        if not is_concordant:
            ref_mapq_discord.append(ref['mapq'])
            query_mapq_discord.append(query['mapq'])

    return ref_mapq_discord, query_mapq_discord


def main():
    parser = argparse.ArgumentParser(description='Generate concordance analysis figure')
    parser.add_argument('output', nargs='?', default='concordance_real.pdf', help='Output file')
    parser.add_argument('--max-reads', type=int, default=None, help='Max reads to process (for testing)')
    args = parser.parse_args()

    # Create figure with 2 rows: bar chart + MAPQ analysis
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    results_summary = {}

    for idx, (dataset_name, dataset_cfg) in enumerate(DATASETS.items()):
        dataset_path = Path(dataset_cfg['path'])

        print(f"\n{'='*60}")
        print(f"Processing: {dataset_cfg['label']}")
        print(f"{'='*60}")

        # Load all aligner mappings
        mappings = {}
        for aligner in ALIGNERS:
            sam_path = dataset_path / f'{aligner}.sam'
            if not sam_path.exists():
                print(f"  {aligner}: SAM not found")
                continue

            print(f"  Loading {aligner}...")
            mappings[aligner] = parse_sam_mappings(sam_path, args.max_reads)
            print(f"    {len(mappings[aligner])} primary mappings")

        if 'minimap2' not in mappings:
            print("  ERROR: minimap2 required as reference")
            continue

        # Calculate concordance vs minimap2 (reference)
        concordance_results = {}
        for aligner in ['freemap']:
            if aligner not in mappings:
                continue

            result = calculate_concordance(mappings['minimap2'], mappings[aligner])
            concordance_results[aligner] = result
            print(f"  {aligner} vs minimap2: {result['rate']:.2f}% concordant ({result['concordant']}/{result['total']})")

        results_summary[dataset_name] = concordance_results

        # Panel A/C: Concordance bar chart
        ax = axes[0, idx]
        aligners_plot = [a for a in ['freemap'] if a in concordance_results]
        rates = [concordance_results[a]['rate'] for a in aligners_plot]
        colors = [ALIGNER_COLORS[a] for a in aligners_plot]

        bars = ax.bar(aligners_plot, rates, color=colors, edgecolor='black', linewidth=1)

        # Add value labels
        for bar, rate in zip(bars, rates):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   f'{rate:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

        ax.set_ylabel('Concordance with minimap2 (%)', fontsize=10)
        ax.set_title(f'{dataset_cfg["label"]}', fontsize=11, fontweight='bold')
        ax.set_ylim(0, 105)
        ax.axhline(y=100, color='gray', linestyle='--', alpha=0.5)

        # Panel B/D: MAPQ distribution for discordant reads
        ax = axes[1, idx]

        for aligner in ['freemap']:
            if aligner not in mappings:
                continue

            ref_mapq, query_mapq = analyze_discordant_mapq(mappings['minimap2'], mappings[aligner])

            if query_mapq:
                # Histogram of query aligner MAPQ for discordant reads
                bins = np.arange(0, 65, 5)
                ax.hist(query_mapq, bins=bins, alpha=0.5, label=f'{aligner} (n={len(query_mapq)})',
                       color=ALIGNER_COLORS[aligner], edgecolor='black', linewidth=0.5)

        ax.set_xlabel('MAPQ of discordant reads', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(f'MAPQ distribution of discordant mappings', fontsize=10)
        ax.legend(loc='upper right', fontsize=9)
        ax.set_xlim(0, 65)

    # Add panel labels
    for idx, ax in enumerate(axes.flat):
        ax.text(-0.1, 1.05, chr(65 + idx), transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top')

    plt.suptitle('Concordance Analysis: Real Data (vs minimap2 reference)', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(args.output, dpi=200, bbox_inches='tight')
    print(f"\nSaved: {args.output}")

    # Print summary table
    print("\n" + "=" * 70)
    print("CONCORDANCE SUMMARY (vs minimap2)")
    print("=" * 70)
    print(f"{'Dataset':<25} {'Aligner':<12} {'Concordant':<12} {'Total':<12} {'Rate':<10}")
    print("-" * 70)
    for dataset_name, results in results_summary.items():
        for aligner, r in results.items():
            print(f"{dataset_name:<25} {aligner:<12} {r['concordant']:<12} {r['total']:<12} {r['rate']:.2f}%")


if __name__ == '__main__':
    main()
