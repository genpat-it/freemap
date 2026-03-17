#!/usr/bin/env python3
"""
mapq_calibration_multi.py - Multi-panel MAPQ calibration figure

Shows MAPQ calibration across organisms (E.coli, Human) and technologies (ONT, HiFi, CLR)

Usage:
    python3 mapq_calibration_multi.py output.pdf
"""

import sys
import os
import re
import gzip
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

POSITION_TOLERANCE = 100  # bp tolerance for correct mapping

_DATA_DIR = os.environ.get('DATA_DIR', str(Path(__file__).resolve().parents[3] / 'data'))

# Dataset configurations: 2x2 matrix (organism x technology)
# Using chr22 subsets for Human to keep processing time reasonable
DATASETS = {
    # Row 1: E.coli (full datasets, ~20k reads each)
    ('E.coli', 'ONT'): {
        'path': f'{_DATA_DIR}/L1_ecoli_ont_50x',
        'label': 'E.coli ONT (sim)',
    },
    ('E.coli', 'HiFi'): {
        'path': f'{_DATA_DIR}/L2_ecoli_hifi_50x',
        'label': 'E.coli HiFi (sim)',
    },
    # Row 2: Human chr22 (simulated, ~30-40k reads each)
    ('Human', 'ONT'): {
        'path': f'{_DATA_DIR}/H4_human_ont_10x_chr22',
        'label': 'Human chr22 ONT (sim)',
    },
    ('Human', 'HiFi'): {
        'path': f'{_DATA_DIR}/H1_human_hifi_5x_chr22',
        'label': 'Human chr22 HiFi (sim)',
    },
}

ORGANISMS = ['E.coli', 'Human']
TECHNOLOGIES = ['ONT', 'HiFi']

ALIGNERS = {
    'freemap': {'color': 'coral', 'marker': 'o'},
    'minimap2': {'color': 'steelblue', 'marker': 's'},
}


def parse_maf_truth(maf_paths):
    """Parse MAF file(s) to get ground truth positions."""
    truth = {}

    if not isinstance(maf_paths, list):
        maf_paths = [maf_paths]

    for maf_path in maf_paths:
        maf_path = Path(maf_path)

        if str(maf_path).endswith('.gz'):
            opener = lambda p: gzip.open(p, 'rt')
        else:
            opener = lambda p: open(p, 'r')

        with opener(maf_path) as f:
            ref_line = None
            for line in f:
                line = line.strip()
                if line.startswith('a'):
                    ref_line = None
                elif line.startswith('s ref') or line.startswith('s NC_') or line.startswith('s chr'):
                    parts = line.split()
                    ref_line = {
                        'chrom': parts[1],
                        'start': int(parts[2]),
                        'strand': parts[4],
                    }
                elif line.startswith('s S') and ref_line:
                    parts = line.split()
                    read_name = parts[1]
                    read_strand = parts[4]
                    true_strand = '+' if ref_line['strand'] == read_strand else '-'
                    truth[read_name] = {
                        'chrom': ref_line['chrom'],
                        'pos': ref_line['start'],
                        'strand': true_strand,
                    }

    return truth


def parse_sam_mappings(sam_path):
    """Parse SAM file to get mapping positions and MAPQ."""
    mappings = {}

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

            if flag & 4:
                continue

            strand = '-' if flag & 16 else '+'

            mappings[qname] = {
                'chrom': rname,
                'pos': pos,
                'strand': strand,
                'mapq': mapq,
            }

    return mappings


def is_correct_mapping(truth, mapping, tolerance=POSITION_TOLERANCE):
    """Check if mapping is correct within tolerance."""
    if truth['strand'] != mapping['strand']:
        return False
    pos_diff = abs(truth['pos'] - mapping['pos'])
    return pos_diff <= tolerance


def calculate_calibration_data(truth, mappings, bin_size=5):
    """Calculate error rate per MAPQ bin."""
    mapq_correct = defaultdict(int)
    mapq_total = defaultdict(int)

    for read_name, mapping in mappings.items():
        if read_name not in truth:
            continue

        mapq = mapping['mapq']
        is_correct = is_correct_mapping(truth[read_name], mapping)

        mapq_total[mapq] += 1
        if is_correct:
            mapq_correct[mapq] += 1

    # Bin results
    bins = defaultdict(lambda: {'correct': 0, 'total': 0})
    for mapq in mapq_total:
        bin_start = (mapq // bin_size) * bin_size
        bins[bin_start]['correct'] += mapq_correct[mapq]
        bins[bin_start]['total'] += mapq_total[mapq]

    # Calculate stats
    total = sum(mapq_total.values())
    correct = sum(mapq_correct.values())
    accuracy = correct / total * 100 if total > 0 else 0

    return bins, total, correct, accuracy


def plot_calibration(ax, bins, color, marker, label, bin_size=5, min_reads=10):
    """Plot calibration data on axis."""
    mapq_values = []
    error_rates = []
    sizes = []

    for mapq_bin in sorted(bins.keys()):
        data = bins[mapq_bin]
        if data['total'] < min_reads:
            continue

        error_count = data['total'] - data['correct']
        error_rate = error_count / data['total']

        if error_rate == 0:
            error_rate = 0.5 / data['total']

        mapq_values.append(mapq_bin + bin_size / 2)
        error_rates.append(error_rate)
        sizes.append(np.sqrt(data['total']) * 3)

    ax.scatter(mapq_values, error_rates,
               s=sizes, c=color, marker=marker,
               label=label, alpha=0.7,
               edgecolors='black', linewidths=0.5)


def main():
    parser = argparse.ArgumentParser(description='Generate multi-panel MAPQ calibration figure')
    parser.add_argument('output', nargs='?', default='mapq_calibration_multi.pdf', help='Output file')
    args = parser.parse_args()

    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    for row_idx, organism in enumerate(ORGANISMS):
        for col_idx, tech in enumerate(TECHNOLOGIES):
            ax = axes[row_idx, col_idx]
            key = (organism, tech)

            if key not in DATASETS:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue

            dataset_cfg = DATASETS[key]
            dataset_path = Path(dataset_cfg['path'])

            print(f"\n{'='*60}")
            print(f"Processing: {dataset_cfg['label']}")
            print(f"{'='*60}")

            # Find and parse MAF files
            maf_files = list(dataset_path.glob('*.maf')) + list(dataset_path.glob('*.maf.gz'))
            if not maf_files:
                print(f"  ERROR: No MAF files found in {dataset_path}")
                ax.text(0.5, 0.5, 'No MAF', ha='center', va='center', transform=ax.transAxes)
                continue

            print(f"  Parsing {len(maf_files)} MAF file(s)...")
            truth = parse_maf_truth(maf_files)
            print(f"  Found {len(truth)} reads with ground truth")

            # Process each aligner
            for aligner_name, aligner_cfg in ALIGNERS.items():
                sam_path = dataset_path / f'{aligner_name}.sam'
                if not sam_path.exists():
                    print(f"  {aligner_name}: SAM not found")
                    continue

                print(f"  Processing {aligner_name}...")
                mappings = parse_sam_mappings(sam_path)

                bins, total, correct, accuracy = calculate_calibration_data(truth, mappings)

                print(f"    {correct}/{total} correct ({accuracy:.2f}%)")

                plot_calibration(ax, bins, aligner_cfg['color'], aligner_cfg['marker'],
                                f"{aligner_name} ({accuracy:.1f}%)")

            # Plot perfect calibration line
            mapq_line = np.arange(0, 65, 1)
            theoretical_error = 10 ** (-mapq_line / 10)
            ax.plot(mapq_line, theoretical_error, 'k--', linewidth=1.5,
                    label='Perfect calibration', alpha=0.7)

            # Format axis
            ax.set_yscale('log')
            ax.set_xlim(-2, 65)
            ax.set_ylim(1e-4, 1.5)
            ax.grid(True, alpha=0.3, which='both')
            ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
            ax.axhline(y=0.01, color='gray', linestyle=':', alpha=0.5)

            # Show legend on each panel
            ax.legend(loc='upper right', fontsize=7)

            # Panel title
            ax.set_title(dataset_cfg['label'], fontsize=10, fontweight='bold')

            # Labels
            if row_idx == 1:
                ax.set_xlabel('MAPQ', fontsize=10)
            if col_idx == 0:
                ax.set_ylabel('Error rate', fontsize=10)

    # Add panel labels (A-D)
    for idx, ax in enumerate(axes.flat):
        ax.text(-0.12, 1.12, chr(65 + idx), transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top')

    plt.tight_layout()
    plt.savefig(args.output, dpi=200, bbox_inches='tight')
    print(f"\nSaved: {args.output}")


if __name__ == '__main__':
    main()
