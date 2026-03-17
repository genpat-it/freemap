#!/usr/bin/env python3
"""
coverage_scatter_paper.py - Multi-panel coverage scatter for paper (no titles)

5 panels showing freemap vs minimap2 coverage agreement:
  Row 1: E.coli ONT, E.coli HiFi, E.coli CLR (simulated)
  Row 2: Human HiFi (sim), Human ONT (sim)

Plus 2 extra panels for GIAB real data showing the decoy effect.

Usage:
    python3 coverage_scatter_paper.py output.pdf
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path

BASE = Path(os.environ.get('DATA_DIR', str(Path(__file__).resolve().parents[3] / 'data')))

# Simulated datasets
SIMULATED = [
    ("L1_ecoli_ont_50x", "freemap.depth", "minimap2.depth"),
    ("L2_ecoli_hifi_50x", "freemap_ccs.depth", "minimap2.depth"),
    ("L3_ecoli_clr_50x", "freemap.depth", "minimap2.depth"),
    ("H1_human_hifi_5x", "freemap_chr22.depth", "minimap2_chr22.depth"),
    ("H4_human_ont_10x", "freemap_chr22.depth", "minimap2_chr22.depth"),
]

# Real GIAB datasets
GIAB = [
    ("D1_giab_hg002_hifi", "freemap_chr1.depth", "minimap2_chr1.depth"),
    ("D2_giab_hg002_ont", "freemap_chr1.depth", "minimap2_chr1.depth"),
]

def parse_depth_binned(filename, bin_size=1000):
    """Parse depth file into binned coverage."""
    sum_cov = defaultdict(float)
    count = defaultdict(int)

    with open(filename) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t', 3)
            if len(parts) >= 3:
                chrom = parts[0]
                pos = int(parts[1])
                depth = int(parts[2])
                bin_key = (chrom, (pos - 1) // bin_size)
                sum_cov[bin_key] += depth
                count[bin_key] += 1

    return {k: sum_cov[k] / count[k] for k in sum_cov}


def make_panel(ax, d1_bin, d2_bin, panel_label, xlabel, ylabel, rng):
    """Create a single scatter panel without title."""
    common_bin = sorted(set(d1_bin) & set(d2_bin))

    if len(common_bin) < 2:
        ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
        return float('nan')

    x_bin = np.array([d1_bin[k] for k in common_bin])
    y_bin = np.array([d2_bin[k] for k in common_bin])

    if np.std(x_bin) == 0 or np.std(y_bin) == 0:
        r_bin = float('nan')
    else:
        r_bin = np.corrcoef(x_bin, y_bin)[0, 1]

    # Subsample for visibility (deterministic)
    n_points = min(10000, len(x_bin))
    idx = rng.choice(len(x_bin), size=n_points, replace=False)

    ax.scatter(x_bin[idx], y_bin[idx], alpha=0.4, s=3, c='steelblue')

    # Axis limits based on 99th percentile
    max_val = max(np.percentile(x_bin, 99), np.percentile(y_bin, 99))
    ax.plot([0, max_val * 1.1], [0, max_val * 1.1], 'r--', lw=1, alpha=0.7)
    ax.set_xlim(0, max_val * 1.1)
    ax.set_ylim(0, max_val * 1.1)
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)

    # Panel label in top-left
    ax.text(0.05, 0.95, panel_label, transform=ax.transAxes,
            ha='left', va='top', fontsize=11, fontweight='bold')

    # Correlation annotation in bottom-right
    ax.text(0.95, 0.05, f'r = {r_bin:.2f}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    return r_bin


def main():
    output_file = sys.argv[1] if len(sys.argv) > 1 else "fig_coverage_scatter.pdf"
    rng = np.random.default_rng(42)

    # Create figure with 3 panels on top, 2 centered on bottom
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(11, 7))
    gs_top = GridSpec(1, 3, figure=fig, top=0.95, bottom=0.55,
                      left=0.06, right=0.98, wspace=0.35)
    gs_bot = GridSpec(1, 2, figure=fig, top=0.45, bottom=0.05,
                      left=0.17, right=0.87, wspace=0.35)

    axes_top = [fig.add_subplot(gs_top[0, i]) for i in range(3)]
    axes_bot = [fig.add_subplot(gs_bot[0, i]) for i in range(2)]
    all_axes = axes_top + axes_bot

    panel_labels = ['A', 'B', 'C', 'D', 'E']
    panel_names = [
        'E. coli ONT 50x',
        'E. coli HiFi 50x',
        'E. coli CLR 50x',
        'Human HiFi 5x',
        'Human ONT 10x',
    ]

    # First 5 panels: simulated data
    for i, (dataset, fm_file, mm2_file) in enumerate(SIMULATED):
        ax = all_axes[i]

        fm_path = BASE / dataset / fm_file
        mm2_path = BASE / dataset / mm2_file

        print(f"Loading {dataset}...")
        try:
            fm_bin = parse_depth_binned(fm_path)
            mm2_bin = parse_depth_binned(mm2_path)
            r = make_panel(ax, fm_bin, mm2_bin, panel_labels[i],
                          "freemap depth", "minimap2 depth", rng)
            print(f"  {panel_labels[i]}: {panel_names[i]}, r = {r:.3f}")
        except Exception as e:
            print(f"  Error: {e}")
            ax.text(0.5, 0.5, f"Error", ha='center', va='center', transform=ax.transAxes)

    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"\nSaved: {output_file}")

    # Now create a separate 1x2 figure for GIAB real data
    output_giab = output_file.replace('.pdf', '_giab.pdf')
    fig2, axes2 = plt.subplots(1, 2, figsize=(8, 4))

    giab_labels = ['G', 'H']
    giab_names = ['GIAB HiFi (real)', 'GIAB ONT (real)']

    for i, (dataset, fm_file, mm2_file) in enumerate(GIAB):
        ax = axes2[i]
        fm_path = BASE / dataset / fm_file
        mm2_path = BASE / dataset / mm2_file

        print(f"Loading {dataset}...")
        try:
            fm_bin = parse_depth_binned(fm_path)
            mm2_bin = parse_depth_binned(mm2_path)
            r = make_panel(ax, fm_bin, mm2_bin, giab_labels[i],
                          "freemap depth", "minimap2 depth", rng)
            print(f"  {giab_labels[i]}: {giab_names[i]}, r = {r:.3f}")
        except Exception as e:
            print(f"  Error: {e}")

    plt.tight_layout()
    plt.savefig(output_giab, dpi=200, bbox_inches='tight')
    print(f"Saved: {output_giab}")


if __name__ == '__main__':
    main()
