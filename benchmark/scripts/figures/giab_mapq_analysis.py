#!/usr/bin/env python3
"""
giab_mapq_analysis.py - MAPQ analysis using GIAB high-confidence regions

For real data, we can validate MAPQ by checking if discordant reads
fall inside or outside GIAB high-confidence regions:
- Inside HC regions: alignment should be reliable
- Outside HC regions: alignment uncertainty is expected
"""

import sys
import os
import numpy as np
from pathlib import Path
from collections import defaultdict
from intervaltree import IntervalTree

# NCBI accession to chromosome name mapping (GRCh38)
NCBI_TO_CHR = {
    'CM000663.2': 'chr1',
    'CM000664.2': 'chr2',
    'CM000665.2': 'chr3',
    'CM000666.2': 'chr4',
    'CM000667.2': 'chr5',
    'CM000668.2': 'chr6',
    'CM000669.2': 'chr7',
    'CM000670.2': 'chr8',
    'CM000671.2': 'chr9',
    'CM000672.2': 'chr10',
    'CM000673.2': 'chr11',
    'CM000674.2': 'chr12',
    'CM000675.2': 'chr13',
    'CM000676.2': 'chr14',
    'CM000677.2': 'chr15',
    'CM000678.2': 'chr16',
    'CM000679.2': 'chr17',
    'CM000680.2': 'chr18',
    'CM000681.2': 'chr19',
    'CM000682.2': 'chr20',
    'CM000683.2': 'chr21',
    'CM000684.2': 'chr22',
    'CM000685.2': 'chrX',
    'CM000686.2': 'chrY',
}

POSITION_TOLERANCE = 500

_DATA_DIR = os.environ.get('DATA_DIR', str(Path(__file__).resolve().parents[3] / 'data'))

DATASETS = {
    'D1 HiFi': f'{_DATA_DIR}/D1_giab_hg002_hifi',
    'D2 ONT': f'{_DATA_DIR}/D2_giab_hg002_ont',
}

GIAB_BED = f'{_DATA_DIR}/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'

ALIGNERS = ['freemap', 'minimap2']


def load_high_confidence_regions(bed_path):
    """Load GIAB high-confidence regions into interval trees."""
    print(f"Loading high-confidence regions from {bed_path}...")
    regions = defaultdict(IntervalTree)
    with open(bed_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            regions[chrom].addi(start, end)

    total = sum(len(tree) for tree in regions.values())
    print(f"  Loaded {total:,} high-confidence regions")
    return regions


def is_in_high_confidence(chrom_ncbi, pos, hc_regions):
    """Check if position is in high-confidence region."""
    chrom = NCBI_TO_CHR.get(chrom_ncbi)
    if not chrom:
        return False  # Unknown chromosome
    if chrom not in hc_regions:
        return False
    return len(hc_regions[chrom].at(pos)) > 0


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
    return mappings


def is_concordant(m1, m2, tolerance=POSITION_TOLERANCE):
    """Check if two mappings are concordant."""
    if m1['chrom'] != m2['chrom']:
        return False
    if m1['strand'] != m2['strand']:
        return False
    return abs(m1['pos'] - m2['pos']) <= tolerance


def main():
    print("=" * 80)
    print("MAPQ Analysis with GIAB High-Confidence Regions")
    print("=" * 80)

    # Load high-confidence regions
    hc_regions = load_high_confidence_regions(GIAB_BED)

    results = {}

    for dataset_name, dataset_path in DATASETS.items():
        print(f"\n{'='*60}")
        print(f"Dataset: {dataset_name}")
        print(f"{'='*60}")

        # Load all aligner mappings
        mappings = {}
        for aligner in ALIGNERS:
            sam_path = Path(dataset_path) / f'{aligner}.sam'
            if not sam_path.exists():
                continue
            print(f"  Loading {aligner}...")
            mappings[aligner] = parse_sam_mappings(sam_path)
            print(f"    {len(mappings[aligner]):,} primary mappings")

        if 'minimap2' not in mappings:
            continue

        # Analyze freemap vs minimap2 (both perspectives)
        if 'freemap' not in mappings:
            continue

        common_reads = set(mappings['freemap'].keys()) & set(mappings['minimap2'].keys())

        # Categorize reads by both aligners' perspectives
        fm_disc_in_hc = []
        fm_disc_out_hc = []
        fm_conc_in_hc = []
        fm_conc_out_hc = []
        mm2_disc_in_hc = []
        mm2_disc_out_hc = []
        mm2_conc_in_hc = []
        mm2_conc_out_hc = []

        # Directional analysis: discordant reads where one maps in HC, other outside
        fm_in_hc_mm2_out = 0
        mm2_in_hc_fm_out = 0

        for read in common_reads:
            m_fm = mappings['freemap'][read]
            m_mm2 = mappings['minimap2'][read]

            fm_in_hc = is_in_high_confidence(m_fm['chrom'], m_fm['pos'], hc_regions)
            mm2_in_hc = is_in_high_confidence(m_mm2['chrom'], m_mm2['pos'], hc_regions)
            concordant = is_concordant(m_fm, m_mm2)

            if concordant:
                if fm_in_hc:
                    fm_conc_in_hc.append(m_fm['mapq'])
                    mm2_conc_in_hc.append(m_mm2['mapq'])
                else:
                    fm_conc_out_hc.append(m_fm['mapq'])
                    mm2_conc_out_hc.append(m_mm2['mapq'])
            else:
                # freemap perspective (use freemap's position for HC check)
                if fm_in_hc:
                    fm_disc_in_hc.append(m_fm['mapq'])
                else:
                    fm_disc_out_hc.append(m_fm['mapq'])
                # minimap2 perspective (use minimap2's position for HC check)
                if mm2_in_hc:
                    mm2_disc_in_hc.append(m_mm2['mapq'])
                else:
                    mm2_disc_out_hc.append(m_mm2['mapq'])
                # Directional: one in HC, other outside
                if fm_in_hc and not mm2_in_hc:
                    fm_in_hc_mm2_out += 1
                elif mm2_in_hc and not fm_in_hc:
                    mm2_in_hc_fm_out += 1

        def mapq_stats(mapq_list, name):
            if not mapq_list:
                return None
            n = len(mapq_list)
            mean = np.mean(mapq_list)
            pct_low = sum(1 for m in mapq_list if m < 30) / n * 100
            print(f"    {name}: n={n:,}, mean MAPQ={mean:.1f}, MAPQ<30={pct_low:.1f}%")
            return {'n': n, 'mean': mean, 'pct_low': pct_low}

        print(f"\n  freemap vs minimap2:")
        print(f"  -- freemap MAPQ --")
        r_fm = {}
        r_fm['conc_in_hc'] = mapq_stats(fm_conc_in_hc, "Concordant in HC")
        r_fm['conc_out_hc'] = mapq_stats(fm_conc_out_hc, "Concordant outside HC")
        r_fm['disc_in_hc'] = mapq_stats(fm_disc_in_hc, "Discordant in HC")
        r_fm['disc_out_hc'] = mapq_stats(fm_disc_out_hc, "Discordant outside HC")
        results[(dataset_name, 'freemap')] = r_fm

        print(f"  -- minimap2 MAPQ --")
        r_mm2 = {}
        r_mm2['conc_in_hc'] = mapq_stats(mm2_conc_in_hc, "Concordant in HC")
        r_mm2['conc_out_hc'] = mapq_stats(mm2_conc_out_hc, "Concordant outside HC")
        r_mm2['disc_in_hc'] = mapq_stats(mm2_disc_in_hc, "Discordant in HC")
        r_mm2['disc_out_hc'] = mapq_stats(mm2_disc_out_hc, "Discordant outside HC")
        results[(dataset_name, 'minimap2')] = r_mm2

        # Directional analysis
        decidable = fm_in_hc_mm2_out + mm2_in_hc_fm_out
        print(f"\n  Directional analysis (discordant, one in HC / one outside):")
        print(f"    Decidable: {decidable:,}")
        if decidable > 0:
            print(f"    freemap in HC:  {fm_in_hc_mm2_out:,} ({fm_in_hc_mm2_out/decidable*100:.1f}%)")
            print(f"    minimap2 in HC: {mm2_in_hc_fm_out:,} ({mm2_in_hc_fm_out/decidable*100:.1f}%)")
        results[(dataset_name, 'directional')] = {
            'decidable': decidable,
            'fm_in_hc': fm_in_hc_mm2_out,
            'mm2_in_hc': mm2_in_hc_fm_out,
        }

    # Summary
    print("\n" + "=" * 80)
    print("KEY INSIGHT")
    print("=" * 80)
    print("""
For well-calibrated MAPQ:
1. Discordant reads INSIDE high-confidence regions should have LOW MAPQ
   (We expect reliable alignment in these regions, so disagreement = uncertainty)
2. Discordant reads OUTSIDE high-confidence regions may have any MAPQ
   (Alignment uncertainty is expected in difficult regions)

A good aligner assigns LOW MAPQ to discordant reads in HC regions.
""")

    # Print summary table
    print("\nSUMMARY: Discordant reads in High-Confidence regions")
    print("-" * 70)
    print(f"{'Dataset':<12} {'Aligner':<12} {'N Disc in HC':<15} {'MAPQ<30':<12} {'Mean MAPQ':<12}")
    print("-" * 70)
    for (dataset, aligner), r in results.items():
        if r.get('disc_in_hc'):
            d = r['disc_in_hc']
            print(f"{dataset:<12} {aligner:<12} {d['n']:<15,} {d['pct_low']:.1f}%        {d['mean']:.1f}")


if __name__ == '__main__':
    main()
