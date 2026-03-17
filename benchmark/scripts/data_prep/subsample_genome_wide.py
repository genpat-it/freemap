#!/usr/bin/env python3
"""
Subsample reads genome-wide from per-chromosome files.

Usage:
    python3 subsample_genome_wide.py --dataset H5 --output reads_100k.fq --target 100000
    python3 subsample_genome_wide.py --dataset H4 --output reads_100k.fq --target 100000
    python3 subsample_genome_wide.py --dataset H1 --output reads_100k.fq --target 100000
"""

import argparse
import gzip
import subprocess
import sys
import os
from pathlib import Path


def iter_fastq(fileobj):
    """Yield FASTQ records (4-line tuples) from a file object."""
    while True:
        header = fileobj.readline()
        if not header:
            break
        seq = fileobj.readline()
        plus = fileobj.readline()
        qual = fileobj.readline()
        if not qual:
            break
        yield (header, seq, plus, qual)


def subsample_plain_fastq(input_path, output_path, target, total_reads):
    """Subsample from a plain FASTQ file with known read count."""
    stride = max(1, total_reads // target)
    count = 0
    written = 0
    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        for record in iter_fastq(fin):
            if count % stride == 0 and written < target:
                fout.writelines(record)
                written += 1
            count += 1
    return written


def subsample_fqgz_files(data_dir, output_path, target):
    """Subsample from multiple .fq.gz files."""
    fqgz_files = sorted(data_dir.glob('reads_*.fq.gz'))
    if not fqgz_files:
        print(f"ERROR: No reads_*.fq.gz files in {data_dir}", file=sys.stderr)
        sys.exit(1)

    # Count total reads across all files
    print(f"Counting reads across {len(fqgz_files)} files...")
    total = 0
    file_counts = []
    for f in fqgz_files:
        n = 0
        with gzip.open(f, 'rt') as gz:
            for line in gz:
                if line.startswith('@'):
                    n += 1
                    # Skip 3 lines (seq, +, qual) - but in counting mode
                    # just count @ lines at position 0 mod 4
        # Actually, simpler: count lines / 4
        # Re-count properly
        pass

    # Simpler: count lines
    print(f"Counting total lines (this may take a few minutes)...")
    total_lines = 0
    for f in fqgz_files:
        with gzip.open(f, 'rt') as gz:
            for _ in gz:
                total_lines += 1
    total = total_lines // 4
    print(f"Total reads: {total:,}")

    stride = max(1, total // target)
    print(f"Using stride: {stride} (expect ~{total // stride:,} reads)")

    # Second pass: write every Nth read
    count = 0
    written = 0
    with open(output_path, 'w') as fout:
        for f in fqgz_files:
            with gzip.open(f, 'rt') as gz:
                for record in iter_fastq(gz):
                    if count % stride == 0 and written < target:
                        fout.writelines(record)
                        written += 1
                    count += 1
            if written >= target:
                break
    return written


def subsample_bam_files(data_dir, output_path, target):
    """Subsample from multiple .bam files using samtools."""
    bam_files = sorted(data_dir.glob('reads_*.bam'))
    if not bam_files:
        print(f"ERROR: No reads_*.bam files in {data_dir}", file=sys.stderr)
        sys.exit(1)

    # Count total reads across all BAM files
    print(f"Counting reads across {len(bam_files)} BAM files...")
    total = 0
    for f in bam_files:
        result = subprocess.run(
            ['samtools', 'view', '-c', str(f)],
            capture_output=True, text=True
        )
        n = int(result.stdout.strip())
        total += n
    print(f"Total reads: {total:,}")

    stride = max(1, total // target)
    print(f"Using stride: {stride} (expect ~{total // stride:,} reads)")

    # Second pass: extract every Nth read
    count = 0
    written = 0
    with open(output_path, 'w') as fout:
        for f in bam_files:
            proc = subprocess.Popen(
                ['samtools', 'fastq', str(f)],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                text=True
            )
            for record in iter_fastq(proc.stdout):
                if count % stride == 0 and written < target:
                    fout.writelines(record)
                    written += 1
                count += 1
            proc.wait()
            if written >= target:
                break
    return written


def main():
    parser = argparse.ArgumentParser(description='Subsample reads genome-wide')
    parser.add_argument('--dataset', required=True, choices=['H1', 'H4', 'H5'],
                        help='Dataset to subsample')
    parser.add_argument('--output', required=True, help='Output FASTQ path')
    parser.add_argument('--target', type=int, default=100000,
                        help='Target number of reads (default: 100000)')
    parser.add_argument('--data-dir', default=None,
                        help='Override data directory')
    args = parser.parse_args()

    data_root = Path(os.environ.get('DATA_DIR',
                                     str(Path(__file__).resolve().parents[3] / 'data')))

    dataset_dirs = {
        'H1': 'H1_human_hifi_5x',
        'H4': 'H4_human_ont_10x',
        'H5': 'H5_human_clr_5x',
    }

    if args.data_dir:
        data_dir = Path(args.data_dir)
    else:
        data_dir = data_root / dataset_dirs[args.dataset]

    if not data_dir.exists():
        # Try local data/ symlink
        data_dir = Path(__file__).resolve().parent.parent.parent.parent / 'data' / dataset_dirs[args.dataset]

    print(f"Dataset: {args.dataset}")
    print(f"Data dir: {data_dir}")
    print(f"Target: {args.target:,} reads")
    print(f"Output: {args.output}")
    print()

    if args.dataset == 'H5':
        # H5 CLR has reads_combined.fq
        combined = data_dir / 'reads_combined.fq'
        if combined.exists():
            # Count reads (already known: ~1.56M)
            print(f"Using combined file: {combined}")
            total_lines = sum(1 for _ in open(combined))
            total = total_lines // 4
            print(f"Total reads: {total:,}")
            written = subsample_plain_fastq(combined, args.output, args.target, total)
        else:
            written = subsample_fqgz_files(data_dir, args.output, args.target)
    elif args.dataset == 'H4':
        written = subsample_fqgz_files(data_dir, args.output, args.target)
    elif args.dataset == 'H1':
        written = subsample_bam_files(data_dir, args.output, args.target)

    print(f"\nWrote {written:,} reads to {args.output}")

    # Verify
    verify_lines = sum(1 for _ in open(args.output))
    verify_reads = verify_lines // 4
    print(f"Verification: {verify_reads:,} reads ({verify_lines:,} lines)")

    # Show chromosome distribution
    from collections import Counter
    chroms = Counter()
    with open(args.output) as f:
        for line in f:
            if line.startswith('@'):
                # Read name format: @S<chr>_<id> or @S<chr>/<id>
                name = line.strip().lstrip('@')
                # Extract chromosome number
                if name.startswith('S'):
                    parts = name[1:].split('_', 1)
                    if not parts[0]:
                        parts = name[1:].split('/', 1)
                    try:
                        chrom = int(parts[0])
                        chroms[chrom] += 1
                    except ValueError:
                        chroms['other'] += 1
    print(f"\nChromosome distribution ({len(chroms)} chromosomes):")
    for chrom in sorted(chroms.keys(), key=lambda x: (isinstance(x, str), x)):
        print(f"  S{chrom}: {chroms[chrom]:,} reads")


if __name__ == '__main__':
    main()
