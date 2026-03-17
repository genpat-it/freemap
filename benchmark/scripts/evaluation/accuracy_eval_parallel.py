#!/usr/bin/env python3
"""
accuracy_eval_parallel.py - Multithreaded accuracy evaluation against MAF ground truth

Uses multiprocessing.Pool with fork() to share the SAM alignment dict (copy-on-write).
Each worker processes one MAF file independently. ~30x faster than single-threaded on 32 cores.

Usage:
    python3 accuracy_eval_parallel.py --maf-dir data/H1_human_hifi_5x \
        --sam data/H1_human_hifi_5x/freemap.sam \
        --fai data/references/GRCh38.fa.fai --metric both --json -t 32
"""

import argparse
import gzip
import glob
import os
import subprocess
import sys
from multiprocessing import Pool
from typing import Dict, Tuple, Optional

# Global alignment dict — shared via fork() copy-on-write
_ALIGNMENTS = {}
_TOLERANCE = 500
_METRIC = 'position'


def _zmw_id(name: str) -> Optional[str]:
    parts = name.rsplit('/', 1)
    if len(parts) == 2 and (parts[1] == 'ccs' or parts[1].isdigit()):
        return parts[0]
    return None


def cigar_ref_length(cigar: str) -> int:
    num = 0
    ref_len = 0
    for c in cigar:
        if c.isdigit():
            num = num * 10 + int(c)
        else:
            if c in ('M', 'D', 'N', '=', 'X'):
                ref_len += num
            num = 0
    return ref_len


def parse_sam(sam_file: str) -> Dict[str, Tuple[str, int, int, int]]:
    alignments = {}
    fh = gzip.open(sam_file, 'rt') if sam_file.endswith('.gz') else open(sam_file, 'r')
    with fh:
        for line in fh:
            if line[0] == '@' or line[0] == '[':
                continue
            fields = line.split('\t', 12)
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 0x904:
                continue
            read_name = fields[0]
            cigar = fields[5]
            ref_len = cigar_ref_length(cigar) if cigar != '*' else 0
            val = (fields[2], int(fields[3]), int(fields[4]), ref_len)
            alignments[read_name] = val
            zmw = _zmw_id(read_name)
            if zmw and zmw not in alignments:
                alignments[zmw] = val
    return alignments


def find_ref_name_for_maf(maf_path: str) -> Optional[str]:
    base = maf_path
    if base.endswith('.gz'):
        base = base[:-3]
    if base.endswith('.maf'):
        base = base[:-4]
    ref_path = base + '.ref'
    if not os.path.exists(ref_path):
        return None
    with open(ref_path, 'r') as f:
        first = f.readline().strip()
    if not first.startswith('>'):
        return None
    return first[1:].split()[0]


def find_aligned_ref_interval(ref_text: str, read_text: str, ref_start: int):
    min_len = min(len(ref_text), len(read_text))
    first_pos = None
    last_pos = None
    ref_offset = 0
    for i in range(min_len):
        r, q = ref_text[i], read_text[i]
        if r != '-' and q != '-':
            current_pos = ref_start + ref_offset
            if first_pos is None:
                first_pos = current_pos
            last_pos = current_pos
        if r != '-':
            ref_offset += 1
    if first_pos is None:
        return ref_start + 1, 1
    return first_pos + 1, last_pos - first_pos + 1


def _open_maf(path):
    """Open MAF file, using pigz for faster gzip decompression if available."""
    if path.endswith('.gz'):
        try:
            proc = subprocess.Popen(['pigz', '-dc', path], stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL)
            return proc.stdout, proc
        except FileNotFoundError:
            return gzip.open(path, 'rt'), None
    return open(path, 'r'), None


def process_one_maf(args):
    """Worker: process a single MAF file. Returns (counters_tuple, n_zmw_deduped)."""
    maf_path, ref_name_override, ref_names_set = args

    total = mapped = correct = correct_overlap = 0
    wrong_chr = wrong_pos = wrong_pos_overlap = errors_low_mapq = 0
    seen_zmws = set()
    tolerance = _TOLERANCE
    metric = _METRIC
    alignments = _ALIGNMENTS

    ref_names_in_maf = ref_names_set.copy() if ref_names_set else {'ref'}
    if ref_name_override:
        ref_names_in_maf.add(ref_name_override)

    fh, proc = _open_maf(maf_path)
    try:
        current_s = []
        for raw_line in fh:
            if isinstance(raw_line, bytes):
                line = raw_line.decode('utf-8', errors='replace').strip()
            else:
                line = raw_line.strip()
            if not line or line[0] == '#':
                continue
            if line[0] == 'a':
                if len(current_s) >= 2:
                    # Process block
                    ref_line = read_line = None
                    for s in current_s:
                        name = s[1]
                        if name in ref_names_in_maf and ref_line is None:
                            ref_line = s
                        elif name not in ref_names_in_maf and read_line is None:
                            read_line = s
                    if ref_line and read_line:
                        ref_start = int(ref_line[2])
                        ref_size = int(ref_line[3])
                        ref_strand = ref_line[4]
                        ref_src_size = int(ref_line[5])
                        ref_text = ref_line[6] if len(ref_line) > 6 else ""
                        read_name = read_line[1]
                        read_text = read_line[6] if len(read_line) > 6 else ""

                        if ref_strand == '-':
                            ref_start = ref_src_size - ref_start - ref_size

                        if ref_text and read_text:
                            true_pos, true_len = find_aligned_ref_interval(ref_text, read_text, ref_start)
                        else:
                            true_pos, true_len = ref_start + 1, ref_size

                        true_chr = ref_name_override if (ref_name_override and ref_line[1] == 'ref') else ref_line[1]

                        # Lookup alignment
                        aln = alignments.get(read_name)
                        if aln is None:
                            zmw = _zmw_id(read_name)
                            if zmw:
                                if zmw in seen_zmws:
                                    current_s = []
                                    continue
                                seen_zmws.add(zmw)
                                aln = alignments.get(zmw)

                        # Evaluate
                        total += 1
                        if aln is not None:
                            mapped += 1
                            aln_chr, aln_pos, mapq, aln_len = aln
                            if aln_chr != true_chr:
                                wrong_chr += 1
                                if mapq < 30:
                                    errors_low_mapq += 1
                            else:
                                ov = max(0, min(aln_pos + aln_len, true_pos + true_len) - max(aln_pos, true_pos))
                                is_pos = abs(aln_pos - true_pos) <= tolerance
                                is_ovl = ov >= 0.1 * max(true_len, aln_len)
                                if metric in ('position', 'both'):
                                    if is_pos:
                                        correct += 1
                                    else:
                                        wrong_pos += 1
                                        if mapq < 30:
                                            errors_low_mapq += 1
                                if metric in ('overlap', 'both'):
                                    if is_ovl:
                                        correct_overlap += 1
                                    else:
                                        wrong_pos_overlap += 1
                current_s = []
                continue
            if line[0] == 's':
                current_s.append(line.split())

        # Last block
        if len(current_s) >= 2:
            ref_line = read_line = None
            for s in current_s:
                name = s[1]
                if name in ref_names_in_maf and ref_line is None:
                    ref_line = s
                elif name not in ref_names_in_maf and read_line is None:
                    read_line = s
            if ref_line and read_line:
                ref_start = int(ref_line[2])
                ref_size = int(ref_line[3])
                ref_strand = ref_line[4]
                ref_src_size = int(ref_line[5])
                ref_text = ref_line[6] if len(ref_line) > 6 else ""
                read_name = read_line[1]
                read_text = read_line[6] if len(read_line) > 6 else ""
                if ref_strand == '-':
                    ref_start = ref_src_size - ref_start - ref_size
                if ref_text and read_text:
                    true_pos, true_len = find_aligned_ref_interval(ref_text, read_text, ref_start)
                else:
                    true_pos, true_len = ref_start + 1, ref_size
                true_chr = ref_name_override if (ref_name_override and ref_line[1] == 'ref') else ref_line[1]
                aln = alignments.get(read_name)
                if aln is None:
                    zmw = _zmw_id(read_name)
                    if zmw:
                        if zmw not in seen_zmws:
                            seen_zmws.add(zmw)
                            aln = alignments.get(zmw)
                        else:
                            aln = None  # skip duplicate
                if aln is not None or read_name not in seen_zmws:
                    total += 1
                    if aln is not None:
                        mapped += 1
                        aln_chr, aln_pos, mapq, aln_len = aln
                        if aln_chr != true_chr:
                            wrong_chr += 1
                            if mapq < 30:
                                errors_low_mapq += 1
                        else:
                            ov = max(0, min(aln_pos + aln_len, true_pos + true_len) - max(aln_pos, true_pos))
                            is_pos = abs(aln_pos - true_pos) <= tolerance
                            is_ovl = ov >= 0.1 * max(true_len, aln_len)
                            if metric in ('position', 'both'):
                                if is_pos:
                                    correct += 1
                                else:
                                    wrong_pos += 1
                                    if mapq < 30:
                                        errors_low_mapq += 1
                            if metric in ('overlap', 'both'):
                                if is_ovl:
                                    correct_overlap += 1
                                else:
                                    wrong_pos_overlap += 1
    finally:
        fh.close()
        if proc:
            proc.wait()

    return (total, mapped, correct, correct_overlap, wrong_chr, wrong_pos,
            wrong_pos_overlap, errors_low_mapq, len(seen_zmws))


def main():
    global _ALIGNMENTS, _TOLERANCE, _METRIC

    ap = argparse.ArgumentParser(description='Parallel accuracy evaluation against MAF ground truth')
    ap.add_argument('--maf', help='Single MAF file')
    ap.add_argument('--maf-dir', help='Directory with per-chunk MAF(+.ref) files')
    ap.add_argument('--sam', required=True, help='Alignment SAM')
    ap.add_argument('--tolerance', type=int, default=500, help='Position tolerance (bp)')
    ap.add_argument('--metric', choices=['position', 'overlap', 'both'], default='position')
    ap.add_argument('--ref-name', help='Reference name override')
    ap.add_argument('--fai', help='FASTA index (.fai)')
    ap.add_argument('--json', action='store_true', help='Output JSON')
    ap.add_argument('-t', '--threads', type=int, default=8, help='Number of worker processes')
    args = ap.parse_args()

    ref_names = None
    if args.fai:
        ref_names = {'ref'}
        with open(args.fai) as f:
            for line in f:
                ref_names.add(line.split('\t')[0])

    # Load SAM into global dict (will be shared via fork COW)
    print("Loading SAM alignments...", file=sys.stderr)
    _ALIGNMENTS = parse_sam(args.sam)
    print(f"  Loaded {len(_ALIGNMENTS):,} entries", file=sys.stderr)
    _TOLERANCE = args.tolerance
    _METRIC = args.metric

    # Build work items
    work = []
    if args.maf_dir:
        maf_files = sorted(
            p for p in glob.glob(os.path.join(args.maf_dir, '*.maf*'))
            if os.path.isfile(p)
        )
        for maf_path in maf_files:
            ref_name = find_ref_name_for_maf(maf_path)
            if ref_name is None:
                print(f"  SKIP {os.path.basename(maf_path)} (no .ref)", file=sys.stderr)
                continue
            work.append((maf_path, ref_name, ref_names))
    elif args.maf:
        work.append((args.maf, args.ref_name, ref_names))
    else:
        ap.error('--maf-dir or --maf required')

    print(f"  Processing {len(work)} MAF files with {args.threads} workers...", file=sys.stderr)

    # Process in parallel
    total = mapped = correct = correct_overlap = 0
    wrong_chr = wrong_pos = wrong_pos_overlap = errors_low_mapq = 0
    total_zmw_deduped = 0

    with Pool(args.threads) as pool:
        for i, result in enumerate(pool.imap_unordered(process_one_maf, work), 1):
            t, m, c, co, wc, wp, wpo, elm, nzmw = result
            total += t
            mapped += m
            correct += c
            correct_overlap += co
            wrong_chr += wc
            wrong_pos += wp
            wrong_pos_overlap += wpo
            errors_low_mapq += elm
            total_zmw_deduped += nzmw
            if i % 50 == 0 or i == len(work):
                print(f"  [{i}/{len(work)}] total={total:,} mapped={mapped:,}", file=sys.stderr)

    total_errors = wrong_chr + wrong_pos
    results = {
        'total': total,
        'mapped': mapped,
        'correct': correct,
        'correct_overlap': correct_overlap,
        'wrong_chr': wrong_chr,
        'wrong_pos': wrong_pos,
        'wrong_pos_overlap': wrong_pos_overlap,
        'unmapped': total - mapped,
        'accuracy': 100.0 * correct / mapped if mapped > 0 else 0.0,
        'accuracy_overlap': 100.0 * correct_overlap / mapped if mapped > 0 else 0.0,
        'mapping_rate': 100.0 * mapped / total if total > 0 else 0.0,
        'error_detection': 100.0 * errors_low_mapq / total_errors if total_errors > 0 else None,
    }

    if args.json:
        import json
        print(json.dumps(results, indent=2))
    else:
        print(f"Total reads:      {results['total']:,}")
        print(f"Mapped:           {results['mapped']:,}")
        if args.metric in ('position', 'both'):
            print(f"Correct (pos):    {results['correct']:,} ({results['accuracy']:.2f}%)")
        if args.metric in ('overlap', 'both'):
            print(f"Correct (ovl):    {results['correct_overlap']:,} ({results['accuracy_overlap']:.2f}%)")
        print(f"Wrong chromosome: {results['wrong_chr']:,}")
        print(f"Unmapped:         {results['unmapped']:,}")
        print(f"Mapping rate:     {results['mapping_rate']:.2f}%")
        if results['error_detection'] is not None:
            print(f"Error detection:  {results['error_detection']:.1f}% (MAPQ<30)")
        if total_zmw_deduped:
            print(f"ZMW deduped:      {total_zmw_deduped:,}")


if __name__ == '__main__':
    main()
