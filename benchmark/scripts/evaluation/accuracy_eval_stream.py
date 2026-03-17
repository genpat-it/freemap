#!/usr/bin/env python3
"""
accuracy_eval_stream.py - Evaluate alignment accuracy against ground truth MAF (streaming)

Memory-efficient: loads SAM (small) first, then streams MAF files one at a time,
evaluating each read on-the-fly without accumulating the full truth dict.

Metrics: position (default) or overlap (mm2-like 10%).
Primary alignments only; requires --ref-name to map MAF 'ref' to SAM names.
"""

import argparse
import gzip
import glob
import os
import sys
from typing import Dict, Tuple, Optional


def find_aligned_ref_interval(ref_text: str, read_text: str, ref_start_0based: int) -> Tuple[int, int]:
    """Return (start_1based, length) of aligned ref interval from MAF texts."""
    min_len = min(len(ref_text), len(read_text))
    first_pos = None
    last_pos = None
    ref_offset = 0

    for i in range(min_len):
        r, q = ref_text[i], read_text[i]
        if r != '-' and q != '-':
            current_pos = ref_start_0based + ref_offset
            if first_pos is None:
                first_pos = current_pos
            last_pos = current_pos
        if r != '-':
            ref_offset += 1

    if first_pos is None:
        return ref_start_0based + 1, 1
    return first_pos + 1, last_pos - first_pos + 1


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


def _zmw_id(name: str) -> Optional[str]:
    """Extract ZMW ID from PacBio-style read names.
    S1/123/ccs -> S1/123,  S1/123/0 -> S1/123,  other -> None."""
    parts = name.rsplit('/', 1)
    if len(parts) == 2 and (parts[1] == 'ccs' or parts[1].isdigit()):
        return parts[0]
    return None


def parse_sam_stream(sam_file: str) -> Dict[str, Tuple[str, int, int, int]]:
    """Primary alignments only: read_name -> (ref_name, pos, mapq, ref_len).
    Also indexes by ZMW ID for HiFi CCS reads (S1/N/ccs -> S1/N)."""
    alignments = {}
    fh = gzip.open(sam_file, 'rt') if sam_file.endswith('.gz') else open(sam_file, 'r')
    with fh:
        for line in fh:
            if line.startswith('@') or line.startswith('['):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                continue
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            if flag & 4:
                continue
            if flag & 0x900:
                continue
            read_name = fields[0]
            ref_name = fields[2]
            try:
                pos = int(fields[3])
                mapq = int(fields[4])
            except ValueError:
                continue
            cigar = fields[5]
            ref_len = cigar_ref_length(cigar) if cigar != '*' else 0
            val = (ref_name, pos, mapq, ref_len)
            alignments[read_name] = val
            # Also index by ZMW ID for CCS reads
            zmw = _zmw_id(read_name)
            if zmw and zmw not in alignments:
                alignments[zmw] = val
    return alignments


def compute_overlap(start1: int, len1: int, start2: int, len2: int) -> int:
    end1 = start1 + len1
    end2 = start2 + len2
    return max(0, min(end1, end2) - max(start1, start2))


def find_ref_name_for_maf(maf_path: str) -> Optional[str]:
    """Return ref name from companion .ref file, or None if missing."""
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


class AccuracyCounters:
    """Accumulates accuracy statistics without storing per-read data."""
    __slots__ = ('total', 'mapped', 'correct', 'correct_overlap',
                 'wrong_chr', 'wrong_pos', 'wrong_pos_overlap', 'errors_low_mapq')

    def __init__(self):
        self.total = 0
        self.mapped = 0
        self.correct = 0
        self.correct_overlap = 0
        self.wrong_chr = 0
        self.wrong_pos = 0
        self.wrong_pos_overlap = 0
        self.errors_low_mapq = 0

    def evaluate_read(self, true_chr, true_pos, true_len, aln, tolerance, metric):
        """Evaluate a single read against its alignment."""
        self.total += 1
        if aln is None:
            return
        self.mapped += 1
        aln_chr, aln_pos, mapq, aln_len = aln
        if aln_chr != true_chr:
            self.wrong_chr += 1
            if mapq < 30:
                self.errors_low_mapq += 1
        else:
            overlap = compute_overlap(true_pos, true_len, aln_pos, aln_len)
            is_correct_pos = abs(aln_pos - true_pos) <= tolerance
            is_correct_ovl = overlap >= 0.1 * max(true_len, aln_len)

            if metric in ('position', 'both'):
                if is_correct_pos:
                    self.correct += 1
                else:
                    self.wrong_pos += 1
                    if mapq < 30:
                        self.errors_low_mapq += 1

            if metric in ('overlap', 'both'):
                if is_correct_ovl:
                    self.correct_overlap += 1
                else:
                    self.wrong_pos_overlap += 1

    def results(self):
        total_errors = self.wrong_chr + self.wrong_pos
        return {
            'total': self.total,
            'mapped': self.mapped,
            'correct': self.correct,
            'correct_overlap': self.correct_overlap,
            'wrong_chr': self.wrong_chr,
            'wrong_pos': self.wrong_pos,
            'wrong_pos_overlap': self.wrong_pos_overlap,
            'unmapped': self.total - self.mapped,
            'accuracy': 100.0 * self.correct / self.mapped if self.mapped > 0 else 0.0,
            'accuracy_overlap': 100.0 * self.correct_overlap / self.mapped if self.mapped > 0 else 0.0,
            'mapping_rate': 100.0 * self.mapped / self.total if self.total > 0 else 0.0,
            'error_detection': 100.0 * self.errors_low_mapq / total_errors if total_errors > 0 else None,
        }


def _process_maf_entry(current_s, ref_names_in_maf, ref_name_override):
    """Extract (read_name, ref_name, true_pos, true_len) from a MAF alignment block."""
    ref_line = None
    read_line = None
    for s in current_s:
        name = s[1]
        if name in ref_names_in_maf and ref_line is None:
            ref_line = s
        elif name not in ref_names_in_maf and read_line is None:
            read_line = s
    if not ref_line or not read_line:
        return None

    maf_ref_name = ref_line[1]
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

    out_ref = ref_name_override if (ref_name_override and maf_ref_name == 'ref') else maf_ref_name
    return read_name, out_ref, true_pos, true_len


def stream_evaluate_maf(maf_file: str, alignments: dict, counters: AccuracyCounters,
                        tolerance: int, metric: str,
                        ref_name_override: str = None, ref_names: set = None,
                        seen_zmws: set = None):
    """Stream-parse a single MAF file, evaluating reads on-the-fly against alignments dict.
    No truth dict is built — each read is evaluated and discarded immediately.
    seen_zmws tracks ZMW IDs already evaluated (for HiFi dedup across subreads)."""
    ref_names_in_maf = ref_names.copy() if ref_names else {'ref'}
    if ref_name_override:
        ref_names_in_maf.add(ref_name_override)

    fh = gzip.open(maf_file, 'rt') if maf_file.endswith('.gz') else open(maf_file, 'r')
    with fh:
        current_s = []
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('a'):
                if len(current_s) >= 2:
                    result = _process_maf_entry(current_s, ref_names_in_maf, ref_name_override)
                    if result:
                        read_name, true_chr, true_pos, true_len = result
                        aln = alignments.get(read_name)
                        if aln is None:
                            zmw = _zmw_id(read_name)
                            if zmw:
                                if seen_zmws is not None:
                                    if zmw in seen_zmws:
                                        current_s = []
                                        continue
                                    seen_zmws.add(zmw)
                                aln = alignments.get(zmw)
                        counters.evaluate_read(true_chr, true_pos, true_len, aln, tolerance, metric)
                current_s = []
                continue
            if line.startswith('s'):
                current_s.append(line.split())

        # Process last entry
        if len(current_s) >= 2:
            result = _process_maf_entry(current_s, ref_names_in_maf, ref_name_override)
            if result:
                read_name, true_chr, true_pos, true_len = result
                aln = alignments.get(read_name)
                if aln is None:
                    zmw = _zmw_id(read_name)
                    if zmw:
                        if seen_zmws is not None:
                            if zmw in seen_zmws:
                                return
                            seen_zmws.add(zmw)
                        aln = alignments.get(zmw)
                counters.evaluate_read(true_chr, true_pos, true_len, aln, tolerance, metric)


def main():
    ap = argparse.ArgumentParser(description='Evaluate alignment accuracy (streaming, memory-efficient)')
    ap.add_argument('--maf', help='Ground truth MAF (gz or plain)')
    ap.add_argument('--maf-dir', help='Directory with per-chunk MAF(+.ref) files')
    ap.add_argument('--sam', required=True, help='Alignment SAM/BAM (gz or plain SAM)')
    ap.add_argument('--tolerance', type=int, default=500, help='Position tolerance (bp)')
    ap.add_argument('--metric', choices=['position', 'overlap', 'both'], default='position',
                    help='Correctness metric: position, overlap>=10%, or both')
    ap.add_argument('--ref-name',
                    help='Reference name to map MAF "ref" to (e.g., NC_000913.3 or chr22)')
    ap.add_argument('--fai', help='FASTA index (.fai) to load all reference names')
    ap.add_argument('--json', action='store_true', help='Output JSON')
    args = ap.parse_args()

    # Load reference names from .fai if provided
    ref_names = None
    if args.fai:
        ref_names = {'ref'}
        with open(args.fai) as f:
            for line in f:
                ref_names.add(line.split('\t')[0])

    # Step 1: Load SAM alignments (memory-efficient: only stores name + 4 values per read)
    print("Loading SAM alignments...", file=sys.stderr)
    alignments = parse_sam_stream(args.sam)
    print(f"  Loaded {len(alignments):,} primary alignments", file=sys.stderr)

    # Step 2: Stream MAF files, evaluating each read on-the-fly
    counters = AccuracyCounters()
    seen_zmws = set()  # Dedup HiFi subreads: only count first subread per ZMW

    if args.maf_dir:
        maf_files = sorted(
            p for p in glob.glob(os.path.join(args.maf_dir, '*.maf*'))
            if os.path.isfile(p)
        )
        if not maf_files:
            ap.error(f'No MAF files found in {args.maf_dir}')
        for i, maf_path in enumerate(maf_files, 1):
            ref_name = find_ref_name_for_maf(maf_path)
            if ref_name is None:
                print(f"  [{i}/{len(maf_files)}] {os.path.basename(maf_path)} — SKIPPED (no .ref)",
                      file=sys.stderr)
                continue
            print(f"  [{i}/{len(maf_files)}] {os.path.basename(maf_path)} (ref={ref_name})",
                  file=sys.stderr)
            stream_evaluate_maf(maf_path, alignments, counters, args.tolerance, args.metric,
                                ref_name, ref_names, seen_zmws)
    else:
        if not args.maf or (not args.ref_name and not args.fai):
            ap.error('Either --maf-dir, or --maf with --ref-name/--fai are required')
        stream_evaluate_maf(args.maf, alignments, counters, args.tolerance, args.metric,
                            args.ref_name, ref_names, seen_zmws)

    results = counters.results()

    if args.json:
        import json
        print(json.dumps(results, indent=2))
    else:
        if args.metric == 'both':
            print("Metric:           pos±{0}bp and overlap>=10%".format(args.tolerance))
        else:
            metric_desc = "overlap>=10%" if args.metric == 'overlap' else f"pos±{args.tolerance}bp"
            print(f"Metric:           {metric_desc}")
        print(f"Total reads:      {results['total']:,}")
        print(f"Mapped:           {results['mapped']:,}")
        if args.metric in ('position', 'both'):
            print(f"Correct (pos):    {results['correct']:,} ({results['accuracy']:.2f}%)")
            print(f"Wrong position:   {results['wrong_pos']:,}")
        if args.metric in ('overlap', 'both'):
            print(f"Correct (ovl):    {results['correct_overlap']:,} ({results['accuracy_overlap']:.2f}%)")
            print(f"Wrong ovl:        {results['wrong_pos_overlap']:,}")
        print(f"Wrong chromosome: {results['wrong_chr']:,}")
        print(f"Unmapped:         {results['unmapped']:,}")
        print(f"Mapping rate:     {results['mapping_rate']:.2f}%")
        if results['error_detection'] is None:
            print("Error detection:  N/A (no errors)")
        else:
            print(f"Error detection:  {results['error_detection']:.1f}% (MAPQ<30)")


if __name__ == '__main__':
    main()
