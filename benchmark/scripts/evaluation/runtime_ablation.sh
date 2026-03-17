#!/bin/bash
# runtime_ablation.sh — Reproduce SM11 runtime ablation tables
#
# Usage:
#   bash benchmark/scripts/evaluation/runtime_ablation.sh [DATA_DIR]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

FREEMAP="$PROJECT_DIR/target/release/freemap"
T="${THREADS:-128}"

if [[ ! -x "$FREEMAP" ]]; then
    log_error "freemap not found — run: cargo build --release"
    exit 1
fi

run_fm() {
    local label="$1"; shift
    echo "  $label"
    "$FREEMAP" -t "$T" "$@" 2>&1 | grep -E 'Align time|Output time|Total time|Throughput'
    echo ""
}

# ── E. coli ONT ─────────────────────────────────────────────────────────────

REF_E="$DATA_DIR/references/ecoli_ref.fa"
READS_E="$DATA_DIR/L1_ecoli_ont_50x/reads_0001.fq"

if [[ -f "$REF_E" && -f "$READS_E" ]]; then
    echo "========================================================="
    echo " SM11a — E. coli ONT (23k reads)"
    echo "========================================================="
    echo ""
    run_fm "1. Seed+chain only (PAF, no CIGAR)" -k 15 -w 10 -f 20 "$REF_E" "$READS_E" /dev/null
    run_fm "2. + Trajectory CIGAR (-g)" -k 15 -w 10 -f 20 -g "$REF_E" "$READS_E" /dev/null
    run_fm "3. + Polish (-g -p)" -k 15 -w 10 -f 20 -g -p "$REF_E" "$READS_E" /dev/null
    run_fm "4. + SAM output (-g -p -a)" -k 15 -w 10 -f 20 -g -p -a "$REF_E" "$READS_E" /dev/null

    if command -v minimap2 &>/dev/null; then
        echo "  5. minimap2 PAF (no CIGAR)"
        /usr/bin/time -f "  Wall: %e s, CPU: %U+%S s" minimap2 -t "$T" -x map-ont "$REF_E" "$READS_E" -o /dev/null 2>&1
        echo ""
        echo "  6. minimap2 SAM (with DP)"
        /usr/bin/time -f "  Wall: %e s, CPU: %U+%S s" minimap2 -a -t "$T" -x map-ont "$REF_E" "$READS_E" -o /dev/null 2>&1
        echo ""
    else
        log_warn "minimap2 not found, skipping comparison"
    fi
else
    log_warn "E. coli data not found, skipping"
fi

# ── Human HiFi ──────────────────────────────────────────────────────────────

REF_H="$DATA_DIR/references/GRCh38.fa"
READS_H="$DATA_DIR/H1_human_hifi_5x/reads_100k.fq"
IDX_FM="$DATA_DIR/references/GRCh38_maphifi.fmi"
IDX_MM2="$DATA_DIR/references/GRCh38_hifi.mmi"

if [[ -f "$REF_H" && -f "$READS_H" ]]; then
    echo "========================================================="
    echo " SM11b — Human HiFi (100k reads, GRCh38)"
    echo "========================================================="
    echo ""

    FM_IDX_ARGS=""
    [[ -f "$IDX_FM" ]] && FM_IDX_ARGS="-i $IDX_FM"

    run_fm "1. Seed+chain only (PAF, no CIGAR)" -k 19 -w 19 -f 20 $FM_IDX_ARGS "$REF_H" "$READS_H" /dev/null
    run_fm "2. + Trajectory CIGAR (-g)" -k 19 -w 19 -f 20 -g $FM_IDX_ARGS "$REF_H" "$READS_H" /dev/null
    run_fm "3. + Polish (-g -p)" -k 19 -w 19 -f 20 -g -p $FM_IDX_ARGS "$REF_H" "$READS_H" /dev/null
    run_fm "4. + SAM output (-g -p -a)" -k 19 -w 19 -f 20 -g -p -a $FM_IDX_ARGS "$REF_H" "$READS_H" /dev/null

    if command -v minimap2 &>/dev/null; then
        MM2_ARGS=""
        [[ -f "$IDX_MM2" ]] && MM2_ARGS="$IDX_MM2" || MM2_ARGS="$REF_H"

        echo "  5. minimap2 PAF (no CIGAR)"
        /usr/bin/time -f "  Wall: %e s, CPU: %U+%S s" minimap2 -t "$T" -x map-hifi "$MM2_ARGS" "$READS_H" -o /dev/null 2>&1
        echo ""
        echo "  6. minimap2 SAM (with DP)"
        /usr/bin/time -f "  Wall: %e s, CPU: %U+%S s" minimap2 -a -t "$T" -x map-hifi "$MM2_ARGS" "$READS_H" -o /dev/null 2>&1
        echo ""
    else
        log_warn "minimap2 not found, skipping comparison"
    fi
else
    log_warn "Human HiFi data not found, skipping"
fi

echo "========================================================="
echo " Done."
echo "========================================================="
