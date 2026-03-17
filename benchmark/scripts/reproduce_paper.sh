#!/bin/bash
# reproduce_paper.sh - Reproduce all freemap paper results
# =========================================================
#
# Single entry point for reproducing every table and figure in the paper.
#
# Usage:
#   ./reproduce_paper.sh [DATA_DIR] [SECTION]
#
#   DATA_DIR  Path to benchmark data (default: ../../data relative to this script,
#             or $DATA_DIR env var)
#   SECTION   Run a specific section instead of all (see list below)
#
# Examples:
#   ./reproduce_paper.sh                          # All results, default data path
#   ./reproduce_paper.sh /mnt/data                # All results, custom data path
#   ./reproduce_paper.sh /mnt/data table1         # Only Table 1
#   ./reproduce_paper.sh . table3                  # Table 3, data in current dir
#
# Sections: table1 table2 table3 table4 table5 figure1 figure3
#           supp_short supp_concordance supp_coverage
#
# Requirements:
#   - Python 3.8+ with numpy, matplotlib, scipy, intervaltree
#   - samtools (for coverage analysis)
#   - Pre-computed alignment files (SAM/BAM, MAF, depth) in DATA_DIR
#
# Data layout expected under DATA_DIR:
#   L1_ecoli_ont_50x/        freemap.sam, minimap2.sam, *.maf[.gz]
#   L2_ecoli_hifi_50x/       freemap_ccs.sam, minimap2.sam, *.maf[.gz]
#   L3_ecoli_clr_50x/        freemap.sam, minimap2.sam, *.maf[.gz]
#   H1_human_hifi_5x/        freemap.sam, minimap2.sam, *.maf[.gz]
#   H1_human_hifi_5x_chr22/  freemap.sam, minimap2.sam, *.maf[.gz]
#   H4_human_ont_10x/        freemap.sam, minimap2.sam, *.maf[.gz]
#   H4_human_ont_10x_chr22/  freemap.sam, minimap2.sam, *.maf[.gz]
#   H5_human_clr_5x/         freemap.sam, minimap2.sam, *.maf[.gz]
#   S2_human_short_10M/      freemap.sam, minimap2.sam, truth.sam
#   D1_giab_hg002_hifi/      freemap.sam, minimap2.sam
#   D2_giab_hg002_ont/       freemap.sam, minimap2.sam
#   giab_truth/               HG002 high-confidence BED

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# Known sections for matching
KNOWN_SECTIONS=" all table1 table2 table3 table4 table5 table6 figure1 figure2 figure3 supp_short supp_concordance supp_coverage "

# Parse positional arguments (order-independent)
SECTION=""
for arg in "$@"; do
    if [[ "$KNOWN_SECTIONS" == *" $arg "* ]]; then
        SECTION="$arg"
    elif [[ -d "$arg" ]]; then
        DATA_DIR="$(cd "$arg" && pwd)"
        export DATA_DIR
    else
        log_error "Unknown argument: $arg (not a known section or existing directory)"
        echo ""
        echo "Available sections:"
        echo "  table1 table2 table3 table4 table5"
        echo "  figure1 figure3"
        echo "  supp_short supp_concordance supp_coverage"
        echo "  all (default)"
        exit 1
    fi
done
SECTION="${SECTION:-all}"

RESULTS_DIR="$PROJECT_DIR/benchmark/results"
mkdir -p "$RESULTS_DIR"

# ──────────────────────────────────────────────
# Prerequisites
# ──────────────────────────────────────────────

check_prerequisites() {
    log_info "Checking prerequisites..."
    log_info "DATA_DIR = $DATA_DIR"
    log_info "PROJECT_DIR = $PROJECT_DIR"

    if [[ ! -d "$DATA_DIR" ]]; then
        log_error "DATA_DIR not found: $DATA_DIR"
        log_error "Set DATA_DIR or pass it as first argument."
        exit 1
    fi

    if ! command -v python3 &>/dev/null; then
        log_error "python3 not found"
        exit 1
    fi

    log_info "python3: $(python3 --version 2>&1)"
    log_info "Prerequisites OK"
    echo ""
}

# ──────────────────────────────────────────────
# Table 1: E. coli simulated long-read accuracy
# ──────────────────────────────────────────────

table1() {
    log_info "=== Table 1: E. coli simulated accuracy ==="

    local datasets=("L1_ecoli_ont_50x" "L2_ecoli_hifi_50x" "L3_ecoli_clr_50x")

    for ds in "${datasets[@]}"; do
        local ds_dir="$DATA_DIR/$ds"
        if [[ ! -d "$ds_dir" ]]; then
            log_warn "Skipping $ds (directory not found)"
            continue
        fi

        echo ""
        log_info "--- $ds ---"

        # Find MAF files
        local maf_args=""
        if ls "$ds_dir"/*.maf.gz &>/dev/null 2>&1; then
            for maf in "$ds_dir"/*.maf.gz; do
                maf_args="--maf $maf"
                break
            done
        elif ls "$ds_dir"/*.maf &>/dev/null 2>&1; then
            for maf in "$ds_dir"/*.maf; do
                maf_args="--maf $maf"
                break
            done
        fi

        if [[ -z "$maf_args" ]]; then
            # Try maf-dir for per-chunk MAFs
            if [[ -d "$ds_dir/maf" ]]; then
                maf_args="--maf-dir $ds_dir/maf"
            else
                log_warn "No MAF files found for $ds"
                continue
            fi
        fi

        for aligner in freemap minimap2; do
            local sam="$ds_dir/${aligner}.sam"
            # Try alternative names
            [[ ! -f "$sam" && -f "$ds_dir/${aligner}_ccs.sam" ]] && sam="$ds_dir/${aligner}_ccs.sam"
            if [[ ! -f "$sam" ]]; then
                log_warn "  $aligner: SAM not found"
                continue
            fi

            echo "  $aligner:"
            python3 "$EVAL_DIR/accuracy_eval_stream.py" \
                $maf_args \
                --sam "$sam" \
                --ref-name NC_000913.3 \
                --metric position
        done
    done
}

# ──────────────────────────────────────────────
# Table 2: Human simulated long-read accuracy
# ──────────────────────────────────────────────

table2() {
    log_info "=== Table 2: Human simulated accuracy ==="
    log_info "NOTE: Full-genome datasets have hundreds of MAF chunks (~95 GB compressed)."
    log_info "      Parsing takes 1-2 hours per dataset. This is expected."

    local datasets=("H1_human_hifi_5x" "H4_human_ont_10x" "H5_human_clr_5x")

    for ds in "${datasets[@]}"; do
        local ds_dir="$DATA_DIR/$ds"
        if [[ ! -d "$ds_dir" ]]; then
            log_warn "Skipping $ds (directory not found)"
            continue
        fi

        echo ""
        log_info "--- $ds ---"

        # Determine MAF source: single truth.maf, maf/ subdir, or chunked reads_*.maf.gz
        local maf_args=""
        if [[ -f "$ds_dir/truth.maf" ]]; then
            maf_args="--maf $ds_dir/truth.maf"
        elif [[ -d "$ds_dir/maf" ]]; then
            maf_args="--maf-dir $ds_dir/maf"
        elif ls "$ds_dir"/reads_*.maf.gz &>/dev/null 2>&1; then
            maf_args="--maf-dir $ds_dir"
        elif ls "$ds_dir"/*.maf.gz &>/dev/null 2>&1; then
            for maf in "$ds_dir"/*.maf.gz; do
                maf_args="--maf $maf"
                break
            done
        fi

        if [[ -z "$maf_args" ]]; then
            log_warn "No MAF files found for $ds"
            continue
        fi

        for aligner in freemap minimap2; do
            local sam="$ds_dir/${aligner}.sam"
            if [[ ! -f "$sam" ]]; then
                log_warn "  $aligner: SAM not found"
                continue
            fi

            echo "  $aligner:"
            python3 "$EVAL_DIR/accuracy_eval_stream.py" \
                $maf_args \
                --sam "$sam" \
                --fai "$DATA_DIR/references/GRCh38.fa.fai" \
                --metric position
        done
    done
}

# ──────────────────────────────────────────────
# Table 3: Concordance (simulated + real)
# ──────────────────────────────────────────────

table3() {
    log_info "=== Table 3: Pairwise concordance ==="

    local datasets=(
        "L1_ecoli_ont_50x"
        "L2_ecoli_hifi_50x"
        "L3_ecoli_clr_50x"
        "H1_human_hifi_5x"
        "H4_human_ont_10x"
        "H5_human_clr_5x"
        "D1_giab_hg002_hifi"
        "D2_giab_hg002_ont"
    )

    for ds in "${datasets[@]}"; do
        local ds_dir="$DATA_DIR/$ds"
        if [[ ! -d "$ds_dir" ]]; then
            log_warn "Skipping $ds"
            continue
        fi

        local fm_sam="$ds_dir/freemap.sam"
        local mm2_sam="$ds_dir/minimap2.sam"
        [[ ! -f "$fm_sam" && -f "$ds_dir/freemap_ccs.sam" ]] && fm_sam="$ds_dir/freemap_ccs.sam"

        if [[ ! -f "$fm_sam" || ! -f "$mm2_sam" ]]; then
            log_warn "Skipping $ds (SAM files missing)"
            continue
        fi

        echo ""
        log_info "--- $ds ---"
        python3 "$EVAL_DIR/concordance_eval.py" \
            --sam1 "$fm_sam" --sam2 "$mm2_sam" \
            --name1 freemap --name2 minimap2 \
            --metric position
    done
}

# ──────────────────────────────────────────────
# Table 4: CIGAR structural concordance
# ──────────────────────────────────────────────

table4() {
    log_info "=== Table 4: CIGAR structural concordance ==="

    local datasets=(
        "L1_ecoli_ont_50x"
        "L2_ecoli_hifi_50x"
        "L3_ecoli_clr_50x"
        "H1_human_hifi_5x"
        "H4_human_ont_10x"
        "H5_human_clr_5x"
    )

    for ds in "${datasets[@]}"; do
        local ds_dir="$DATA_DIR/$ds"
        if [[ ! -d "$ds_dir" ]]; then
            log_warn "Skipping $ds"
            continue
        fi

        local fm_sam="$ds_dir/freemap.sam"
        local mm2_sam="$ds_dir/minimap2.sam"
        [[ ! -f "$fm_sam" && -f "$ds_dir/freemap_ccs.sam" ]] && fm_sam="$ds_dir/freemap_ccs.sam"

        if [[ ! -f "$fm_sam" || ! -f "$mm2_sam" ]]; then
            log_warn "Skipping $ds (SAM files missing)"
            continue
        fi

        echo ""
        log_info "--- $ds ---"
        python3 "$EVAL_DIR/cigar_stats.py" \
            --compare "$fm_sam" "$mm2_sam" \
            --name1 freemap --name2 minimap2 \
            --concordant
    done
}

# ──────────────────────────────────────────────
# Tables 5 & 6: GIAB high-confidence MAPQ analysis
# ──────────────────────────────────────────────

table5() {
    log_info "=== Tables 5 & 6: GIAB MAPQ analysis ==="
    python3 "$FIG_DIR/giab_mapq_analysis.py"
}

# ──────────────────────────────────────────────
# Figures 1 & 2: Coverage scatter
# ──────────────────────────────────────────────

figure1() {
    log_info "=== Figures 1 & 2: Coverage scatter ==="
    python3 "$FIG_DIR/coverage_scatter_paper.py" "$RESULTS_DIR/fig_coverage_scatter.pdf"
    log_info "Saved: $RESULTS_DIR/fig_coverage_scatter.pdf"
    log_info "Saved: $RESULTS_DIR/fig_coverage_scatter_giab.pdf"
}

# ──────────────────────────────────────────────
# Figure 3: MAPQ calibration
# ──────────────────────────────────────────────

figure3() {
    log_info "=== Figure 3: MAPQ calibration ==="
    python3 "$FIG_DIR/mapq_calibration_multi.py" "$RESULTS_DIR/fig_mapq_calibration.pdf"
    log_info "Saved: $RESULTS_DIR/fig_mapq_calibration.pdf"
}

# ──────────────────────────────────────────────
# Supplementary: Short-read accuracy (Tables S2, S3)
# ──────────────────────────────────────────────

supp_short() {
    log_info "=== Supplementary: Short-read accuracy (Tables S2, S3) ==="

    local ds_dir="$DATA_DIR/S2_human_short_10M"
    if [[ ! -d "$ds_dir" ]]; then
        log_warn "Skipping: $ds_dir not found"
        return
    fi

    local truth="$ds_dir/truth.sam"
    if [[ ! -f "$truth" ]]; then
        log_warn "truth.sam not found in $ds_dir"
        return
    fi

    for aligner in freemap minimap2 bwa-mem2; do
        local sam="$ds_dir/${aligner}.sam"
        if [[ ! -f "$sam" ]]; then
            log_warn "  $aligner: SAM not found"
            continue
        fi

        echo ""
        echo "  $aligner:"
        python3 "$EVAL_DIR/accuracy_eval_short.py" \
            --truth "$truth" --sam "$sam" --metric position
    done
}

# ──────────────────────────────────────────────
# Supplementary: Real-data concordance (Table S4, S6)
# ──────────────────────────────────────────────

supp_concordance() {
    log_info "=== Supplementary: Real-data concordance ==="
    python3 "$EVAL_DIR/concordance_real.py" "$RESULTS_DIR/fig_concordance_real.pdf"
    log_info "Saved: $RESULTS_DIR/fig_concordance_real.pdf"
}

# ──────────────────────────────────────────────
# Supplementary: Coverage correlation (Table S12)
# ──────────────────────────────────────────────

supp_coverage() {
    log_info "=== Supplementary: Coverage correlation (Table S12) ==="

    local datasets=(
        "L1_ecoli_ont_50x"
        "L2_ecoli_hifi_50x"
        "L3_ecoli_clr_50x"
        "H1_human_hifi_5x"
        "H4_human_ont_10x"
    )

    for ds in "${datasets[@]}"; do
        local ds_dir="$DATA_DIR/$ds"
        if [[ ! -d "$ds_dir" ]]; then
            log_warn "Skipping $ds"
            continue
        fi

        # Find depth files
        local fm_depth="" mm2_depth=""
        for f in freemap.depth freemap_chr22.depth freemap_ccs.depth; do
            [[ -f "$ds_dir/$f" ]] && fm_depth="$ds_dir/$f" && break
        done
        for f in minimap2.depth minimap2_chr22.depth; do
            [[ -f "$ds_dir/$f" ]] && mm2_depth="$ds_dir/$f" && break
        done

        if [[ -z "$fm_depth" || -z "$mm2_depth" ]]; then
            log_warn "Depth files missing for $ds"
            continue
        fi

        echo ""
        log_info "--- $ds (per-base) ---"
        python3 "$FIG_DIR/coverage_corr.py" "$fm_depth" "$mm2_depth"

        echo ""
        log_info "--- $ds (1kb bins) ---"
        python3 "$FIG_DIR/coverage_corr_binned.py" "$fm_depth" "$mm2_depth" 1000
    done
}

# ──────────────────────────────────────────────
# Run all
# ──────────────────────────────────────────────

run_all() {
    table1
    echo ""
    table2
    echo ""
    table3
    echo ""
    table4
    echo ""
    table5
    echo ""
    figure1
    echo ""
    figure3
    echo ""
    supp_short
    echo ""
    supp_concordance
    echo ""
    supp_coverage
}

# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

echo "========================================"
echo " freemap — Paper Reproduction"
echo "========================================"
echo " Data:    $DATA_DIR"
echo " Results: $RESULTS_DIR"
echo " Section: $SECTION"
echo "========================================"
echo ""

check_prerequisites

case "$SECTION" in
    all)              run_all ;;
    table1)           table1 ;;
    table2)           table2 ;;
    table3)           table3 ;;
    table4)           table4 ;;
    table5|table6)    table5 ;;
    figure1|figure2)  figure1 ;;
    figure3)          figure3 ;;
    supp_short)       supp_short ;;
    supp_concordance) supp_concordance ;;
    supp_coverage)    supp_coverage ;;
    *)
        log_error "Unknown section: $SECTION"
        echo ""
        echo "Available sections:"
        echo "  table1 table2 table3 table4 table5"
        echo "  figure1 figure3"
        echo "  supp_short supp_concordance supp_coverage"
        echo "  all (default)"
        exit 1
        ;;
esac

echo ""
echo "========================================"
log_info "Done."
echo "========================================"
