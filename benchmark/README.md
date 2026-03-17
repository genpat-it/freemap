# freemap Benchmark Suite

Reproducibility scripts for the freemap paper. One script reproduces all tables and figures.

## Quick Start

```bash
# 1. Build freemap
cargo build --release

# 2. Reproduce all paper results
bash benchmark/scripts/reproduce_paper.sh /path/to/data

# 3. Reproduce a single section
bash benchmark/scripts/reproduce_paper.sh /path/to/data table1
```

`/path/to/data` is wherever your benchmark datasets live (default: `data/` symlink in project root). You can also set `DATA_DIR` as an environment variable.

## Available Sections

| Section | Paper Element |
|---------|---------------|
| `table1` | Table 1 — E. coli simulated long-read accuracy |
| `table2` | Table 2 — Human simulated long-read accuracy |
| `table3` | Table 3 — Pairwise concordance (simulated + real) |
| `table4` | Table 4 — CIGAR structural concordance |
| `table5` | Tables 5 & 6 — GIAB high-confidence MAPQ analysis |
| `figure1` | Figures 1 & 2 — Coverage scatter plots |
| `figure3` | Figure 3 — MAPQ calibration |
| `supp_short` | Tables S2, S3 — Short-read accuracy |
| `supp_concordance` | Tables S4, S6 — Real-data concordance |
| `supp_coverage` | Table S12 — Coverage correlation |
| `all` | Everything (default) |

## Directory Structure

```
benchmark/scripts/
  reproduce_paper.sh        # Single entry point
  config.sh                 # Shared configuration (paths, logging)
  evaluation/               # Core analysis scripts
    accuracy_eval_stream.py  # Simulated long-read accuracy (Tables 1, 2)
    accuracy_eval_short.py   # Simulated short-read accuracy (Tables S2, S3)
    concordance_eval.py      # Pairwise concordance (Tables 3, S4, S6)
    concordance_real.py      # Real-data concordance helper
    cigar_stats.py           # CIGAR structural analysis (Tables 4, S5)
    concordant_coverage_v2.py  # Concordant-only coverage (text values)
  figures/                   # Figure and table generation
    coverage_scatter_paper.py  # Figures 1 & 2
    coverage_corr.py           # Table S12 (per-base)
    coverage_corr_binned.py    # Table S12 (1kb bins)
    mapq_calibration_multi.py  # Figure 3
    giab_mapq_analysis.py      # Tables 5 & 6
  data_prep/                 # Data preparation
    subsample_genome_wide.py   # Generate subsampled read files
```

## Expected Data Layout

All scripts expect the following directories under `DATA_DIR`:

| Dataset | Contents |
|---------|----------|
| `L1_ecoli_ont_50x/` | `freemap.sam`, `minimap2.sam`, `*.maf[.gz]`, `*.depth` |
| `L2_ecoli_hifi_50x/` | `freemap_ccs.sam`, `minimap2.sam`, `*.maf[.gz]`, `*.depth` |
| `L3_ecoli_clr_50x/` | `freemap.sam`, `minimap2.sam`, `*.maf[.gz]`, `*.depth` |
| `H1_human_hifi_5x/` | `freemap.sam`, `minimap2.sam`, MAF dir, `*_chr22.depth` |
| `H4_human_ont_10x/` | `freemap.sam`, `minimap2.sam`, MAF dir, `*_chr22.depth` |
| `H5_human_clr_5x/` | `freemap.sam`, `minimap2.sam`, MAF dir |
| `S2_human_short_10M/` | `freemap.sam`, `minimap2.sam`, `bwa-mem2.sam`, `truth.sam` |
| `D1_giab_hg002_hifi/` | `freemap.sam`, `minimap2.sam`, `*_chr1.depth` |
| `D2_giab_hg002_ont/` | `freemap.sam`, `minimap2.sam`, `*_chr1.depth` |
| `giab_truth/` | `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |

### Data Sources

- **GIAB HG002**: [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/)
- **GRCh38 Reference**: [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40)
- **Simulated reads**: Generated with mason2 (Illumina) and pbsim3 (long reads)

## Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| freemap | 0.0.1+ | Primary aligner |
| Python | 3.8+ | Analysis scripts |
| numpy, matplotlib, scipy | latest | Python dependencies |
| intervaltree | latest | GIAB region lookup |
| samtools | 1.15+ | Depth file generation |
