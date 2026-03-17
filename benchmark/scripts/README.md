# Benchmark Scripts

All scripts for reproducing the paper's tables and figures.

## Quick Start

```bash
# Reproduce everything
bash reproduce_paper.sh /path/to/data

# Reproduce one section
bash reproduce_paper.sh /path/to/data table1
```

See `../README.md` for full documentation, available sections, and expected data layout.

## Structure

| Path | Purpose |
|------|---------|
| `reproduce_paper.sh` | Single entry point |
| `config.sh` | Shared paths and logging |
| `evaluation/` | Core analysis (accuracy, concordance, CIGAR) |
| `figures/` | Figures and table generation |
| `data_prep/` | Simulated data preparation |

## Correctness Metrics

| Metric | Description |
|--------|-------------|
| `position` (default) | Correct if mapped within 500 bp of truth |
| `overlap` | Correct if overlap >= 10% of max(truth_len, aln_len) |

Primary alignments only -- secondary (0x100) and supplementary (0x800) are filtered.
