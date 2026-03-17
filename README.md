<p align="center">
  <img src="logo.svg?v=6" alt="freemap logo" width="200"/>
</p>

# freemap

[![CI](https://github.com/genpat-it/freemap/actions/workflows/ci.yml/badge.svg)](https://github.com/genpat-it/freemap/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![GitHub release](https://img.shields.io/github/v/release/genpat-it/freemap)](https://github.com/genpat-it/freemap/releases)
[![Docker](https://img.shields.io/badge/Docker-ghcr.io%2Fgenpat--it%2Ffreemap-blue?logo=docker)](https://github.com/genpat-it/freemap/pkgs/container/freemap)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18671716.svg)](https://doi.org/10.5281/zenodo.18671716)

A fast trajectory-based sequence aligner for long reads (ONT, PacBio) and short reads (Illumina). freemap uses a novel DP-free CIGAR generation approach that adds only ~2-3% overhead versus ~200%+ for traditional dynamic programming methods.

## Installation

```bash
cargo build --release
```

For maximum performance on your CPU:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

The binary is at `target/release/freemap`.

## Usage

```bash
# ONT long reads
freemap -x map-ont -a ref.fa reads.fq out.sam

# PacBio HiFi
freemap -x map-hifi -a ref.fa reads.fq out.sam

# PacBio CLR
freemap -x map-pb -a ref.fa reads.fq out.sam

# Illumina short reads (single-end)
freemap -x sr -a ref.fa reads.fq out.sam

# Illumina short reads (paired-end)
freemap -x sr -a -1 R1.fq -2 R2.fq ref.fa out.sam

# PAF output (default)
freemap -x map-ont ref.fa reads.fq out.paf
```

### Presets

| Preset | Technology | Description |
|--------|-----------|-------------|
| `map-ont` | ONT | k=15, w=10, trajectory CIGAR |
| `map-hifi` | PacBio HiFi | k=19, w=19, trajectory CIGAR |
| `map-pb` | PacBio CLR | k=19, w=10, trajectory CIGAR |
| `sr` | Illumina | k=21, w=11, polish mode |

### Paired-end options

| Flag | Description |
|------|-------------|
| `-1 FILE` | First read file (R1) |
| `-2 FILE` | Second read file (R2) |
| `-I MIN:MAX` | Expected insert size range [0:1000] |

### Indexing and algorithm

| Flag | Description |
|------|-------------|
| `-k INT` | k-mer size [19] |
| `-w INT` | Minimizer window size [19] |
| `-f INT` | Max k-mer frequency [200] |
| `-L INT` | Chaining lookback limit [16] |
| `-G INT` | Max gap difference for penalty [50] |
| `-S INT` | Gap penalty scaling factor [5] |
| `-t INT` | Threads [all] |

### Alignment modes

| Flag | Description |
|------|-------------|
| `-g` | Trajectory mode: CIGAR from geometry (no DP) |
| `-p, --polish` | Polish mode: base-level indel detection (DP-free) |
| `-H` | Homopolymer compression (recommended for ONT) |
| `-r` | Refine boundaries with micro-anchors |
| `-R` | Short-read mode |
| `-u` | Ultralong mode: relaxed chaining for ONT ultralong reads |

### Output

| Flag | Description |
|------|-------------|
| `-a` | SAM output (default: PAF) |
| `-c` | Generate detailed CIGAR (heuristic gap alignment) |
| `--multi` | Output secondary and supplementary alignments |
| `--max-secondary N` | Max secondary alignments per read [5] |
| `-q` | Quiet mode |

### Index I/O

| Flag | Description |
|------|-------------|
| `-d FILE` | Save pre-built index to disk |
| `-i FILE` | Load pre-built index from disk |

### Advanced

| Flag | Description |
|------|-------------|
| `--ransac-threshold FLOAT` | RANSAC inlier threshold for trajectory regression [25.0] |
| `--polish-max-indel INT` | Max indel size for CIGAR polishing [4] |
| `--band INT` | Custom chaining band width |
| `--lookback INT` | Override chaining lookback limit |
| `--no-tiebreaker` | Disable primary-chromosome tiebreaker |

Run `freemap -h` for the full flag reference.

## Quick start with test data

To verify your installation, you can run freemap on publicly available *E. coli* K-12 data:

```bash
# Download E. coli K-12 MG1655 reference
wget -O ecoli.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
gunzip ecoli.fa.gz

# Simulate 1000 HiFi reads with pbsim3 (optional — or use your own reads)
pbsim --strategy wgs --method qshmm --qshmm QSHMM-RSII.model \
      --depth 5 --genome ecoli.fa --prefix ecoli_test

# Build index and align
freemap -x map-hifi -a ecoli.fa ecoli_test_0001.fastq out.sam
```

**Expected output:** `out.sam` — a standard SAM file with header lines (`@SQ`, `@PG`) followed by one alignment record per read, including trajectory-based CIGAR strings. You can verify with:

```bash
# Check mapped reads
samtools flagstat out.sam

# View coverage
samtools sort out.sam -o out.bam && samtools index out.bam
samtools depth -a out.bam | head
```

### Test datasets used in the paper

| Dataset | Source | Accession / URL |
|---------|--------|-----------------|
| *E. coli* K-12 MG1655 | NCBI | [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3) |
| Human GRCh38 | NCBI | [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) |
| GIAB HG002 (HiFi, ONT) | GIAB | [HG002 data](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/) |
| GIAB HC regions v4.2.1 | GIAB | [v4.2.1 benchmark](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/) |

Simulated reads were generated with [pbsim3](https://github.com/yukiteruono/pbsim3) v3.0.0 using fixed seeds for reproducibility. Generation scripts and parameters are documented in `benchmark/scripts/`.

## Reproducing paper results

```bash
# Reproduce all tables and figures
bash benchmark/scripts/reproduce_paper.sh /path/to/data

# Reproduce a single section
bash benchmark/scripts/reproduce_paper.sh /path/to/data table1
```

See [benchmark/README.md](benchmark/README.md) for data layout and requirements.

## Citation

If you use freemap in your research, please cite:

> de Ruvo A., Radomski N., Flammini M., Di Pasquale A. (2026). freemap: DP-free CIGAR generation for long reads via trajectory inference. *Bioinformatics* (submitted).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18671716.svg)](https://doi.org/10.5281/zenodo.18671716)

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

[MIT](LICENSE)
