# tagbam

[![CI](https://github.com/biobenkj/tagbam/actions/workflows/ci.yml/badge.svg)](https://github.com/biobenkj/tagbam/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/biobenkj/tagbam?display_name=tag&sort=semver)](https://github.com/biobenkj/tagbam/releases)

### Re-tag BAM files by parsing cell barcodes and UMIs from read names

`tagbam` extracts cell barcodes and UMIs from read names and adds them as standard BAM tags for single-cell analysis pipelines.

Reads from read names in format: `{uuid}_{i7}-{i5}-{CBC}_{UMI}`

Creates BAM tags:
- `CB:Z` - Cell barcode (concatenated i7+i5+CBC)
- `CY:Z` - Cell barcode quality (from FASTQ `|BQ:` header if provided, otherwise perfect quality)
- `UB:Z` - UMI sequence
- `UY:Z` - UMI quality (from FASTQ `|BQ:` header if provided, otherwise perfect quality)

## Install

### From source
```bash
cargo install --git https://github.com/biobenkj/tagbam
```

### Download binaries
Grab macOS (universal) and Linux builds from [Releases](https://github.com/biobenkj/tagbam/releases).

## Usage

### Basic usage with output file

```bash
tagbam --input input.bam --output tagged.bam
```

### In-place tagging (modifies input file directly)

```bash
tagbam --input input.bam --in-place
```

### Skip unparseable read names

```bash
tagbam --input input.bam --output tagged.bam --skip-unparseable
```

By default, the tool errors on unparseable read names. Use `--skip-unparseable` to continue processing and skip those reads with a warning.

### Multi-threaded compression/decompression

```bash
tagbam --input input.bam --output tagged.bam --threads 8
```

Control the number of threads used for BAM compression and decompression (default: 4). Increasing threads can significantly improve performance on large files. A good starting point is to match your CPU core count.

## Read name format

Read names must follow this format: `{uuid}_{i7}-{i5}-{CBC}_{UMI}`

### Example

```
2efc6b85-aa0d-4c1d-ab33-bf5f442fe47c_TTGGCTCC-GGTCGGCG-ACTTGA_GAAGCAGT
```

This produces:
- `CB:Z:TTGGCTCCGGTCGGCGACTTGA` (22 bases: i7+i5+CBC concatenated)
- `CY:Z:IIIIIIIIIIIIIIIIIIIIII` (22 I's for perfect quality Q40 if no `--fastq-bq`)
- `UB:Z:GAAGCAGT` (8 bases)
- `UY:Z:IIIIIIII` (8 I's for perfect quality Q40 if no `--fastq-bq`)

### Supplying barcode/UMI qualities from FASTQ (`--fastq-bq`)

If your FASTQ headers include a `|BQ:` token (e.g., `|BQ:i7:<qual>;i5:<qual>;CBC:<qual>;UMI:<qual>`), you can supply that FASTQ to reuse the barcode/UMI qualities when tagging the BAM:

```bash
tagbam --input input.bam --output tagged.bam --fastq-bq demuxed.fastq
```

- `CY` is populated from concatenated i7+i5+CBC qualities in the `|BQ:` token.
- `UY` is populated from the `UMI` quality in the `|BQ:` token if present; otherwise it falls back to perfect quality.
- Reads without a `|BQ:` token (or absent in the FASTQ map) still receive perfect-quality tags.

## Behavior

- **Existing tags**: Reads that already have any of the four tags (`CB`, `CY`, `UB`, `UY`) are skipped with a warning. Existing tags are preserved.
- **Invalid read names**: By default, the tool exits with an error. Use `--skip-unparseable` to skip these reads and continue.
- **Input/Output**: Either `--output` or `--in-place` must be specified (they are mutually exclusive).

## MSRV

- Minimum Supported Rust Version: **1.83**

## License

Dual-licensed under **MIT** or **Apache-2.0**.
You may choose either license.
See [`LICENSE-MIT`](LICENSE-MIT) and [`LICENSE-APACHE`](LICENSE-APACHE) for details.

## Security

If you discover a security vulnerability in `tagbam`, please **do not** open a public issue.
Instead, email **Ben Johnson** at the address listed on the GitHub profile for [@biobenkj](https://github.com/biobenkj).
We will coordinate a fix and disclosure.
