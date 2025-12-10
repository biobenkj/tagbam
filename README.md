# tagbam

[![CI](https://github.com/biobenkj/tagbam/actions/workflows/ci.yml/badge.svg)](https://github.com/biobenkj/tagbam/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/biobenkj/tagbam?display_name=tag&sort=semver)](https://github.com/biobenkj/tagbam/releases)

### Re-tag BAM files by parsing cell barcodes and UMIs from read names

`tagbam` extracts cell barcodes and UMIs from read names and adds them as standard BAM tags for single-cell analysis pipelines.

Reads from read names in format: `{uuid}_{i7}-{i5}-{CBC}_{UMI}`

Creates BAM tags:
- `CB:Z` - Cell barcode (concatenated i7+i5+CBC)
- `CY:Z` - Cell barcode quality (perfect quality: all 'I')
- `UB:Z` - UMI sequence
- `UY:Z` - UMI quality (perfect quality: all 'I')

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

## Read name format

Read names must follow this format: `{uuid}_{i7}-{i5}-{CBC}_{UMI}`

### Example

```
2efc6b85-aa0d-4c1d-ab33-bf5f442fe47c_TTGGCTCC-GGTCGGCG-ACTTGA_GAAGCAGT
```

This produces:
- `CB:Z:TTGGCTCCGGTCGGCGACTTGA` (22 bases: i7+i5+CBC concatenated)
- `CY:Z:IIIIIIIIIIIIIIIIIIIIII` (22 I's for perfect quality Q40)
- `UB:Z:GAAGCAGT` (8 bases)
- `UY:Z:IIIIIIII` (8 I's for perfect quality Q40)

## Behavior

- **Existing tags**: Reads that already have any of the four tags (`CB`, `CY`, `UB`, `UY`) are skipped with a warning. Existing tags are preserved.
- **Invalid read names**: By default, the tool exits with an error. Use `--skip-unparseable` to skip these reads and continue.
- **Input/Output**: Either `--output` or `--in-place` must be specified (they are mutually exclusive).

## MSRV

- Minimum Supported Rust Version: **1.82**

## License

Dual-licensed under **MIT** or **Apache-2.0**.
You may choose either license.
See [`LICENSE-MIT`](LICENSE-MIT) and [`LICENSE-APACHE`](LICENSE-APACHE) for details.

## Security

If you discover a security vulnerability in `tagbam`, please **do not** open a public issue.
Instead, email **Ben Johnson** at the address listed on the GitHub profile for [@biobenkj](https://github.com/biobenkj).
We will coordinate a fix and disclosure.
