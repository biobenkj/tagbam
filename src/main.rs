use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "tagbam",
    version,
    about = "Re-tag BAM files by parsing cell barcodes and UMIs from read names",
    long_about = "Parses read names in format {uuid}_{i7}-{i5}-{CBC}_{UMI} and adds BAM tags:\n\
                  - CB:Z (cell barcode: i7+i5+CBC concatenated)\n\
                  - CY:Z (cell barcode quality: all 'I' for perfect quality)\n\
                  - UB:Z (UMI sequence)\n\
                  - UY:Z (UMI quality: all 'I' for perfect quality)"
)]
struct Cli {
    /// Input BAM file
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,

    /// Output BAM file (required unless --in-place is used)
    #[arg(short, long, value_name = "FILE", conflicts_with = "in_place")]
    output: Option<PathBuf>,

    /// Modify the input BAM file in-place
    #[arg(long, conflicts_with = "output")]
    in_place: bool,

    /// Skip reads with unparseable names instead of erroring
    #[arg(long)]
    skip_unparseable: bool,
}

/// Parsed components from read name: {uuid}_{i7}-{i5}-{CBC}_{UMI}
#[derive(Debug, PartialEq)]
struct ReadNameComponents {
    i7: String,
    i5: String,
    cbc: String,
    umi: String,
}

/// Parse read name in format: {uuid}_{i7}-{i5}-{CBC}_{UMI}
///
/// Example: 2efc6b85-aa0d-4c1d-ab33-bf5f442fe47c_TTGGCTCC-GGTCGGCG-ACTTGA_GAAGCAGT
/// Returns: ReadNameComponents { i7: "TTGGCTCC", i5: "GGTCGGCG", cbc: "ACTTGA", umi: "GAAGCAGT" }
fn parse_read_name(name: &str) -> Result<ReadNameComponents> {
    let parts: Vec<&str> = name.split('_').collect();

    if parts.len() != 3 {
        anyhow::bail!(
            "Expected 3 underscore-separated parts in read name, found {}: '{}'",
            parts.len(),
            name
        );
    }

    // parts[0] = uuid (ignored)
    // parts[1] = i7-i5-CBC
    // parts[2] = UMI

    let barcode_parts: Vec<&str> = parts[1].split('-').collect();
    if barcode_parts.len() != 3 {
        anyhow::bail!(
            "Expected 3 hyphen-separated barcode parts, found {}: '{}'",
            barcode_parts.len(),
            parts[1]
        );
    }

    Ok(ReadNameComponents {
        i7: barcode_parts[0].to_string(),
        i5: barcode_parts[1].to_string(),
        cbc: barcode_parts[2].to_string(),
        umi: parts[2].to_string(),
    })
}

/// Create perfect quality string: 'I' repeated n times (Phred Q40, ASCII 73)
fn perfect_quality(length: usize) -> Vec<u8> {
    vec![b'I'; length]
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Validate that either output or in_place is specified
    if cli.output.is_none() && !cli.in_place {
        anyhow::bail!("Either --output or --in-place must be specified");
    }

    let mut reader = bam::Reader::from_path(&cli.input)
        .with_context(|| format!("Failed to open input BAM: {:?}", cli.input))?;

    let header = bam::Header::from_template(reader.header());

    // Determine output path: either specified output, or a temp file for in-place mode
    let output_path = if let Some(ref out) = cli.output {
        out.clone()
    } else {
        // In-place mode: create a temp file in the same directory
        let input_dir = cli
            .input
            .parent()
            .ok_or_else(|| anyhow::anyhow!("Cannot determine parent directory of input file"))?;
        input_dir.join(format!(
            ".{}.tmp",
            cli.input
                .file_name()
                .ok_or_else(|| anyhow::anyhow!("Cannot determine input file name"))?
                .to_string_lossy()
        ))
    };

    let mut writer = bam::Writer::from_path(&output_path, &header, bam::Format::Bam)
        .with_context(|| format!("Failed to create output BAM: {:?}", output_path))?;

    let mut n_total: u64 = 0;
    let mut n_tagged: u64 = 0;
    let mut n_skipped: u64 = 0;

    for result in reader.records() {
        let mut record = result.context("Failed to read BAM record")?;
        n_total += 1;

        let qname = std::str::from_utf8(record.qname()).context("Read name is not valid UTF-8")?;

        match parse_read_name(qname) {
            Ok(components) => {
                // Check if any of our tags already exist
                let has_existing_tags = record.aux(b"CB").is_ok()
                    || record.aux(b"CY").is_ok()
                    || record.aux(b"UB").is_ok()
                    || record.aux(b"UY").is_ok();

                if has_existing_tags {
                    eprintln!(
                        "Warning: Read '{}' already has CB/CY/UB/UY tags, skipping",
                        qname
                    );
                    n_skipped += 1;
                } else {
                    // Concatenate cell barcode: i7 + i5 + CBC
                    let cell_barcode =
                        format!("{}{}{}", components.i7, components.i5, components.cbc);
                    let cell_barcode_qual = perfect_quality(cell_barcode.len());

                    let umi_qual = perfect_quality(components.umi.len());

                    // Add tags to BAM record
                    record.push_aux(b"CB", bam::record::Aux::String(&cell_barcode))?;
                    record.push_aux(
                        b"CY",
                        bam::record::Aux::String(std::str::from_utf8(&cell_barcode_qual).unwrap()),
                    )?;
                    record.push_aux(b"UB", bam::record::Aux::String(&components.umi))?;
                    record.push_aux(
                        b"UY",
                        bam::record::Aux::String(std::str::from_utf8(&umi_qual).unwrap()),
                    )?;

                    n_tagged += 1;
                }
            }
            Err(e) => {
                if cli.skip_unparseable {
                    eprintln!("Warning: Skipping unparseable read name '{}': {}", qname, e);
                    n_skipped += 1;
                } else {
                    return Err(e).context(format!("Failed to parse read name '{}'", qname));
                }
            }
        }

        writer
            .write(&record)
            .context("Failed to write BAM record")?;
    }

    // Ensure writer is flushed and closed before moving the file
    drop(writer);

    // If in-place mode, replace the original file with the temp file
    if cli.in_place {
        std::fs::rename(&output_path, &cli.input)
            .with_context(|| "Failed to replace input file with tagged version".to_string())?;
        eprintln!(
            "In-place tagging complete: {} reads processed, {} tagged, {} skipped",
            n_total, n_tagged, n_skipped
        );
    } else {
        eprintln!(
            "Processed {} reads: {} tagged, {} skipped",
            n_total, n_tagged, n_skipped
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_valid_read_name() {
        let name = "2efc6b85-aa0d-4c1d-ab33-bf5f442fe47c_TTGGCTCC-GGTCGGCG-ACTTGA_GAAGCAGT";
        let result = parse_read_name(name).unwrap();

        assert_eq!(
            result,
            ReadNameComponents {
                i7: "TTGGCTCC".to_string(),
                i5: "GGTCGGCG".to_string(),
                cbc: "ACTTGA".to_string(),
                umi: "GAAGCAGT".to_string(),
            }
        );
    }

    #[test]
    fn parse_different_lengths() {
        let name = "uuid_AAA-BB-CCCCCC_UUUU";
        let result = parse_read_name(name).unwrap();

        assert_eq!(result.i7, "AAA");
        assert_eq!(result.i5, "BB");
        assert_eq!(result.cbc, "CCCCCC");
        assert_eq!(result.umi, "UUUU");
    }

    #[test]
    fn parse_missing_underscore() {
        let name = "uuid_TTGGCTCC-GGTCGGCG-ACTTGAGAAGCAGT"; // Missing underscore before UMI
        assert!(parse_read_name(name).is_err());
    }

    #[test]
    fn parse_missing_hyphen() {
        let name = "uuid_TTGGCTCC-GGTCGGCGACTTGA_GAAGCAGT"; // Missing hyphen in barcodes
        assert!(parse_read_name(name).is_err());
    }

    #[test]
    fn perfect_quality_length() {
        let qual = perfect_quality(8);
        assert_eq!(qual.len(), 8);
        assert!(qual.iter().all(|&b| b == b'I'));
    }

    #[test]
    fn perfect_quality_ascii() {
        let qual = perfect_quality(5);
        assert_eq!(std::str::from_utf8(&qual).unwrap(), "IIIII");
    }
}
