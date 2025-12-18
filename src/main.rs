use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bgzf;
use rust_htslib::tpool::ThreadPool;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Read as IoRead, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::{fs::File, str};

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

    /// Optional FASTQ (plain, gzip, or bgzip) with BQ tag in header for barcode qualities (loads into memory)
    #[arg(long, value_name = "FASTQ")]
    fastq_bq: Option<PathBuf>,

    /// Optional cache file for --fastq-bq (loads if present, otherwise created)
    #[arg(long, value_name = "CACHE", requires = "fastq_bq")]
    fastq_bq_cache: Option<PathBuf>,

    /// Number of threads for BAM compression/decompression
    #[arg(short = 't', long, default_value = "4")]
    threads: usize,
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

/// Parsed barcode/UMI qualities from a BQ token.
#[derive(Debug, Clone)]
struct BqQuals {
    cb: Vec<u8>,
    umi: Option<Vec<u8>>,
}

const BQ_CACHE_MAGIC: &[u8; 8] = b"TBQMAP01";

struct FastqReader {
    reader: bgzf::Reader,
    _tpool: Option<ThreadPool>,
}

impl FastqReader {
    fn from_path(path: &Path, threads: usize) -> Result<Self> {
        let mut reader = bgzf::Reader::from_path(path)
            .with_context(|| format!("Failed to open FASTQ: {:?}", path))?;
        let tpool = if threads > 1 {
            let is_bgzip = bgzf::is_bgzip(path)
                .with_context(|| format!("Failed to detect bgzip FASTQ: {:?}", path))?;
            if is_bgzip {
                let thread_count =
                    u32::try_from(threads).context("FASTQ thread count exceeds u32")?;
                let tpool =
                    ThreadPool::new(thread_count).context("Failed to create FASTQ thread pool")?;
                reader
                    .set_thread_pool(&tpool)
                    .context("Failed to set FASTQ thread pool")?;
                Some(tpool)
            } else {
                None
            }
        } else {
            None
        };
        Ok(Self {
            reader,
            _tpool: tpool,
        })
    }
}

impl IoRead for FastqReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.reader.read(buf)
    }
}

/// Parse `|BQ:` token from a FASTQ header and return concatenated i7+i5+CBC qualities and UMI qualities if present.
fn parse_bq_token(header: &str) -> Option<BqQuals> {
    let bq_start = header.find("|BQ:")?;
    let token = &header[bq_start + 4..];
    let token = token.split_whitespace().next().unwrap_or(token);

    let mut i7 = None;
    let mut i5 = None;
    let mut cbc = None;
    let mut umi = None;

    for part in token.split(';') {
        if let Some((label, qual)) = part.split_once(':') {
            match label {
                "i7" => i7 = Some(qual.to_string()),
                "i5" => i5 = Some(qual.to_string()),
                "CBC" => cbc = Some(qual.to_string()),
                "UMI" => umi = Some(qual.to_string()),
                _ => {}
            }
        }
    }

    match (i7, i5, cbc) {
        (Some(a), Some(b), Some(c)) => Some(BqQuals {
            cb: [a, b, c].concat().into_bytes(),
            umi: umi.map(|u| u.into_bytes()),
        }),
        _ => None,
    }
}

/// Load FASTQ into memory and build a map: read name -> parsed barcode/UMI qualities.
fn load_bq_map(fastq_path: &Path, threads: usize) -> Result<HashMap<String, BqQuals>> {
    let mut reader = open_fastq_reader(fastq_path, threads)?.lines();
    let mut map = HashMap::new();

    loop {
        let header = match reader.next() {
            Some(Ok(h)) => h,
            Some(Err(e)) => return Err(e).context("Failed reading FASTQ header"),
            None => break,
        };

        // Skip seq, plus, qual
        let _ = reader.next();
        let _ = reader.next();
        let _ = reader.next();

        let name = header
            .strip_prefix('@')
            .unwrap_or(&header)
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        if let Some(quals) = parse_bq_token(&header) {
            map.insert(name, quals);
        }
    }

    Ok(map)
}

fn open_fastq_reader(fastq_path: &Path, threads: usize) -> Result<Box<dyn BufRead>> {
    let reader = FastqReader::from_path(fastq_path, threads)?;
    Ok(Box::new(BufReader::new(reader)))
}

fn load_bq_map_with_cache(
    fastq_path: &Path,
    cache_path: Option<&Path>,
    threads: usize,
) -> Result<HashMap<String, BqQuals>> {
    if let Some(cache_path) = cache_path {
        if cache_path.exists() {
            return read_bq_cache(cache_path)
                .with_context(|| format!("Failed to read BQ cache: {:?}", cache_path));
        }
    }

    let map = load_bq_map(fastq_path, threads)?;
    if let Some(cache_path) = cache_path {
        write_bq_cache(cache_path, &map)
            .with_context(|| format!("Failed to write BQ cache: {:?}", cache_path))?;
    }
    Ok(map)
}

fn read_bq_cache(cache_path: &Path) -> Result<HashMap<String, BqQuals>> {
    let file = File::open(cache_path)
        .with_context(|| format!("Failed to open BQ cache: {:?}", cache_path))?;
    let mut reader = BufReader::new(file);

    let mut magic = [0u8; 8];
    reader
        .read_exact(&mut magic)
        .context("Failed to read BQ cache header")?;
    if &magic != BQ_CACHE_MAGIC {
        anyhow::bail!("BQ cache has invalid header");
    }

    let count = read_u64(&mut reader).context("Failed to read BQ cache entry count")?;
    let count_usize = usize::try_from(count).context("BQ cache entry count exceeds usize")?;
    let mut map = HashMap::with_capacity(count_usize);

    for _ in 0..count_usize {
        let name_bytes = read_bytes(&mut reader).context("Failed to read BQ cache name")?;
        let cb = read_bytes(&mut reader).context("Failed to read BQ cache CB qualities")?;
        let umi_present = read_u8(&mut reader).context("Failed to read BQ cache UMI presence")?;
        let umi = if umi_present == 1 {
            Some(read_bytes(&mut reader).context("Failed to read BQ cache UMI qualities")?)
        } else {
            None
        };

        let name = String::from_utf8(name_bytes).context("BQ cache contains non-UTF8 read name")?;
        map.insert(name, BqQuals { cb, umi });
    }

    Ok(map)
}

fn write_bq_cache(cache_path: &Path, map: &HashMap<String, BqQuals>) -> Result<()> {
    let file = File::create(cache_path)
        .with_context(|| format!("Failed to create BQ cache: {:?}", cache_path))?;
    let mut writer = BufWriter::new(file);

    writer
        .write_all(BQ_CACHE_MAGIC)
        .context("Failed to write BQ cache header")?;
    let count = u64::try_from(map.len()).context("BQ cache entry count exceeds u64")?;
    write_u64(&mut writer, count).context("Failed to write BQ cache entry count")?;

    for (name, quals) in map {
        write_bytes(&mut writer, name.as_bytes()).context("Failed to write BQ cache name")?;
        write_bytes(&mut writer, &quals.cb).context("Failed to write BQ cache CB qualities")?;
        match &quals.umi {
            Some(umi) => {
                write_u8(&mut writer, 1).context("Failed to write BQ cache UMI presence")?;
                write_bytes(&mut writer, umi).context("Failed to write BQ cache UMI qualities")?;
            }
            None => {
                write_u8(&mut writer, 0).context("Failed to write BQ cache UMI presence")?;
            }
        }
    }

    writer.flush().context("Failed to flush BQ cache")?;
    Ok(())
}

fn read_u64<R: IoRead>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn write_u64<W: IoWrite>(writer: &mut W, value: u64) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn read_u8<R: IoRead>(reader: &mut R) -> Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn write_u8<W: IoWrite>(writer: &mut W, value: u8) -> Result<()> {
    writer.write_all(&[value])?;
    Ok(())
}

fn read_bytes<R: IoRead>(reader: &mut R) -> Result<Vec<u8>> {
    let len = read_u64(reader)?;
    let len_usize = usize::try_from(len).context("BQ cache length exceeds usize")?;
    let mut buf = vec![0u8; len_usize];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

fn write_bytes<W: IoWrite>(writer: &mut W, bytes: &[u8]) -> Result<()> {
    let len = u64::try_from(bytes.len()).context("BQ cache length exceeds u64")?;
    write_u64(writer, len)?;
    writer.write_all(bytes)?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Validate that either output or in_place is specified
    if cli.output.is_none() && !cli.in_place {
        anyhow::bail!("Either --output or --in-place must be specified");
    }

    let bq_map = if let Some(ref fastq) = cli.fastq_bq {
        Some(load_bq_map_with_cache(
            fastq,
            cli.fastq_bq_cache.as_deref(),
            cli.threads,
        )?)
    } else {
        None
    };

    let mut reader = bam::Reader::from_path(&cli.input)
        .with_context(|| format!("Failed to open input BAM: {:?}", cli.input))?;

    // Enable multi-threaded decompression
    reader.set_threads(cli.threads)?;

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

    // Enable multi-threaded compression
    writer.set_threads(cli.threads)?;

    let mut n_total: u64 = 0;
    let mut n_tagged: u64 = 0;
    let mut n_skipped: u64 = 0;

    for result in reader.records() {
        let mut record = result.context("Failed to read BAM record")?;
        n_total += 1;

        let qname = str::from_utf8(record.qname()).context("Read name is not valid UTF-8")?;

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
                    let (cell_barcode_qual, umi_qual) = if let Some(bq_map) = bq_map.as_ref() {
                        if let Some(quals) = bq_map.get(qname) {
                            let cbq = quals.cb.clone();
                            let umi_q = quals
                                .umi
                                .clone()
                                .unwrap_or_else(|| perfect_quality(components.umi.len()));
                            (cbq, umi_q)
                        } else {
                            (
                                perfect_quality(cell_barcode.len()),
                                perfect_quality(components.umi.len()),
                            )
                        }
                    } else {
                        (
                            perfect_quality(cell_barcode.len()),
                            perfect_quality(components.umi.len()),
                        )
                    };

                    // Add tags to BAM record
                    record.push_aux(b"CB", bam::record::Aux::String(&cell_barcode))?;
                    record.push_aux(
                        b"CY",
                        bam::record::Aux::String(str::from_utf8(&cell_barcode_qual).unwrap()),
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
    use tempfile::NamedTempFile;

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

    #[test]
    fn bq_cache_roundtrip() {
        let mut map = HashMap::new();
        map.insert(
            "read1".to_string(),
            BqQuals {
                cb: b"ABC".to_vec(),
                umi: Some(b"XYZ".to_vec()),
            },
        );
        map.insert(
            "read2".to_string(),
            BqQuals {
                cb: b"QQ".to_vec(),
                umi: None,
            },
        );

        let file = NamedTempFile::new().unwrap();
        write_bq_cache(file.path(), &map).unwrap();
        let loaded = read_bq_cache(file.path()).unwrap();

        assert_eq!(loaded.len(), map.len());
        assert_eq!(loaded.get("read1").unwrap().cb, b"ABC".to_vec());
        assert_eq!(loaded.get("read1").unwrap().umi.as_ref().unwrap(), b"XYZ");
        assert_eq!(loaded.get("read2").unwrap().cb, b"QQ".to_vec());
        assert!(loaded.get("read2").unwrap().umi.is_none());
    }
}
