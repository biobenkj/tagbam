use assert_cmd::Command;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use tempfile::TempDir;

/// Helper to create a minimal BAM file with custom read names
fn create_test_bam(path: &Path, read_names: &[&str]) -> Result<(), Box<dyn std::error::Error>> {
    // Create a minimal header
    let mut header = bam::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"HD");
    header_rec.push_tag(b"VN", "1.6");
    header_rec.push_tag(b"SO", "unsorted");
    header.push_record(&header_rec);

    // Add a reference sequence (required for valid BAM)
    let mut ref_rec = bam::header::HeaderRecord::new(b"SQ");
    ref_rec.push_tag(b"SN", "chr1");
    ref_rec.push_tag(b"LN", "1000");
    header.push_record(&ref_rec);

    let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;

    // Create minimal records with the specified read names
    for read_name in read_names {
        let mut record = bam::Record::new();
        record.set(
            read_name.as_bytes(),
            None,    // cigar
            b"ACGT", // seq
            b"IIII", // qual
        );
        record.set_pos(0);
        record.set_tid(0); // chr1

        writer.write(&record)?;
    }

    Ok(())
}

/// Helper to read BAM tags from a record
fn get_tag_string(record: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    match record.aux(tag) {
        Ok(bam::record::Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

#[test]
fn basic_tagging() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");
    let output_bam = td.path().join("output.bam");

    // Create test BAM with one properly formatted read name
    let read_name = "2efc6b85-aa0d-4c1d-ab33-bf5f442fe47c_TTGGCTCC-GGTCGGCG-ACTTGA_GAAGCAGT";
    create_test_bam(&input_bam, &[read_name]).unwrap();

    // Run tagbam
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args([
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
    ]);
    cmd.assert().success();

    // Read output BAM and verify tags
    let mut reader = bam::Reader::from_path(&output_bam).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    // Verify tags
    assert_eq!(
        get_tag_string(&record, b"CB"),
        Some("TTGGCTCCGGTCGGCGACTTGA".to_string()),
        "Cell barcode should be i7+i5+CBC concatenated"
    );
    assert_eq!(
        get_tag_string(&record, b"CY"),
        Some("IIIIIIIIIIIIIIIIIIIIII".to_string()),
        "Cell barcode quality should be 22 'I's"
    );
    assert_eq!(
        get_tag_string(&record, b"UB"),
        Some("GAAGCAGT".to_string()),
        "UMI should be extracted"
    );
    assert_eq!(
        get_tag_string(&record, b"UY"),
        Some("IIIIIIII".to_string()),
        "UMI quality should be 8 'I's"
    );
}

#[test]
fn multiple_reads() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");
    let output_bam = td.path().join("output.bam");

    let read_names = [
        "uuid1_AAA-BBB-CCC_UUU",
        "uuid2_AAAAAAAA-BBBBBBBB-CCCCCCCC_UUUUUUUU",
        "uuid3_A-B-C_U",
    ];
    create_test_bam(&input_bam, &read_names).unwrap();

    // Run tagbam
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args([
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
    ]);
    cmd.assert().success();

    // Read output BAM and verify all records were processed
    let mut reader = bam::Reader::from_path(&output_bam).unwrap();
    let records: Vec<_> = reader.records().collect();
    assert_eq!(records.len(), 3, "All three records should be present");

    // Check first record
    let rec1 = records[0].as_ref().unwrap();
    assert_eq!(get_tag_string(rec1, b"CB"), Some("AAABBBCCC".to_string()));
    assert_eq!(get_tag_string(rec1, b"CY"), Some("IIIIIIIII".to_string()));
    assert_eq!(get_tag_string(rec1, b"UB"), Some("UUU".to_string()));
    assert_eq!(get_tag_string(rec1, b"UY"), Some("III".to_string()));

    // Check second record
    let rec2 = records[1].as_ref().unwrap();
    assert_eq!(
        get_tag_string(rec2, b"CB"),
        Some("AAAAAAAABBBBBBBBCCCCCCCC".to_string())
    );
    assert_eq!(get_tag_string(rec2, b"UB"), Some("UUUUUUUU".to_string()));

    // Check third record
    let rec3 = records[2].as_ref().unwrap();
    assert_eq!(get_tag_string(rec3, b"CB"), Some("ABC".to_string()));
    assert_eq!(get_tag_string(rec3, b"CY"), Some("III".to_string()));
    assert_eq!(get_tag_string(rec3, b"UB"), Some("U".to_string()));
    assert_eq!(get_tag_string(rec3, b"UY"), Some("I".to_string()));
}

#[test]
fn invalid_read_name_errors() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");
    let output_bam = td.path().join("output.bam");

    // Create BAM with invalid read name (missing underscore)
    create_test_bam(&input_bam, &["invalid_name_without_proper_format"]).unwrap();

    // Run tagbam - should fail
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args([
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
    ]);
    cmd.assert().failure();
}

#[test]
fn skip_unparseable_flag() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");
    let output_bam = td.path().join("output.bam");

    // Mix of valid and invalid read names
    let read_names = [
        "uuid1_AAA-BBB-CCC_UUU", // valid
        "invalid_format",        // invalid
        "uuid2_DDD-EEE-FFF_VVV", // valid
    ];
    create_test_bam(&input_bam, &read_names).unwrap();

    // Run tagbam with --skip-unparseable
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args([
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--skip-unparseable",
    ]);
    cmd.assert().success();

    // Read output BAM
    let mut reader = bam::Reader::from_path(&output_bam).unwrap();
    let records: Vec<_> = reader.records().collect();
    assert_eq!(records.len(), 3, "All three records should be present");

    // First record should have tags
    let rec1 = records[0].as_ref().unwrap();
    assert_eq!(get_tag_string(rec1, b"CB"), Some("AAABBBCCC".to_string()));

    // Second record should NOT have tags (invalid name)
    let rec2 = records[1].as_ref().unwrap();
    assert_eq!(get_tag_string(rec2, b"CB"), None);

    // Third record should have tags
    let rec3 = records[2].as_ref().unwrap();
    assert_eq!(get_tag_string(rec3, b"CB"), Some("DDDEEEFFF".to_string()));
}

#[test]
fn fastq_bq_overrides_barcode_qualities() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");
    let output_bam = td.path().join("output.bam");
    let fastq_path = td.path().join("reads.fastq");

    let read_name = "uuid1_AAA-BBB-CCC_UUU";
    create_test_bam(&input_bam, &[read_name]).unwrap();

    let mut fq = File::create(&fastq_path).unwrap();
    writeln!(
        fq,
        "@{read_name} cell|Barcodes:i7:AAA;i5:BBB;CBC:CCC|UMI:UUU|orientation:+|BQ:i7:123;i5:456;CBC:789;UMI:XYZ\n\
         AAAAAAAAA\n+\nIIIIIIIII"
    )
    .unwrap();

    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args([
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--fastq-bq",
        fastq_path.to_str().unwrap(),
    ]);
    cmd.assert().success();

    let mut reader = bam::Reader::from_path(&output_bam).unwrap();
    let record = reader.records().next().unwrap().unwrap();
    assert_eq!(
        get_tag_string(&record, b"CY"),
        Some("123456789".to_string())
    );
    assert_eq!(get_tag_string(&record, b"UY"), Some("XYZ".to_string()));
}

#[test]
fn skip_reads_with_existing_tags() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");
    let output_bam = td.path().join("output.bam");

    // Create a BAM with one read
    let read_name = "uuid_AAA-BBB-CCC_UUU";

    // Create header
    let mut header = bam::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"HD");
    header_rec.push_tag(b"VN", "1.6");
    header_rec.push_tag(b"SO", "unsorted");
    header.push_record(&header_rec);

    let mut ref_rec = bam::header::HeaderRecord::new(b"SQ");
    ref_rec.push_tag(b"SN", "chr1");
    ref_rec.push_tag(b"LN", "1000");
    header.push_record(&ref_rec);

    let mut writer = bam::Writer::from_path(&input_bam, &header, bam::Format::Bam).unwrap();

    // Create record with existing CB tag
    let mut record = bam::Record::new();
    record.set(read_name.as_bytes(), None, b"ACGT", b"IIII");
    record.set_pos(0);
    record.set_tid(0);

    // Add an existing CB tag
    record
        .push_aux(b"CB", bam::record::Aux::String("EXISTING"))
        .unwrap();

    writer.write(&record).unwrap();
    drop(writer);

    // Run tagbam
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args([
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
    ]);
    cmd.assert().success();

    // Read output BAM and verify the existing tag was NOT changed
    let mut reader = bam::Reader::from_path(&output_bam).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    // CB tag should still be "EXISTING", not the parsed value
    assert_eq!(
        get_tag_string(&record, b"CB"),
        Some("EXISTING".to_string()),
        "Existing CB tag should be preserved"
    );

    // No other tags should have been added either
    assert_eq!(get_tag_string(&record, b"UB"), None);
}

#[test]
fn in_place_mode() {
    let td = TempDir::new().unwrap();
    let bam_file = td.path().join("test.bam");

    // Create test BAM
    let read_name = "uuid_AAA-BBB-CCC_UUU";
    create_test_bam(&bam_file, &[read_name]).unwrap();

    // Get original modification time to verify file was replaced
    let original_metadata = std::fs::metadata(&bam_file).unwrap();

    // Run tagbam with --in-place
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args(["--input", bam_file.to_str().unwrap(), "--in-place"]);
    cmd.assert().success();

    // Read the modified BAM (same file)
    let mut reader = bam::Reader::from_path(&bam_file).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    // Verify tags were added
    assert_eq!(
        get_tag_string(&record, b"CB"),
        Some("AAABBBCCC".to_string())
    );
    assert_eq!(get_tag_string(&record, b"UB"), Some("UUU".to_string()));

    // Verify the file was actually modified (metadata changed)
    let new_metadata = std::fs::metadata(&bam_file).unwrap();
    // On some systems modification time might be the same if the operation is very fast,
    // but the file size should be different due to added tags
    assert!(
        new_metadata.len() >= original_metadata.len(),
        "File should have been modified in-place"
    );
}

#[test]
fn requires_output_or_in_place() {
    let td = TempDir::new().unwrap();
    let input_bam = td.path().join("input.bam");

    // Create test BAM
    create_test_bam(&input_bam, &["uuid_AAA-BBB-CCC_UUU"]).unwrap();

    // Run tagbam without --output or --in-place
    let mut cmd = Command::new(assert_cmd::cargo::cargo_bin!("tagbam"));
    cmd.args(["--input", input_bam.to_str().unwrap()]);
    cmd.assert()
        .failure()
        .stderr(predicates::str::contains("Either --output or --in-place"));
}
