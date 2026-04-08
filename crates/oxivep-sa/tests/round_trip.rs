//! End-to-end round-trip tests for SaWriter -> SaReader.

use oxivep_cache::annotation::AnnotationProvider;
use oxivep_sa::common::AnnotationRecord;
use oxivep_sa::index::IndexHeader;
use oxivep_sa::reader::SaReader;
use oxivep_sa::writer::SaWriter;

#[test]
fn test_writer_reader_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("test_clinvar");

    let header = IndexHeader {
        schema_version: oxivep_sa::common::SCHEMA_VERSION,
        json_key: "clinvar".into(),
        name: "ClinVar".into(),
        version: "2024-12".into(),
        description: "Test ClinVar data".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    // Create records sorted by (chrom_idx, position)
    let chrom_map = vec!["chr1".to_string(), "chr2".to_string()];
    let records: Vec<AnnotationRecord> = (0..1000)
        .map(|i| AnnotationRecord {
            chrom_idx: if i < 700 { 0 } else { 1 },
            position: if i < 700 { 1000 + i } else { 1000 + (i - 700) },
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: format!(r#"{{"sig":"pathogenic","id":{}}}"#, i),
        })
        .collect();

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();

    // Verify files exist
    assert!(base.with_extension("osa").exists());
    assert!(base.with_extension("osa.idx").exists());

    // Open reader
    let reader = SaReader::open(&base.with_extension("osa")).unwrap();
    assert_eq!(reader.name(), "ClinVar");
    assert_eq!(reader.json_key(), "clinvar");

    // Query a known record on chr1
    let result = reader
        .annotate_position("chr1", 1050, "A", "G")
        .unwrap();
    assert!(result.is_some());
    let json = match result.unwrap() {
        oxivep_cache::annotation::AnnotationValue::Json(j) => j,
        other => panic!("Expected Json, got {:?}", other),
    };
    assert!(json.contains(r#""sig":"pathogenic""#));
    assert!(json.contains(r#""id":50"#));

    // Query a known record on chr2
    let result = reader
        .annotate_position("chr2", 1100, "A", "G")
        .unwrap();
    assert!(result.is_some());

    // Query with wrong allele (should not match when match_by_allele=true)
    let result = reader
        .annotate_position("chr1", 1050, "A", "T")
        .unwrap();
    assert!(result.is_none());

    // Query non-existent position
    let result = reader
        .annotate_position("chr1", 99999, "A", "G")
        .unwrap();
    assert!(result.is_none());

    // Query non-existent chromosome
    let result = reader
        .annotate_position("chr99", 1000, "A", "G")
        .unwrap();
    assert!(result.is_none());
}

#[test]
fn test_preload_then_query() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("test_gnomad");

    let header = IndexHeader {
        schema_version: oxivep_sa::common::SCHEMA_VERSION,
        json_key: "gnomad".into(),
        name: "gnomAD".into(),
        version: "4.0".into(),
        description: "Test gnomAD data".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    let chrom_map = vec!["chr1".to_string()];
    let records: Vec<AnnotationRecord> = (0..500)
        .map(|i| AnnotationRecord {
            chrom_idx: 0,
            position: 10000 + i * 10,
            ref_allele: "C".into(),
            alt_allele: "T".into(),
            json: format!(r#"{{"af":{:.6}}}"#, i as f64 / 10000.0),
        })
        .collect();

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();

    let reader = SaReader::open(&base.with_extension("osa")).unwrap();

    // Preload a range of positions
    let positions: Vec<u64> = (0..50).map(|i| 10000 + i * 10).collect();
    reader.preload("chr1", &positions).unwrap();

    // Query preloaded positions
    let result = reader
        .annotate_position("chr1", 10100, "C", "T")
        .unwrap();
    assert!(result.is_some());
}

#[test]
fn test_positional_annotation() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("test_phylop");

    let header = IndexHeader {
        schema_version: oxivep_sa::common::SCHEMA_VERSION,
        json_key: "phylopScore".into(),
        name: "PhyloP".into(),
        version: "1.0".into(),
        description: "PhyloP conservation scores".into(),
        assembly: "GRCh38".into(),
        match_by_allele: false,
        is_array: false,
        is_positional: true,
    };

    let chrom_map = vec!["chr1".to_string()];
    let records: Vec<AnnotationRecord> = (0..100)
        .map(|i| AnnotationRecord {
            chrom_idx: 0,
            position: 1000 + i,
            ref_allele: String::new(),
            alt_allele: String::new(),
            json: format!("{:.3}", (i as f64 - 50.0) / 10.0),
        })
        .collect();

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();

    let reader = SaReader::open(&base.with_extension("osa")).unwrap();

    // Positional: alleles don't matter
    let result = reader
        .annotate_position("chr1", 1025, "A", "G")
        .unwrap();
    assert!(result.is_some());
    match result.unwrap() {
        oxivep_cache::annotation::AnnotationValue::Positional(s) => {
            assert_eq!(s, "-2.500");
        }
        other => panic!("Expected Positional, got {:?}", other),
    }
}
