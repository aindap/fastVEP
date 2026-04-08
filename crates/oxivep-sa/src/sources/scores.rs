//! Generic position-level score parsers (PhyloP, GERP, DANN, etc.).
//!
//! These parsers produce positional AnnotationRecords where the JSON is
//! just the numeric score value as a string. The SaWriter stores them
//! as positional annotations (match_by_allele=false, is_positional=true).

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a tab-separated score file (chrom, pos, score) into AnnotationRecords.
///
/// Supports formats like:
/// - BED-like: `chr1\t12345\t12346\t2.345`  (4 columns: chrom, start, end, score)
/// - Simple:   `chr1\t12345\t2.345`          (3 columns: chrom, pos, score)
///
/// Positions are 0-based in BED format (converted to 1-based internally)
/// or 1-based in simple format. Set `zero_based` accordingly.
pub fn parse_score_tsv<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    zero_based: bool,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading score TSV line")?;
        if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        let (chrom_str, pos_str, score_str) = match fields.len() {
            // BED 4-column: chrom, start (0-based), end, score
            4.. => (fields[0], fields[1], fields[3]),
            // Simple 3-column: chrom, pos, score
            3 => (fields[0], fields[1], fields[2]),
            _ => continue,
        };

        let chrom = normalize_chrom(chrom_str);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match pos_str.parse::<u32>() {
            Ok(p) => {
                if zero_based { p + 1 } else { p }
            }
            Err(_) => continue,
        };

        // Validate score is a number
        let score: f64 = match score_str.trim().parse() {
            Ok(s) => s,
            Err(_) => continue,
        };

        // Format score compactly
        let json = format_score(score);

        records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele: String::new(),
            alt_allele: String::new(),
            json,
        });
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Parse a UCSC wiggle fixed-step (wigFix) file into AnnotationRecords.
///
/// Format:
/// ```text
/// fixedStep chrom=chr1 start=10001 step=1
/// 0.064
/// -0.002
/// ...
/// ```
pub fn parse_wigfix<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mut current_chrom_idx: Option<u16> = None;
    let mut current_pos: u32 = 0;
    let mut step: u32 = 1;

    for line in reader.lines() {
        let line = line.context("Reading wigFix line")?;

        if line.starts_with("fixedStep") {
            // Parse header: fixedStep chrom=chr1 start=10001 step=1
            let mut chrom = None;
            let mut start = None;
            step = 1;

            for part in line.split_whitespace().skip(1) {
                if let Some((key, val)) = part.split_once('=') {
                    match key {
                        "chrom" => chrom = Some(normalize_chrom(val)),
                        "start" => start = val.parse().ok(),
                        "step" => step = val.parse().unwrap_or(1),
                        _ => {}
                    }
                }
            }

            current_chrom_idx = chrom.as_ref().and_then(|c| chrom_to_idx.get(c)).copied();
            current_pos = start.unwrap_or(1);
            continue;
        }

        if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
            continue;
        }

        if let Some(chrom_idx) = current_chrom_idx {
            if let Ok(score) = line.trim().parse::<f64>() {
                records.push(AnnotationRecord {
                    chrom_idx,
                    position: current_pos,
                    ref_allele: String::new(),
                    alt_allele: String::new(),
                    json: format_score(score),
                });
            }
            current_pos += step;
        }
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Format a score compactly: drop trailing zeros, max 4 decimal places.
fn format_score(score: f64) -> String {
    if score == 0.0 {
        return "0".into();
    }
    // Use up to 4 decimal places, strip trailing zeros
    let s = format!("{:.4}", score);
    let s = s.trim_end_matches('0');
    let s = s.trim_end_matches('.');
    s.to_string()
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_chrom_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0);
        m.insert("chr2".into(), 1);
        m
    }

    #[test]
    fn test_parse_score_tsv_bed4() {
        let data = "\
chr1\t99\t100\t2.345
chr1\t199\t200\t-1.5
chr2\t49\t50\t0.001
";
        let records = parse_score_tsv(data.as_bytes(), &test_chrom_map(), true).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 100); // 0-based 99 -> 1-based 100
        assert_eq!(records[0].json, "2.345");
        assert_eq!(records[1].json, "-1.5");
        assert_eq!(records[2].chrom_idx, 1);
        assert_eq!(records[2].json, "0.001");
    }

    #[test]
    fn test_parse_score_tsv_simple() {
        let data = "chr1\t100\t3.14\nchr1\t200\t-0.5\n";
        let records = parse_score_tsv(data.as_bytes(), &test_chrom_map(), false).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].position, 100);
        assert_eq!(records[0].json, "3.14");
    }

    #[test]
    fn test_parse_wigfix() {
        let data = "\
fixedStep chrom=chr1 start=10001 step=1
0.064
-0.002
1.500
fixedStep chrom=chr2 start=5000 step=1
-3.210
";
        let records = parse_wigfix(data.as_bytes(), &test_chrom_map()).unwrap();
        assert_eq!(records.len(), 4);
        // chr1 records
        assert_eq!(records[0].chrom_idx, 0);
        assert_eq!(records[0].position, 10001);
        assert_eq!(records[0].json, "0.064");
        assert_eq!(records[1].position, 10002);
        assert_eq!(records[1].json, "-0.002");
        assert_eq!(records[2].position, 10003);
        assert_eq!(records[2].json, "1.5");
        // chr2 record
        assert_eq!(records[3].chrom_idx, 1);
        assert_eq!(records[3].position, 5000);
        assert_eq!(records[3].json, "-3.21");
    }

    #[test]
    fn test_format_score() {
        assert_eq!(format_score(2.3450), "2.345");
        assert_eq!(format_score(0.0), "0");
        assert_eq!(format_score(-1.0), "-1");
        assert_eq!(format_score(0.12345), "0.1235"); // rounds to 4dp
        assert_eq!(format_score(3.0001), "3.0001");
    }
}
