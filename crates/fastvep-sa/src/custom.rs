//! Custom annotation providers.
//!
//! Allows users to provide their own annotation files (VCF, BED, TSV)
//! that get indexed and queried at annotation time.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a custom VCF annotation file.
///
/// Extracts specified INFO fields as JSON annotations.
/// If `info_fields` is empty, all INFO fields are included.
pub fn parse_custom_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    name: &str,
    info_fields: &[String],
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading custom VCF")?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = normalize_chrom(fields[0]);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];
        let info = fields[7];

        let info_map = parse_info(info);

        // Build JSON from requested INFO fields
        let mut parts = Vec::new();
        if info_fields.is_empty() {
            // Include all INFO fields
            for (key, val) in &info_map {
                parts.push(format!("\"{}\":\"{}\"", key, escape_json(val)));
            }
        } else {
            for field in info_fields {
                if let Some(val) = info_map.get(field.as_str()) {
                    parts.push(format!("\"{}\":\"{}\"", field, escape_json(val)));
                }
            }
        }

        if parts.is_empty() {
            parts.push(format!("\"source\":\"{}\"", name));
        }

        let json = format!("{{{}}}", parts.join(","));

        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }
            records.push(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.to_string(),
                json: json.clone(),
            });
        }
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Parse a custom BED annotation file into interval records.
///
/// Format: chrom, start (0-based), end, name, [score], [additional fields]
pub fn parse_custom_bed<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<crate::common::IntervalRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading custom BED")?;
        if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chrom = normalize_chrom(fields[0]);
        if !chrom_to_idx.contains_key(&chrom) {
            continue;
        }

        let start: u32 = match fields[1].parse::<u32>() {
            Ok(s) => s + 1, // BED is 0-based, convert to 1-based
            Err(_) => continue,
        };
        let end: u32 = match fields[2].parse() {
            Ok(e) => e,
            Err(_) => continue,
        };

        let name = fields.get(3).unwrap_or(&".").to_string();
        let score = fields.get(4).unwrap_or(&".").to_string();

        let mut parts = Vec::new();
        if name != "." {
            parts.push(format!("\"name\":\"{}\"", escape_json(&name)));
        }
        if score != "." {
            if let Ok(s) = score.parse::<f64>() {
                parts.push(format!("\"score\":{}", s));
            }
        }

        records.push(crate::common::IntervalRecord {
            chrom,
            start,
            end,
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    Ok(records)
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in info.split(';') {
        if let Some((k, v)) = pair.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    map
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") { chrom.to_string() } else { format!("chr{}", chrom) }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_custom_vcf() {
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tMY_SCORE=0.95;MY_FLAG=true\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);

        // All fields
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "test", &[]).unwrap();
        assert_eq!(recs.len(), 1);
        assert!(recs[0].json.contains("MY_SCORE"));

        // Specific fields
        let recs = parse_custom_vcf(
            vcf.as_bytes(), &m, "test",
            &["MY_SCORE".to_string()],
        ).unwrap();
        assert!(recs[0].json.contains("MY_SCORE"));
        assert!(!recs[0].json.contains("MY_FLAG"));
    }

    #[test]
    fn test_parse_custom_bed() {
        let bed = "chr1\t99\t200\tregion1\t0.5\nchr1\t499\t600\tregion2\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_bed(bed.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].start, 100); // 0-based -> 1-based
        assert_eq!(recs[0].end, 200);
        assert!(recs[0].json.contains("region1"));
        assert!(recs[0].json.contains("0.5"));
    }
}
