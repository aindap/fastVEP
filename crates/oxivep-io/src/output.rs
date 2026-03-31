use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use oxivep_core::Consequence;

/// Format a VCF CSQ INFO field value from a VariationFeature.
///
/// Fields match the standard VEP CSQ format:
/// Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
/// EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
/// Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS
pub fn format_csq(vf: &VariationFeature, fields: &[&str]) -> String {
    let mut entries = Vec::new();

    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let entry = format_csq_entry(tv, aa, fields);
            entries.push(entry);
        }
    }

    entries.join(",")
}

fn format_csq_entry(tv: &TranscriptVariation, aa: &AlleleAnnotation, fields: &[&str]) -> String {
    let mut parts = Vec::new();

    for field in fields {
        let value = match *field {
            "Allele" => aa.allele.to_string(),
            "Consequence" => aa
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join("&"),
            "IMPACT" => format!("{:?}", aa.impact).to_uppercase(),
            "SYMBOL" => tv.gene_symbol.clone().unwrap_or_default(),
            "Gene" => tv.gene_id.clone(),
            "Feature_type" => "Transcript".to_string(),
            "Feature" => tv.transcript_id.clone(),
            "BIOTYPE" => tv.biotype.clone(),
            "EXON" => aa
                .exon
                .map(|(n, t)| format!("{}/{}", n, t))
                .unwrap_or_default(),
            "INTRON" => aa
                .intron
                .map(|(n, t)| format!("{}/{}", n, t))
                .unwrap_or_default(),
            "HGVSc" => aa.hgvsc.clone().unwrap_or_default(),
            "HGVSp" => aa.hgvsp.clone().unwrap_or_default(),
            "cDNA_position" => format_position_range(aa.cdna_position),
            "CDS_position" => format_position_range(aa.cds_position),
            "Protein_position" => format_position_range(aa.protein_position),
            "Amino_acids" => aa
                .amino_acids
                .as_ref()
                .map(|(r, a)| format!("{}/{}", r, a))
                .unwrap_or_default(),
            "Codons" => aa
                .codons
                .as_ref()
                .map(|(r, a)| format!("{}/{}", r, a))
                .unwrap_or_default(),
            "Existing_variation" => aa.existing_variation.join("&"),
            "DISTANCE" => aa
                .distance
                .map(|d| d.to_string())
                .unwrap_or_default(),
            "STRAND" => format!("{}", tv.strand.as_int()),
            "CANONICAL" => {
                if tv.canonical {
                    "YES".to_string()
                } else {
                    String::new()
                }
            }
            "MANE_SELECT" => tv.mane_select.clone().unwrap_or_default(),
            "MANE_PLUS_CLINICAL" => tv.mane_plus_clinical.clone().unwrap_or_default(),
            "TSL" => tv.tsl.map(|t| t.to_string()).unwrap_or_default(),
            "APPRIS" => tv.appris.clone().unwrap_or_default(),
            "CCDS" => tv.ccds.clone().unwrap_or_default(),
            "ENSP" => tv.protein_id.clone().unwrap_or_default(),
            "SIFT" => aa.sift.clone().unwrap_or_default(),
            "PolyPhen" => aa.polyphen.clone().unwrap_or_default(),
            "SOURCE" => tv.source.clone().unwrap_or_default(),
            _ => String::new(),
        };
        parts.push(escape_csq_value(&value));
    }

    parts.join("|")
}

/// Escape special characters in CSQ field values.
fn escape_csq_value(value: &str) -> String {
    value
        .replace(',', "&")
        .replace(';', "%3B")
        .replace('|', "&")
        .replace(' ', "_")
}

fn format_position_range(pos: Option<(u64, u64)>) -> String {
    match pos {
        Some((start, end)) if start == end => start.to_string(),
        Some((start, end)) => format!("{}-{}", start, end),
        None => String::new(),
    }
}

/// Default CSQ fields used by VEP.
pub const DEFAULT_CSQ_FIELDS: &[&str] = &[
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
];

/// Generate the VCF INFO header line for CSQ.
pub fn csq_header_line(fields: &[&str]) -> String {
    format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from OxiVEP. Format: {}\">",
        fields.join("|")
    )
}

/// Format a VariationFeature as a tab-delimited VEP output line.
pub fn format_tab_line(vf: &VariationFeature) -> Vec<String> {
    let mut lines = Vec::new();

    let location = if vf.position.start == vf.position.end {
        format!("{}:{}", vf.position.chromosome, vf.position.start)
    } else {
        format!(
            "{}:{}-{}",
            vf.position.chromosome, vf.position.start, vf.position.end
        )
    };

    let uploaded_variation = vf
        .variation_name
        .clone()
        .unwrap_or_else(|| format!("{}_{}", location, vf.allele_string));

    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let consequence_str = aa
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join(",");

            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                uploaded_variation,
                location,
                aa.allele,
                tv.gene_id,
                tv.transcript_id,
                "Transcript",
                consequence_str,
                format_position_range(aa.cdna_position),
                format_position_range(aa.cds_position),
                format_position_range(aa.protein_position),
                aa.amino_acids
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or_default(),
                aa.codons
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or_default(),
                aa.existing_variation.join(","),
            );
            lines.push(line);
        }
    }

    // If no transcript annotations, still output the variant with intergenic
    if vf.transcript_variations.is_empty() {
        for alt in &vf.alt_alleles {
            let line = format!(
                "{}\t{}\t{}\t-\t-\t-\t{}\t-\t-\t-\t-\t-\t-",
                uploaded_variation,
                location,
                alt,
                Consequence::IntergenicVariant.so_term(),
            );
            lines.push(line);
        }
    }

    lines
}

/// Format a VariationFeature as JSON.
pub fn format_json(vf: &VariationFeature) -> serde_json::Value {
    let mut obj = serde_json::Map::new();

    obj.insert("id".into(), json_str(&vf.variation_name));
    obj.insert(
        "seq_region_name".into(),
        serde_json::Value::String(vf.position.chromosome.clone()),
    );
    obj.insert("start".into(), serde_json::Value::Number(vf.position.start.into()));
    obj.insert("end".into(), serde_json::Value::Number(vf.position.end.into()));
    obj.insert(
        "allele_string".into(),
        serde_json::Value::String(vf.allele_string.clone()),
    );
    obj.insert("strand".into(), serde_json::Value::Number(1.into()));

    if let Some(ref msq) = vf.most_severe_consequence {
        obj.insert(
            "most_severe_consequence".into(),
            serde_json::Value::String(msq.so_term().to_string()),
        );
    }

    let transcript_consequences: Vec<serde_json::Value> = vf
        .transcript_variations
        .iter()
        .flat_map(|tv| {
            tv.allele_annotations.iter().map(move |aa| {
                let mut tc = serde_json::Map::new();
                tc.insert(
                    "gene_id".into(),
                    serde_json::Value::String(tv.gene_id.clone()),
                );
                tc.insert(
                    "transcript_id".into(),
                    serde_json::Value::String(tv.transcript_id.clone()),
                );
                tc.insert(
                    "biotype".into(),
                    serde_json::Value::String(tv.biotype.clone()),
                );
                if let Some(ref sym) = tv.gene_symbol {
                    tc.insert(
                        "gene_symbol".into(),
                        serde_json::Value::String(sym.clone()),
                    );
                }
                tc.insert(
                    "consequence_terms".into(),
                    serde_json::Value::Array(
                        aa.consequences
                            .iter()
                            .map(|c| serde_json::Value::String(c.so_term().to_string()))
                            .collect(),
                    ),
                );
                tc.insert(
                    "impact".into(),
                    serde_json::Value::String(format!("{:?}", aa.impact).to_uppercase()),
                );
                tc.insert(
                    "variant_allele".into(),
                    serde_json::Value::String(aa.allele.to_string()),
                );
                serde_json::Value::Object(tc)
            })
        })
        .collect();

    obj.insert(
        "transcript_consequences".into(),
        serde_json::Value::Array(transcript_consequences),
    );

    serde_json::Value::Object(obj)
}

fn json_str(opt: &Option<String>) -> serde_json::Value {
    match opt {
        Some(s) => serde_json::Value::String(s.clone()),
        None => serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape_csq_value() {
        assert_eq!(escape_csq_value("hello,world"), "hello&world");
        assert_eq!(escape_csq_value("a;b"), "a%3Bb");
        assert_eq!(escape_csq_value("a|b"), "a&b");
        assert_eq!(escape_csq_value("a b"), "a_b");
    }

    #[test]
    fn test_csq_header() {
        let header = csq_header_line(&["Allele", "Consequence"]);
        assert!(header.contains("Format: Allele|Consequence"));
    }

    #[test]
    fn test_format_position_range() {
        assert_eq!(format_position_range(Some((100, 100))), "100");
        assert_eq!(format_position_range(Some((100, 200))), "100-200");
        assert_eq!(format_position_range(None), "");
    }
}
