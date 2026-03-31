use oxivep_core::Allele;

/// Generate HGVSc (coding DNA) notation.
///
/// Uses the transcript ID as the reference and cDNA position numbering.
/// Format: ENST00000001.1:c.123A>G
pub fn hgvsc(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
) -> Option<String> {
    let prefix = format!("{}:c.", transcript_id);

    // Convert cDNA position to CDS position (relative to ATG)
    let cds_pos_start = cdna_start as i64 - coding_start as i64 + 1;
    let cds_pos_end = cdna_end as i64 - coding_start as i64 + 1;

    let pos_str = if cds_pos_start < 0 {
        // 5' UTR: negative positions
        format!("{}", cds_pos_start)
    } else if cds_pos_start == cds_pos_end {
        format!("{}", cds_pos_start)
    } else {
        format!("{}_{}", cds_pos_start, cds_pos_end)
    };

    let notation = match (ref_allele, alt_allele) {
        // SNV
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix,
                pos_str,
                ref_bases[0] as char,
                alt_bases[0] as char
            )
        }
        // Deletion
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            if ref_bases.len() == 1 {
                format!("{}{}del", prefix, pos_str)
            } else {
                format!("{}{}del", prefix, pos_str)
            }
        }
        // Insertion
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_pos = format!("{}_{}", cds_pos_end, cds_pos_end + 1);
            format!(
                "{}{}ins{}",
                prefix,
                ins_pos,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        // MNV or complex
        (Allele::Sequence(_), Allele::Sequence(alt_bases)) => {
            format!(
                "{}{}delins{}",
                prefix,
                pos_str,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        _ => return None,
    };

    Some(notation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvsc_snv() {
        let result = hgvsc(
            "ENST00000001",
            151, 151,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51, // coding_start in cDNA
        );
        // cds_pos = 151 - 51 + 1 = 101
        assert_eq!(result, Some("ENST00000001:c.101A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion() {
        let result = hgvsc(
            "ENST00000001",
            54, 56,
            &Allele::Sequence(b"ACG".to_vec()),
            &Allele::Deletion,
            51,
        );
        // cds_pos_start = 54-51+1=4, end=56-51+1=6
        assert_eq!(result, Some("ENST00000001:c.4_6del".to_string()));
    }

    #[test]
    fn test_hgvsc_insertion() {
        let result = hgvsc(
            "ENST00000001",
            54, 53, // insertion between 53 and 54 in cDNA
            &Allele::Deletion,
            &Allele::Sequence(b"TTT".to_vec()),
            51,
        );
        // cds_pos_end = 53-51+1=3
        assert_eq!(result, Some("ENST00000001:c.3_4insTTT".to_string()));
    }

    #[test]
    fn test_hgvsc_5_utr() {
        let result = hgvsc(
            "ENST00000001",
            10, 10,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51,
        );
        // cds_pos = 10 - 51 + 1 = -40
        assert_eq!(result, Some("ENST00000001:c.-40A>G".to_string()));
    }
}
