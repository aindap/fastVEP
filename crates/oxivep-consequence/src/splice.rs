use oxivep_core::Strand;
use oxivep_genome::Transcript;

/// Splice site boundaries relative to the exon-intron junction.
///
/// Donor site (5' end of intron on forward strand):
///   exon ...XY | GT...  intron
///   splice_donor: positions 1-2 into intron (GT)
///   splice_donor_5th_base: position 5 into intron
///   splice_donor_region: positions 3-6 into intron
///   splice_region (exonic side): 3 bases into exon from boundary
///   splice_region (intronic side): 3-8 bases into intron
///
/// Acceptor site (3' end of intron on forward strand):
///   intron  ...AG | XY...  exon
///   splice_acceptor: last 2 bases of intron (AG)
///   splice_polypyrimidine: 3-15 bases from end of intron
///   splice_region (exonic side): 1-3 bases into exon
///   splice_region (intronic side): 3-8 bases into intron

/// Check if a genomic position is in a splice donor site (first 2 intronic bases at 5' of intron).
pub fn is_splice_donor(transcript: &Transcript, genomic_pos: u64) -> bool {
    for_each_intron_boundary(transcript, |donor_start, donor_end, _acc_start, _acc_end| {
        genomic_pos >= donor_start && genomic_pos <= donor_end
    })
}

/// Check if a genomic position is in a splice acceptor site (last 2 intronic bases at 3' of intron).
pub fn is_splice_acceptor(transcript: &Transcript, genomic_pos: u64) -> bool {
    for_each_intron_boundary(transcript, |_donor_start, _donor_end, acc_start, acc_end| {
        genomic_pos >= acc_start && genomic_pos <= acc_end
    })
}

/// Check if position is the 5th base of the donor site.
pub fn is_splice_donor_5th_base(transcript: &Transcript, genomic_pos: u64) -> bool {
    for_each_intron_boundary_extended(
        transcript,
        |intron_start, intron_end, is_donor_at_start| {
            if is_donor_at_start {
                genomic_pos == intron_start + 4
            } else {
                genomic_pos == intron_end - 4
            }
        },
    )
}

/// Check if position is in the splice donor region (positions 3-6 of intron).
pub fn is_splice_donor_region(transcript: &Transcript, genomic_pos: u64) -> bool {
    for_each_intron_boundary_extended(
        transcript,
        |intron_start, intron_end, is_donor_at_start| {
            if is_donor_at_start {
                genomic_pos >= intron_start + 2 && genomic_pos <= intron_start + 5
            } else {
                genomic_pos >= intron_end - 5 && genomic_pos <= intron_end - 2
            }
        },
    )
}

/// Check if position is in the splice polypyrimidine tract (3-15 bases from acceptor).
pub fn is_splice_polypyrimidine_tract(transcript: &Transcript, genomic_pos: u64) -> bool {
    for_each_intron_boundary_extended(
        transcript,
        |intron_start, intron_end, is_donor_at_start| {
            // Polypyrimidine tract is near the acceptor end
            if is_donor_at_start {
                // Acceptor is at intron_end
                let acc_region_start = if intron_end >= 14 { intron_end - 14 } else { intron_start };
                genomic_pos >= acc_region_start && genomic_pos <= intron_end.saturating_sub(2)
            } else {
                // Acceptor is at intron_start
                let acc_region_end = (intron_start + 14).min(intron_end);
                genomic_pos >= intron_start + 2 && genomic_pos <= acc_region_end
            }
        },
    )
}

/// Check if position is in a splice region (3-8 bases into intron from either end,
/// or 1-3 bases into exon from the boundary).
pub fn is_splice_region(transcript: &Transcript, genomic_pos: u64) -> bool {
    let sorted_exons = sorted_exons(transcript);
    let n = sorted_exons.len();
    if n < 2 {
        return false;
    }

    for i in 0..n - 1 {
        let exon1 = sorted_exons[i];
        let exon2 = sorted_exons[i + 1];

        let intron_start = exon1.end + 1;
        let intron_end = exon2.start - 1;

        if intron_start > intron_end {
            continue;
        }

        // Intronic splice region: 3-8 bases from each intron boundary
        let donor_region_end = (intron_start + 7).min(intron_end);
        if genomic_pos >= intron_start + 2 && genomic_pos <= donor_region_end {
            return true;
        }
        let acc_region_start = if intron_end >= 7 {
            intron_end - 7
        } else {
            intron_start
        };
        if genomic_pos >= acc_region_start && genomic_pos <= intron_end.saturating_sub(2) {
            return true;
        }

        // Exonic splice region: 1-3 bases at exon boundary
        // End of exon1 (3 bases from end, towards intron)
        let exon1_region_start = if exon1.end >= 2 { exon1.end - 2 } else { exon1.start };
        if genomic_pos >= exon1_region_start && genomic_pos <= exon1.end {
            return true;
        }
        // Start of exon2 (3 bases from start, towards intron)
        let exon2_region_end = (exon2.start + 2).min(exon2.end);
        if genomic_pos >= exon2.start && genomic_pos <= exon2_region_end {
            return true;
        }
    }

    false
}

/// Helper: iterate intron boundaries and check a condition.
/// Calls `check(donor_start, donor_end, acceptor_start, acceptor_end)`.
fn for_each_intron_boundary<F>(transcript: &Transcript, check: F) -> bool
where
    F: Fn(u64, u64, u64, u64) -> bool,
{
    let sorted = sorted_exons(transcript);
    let n = sorted.len();
    if n < 2 {
        return false;
    }

    for i in 0..n - 1 {
        let exon1 = sorted[i];
        let exon2 = sorted[i + 1];

        let intron_start = exon1.end + 1;
        let intron_end = exon2.start - 1;

        if intron_start > intron_end {
            continue;
        }

        // On forward strand: donor at intron_start, acceptor at intron_end
        // On reverse strand: donor at intron_end, acceptor at intron_start
        // Since we sort exons in transcript order, the first exon boundary
        // is always the donor side in transcript terms
        let (donor_start, donor_end, acc_start, acc_end) = match transcript.strand {
            Strand::Forward => (
                intron_start,
                (intron_start + 1).min(intron_end),
                if intron_end >= 1 { intron_end - 1 } else { intron_start },
                intron_end,
            ),
            Strand::Reverse => (
                if intron_end >= 1 { intron_end - 1 } else { intron_start },
                intron_end,
                intron_start,
                (intron_start + 1).min(intron_end),
            ),
        };

        if check(donor_start, donor_end, acc_start, acc_end) {
            return true;
        }
    }

    false
}

/// Extended intron boundary helper that provides full intron coords.
fn for_each_intron_boundary_extended<F>(transcript: &Transcript, check: F) -> bool
where
    F: Fn(u64, u64, bool) -> bool,
{
    let sorted = sorted_exons(transcript);
    let n = sorted.len();
    if n < 2 {
        return false;
    }

    for i in 0..n - 1 {
        let intron_start = sorted[i].end + 1;
        let intron_end = sorted[i + 1].start - 1;

        if intron_start > intron_end {
            continue;
        }

        // is_donor_at_start: true for forward strand (donor=5' end=start of intron)
        let is_donor_at_start = transcript.strand == Strand::Forward;

        if check(intron_start, intron_end, is_donor_at_start) {
            return true;
        }
    }

    false
}

fn sorted_exons(transcript: &Transcript) -> Vec<&oxivep_genome::Exon> {
    let mut exons: Vec<&oxivep_genome::Exon> = transcript.exons.iter().collect();
    match transcript.strand {
        Strand::Forward => exons.sort_by_key(|e| e.start),
        Strand::Reverse => exons.sort_by(|a, b| b.start.cmp(&a.start)),
    }
    exons
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxivep_genome::{Exon, Gene, Transcript, Translation};

    fn make_forward_transcript() -> Transcript {
        // Exon1: 1000-1200, Intron1: 1201-1999, Exon2: 2000-2300
        Transcript {
            stable_id: "ENST_TEST".into(),
            gene: Gene {
                stable_id: "ENSG_TEST".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 2300,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 2300,
            strand: Strand::Forward,
            exons: vec![
                Exon { stable_id: "E1".into(), start: 1000, end: 1200, strand: Strand::Forward, phase: 0, end_phase: 0, rank: 1 },
                Exon { stable_id: "E2".into(), start: 2000, end: 2300, strand: Strand::Forward, phase: 0, end_phase: 0, rank: 2 },
            ],
            translation: Some(Translation { stable_id: "P1".into(), genomic_start: 1000, genomic_end: 2300, start_exon_rank: 1, start_exon_offset: 0, end_exon_rank: 2, end_exon_offset: 300 }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(502),
            coding_region_start: Some(1000),
            coding_region_end: Some(2300),
            spliced_seq: None, translateable_seq: None, peptide: None,
            canonical: false, mane_select: None, mane_plus_clinical: None,
            tsl: None, appris: None, ccds: None, protein_id: None,
            swissprot: vec![], trembl: vec![], uniparc: vec![],
            refseq_id: None, source: None, gencode_primary: false,
        }
    }

    #[test]
    fn test_splice_donor() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Donor = first 2 bases: 1201, 1202
        assert!(is_splice_donor(&tr, 1201));
        assert!(is_splice_donor(&tr, 1202));
        assert!(!is_splice_donor(&tr, 1203));
        assert!(!is_splice_donor(&tr, 1200)); // exonic
    }

    #[test]
    fn test_splice_acceptor() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Acceptor = last 2 bases: 1998, 1999
        assert!(is_splice_acceptor(&tr, 1998));
        assert!(is_splice_acceptor(&tr, 1999));
        assert!(!is_splice_acceptor(&tr, 1997));
        assert!(!is_splice_acceptor(&tr, 2000)); // exonic
    }

    #[test]
    fn test_splice_region() {
        let tr = make_forward_transcript();
        // Exonic splice region: last 3 bases of exon1 (1198, 1199, 1200)
        assert!(is_splice_region(&tr, 1198));
        assert!(is_splice_region(&tr, 1200));
        // Exonic splice region: first 3 bases of exon2 (2000, 2001, 2002)
        assert!(is_splice_region(&tr, 2000));
        assert!(is_splice_region(&tr, 2002));
        // Intronic splice region: 3-8 bases from donor (1203-1208)
        assert!(is_splice_region(&tr, 1203));
        assert!(is_splice_region(&tr, 1208));
        assert!(!is_splice_region(&tr, 1209));
        // Mid-intron: not splice region
        assert!(!is_splice_region(&tr, 1500));
    }
}
