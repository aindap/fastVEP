use oxivep_core::Strand;
use serde::{Deserialize, Serialize};

/// A gene model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub stable_id: String,
    pub symbol: Option<String>,
    pub symbol_source: Option<String>,
    pub hgnc_id: Option<String>,
    pub biotype: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

/// A transcript model with all data needed for consequence prediction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcript {
    pub stable_id: String,
    pub gene: Gene,
    pub biotype: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub translation: Option<Translation>,

    // Pre-computed fields for consequence prediction
    /// Start of the coding region in cDNA coordinates (1-based).
    pub cdna_coding_start: Option<u64>,
    /// End of the coding region in cDNA coordinates (1-based).
    pub cdna_coding_end: Option<u64>,
    /// Start of the coding region in genomic coordinates.
    pub coding_region_start: Option<u64>,
    /// End of the coding region in genomic coordinates.
    pub coding_region_end: Option<u64>,

    // Spliced and translated sequences
    pub spliced_seq: Option<String>,
    pub translateable_seq: Option<String>,
    pub peptide: Option<String>,

    // Annotation metadata
    pub canonical: bool,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    pub ccds: Option<String>,
    pub protein_id: Option<String>,
    pub swissprot: Vec<String>,
    pub trembl: Vec<String>,
    pub uniparc: Vec<String>,
    pub refseq_id: Option<String>,
    pub source: Option<String>,
    pub gencode_primary: bool,
}

impl Transcript {
    /// Whether this transcript is protein-coding.
    pub fn is_coding(&self) -> bool {
        self.translation.is_some()
    }

    /// Total number of exons.
    pub fn exon_count(&self) -> usize {
        self.exons.len()
    }

    /// Total number of introns (exons - 1, minimum 0).
    pub fn intron_count(&self) -> usize {
        self.exons.len().saturating_sub(1)
    }

    /// Get the cDNA length (sum of exon lengths).
    pub fn cdna_length(&self) -> u64 {
        self.exons.iter().map(|e| e.end - e.start + 1).sum()
    }

    /// Map a genomic position to a cDNA position.
    /// Returns None if the position is not in an exon.
    pub fn genomic_to_cdna(&self, genomic_pos: u64) -> Option<u64> {
        let mut cdna_pos = 0u64;
        let sorted_exons = self.sorted_exons();

        for exon in &sorted_exons {
            if genomic_pos >= exon.start && genomic_pos <= exon.end {
                return match self.strand {
                    Strand::Forward => Some(cdna_pos + (genomic_pos - exon.start) + 1),
                    Strand::Reverse => Some(cdna_pos + (exon.end - genomic_pos) + 1),
                };
            }
            cdna_pos += exon.end - exon.start + 1;
        }
        None
    }

    /// Map a cDNA position to a CDS position.
    /// Returns None if the position is not in the coding region.
    pub fn cdna_to_cds(&self, cdna_pos: u64) -> Option<u64> {
        let coding_start = self.cdna_coding_start?;
        let coding_end = self.cdna_coding_end?;
        if cdna_pos >= coding_start && cdna_pos <= coding_end {
            Some(cdna_pos - coding_start + 1)
        } else {
            None
        }
    }

    /// Map a CDS position (1-based) to a protein position (1-based).
    pub fn cds_to_protein(cds_pos: u64) -> u64 {
        (cds_pos - 1) / 3 + 1
    }

    /// Get exons sorted by transcript order (forward for +, reverse for -).
    fn sorted_exons(&self) -> Vec<&Exon> {
        let mut exons: Vec<&Exon> = self.exons.iter().collect();
        match self.strand {
            Strand::Forward => exons.sort_by_key(|e| e.start),
            Strand::Reverse => exons.sort_by(|a, b| b.start.cmp(&a.start)),
        }
        exons
    }

    /// Determine which exon (0-indexed rank) a genomic position falls in.
    /// Returns (exon_index, total_exons) or None if intronic/outside.
    pub fn exon_at(&self, genomic_pos: u64) -> Option<(usize, usize)> {
        let sorted = self.sorted_exons();
        for (i, exon) in sorted.iter().enumerate() {
            if genomic_pos >= exon.start && genomic_pos <= exon.end {
                return Some((i, sorted.len()));
            }
        }
        None
    }

    /// Determine which intron (0-indexed rank) a genomic position falls in.
    /// Returns (intron_index, total_introns) or None if exonic/outside.
    pub fn intron_at(&self, genomic_pos: u64) -> Option<(usize, usize)> {
        let sorted = self.sorted_exons();
        let n_introns = sorted.len().saturating_sub(1);
        for i in 0..n_introns {
            let intron_start = sorted[i].end + 1;
            let intron_end = sorted[i + 1].start - 1;
            if genomic_pos >= intron_start && genomic_pos <= intron_end {
                return Some((i, n_introns));
            }
        }
        None
    }
}

/// An exon within a transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Exon {
    pub stable_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    /// Reading frame phase at start of exon (-1, 0, 1, 2).
    pub phase: i8,
    /// Reading frame phase at end of exon.
    pub end_phase: i8,
    /// 1-based rank of the exon in the transcript.
    pub rank: u32,
}

impl Exon {
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }
}

/// Translation metadata for a protein-coding transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Translation {
    pub stable_id: String,
    /// Genomic start of the translation (start codon first base).
    pub genomic_start: u64,
    /// Genomic end of the translation (stop codon last base).
    pub genomic_end: u64,
    /// The exon rank where translation starts.
    pub start_exon_rank: u32,
    /// Offset within the start exon (0-based).
    pub start_exon_offset: u64,
    /// The exon rank where translation ends.
    pub end_exon_rank: u32,
    /// Offset within the end exon (0-based).
    pub end_exon_offset: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST00000001".into(),
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TESTGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: Some("HGNC:1".into()),
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE00000001".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "ENSE00000002".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 1,
                    rank: 2,
                },
                Exon {
                    stable_id: "ENSE00000003".into(),
                    start: 4000,
                    end: 5000,
                    strand: Strand::Forward,
                    phase: 1,
                    end_phase: -1,
                    rank: 3,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP00000001".into(),
                genomic_start: 1050,
                genomic_end: 4500,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 3,
                end_exon_offset: 500,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(952),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP00000001".into()),
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
        }
    }

    #[test]
    fn test_is_coding() {
        let tr = make_test_transcript();
        assert!(tr.is_coding());
    }

    #[test]
    fn test_exon_count() {
        let tr = make_test_transcript();
        assert_eq!(tr.exon_count(), 3);
        assert_eq!(tr.intron_count(), 2);
    }

    #[test]
    fn test_genomic_to_cdna_forward() {
        let tr = make_test_transcript();
        // Position in first exon
        assert_eq!(tr.genomic_to_cdna(1000), Some(1));
        assert_eq!(tr.genomic_to_cdna(1200), Some(201));
        // Position in second exon
        assert_eq!(tr.genomic_to_cdna(2000), Some(202));
        // Position in intron
        assert_eq!(tr.genomic_to_cdna(1500), None);
    }

    #[test]
    fn test_exon_at() {
        let tr = make_test_transcript();
        assert_eq!(tr.exon_at(1100), Some((0, 3)));
        assert_eq!(tr.exon_at(2100), Some((1, 3)));
        assert_eq!(tr.exon_at(4500), Some((2, 3)));
        assert_eq!(tr.exon_at(1500), None); // intron
    }

    #[test]
    fn test_intron_at() {
        let tr = make_test_transcript();
        assert_eq!(tr.intron_at(1500), Some((0, 2)));
        assert_eq!(tr.intron_at(3000), Some((1, 2)));
        assert_eq!(tr.intron_at(1100), None); // exon
    }

    #[test]
    fn test_cds_to_protein() {
        assert_eq!(Transcript::cds_to_protein(1), 1);
        assert_eq!(Transcript::cds_to_protein(3), 1);
        assert_eq!(Transcript::cds_to_protein(4), 2);
        assert_eq!(Transcript::cds_to_protein(6), 2);
        assert_eq!(Transcript::cds_to_protein(7), 3);
    }
}
