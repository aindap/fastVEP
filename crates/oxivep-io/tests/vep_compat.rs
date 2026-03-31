/// Integration tests validating OxiVEP against Ensembl VEP test patterns.
///
/// These tests verify that OxiVEP's VCF parsing, allele normalization,
/// and variant representation matches the behavior documented in
/// ensembl-vep's Parser_VCF.t test suite.
use oxivep_core::{Allele, Strand};
use oxivep_io::vcf::{parse_vcf_line, VcfParser};

// =============================================================================
// VCF Parsing — matches ensembl-vep Parser_VCF.t assertions
// =============================================================================

#[test]
fn test_vep_snv_basic_parsing() {
    // From ensembl-vep test.vcf line 1: 21 25585733 rs142513484 C T
    // VEP Parser_VCF.t expects: chr=21, start=25585733, end=25585733, allele_string=C/T
    let line = "21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "21");
    assert_eq!(vf.position.start, 25585733);
    assert_eq!(vf.position.end, 25585733);
    assert_eq!(vf.allele_string, "C/T");
    assert_eq!(vf.variation_name, Some("rs142513484".to_string()));
    assert_eq!(vf.position.strand, Strand::Forward);
    assert!(!vf.is_indel());
}

#[test]
fn test_vep_insertion_parsing() {
    // From test_not_ordered.vcf: 21 30000016 rs202173120 T TCA
    // VEP behavior: strip shared first base T → ref="-", alt="CA", start=30000017, end=30000016
    let line = "21\t30000016\trs202173120\tT\tTCA\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 30000017);
    assert_eq!(vf.position.end, 30000016); // insertion: end < start
    assert_eq!(vf.allele_string, "-/CA");
    assert_eq!(vf.ref_allele, Allele::Deletion);
    assert_eq!(vf.alt_alleles[0], Allele::Sequence(b"CA".to_vec()));
    assert!(vf.is_insertion());
}

#[test]
fn test_vep_insertion_single_base() {
    // 21 30000105 rs34788396 C CC → insertion of 1 base
    let line = "21\t30000105\trs34788396\tC\tCC\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 30000106);
    assert_eq!(vf.position.end, 30000105);
    assert_eq!(vf.allele_string, "-/C");
    assert!(vf.is_insertion());
}

#[test]
fn test_vep_large_insertion() {
    // 21 30000094 rs960445955 C CATATTCTCCCCTATT → large insertion
    let line = "21\t30000094\trs960445955\tC\tCATATTCTCCCCTATT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 30000095);
    assert_eq!(vf.position.end, 30000094);
    assert_eq!(vf.allele_string, "-/ATATTCTCCCCTATT");
    assert!(vf.is_insertion());
}

#[test]
fn test_vep_deletion_parsing() {
    // Simulated deletion: ref=GGA, alt=G at chr22:19353532
    // VEP strips shared G: ref=GA, alt=-, start=19353533
    let line = "22\t19353532\trs1162718428\tGGA\tG\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "22");
    assert_eq!(vf.position.start, 19353533);
    assert_eq!(vf.position.end, 19353534);
    assert_eq!(vf.allele_string, "GA/-");
    assert_eq!(vf.ref_allele, Allele::Sequence(b"GA".to_vec()));
    assert_eq!(vf.alt_alleles[0], Allele::Deletion);
    assert!(vf.is_deletion());
}

#[test]
fn test_vep_multi_allelic_snv() {
    // Multi-allelic SNV: two alt alleles
    let line = "21\t25585733\t.\tC\tT,G\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.allele_string, "C/T/G");
    assert_eq!(vf.alt_alleles.len(), 2);
    assert_eq!(vf.alt_alleles[0], Allele::Sequence(b"T".to_vec()));
    assert_eq!(vf.alt_alleles[1], Allele::Sequence(b"G".to_vec()));
}

#[test]
fn test_vep_multi_allelic_indel() {
    // Multi-allelic with indel: ref=ACG, alt=A,ACGT
    // All share first base A, so strip: ref=CG, alt=-,CGT
    let line = "21\t25585733\t.\tACG\tA,ACGT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 25585734);
    assert_eq!(vf.allele_string, "CG/-/CGT");
}

#[test]
fn test_vep_star_allele_handling() {
    // Star allele should be preserved during normalization
    let line = "21\t25585733\t.\tACG\tA,*\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 25585734);
    assert!(vf.allele_string.contains("*"));
}

#[test]
fn test_vep_non_variant_site() {
    // REF-only site: alt is "."
    let line = "21\t25585733\t.\tC\t.\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.allele_string, "C");
}

#[test]
fn test_vep_chr_prefix_handling() {
    // Chromosome with "chr" prefix should be preserved
    let line = "chr21\t25585733\trs1\tC\tT\t.\tPASS\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "chr21");
}

#[test]
fn test_vep_mt_chromosome() {
    // Mitochondrial chromosome
    let line = "MT\t4472\t.\tT\tA\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "MT");
    assert_eq!(vf.position.start, 4472);
    assert_eq!(vf.allele_string, "T/A");
}

#[test]
fn test_vep_mnv_parsing() {
    // MNV: multi-nucleotide variant (same ref/alt length >1)
    let line = "21\t25585733\t.\tAC\tGT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 25585733);
    assert_eq!(vf.position.end, 25585734);
    assert_eq!(vf.allele_string, "AC/GT");
    assert!(!vf.is_indel());
}

#[test]
fn test_vep_complex_indel() {
    // Complex indel: different ref/alt lengths, not simple ins/del
    let line = "21\t25585733\t.\tACG\tTT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    // No shared first base, so no stripping
    assert_eq!(vf.position.start, 25585733);
    assert_eq!(vf.allele_string, "ACG/TT");
}

// =============================================================================
// VCF Parser — multi-line / header handling
// =============================================================================

#[test]
fn test_vep_full_vcf_parser() {
    // Simulates the structure of ensembl-vep test.vcf
    let vcf = "##fileformat=VCFv4.1\n\
##contig=<ID=21,assembly=GCF_000001405.26,length=46709983>\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"SV length\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\n\
21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0\n\
21\t25587701\trs187353664\tT\tC\t.\t.\t.\tGT\t0|0\n\
21\t25587758\trs116645811\tG\tA\t.\t.\t.\tGT\t0|0\n";

    let mut parser = VcfParser::new(vcf.as_bytes()).unwrap();

    // Should have 4 header lines
    assert_eq!(parser.header_lines().len(), 4);
    assert!(parser.header_lines()[0].starts_with("##fileformat"));
    assert!(parser.header_lines()[3].starts_with("#CHROM"));

    let variants = parser.read_all().unwrap();
    assert_eq!(variants.len(), 3);

    assert_eq!(variants[0].variation_name.as_deref(), Some("rs142513484"));
    assert_eq!(variants[0].position.start, 25585733);
    assert_eq!(variants[0].allele_string, "C/T");

    assert_eq!(variants[1].variation_name.as_deref(), Some("rs187353664"));
    assert_eq!(variants[2].variation_name.as_deref(), Some("rs116645811"));
}

#[test]
fn test_vep_unordered_variants() {
    // Variants from test_not_ordered.vcf — mixed chromosomes and positions
    let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t230710034\trs894744940\tC\tT\t.\t.\t.\n\
21\t25587758\trs116645811\tG\tA\t.\t.\t.\n\
21\t25592836\trs1135638\tG\tA\t.\t.\t.\n\
1\t230710045\trs754176245\tC\tT\t.\t.\t.\n";

    let mut parser = VcfParser::new(vcf.as_bytes()).unwrap();
    let variants = parser.read_all().unwrap();

    assert_eq!(variants.len(), 4);
    // Verify all parsed regardless of order
    assert_eq!(variants[0].position.chromosome, "1");
    assert_eq!(variants[1].position.chromosome, "21");
    assert_eq!(variants[2].position.chromosome, "21");
    assert_eq!(variants[3].position.chromosome, "1");
}

#[test]
fn test_vep_indel_variety() {
    // Various indel types from test_not_ordered.vcf
    let test_cases = vec![
        // (line, expected_start, expected_allele_string)
        ("21\t30000016\trs202173120\tT\tTCA\t.\t.\t.", 30000017, "-/CA"),
        ("21\t30000018\trs144102269\tA\tAAT\t.\t.\t.", 30000019, "-/AT"),
        ("21\t30000020\trs35748481\tG\tGAT\t.\t.\t.", 30000021, "-/AT"),
        ("21\t30000105\trs34788396\tC\tCC\t.\t.\t.", 30000106, "-/C"),
    ];

    for (line, expected_start, expected_alleles) in test_cases {
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(
            vf.position.start, expected_start,
            "Failed for: {} — got start={}", line, vf.position.start
        );
        assert_eq!(
            vf.allele_string, expected_alleles,
            "Failed for: {} — got alleles={}", line, vf.allele_string
        );
    }
}

#[test]
fn test_vep_filter_status_preserved() {
    let line = "21\t25585733\trs1\tC\tT\t30\tPASS\tDP=50";
    let vf = parse_vcf_line(line).unwrap();

    let fields = vf.vcf_fields.unwrap();
    assert_eq!(fields.filter, "PASS");
    assert_eq!(fields.qual, "30");
    assert_eq!(fields.info, "DP=50");
}

#[test]
fn test_vep_genotype_fields_preserved() {
    let line = "21\t25585733\trs1\tC\tT\t.\t.\t.\tGT\t0|0";
    let vf = parse_vcf_line(line).unwrap();

    let fields = vf.vcf_fields.unwrap();
    assert_eq!(fields.rest.len(), 2);
    assert_eq!(fields.rest[0], "GT");
    assert_eq!(fields.rest[1], "0|0");
}

// =============================================================================
// Consequence prediction validation
// =============================================================================

#[test]
fn test_all_consequence_types_have_valid_so_terms() {
    use oxivep_core::Consequence;

    // Verify all SO terms are valid and round-trip
    let all = [
        Consequence::TranscriptAblation,
        Consequence::SpliceAcceptorVariant,
        Consequence::SpliceDonorVariant,
        Consequence::StopGained,
        Consequence::FrameshiftVariant,
        Consequence::StopLost,
        Consequence::StartLost,
        Consequence::TranscriptAmplification,
        Consequence::FeatureElongation,
        Consequence::FeatureTruncation,
        Consequence::InframeInsertion,
        Consequence::InframeDeletion,
        Consequence::MissenseVariant,
        Consequence::ProteinAlteringVariant,
        Consequence::SpliceRegionVariant,
        Consequence::SpliceDonorFifthBaseVariant,
        Consequence::SpliceDonorRegionVariant,
        Consequence::SplicePolypyrimidineTractVariant,
        Consequence::IncompleteTerminalCodonVariant,
        Consequence::StartRetainedVariant,
        Consequence::StopRetainedVariant,
        Consequence::SynonymousVariant,
        Consequence::CodingSequenceVariant,
        Consequence::MatureMirnaVariant,
        Consequence::FivePrimeUtrVariant,
        Consequence::ThreePrimeUtrVariant,
        Consequence::NonCodingTranscriptExonVariant,
        Consequence::IntronVariant,
        Consequence::NmdTranscriptVariant,
        Consequence::NonCodingTranscriptVariant,
        Consequence::CodingTranscriptVariant,
        Consequence::UpstreamGeneVariant,
        Consequence::DownstreamGeneVariant,
        Consequence::TfbsAblation,
        Consequence::TfbsAmplification,
        Consequence::TfBindingSiteVariant,
        Consequence::RegulatoryRegionAblation,
        Consequence::RegulatoryRegionAmplification,
        Consequence::RegulatoryRegionVariant,
        Consequence::IntergenicVariant,
        Consequence::SequenceVariant,
    ];

    for c in &all {
        let term = c.so_term();
        let parsed = Consequence::from_so_term(term);
        assert_eq!(parsed, Some(*c), "SO term round-trip failed for {:?} ({})", c, term);
    }

    // Verify strict ordering by rank
    for i in 0..all.len() - 1 {
        assert!(
            all[i].rank() < all[i + 1].rank(),
            "{:?} (rank {}) should be more severe than {:?} (rank {})",
            all[i], all[i].rank(), all[i + 1], all[i + 1].rank()
        );
    }
}

#[test]
fn test_vep_impact_classification() {
    use oxivep_core::{Consequence, Impact};

    // Verify impact matches VEP's classification
    let high = [
        Consequence::TranscriptAblation,
        Consequence::SpliceAcceptorVariant,
        Consequence::SpliceDonorVariant,
        Consequence::StopGained,
        Consequence::FrameshiftVariant,
        Consequence::StopLost,
        Consequence::StartLost,
    ];
    for c in &high {
        assert_eq!(c.impact(), Impact::High, "{:?} should be HIGH", c);
    }

    let moderate = [
        Consequence::InframeInsertion,
        Consequence::InframeDeletion,
        Consequence::MissenseVariant,
        Consequence::ProteinAlteringVariant,
    ];
    for c in &moderate {
        assert_eq!(c.impact(), Impact::Moderate, "{:?} should be MODERATE", c);
    }

    let low = [
        Consequence::SpliceRegionVariant,
        Consequence::SynonymousVariant,
        Consequence::StopRetainedVariant,
        Consequence::StartRetainedVariant,
    ];
    for c in &low {
        assert_eq!(c.impact(), Impact::Low, "{:?} should be LOW", c);
    }

    let modifier = [
        Consequence::IntronVariant,
        Consequence::UpstreamGeneVariant,
        Consequence::DownstreamGeneVariant,
        Consequence::IntergenicVariant,
        Consequence::FivePrimeUtrVariant,
        Consequence::ThreePrimeUtrVariant,
    ];
    for c in &modifier {
        assert_eq!(c.impact(), Impact::Modifier, "{:?} should be MODIFIER", c);
    }
}

// =============================================================================
// Codon table validation
// =============================================================================

#[test]
fn test_all_64_codons_translate() {
    use oxivep_genome::CodonTable;

    let table = CodonTable::standard();
    let bases = [b'A', b'C', b'G', b'T'];

    let mut count = 0;
    for &b1 in &bases {
        for &b2 in &bases {
            for &b3 in &bases {
                let codon = [b1, b2, b3];
                let aa = table.translate(&codon);
                assert_ne!(aa, b'X', "Codon {:?} should translate", std::str::from_utf8(&codon));
                count += 1;
            }
        }
    }
    assert_eq!(count, 64);
}

#[test]
fn test_codon_table_matches_ncbi() {
    use oxivep_genome::CodonTable;

    let table = CodonTable::standard();

    // NCBI standard genetic code verification
    assert_eq!(table.translate(b"ATG"), b'M'); // Met (start)
    assert_eq!(table.translate(b"TAA"), b'*'); // Stop
    assert_eq!(table.translate(b"TAG"), b'*'); // Stop
    assert_eq!(table.translate(b"TGA"), b'*'); // Stop

    // Verify a selection of other codons
    assert_eq!(table.translate(b"TTT"), b'F'); // Phe
    assert_eq!(table.translate(b"TTC"), b'F'); // Phe
    assert_eq!(table.translate(b"CTG"), b'L'); // Leu (most common)
    assert_eq!(table.translate(b"GAT"), b'D'); // Asp
    assert_eq!(table.translate(b"GAC"), b'D'); // Asp
    assert_eq!(table.translate(b"TGG"), b'W'); // Trp (only codon)
    assert_eq!(table.translate(b"CGT"), b'R'); // Arg
    assert_eq!(table.translate(b"AGA"), b'R'); // Arg
    assert_eq!(table.translate(b"GGG"), b'G'); // Gly
}
