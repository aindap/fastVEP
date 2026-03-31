# OxiVEP

A high-performance Variant Effect Predictor written in Rust. OxiVEP predicts the functional consequences of genomic variants (SNPs, insertions, deletions) on genes, transcripts, and protein sequences.

OxiVEP is inspired by and aims to be compatible with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), while delivering significantly better performance through Rust's zero-cost abstractions and native parallelism.

## Features

- **Variant Consequence Prediction** — Classifies variants using [Sequence Ontology](http://www.sequenceontology.org/) terms (missense, frameshift, splice donor, etc.)
- **HGVS Nomenclature** — Generates HGVSg, HGVSc, and HGVSp notations
- **Multiple Output Formats** — VCF (with CSQ field), tab-delimited, and JSON
- **GFF3 Annotation Support** — Load gene models from standard GFF3 files
- **Splice Site Detection** — Identifies splice donor, acceptor, region, 5th base, donor region, and polypyrimidine tract variants
- **Multi-allelic Support** — Handles multi-allelic VCF records with proper allele normalization
- **Web Interface** — Built-in web GUI for interactive variant annotation

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/OxiVEP.git
cd OxiVEP

# Build (requires Rust 1.70+)
cargo build --release

# The binary is at target/release/oxivep
```

### Basic Usage

```bash
# Annotate a VCF file using a GFF3 gene model
oxivep annotate \
  -i variants.vcf \
  -o annotated.vcf \
  --gff3 genes.gff3 \
  --hgvs

# Output as tab-delimited
oxivep annotate \
  -i variants.vcf \
  --gff3 genes.gff3 \
  --output-format tab

# Output as JSON
oxivep annotate \
  -i variants.vcf \
  --gff3 genes.gff3 \
  --output-format json

# Launch the web interface
oxivep web --gff3 genes.gff3 --port 8080
```

### Reading from stdin

```bash
cat variants.vcf | oxivep annotate -i - --gff3 genes.gff3
```

## Command Reference

### `oxivep annotate`

Annotate variants with predicted functional consequences.

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VCF file (`-` for stdin) | *required* |
| `-o, --output` | Output file (`-` for stdout) | `-` |
| `--gff3` | GFF3 gene annotation file | — |
| `--fasta` | Reference FASTA file | — |
| `--output-format` | `vcf`, `tab`, or `json` | `vcf` |
| `--hgvs` | Include HGVS notations | off |
| `--pick` | Report only the most severe consequence per variant | off |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--everything` | Enable all annotation flags | off |
| `--symbol` | Include gene symbol | off |
| `--canonical` | Flag canonical transcripts | off |

### `oxivep web`

Launch the interactive web interface.

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 gene annotation file | — |
| `--fasta` | Reference FASTA file | — |
| `--port` | HTTP port | `8080` |

### `oxivep filter`

Filter annotated VEP output (coming soon).

## Output Formats

### VCF Output

Annotations are added as a `CSQ` field in the INFO column:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from OxiVEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND">
```

Example:
```
chr1  1100  rs003  A  G  30  PASS  CSQ=G|missense_variant|MODERATE|BRCA1|ENSG001|Transcript|ENST001|protein_coding|1/3||ENST001:c.51A>G||101|51|17|||||1
```

### Tab Output

One line per variant-transcript-allele combination:

```
#Uploaded_variation  Location  Allele  Gene  Feature  Feature_type  Consequence  cDNA_position  CDS_position  Protein_position  Amino_acids  Codons  Existing_variation
rs003  chr1:1100  G  ENSG001  ENST001  Transcript  missense_variant  101  51  17
```

### JSON Output

```json
{
  "id": "rs003",
  "seq_region_name": "chr1",
  "start": 1100,
  "end": 1100,
  "allele_string": "A/G",
  "most_severe_consequence": "missense_variant",
  "transcript_consequences": [
    {
      "gene_id": "ENSG001",
      "transcript_id": "ENST001",
      "consequence_terms": ["missense_variant"],
      "impact": "MODERATE",
      "variant_allele": "G"
    }
  ]
}
```

## Consequence Types

OxiVEP predicts 41 consequence types organized by impact:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant |
| **LOW** | splice_region_variant, splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, synonymous_variant, start_retained_variant, stop_retained_variant, incomplete_terminal_codon_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant, and others |

## Input Requirements

### VCF Files

Standard VCF 4.1+ format. OxiVEP handles:
- SNVs, MNVs, insertions, deletions
- Multi-allelic sites (comma-separated ALT)
- Star alleles (`*`)
- Chromosome prefixes (`chr1` or `1`)
- Phased and unphased genotypes

### GFF3 Files

Standard GFF3 format from Ensembl, NCBI, or GENCODE. Required feature types:
- `gene` — Gene boundaries and symbols
- `mRNA` / `transcript` — Transcript models
- `exon` — Exon coordinates
- `CDS` — Coding sequence regions

Download from Ensembl:
```bash
# Human GRCh38
wget https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz
gunzip Homo_sapiens.GRCh38.112.gff3.gz
```

## Architecture

OxiVEP is organized as a Cargo workspace with modular crates:

```
crates/
  oxivep-core/         # Core types: Consequence, GenomicPosition, Allele, Impact
  oxivep-genome/       # Transcript, Exon, Gene, CodonTable, coordinate mapping
  oxivep-cache/        # GFF3 parser, FASTA reader, annotation providers
  oxivep-consequence/  # Consequence prediction engine, splice site detection
  oxivep-hgvs/         # HGVS nomenclature generation (c., p., g.)
  oxivep-io/           # VCF parser, output formatters (CSQ, tab, JSON)
  oxivep-filter/       # Variant filtering
  oxivep-cli/          # CLI binary and annotation pipeline
web/                   # Web GUI (HTML/CSS/JS)
```

## Running Tests

```bash
# Run all tests
cargo test --workspace

# Run tests for a specific crate
cargo test -p oxivep-consequence

# Run with output
cargo test --workspace -- --nocapture
```

## Comparison with Ensembl VEP

| Feature | Ensembl VEP | OxiVEP |
|---------|-------------|--------|
| Language | Perl | Rust |
| Annotation source | Cache (Storable), Database, GFF3 | GFF3, FASTA |
| Input formats | VCF, HGVS, IDs, SPDI, Region | VCF |
| Output formats | VCF, Tab, JSON, VEP native | VCF, Tab, JSON |
| Parallelism | Fork-based | Thread-based (rayon) |
| Consequence types | ~35 SO terms | 41 SO terms |
| HGVS | Full (c., p., g.) | Full (c., p., g.) |
| Splice detection | Donor, acceptor, region | Donor, acceptor, region, 5th base, polypyrimidine |
| Web interface | REST API | Built-in GUI |
| Dependencies | Perl 5.22+, DBI, many CPAN modules | None (static binary) |

## Performance Benchmarks

Benchmarked on Apple M-series (ARM64), single-threaded, release build.

### Throughput

| Dataset | Variants | Time | Throughput |
|---------|----------|------|------------|
| Human chr21 (20 genes) | 1,000 | 0.037s | **27,000 variants/sec** |
| Human chr21 (20 genes) | 10,000 | 0.057s | **175,000 variants/sec** |
| Human chr21 (20 genes) | 50,000 | 0.128s | **391,000 variants/sec** |
| Mouse chr19 (10 genes) | 500 | 0.033s | **15,200 variants/sec** |
| Zebrafish chr5 (8 genes) | 300 | 0.031s | **9,800 variants/sec** |

### Comparison with Ensembl VEP

| Metric | Ensembl VEP (Perl) | OxiVEP (Rust) |
|--------|-------------------|---------------|
| Startup time | 5-15s (loading modules) | <0.01s |
| 1,000 SNVs (offline, cache) | ~3-10s | **0.04s** |
| Memory (1K variants) | ~500MB | **~5MB** |
| Dependencies | Perl 5.22+, DBI, 10+ CPAN modules | **None (static binary)** |

### Annotation Accuracy

Tested against Ensembl VEP test suite patterns (101 tests, 22 VEP compatibility tests):

| Category | Tests | Status |
|----------|-------|--------|
| VCF parsing (SNV, indel, multi-allelic, MNV) | 22 | All pass |
| Allele normalization (VEP-compatible) | 8 | All pass |
| Consequence prediction (all SO terms) | 21 | All pass |
| HGVS nomenclature (c., p., g.) | 14 | All pass |
| Codon table (NCBI standard genetic code) | 6 | All pass |
| GFF3 parsing (forward/reverse strand) | 6 | All pass |
| Cache/FASTA reading | 12 | All pass |
| Output formatting (CSQ, tab, JSON) | 12 | All pass |

Consequence distribution from 1,000 realistic human chr21 variants:
```
859 intergenic_variant
 82 intron_variant
 66 missense_variant
 15 3_prime_UTR_variant
 14 splice_donor_variant
 14 5_prime_UTR_variant
 10 splice_acceptor_variant
  7 upstream_gene_variant
  4 non_coding_transcript_variant
  2 downstream_gene_variant
```

## License

Apache License 2.0

## Acknowledgements

OxiVEP is inspired by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) by EMBL-EBI. The consequence prediction logic follows the Sequence Ontology term definitions and the Ensembl variant annotation framework.
