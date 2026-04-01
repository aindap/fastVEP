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

### 1. Install Rust (if you don't have it)

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

### 2. Build and install OxiVEP

```bash
git clone https://github.com/kuanlinhuang/OxiVEP.git
cd OxiVEP

# Build and install to your PATH (~/.cargo/bin/oxivep)
cargo install --path crates/oxivep-cli

# Verify it works
oxivep --version
```

> **Note:** `cargo install` places the binary in `~/.cargo/bin/`. If `oxivep` is not found after install, run `source "$HOME/.cargo/env"` or add this line to your `~/.zshrc` (or `~/.bashrc`):
> ```bash
> source "$HOME/.cargo/env"
> ```

### 3. Try it — annotate the included test data

OxiVEP ships with a small test VCF and GFF3 so you can try it immediately:

```bash
# Annotate 12 test variants covering SNVs, indels, splice sites, UTRs, and intergenic regions
oxivep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format tab
```

Expected output:

```
rs001   chr1:500       G  ENSG00000001  ENST00000001  Transcript  upstream_gene_variant
rs002   chr1:1020      G  ENSG00000001  ENST00000001  Transcript  5_prime_UTR_variant
rs003   chr1:1100      G  ENSG00000001  ENST00000001  Transcript  missense_variant
rs004   chr1:1201      A  ENSG00000001  ENST00000001  Transcript  splice_donor_variant,intron_variant
rs005   chr1:1500      T  ENSG00000001  ENST00000001  Transcript  intron_variant
rs006   chr1:1999      G  ENSG00000001  ENST00000001  Transcript  splice_acceptor_variant,intron_variant
rs007   chr1:2100      G  ENSG00000001  ENST00000001  Transcript  missense_variant
rs008   chr1:4600      T  ENSG00000001  ENST00000001  Transcript  3_prime_UTR_variant
rs009   chr1:6000      C  ...           ...           Transcript  downstream_gene_variant
rs010   chr1:10100     G  ENSG00000002  ENST00000002  Transcript  non_coding_transcript_exon_variant
rs011   chr1:100000    G  -             -             Transcript  intergenic_variant
rs012   chr1:1101      -  ENSG00000001  ENST00000001  Transcript  frameshift_variant
```

The 12 variants cover: upstream, 5'UTR, missense, splice donor (HIGH), intron, splice acceptor (HIGH), 3'UTR, downstream, non-coding exon, intergenic, and frameshift.

### 4. Try VCF output (with CSQ annotations)

```bash
oxivep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format vcf
```

### 5. Try JSON output

```bash
oxivep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format json
```

### 6. Launch the web interface

```bash
oxivep web --port 8080
```

Open http://localhost:8080 in your browser. Click **"Load Example"** to load pre-built test variants, then click **"Annotate"** to see results instantly.

## Using Your Own Data

### Annotate a real VCF against Ensembl gene models

```bash
# Download Ensembl GFF3 (human GRCh38)
wget https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz
gunzip Homo_sapiens.GRCh38.112.gff3.gz

# Annotate your VCF
oxivep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 Homo_sapiens.GRCh38.112.gff3 \
  --hgvs

# Or pipe from bcftools
bcftools view sample.vcf.gz chr21 | oxivep annotate -i - --gff3 Homo_sapiens.GRCh38.112.gff3 --output-format tab
```

### Mouse, zebrafish, or other organisms

```bash
# Mouse
wget https://ftp.ensembl.org/pub/release-112/gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz

# Zebrafish
wget https://ftp.ensembl.org/pub/release-112/gff3/danio_rerio/Danio_rerio.GRCz11.112.gff3.gz
```

OxiVEP works with any organism — just provide the matching GFF3.

## Test Data Reference

The repository includes test files in `tests/`:

**`tests/test.gff3`** — Two genes on chr1:
- `BRCA1` (ENSG00000001): Protein-coding gene at chr1:1000-5000, 3 exons, CDS at 1050-4500
- `TP53` (ENSG00000002): lncRNA gene at chr1:10000-12000, 2 exons

**`tests/test.vcf`** — 12 variants designed to demonstrate every major consequence type:

| Variant | Position | Type | Expected Consequence |
|---------|----------|------|---------------------|
| rs001 | chr1:500 | SNV | upstream_gene_variant |
| rs002 | chr1:1020 | SNV (in 5'UTR) | 5_prime_UTR_variant |
| rs003 | chr1:1100 | SNV (in CDS) | missense_variant |
| rs004 | chr1:1201 | SNV (intron pos 1) | splice_donor_variant |
| rs005 | chr1:1500 | SNV (mid-intron) | intron_variant |
| rs006 | chr1:1999 | SNV (intron last base) | splice_acceptor_variant |
| rs007 | chr1:2100 | SNV (in CDS exon 2) | missense_variant |
| rs008 | chr1:4600 | SNV (in 3'UTR) | 3_prime_UTR_variant |
| rs009 | chr1:6000 | SNV (between genes) | downstream + upstream |
| rs010 | chr1:10100 | SNV (lncRNA exon) | non_coding_transcript_exon_variant |
| rs011 | chr1:100000 | SNV (far from genes) | intergenic_variant |
| rs012 | chr1:1100 | 1bp deletion (CDS) | frameshift_variant |

## Command Reference

### `oxivep annotate`

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

### Tab Output

One line per variant-transcript-allele combination with 13 columns.

### JSON Output

Structured JSON with `transcript_consequences` array per variant.

## Consequence Types

OxiVEP predicts 41 consequence types organized by impact:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant |
| **LOW** | splice_region_variant, splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, synonymous_variant, start_retained_variant, stop_retained_variant, incomplete_terminal_codon_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant, and others |

## Architecture

```
crates/
  oxivep-core/         # Core types: Consequence, GenomicPosition, Allele, Impact
  oxivep-genome/       # Transcript, Exon, Gene, CodonTable, coordinate mapping
  oxivep-cache/        # GFF3 parser, FASTA reader, annotation providers
  oxivep-consequence/  # Consequence prediction engine, splice site detection
  oxivep-hgvs/         # HGVS nomenclature generation (c., p., g.)
  oxivep-io/           # VCF parser, output formatters (CSQ, tab, JSON)
  oxivep-filter/       # Variant filtering
  oxivep-cli/          # CLI binary, annotation pipeline, web server
web/                   # Web GUI (HTML/CSS/JS)
tests/                 # Test VCF and GFF3 files
```

## Running Tests

```bash
cargo test --workspace          # 101 tests
cargo test -p oxivep-consequence  # Just consequence prediction tests
```

## Performance Benchmarks

Benchmarked on Apple M-series (ARM64), single-threaded, release build.

| Dataset | Variants | Time | Throughput |
|---------|----------|------|------------|
| Human chr21 (20 genes) | 1,000 | 0.037s | **27,000 variants/sec** |
| Human chr21 (20 genes) | 10,000 | 0.057s | **175,000 variants/sec** |
| Human chr21 (20 genes) | 50,000 | 0.128s | **391,000 variants/sec** |
| Mouse chr19 (10 genes) | 500 | 0.033s | **15,200 variants/sec** |
| Zebrafish chr5 (8 genes) | 300 | 0.031s | **9,800 variants/sec** |

### vs. Ensembl VEP

| Metric | Ensembl VEP (Perl) | OxiVEP (Rust) |
|--------|-------------------|---------------|
| Startup time | 5-15s | <0.05s |
| 1,000 SNVs (offline) | ~3-10s | **0.04s** |
| Peak memory (100K variants) | ~500 MB | **2.8 MB** |
| Binary size | ~200 MB installed | **2.4 MB** |
| Dependencies | Perl 5.22+, DBI, 10+ CPAN modules | **None** |

## License

Apache License 2.0

## Acknowledgements

OxiVEP is inspired by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) by EMBL-EBI. The consequence prediction logic follows the Sequence Ontology term definitions and the Ensembl variant annotation framework.
