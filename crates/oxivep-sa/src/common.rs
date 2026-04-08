//! Constants and common types for the OxiSA binary annotation format.

/// Magic bytes for position/allele-level annotation files (.osa).
pub const OSA_MAGIC: &[u8; 8] = b"OXISA_01";

/// Magic bytes for interval-level annotation files (.osi).
pub const OSI_MAGIC: &[u8; 8] = b"OXISI_01";

/// Magic bytes for gene-level annotation files (.oga).
pub const OGA_MAGIC: &[u8; 8] = b"OXIGA_01";

/// Current schema version. Bump when the binary format changes.
pub const SCHEMA_VERSION: u16 = 1;

/// Default block size for compression (8 MiB).
pub const DEFAULT_BLOCK_SIZE: usize = 8 * 1024 * 1024;

/// Default zstd compression level (3 is a good speed/ratio tradeoff).
pub const ZSTD_LEVEL: i32 = 3;

/// File extension for position/allele-level annotations.
pub const OSA_EXT: &str = "osa";

/// File extension for the index file.
pub const IDX_EXT: &str = "osa.idx";

/// File extension for interval-level annotations.
pub const OSI_EXT: &str = "osi";

/// File extension for gene-level annotations.
pub const OGA_EXT: &str = "oga";

/// A single annotation record ready for writing.
#[derive(Debug, Clone)]
pub struct AnnotationRecord {
    /// Chromosome index (numeric, mapped externally).
    pub chrom_idx: u16,
    /// 1-based genomic position.
    pub position: u32,
    /// Reference allele (empty string for positional annotations).
    pub ref_allele: String,
    /// Alternate allele (empty string for positional annotations).
    pub alt_allele: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// A single interval annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IntervalRecord {
    /// Chromosome name.
    pub chrom: String,
    /// 1-based start position (inclusive).
    pub start: u32,
    /// 1-based end position (inclusive).
    pub end: u32,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// A single gene annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GeneRecord {
    /// Gene symbol (e.g., "BRCA1").
    pub gene_symbol: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}
