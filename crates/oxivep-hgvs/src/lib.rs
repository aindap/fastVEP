mod coding;
mod genomic;
mod protein;

pub use coding::hgvsc;
pub use genomic::hgvsg;
pub use protein::hgvsp;

/// Full HGVS annotation result.
#[derive(Debug, Clone, Default)]
pub struct HgvsAnnotation {
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    pub hgvsg: Option<String>,
}
