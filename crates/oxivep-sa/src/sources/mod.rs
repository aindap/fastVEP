//! Source-specific parsers for building annotation databases.
//!
//! Each submodule implements a parser for a specific data source
//! that produces `AnnotationRecord`s for the `SaWriter`.

pub mod clinvar;
pub mod dbnsfp;
pub mod dbsnp;
pub mod gnomad;
pub mod primateai;
pub mod revel;
pub mod scores;
pub mod spliceai;
