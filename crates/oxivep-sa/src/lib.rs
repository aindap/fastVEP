//! OxiSA: Supplementary annotation format for OxiVEP.
//!
//! This crate provides binary reader/writer for three annotation file types:
//!
//! - **`.osa`** — Position/allele-level annotations (ClinVar, gnomAD, dbSNP, etc.)
//! - **`.osi`** — Interval-level annotations (SV databases, regulatory regions)
//! - **`.oga`** — Gene-level annotations (OMIM, pLI scores, ClinGen)
//!
//! All formats use zstd block compression with per-chromosome indexing for
//! fast random access during parallel annotation.

pub mod block;
pub mod common;
pub mod custom;
pub mod gene;
pub mod index;
pub mod interval;
pub mod kmer16;
pub mod reader;
pub mod sources;
pub mod var32;
pub mod writer;
pub mod zigzag;
