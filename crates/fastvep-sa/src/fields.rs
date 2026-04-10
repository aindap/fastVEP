//! Field type definitions for the .osa2 format.
//!
//! Each annotation field (e.g., gnomAD AF, ClinVar significance) has a type
//! that determines how it's stored as a u32 integer and reconstructed for output.

use serde::{Deserialize, Serialize};

/// How a field's values are stored internally.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FieldType {
    /// Raw integer (u32). Optional zigzag for signed values.
    Integer,
    /// Float scaled to integer: `stored = (value * multiplier) as u32`.
    Float,
    /// String mapped to u32 index via a lookup table.
    Categorical,
    /// Boolean flag: 0 or 1.
    Flag,
    /// Opaque JSON string stored separately (not in u32 arrays).
    JsonBlob,
}

/// Definition of a single annotation field within an .osa2 file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Field {
    /// Original VCF INFO field name (e.g., "AF").
    pub field: String,
    /// Output alias (e.g., "gnomad_af"). Used as the JSON key.
    pub alias: String,
    /// How the field is stored.
    pub ftype: FieldType,
    /// For Float fields: multiply by this before casting to u32.
    /// E.g., 2_000_000 for allele frequencies (6 significant digits).
    #[serde(default = "default_multiplier")]
    pub multiplier: u32,
    /// Whether to apply zigzag encoding (for signed values).
    #[serde(default)]
    pub zigzag: bool,
    /// Sentinel value stored when the field is missing. Default: u32::MAX.
    #[serde(default = "default_missing")]
    pub missing_value: u32,
    /// For Categorical fields: the string to use for missing values.
    #[serde(default = "default_missing_string")]
    pub missing_string: String,
    /// Human-readable description.
    #[serde(default)]
    pub description: String,
}

fn default_multiplier() -> u32 { 1 }
fn default_missing() -> u32 { u32::MAX }
fn default_missing_string() -> String { ".".into() }

impl Field {
    /// Encode a float value as a u32 using this field's multiplier.
    #[inline]
    pub fn encode_float(&self, value: f64) -> u32 {
        let scaled = value * self.multiplier as f64;
        if self.zigzag {
            crate::zigzag::encode(scaled as i32)
        } else {
            scaled as u32
        }
    }

    /// Decode a stored u32 back to a float.
    #[inline]
    pub fn decode_float(&self, stored: u32) -> f64 {
        if stored == self.missing_value {
            return f64::NAN;
        }
        let raw = if self.zigzag {
            crate::zigzag::decode(stored) as f64
        } else {
            stored as f64
        };
        raw / self.multiplier as f64
    }

    /// Encode an integer value, optionally with zigzag.
    #[inline]
    pub fn encode_int(&self, value: i64) -> u32 {
        if self.zigzag {
            crate::zigzag::encode(value as i32)
        } else {
            value as u32
        }
    }

    /// Decode a stored u32 back to an integer.
    #[inline]
    pub fn decode_int(&self, stored: u32) -> i64 {
        if stored == self.missing_value {
            return 0;
        }
        if self.zigzag {
            crate::zigzag::decode(stored) as i64
        } else {
            stored as i64
        }
    }
}

/// Format a decoded value as a JSON fragment string.
pub fn format_value(field: &Field, stored: u32, strings: Option<&[String]>) -> String {
    if stored == field.missing_value {
        return "null".into();
    }
    match field.ftype {
        FieldType::Float => {
            let val = field.decode_float(stored);
            if val.abs() < 1e-10 { "0".into() }
            else { format!("{:.6e}", val) }
        }
        FieldType::Integer => {
            let val = field.decode_int(stored);
            format!("{}", val)
        }
        FieldType::Categorical => {
            if let Some(strs) = strings {
                let idx = stored as usize;
                if idx < strs.len() {
                    format!("\"{}\"", strs[idx])
                } else {
                    "null".into()
                }
            } else {
                format!("{}", stored)
            }
        }
        FieldType::Flag => {
            if stored != 0 { "true".into() } else { "false".into() }
        }
        FieldType::JsonBlob => {
            // JsonBlob values are stored separately, not in u32 arrays
            "null".into()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_float_round_trip() {
        let field = Field {
            field: "AF".into(), alias: "af".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let original = 0.001234;
        let stored = field.encode_float(original);
        let decoded = field.decode_float(stored);
        assert!((decoded - original).abs() < 1e-6);
    }

    #[test]
    fn test_zigzag_float() {
        let field = Field {
            field: "score".into(), alias: "score".into(), ftype: FieldType::Float,
            multiplier: 10_000, zigzag: true, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let original = -2.5;
        let stored = field.encode_float(original);
        let decoded = field.decode_float(stored);
        assert!((decoded - original).abs() < 0.001);
    }

    #[test]
    fn test_missing_value() {
        let field = Field {
            field: "AF".into(), alias: "af".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        assert!(field.decode_float(u32::MAX).is_nan());
        assert_eq!(format_value(&field, u32::MAX, None), "null");
    }

    #[test]
    fn test_categorical() {
        let field = Field {
            field: "sig".into(), alias: "sig".into(), ftype: FieldType::Categorical,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let strings = vec!["Benign".to_string(), "Pathogenic".to_string()];
        assert_eq!(format_value(&field, 0, Some(&strings)), "\"Benign\"");
        assert_eq!(format_value(&field, 1, Some(&strings)), "\"Pathogenic\"");
    }
}
