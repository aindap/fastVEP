//! ZigZag encoding for signed integers.
//!
//! Maps signed integers to unsigned integers so that values with small
//! absolute value have small encoded values, improving compression.
//! Used for signed scores (e.g., PhyloP: -20 to +30).
//!
//! Mapping: 0 → 0, -1 → 1, 1 → 2, -2 → 3, 2 → 4, ...

/// Encode a signed i32 as an unsigned u32 using zigzag encoding.
#[inline]
pub fn encode(val: i32) -> u32 {
    ((val << 1) ^ (val >> 31)) as u32
}

/// Decode a zigzag-encoded u32 back to i32.
#[inline]
pub fn decode(val: u32) -> i32 {
    ((val >> 1) as i32) ^ -((val & 1) as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round_trip() {
        for v in [-1000, -100, -1, 0, 1, 100, 1000, i32::MIN, i32::MAX] {
            assert_eq!(decode(encode(v)), v, "failed for {}", v);
        }
    }

    #[test]
    fn test_small_values_encode_small() {
        assert_eq!(encode(0), 0);
        assert_eq!(encode(-1), 1);
        assert_eq!(encode(1), 2);
        assert_eq!(encode(-2), 3);
        assert_eq!(encode(2), 4);
    }

    #[test]
    fn test_negative_values() {
        assert_eq!(decode(encode(-20)), -20);
        assert_eq!(decode(encode(-12345)), -12345);
    }
}
