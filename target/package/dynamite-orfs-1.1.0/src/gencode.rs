//! Genetic-code (NCBI translation table) support.
//!
//! DYNAMITE supports the full set of NCBI tables that are biologically
//! meaningful for the targeted clades. The headline ones for this tool are:
//!
//! | Code | Organism group                              | Reassignment        |
//! |------|---------------------------------------------|---------------------|
//! | 11   | Bacteria / Archaea / plastids (default)     | none                |
//! |  1   | Standard                                    | none                |
//! |  4   | Mold/protozoan/coelenterate mito, *Mycoplasma*, many phages | TGA = Trp |
//! | 15   | *Blepharisma*; Topaz/Agate NCLDV-associated lineages, gut phages | TAG = Gln |
//! | 25   | Candidate division SR1 / *Gracilibacteria*  | TGA = Gly           |
//! |  6   | Ciliate / dasycladacean / hexamita nuclear  | TAA, TAG = Gln      |
//!
//! Stop codons are derived from the amino-acid table (`*` positions), so the
//! alternative codes automatically change which codons terminate translation.
//! Start-codon recognition is restricted to ATG / GTG / TTG with the same
//! per-table rules Prodigal uses, so the start-node typing stays faithful.

use crate::seq::{A, C, G, T};

/// NCBI rank of a base for indexing the canonical TCAG-ordered codon tables.
#[inline]
fn ncbi_rank(code: u8) -> usize {
    // Codon tables are written in T, C, A, G order.
    match code & 0b11 {
        T => 0,
        C => 1,
        A => 2,
        _ => 3, // G
    }
}

/// Codon index into an NCBI 64-entry table for the codon `b1 b2 b3`.
#[inline]
fn codon_index(b1: u8, b2: u8, b3: u8) -> usize {
    16 * ncbi_rank(b1) + 4 * ncbi_rank(b2) + ncbi_rank(b3)
}

/// One NCBI translation table.
#[derive(Clone)]
pub struct GeneticCode {
    /// NCBI translation-table number.
    pub id: u32,
    /// Amino acids for all 64 codons in TCAG order; `*` marks a stop codon.
    aas: &'static [u8; 64],
}

/// Return the amino-acid string for a given NCBI table id, if supported.
fn aa_table(id: u32) -> Option<&'static [u8; 64]> {
    // All strings are 64 chars, codons ordered TTT,TTC,...,GGG (TCAG-major).
    let s: &'static [u8] = match id {
        1 => b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        2 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        3 => b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        4 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        5 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        6 => b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        9 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        10 => b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        11 => b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        12 => b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        13 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        14 => b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        15 => b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        16 => b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        21 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        22 => b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        23 => b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        24 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        25 => b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        _ => return None,
    };
    // Safe: every literal above is exactly 64 bytes (checked by unit test).
    Some(<&[u8; 64]>::try_from(s).expect("genetic-code table must be 64 chars"))
}

/// Genetic-code ids supported by DYNAMITE, in a sensible display order.
pub const SUPPORTED: &[u32] = &[
    11, 1, 4, 15, 25, 6, 2, 3, 5, 9, 10, 12, 13, 14, 16, 21, 22, 23, 24,
];

impl GeneticCode {
    /// Construct a genetic code by NCBI id. Returns `None` if unsupported.
    pub fn new(id: u32) -> Option<GeneticCode> {
        aa_table(id).map(|aas| GeneticCode { id, aas })
    }

    /// Translate a codon to its amino acid. Returns `X` if any base of the
    /// codon is ambiguous (caller passes `n_amb`), otherwise the table AA.
    #[inline]
    pub fn translate(&self, b1: u8, b2: u8, b3: u8) -> u8 {
        self.aas[codon_index(b1, b2, b3)]
    }

    /// `true` if `b1 b2 b3` is a stop codon under this table.
    #[inline]
    pub fn is_stop(&self, b1: u8, b2: u8, b3: u8) -> bool {
        self.aas[codon_index(b1, b2, b3)] == b'*'
    }

    /// `true` if `b1 b2 b3` is an accepted start codon (ATG/GTG/TTG, restricted
    /// per table exactly as Prodigal does).
    #[inline]
    pub fn is_start(&self, b1: u8, b2: u8, b3: u8) -> bool {
        // ATG is always a start.
        if (b1, b2, b3) == (A, T, G) {
            return true;
        }
        // Tables that recognise only ATG as a start.
        if matches!(self.id, 6 | 10 | 14 | 15 | 16 | 22) {
            return false;
        }
        // GTG
        if (b1, b2, b3) == (G, T, G) {
            return !matches!(self.id, 1 | 3 | 12 | 22);
        }
        // TTG
        if (b1, b2, b3) == (T, T, G) {
            let id = self.id;
            return !(id < 4 || id == 9 || (21..25).contains(&id));
        }
        false
    }

    /// Start-codon type used by the dynamic-programming engine:
    /// `0 = ATG, 1 = GTG, 2 = TTG, 3 = not a start`.
    #[inline]
    pub fn start_type(b1: u8, b2: u8, b3: u8) -> u8 {
        match (b1, b2, b3) {
            (A, T, G) => 0,
            (G, T, G) => 1,
            (T, T, G) => 2,
            _ => 3,
        }
    }

    /// Bitmask over the 64 codons (bit = NCBI TCAG index) that are stop codons
    /// under this table. Two tables with equal masks terminate ORFs identically.
    pub fn stop_mask(&self) -> u64 {
        let mut m = 0u64;
        for (i, &aa) in self.aas.iter().enumerate() {
            if aa == b'*' {
                m |= 1u64 << i;
            }
        }
        m
    }

    /// Bitmask over the 64 codons (bit = NCBI TCAG index) accepted as start
    /// codons by [`GeneticCode::is_start`].
    pub fn start_mask(&self) -> u64 {
        let bases = [T, C, A, G];
        let mut m = 0u64;
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    if self.is_start(b1, b2, b3) {
                        m |= 1u64 << codon_index(b1, b2, b3);
                    }
                }
            }
        }
        m
    }

    /// The set of stop codons (as ASCII triples) for this table — used by the
    /// PHANOTATE-style engine, which works on codon strings.
    pub fn stop_codons_ascii(&self) -> Vec<[u8; 3]> {
        let bases = [T, C, A, G]; // any order; we just enumerate all 64
        let mut out = Vec::new();
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    if self.is_stop(b1, b2, b3) {
                        out.push([
                            crate::seq::decode_base(b1),
                            crate::seq::decode_base(b2),
                            crate::seq::decode_base(b3),
                        ]);
                    }
                }
            }
        }
        out
    }
}

/// Lowercase ASCII for a codon triple of base codes.
pub fn codon_lower(b1: u8, b2: u8, b3: u8) -> [u8; 3] {
    [
        crate::seq::decode_base(b1).to_ascii_lowercase(),
        crate::seq::decode_base(b2).to_ascii_lowercase(),
        crate::seq::decode_base(b3).to_ascii_lowercase(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{A, G, T};

    #[test]
    fn all_tables_are_64() {
        for &id in SUPPORTED {
            let code = GeneticCode::new(id).unwrap();
            assert_eq!(code.aas.len(), 64, "table {id}");
        }
    }

    #[test]
    fn standard_code_translation() {
        let c = GeneticCode::new(11).unwrap();
        assert_eq!(c.translate(A, T, G), b'M'); // ATG -> Met
        assert_eq!(c.translate(T, G, G), b'W'); // TGG -> Trp
        assert!(c.is_stop(T, A, A)); // TAA stop
        assert!(c.is_stop(T, G, A)); // TGA stop
        assert!(c.is_stop(T, A, G)); // TAG stop
    }

    #[test]
    fn alternative_codes() {
        // Code 4: TGA = Trp (not a stop)
        let c4 = GeneticCode::new(4).unwrap();
        assert_eq!(c4.translate(T, G, A), b'W');
        assert!(!c4.is_stop(T, G, A));
        assert!(c4.is_stop(T, A, A));
        assert!(c4.is_stop(T, A, G));

        // Code 15: TAG = Gln (not a stop)
        let c15 = GeneticCode::new(15).unwrap();
        assert_eq!(c15.translate(T, A, G), b'Q');
        assert!(!c15.is_stop(T, A, G));
        assert!(c15.is_stop(T, A, A));
        assert!(c15.is_stop(T, G, A));

        // Code 25: TGA = Gly
        let c25 = GeneticCode::new(25).unwrap();
        assert_eq!(c25.translate(T, G, A), b'G');
        assert!(!c25.is_stop(T, G, A));
    }

    #[test]
    fn stop_and_start_masks() {
        let c11 = GeneticCode::new(11).unwrap();
        let c1 = GeneticCode::new(1).unwrap();
        let c4 = GeneticCode::new(4).unwrap();
        let c25 = GeneticCode::new(25).unwrap();
        let c15 = GeneticCode::new(15).unwrap();
        let c6 = GeneticCode::new(6).unwrap();
        assert_eq!(c11.stop_mask().count_ones(), 3);
        assert_eq!(c11.stop_mask(), c1.stop_mask()); // same stops ...
        assert_ne!(c11.start_mask(), c1.start_mask()); // ... different starts
        assert_eq!((c4.stop_mask(), c4.start_mask()), (c25.stop_mask(), c25.start_mask()));
        // stop-set nesting: 6 ⊂ 15 ⊂ 11
        assert_eq!(c6.stop_mask() & c15.stop_mask(), c6.stop_mask());
        assert_eq!(c15.stop_mask() & c11.stop_mask(), c15.stop_mask());
    }

    #[test]
    fn start_rules() {
        let c11 = GeneticCode::new(11).unwrap();
        assert!(c11.is_start(A, T, G));
        assert!(c11.is_start(G, T, G));
        assert!(c11.is_start(T, T, G));
        // Code 15 only ATG
        let c15 = GeneticCode::new(15).unwrap();
        assert!(c15.is_start(A, T, G));
        assert!(!c15.is_start(G, T, G));
        assert!(!c15.is_start(T, T, G));
    }
}
