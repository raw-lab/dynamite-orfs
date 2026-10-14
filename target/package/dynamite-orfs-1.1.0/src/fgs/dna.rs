use std::cmp::{max, min};

pub const CG_MIN: usize = 26;
pub const CG_MAX: usize = 70;

pub const ACGT: usize = 4;
pub const BI_ACGT: usize = 4 * 4;
pub const TRI_ACGT: usize = 4 * 4 * 4;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Nuc {
    A,
    C,
    G,
    T,
    N,
    Ai,
    Ci,
    Gi,
    Ti,
    Ni,
    Xi,
}

impl Nuc {
    pub fn to_int(&self) -> Option<usize> {
        match self {
            Nuc::A | Nuc::Ai => Some(0),
            Nuc::C | Nuc::Ci => Some(1),
            Nuc::G | Nuc::Gi => Some(2),
            Nuc::T | Nuc::Ti => Some(3),
            Nuc::N | Nuc::Ni => None,
            Nuc::Xi => None,
        }
    }

    pub fn to_lower(&self) -> Nuc {
        match self {
            Nuc::A | Nuc::Ai => Nuc::Ai,
            Nuc::C | Nuc::Ci => Nuc::Ci,
            Nuc::G | Nuc::Gi => Nuc::Gi,
            Nuc::T | Nuc::Ti => Nuc::Ti,
            Nuc::N | Nuc::Ni => Nuc::Ni,
            Nuc::Xi => Nuc::Xi,
        }
    }

    pub fn is_insertion(&self) -> bool {
        match self {
            Nuc::Ai | Nuc::Ci | Nuc::Gi | Nuc::Ti | Nuc::Ni => true,
            _ => false,
        }
    }

    pub fn rc(&self) -> Nuc {
        match self {
            Nuc::A => Nuc::T,
            Nuc::C => Nuc::G,
            Nuc::G => Nuc::C,
            Nuc::T => Nuc::A,
            Nuc::N => Nuc::N,
            Nuc::Ai => Nuc::Ti,
            Nuc::Ci => Nuc::Gi,
            Nuc::Gi => Nuc::Ci,
            Nuc::Ti => Nuc::Ai,
            Nuc::Ni => Nuc::Ni,
            Nuc::Xi => Nuc::Xi,
        }
    }
}

impl From<u8> for Nuc {
    fn from(nt: u8) -> Nuc {
        match nt {
            b'A' => Nuc::A,
            b'C' => Nuc::C,
            b'G' => Nuc::G,
            b'T' => Nuc::T,
            b'a' => Nuc::Ai,
            b'c' => Nuc::Ci,
            b'g' => Nuc::Gi,
            b't' => Nuc::Ti,
            b'x' => Nuc::Xi,
            b'n' => Nuc::Ni,
            _ => Nuc::N,
        }
    }
}

impl From<Nuc> for u8 {
    fn from(nt: Nuc) -> u8 {
        match nt {
            Nuc::A => b'A',
            Nuc::C => b'C',
            Nuc::G => b'G',
            Nuc::T => b'T',
            Nuc::N => b'N',
            Nuc::Ai => b'a',
            Nuc::Ci => b'c',
            Nuc::Gi => b'g',
            Nuc::Ti => b't',
            Nuc::Ni => b'n',
            Nuc::Xi => b'N',
        }
    }
}

pub fn dna(seq: &str) -> Vec<Nuc> {
    seq.as_bytes().into_iter().copied().map(Nuc::from).collect()
}

pub const CODON_CODE: [u8; TRI_ACGT] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S', b'I', b'I', b'M', b'I',
    b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P', b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L',
    b'E', b'D', b'E', b'D', b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C', b'L', b'F', b'L', b'F',
];

pub const ANTI_CODON_CODE: [u8; TRI_ACGT] = [
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'M', b'W', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'I', b'*', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
];

pub fn trinucleotide(n: &[Nuc]) -> Option<usize> {
    match n.get(0..3).unwrap() {
        [Nuc::A, Nuc::A, Nuc::A]
        | [Nuc::A, Nuc::A, Nuc::Ai]
        | [Nuc::A, Nuc::Ai, Nuc::A]
        | [Nuc::A, Nuc::Ai, Nuc::Ai]
        | [Nuc::Ai, Nuc::A, Nuc::A]
        | [Nuc::Ai, Nuc::A, Nuc::Ai]
        | [Nuc::Ai, Nuc::Ai, Nuc::A]
        | [Nuc::Ai, Nuc::Ai, Nuc::Ai] => Some(0),
        [Nuc::A, Nuc::A, Nuc::C]
        | [Nuc::A, Nuc::A, Nuc::Ci]
        | [Nuc::A, Nuc::Ai, Nuc::C]
        | [Nuc::A, Nuc::Ai, Nuc::Ci]
        | [Nuc::Ai, Nuc::A, Nuc::C]
        | [Nuc::Ai, Nuc::A, Nuc::Ci]
        | [Nuc::Ai, Nuc::Ai, Nuc::C]
        | [Nuc::Ai, Nuc::Ai, Nuc::Ci] => Some(1),
        [Nuc::A, Nuc::A, Nuc::G]
        | [Nuc::A, Nuc::A, Nuc::Gi]
        | [Nuc::A, Nuc::Ai, Nuc::G]
        | [Nuc::A, Nuc::Ai, Nuc::Gi]
        | [Nuc::Ai, Nuc::A, Nuc::G]
        | [Nuc::Ai, Nuc::A, Nuc::Gi]
        | [Nuc::Ai, Nuc::Ai, Nuc::G]
        | [Nuc::Ai, Nuc::Ai, Nuc::Gi] => Some(2),
        [Nuc::A, Nuc::A, Nuc::T]
        | [Nuc::A, Nuc::A, Nuc::Ti]
        | [Nuc::A, Nuc::Ai, Nuc::T]
        | [Nuc::A, Nuc::Ai, Nuc::Ti]
        | [Nuc::Ai, Nuc::A, Nuc::T]
        | [Nuc::Ai, Nuc::A, Nuc::Ti]
        | [Nuc::Ai, Nuc::Ai, Nuc::T]
        | [Nuc::Ai, Nuc::Ai, Nuc::Ti] => Some(3),
        [Nuc::A, Nuc::C, Nuc::A]
        | [Nuc::A, Nuc::C, Nuc::Ai]
        | [Nuc::A, Nuc::Ci, Nuc::A]
        | [Nuc::A, Nuc::Ci, Nuc::Ai]
        | [Nuc::Ai, Nuc::C, Nuc::A]
        | [Nuc::Ai, Nuc::C, Nuc::Ai]
        | [Nuc::Ai, Nuc::Ci, Nuc::A]
        | [Nuc::Ai, Nuc::Ci, Nuc::Ai] => Some(4),
        [Nuc::A, Nuc::C, Nuc::C]
        | [Nuc::A, Nuc::C, Nuc::Ci]
        | [Nuc::A, Nuc::Ci, Nuc::C]
        | [Nuc::A, Nuc::Ci, Nuc::Ci]
        | [Nuc::Ai, Nuc::C, Nuc::C]
        | [Nuc::Ai, Nuc::C, Nuc::Ci]
        | [Nuc::Ai, Nuc::Ci, Nuc::C]
        | [Nuc::Ai, Nuc::Ci, Nuc::Ci] => Some(5),
        [Nuc::A, Nuc::C, Nuc::G]
        | [Nuc::A, Nuc::C, Nuc::Gi]
        | [Nuc::A, Nuc::Ci, Nuc::G]
        | [Nuc::A, Nuc::Ci, Nuc::Gi]
        | [Nuc::Ai, Nuc::C, Nuc::G]
        | [Nuc::Ai, Nuc::C, Nuc::Gi]
        | [Nuc::Ai, Nuc::Ci, Nuc::G]
        | [Nuc::Ai, Nuc::Ci, Nuc::Gi] => Some(6),
        [Nuc::A, Nuc::C, Nuc::T]
        | [Nuc::A, Nuc::C, Nuc::Ti]
        | [Nuc::A, Nuc::Ci, Nuc::T]
        | [Nuc::A, Nuc::Ci, Nuc::Ti]
        | [Nuc::Ai, Nuc::C, Nuc::T]
        | [Nuc::Ai, Nuc::C, Nuc::Ti]
        | [Nuc::Ai, Nuc::Ci, Nuc::T]
        | [Nuc::Ai, Nuc::Ci, Nuc::Ti] => Some(7),
        [Nuc::A, Nuc::G, Nuc::A]
        | [Nuc::A, Nuc::G, Nuc::Ai]
        | [Nuc::A, Nuc::Gi, Nuc::A]
        | [Nuc::A, Nuc::Gi, Nuc::Ai]
        | [Nuc::Ai, Nuc::G, Nuc::A]
        | [Nuc::Ai, Nuc::G, Nuc::Ai]
        | [Nuc::Ai, Nuc::Gi, Nuc::A]
        | [Nuc::Ai, Nuc::Gi, Nuc::Ai] => Some(8),
        [Nuc::A, Nuc::G, Nuc::C]
        | [Nuc::A, Nuc::G, Nuc::Ci]
        | [Nuc::A, Nuc::Gi, Nuc::C]
        | [Nuc::A, Nuc::Gi, Nuc::Ci]
        | [Nuc::Ai, Nuc::G, Nuc::C]
        | [Nuc::Ai, Nuc::G, Nuc::Ci]
        | [Nuc::Ai, Nuc::Gi, Nuc::C]
        | [Nuc::Ai, Nuc::Gi, Nuc::Ci] => Some(9),
        [Nuc::A, Nuc::G, Nuc::G]
        | [Nuc::A, Nuc::G, Nuc::Gi]
        | [Nuc::A, Nuc::Gi, Nuc::G]
        | [Nuc::A, Nuc::Gi, Nuc::Gi]
        | [Nuc::Ai, Nuc::G, Nuc::G]
        | [Nuc::Ai, Nuc::G, Nuc::Gi]
        | [Nuc::Ai, Nuc::Gi, Nuc::G]
        | [Nuc::Ai, Nuc::Gi, Nuc::Gi] => Some(10),
        [Nuc::A, Nuc::G, Nuc::T]
        | [Nuc::A, Nuc::G, Nuc::Ti]
        | [Nuc::A, Nuc::Gi, Nuc::T]
        | [Nuc::A, Nuc::Gi, Nuc::Ti]
        | [Nuc::Ai, Nuc::G, Nuc::T]
        | [Nuc::Ai, Nuc::G, Nuc::Ti]
        | [Nuc::Ai, Nuc::Gi, Nuc::T]
        | [Nuc::Ai, Nuc::Gi, Nuc::Ti] => Some(11),
        [Nuc::A, Nuc::T, Nuc::A]
        | [Nuc::A, Nuc::T, Nuc::Ai]
        | [Nuc::A, Nuc::Ti, Nuc::A]
        | [Nuc::A, Nuc::Ti, Nuc::Ai]
        | [Nuc::Ai, Nuc::T, Nuc::A]
        | [Nuc::Ai, Nuc::T, Nuc::Ai]
        | [Nuc::Ai, Nuc::Ti, Nuc::A]
        | [Nuc::Ai, Nuc::Ti, Nuc::Ai] => Some(12),
        [Nuc::A, Nuc::T, Nuc::C]
        | [Nuc::A, Nuc::T, Nuc::Ci]
        | [Nuc::A, Nuc::Ti, Nuc::C]
        | [Nuc::A, Nuc::Ti, Nuc::Ci]
        | [Nuc::Ai, Nuc::T, Nuc::C]
        | [Nuc::Ai, Nuc::T, Nuc::Ci]
        | [Nuc::Ai, Nuc::Ti, Nuc::C]
        | [Nuc::Ai, Nuc::Ti, Nuc::Ci] => Some(13),
        [Nuc::A, Nuc::T, Nuc::G]
        | [Nuc::A, Nuc::T, Nuc::Gi]
        | [Nuc::A, Nuc::Ti, Nuc::G]
        | [Nuc::A, Nuc::Ti, Nuc::Gi]
        | [Nuc::Ai, Nuc::T, Nuc::G]
        | [Nuc::Ai, Nuc::T, Nuc::Gi]
        | [Nuc::Ai, Nuc::Ti, Nuc::G]
        | [Nuc::Ai, Nuc::Ti, Nuc::Gi] => Some(14),
        [Nuc::A, Nuc::T, Nuc::T]
        | [Nuc::A, Nuc::T, Nuc::Ti]
        | [Nuc::A, Nuc::Ti, Nuc::T]
        | [Nuc::A, Nuc::Ti, Nuc::Ti]
        | [Nuc::Ai, Nuc::T, Nuc::T]
        | [Nuc::Ai, Nuc::T, Nuc::Ti]
        | [Nuc::Ai, Nuc::Ti, Nuc::T]
        | [Nuc::Ai, Nuc::Ti, Nuc::Ti] => Some(15),
        [Nuc::C, Nuc::A, Nuc::A]
        | [Nuc::C, Nuc::A, Nuc::Ai]
        | [Nuc::C, Nuc::Ai, Nuc::A]
        | [Nuc::C, Nuc::Ai, Nuc::Ai]
        | [Nuc::Ci, Nuc::A, Nuc::A]
        | [Nuc::Ci, Nuc::A, Nuc::Ai]
        | [Nuc::Ci, Nuc::Ai, Nuc::A]
        | [Nuc::Ci, Nuc::Ai, Nuc::Ai] => Some(16),
        [Nuc::C, Nuc::A, Nuc::C]
        | [Nuc::C, Nuc::A, Nuc::Ci]
        | [Nuc::C, Nuc::Ai, Nuc::C]
        | [Nuc::C, Nuc::Ai, Nuc::Ci]
        | [Nuc::Ci, Nuc::A, Nuc::C]
        | [Nuc::Ci, Nuc::A, Nuc::Ci]
        | [Nuc::Ci, Nuc::Ai, Nuc::C]
        | [Nuc::Ci, Nuc::Ai, Nuc::Ci] => Some(17),
        [Nuc::C, Nuc::A, Nuc::G]
        | [Nuc::C, Nuc::A, Nuc::Gi]
        | [Nuc::C, Nuc::Ai, Nuc::G]
        | [Nuc::C, Nuc::Ai, Nuc::Gi]
        | [Nuc::Ci, Nuc::A, Nuc::G]
        | [Nuc::Ci, Nuc::A, Nuc::Gi]
        | [Nuc::Ci, Nuc::Ai, Nuc::G]
        | [Nuc::Ci, Nuc::Ai, Nuc::Gi] => Some(18),
        [Nuc::C, Nuc::A, Nuc::T]
        | [Nuc::C, Nuc::A, Nuc::Ti]
        | [Nuc::C, Nuc::Ai, Nuc::T]
        | [Nuc::C, Nuc::Ai, Nuc::Ti]
        | [Nuc::Ci, Nuc::A, Nuc::T]
        | [Nuc::Ci, Nuc::A, Nuc::Ti]
        | [Nuc::Ci, Nuc::Ai, Nuc::T]
        | [Nuc::Ci, Nuc::Ai, Nuc::Ti] => Some(19),
        [Nuc::C, Nuc::C, Nuc::A]
        | [Nuc::C, Nuc::C, Nuc::Ai]
        | [Nuc::C, Nuc::Ci, Nuc::A]
        | [Nuc::C, Nuc::Ci, Nuc::Ai]
        | [Nuc::Ci, Nuc::C, Nuc::A]
        | [Nuc::Ci, Nuc::C, Nuc::Ai]
        | [Nuc::Ci, Nuc::Ci, Nuc::A]
        | [Nuc::Ci, Nuc::Ci, Nuc::Ai] => Some(20),
        [Nuc::C, Nuc::C, Nuc::C]
        | [Nuc::C, Nuc::C, Nuc::Ci]
        | [Nuc::C, Nuc::Ci, Nuc::C]
        | [Nuc::C, Nuc::Ci, Nuc::Ci]
        | [Nuc::Ci, Nuc::C, Nuc::C]
        | [Nuc::Ci, Nuc::C, Nuc::Ci]
        | [Nuc::Ci, Nuc::Ci, Nuc::C]
        | [Nuc::Ci, Nuc::Ci, Nuc::Ci] => Some(21),
        [Nuc::C, Nuc::C, Nuc::G]
        | [Nuc::C, Nuc::C, Nuc::Gi]
        | [Nuc::C, Nuc::Ci, Nuc::G]
        | [Nuc::C, Nuc::Ci, Nuc::Gi]
        | [Nuc::Ci, Nuc::C, Nuc::G]
        | [Nuc::Ci, Nuc::C, Nuc::Gi]
        | [Nuc::Ci, Nuc::Ci, Nuc::G]
        | [Nuc::Ci, Nuc::Ci, Nuc::Gi] => Some(22),
        [Nuc::C, Nuc::C, Nuc::T]
        | [Nuc::C, Nuc::C, Nuc::Ti]
        | [Nuc::C, Nuc::Ci, Nuc::T]
        | [Nuc::C, Nuc::Ci, Nuc::Ti]
        | [Nuc::Ci, Nuc::C, Nuc::T]
        | [Nuc::Ci, Nuc::C, Nuc::Ti]
        | [Nuc::Ci, Nuc::Ci, Nuc::T]
        | [Nuc::Ci, Nuc::Ci, Nuc::Ti] => Some(23),
        [Nuc::C, Nuc::G, Nuc::A]
        | [Nuc::C, Nuc::G, Nuc::Ai]
        | [Nuc::C, Nuc::Gi, Nuc::A]
        | [Nuc::C, Nuc::Gi, Nuc::Ai]
        | [Nuc::Ci, Nuc::G, Nuc::A]
        | [Nuc::Ci, Nuc::G, Nuc::Ai]
        | [Nuc::Ci, Nuc::Gi, Nuc::A]
        | [Nuc::Ci, Nuc::Gi, Nuc::Ai] => Some(24),
        [Nuc::C, Nuc::G, Nuc::C]
        | [Nuc::C, Nuc::G, Nuc::Ci]
        | [Nuc::C, Nuc::Gi, Nuc::C]
        | [Nuc::C, Nuc::Gi, Nuc::Ci]
        | [Nuc::Ci, Nuc::G, Nuc::C]
        | [Nuc::Ci, Nuc::G, Nuc::Ci]
        | [Nuc::Ci, Nuc::Gi, Nuc::C]
        | [Nuc::Ci, Nuc::Gi, Nuc::Ci] => Some(25),
        [Nuc::C, Nuc::G, Nuc::G]
        | [Nuc::C, Nuc::G, Nuc::Gi]
        | [Nuc::C, Nuc::Gi, Nuc::G]
        | [Nuc::C, Nuc::Gi, Nuc::Gi]
        | [Nuc::Ci, Nuc::G, Nuc::G]
        | [Nuc::Ci, Nuc::G, Nuc::Gi]
        | [Nuc::Ci, Nuc::Gi, Nuc::G]
        | [Nuc::Ci, Nuc::Gi, Nuc::Gi] => Some(26),
        [Nuc::C, Nuc::G, Nuc::T]
        | [Nuc::C, Nuc::G, Nuc::Ti]
        | [Nuc::C, Nuc::Gi, Nuc::T]
        | [Nuc::C, Nuc::Gi, Nuc::Ti]
        | [Nuc::Ci, Nuc::G, Nuc::T]
        | [Nuc::Ci, Nuc::G, Nuc::Ti]
        | [Nuc::Ci, Nuc::Gi, Nuc::T]
        | [Nuc::Ci, Nuc::Gi, Nuc::Ti] => Some(27),
        [Nuc::C, Nuc::T, Nuc::A]
        | [Nuc::C, Nuc::T, Nuc::Ai]
        | [Nuc::C, Nuc::Ti, Nuc::A]
        | [Nuc::C, Nuc::Ti, Nuc::Ai]
        | [Nuc::Ci, Nuc::T, Nuc::A]
        | [Nuc::Ci, Nuc::T, Nuc::Ai]
        | [Nuc::Ci, Nuc::Ti, Nuc::A]
        | [Nuc::Ci, Nuc::Ti, Nuc::Ai] => Some(28),
        [Nuc::C, Nuc::T, Nuc::C]
        | [Nuc::C, Nuc::T, Nuc::Ci]
        | [Nuc::C, Nuc::Ti, Nuc::C]
        | [Nuc::C, Nuc::Ti, Nuc::Ci]
        | [Nuc::Ci, Nuc::T, Nuc::C]
        | [Nuc::Ci, Nuc::T, Nuc::Ci]
        | [Nuc::Ci, Nuc::Ti, Nuc::C]
        | [Nuc::Ci, Nuc::Ti, Nuc::Ci] => Some(29),
        [Nuc::C, Nuc::T, Nuc::G]
        | [Nuc::C, Nuc::T, Nuc::Gi]
        | [Nuc::C, Nuc::Ti, Nuc::G]
        | [Nuc::C, Nuc::Ti, Nuc::Gi]
        | [Nuc::Ci, Nuc::T, Nuc::G]
        | [Nuc::Ci, Nuc::T, Nuc::Gi]
        | [Nuc::Ci, Nuc::Ti, Nuc::G]
        | [Nuc::Ci, Nuc::Ti, Nuc::Gi] => Some(30),
        [Nuc::C, Nuc::T, Nuc::T]
        | [Nuc::C, Nuc::T, Nuc::Ti]
        | [Nuc::C, Nuc::Ti, Nuc::T]
        | [Nuc::C, Nuc::Ti, Nuc::Ti]
        | [Nuc::Ci, Nuc::T, Nuc::T]
        | [Nuc::Ci, Nuc::T, Nuc::Ti]
        | [Nuc::Ci, Nuc::Ti, Nuc::T]
        | [Nuc::Ci, Nuc::Ti, Nuc::Ti] => Some(31),
        [Nuc::G, Nuc::A, Nuc::A]
        | [Nuc::G, Nuc::A, Nuc::Ai]
        | [Nuc::G, Nuc::Ai, Nuc::A]
        | [Nuc::G, Nuc::Ai, Nuc::Ai]
        | [Nuc::Gi, Nuc::A, Nuc::A]
        | [Nuc::Gi, Nuc::A, Nuc::Ai]
        | [Nuc::Gi, Nuc::Ai, Nuc::A]
        | [Nuc::Gi, Nuc::Ai, Nuc::Ai] => Some(32),
        [Nuc::G, Nuc::A, Nuc::C]
        | [Nuc::G, Nuc::A, Nuc::Ci]
        | [Nuc::G, Nuc::Ai, Nuc::C]
        | [Nuc::G, Nuc::Ai, Nuc::Ci]
        | [Nuc::Gi, Nuc::A, Nuc::C]
        | [Nuc::Gi, Nuc::A, Nuc::Ci]
        | [Nuc::Gi, Nuc::Ai, Nuc::C]
        | [Nuc::Gi, Nuc::Ai, Nuc::Ci] => Some(33),
        [Nuc::G, Nuc::A, Nuc::G]
        | [Nuc::G, Nuc::A, Nuc::Gi]
        | [Nuc::G, Nuc::Ai, Nuc::G]
        | [Nuc::G, Nuc::Ai, Nuc::Gi]
        | [Nuc::Gi, Nuc::A, Nuc::G]
        | [Nuc::Gi, Nuc::A, Nuc::Gi]
        | [Nuc::Gi, Nuc::Ai, Nuc::G]
        | [Nuc::Gi, Nuc::Ai, Nuc::Gi] => Some(34),
        [Nuc::G, Nuc::A, Nuc::T]
        | [Nuc::G, Nuc::A, Nuc::Ti]
        | [Nuc::G, Nuc::Ai, Nuc::T]
        | [Nuc::G, Nuc::Ai, Nuc::Ti]
        | [Nuc::Gi, Nuc::A, Nuc::T]
        | [Nuc::Gi, Nuc::A, Nuc::Ti]
        | [Nuc::Gi, Nuc::Ai, Nuc::T]
        | [Nuc::Gi, Nuc::Ai, Nuc::Ti] => Some(35),
        [Nuc::G, Nuc::C, Nuc::A]
        | [Nuc::G, Nuc::C, Nuc::Ai]
        | [Nuc::G, Nuc::Ci, Nuc::A]
        | [Nuc::G, Nuc::Ci, Nuc::Ai]
        | [Nuc::Gi, Nuc::C, Nuc::A]
        | [Nuc::Gi, Nuc::C, Nuc::Ai]
        | [Nuc::Gi, Nuc::Ci, Nuc::A]
        | [Nuc::Gi, Nuc::Ci, Nuc::Ai] => Some(36),
        [Nuc::G, Nuc::C, Nuc::C]
        | [Nuc::G, Nuc::C, Nuc::Ci]
        | [Nuc::G, Nuc::Ci, Nuc::C]
        | [Nuc::G, Nuc::Ci, Nuc::Ci]
        | [Nuc::Gi, Nuc::C, Nuc::C]
        | [Nuc::Gi, Nuc::C, Nuc::Ci]
        | [Nuc::Gi, Nuc::Ci, Nuc::C]
        | [Nuc::Gi, Nuc::Ci, Nuc::Ci] => Some(37),
        [Nuc::G, Nuc::C, Nuc::G]
        | [Nuc::G, Nuc::C, Nuc::Gi]
        | [Nuc::G, Nuc::Ci, Nuc::G]
        | [Nuc::G, Nuc::Ci, Nuc::Gi]
        | [Nuc::Gi, Nuc::C, Nuc::G]
        | [Nuc::Gi, Nuc::C, Nuc::Gi]
        | [Nuc::Gi, Nuc::Ci, Nuc::G]
        | [Nuc::Gi, Nuc::Ci, Nuc::Gi] => Some(38),
        [Nuc::G, Nuc::C, Nuc::T]
        | [Nuc::G, Nuc::C, Nuc::Ti]
        | [Nuc::G, Nuc::Ci, Nuc::T]
        | [Nuc::G, Nuc::Ci, Nuc::Ti]
        | [Nuc::Gi, Nuc::C, Nuc::T]
        | [Nuc::Gi, Nuc::C, Nuc::Ti]
        | [Nuc::Gi, Nuc::Ci, Nuc::T]
        | [Nuc::Gi, Nuc::Ci, Nuc::Ti] => Some(39),
        [Nuc::G, Nuc::G, Nuc::A]
        | [Nuc::G, Nuc::G, Nuc::Ai]
        | [Nuc::G, Nuc::Gi, Nuc::A]
        | [Nuc::G, Nuc::Gi, Nuc::Ai]
        | [Nuc::Gi, Nuc::G, Nuc::A]
        | [Nuc::Gi, Nuc::G, Nuc::Ai]
        | [Nuc::Gi, Nuc::Gi, Nuc::A]
        | [Nuc::Gi, Nuc::Gi, Nuc::Ai] => Some(40),
        [Nuc::G, Nuc::G, Nuc::C]
        | [Nuc::G, Nuc::G, Nuc::Ci]
        | [Nuc::G, Nuc::Gi, Nuc::C]
        | [Nuc::G, Nuc::Gi, Nuc::Ci]
        | [Nuc::Gi, Nuc::G, Nuc::C]
        | [Nuc::Gi, Nuc::G, Nuc::Ci]
        | [Nuc::Gi, Nuc::Gi, Nuc::C]
        | [Nuc::Gi, Nuc::Gi, Nuc::Ci] => Some(41),
        [Nuc::G, Nuc::G, Nuc::G]
        | [Nuc::G, Nuc::G, Nuc::Gi]
        | [Nuc::G, Nuc::Gi, Nuc::G]
        | [Nuc::G, Nuc::Gi, Nuc::Gi]
        | [Nuc::Gi, Nuc::G, Nuc::G]
        | [Nuc::Gi, Nuc::G, Nuc::Gi]
        | [Nuc::Gi, Nuc::Gi, Nuc::G]
        | [Nuc::Gi, Nuc::Gi, Nuc::Gi] => Some(42),
        [Nuc::G, Nuc::G, Nuc::T]
        | [Nuc::G, Nuc::G, Nuc::Ti]
        | [Nuc::G, Nuc::Gi, Nuc::T]
        | [Nuc::G, Nuc::Gi, Nuc::Ti]
        | [Nuc::Gi, Nuc::G, Nuc::T]
        | [Nuc::Gi, Nuc::G, Nuc::Ti]
        | [Nuc::Gi, Nuc::Gi, Nuc::T]
        | [Nuc::Gi, Nuc::Gi, Nuc::Ti] => Some(43),
        [Nuc::G, Nuc::T, Nuc::A]
        | [Nuc::G, Nuc::T, Nuc::Ai]
        | [Nuc::G, Nuc::Ti, Nuc::A]
        | [Nuc::G, Nuc::Ti, Nuc::Ai]
        | [Nuc::Gi, Nuc::T, Nuc::A]
        | [Nuc::Gi, Nuc::T, Nuc::Ai]
        | [Nuc::Gi, Nuc::Ti, Nuc::A]
        | [Nuc::Gi, Nuc::Ti, Nuc::Ai] => Some(44),
        [Nuc::G, Nuc::T, Nuc::C]
        | [Nuc::G, Nuc::T, Nuc::Ci]
        | [Nuc::G, Nuc::Ti, Nuc::C]
        | [Nuc::G, Nuc::Ti, Nuc::Ci]
        | [Nuc::Gi, Nuc::T, Nuc::C]
        | [Nuc::Gi, Nuc::T, Nuc::Ci]
        | [Nuc::Gi, Nuc::Ti, Nuc::C]
        | [Nuc::Gi, Nuc::Ti, Nuc::Ci] => Some(45),
        [Nuc::G, Nuc::T, Nuc::G]
        | [Nuc::G, Nuc::T, Nuc::Gi]
        | [Nuc::G, Nuc::Ti, Nuc::G]
        | [Nuc::G, Nuc::Ti, Nuc::Gi]
        | [Nuc::Gi, Nuc::T, Nuc::G]
        | [Nuc::Gi, Nuc::T, Nuc::Gi]
        | [Nuc::Gi, Nuc::Ti, Nuc::G]
        | [Nuc::Gi, Nuc::Ti, Nuc::Gi] => Some(46),
        [Nuc::G, Nuc::T, Nuc::T]
        | [Nuc::G, Nuc::T, Nuc::Ti]
        | [Nuc::G, Nuc::Ti, Nuc::T]
        | [Nuc::G, Nuc::Ti, Nuc::Ti]
        | [Nuc::Gi, Nuc::T, Nuc::T]
        | [Nuc::Gi, Nuc::T, Nuc::Ti]
        | [Nuc::Gi, Nuc::Ti, Nuc::T]
        | [Nuc::Gi, Nuc::Ti, Nuc::Ti] => Some(47),
        [Nuc::T, Nuc::A, Nuc::A]
        | [Nuc::T, Nuc::A, Nuc::Ai]
        | [Nuc::T, Nuc::Ai, Nuc::A]
        | [Nuc::T, Nuc::Ai, Nuc::Ai]
        | [Nuc::Ti, Nuc::A, Nuc::A]
        | [Nuc::Ti, Nuc::A, Nuc::Ai]
        | [Nuc::Ti, Nuc::Ai, Nuc::A]
        | [Nuc::Ti, Nuc::Ai, Nuc::Ai] => Some(48),
        [Nuc::T, Nuc::A, Nuc::C]
        | [Nuc::T, Nuc::A, Nuc::Ci]
        | [Nuc::T, Nuc::Ai, Nuc::C]
        | [Nuc::T, Nuc::Ai, Nuc::Ci]
        | [Nuc::Ti, Nuc::A, Nuc::C]
        | [Nuc::Ti, Nuc::A, Nuc::Ci]
        | [Nuc::Ti, Nuc::Ai, Nuc::C]
        | [Nuc::Ti, Nuc::Ai, Nuc::Ci] => Some(49),
        [Nuc::T, Nuc::A, Nuc::G]
        | [Nuc::T, Nuc::A, Nuc::Gi]
        | [Nuc::T, Nuc::Ai, Nuc::G]
        | [Nuc::T, Nuc::Ai, Nuc::Gi]
        | [Nuc::Ti, Nuc::A, Nuc::G]
        | [Nuc::Ti, Nuc::A, Nuc::Gi]
        | [Nuc::Ti, Nuc::Ai, Nuc::G]
        | [Nuc::Ti, Nuc::Ai, Nuc::Gi] => Some(50),
        [Nuc::T, Nuc::A, Nuc::T]
        | [Nuc::T, Nuc::A, Nuc::Ti]
        | [Nuc::T, Nuc::Ai, Nuc::T]
        | [Nuc::T, Nuc::Ai, Nuc::Ti]
        | [Nuc::Ti, Nuc::A, Nuc::T]
        | [Nuc::Ti, Nuc::A, Nuc::Ti]
        | [Nuc::Ti, Nuc::Ai, Nuc::T]
        | [Nuc::Ti, Nuc::Ai, Nuc::Ti] => Some(51),
        [Nuc::T, Nuc::C, Nuc::A]
        | [Nuc::T, Nuc::C, Nuc::Ai]
        | [Nuc::T, Nuc::Ci, Nuc::A]
        | [Nuc::T, Nuc::Ci, Nuc::Ai]
        | [Nuc::Ti, Nuc::C, Nuc::A]
        | [Nuc::Ti, Nuc::C, Nuc::Ai]
        | [Nuc::Ti, Nuc::Ci, Nuc::A]
        | [Nuc::Ti, Nuc::Ci, Nuc::Ai] => Some(52),
        [Nuc::T, Nuc::C, Nuc::C]
        | [Nuc::T, Nuc::C, Nuc::Ci]
        | [Nuc::T, Nuc::Ci, Nuc::C]
        | [Nuc::T, Nuc::Ci, Nuc::Ci]
        | [Nuc::Ti, Nuc::C, Nuc::C]
        | [Nuc::Ti, Nuc::C, Nuc::Ci]
        | [Nuc::Ti, Nuc::Ci, Nuc::C]
        | [Nuc::Ti, Nuc::Ci, Nuc::Ci] => Some(53),
        [Nuc::T, Nuc::C, Nuc::G]
        | [Nuc::T, Nuc::C, Nuc::Gi]
        | [Nuc::T, Nuc::Ci, Nuc::G]
        | [Nuc::T, Nuc::Ci, Nuc::Gi]
        | [Nuc::Ti, Nuc::C, Nuc::G]
        | [Nuc::Ti, Nuc::C, Nuc::Gi]
        | [Nuc::Ti, Nuc::Ci, Nuc::G]
        | [Nuc::Ti, Nuc::Ci, Nuc::Gi] => Some(54),
        [Nuc::T, Nuc::C, Nuc::T]
        | [Nuc::T, Nuc::C, Nuc::Ti]
        | [Nuc::T, Nuc::Ci, Nuc::T]
        | [Nuc::T, Nuc::Ci, Nuc::Ti]
        | [Nuc::Ti, Nuc::C, Nuc::T]
        | [Nuc::Ti, Nuc::C, Nuc::Ti]
        | [Nuc::Ti, Nuc::Ci, Nuc::T]
        | [Nuc::Ti, Nuc::Ci, Nuc::Ti] => Some(55),
        [Nuc::T, Nuc::G, Nuc::A]
        | [Nuc::T, Nuc::G, Nuc::Ai]
        | [Nuc::T, Nuc::Gi, Nuc::A]
        | [Nuc::T, Nuc::Gi, Nuc::Ai]
        | [Nuc::Ti, Nuc::G, Nuc::A]
        | [Nuc::Ti, Nuc::G, Nuc::Ai]
        | [Nuc::Ti, Nuc::Gi, Nuc::A]
        | [Nuc::Ti, Nuc::Gi, Nuc::Ai] => Some(56),
        [Nuc::T, Nuc::G, Nuc::C]
        | [Nuc::T, Nuc::G, Nuc::Ci]
        | [Nuc::T, Nuc::Gi, Nuc::C]
        | [Nuc::T, Nuc::Gi, Nuc::Ci]
        | [Nuc::Ti, Nuc::G, Nuc::C]
        | [Nuc::Ti, Nuc::G, Nuc::Ci]
        | [Nuc::Ti, Nuc::Gi, Nuc::C]
        | [Nuc::Ti, Nuc::Gi, Nuc::Ci] => Some(57),
        [Nuc::T, Nuc::G, Nuc::G]
        | [Nuc::T, Nuc::G, Nuc::Gi]
        | [Nuc::T, Nuc::Gi, Nuc::G]
        | [Nuc::T, Nuc::Gi, Nuc::Gi]
        | [Nuc::Ti, Nuc::G, Nuc::G]
        | [Nuc::Ti, Nuc::G, Nuc::Gi]
        | [Nuc::Ti, Nuc::Gi, Nuc::G]
        | [Nuc::Ti, Nuc::Gi, Nuc::Gi] => Some(58),
        [Nuc::T, Nuc::G, Nuc::T]
        | [Nuc::T, Nuc::G, Nuc::Ti]
        | [Nuc::T, Nuc::Gi, Nuc::T]
        | [Nuc::T, Nuc::Gi, Nuc::Ti]
        | [Nuc::Ti, Nuc::G, Nuc::T]
        | [Nuc::Ti, Nuc::G, Nuc::Ti]
        | [Nuc::Ti, Nuc::Gi, Nuc::T]
        | [Nuc::Ti, Nuc::Gi, Nuc::Ti] => Some(59),
        [Nuc::T, Nuc::T, Nuc::A]
        | [Nuc::T, Nuc::T, Nuc::Ai]
        | [Nuc::T, Nuc::Ti, Nuc::A]
        | [Nuc::T, Nuc::Ti, Nuc::Ai]
        | [Nuc::Ti, Nuc::T, Nuc::A]
        | [Nuc::Ti, Nuc::T, Nuc::Ai]
        | [Nuc::Ti, Nuc::Ti, Nuc::A]
        | [Nuc::Ti, Nuc::Ti, Nuc::Ai] => Some(60),
        [Nuc::T, Nuc::T, Nuc::C]
        | [Nuc::T, Nuc::T, Nuc::Ci]
        | [Nuc::T, Nuc::Ti, Nuc::C]
        | [Nuc::T, Nuc::Ti, Nuc::Ci]
        | [Nuc::Ti, Nuc::T, Nuc::C]
        | [Nuc::Ti, Nuc::T, Nuc::Ci]
        | [Nuc::Ti, Nuc::Ti, Nuc::C]
        | [Nuc::Ti, Nuc::Ti, Nuc::Ci] => Some(61),
        [Nuc::T, Nuc::T, Nuc::G]
        | [Nuc::T, Nuc::T, Nuc::Gi]
        | [Nuc::T, Nuc::Ti, Nuc::G]
        | [Nuc::T, Nuc::Ti, Nuc::Gi]
        | [Nuc::Ti, Nuc::T, Nuc::G]
        | [Nuc::Ti, Nuc::T, Nuc::Gi]
        | [Nuc::Ti, Nuc::Ti, Nuc::G]
        | [Nuc::Ti, Nuc::Ti, Nuc::Gi] => Some(62),
        [Nuc::T, Nuc::T, Nuc::T]
        | [Nuc::T, Nuc::T, Nuc::Ti]
        | [Nuc::T, Nuc::Ti, Nuc::T]
        | [Nuc::T, Nuc::Ti, Nuc::Ti]
        | [Nuc::Ti, Nuc::T, Nuc::T]
        | [Nuc::Ti, Nuc::T, Nuc::Ti]
        | [Nuc::Ti, Nuc::Ti, Nuc::T]
        | [Nuc::Ti, Nuc::Ti, Nuc::Ti] => Some(63),
        _ => None,
    }
}

pub fn count_cg_content(seq: &[Nuc]) -> usize {
    let mut count = 0;
    for &l in seq.iter() {
        if l == Nuc::C || l == Nuc::G {
            count += 1;
        }
    }
    min(
        CG_MAX - CG_MIN - 1,
        max(CG_MIN, count * 100 / seq.len()) - CG_MIN,
    )
}
