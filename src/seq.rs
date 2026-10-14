//! Encoded nucleotide sequence and the low-level sequence statistics shared by
//! both gene-calling engines.
//!
//! Bases are packed into small integers so that k-mer indices, GC tests and
//! translation are all cheap and mutually consistent:
//!
//! ```text
//!   A = 0b00   C = 0b01   G = 0b10   T = 0b11
//! ```
//!
//! (Prodigal packs the same four bases as A=00, G=01, C=10, T=11. The two
//! layouts index the 4096 hexamers in a different *order*, but every statistic
//! is trained and scored inside one layout, so the models and gene calls are
//! identical — only a raw dump of the hexamer table would need re-indexing.)
//! `is_gc` is simply "the code is C or G".
//!
//! ## Ambiguous bases (N, IUPAC codes)
//!
//! Ambiguous / non-ACGT characters must never masquerade as a real base in a
//! *decision*. Each one is stored as a placeholder code (C — exactly what
//! Prodigal does internally, so the k-mer / GC statistics stay Prodigal-
//! identical) **and** flagged in the [`Seq::n`] mask. Every codon-level
//! decision consults that mask:
//!
//! * a codon containing an ambiguous base is never a start or a stop codon,
//!   under any genetic code ([`Seq::codon_has_amb_fwd`] / [`Seq::codon_has_amb_rev`]);
//! * it translates to `X`, and the nucleotide is written back out as `N`.
//!
//! Runs of at least [`GAP_MIN_RUN`] ambiguous bases are treated as **assembly
//! gaps**: callers split contigs at them ([`gap_free_segments`]) so that no gene
//! can span, start in, or end in a gap — genes that run into one are reported
//! as partial at that end, exactly like genes that run off a contig edge.

/// Numeric code for adenine.
pub const A: u8 = 0;
/// Numeric code for cytosine.
pub const C: u8 = 1;
/// Numeric code for guanine.
pub const G: u8 = 2;
/// Numeric code for thymine.
pub const T: u8 = 3;

/// Placeholder code stored for an ambiguous base (always paired with a `true`
/// entry in the ambiguity mask). C is Prodigal's choice: it keeps the hexamer
/// and GC-frame statistics identical to Prodigal's on sequences containing N.
pub const AMBIGUOUS_CODE: u8 = C;

/// Minimum run of consecutive ambiguous bases (N or other IUPAC codes) that is
/// treated as an assembly gap. Genes never span, start in, or end in a gap;
/// shorter runs are kept in-frame (translated as `X`, never a start/stop).
pub const GAP_MIN_RUN: usize = 10;

/// Minimum run of consecutive ambiguous bases masked by Prodigal 2.x's `-m`
/// option (Prodigal's `MASK_SIZE`; pyrodigal's `min_mask` default).
pub const MASK_MIN_RUN: usize = 50;

/// Map an input byte to a 2-bit base code, returning `(code, is_ambiguous)`.
///
/// `U`/`u` is treated as `T`. Every other non-ACGT character (N and the IUPAC
/// ambiguity codes, gaps, digits …) is flagged ambiguous and stored as the
/// [`AMBIGUOUS_CODE`] placeholder; callers must check the flag before using
/// the code in any codon-level decision.
#[inline]
pub fn encode_base(b: u8) -> (u8, bool) {
    match b {
        b'A' | b'a' => (A, false),
        b'C' | b'c' => (C, false),
        b'G' | b'g' => (G, false),
        b'T' | b't' | b'U' | b'u' => (T, false),
        _ => (AMBIGUOUS_CODE, true),
    }
}

/// `true` if the raw input byte is not an unambiguous A/C/G/T/U.
#[inline]
pub fn is_ambiguous_byte(b: u8) -> bool {
    encode_base(b).1
}

/// Decode a 2-bit base code back to an upper-case ASCII nucleotide.
#[inline]
pub fn decode_base(code: u8) -> u8 {
    match code & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// An encoded contig: forward strand, its reverse complement, and the
/// ambiguity bitmap. Both strands are pre-computed because the engines scan
/// each strand independently many times.
#[derive(Clone)]
pub struct Seq {
    /// Forward-strand base codes (`0..=3`), 5' -> 3'.
    pub fwd: Vec<u8>,
    /// Reverse-complement base codes (`0..=3`), 5' -> 3' of the minus strand.
    pub rev: Vec<u8>,
    /// `true` where the original input base was ambiguous / non-ACGT (forward
    /// coordinates). Reverse-strand position `i` maps to `n[len - 1 - i]`.
    pub n: Vec<bool>,
    /// GC fraction over unambiguous bases (ambiguous bases are excluded from
    /// both numerator and denominator).
    pub gc: f64,
}

impl Seq {
    /// Build an encoded sequence from raw ASCII nucleotides.
    pub fn from_ascii(bytes: &[u8]) -> Seq {
        let len = bytes.len();
        let mut fwd = Vec::with_capacity(len);
        let mut n = Vec::with_capacity(len);
        let mut gc_count = 0u64;
        let mut acgt = 0u64;
        for &b in bytes {
            let (code, amb) = encode_base(b);
            if !amb {
                acgt += 1;
                if code == C || code == G {
                    gc_count += 1;
                }
            }
            fwd.push(code);
            n.push(amb);
        }
        let mut rev = revcomp(&fwd);
        // Prodigal's `rcom_seq` keeps an ambiguous base as the placeholder on
        // the reverse strand too (complementing it would turn C into G and
        // change every reverse-strand statistic on input containing N).
        let len = n.len();
        for (i, &amb) in n.iter().enumerate() {
            if amb {
                rev[len - 1 - i] = AMBIGUOUS_CODE;
            }
        }
        let gc = if acgt > 0 {
            gc_count as f64 / acgt as f64
        } else {
            0.5
        };
        Seq { fwd, rev, n, gc }
    }

    /// Length of the contig in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.fwd.len()
    }

    /// Whether the contig is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.fwd.is_empty()
    }

    /// `true` if the base at forward position `i` is ambiguous (out of range
    /// counts as ambiguous).
    #[inline]
    pub fn amb_fwd(&self, i: usize) -> bool {
        self.n.get(i).copied().unwrap_or(true)
    }

    /// `true` if the base at reverse-strand position `i` (5'->3' of the minus
    /// strand, i.e. an index into [`Seq::rev`]) is ambiguous.
    #[inline]
    pub fn amb_rev(&self, i: usize) -> bool {
        let len = self.n.len();
        if i >= len {
            return true;
        }
        self.n[len - 1 - i]
    }

    /// `true` if any base of the forward codon starting at `i` is ambiguous
    /// (or the codon runs off the sequence).
    #[inline]
    pub fn codon_has_amb_fwd(&self, i: usize) -> bool {
        self.amb_fwd(i) || self.amb_fwd(i + 1) || self.amb_fwd(i + 2)
    }

    /// `true` if any base of the reverse-strand codon starting at `i` (index
    /// into [`Seq::rev`]) is ambiguous (or the codon runs off the sequence).
    #[inline]
    pub fn codon_has_amb_rev(&self, i: usize) -> bool {
        self.amb_rev(i) || self.amb_rev(i + 1) || self.amb_rev(i + 2)
    }

    /// The ambiguity mask in reverse-strand coordinates (aligned with
    /// [`Seq::rev`]).
    pub fn n_rev(&self) -> Vec<bool> {
        self.n.iter().rev().copied().collect()
    }

    /// `true` if the sequence contains at least one ambiguous base.
    pub fn has_ambiguity(&self) -> bool {
        self.n.iter().any(|&b| b)
    }

    /// GC fraction exactly as Prodigal computes it for training:
    /// unambiguous G+C over the *total* length (ambiguous bases count in the
    /// denominator only).
    pub fn gc_prodigal(&self) -> f64 {
        if self.fwd.is_empty() {
            return 0.0;
        }
        let gc = self
            .fwd
            .iter()
            .zip(self.n.iter())
            .filter(|(&c, &amb)| !amb && (c == C || c == G))
            .count();
        gc as f64 / self.fwd.len() as f64
    }

    /// Reconstruct the forward strand as uppercase ASCII (`A/C/G/T`), with
    /// ambiguous positions rendered as `N`. Used to feed sequence to vendored
    /// engines that parse raw ASCII (e.g. the FragGeneScan HMM).
    pub fn to_ascii_upper(&self) -> Vec<u8> {
        self.fwd
            .iter()
            .zip(self.n.iter())
            .map(|(&c, &amb)| if amb { b'N' } else { decode_base(c) })
            .collect()
    }

    /// Copy of the sub-sequence `[start, start + len)` (forward coordinates)
    /// as an independent [`Seq`] (ambiguity preserved).
    pub fn subseq(&self, start: usize, len: usize) -> Seq {
        let end = (start + len).min(self.len());
        let start = start.min(end);
        let ascii: Vec<u8> = self.fwd[start..end]
            .iter()
            .zip(self.n[start..end].iter())
            .map(|(&c, &amb)| if amb { b'N' } else { decode_base(c) })
            .collect();
        Seq::from_ascii(&ascii)
    }
}

/// Runs of at least `min_run` consecutive ambiguous bases in an encoded
/// sequence, as half-open `[begin, end)` 0-based intervals in increasing order
/// (a trailing run is included). These are Prodigal's masked regions (`-m`,
/// with `min_run` = [`MASK_MIN_RUN`]) and DYNAMITE's assembly gaps (with
/// [`GAP_MIN_RUN`]).
pub fn ambiguous_runs(seq: &Seq, min_run: usize) -> Vec<(usize, usize)> {
    let min_run = min_run.max(1);
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < seq.n.len() {
        if seq.n[i] {
            let start = i;
            while i < seq.n.len() && seq.n[i] {
                i += 1;
            }
            if i - start >= min_run {
                runs.push((start, i));
            }
        } else {
            i += 1;
        }
    }
    runs
}

/// Split a sequence into maximal stretches that contain no assembly gap
/// (a run of at least `min_run` ambiguous bases). `is_amb(i)` reports whether
/// position `i` (of `len`) is ambiguous. Returns `(offset, length)` pairs in
/// order; gap bases belong to no segment. A sequence without gaps yields a
/// single segment covering all of it (an empty sequence yields none).
pub fn gap_free_segments_by<F: Fn(usize) -> bool>(
    len: usize,
    min_run: usize,
    is_amb: F,
) -> Vec<(usize, usize)> {
    let min_run = min_run.max(1);
    let mut segs = Vec::new();
    let mut seg_start = 0usize;
    let mut i = 0usize;
    while i < len {
        if is_amb(i) {
            let run_start = i;
            while i < len && is_amb(i) {
                i += 1;
            }
            if i - run_start >= min_run {
                if run_start > seg_start {
                    segs.push((seg_start, run_start - seg_start));
                }
                seg_start = i;
            }
        } else {
            i += 1;
        }
    }
    if len > seg_start {
        segs.push((seg_start, len - seg_start));
    }
    segs
}

/// [`gap_free_segments_by`] over raw ASCII input.
pub fn gap_free_segments(ascii: &[u8], min_run: usize) -> Vec<(usize, usize)> {
    gap_free_segments_by(ascii.len(), min_run, |i| is_ambiguous_byte(ascii[i]))
}

/// Reverse-complement a slice of base codes.
///
/// Complementing in this encoding is simply `3 - code` (A<->T, C<->G); the
/// sequence is then reversed.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| 3 - (c & 0b11)).collect()
}

/// `true` if the base at `i` is G or C.
#[inline]
pub fn is_gc(seq: &[u8], i: usize) -> bool {
    let c = seq[i];
    c == C || c == G
}

/// GC fraction over the inclusive range `[a, b]`.
pub fn gc_content(seq: &[u8], a: usize, b: usize) -> f64 {
    if b < a {
        return 0.0;
    }
    let mut gc = 0.0f64;
    let mut tot = 0.0f64;
    for i in a..=b {
        if is_gc(seq, i) {
            gc += 1.0;
        }
        tot += 1.0;
    }
    if tot == 0.0 {
        0.0
    } else {
        gc / tot
    }
}

/// Index of the `len`-mer starting at `pos`, base-4 little-endian.
///
/// Same construction as Prodigal's `mer_ndx`: each base contributes
/// `code << (2*i)`. For `len == 6` the result is in `0..4096`. Positions that
/// would run off the end of the sequence are treated as `A` (code 0).
#[inline]
pub fn mer_ndx(len: usize, seq: &[u8], pos: usize) -> usize {
    let mut ndx = 0usize;
    let slen = seq.len();
    for i in 0..len {
        let p = pos + i;
        let code = if p < slen { seq[p] as usize } else { 0 };
        ndx |= code << (2 * i);
    }
    ndx
}

/// Render a `len`-mer index back to ASCII text (handy for motif reporting).
pub fn mer_text(len: usize, ndx: usize) -> String {
    if len == 0 {
        return "None".to_string();
    }
    let mut s = String::with_capacity(len);
    for i in 0..len {
        let code = ((ndx >> (2 * i)) & 0b11) as u8;
        s.push(decode_base(code) as char);
    }
    s
}

/// Compute the background frequency of every `len`-mer across both strands,
/// matching Prodigal's `calc_mer_bg`. The returned vector has `4^len` entries.
pub fn calc_mer_bg(len: usize, fwd: &[u8], rev: &[u8]) -> Vec<f64> {
    let size = 1usize << (2 * len);
    let mut counts = vec![0u64; size];
    let slen = fwd.len();
    if slen < len {
        return vec![0.0; size];
    }
    let mut glob = 0u64;
    for i in 0..=(slen - len) {
        counts[mer_ndx(len, fwd, i)] += 1;
        counts[mer_ndx(len, rev, i)] += 1;
        glob += 2;
    }
    let denom = glob as f64;
    counts
        .into_iter()
        .map(|c| if denom > 0.0 { c as f64 / denom } else { 0.0 })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_roundtrip() {
        let s = Seq::from_ascii(b"ATGCAAGTAG");
        // revcomp of revcomp == original forward
        assert_eq!(revcomp(&s.rev), s.fwd);
        // ambiguous bases keep the placeholder on both strands (Prodigal)
        let s = Seq::from_ascii(b"ATGCAANTAG");
        assert_eq!(s.fwd[6], AMBIGUOUS_CODE);
        assert_eq!(s.rev[3], AMBIGUOUS_CODE);
        assert!(s.amb_rev(3));
    }

    #[test]
    fn gc_test() {
        // GCGC -> all GC
        let s = Seq::from_ascii(b"GCGC");
        assert!((s.gc - 1.0).abs() < 1e-9);
        for i in 0..4 {
            assert!(is_gc(&s.fwd, i));
        }
        // ATAT -> no GC
        let s2 = Seq::from_ascii(b"ATAT");
        assert!(s2.gc.abs() < 1e-9);
    }

    #[test]
    fn gc_ignores_ambiguous_codes() {
        // Regression: S/B/V used to count as GC in the numerator but not the
        // denominator, so "ASSS" reported GC = 3.0.
        let s = Seq::from_ascii(b"ASSS");
        assert!(s.gc <= 1.0 && s.gc.abs() < 1e-9, "gc = {}", s.gc);
        let s = Seq::from_ascii(b"ACNNNNGT");
        assert!((s.gc - 0.5).abs() < 1e-9);
        assert!((s.gc_prodigal() - 0.25).abs() < 1e-9);
    }

    #[test]
    fn ambiguous_runs_are_half_open_and_include_trailing_runs() {
        // runs of 3 (leading), 2 (too short), 4 (internal), 3 (trailing)
        let s = Seq::from_ascii(b"NNNACNNGTRYKSAAANNN");
        assert_eq!(ambiguous_runs(&s, 3), vec![(0, 3), (9, 13), (16, 19)]);
        assert_eq!(ambiguous_runs(&s, 4), vec![(9, 13)]);
        assert!(ambiguous_runs(&Seq::from_ascii(b"ACGT"), 1).is_empty());
    }

    #[test]
    fn ambiguous_bases_are_flagged_not_read_as_a() {
        let s = Seq::from_ascii(b"TNNACG");
        assert!(s.n[1] && s.n[2]);
        assert_ne!(s.fwd[1], A, "N must not be stored as A");
        assert!(s.codon_has_amb_fwd(0));
        assert!(!s.codon_has_amb_fwd(3));
        // reverse strand: rev[i] pairs with fwd[len-1-i]
        assert!(s.amb_rev(3) && s.amb_rev(4));
        assert!(!s.amb_rev(0));
        assert_eq!(s.to_ascii_upper(), b"TNNACG".to_vec());
    }

    #[test]
    fn segments_split_only_at_long_runs() {
        let mut v = b"ACGTACGT".to_vec();
        v.extend_from_slice(&[b'N'; GAP_MIN_RUN]);
        v.extend_from_slice(b"ACGNNACG"); // short run stays inline
        let segs = gap_free_segments(&v, GAP_MIN_RUN);
        assert_eq!(segs, vec![(0, 8), (8 + GAP_MIN_RUN, 8)]);
        assert_eq!(gap_free_segments(b"NNNNNNNNNNNN", GAP_MIN_RUN), vec![]);
        assert_eq!(gap_free_segments(b"ACGT", GAP_MIN_RUN), vec![(0, 4)]);
    }

    #[test]
    fn mer_index_matches_base4() {
        // ATG = A(0) T(3) G(2) -> 0<<0 | 3<<2 | 2<<4 = 12 + 32 = 44
        let s = Seq::from_ascii(b"ATG");
        assert_eq!(mer_ndx(3, &s.fwd, 0), (3 << 2) | (2 << 4));
        assert_eq!(&mer_text(3, mer_ndx(3, &s.fwd, 0)), "ATG");
    }
}
