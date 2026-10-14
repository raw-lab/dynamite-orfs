//! A small, engine-independent open-reading-frame enumerator.
//!
//! Several callers (GeneMark, AUGUSTUS, the six-frame translator) need the
//! same primitive: walk all six frames of a contig and emit every maximal
//! stop-to-stop region, together with the best (most upstream) start codon
//! inside it. This module provides that, parameterised by genetic code,
//! minimum length, and whether ORFs running off the contig ends are allowed
//! (partial genes).
//!
//! Partial flags follow DYNAMITE's (= Prodigal's) **genomic** convention:
//! `partial_left` / `partial_right` refer to the lower / higher coordinate end,
//! on both strands. Codons containing an ambiguous base (N / IUPAC) are never
//! starts or stops.

use crate::gencode::GeneticCode;
use crate::seq::Seq;

/// A candidate ORF in 1-based inclusive coordinates (stop codon included for
/// complete genes). `strand` is +1 / -1.
#[derive(Clone, Debug)]
pub struct OrfCand {
    pub begin: usize,
    pub end: usize,
    pub strand: i8,
    /// 0=ATG, 1=GTG, 2=TTG, 3=edge/none.
    pub start_type: u8,
    /// Incomplete at the lower-coordinate end (genomic convention).
    pub partial_left: bool,
    /// Incomplete at the higher-coordinate end (genomic convention).
    pub partial_right: bool,
}

/// Enumerate ORFs over both strands — the classic six-frame ORF listing.
///
/// * `min_len` — minimum ORF length in nucleotides (stop included).
/// * `require_start` — if true, only ORFs that begin with a start codon are
///   reported.
/// * `allow_partial_ends` — if true (and `require_start` is false) every
///   stop-to-stop region is listed, including those without a start codon
///   (reported 5'-partial) and the fragments that run off either contig end.
///
/// Gene callers should use [`find_gene_orfs`], which never treats an interior
/// start-less region as a gene.
pub fn find_orfs(
    seq: &Seq,
    code: &GeneticCode,
    min_len: usize,
    require_start: bool,
    allow_partial_ends: bool,
) -> Vec<OrfCand> {
    scan_both(seq, code, min_len, require_start, allow_partial_ends, true)
}

/// Enumerate candidate **genes** over both strands. Like [`find_orfs`], but a
/// start-less ORF is only allowed where it is genuinely truncated: between the
/// strand's 5' contig edge and the first in-frame stop (and only when
/// `allow_partial_ends`). Interior stop-to-stop regions without a start codon
/// are never reported — the upstream in-frame stop proves they cannot be a
/// truncated gene (they used to be emitted as partial genes in mid-contig).
pub fn find_gene_orfs(
    seq: &Seq,
    code: &GeneticCode,
    min_len: usize,
    require_start: bool,
    allow_partial_ends: bool,
) -> Vec<OrfCand> {
    scan_both(seq, code, min_len, require_start, allow_partial_ends, false)
}

fn scan_both(
    seq: &Seq,
    code: &GeneticCode,
    min_len: usize,
    require_start: bool,
    allow_partial_ends: bool,
    interior_startless: bool,
) -> Vec<OrfCand> {
    let mut out = Vec::new();
    let n = seq.len();
    if n < 3 {
        return out;
    }
    let ramb = seq.n_rev();
    let opts = ScanOpts {
        min_len,
        require_start,
        allow_partial_ends,
        interior_startless,
    };
    // Forward strand uses seq.fwd (2-bit codes); reverse uses seq.rev.
    scan_strand(&seq.fwd, &seq.n, n, 1, code, opts, &mut out);
    scan_strand(&seq.rev, &ramb, n, -1, code, opts, &mut out);
    out
}

#[derive(Clone, Copy)]
struct ScanOpts {
    min_len: usize,
    require_start: bool,
    allow_partial_ends: bool,
    interior_startless: bool,
}

fn scan_strand(
    s: &[u8],
    amb: &[bool],
    n: usize,
    strand: i8,
    code: &GeneticCode,
    o: ScanOpts,
    out: &mut Vec<OrfCand>,
) {
    let clean = |i: usize| !(amb[i] || amb[i + 1] || amb[i + 2]);
    for frame in 0..3usize {
        // `region_start` is the first codon of the current stop-to-stop region;
        // `edge_region` is true while that region still begins at the strand's
        // 5' edge (no in-frame stop seen yet) — only there is a start-less ORF
        // genuinely truncated rather than simply non-coding.
        let mut region_start = frame;
        let mut edge_region = true;
        let mut first_start: Option<usize> = None;
        let startless_ok =
            |edge: bool| o.allow_partial_ends && !o.require_start && (edge || o.interior_startless);
        let mut i = frame;
        while i + 3 <= n {
            let ok = clean(i);
            let (b0, b1, b2) = (s[i], s[i + 1], s[i + 2]);
            let is_start = ok && code.is_start(b0, b1, b2);
            let is_stop = ok && code.is_stop(b0, b1, b2);
            if is_start && first_start.is_none() {
                first_start = Some(i);
            }
            if is_stop {
                let begin0 = match first_start {
                    Some(p) => Some(p),
                    None if startless_ok(edge_region) => Some(region_start),
                    None => None,
                };
                if let Some(b0pos) = begin0 {
                    let stop_end = i + 2; // 0-based inclusive end of stop codon
                    let len = stop_end + 1 - b0pos;
                    if len >= o.min_len {
                        push_orf(s, n, strand, b0pos, stop_end, first_start.is_none(), false, out);
                    }
                }
                first_start = None;
                edge_region = false;
                region_start = i + 3;
            }
            i += 3;
        }
        // Trailing region with no terminal stop (3'-partial), if allowed.
        if o.allow_partial_ends {
            let begin0 = match first_start {
                Some(p) => Some(p),
                None if startless_ok(edge_region) => Some(region_start),
                None => None,
            };
            if let Some(b0pos) = begin0 {
                // Extend to the last full codon in this frame.
                let last = frame + ((n - frame) / 3) * 3;
                if last >= b0pos + 3 {
                    let stop_end = last - 1;
                    let len = stop_end + 1 - b0pos;
                    if len >= o.min_len {
                        push_orf(s, n, strand, b0pos, stop_end, first_start.is_none(), true, out);
                    }
                }
            }
        }
    }
}

/// Record one ORF given strand-local 0-based `[begin0, end0]` and its
/// strand-relative 5' / 3' partial status.
#[allow(clippy::too_many_arguments)]
fn push_orf(
    s: &[u8],
    n: usize,
    strand: i8,
    begin0: usize,
    end0: usize,
    partial_5p: bool,
    partial_3p: bool,
    out: &mut Vec<OrfCand>,
) {
    // Determine start type from the first codon (on the strand's coordinates).
    let st = if partial_5p {
        3
    } else {
        GeneticCode::start_type(s[begin0], s[begin0 + 1], s[begin0 + 2])
    };
    // Convert strand-local 0-based [begin0,end0] to global 1-based [begin,end]
    // and the 5'/3' flags to the genomic left/right convention.
    let (begin, end, pl, pr) = if strand == 1 {
        (begin0 + 1, end0 + 1, partial_5p, partial_3p)
    } else {
        // reverse strand: position p on rev corresponds to global n-1-p; the
        // 5' end of a minus-strand ORF is its right (higher) coordinate.
        let g_begin = n - 1 - end0; // smaller global coordinate
        let g_end = n - 1 - begin0;
        (g_begin + 1, g_end + 1, partial_3p, partial_5p)
    };
    out.push(OrfCand {
        begin,
        end,
        strand,
        start_type: st,
        partial_left: pl,
        partial_right: pr,
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code() -> GeneticCode {
        GeneticCode::new(11).unwrap()
    }

    #[test]
    fn no_startless_orfs_in_the_interior() {
        // stop | 40 start-less codons | stop | ATG ... stop
        let mut v = b"AAATAA".to_vec();
        v.extend(std::iter::repeat(*b"GCA").take(40).flatten());
        v.extend_from_slice(b"TAA");
        v.extend_from_slice(b"ATG");
        v.extend(std::iter::repeat(*b"GCA").take(40).flatten());
        v.extend_from_slice(b"TAA");
        let s = Seq::from_ascii(&v);
        let orfs = find_gene_orfs(&s, &code(), 60, false, true);
        let fwd0: Vec<&OrfCand> = orfs
            .iter()
            .filter(|o| o.strand == 1 && (o.begin - 1) % 3 == 0)
            .collect();
        // Only the ATG ORF in frame 0; the start-less interior region is not a gene.
        assert_eq!(fwd0.len(), 1, "{:?}", fwd0);
        assert!(!fwd0[0].partial_left && !fwd0[0].partial_right);
    }

    #[test]
    fn orf_listing_keeps_every_stop_to_stop_region() {
        // `find_orfs` is the six-frame ORF lister: interior start-less regions
        // are listed (5'-partial), as `translate --orfs` documents.
        let mut v = b"AAATAA".to_vec();
        v.extend(std::iter::repeat(*b"GCA").take(40).flatten());
        v.extend_from_slice(b"TAA");
        let s = Seq::from_ascii(&v);
        let orfs = find_orfs(&s, &code(), 60, false, true);
        assert!(orfs.iter().any(|o| o.strand == 1 && o.begin == 7 && o.end == 129 && o.partial_left));
        // ...but never with require_start
        let orfs = find_orfs(&s, &code(), 60, true, false);
        assert!(!orfs.iter().any(|o| o.strand == 1 && o.begin == 7));
    }

    #[test]
    fn reverse_edge_orf_flags_are_genomic() {
        // A minus-strand ORF whose start runs off the right edge must be
        // flagged on the right (higher-coordinate) end.
        let mut fwd = b"TTA".to_vec(); // revcomp TAA -> stop at the left end of the minus ORF
        fwd.extend(std::iter::repeat(*b"TGC").take(40).flatten());
        let s = Seq::from_ascii(&fwd);
        let orfs = find_orfs(&s, &code(), 60, false, true);
        let o = orfs
            .iter()
            .find(|o| o.strand == -1 && o.begin == 1)
            .expect("minus-strand edge ORF");
        assert!(!o.partial_left, "stop present at the left end");
        assert!(o.partial_right, "start missing at the right (5') end");
        assert_eq!(o.start_type, 3);
    }

    #[test]
    fn ambiguous_codons_are_not_stops() {
        // ATG ... T NN ... TAA : `TNN` must not terminate the ORF.
        let mut v = b"ATG".to_vec();
        v.extend(std::iter::repeat(*b"GCA").take(30).flatten());
        v.extend_from_slice(b"TNN");
        v.extend(std::iter::repeat(*b"GCA").take(30).flatten());
        v.extend_from_slice(b"TAA");
        let s = Seq::from_ascii(&v);
        let orfs = find_orfs(&s, &code(), 60, true, false);
        assert!(orfs
            .iter()
            .any(|o| o.strand == 1 && o.begin == 1 && o.end == v.len()));
    }
}
