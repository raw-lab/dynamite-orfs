//! Engine-independent post-processing applied to every gene call:
//! **assembly-gap segmentation**, **gap-codon trimming** and **codon-based
//! inference of partial ends**.
//!
//! * Contigs are split at assembly gaps (runs of at least
//!   [`crate::seq::GAP_MIN_RUN`] ambiguous bases) before calling, so no gene can
//!   span a gap or start/end inside one; a gene that runs into a gap is
//!   truncated there and reported partial at that end — exactly like a gene
//!   that runs off a contig edge.
//! * Codons containing an ambiguous base at either end of a gene are trimmed
//!   (in frame) and that end is marked partial.
//! * Partial ends are inferred from the sequence itself: a gene is 5'-partial
//!   unless it begins with a start codon of the active genetic code, and
//!   3'-partial unless it ends with a stop codon. Following Prodigal, a start
//!   lying in the first codon of an *open* edge (contig end or gap boundary)
//!   is also treated as partial — the real start may lie upstream, outside the
//!   sequence.
//! * Partial flags use Prodigal's genomic `partial=XY` convention (left/right
//!   end), the start type becomes `Edge` for 5'-partial genes, and the protein
//!   is re-translated (`M` only for a complete 5' end, `X` for ambiguous
//!   codons, no stray `*`).

use crate::gencode::GeneticCode;
use crate::gene::Gene;
use crate::seq::{gap_free_segments, gap_free_segments_by, Seq, GAP_MIN_RUN};

/// Minimum length (nt) a gene must keep after terminal gap codons are trimmed
/// (Prodigal's minimum edge-gene length). Shorter remnants are dropped.
pub const MIN_TRIMMED_LEN: usize = 60;

/// Where a gene's nucleotide / protein strings come from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SeqSource {
    /// The gene is an exact, codon-aligned genome interval: trim terminal gap
    /// codons, infer partial ends, and re-extract / re-translate from the
    /// genome.
    Genome,
    /// The engine supplied its own (frameshift-corrected) sequences —
    /// FragGeneScan. Partial ends are inferred from the genome codons at the
    /// gene's coordinates, but nothing is trimmed or re-translated.
    Engine,
}

/// Which codons count as a start codon when inferring a gene's 5' end.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum StartRule {
    /// The genetic code's own start codons (Prodigal's per-table rules) —
    /// for engines that only ever open genes at those codons.
    Table,
    /// ATG / GTG / TTG under every table, plus the table's own starts — for
    /// engines whose start set does not depend on the table (PHANOTATE,
    /// FragGeneScan).
    Bacterial,
}

/// `(offset, len)` of each gap-free segment of a raw ASCII contig.
pub fn segments(ascii: &[u8]) -> Vec<(usize, usize)> {
    gap_free_segments(ascii, GAP_MIN_RUN)
}

/// `(offset, len)` of each gap-free segment of an encoded contig.
pub fn segments_of_seq(seq: &Seq) -> Vec<(usize, usize)> {
    gap_free_segments_by(seq.len(), GAP_MIN_RUN, |i| seq.n[i])
}

/// The codon at coding position `k` (0-based, 5'->3') of a gene, as
/// `(b1, b2, b3, any_ambiguous)`. `None` if it falls outside the sequence.
fn coding_codon(seq: &Seq, g: &Gene, k: usize) -> Option<(u8, u8, u8, bool)> {
    let n = seq.len();
    if g.begin == 0 || g.end > n || g.begin > g.end {
        return None;
    }
    if g.strand >= 0 {
        let p = g.begin - 1 + 3 * k;
        if p + 3 > g.end {
            return None;
        }
        Some((
            seq.fwd[p],
            seq.fwd[p + 1],
            seq.fwd[p + 2],
            seq.codon_has_amb_fwd(p),
        ))
    } else {
        // 5' end of a minus-strand gene is `end`; rev index of that base:
        let r = n - g.end + 3 * k;
        if r + 3 > n - g.begin + 1 {
            return None;
        }
        Some((
            seq.rev[r],
            seq.rev[r + 1],
            seq.rev[r + 2],
            seq.codon_has_amb_rev(r),
        ))
    }
}

/// The codon ending at the gene's 3' coordinate (genome orientation-aware).
fn last_codon(seq: &Seq, g: &Gene) -> Option<(u8, u8, u8, bool)> {
    let n = seq.len();
    if g.begin == 0 || g.end > n || g.end < g.begin + 2 {
        return None;
    }
    if g.strand >= 0 {
        let p = g.end - 3;
        Some((
            seq.fwd[p],
            seq.fwd[p + 1],
            seq.fwd[p + 2],
            seq.codon_has_amb_fwd(p),
        ))
    } else {
        let r = n - g.begin - 2; // rev index of the 3' codon's first base
        Some((
            seq.rev[r],
            seq.rev[r + 1],
            seq.rev[r + 2],
            seq.codon_has_amb_rev(r),
        ))
    }
}

#[inline]
fn is_start(c: Option<(u8, u8, u8, bool)>, code: &GeneticCode, rule: StartRule) -> bool {
    match c {
        Some((a, b, d, false)) => {
            code.is_start(a, b, d)
                || (rule == StartRule::Bacterial && GeneticCode::start_type(a, b, d) < 3)
        }
        _ => false,
    }
}

#[inline]
fn is_stop(c: Option<(u8, u8, u8, bool)>, code: &GeneticCode) -> bool {
    matches!(c, Some((a, b, d, false)) if code.is_stop(a, b, d))
}

/// `true` if the gene's 5' codon begins within the first 3 bases of its
/// strand (an open edge: the real start may lie upstream of the sequence).
#[inline]
fn starts_at_edge(seq: &Seq, g: &Gene) -> bool {
    if g.strand >= 0 {
        g.begin <= 3
    } else {
        g.end + 2 >= seq.len()
    }
}

fn gc_of_interval(seq: &Seq, g: &Gene) -> f64 {
    let b = g.begin.saturating_sub(1);
    let e = g.end.min(seq.len());
    let mut gc = 0usize;
    let mut tot = 0usize;
    for i in b..e {
        if !seq.n[i] {
            tot += 1;
            if seq.fwd[i] == crate::seq::C || seq.fwd[i] == crate::seq::G {
                gc += 1;
            }
        }
    }
    if tot == 0 {
        0.0
    } else {
        gc as f64 / tot as f64
    }
}

/// Finalise one gene in place. Returns `false` if it should be discarded
/// (nothing but gap codons, or too short once they are trimmed).
pub fn finalize_gene(
    g: &mut Gene,
    seq: &Seq,
    code: &GeneticCode,
    closed: bool,
    source: SeqSource,
    starts: StartRule,
) -> bool {
    finalize_gene_ends(g, seq, code, closed, closed, source, starts)
}

/// [`finalize_gene`] with each end of the sequence open or closed on its own:
/// Prodigal's open-edge rule (a start in the first codon of an open edge is
/// partial) applies only to an open end.
pub fn finalize_gene_ends(
    g: &mut Gene,
    seq: &Seq,
    code: &GeneticCode,
    closed_left: bool,
    closed_right: bool,
    source: SeqSource,
    starts: StartRule,
) -> bool {
    if g.begin == 0 || g.end > seq.len() || g.begin > g.end {
        return false;
    }
    let len = g.end - g.begin + 1;
    let codon_aligned = len % 3 == 0;
    let mut p5 = g.partial_5p();
    let mut trimmed = false;
    let mut trimmed_3p = false;

    if source == SeqSource::Genome && codon_aligned {
        // ---- trim terminal codons that touch an ambiguous base ----
        let m = len / 3;
        let mut lead = 0usize;
        while lead < m && coding_codon(seq, g, lead).map(|c| c.3).unwrap_or(true) {
            lead += 1;
        }
        if lead == m {
            return false; // nothing but gap codons
        }
        let mut trail = 0usize;
        while trail < m - lead
            && coding_codon(seq, g, m - 1 - trail)
                .map(|c| c.3)
                .unwrap_or(true)
        {
            trail += 1;
        }
        if lead > 0 || trail > 0 {
            if (m - lead - trail) * 3 < MIN_TRIMMED_LEN {
                return false;
            }
            if g.strand >= 0 {
                g.begin += 3 * lead;
                g.end -= 3 * trail;
            } else {
                g.end -= 3 * lead;
                g.begin += 3 * trail;
            }
            trimmed = true;
            p5 |= lead > 0;
            trimmed_3p = trail > 0;
        }
    }

    // ---- infer partial ends from the codons themselves ----
    let first = coding_codon(seq, g, 0);
    let has_start = is_start(first, code, starts);
    let has_stop = is_stop(last_codon(seq, g), code);
    // 5': the engine's verdict, a trimmed gap, a missing start codon, or a
    // start sitting in the first codon of an open edge (Prodigal's rule).
    let open_5p_edge = if g.strand >= 0 { !closed_left } else { !closed_right };
    let p5 = p5 || !has_start || (open_5p_edge && starts_at_edge(seq, g));
    // 3': a real stop codon is a definitive end; anything else is truncated.
    let p3 = !has_stop || trimmed_3p;
    g.set_partial_5p_3p(p5, p3);
    g.start_type = if p5 {
        3
    } else {
        match first {
            Some((a, b, d, false)) => GeneticCode::start_type(a, b, d),
            _ => 3,
        }
    };

    if source == SeqSource::Genome && codon_aligned {
        g.fill_sequences(seq, code);
        if trimmed {
            g.gc = gc_of_interval(seq, g);
        }
    }
    true
}

/// Finalise every gene of one contig / segment (see the module docs). Genes
/// that are reduced to nothing by gap trimming are removed.
pub fn finalize_genes(
    genes: Vec<Gene>,
    seq: &Seq,
    code: &GeneticCode,
    closed: bool,
    source: SeqSource,
    starts: StartRule,
) -> Vec<Gene> {
    finalize_genes_ends(genes, seq, code, closed, closed, source, starts)
}

/// [`finalize_genes`] with each end of the sequence open or closed on its own
/// (see [`finalize_gene_ends`]).
pub fn finalize_genes_ends(
    genes: Vec<Gene>,
    seq: &Seq,
    code: &GeneticCode,
    closed_left: bool,
    closed_right: bool,
    source: SeqSource,
    starts: StartRule,
) -> Vec<Gene> {
    genes
        .into_iter()
        .filter_map(|mut g| {
            if finalize_gene_ends(&mut g, seq, code, closed_left, closed_right, source, starts) {
                Some(g)
            } else {
                None
            }
        })
        .collect()
}

/// Shift a segment's genes into the parent contig's coordinates.
pub fn shift_genes(genes: &mut [Gene], offset: usize, contig: &str) {
    for g in genes.iter_mut() {
        g.begin += offset;
        g.end += offset;
        if g.contig != contig {
            g.contig = contig.to_string();
        }
    }
}

/// Order genes tagged with the index of the input record they came from:
/// by record (input order), then coordinate, renumbering the 1-based `id`
/// within each record. Unlike [`order_and_renumber`] this keeps records that
/// share an id apart, as 1.0.0 did.
pub fn order_by_record(mut tagged: Vec<(usize, Gene)>) -> Vec<Gene> {
    tagged.sort_by(|(ra, a), (rb, b)| {
        ra.cmp(rb)
            .then(a.begin.cmp(&b.begin))
            .then(a.end.cmp(&b.end))
            .then(b.strand.cmp(&a.strand))
    });
    let mut cur: Option<usize> = None;
    let mut k = 0usize;
    tagged
        .into_iter()
        .map(|(r, mut g)| {
            if cur != Some(r) {
                cur = Some(r);
                k = 0;
            }
            k += 1;
            g.id = k;
            g
        })
        .collect()
}

/// Sort genes by contig (in the given order), then coordinate, and renumber
/// the 1-based `id` within each contig. Genes of records that share an id are
/// merged; [`order_by_record`] keeps them apart.
pub fn order_and_renumber(genes: &mut [Gene], contig_order: &[String]) {
    use std::collections::HashMap;
    let mut rank: HashMap<&str, usize> = HashMap::new();
    for (i, c) in contig_order.iter().enumerate() {
        rank.entry(c.as_str()).or_insert(i);
    }
    genes.sort_by(|a, b| {
        let ra = rank.get(a.contig.as_str()).copied().unwrap_or(usize::MAX);
        let rb = rank.get(b.contig.as_str()).copied().unwrap_or(usize::MAX);
        ra.cmp(&rb)
            .then(a.contig.cmp(&b.contig))
            .then(a.begin.cmp(&b.begin))
            .then(a.end.cmp(&b.end))
            .then(b.strand.cmp(&a.strand))
    });
    let mut cur: Option<String> = None;
    let mut k = 0usize;
    for g in genes.iter_mut() {
        if cur.as_deref() != Some(g.contig.as_str()) {
            cur = Some(g.contig.clone());
            k = 0;
        }
        k += 1;
        g.id = k;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene::Engine;

    fn gene(begin: usize, end: usize, strand: i8) -> Gene {
        Gene {
            contig: "c".into(),
            id: 1,
            begin,
            end,
            strand,
            start_type: 0,
            rbs_motif: "None".into(),
            rbs_spacer: "None".into(),
            partial_left: false,
            partial_right: false,
            score: 0.0,
            cscore: 0.0,
            sscore: 0.0,
            conf: 0.0,
            gc: 0.0,
            trans_table: 11,
            engine: Engine::Dynamic,
            aa: String::new(),
            nt: String::new(),
        }
    }

    fn body(n: usize) -> Vec<u8> {
        std::iter::repeat(*b"GCA").take(n).flatten().collect()
    }

    #[test]
    fn complete_gene_stays_complete() {
        let mut v = b"CCCCCC".to_vec();
        v.extend_from_slice(b"ATG");
        v.extend(body(30));
        v.extend_from_slice(b"TAA");
        v.extend_from_slice(b"CCCCCC");
        let s = Seq::from_ascii(&v);
        let code = GeneticCode::new(11).unwrap();
        let mut g = gene(7, 7 + 3 + 90 + 3 - 1, 1);
        assert!(finalize_gene(&mut g, &s, &code, false, SeqSource::Genome, StartRule::Table));
        assert_eq!(g.partial_str(), "00");
        assert_eq!(g.start_codon_str(), "ATG");
        assert!(g.aa.starts_with('M') && !g.aa.contains('*'));
    }

    #[test]
    fn missing_start_and_stop_are_inferred() {
        // forward ORF claimed complete by an engine but lacking both signals
        let mut v = b"CCCCCC".to_vec();
        v.extend(body(40));
        v.extend_from_slice(b"CCCCCC");
        let s = Seq::from_ascii(&v);
        let code = GeneticCode::new(11).unwrap();
        let mut g = gene(7, 7 + 120 - 1, 1);
        assert!(finalize_gene(&mut g, &s, &code, false, SeqSource::Genome, StartRule::Table));
        assert_eq!(g.partial_str(), "11");
        assert_eq!(g.start_codon_str(), "Edge");
    }

    #[test]
    fn reverse_gene_partial_is_genomic() {
        // minus-strand gene: stop at the left end, start missing at the right
        // (5') end -> partial=01 (left complete, right partial).
        let mut coding = body(40); // no start codon
        coding.extend_from_slice(b"TAA");
        let fwd: Vec<u8> = coding
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        let mut v = b"CCC".to_vec();
        v.extend_from_slice(&fwd);
        v.extend_from_slice(b"CCC");
        let s = Seq::from_ascii(&v);
        let code = GeneticCode::new(11).unwrap();
        let mut g = gene(4, 3 + fwd.len(), -1);
        assert!(finalize_gene(&mut g, &s, &code, false, SeqSource::Genome, StartRule::Table));
        assert_eq!(g.partial_str(), "01");
        assert!(g.partial_5p() && !g.partial_3p());
        assert_eq!(g.start_codon_str(), "Edge");
        assert!(!g.aa.ends_with('*'), "stop must be stripped: {}", g.aa);
    }

    #[test]
    fn gap_codons_are_trimmed_and_flagged() {
        // ATG + 30 codons + 'T' + NNNNN...: the engine's ORF ran into the gap.
        let mut v = b"ATG".to_vec();
        v.extend(body(30));
        v.extend_from_slice(b"TNNNNNNNN");
        let s = Seq::from_ascii(&v);
        let code = GeneticCode::new(11).unwrap();
        let mut g = gene(1, 3 + 90 + 9, 1);
        assert!(finalize_gene(&mut g, &s, &code, true, SeqSource::Genome, StartRule::Table));
        assert_eq!((g.begin, g.end), (1, 93), "gap codons trimmed in frame");
        assert_eq!(g.partial_str(), "01");
        assert!(!g.nt.contains('N') && !g.aa.contains('X'));
    }

    #[test]
    fn edge_start_is_partial_unless_closed() {
        let mut v = b"ATG".to_vec();
        v.extend(body(30));
        v.extend_from_slice(b"TAA");
        let s = Seq::from_ascii(&v);
        let code = GeneticCode::new(11).unwrap();
        let mut open = gene(1, v.len(), 1);
        assert!(finalize_gene(&mut open, &s, &code, false, SeqSource::Genome, StartRule::Table));
        assert_eq!(open.partial_str(), "10");
        let mut closed = gene(1, v.len(), 1);
        assert!(finalize_gene(&mut closed, &s, &code, true, SeqSource::Genome, StartRule::Table));
        assert_eq!(closed.partial_str(), "00");
    }

    #[test]
    fn bacterial_start_rule_accepts_gtg_under_table_1() {
        // PHANOTATE opens genes at GTG under any table; under table 1 (ATG
        // only for Prodigal) such a gene is still complete for that engine.
        let mut v = b"CCCCCC".to_vec();
        v.extend_from_slice(b"GTG");
        v.extend(body(30));
        v.extend_from_slice(b"TAA");
        v.extend_from_slice(b"CCCCCC");
        let s = Seq::from_ascii(&v);
        let code = GeneticCode::new(1).unwrap();
        let mut g = gene(7, 7 + 3 + 90 + 3 - 1, 1);
        assert!(finalize_gene(&mut g, &s, &code, false, SeqSource::Genome, StartRule::Bacterial));
        assert_eq!((g.partial_str().as_str(), g.start_codon_str()), ("00", "GTG"));
        assert!(g.aa.starts_with('M'));
        let mut g = gene(7, 7 + 3 + 90 + 3 - 1, 1);
        assert!(finalize_gene(&mut g, &s, &code, false, SeqSource::Genome, StartRule::Table));
        assert_eq!(g.partial_str(), "10");
    }

    #[test]
    fn renumber_keeps_input_contig_order() {
        let mut gs = vec![gene(50, 200, 1), gene(10, 100, 1), gene(5, 90, 1)];
        gs[0].contig = "contig10".into();
        gs[1].contig = "contig2".into();
        gs[2].contig = "contig10".into();
        order_and_renumber(&mut gs, &["contig2".to_string(), "contig10".to_string()]);
        let order: Vec<(&str, usize, usize)> = gs
            .iter()
            .map(|g| (g.contig.as_str(), g.begin, g.id))
            .collect();
        assert_eq!(
            order,
            vec![("contig2", 10, 1), ("contig10", 5, 1), ("contig10", 50, 2)]
        );
    }
}
