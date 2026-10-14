//! Greedy overlap resolution for the ORF-enumeration engines.
//!
//! GeneMark and AUGUSTUS enumerate every positive-scoring open reading frame
//! in all six frames. Real implementations of those tools run a hidden-state
//! Viterbi decode that yields **one** consistent gene parse and resolves
//! overlaps; DYNAMITE approximates that with a greedy maximum-score
//! selection: take ORFs in descending coding score and keep one only if it does
//! not overlap an already-kept gene (on the same contig) by more than
//! `max_overlap` nucleotides. This collapses nested / duplicate / wrong-strand
//! ORFs while still permitting the short operon-style overlaps seen in compact
//! genomes. It is a heuristic, not a true Viterbi.

use crate::gene::Gene;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};

/// Resolve overlaps greedily by descending coding score, then restore genomic
/// order (contigs in the order they first appear in `genes`, which is the
/// input order) and renumber gene IDs per contig.
pub fn resolve_overlaps(genes: Vec<Gene>, max_overlap: i64) -> Vec<Gene> {
    if genes.len() < 2 {
        return genes;
    }
    if max_overlap < 0 {
        // Degenerate limit: every pair "clashes" (even disjoint genes), so the
        // interval search below does not apply. Keep the original behaviour.
        return resolve_overlaps_naive(genes, max_overlap);
    }
    // Contig order of first appearance (the input order). Sorting by contig
    // *name* used to put "contig10" before "contig2".
    let mut contig_rank: HashMap<String, usize> = HashMap::new();
    for g in &genes {
        let next = contig_rank.len();
        contig_rank.entry(g.contig.clone()).or_insert(next);
    }

    // Highest coding score first; break ties by longer gene, then position
    // (deterministic).
    let mut order: Vec<usize> = (0..genes.len()).collect();
    order.sort_by(|&a, &b| {
        let (ga, gb) = (&genes[a], &genes[b]);
        gb.cscore
            .partial_cmp(&ga.cscore)
            .unwrap_or(Ordering::Equal)
            .then((gb.end - gb.begin).cmp(&(ga.end - ga.begin)))
            .then(a.cmp(&b))
    });

    // Kept intervals per contig, keyed by begin (several genes may share a
    // begin, so the value is a list of ends). `max_len` bounds how far back a
    // kept interval can start and still reach a candidate.
    let mut kept_iv: HashMap<usize, BTreeMap<usize, Vec<usize>>> = HashMap::new();
    let mut max_len: HashMap<usize, usize> = HashMap::new();
    let mut keep = vec![false; genes.len()];
    for &gi in &order {
        let g = &genes[gi];
        let cr = contig_rank[&g.contig];
        let tree = kept_iv.entry(cr).or_default();
        let ml = *max_len.get(&cr).unwrap_or(&0);
        let lo_begin = g.begin.saturating_sub(ml);
        let mut clash = false;
        for (&kb, ends) in tree.range(lo_begin..=g.end) {
            for &ke in ends {
                let lo = kb.max(g.begin) as i64;
                let hi = ke.min(g.end) as i64;
                if hi - lo + 1 > max_overlap {
                    clash = true;
                    break;
                }
            }
            if clash {
                break;
            }
        }
        if !clash {
            tree.entry(g.begin).or_default().push(g.end);
            let len = g.end - g.begin + 1;
            let e = max_len.entry(cr).or_insert(0);
            if len > *e {
                *e = len;
            }
            keep[gi] = true;
        }
    }

    let mut kept: Vec<Gene> = genes
        .into_iter()
        .zip(keep)
        .filter_map(|(g, k)| if k { Some(g) } else { None })
        .collect();

    // Restore genomic order (input contig order) and renumber per contig.
    kept.sort_by(|a, b| {
        contig_rank[&a.contig]
            .cmp(&contig_rank[&b.contig])
            .then(a.begin.cmp(&b.begin))
            .then(a.end.cmp(&b.end))
    });
    let mut counts: HashMap<String, usize> = HashMap::new();
    for g in &mut kept {
        let c = counts.entry(g.contig.clone()).or_insert(0);
        *c += 1;
        g.id = *c;
    }
    kept
}

/// Straightforward O(n²) reference implementation (also used for a negative
/// `max_overlap`).
fn resolve_overlaps_naive(mut genes: Vec<Gene>, max_overlap: i64) -> Vec<Gene> {
    let mut contig_rank: HashMap<String, usize> = HashMap::new();
    for g in &genes {
        let next = contig_rank.len();
        contig_rank.entry(g.contig.clone()).or_insert(next);
    }
    genes.sort_by(|a, b| {
        b.cscore
            .partial_cmp(&a.cscore)
            .unwrap_or(Ordering::Equal)
            .then((b.end - b.begin).cmp(&(a.end - a.begin)))
    });
    let mut kept: Vec<Gene> = Vec::with_capacity(genes.len());
    for g in genes {
        let clashes = kept.iter().any(|k| {
            k.contig == g.contig && {
                let lo = k.begin.max(g.begin) as i64;
                let hi = k.end.min(g.end) as i64;
                (hi - lo + 1).max(0) > max_overlap
            }
        });
        if !clashes {
            kept.push(g);
        }
    }
    kept.sort_by(|a, b| {
        contig_rank[&a.contig]
            .cmp(&contig_rank[&b.contig])
            .then(a.begin.cmp(&b.begin))
            .then(a.end.cmp(&b.end))
    });
    let mut counts: HashMap<String, usize> = HashMap::new();
    for g in &mut kept {
        let c = counts.entry(g.contig.clone()).or_insert(0);
        *c += 1;
        g.id = *c;
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene::Engine;

    fn g(contig: &str, begin: usize, end: usize, score: f64) -> Gene {
        Gene {
            contig: contig.into(),
            id: 0,
            begin,
            end,
            strand: 1,
            start_type: 0,
            rbs_motif: "None".into(),
            rbs_spacer: "None".into(),
            partial_left: false,
            partial_right: false,
            score,
            cscore: score,
            sscore: 0.0,
            conf: 0.0,
            gc: 0.0,
            trans_table: 11,
            engine: Engine::GeneMark,
            aa: String::new(),
            nt: String::new(),
        }
    }

    #[test]
    fn keeps_best_and_input_contig_order() {
        let genes = vec![
            g("contig2", 100, 400, 5.0),
            g("contig10", 1, 300, 9.0),
            g("contig2", 150, 600, 7.0), // overlaps the first by 251 > 60
            g("contig2", 380, 900, 1.0), // overlaps the 150..600 winner by 221
            g("contig10", 290, 800, 3.0), // overlaps 1..300 by 11 <= 60
        ];
        let out = resolve_overlaps(genes, 60);
        let v: Vec<(&str, usize, usize, usize)> = out
            .iter()
            .map(|x| (x.contig.as_str(), x.begin, x.end, x.id))
            .collect();
        assert_eq!(
            v,
            vec![
                ("contig2", 150, 600, 1),
                ("contig10", 1, 300, 1),
                ("contig10", 290, 800, 2)
            ]
        );
    }

    #[test]
    fn interval_search_matches_naive_reference() {
        // deterministic pseudo-random gene sets, several overlap limits
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for round in 0..200 {
            let n = 2 + (next() % 40) as usize;
            let genes: Vec<Gene> = (0..n)
                .map(|_| {
                    let b = 1 + (next() % 5000) as usize;
                    let len = 30 + (next() % 3000) as usize;
                    let c = if next() % 3 == 0 { "b" } else { "a" };
                    g(c, b, b + len, (next() % 1000) as f64 / 10.0)
                })
                .collect();
            for &mo in &[-1i64, 0, 1, 60, 300] {
                let fast = resolve_overlaps(genes.clone(), mo);
                let slow = resolve_overlaps_naive(genes.clone(), mo);
                let key = |v: &[Gene]| v.iter().map(|x| (x.contig.clone(), x.begin, x.end, x.id)).collect::<Vec<_>>();
                assert_eq!(key(&fast), key(&slow), "round {round} max_overlap {mo}");
            }
        }
    }

    #[test]
    fn long_kept_gene_still_blocks_later_candidates() {
        // a long kept gene starting far upstream must still be found
        let genes = vec![g("c", 1, 10_000, 10.0), g("c", 9_000, 9_500, 1.0)];
        let out = resolve_overlaps(genes, 60);
        assert_eq!(out.len(), 1);
    }
}
