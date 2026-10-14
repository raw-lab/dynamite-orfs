//! **GeneMark** — a self-training, 3-periodic inhomogeneous Markov gene finder
//! in the spirit of GeneMark / GeneMark.hmm (Borodovsky & McIninch 1993;
//! Besemer et al. 2001).
//!
//! Algorithm: enumerate start→stop ORFs in all six frames, bootstrap a coding
//! model from the long ORFs (presumed coding), then score every ORF with a
//! coding/non-coding log-likelihood ratio and keep those that look coding.
//!
//! Honest scope: this captures GeneMark's *statistical core* (codon-position
//! Markov chains, self-training, LLR scoring) but is **not** the full
//! GeneMark.hmm — there is no hidden-state Viterbi over the whole sequence with
//! duration/RBS submodels, and no GeneMarkS iterative motif refinement.

use crate::dedup::resolve_overlaps;
use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::markov::MarkovModel;
use crate::orf::{find_gene_orfs, OrfCand};
use crate::seq::Seq;

const MIN_ORF: usize = 90;
const MAX_OVERLAP: i64 = 60; // tolerated overlap between kept genes (nt)
const TRAIN_MIN: usize = 300; // long ORFs used to bootstrap the coding model

/// Call genes on a set of contigs with GeneMark-style self-training.
///
/// The coding model is trained once on the pooled long ORFs across all contigs,
/// then applied to every contig. `code` defaults to 11 if `None`.
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    let tt = fixed.unwrap_or(11);
    let code = match GeneticCode::new(tt) {
        Some(c) => c,
        None => return (Vec::new(), tt),
    };
    let require_start = closed; // closed → starts only; open → keep edge ORFs too

    // 1) Enumerate ORFs per contig; remember which contig each training ORF
    //    came from (coordinates are contig-local).
    let mut per_contig: Vec<Vec<OrfCand>> = Vec::with_capacity(contigs.len());
    let mut training: Vec<Vec<(usize, usize, i8)>> = Vec::with_capacity(contigs.len());
    for (_, s) in contigs {
        let cands = find_gene_orfs(s, &code, MIN_ORF, require_start, !require_start);
        let local: Vec<(usize, usize, i8)> = cands
            .iter()
            .filter(|c| c.end - c.begin + 1 >= TRAIN_MIN && !c.partial_left && !c.partial_right)
            .map(|c| (c.begin, c.end, c.strand))
            .collect();
        training.push(local);
        per_contig.push(cands);
    }

    // 2) Train the coding model. Fall back to a lower length bar if too sparse.
    if training.iter().map(|t| t.len()).sum::<usize>() < 10 {
        for (t, cands) in training.iter_mut().zip(&per_contig) {
            *t = cands
                .iter()
                .filter(|c| c.end - c.begin + 1 >= MIN_ORF + 60)
                .map(|c| (c.begin, c.end, c.strand))
                .collect();
        }
    }

    // Train one model on the long ORFs pooled across every contig.
    let sets: Vec<crate::markov::TrainingSet<'_>> = contigs
        .iter()
        .zip(&training)
        .map(|((_, s), t)| (s, t.as_slice()))
        .collect();
    let model = MarkovModel::train_pooled(&sets);

    // 3) Score every ORF; keep coding ones (LLR > 0). IDs are per-contig.
    let mut genes = Vec::new();
    for ((name, s), cands) in contigs.iter().zip(&per_contig) {
        let mut local = Vec::new();
        for c in cands {
            let llr = model.llr(s, c.begin, c.end, c.strand);
            if llr <= 0.0 {
                continue;
            }
            let mut g = build_gene(s, &code, c, llr, local.len() + 1);
            g.contig = name.clone();
            local.push(g);
        }
        genes.extend(crate::finalize::finalize_genes(local, s, &code, closed, crate::finalize::SeqSource::Genome, crate::finalize::StartRule::Table));
    }
    (resolve_overlaps(genes, MAX_OVERLAP), tt)
}

fn build_gene(seq: &Seq, code: &GeneticCode, c: &OrfCand, llr: f64, id: usize) -> Gene {
    let mut g = Gene {
        contig: String::new(),
        id,
        begin: c.begin,
        end: c.end,
        strand: c.strand,
        start_type: c.start_type,
        rbs_motif: "None".into(),
        rbs_spacer: "None".into(),
        partial_left: c.partial_left,
        partial_right: c.partial_right,
        score: llr,
        cscore: llr,
        sscore: 0.0,
        conf: 100.0 / (1.0 + (-llr / 20.0).exp()),
        gc: 0.0,
        trans_table: code.id,
        engine: Engine::GeneMark,
        aa: String::new(),
        nt: String::new(),
    };
    g.fill_sequences(seq, code);
    g.gc = gc_of(&g.nt);
    g
}

fn gc_of(nt: &str) -> f64 {
    // GC over unambiguous bases only (N is neither GC nor AT).
    let acgt = nt
        .bytes()
        .filter(|b| matches!(b, b'a' | b'A' | b'c' | b'C' | b'g' | b'G' | b't' | b'T'))
        .count();
    if acgt == 0 {
        return 0.0;
    }
    let gc = nt
        .bytes()
        .filter(|b| matches!(b, b'g' | b'G' | b'c' | b'C'))
        .count();
    gc as f64 / acgt as f64
}
