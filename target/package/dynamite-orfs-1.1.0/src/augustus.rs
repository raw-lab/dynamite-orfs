//! **AUGUSTUS** — a eukaryotic ORF/CDS caller in the spirit of AUGUSTUS
//! (Stanke & Waack 2003).
//!
//! Eukaryotic genes default to the standard code (table 1) and almost always
//! initiate at ATG within a Kozak context. This module finds coding ORFs with
//! the shared 3-periodic Markov model, prefers ATG starts in a strong Kozak
//! context (purine at −3, G at +4), and reports single-exon CDS.
//!
//! Honest scope: this is a **single-exon** eukaryotic predictor. It does *not*
//! implement AUGUSTUS's full generalized HMM with intron/exon states, trained
//! splice-site (GT–AG) and branch-point models, or length distributions, so
//! multi-exon genes are reported by their individual coding exons rather than
//! stitched transcripts. For production eukaryotic annotation use AUGUSTUS,
//! BRAKER, or Helixer.

use crate::dedup::resolve_overlaps;
use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::markov::MarkovModel;
use crate::orf::find_gene_orfs;
use crate::seq::Seq;

const MIN_ORF: usize = 150; // eukaryotic single-exon CDS bar (50 codons)
const MAX_OVERLAP: i64 = 60; // tolerated overlap between kept genes (nt)
const TRAIN_MIN: usize = 300;

/// Call eukaryotic CDS. `code` defaults to 1 (standard) if `None`.
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    let tt = fixed.unwrap_or(1);
    let code = match GeneticCode::new(tt) {
        Some(c) => c,
        None => return (Vec::new(), tt),
    };

    // Enumerate ATG-initiated ORFs (start required); keep edge ORFs only if open.
    // Training ORFs stay paired with the contig they were found on.
    let mut per_contig = Vec::with_capacity(contigs.len());
    let mut training: Vec<Vec<(usize, usize, i8)>> = Vec::with_capacity(contigs.len());
    for (_, s) in contigs {
        let cands = find_gene_orfs(s, &code, MIN_ORF, true, !closed);
        training.push(
            cands
                .iter()
                .filter(|c| c.end - c.begin + 1 >= TRAIN_MIN && !c.partial_left && !c.partial_right)
                .map(|c| (c.begin, c.end, c.strand))
                .collect(),
        );
        per_contig.push(cands);
    }
    // One model pooled over every contig (it used to train on contig 0 with
    // the coordinates of *all* contigs: a wrong model, and an index-out-of-
    // bounds abort whenever a later contig was longer than the first).
    let sets: Vec<crate::markov::TrainingSet<'_>> = contigs
        .iter()
        .zip(&training)
        .map(|((_, s), t)| (s, t.as_slice()))
        .collect();
    let model = MarkovModel::train_pooled(&sets);

    let mut genes = Vec::new();
    for ((name, s), cands) in contigs.iter().zip(&per_contig) {
        let mut id = 0usize;
        let mut local = Vec::new();
        for c in cands {
            let mut llr = model.llr(s, c.begin, c.end, c.strand);
            // Kozak preference: small bonus for ATG with a strong context.
            if c.start_type == 0 {
                llr += kozak_bonus(s, c.begin, c.end, c.strand);
            } else {
                // non-ATG eukaryotic starts are rare — penalise lightly
                llr -= 2.0;
            }
            if llr <= 0.0 {
                continue;
            }
            id += 1;
            let mut g = Gene {
                contig: name.clone(),
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
                conf: 100.0 / (1.0 + (-llr / 25.0).exp()),
                gc: 0.0,
                trans_table: code.id,
                engine: Engine::Augustus,
                aa: String::new(),
                nt: String::new(),
            };
            g.fill_sequences(s, &code);
            g.gc = gc_of(&g.nt);
            local.push(g);
        }
        genes.extend(crate::finalize::finalize_genes(local, s, &code, closed, crate::finalize::SeqSource::Genome, crate::finalize::StartRule::Table));
    }
    (resolve_overlaps(genes, MAX_OVERLAP), tt)
}

/// Kozak context bonus: +1.5 for a purine (A/G) at −3, +1.5 for G at +4.
fn kozak_bonus(seq: &Seq, begin: usize, end: usize, strand: i8) -> f64 {
    let n = seq.len();
    let mut bonus = 0.0;
    // Work on the coding strand 5'->3'. Start codon occupies the first 3 bases.
    let (arr, start0): (&[u8], i64) = if strand == 1 {
        (&seq.fwd, begin as i64 - 1)
    } else {
        (&seq.rev, (n - end) as i64)
    };
    let at = |p: i64| -> Option<u8> {
        if p >= 0 && (p as usize) < arr.len() {
            Some(arr[p as usize])
        } else {
            None
        }
    };
    if let Some(b) = at(start0 - 3) {
        if b == 0 || b == 2 {
            bonus += 1.5; // A or G (purine) at -3
        }
    }
    if let Some(b) = at(start0 + 3) {
        if b == 2 {
            bonus += 1.5; // G at +4
        }
    }
    bonus
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
