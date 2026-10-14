//! **FragGeneScan** engine — a thin adapter over the vendored, verbatim
//! [`crate::fgs`] core of FragGeneScanRs (Van der Jeugt et al. 2022, the
//! maintained Rust port of FragGeneScan; Rho et al. 2010).
//!
//! Unlike the other engines this is **not** a from-scratch reimplementation: it
//! drives the real FragGeneScan forward/Viterbi HMM with its original trained
//! models, then maps each predicted gene into DYNAMITE's [`Gene`] using
//! FragGeneScan's own nucleotide/protein extraction. Output is therefore
//! faithful to FragGeneScanRs run with the same model.
//!
//! The model is always FragGeneScan's `complete` one (FragGeneScanRs
//! `-t complete`), whose insertion/deletion transitions are 1e-18: genes are
//! not corrected for sequencing-error frameshifts. The error models
//! (`illumina_*`, `454_*`, `sanger_*`) are embedded but not selected.
//!
//! `whole_genome = false` runs the short-read / fragment decode (genes may run
//! off either end) — used for `reads` and `long-reads` modes, which therefore
//! give identical calls. `whole_genome = true` runs the complete-sequence
//! decode — used when FragGeneScan is forced on assembled contigs.

use crate::fgs::dna::Nuc;
use crate::fgs::hmm;
use crate::fgs::viterbi::viterbi;
use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::seq::Seq;
use std::path::PathBuf;

/// Pull the sequence line out of FragGeneScan's `>header\nSEQUENCE\n` buffer.
fn seq_line(buf: &[u8]) -> String {
    let s = String::from_utf8_lossy(buf);
    s.lines().nth(1).unwrap_or("").to_string()
}

fn gc_of(nt: &str) -> f64 {
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

/// Run the FragGeneScan HMM over every contig/read. `whole_genome` selects the
/// complete-sequence (`true`) vs short-read/fragment (`false`) decode.
///
/// Records are decoded in parallel on the current rayon pool (so `--threads`
/// applies). FragGeneScan reports no partial flags of its own; they are
/// inferred from the genome codons at each gene's ends (no start codon → 5'
/// partial, no stop codon → 3' partial — read-mode genes are frequently
/// fragments).
///
/// Coordinates are FragGeneScanRs's own. For genes **without** frameshifts the
/// nucleotide / protein strings are taken from the genome at those
/// coordinates: upstream's minus-strand start refinement moves the coordinate
/// but appends only the new start codon, so its strings silently drop the
/// codons in between (and read mode omits the stop codon). Frameshift-corrected
/// genes keep FragGeneScan's own (indel-aware) strings.
pub fn call(contigs: &[(String, Seq)], whole_genome: bool) -> (Vec<Gene>, u32) {
    call_closed(contigs, whole_genome, false)
}

/// [`call`] honouring `closed`: with `closed = true` a start codon at the very
/// edge of a sequence counts as a real start (no open-edge 5'-partial rule).
pub fn call_closed(contigs: &[(String, Seq)], whole_genome: bool, closed: bool) -> (Vec<Gene>, u32) {
    use rayon::prelude::*;
    let tt = 11u32; // FragGeneScan's models assume the standard/bacterial code
    let code = match GeneticCode::new(tt) {
        Some(c) => c,
        None => return (Vec::new(), tt),
    };
    // A path that cannot exist forces the embedded (`include_bytes!`) models.
    let (global, locals) = match hmm::get_train_from_file(PathBuf::from("/nonexistent/dynamite-fgs-models"), PathBuf::from("complete")) {
        Ok(m) => m,
        Err(_) => return (Vec::new(), tt),
    };

    let genes: Vec<Gene> = contigs
        .par_iter()
        .flat_map_iter(|(name, s)| {
            let local = call_one(&global, &locals, name, s, whole_genome, &code, tt);
            let mut out = Vec::with_capacity(local.len());
            for (mut g, frameshifted) in local {
                let src = if frameshifted {
                    crate::finalize::SeqSource::Engine
                } else {
                    crate::finalize::SeqSource::Genome
                };
                if crate::finalize::finalize_gene(&mut g, s, &code, closed, src, crate::finalize::StartRule::Bacterial) {
                    if !frameshifted {
                        g.gc = gc_of(&g.nt);
                    }
                    out.push(g);
                }
            }
            out.into_iter()
        })
        .collect();
    (genes, tt)
}

fn call_one(
    global: &hmm::Global,
    locals: &Vec<hmm::Local>,
    name: &str,
    s: &Seq,
    whole_genome: bool,
    code: &GeneticCode,
    tt: u32,
) -> Vec<(Gene, bool)> {
    let mut genes = Vec::new();
    // FragGeneScan's decoder indexes the first codon unconditionally and
    // aborts on sequences shorter than 3 bp; it cannot emit a gene shorter
    // than 60 bp anyway, so such sequences are skipped (output unchanged).
    if s.len() < 3 {
        return genes;
    }
    let dummy = b"s".to_vec();
    let ascii = s.to_ascii_upper();
    let nseq: Vec<Nuc> = ascii.iter().map(|&b| Nuc::from(b)).collect();
    let head: Vec<u8> = name.bytes().collect();
    let pred = viterbi(global, locals, head, nseq, whole_genome);

    for (i, fg) in pred.genes.iter().enumerate() {
        // FragGeneScan's own (insertion/deletion-aware) sequence extraction.
        let mut ntbuf = Vec::new();
        let _ = fg.dna(&mut ntbuf, &dummy, false);
        let mut aabuf = Vec::new();
        let _ = fg.protein(&mut aabuf, &dummy, whole_genome);
        let nt = seq_line(&ntbuf).to_ascii_lowercase();
        let aa = seq_line(&aabuf);
        if aa.is_empty() {
            continue;
        }

        let strand = if fg.forward_strand { 1i8 } else { -1 };
        // Classify the start codon from the coding-strand 5' triplet (the
        // finaliser re-checks it against the genome and the partial flags).
        let b = nt.as_bytes();
        let start_type = if b.len() >= 3 {
            let (c1, a1) = crate::seq::encode_base(b[0]);
            let (c2, a2) = crate::seq::encode_base(b[1]);
            let (c3, a3) = crate::seq::encode_base(b[2]);
            if a1 || a2 || a3 || !code.is_start(c1, c2, c3) {
                3
            } else {
                GeneticCode::start_type(c1, c2, c3)
            }
        } else {
            3
        };

        let frameshifted = !fg.inserted.is_empty() || !fg.deleted.is_empty();
        genes.push((Gene {
            contig: name.to_string(),
            id: i + 1,
            begin: fg.start,
            end: fg.end,
            strand,
            start_type,
            rbs_motif: "None".into(),
            rbs_spacer: "None".into(),
            partial_left: false,
            partial_right: false,
            score: fg.score,
            cscore: fg.score,
            sscore: 0.0,
            conf: 100.0 / (1.0 + (-fg.score / 20.0).exp()),
            gc: gc_of(&nt),
            trans_table: tt,
            engine: Engine::FragGeneScan,
            aa,
            nt,
        }, frameshifted));
    }
    genes
}
