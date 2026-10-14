//! Plain six-frame translation and six-frame ORF enumeration.
//!
//! Unlike the statistical callers this module makes no model-based decisions:
//! it simply translates all three forward and three reverse frames, and (on
//! request) lists every ORF between stop codons. Useful for quick inspection,
//! primer design, and sanity checks.

use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::orf::find_orfs;
use crate::seq::Seq;

/// One translated reading frame.
pub struct Frame {
    /// `+1,+2,+3,-1,-2,-3`.
    pub frame: i8,
    pub protein: String,
}

/// Translate all six frames of a contig (stops shown as `*`, codons with an
/// ambiguous base as `X`).
pub fn translate_six(seq: &Seq, code: &GeneticCode) -> Vec<Frame> {
    let n = seq.len();
    let mut frames = Vec::with_capacity(6);
    for (strand, arr) in [(1i8, &seq.fwd), (-1i8, &seq.rev)] {
        for f in 0..3usize {
            let mut protein = String::new();
            let mut i = f;
            while i + 3 <= n {
                let amb = if strand > 0 { seq.codon_has_amb_fwd(i) } else { seq.codon_has_amb_rev(i) };
                let aa = if amb { b'X' } else { code.translate(arr[i], arr[i + 1], arr[i + 2]) };
                protein.push(aa as char);
                i += 3;
            }
            frames.push(Frame {
                frame: strand * (f as i8 + 1),
                protein,
            });
        }
    }
    frames
}

/// Render the six frames as a protein FASTA string.
pub fn to_fasta(contig: &str, frames: &[Frame]) -> String {
    let mut out = String::new();
    for fr in frames {
        out.push('>');
        out.push_str(contig);
        out.push_str(&format!("_frame{}{}\n", if fr.frame > 0 { "+" } else { "" }, fr.frame));
        // wrap at 60
        let bytes = fr.protein.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let end = (i + 60).min(bytes.len());
            out.push_str(std::str::from_utf8(&bytes[i..end]).unwrap());
            out.push('\n');
            i = end;
        }
        if fr.protein.is_empty() {
            out.push('\n');
        }
    }
    out
}

/// Enumerate six-frame ORFs as [`Gene`] records.
///
/// `require_start` keeps only ORFs that open with a start codon; otherwise
/// contig-edge fragments are reported too (interior stop-to-stop regions
/// without a start codon are never ORFs). The contig is split at assembly
/// gaps first, so no ORF spans a run of N; ORFs that run into a gap are
/// reported partial at that end, with gap codons trimmed. With
/// `require_start` the contig ends are treated as closed (a start codon at the
/// very edge counts as a real start).
pub fn orfs(
    contig: &str,
    seq: &Seq,
    code: &GeneticCode,
    min_len: usize,
    require_start: bool,
) -> Vec<Gene> {
    let mut genes = Vec::new();
    for (off, len) in crate::finalize::segments_of_seq(seq) {
        let sub;
        let s: &Seq = if off == 0 && len == seq.len() {
            seq
        } else {
            sub = seq.subseq(off, len);
            &sub
        };
        let cands = find_orfs(s, code, min_len, require_start, !require_start);
        let mut local = Vec::with_capacity(cands.len());
        for c in cands.into_iter() {
            let mut g = Gene {
                contig: contig.to_string(),
                id: 0,
                begin: c.begin,
                end: c.end,
                strand: c.strand,
                start_type: c.start_type,
                rbs_motif: "None".into(),
                rbs_spacer: "None".into(),
                partial_left: c.partial_left,
                partial_right: c.partial_right,
                score: 0.0,
                cscore: 0.0,
                sscore: 0.0,
                conf: 0.0,
                gc: 0.0,
                trans_table: code.id,
                engine: Engine::SixFrame,
                aa: String::new(),
                nt: String::new(),
            };
            g.fill_sequences(s, code);
            local.push(g);
        }
        let mut local = crate::finalize::finalize_genes(local, s, code, require_start, crate::finalize::SeqSource::Genome, crate::finalize::StartRule::Table);
        crate::finalize::shift_genes(&mut local, off, contig);
        genes.extend(local);
    }
    // Same order as before: forward strand first, then reverse; ids 1..n.
    for (i, g) in genes.iter_mut().enumerate() {
        g.id = i + 1;
    }
    genes
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gencode::GeneticCode;
    use crate::seq::Seq;

    // Regression test for the minus-strand partial-flag / stop-strip bug:
    // no protein may carry an internal or trailing stop codon, and complete
    // genes (start present) must begin with M — on either strand.
    #[test]
    fn proteins_have_no_stray_stops() {
        let dna = b"ATGAAACGCATTGCTGATGCAGGTTCACCGCTGATTGAAGCTGGTGCTGAAGCGTTTAAAGCATAA\
                    GGGCCCTTAATGCATGCATGCATGCATGCATGCATGCATGCATAATTTGGGCCCATGAAATTTTAA";
        let s = Seq::from_ascii(dna);
        let code = GeneticCode::new(11).unwrap();
        let genes = orfs("t", &s, &code, 30, false);
        assert!(!genes.is_empty());
        let mut saw_minus = false;
        for g in &genes {
            assert!(!g.aa.ends_with('*'), "trailing stop {}..{}", g.begin, g.end);
            assert!(!g.aa.trim_end_matches('*').contains('*'), "internal stop {}..{}", g.begin, g.end);
            if !g.partial_5p() {
                assert!(g.aa.starts_with('M'), "complete gene not M-started {}..{}", g.begin, g.end);
            }
            if g.strand == -1 {
                saw_minus = true;
            }
        }
        assert!(saw_minus, "test sequence should contain a minus-strand ORF");
    }
}
