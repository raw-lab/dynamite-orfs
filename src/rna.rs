//! Lightweight, dependency-free non-coding-RNA detection: a structural tRNA
//! scanner and a motif-anchored rRNA locator.
//!
//! **This is heuristic, not a covariance-model search.** The reference tools —
//! tRNAscan-SE (Lowe & Eddy 1997; Chan & Lowe 2019) for tRNAs and barrnap /
//! Infernal+Rfam for rRNAs — use trained covariance models and are far more
//! sensitive and specific. DYNAMITE's built-in scanner exists for quick,
//! offline first-pass screening and for the `doctor` self-test; for publication
//! annotation, run tRNAscan-SE and barrnap.
//!
//! tRNA method: slide canonical-length windows and test for a cloverleaf —
//! a 7-bp acceptor stem (Watson–Crick + G·U wobble) closing the molecule and a
//! 5-bp anticodon stem flanking a 7-nt loop — then read the anticodon and map
//! it to an amino acid via the genetic code. rRNA method: anchor on short
//! universally conserved SSU/LSU motifs and report candidate spans of the
//! expected length.

use crate::gencode::GeneticCode;
use crate::seq::{decode_base, Seq};

/// A predicted non-coding RNA feature (1-based inclusive).
#[derive(Clone, Debug)]
pub struct RnaFeature {
    pub kind: String, // "tRNA", "rRNA_16S", "rRNA_23S", "rRNA_5S", ...
    pub begin: usize,
    pub end: usize,
    pub strand: i8,
    pub score: f64,
    pub note: String,
}

#[inline]
fn pairs(a: u8, b: u8) -> bool {
    // Watson–Crick (b == complement) or G·U / U·G wobble. 2-bit A0 C1 G2 T3.
    b == 3 - a || (a == 2 && b == 3) || (a == 3 && b == 2)
}

fn stem_score(arr: &[u8], i5: usize, i3: usize, len: usize) -> usize {
    // Count pairs between arr[i5..i5+len] and arr[i3-len+1..=i3] read antiparallel.
    let mut s = 0;
    for k in 0..len {
        if i5 + k >= arr.len() || i3 < k {
            break;
        }
        if pairs(arr[i5 + k], arr[i3 - k]) {
            s += 1;
        }
    }
    s
}

/// Watson–Crick-only stem score (no wobble) — used for the acceptor stem so the
/// gate is stricter and random sequence rarely passes.
fn stem_score_wc(arr: &[u8], i5: usize, i3: usize, len: usize) -> usize {
    let mut s = 0;
    for k in 0..len {
        if i5 + k >= arr.len() || i3 < k {
            break;
        }
        if arr[i3 - k] == 3 - arr[i5 + k] {
            s += 1;
        }
    }
    s
}

/// Scan both strands for tRNA-like cloverleaves.
///
/// Gating (all required): a 7-bp acceptor stem (≥6/7 Watson–Crick), a 5-bp
/// anticodon stem (≥4/5) enclosing a 7-nt loop, a 5-bp T-arm stem (≥4/5)
/// enclosing a loop carrying the invariant TΨC (…T-T-C…), and a sensible
/// length. These conserved-structure requirements are what keep the
/// false-positive rate low on random sequence; sensitivity is still well below
/// tRNAscan-SE.
pub fn find_trnas(seq: &Seq, code: &GeneticCode) -> Vec<RnaFeature> {
    let mut out = Vec::new();
    let n = seq.len();
    for (strand, arr) in [(1i8, &seq.fwd), (-1i8, &seq.rev)] {
        for len in [76usize, 75, 77, 74] {
            if len > n {
                continue;
            }
            let mut i = 0;
            while i + len <= n {
                // Acceptor stem: 7 bp, all paired, at most one G·U wobble.
                // Try several 3' offsets so we catch genomic genes (no CCA) and
                // mature tRNAs (discriminator + CCA at the 3' end).
                let mut acc_wc = 0usize;
                let mut acc_all = 0usize;
                for d in 1..=4usize {
                    if i + len < d + 1 {
                        continue;
                    }
                    let end3 = i + len - 1 - d;
                    let w = stem_score_wc(arr, i, end3, 7);
                    let a = stem_score(arr, i, end3, 7);
                    if a > acc_all || (a == acc_all && w > acc_wc) {
                        acc_all = a;
                        acc_wc = w;
                    }
                }
                if acc_wc < 6 || acc_all < 7 {
                    i += 1;
                    continue;
                }
                // Anticodon arm (5/5) in the central region.
                let ac = find_arm(arr, i + len * 28 / 100, i + len * 52 / 100, i + len);
                // T-arm (5/5) in the 3' region; loop carries the invariant T-T-C.
                let tarm = find_arm(arr, i + len * 60 / 100, i + len * 88 / 100, i + len);
                if let (Some((al, acs)), Some((tl, tss))) = (ac, tarm) {
                    if acs < 5 || tss < 5 {
                        i += 1;
                        continue;
                    }
                    // Conserved TΨC: loop[1]=U(T) and loop[2]=C (the "TC" of TΨC).
                    let tloop_ok = tl + 3 < arr.len() && arr[tl + 1] == 3 && arr[tl + 2] == 1;
                    if !tloop_ok {
                        i += 1;
                        continue;
                    }
                    let anticodon = [arr[al + 2], arr[al + 3], arr[al + 4]];
                    // Codon recognised = reverse complement of the anticodon.
                    let codon = [3 - anticodon[2], 3 - anticodon[1], 3 - anticodon[0]];
                    let aa = code.translate(codon[0], codon[1], codon[2]) as char;
                    let score = acc_all as f64 + acs as f64 + tss as f64;
                    if aa != '*' {
                        let (begin, end) = to_global(strand, n, i, len);
                        let ac_txt: String = anticodon
                            .iter()
                            .map(|&b| decode_base(b).to_ascii_uppercase() as char)
                            .collect();
                        out.push(RnaFeature {
                            kind: "tRNA".into(),
                            begin,
                            end,
                            strand,
                            score,
                            note: format!("anticodon={ac_txt};aa={aa};model=heuristic"),
                        });
                        i += len; // skip past this hit
                        continue;
                    }
                }
                i += 1;
            }
        }
    }
    dedup_overlaps(&mut out);
    out
}

/// Find a 5-bp stem enclosing a 7-nt loop within [lo, hi]. Returns
/// (loop_start_index, stem_score) for the best stem found (≥3 pairs).
fn find_arm(arr: &[u8], lo: usize, hi: usize, end: usize) -> Option<(usize, usize)> {
    let mut best: Option<(usize, usize)> = None;
    let mut al = lo.max(5);
    while al + 12 <= end && al <= hi {
        let s = stem_score(arr, al - 5, al + 11, 5);
        if s >= 3 && best.map(|(_, bs)| s > bs).unwrap_or(true) {
            best = Some((al, s));
        }
        al += 1;
    }
    best
}

const SSU_MOTIFS: &[&str] = &[
    "gtgccagcagccgcggtaa", // ~515F universal SSU
    "ggattagataccc",        // ~787 region
    "aaactcaaaggaattgacgg", // ~907 region
];
const LSU_MOTIFS: &[&str] = &[
    "agtcgggacctaaggcgag",
    "ccgtaacttcgggataagg",
];

/// Heuristic rRNA locator: anchor on conserved SSU/LSU motifs and report a
/// candidate span of the expected length.
pub fn find_rrnas(seq: &Seq) -> Vec<RnaFeature> {
    let mut out = Vec::new();
    let n = seq.len();
    // Ambiguous bases stay `N` so they can never complete a motif.
    let fwd: Vec<u8> = seq.to_ascii_upper();
    let rev: Vec<u8> = fwd
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect();
    for (strand, hay) in [(1i8, &fwd), (-1i8, &rev)] {
        let s = String::from_utf8_lossy(hay).to_lowercase();
        for (kind, motifs, span) in [
            ("rRNA_SSU(16S/18S)", SSU_MOTIFS, 1500usize),
            ("rRNA_LSU(23S/28S)", LSU_MOTIFS, 2900usize),
        ] {
            for m in motifs {
                if let Some(pos) = s.find(m) {
                    // Centre an expected-length window on the motif.
                    let start = pos.saturating_sub(span / 3);
                    let end = (start + span).min(n);
                    let (begin, gend) = to_global(strand, n, start, end - start);
                    out.push(RnaFeature {
                        kind: kind.into(),
                        begin,
                        end: gend,
                        strand,
                        score: 1.0,
                        note: format!("anchor_motif={m};model=heuristic;verify_with=barrnap/Infernal"),
                    });
                    break; // one hit per class per strand is enough as a flag
                }
            }
        }
    }
    out
}

/// All ncRNA features (tRNA + rRNA).
pub fn find_all(seq: &Seq, code: &GeneticCode) -> Vec<RnaFeature> {
    let mut v = find_trnas(seq, code);
    v.extend(find_rrnas(seq));
    v.sort_by_key(|f| (f.begin, f.end));
    v
}

fn to_global(strand: i8, n: usize, i0: usize, len: usize) -> (usize, usize) {
    if strand == 1 {
        (i0 + 1, i0 + len)
    } else {
        // window [i0, i0+len) on the reverse array → global coordinates
        let g_begin = n - (i0 + len) + 1;
        let g_end = n - i0;
        (g_begin, g_end)
    }
}

fn dedup_overlaps(feats: &mut Vec<RnaFeature>) {
    feats.sort_by(|a, b| {
        a.begin
            .cmp(&b.begin)
            .then(b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut kept: Vec<RnaFeature> = Vec::new();
    for f in feats.drain(..) {
        if let Some(last) = kept.last() {
            if f.begin <= last.end && f.strand == last.strand && f.kind == last.kind {
                continue; // overlapping same-kind call already kept (higher score)
            }
        }
        kept.push(f);
    }
    *feats = kept;
}

/// Render features as a complete GFF3 document (header + records).
pub fn to_gff(contig: &str, feats: &[RnaFeature]) -> String {
    let mut out = String::from("##gff-version 3\n");
    out.push_str(&to_gff_records(contig, feats));
    out
}

/// Render features as GFF3 records only (no `##gff-version` header), for
/// concatenating several contigs into one valid GFF3 file.
pub fn to_gff_records(contig: &str, feats: &[RnaFeature]) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    for (i, f) in feats.iter().enumerate() {
        let _ = writeln!(
            out,
            "{}\tDYNAMITE\t{}\t{}\t{}\t{:.1}\t{}\t.\tID={}_{}_{};{}",
            contig,
            if f.kind == "tRNA" { "tRNA" } else { "rRNA" },
            f.begin,
            f.end,
            f.score,
            if f.strand == 1 { '+' } else { '-' },
            contig,
            f.kind,
            i + 1,
            f.note,
        );
    }
    out
}
