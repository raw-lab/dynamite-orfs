//! The **PHANOTATE** engine: a faithful Rust port (in spirit and arithmetic)
//! of McNair *et al.* 2019's phage gene caller.
//!
//! PHANOTATE treats a phage genome as a weighted graph over the six reading
//! frames. Open reading frames are *rewarding* edges (negative weight); the
//! gaps and overlaps between them are *penalising* edges (positive weight).
//! The optimal gene set is the shortest source -> target path through that
//! graph.
//!
//! ## Differences from the reference implementation
//!
//! * **Arithmetic.** The reference uses Python `Decimal` (arbitrary precision)
//!   to keep the per-ORF `hold` product from underflowing. We use `f64` with a
//!   floor of `1e-300` on `hold` (so `1/hold` stays finite). For ordinary genes
//!   this is indistinguishable; only pathologically long ORFs *saturate* their
//!   reward instead of growing without bound — they are still selected.
//! * **Shortest path.** The reference shells out to `fastpathz`
//!   (Bellman-Ford). We run an in-process SPFA (queue-based Bellman-Ford-Moore)
//!   with a relaxation guard. The reference multiplies every edge by 1000
//!   before handing it to `fastpathz`; that is a uniform positive scaling and
//!   does not change which path is shortest, so we omit it.
//! * **tRNA masking.** The reference optionally calls `aragorn` / `tRNAscan-SE`.
//!   DYNAMITE never shells out, so tRNA edges are skipped (documented in the
//!   README). Gene calls are unaffected aside from the rare case where a tRNA
//!   would have masked a spurious ORF.
//! * **Stop codons** are taken from the selected NCBI genetic code rather than
//!   hard-coded to `taa,tga,tag`, so alternative-code genomes (e.g. code 15
//!   reading `tga` through as W) get correct ORF boundaries. Start-codon
//!   weights remain PHANOTATE's defaults (`atg`/`gtg`/`ttg`).

use std::collections::{HashMap, HashSet, VecDeque};

use crate::gcframe::{max_idx, min_idx, GcFrame};
use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::rbs::phanotate_score_rbs;
use crate::seq::Seq;

/// Default PHANOTATE minimum ORF length (nucleotides, includes the stop codon).
pub const MIN_ORF_LEN: i64 = 90;

const HOLD_FLOOR: f64 = 1e-300;
const S_CAP: f64 = 1e305;

// ---------------------------------------------------------------------------
// Sequence helpers (operate on lower-case ASCII, matching the reference)
// ---------------------------------------------------------------------------

#[inline]
fn rc_base(b: u8) -> u8 {
    match b {
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        other => other, // n, s, w map to themselves
    }
}

fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| rc_base(b)).collect()
}

// ---------------------------------------------------------------------------
// Edge scoring (faithful to functions.py)
// ---------------------------------------------------------------------------

/// `score_overlap`: penalty for two genes overlapping by `length` bp.
fn score_overlap(length: i64, diff: bool, pstop: f64) -> f64 {
    let o = 1.0 - pstop;
    let mut score = o.powf(length as f64);
    score = 1.0 / score;
    if diff {
        score += 1.0 / 0.05;
    }
    score.min(S_CAP)
}

/// `score_gap`: penalty for a `length` bp gap between genes.
fn score_gap(length: i64, diff: bool, pgap: f64) -> f64 {
    let g = 1.0 - pgap;
    if length > 300 {
        return g.powf(100.0) + length as f64;
    }
    let mut score = g.powf(length as f64 / 3.0);
    score = 1.0 / score;
    if diff {
        score += 1.0 / 0.05;
    }
    score.min(S_CAP)
}

// ---------------------------------------------------------------------------
// ORF store
// ---------------------------------------------------------------------------

struct Orf {
    start: i64,
    stop: i64,
    frame: i32,
    rbs_score: usize,
    seq: Vec<u8>,
    pstop: f64,
    hold: f64,
    weight_rbs: f64,
    weight: f64,
}

impl Orf {
    fn start_codon(&self) -> &[u8] {
        if self.seq.len() >= 3 {
            &self.seq[0..3]
        } else {
            &self.seq[..]
        }
    }
    fn stop_codon(&self) -> &[u8] {
        let n = self.seq.len();
        if n >= 3 {
            &self.seq[n - 3..]
        } else {
            &self.seq[..]
        }
    }
}

/// Per-ORF stop probability from its own base composition (`p_stop`).
fn p_stop(seq: &[u8]) -> f64 {
    let (mut a, mut t, mut g, mut c) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for &b in seq {
        match b {
            b'a' => a += 1.0,
            b't' => t += 1.0,
            b'g' => g += 1.0,
            b'c' => c += 1.0,
            _ => {}
        }
    }
    let length = seq.len() as f64;
    if length == 0.0 {
        return 0.0;
    }
    let (pa, pt, pg, _pc) = (a / length, t / length, g / length, c / length);
    pt * pa * pa + pt * pg * pa + pt * pa * pg
}

struct OrfStore {
    all: Vec<Orf>,
    by_stop: HashMap<i64, Vec<usize>>,
    other_end: HashMap<i64, i64>,
    orf_at: HashMap<(i64, i64), usize>,
    min_orf_len: i64,
}

impl OrfStore {
    fn new(min_orf_len: i64) -> Self {
        OrfStore {
            all: Vec::new(),
            by_stop: HashMap::new(),
            other_end: HashMap::new(),
            orf_at: HashMap::new(),
            min_orf_len,
        }
    }

    /// Mirror of `Orfs.add_orf`.
    fn add_orf(&mut self, start: i64, stop: i64, frame: i32, seq: Vec<u8>, rbs_score: usize) {
        if self.orf_at.contains_key(&(start, stop)) {
            return; // duplicate; reference would raise, we ignore
        }
        let first_for_stop = !self.by_stop.contains_key(&stop);
        let pstop = p_stop(&seq);
        let idx = self.all.len();
        self.all.push(Orf {
            start,
            stop,
            frame,
            rbs_score,
            seq,
            pstop,
            hold: 1.0,
            weight_rbs: 1.0,
            weight: 1.0,
        });
        self.by_stop.entry(stop).or_default().push(idx);
        self.orf_at.insert((start, stop), idx);

        if first_for_stop {
            self.other_end.insert(stop, start);
            self.other_end.insert(start, stop);
        } else {
            self.other_end.insert(start, stop);
            let cur = *self.other_end.get(&stop).unwrap_or(&start);
            if frame > 0 && start < cur {
                self.other_end.insert(stop, start);
            } else if frame < 0 && start > cur {
                self.other_end.insert(stop, start);
            }
        }
    }

    fn is_stop(&self, pos: i64) -> bool {
        self.by_stop.contains_key(&pos)
    }

    fn get_orf(&self, start: i64, stop: i64) -> Option<&Orf> {
        self.orf_at.get(&(start, stop)).map(|&i| &self.all[i])
    }

    /// `iter_in`: for each stop, ORFs ordered by start (asc if +frame, desc if
    /// -frame). Returns a vector of (stop, ordered indices).
    fn iter_in(&self) -> Vec<Vec<usize>> {
        let mut out = Vec::with_capacity(self.by_stop.len());
        for (_stop, idxs) in self.by_stop.iter() {
            let mut v = idxs.clone();
            let plus = self.all[v[0]].frame > 0;
            v.sort_by(|&x, &y| {
                let (sx, sy) = (self.all[x].start, self.all[y].start);
                if plus {
                    sx.cmp(&sy)
                } else {
                    sy.cmp(&sx)
                }
            });
            out.push(v);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Graph
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum NType {
    Start,
    Stop,
    Source,
    Target,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct PNode {
    typ: NType,
    frame: i32,
    position: i64,
}

struct Graph {
    nodes: Vec<PNode>,
    index: HashMap<PNode, usize>,
    // adjacency: out edges (to, weight)
    adj: Vec<Vec<(usize, f64)>>,
}

impl Graph {
    fn new() -> Self {
        Graph {
            nodes: Vec::new(),
            index: HashMap::new(),
            adj: Vec::new(),
        }
    }

    fn node(&mut self, n: PNode) -> usize {
        if let Some(&i) = self.index.get(&n) {
            return i;
        }
        let i = self.nodes.len();
        self.nodes.push(n);
        self.index.insert(n, i);
        self.adj.push(Vec::new());
        i
    }

    fn add_edge(&mut self, from: PNode, to: PNode, weight: f64) {
        let a = self.node(from);
        let b = self.node(to);
        self.adj[a].push((b, weight));
    }
}

// ---------------------------------------------------------------------------
// Engine entry point
// ---------------------------------------------------------------------------

/// Run PHANOTATE on a single contig. `dna_ascii` is the raw (any-case)
/// nucleotide sequence; `code` supplies the stop-codon set.
///
/// Returns the called genes (sorted by left coordinate) plus a pseudo-score
/// (the negated total path weight; higher is better) for optional model
/// comparison.
pub fn predict(
    contig_id: &str,
    dna_ascii: &[u8],
    code: &GeneticCode,
    min_orf_len: i64,
    closed: bool,
) -> (Vec<Gene>, f64) {
    // (U is read as T, as by every other engine.)
    let dna: Vec<u8> = dna_ascii
        .iter()
        .map(|b| match b.to_ascii_lowercase() {
            b'u' => b't',
            c => c,
        })
        .collect();
    let n = dna.len() as i64;
    if n < 3 {
        return (Vec::new(), 0.0);
    }

    // Start / stop codon sets ------------------------------------------------
    let start_weights: Vec<(&[u8], f64)> = vec![
        (b"atg".as_slice(), 1.0),
        (b"gtg".as_slice(), 0.10 / 0.85),
        (b"ttg".as_slice(), 0.05 / 0.85),
    ];
    let start_set: HashSet<[u8; 3]> = start_weights
        .iter()
        .map(|(c, _)| [c[0], c[1], c[2]])
        .collect();
    let stop_set: HashSet<[u8; 3]> = code
        .stop_codons_ascii()
        .into_iter()
        .map(|c| [c[0].to_ascii_lowercase(), c[1].to_ascii_lowercase(), c[2].to_ascii_lowercase()])
        .collect();

    let is_start = |c: &[u8]| c.len() == 3 && start_set.contains(&[c[0], c[1], c[2]]);
    let is_stop = |c: &[u8]| c.len() == 3 && stop_set.contains(&[c[0], c[1], c[2]]);
    let start_weight = |c: &[u8]| -> Option<f64> {
        for (sc, w) in &start_weights {
            if *sc == c {
                return Some(*w);
            }
        }
        None
    };

    // ----------------------------------------------------------------------
    // Pass 1: nucleotide frequency, background RBS, GC-frame plot
    // ----------------------------------------------------------------------
    let mut freq = [0.0f64; 4]; // a,c,g,t   (index by code)
    let mut frame_plot = GcFrame::default();
    let mut background_rbs = [1.0f64; 28];
    let mut training_rbs = [1.0f64; 28];

    for i in 0..dna.len() {
        // ambiguity mapping (only for frequency + GC plot, matching reference)
        let mapped = match dna[i] {
            b'a' | b'c' | b'g' | b't' => dna[i],
            b's' | b'b' | b'v' => b'g',
            _ => b'a',
        };
        match mapped {
            b'a' => {
                freq[0] += 1.0;
                freq[3] += 1.0;
            } // a + revcomp(t)
            b'c' => {
                freq[1] += 1.0;
                freq[2] += 1.0;
            } // c + revcomp(g)
            b'g' => {
                freq[2] += 1.0;
                freq[1] += 1.0;
            }
            b't' => {
                freq[3] += 1.0;
                freq[0] += 1.0;
            }
            _ => {}
        }
        // background RBS over both strands of the 21-mer at i
        let end = (i + 21).min(dna.len());
        let win = &dna[i..end];
        background_rbs[phanotate_score_rbs(win)] += 1.0;
        let rc = rev_comp(win);
        background_rbs[phanotate_score_rbs(&rc)] += 1.0;
        // GC frame plot: feed the encoded base (GC-ness matches the mapping)
        frame_plot.add_base(crate::seq::encode_base(mapped).0);
    }
    let gc_pos_freq = frame_plot.get();

    let pa = freq[0] / (n as f64 * 2.0);
    let pt = freq[3] / (n as f64 * 2.0);
    let pg = freq[2] / (n as f64 * 2.0);
    let pgap = pt * pa * pa + pt * pg * pa + pt * pa * pg; // genome-wide pstop
    let pstop_genome = pgap;

    let bg_sum: f64 = background_rbs.iter().sum();
    for x in background_rbs.iter_mut() {
        *x /= bg_sum;
    }

    // ----------------------------------------------------------------------
    // Pass 2: enumerate ORFs (faithful to get_orfs), one assembly-gap-free
    // segment at a time so that no ORF crosses a run of N. Coordinates stay
    // genome-wide (1-based); statistics and the graph remain whole-genome, as
    // in the reference. A contig without gaps is a single segment.
    // ----------------------------------------------------------------------
    let mut orfs = OrfStore::new(min_orf_len);

    // codon at 1-based position p (bytes p, p+1, p+2)
    let codon = |p: i64| -> [u8; 3] {
        let i = (p - 1) as usize; // 0-based first base
        [dna[i], dna[i + 1], dna[i + 2]]
    };
    // genome-wide reading frame (1..=3) of a 1-based position
    let frame_of = |p: i64| -> i32 { ((p - 1) % 3 + 1) as i32 };

    let rbs_fwd = |start: i64| -> usize {
        let lo = (start - 21).max(0) as usize;
        let hi = (start as usize).min(dna.len());
        if lo >= hi {
            return phanotate_score_rbs(&[]);
        }
        phanotate_score_rbs(&dna[lo..hi])
    };
    let rbs_rev = |start: i64| -> usize {
        let lo = start as usize;
        let hi = ((start + 21) as usize).min(dna.len());
        if lo >= hi {
            return phanotate_score_rbs(&[]);
        }
        phanotate_score_rbs(&rev_comp(&dna[lo..hi]))
    };

    let segments = crate::finalize::segments(dna_ascii);
    for &(off, seg_len) in &segments {
        if seg_len < 3 {
            continue;
        }
        let seg_lo = off as i64 + 1; // 1-based first base of the segment
        let seg_hi = (off + seg_len) as i64; // 1-based last base

        // starts[frame] for frames keyed -3..=-1, 1..=3
        let mut starts: HashMap<i32, Vec<i64>> = HashMap::new();
        let mut stops: HashMap<i32, i64> = HashMap::new();
        // Whether the left boundary held in `stops[-f]` is a real stop codon
        // (it starts out as the segment edge). With --closed, minus-strand
        // ORFs that would run off the left edge are not created.
        let mut rev_stop_real: HashMap<i32, bool> = HashMap::new();
        for f in [1, 2, 3] {
            starts.insert(f, Vec::new());
            starts.insert(-f, Vec::new());
            stops.insert(f, 0);
            rev_stop_real.insert(-f, false);
        }
        // left boundary of each minus frame = first segment position in it
        for k in 0..3i64 {
            stops.insert(-frame_of(seg_lo + k), seg_lo + k);
        }

        // edge starts at the segment's first codons (partial 5' ORFs) — not
        // with --closed
        for k in 0..3i64 {
            let p = seg_lo + k;
            if !closed && p + 2 <= seg_hi && !is_start(&codon(p)) {
                starts.get_mut(&frame_of(p)).unwrap().push(p);
            }
        }

        // main scan over every full codon of the segment
        let mut i = seg_lo;
        while i <= seg_hi - 2 {
            let cod = codon(i);
            let frame = frame_of(i);
            let rc = [rc_base(cod[2]), rc_base(cod[1]), rc_base(cod[0])];

            if is_start(&cod) {
                starts.get_mut(&frame).unwrap().push(i);
            } else if is_start(&rc) {
                starts.get_mut(&(-frame)).unwrap().push(i + 2);
            } else if is_stop(&cod) {
                let stop = i + 2;
                let pend: Vec<i64> = starts.get(&frame).unwrap().iter().rev().cloned().collect();
                for s in pend {
                    let length = stop - s + 1;
                    if length >= orfs.min_orf_len {
                        let seq = dna[(s - 1) as usize..stop as usize].to_vec();
                        let rs = rbs_fwd(s);
                        orfs.add_orf(s, stop - 2, frame, seq, rs);
                        training_rbs[rs] += 1.0;
                    }
                }
                starts.get_mut(&frame).unwrap().clear();
                stops.insert(frame, stop);
            } else if is_stop(&rc) {
                let stop = *stops.get(&(-frame)).unwrap();
                let left_ok = !closed || rev_stop_real[&(-frame)];
                let pend: Vec<i64> = if left_ok { starts.get(&(-frame)).unwrap().clone() } else { Vec::new() };
                for s in pend {
                    let length = s - stop + 1;
                    if length >= orfs.min_orf_len {
                        let lo = (stop - 1).max(0) as usize;
                        let seq = rev_comp(&dna[lo..s as usize]);
                        let rs = rbs_rev(s);
                        orfs.add_orf(s - 2, stop, -frame, seq, rs);
                        training_rbs[rs] += 1.0;
                    }
                }
                starts.get_mut(&(-frame)).unwrap().clear();
                stops.insert(-frame, i);
                rev_stop_real.insert(-frame, true);
            }
            i += 1;
        }

        // fragment ORFs at the end of the segment
        for frame in [1i32, 2, 3] {
            // end of the last full codon of this frame inside the segment
            let mut end_stop = seg_hi;
            while end_stop - 2 >= seg_lo && frame_of(end_stop - 2) != frame {
                end_stop -= 1;
            }
            if end_stop - 2 < seg_lo {
                continue;
            }
            // (3'-partial forward fragments are not created with --closed)
            let pend: Vec<i64> = if closed {
                Vec::new()
            } else {
                starts.get(&frame).unwrap().iter().rev().cloned().collect()
            };
            for s in pend {
                let length = end_stop - s + 1;
                if length >= orfs.min_orf_len && end_stop >= s {
                    let lo = (s - 1).max(0) as usize;
                    let hi = (end_stop as usize).min(dna.len());
                    if lo < hi {
                        let seq = dna[lo..hi].to_vec();
                        let rs = rbs_fwd(s);
                        orfs.add_orf(s, end_stop - 2, frame, seq, rs);
                        training_rbs[rs] += 1.0;
                    }
                }
            }
            // possible partial 3' minus ORF
            let s_edge = end_stop;
            let lo3 = (s_edge - 3).max(0) as usize;
            let hi3 = (s_edge as usize).min(dna.len());
            let edge_codon = if lo3 < hi3 { dna[lo3..hi3].to_vec() } else { Vec::new() };
            if !closed && !is_start(&rev_comp(&edge_codon)) {
                starts.get_mut(&(-frame)).unwrap().push(end_stop);
            }
            let left_ok = !closed || rev_stop_real[&(-frame)];
            let pend: Vec<i64> = if left_ok { starts.get(&(-frame)).unwrap().clone() } else { Vec::new() };
            for s in pend {
                let stop = *stops.get(&(-frame)).unwrap();
                let length = s - stop + 1;
                if length >= orfs.min_orf_len && s >= stop {
                    let lo = (stop - 1).max(0) as usize;
                    let hi = (s as usize).min(dna.len());
                    if lo < hi {
                        let seq = rev_comp(&dna[lo..hi]);
                        let rs = rbs_rev(s);
                        orfs.add_orf(s - 2, stop, -frame, seq, rs);
                        training_rbs[rs] += 1.0;
                    }
                }
            }
        }
    }

    if orfs.all.is_empty() {
        return (Vec::new(), 0.0);
    }

    // ----------------------------------------------------------------------
    // Score ORFs by RBS motif
    // ----------------------------------------------------------------------
    let tr_sum: f64 = training_rbs.iter().sum();
    for x in training_rbs.iter_mut() {
        *x /= tr_sum;
    }
    for orf in orfs.all.iter_mut() {
        orf.weight_rbs = training_rbs[orf.rbs_score] / background_rbs[orf.rbs_score];
    }

    // ----------------------------------------------------------------------
    // Score ORFs by GC frame plot
    // ----------------------------------------------------------------------
    let gc = |base: i64, k: usize| -> i64 {
        if base >= 0 && (base as usize) < gc_pos_freq.len() {
            gc_pos_freq[base as usize][k]
        } else {
            0
        }
    };
    let mut pos_max = [1.0f64; 4];
    let mut pos_min = [1.0f64; 4];

    for group in orfs.iter_in() {
        for &gi in &group {
            let (start, stop, sc) = {
                let o = &orfs.all[gi];
                (o.start, o.stop, o.start_codon().to_vec())
            };
            if sc == b"atg" {
                if start < stop {
                    let nn = ((stop - start) / 8) * 3;
                    let mut base = start + nn;
                    while base < stop - 36 {
                        pos_max[max_idx(gc(base, 0), gc(base, 1), gc(base, 2))] += 1.0;
                        pos_min[min_idx(gc(base, 0), gc(base, 1), gc(base, 2))] += 1.0;
                        base += 3;
                    }
                } else if stop < start {
                    let nn = ((start - stop) / 8) * 3;
                    let mut base = start - nn;
                    while base > stop + 36 {
                        pos_max[max_idx(gc(base, 2), gc(base, 1), gc(base, 0))] += 1.0;
                        pos_min[min_idx(gc(base, 2), gc(base, 1), gc(base, 0))] += 1.0;
                        base -= 3;
                    }
                }
                break; // first ATG ORF for this stop
            }
        }
    }
    let mx = pos_max.iter().cloned().fold(f64::MIN, f64::max);
    for x in pos_max.iter_mut() {
        *x /= mx;
    }
    let mn = pos_min.iter().cloned().fold(f64::MIN, f64::max);
    for x in pos_min.iter_mut() {
        *x /= mn;
    }

    for orf in orfs.all.iter_mut() {
        let (start, stop) = (orf.start, orf.stop);
        if orf.frame > 0 {
            let mut base = start;
            while base < stop {
                let im = max_idx(gc(base, 0), gc(base, 1), gc(base, 2));
                let in_ = min_idx(gc(base, 0), gc(base, 1), gc(base, 2));
                let factor = (1.0 - orf.pstop).powf(pos_max[im]).powf(pos_min[in_]);
                orf.hold = (orf.hold * factor).max(HOLD_FLOOR);
                base += 3;
            }
        } else {
            let mut base = start;
            while base > stop {
                let im = max_idx(gc(base, 2), gc(base, 1), gc(base, 0));
                let in_ = min_idx(gc(base, 2), gc(base, 1), gc(base, 0));
                let factor = (1.0 - orf.pstop).powf(pos_max[im]).powf(pos_min[in_]);
                orf.hold = (orf.hold * factor).max(HOLD_FLOOR);
                base -= 3;
            }
        }
    }

    // finalise per-ORF weight (Orf.score)
    for orf in orfs.all.iter_mut() {
        let mut s = 1.0 / orf.hold;
        if let Some(w) = start_weight(orf.start_codon()) {
            s *= w;
        }
        s *= orf.weight_rbs;
        if s > S_CAP {
            s = S_CAP;
        }
        orf.weight = -s;
    }

    // ----------------------------------------------------------------------
    // Build the graph (faithful to get_graph)
    // ----------------------------------------------------------------------
    let mut g = Graph::new();

    // ORF edges
    for orf in &orfs.all {
        if orf.frame > 0 {
            let src = PNode { typ: NType::Start, frame: orf.frame, position: orf.start };
            let tgt = PNode { typ: NType::Stop, frame: orf.frame, position: orf.stop };
            g.add_edge(src, tgt, orf.weight);
        } else {
            let src = PNode { typ: NType::Stop, frame: orf.frame, position: orf.stop };
            let tgt = PNode { typ: NType::Start, frame: orf.frame, position: orf.start };
            g.add_edge(src, tgt, orf.weight);
        }
    }

    // Long non-coding bridge edges
    {
        let mut bases = vec![false; (n as usize).max(1)];
        for group in orfs.iter_in() {
            if let Some(&gi) = group.first() {
                let o = &orfs.all[gi];
                let mi = o.start.min(o.stop);
                let ma = o.start.max(o.stop);
                let hi = ma.min(n - 1);
                let mut k = mi;
                while k < hi {
                    if k >= 0 && (k as usize) < bases.len() {
                        bases[k as usize] = true;
                    }
                    k += 1;
                }
            }
        }
        // snapshot node list (so we can mutate the graph while iterating)
        let snapshot: Vec<PNode> = g.nodes.clone();
        let mut last = 0i64;
        for pos in 0..bases.len() as i64 {
            if bases[pos as usize] {
                let base = pos;
                if base - last > 500 {
                    for right in &snapshot {
                        for left in &snapshot {
                            let l = left.position;
                            let r = right.position;
                            if last + 1 >= l && l > last - 500 && base - 1 <= r && r < base + 500 {
                                if left.frame * right.frame > 0 {
                                    if left.typ == NType::Stop
                                        && right.typ == NType::Start
                                        && left.frame > 0
                                    {
                                        let s = score_gap(r - l - 3, false, pgap);
                                        g.add_edge(*left, *right, s);
                                    } else if left.typ == NType::Start
                                        && right.typ == NType::Stop
                                        && left.frame < 0
                                    {
                                        let s = score_gap(r - l - 3, false, pgap);
                                        g.add_edge(*left, *right, s);
                                    }
                                } else {
                                    if left.typ == NType::Stop
                                        && right.typ == NType::Stop
                                        && left.frame > 0
                                    {
                                        let s = score_gap(r - l - 3, true, pgap);
                                        g.add_edge(*left, *right, s);
                                    } else if left.typ == NType::Start
                                        && right.typ == NType::Start
                                        && left.frame < 0
                                    {
                                        let s = score_gap(r - l - 3, true, pgap);
                                        g.add_edge(*left, *right, s);
                                    }
                                }
                            }
                        }
                    }
                }
                last = base;
            }
        }
    }

    // Connect ORFs to each other (the O(N^2) connection loop, windowed by 500)
    {
        // sort nodes by position for the sliding window
        let snapshot: Vec<PNode> = g.nodes.clone();
        let mut order: Vec<usize> = (0..snapshot.len()).collect();
        order.sort_by_key(|&k| snapshot[k].position);

        let pstop_of = |pos: i64| -> f64 {
            // o1/o2 lookup as in get_graph
            if orfs.is_stop(pos) {
                let oe = *orfs.other_end.get(&pos).unwrap_or(&pos);
                if let Some(o) = orfs.get_orf(oe, pos) {
                    return o.pstop;
                }
                if let Some(&oe2) = orfs.other_end.get(&pos) {
                    if let Some(o) = orfs.get_orf(pos, oe2) {
                        return o.pstop;
                    }
                }
                pgap
            } else {
                pgap
            }
        };

        for oi in 0..order.len() {
            let right = snapshot[order[oi]];
            let r = right.position;
            let r_other = *orfs.other_end.get(&r).unwrap_or(&r);
            // left nodes have position in (r-500, r); scan backwards
            let mut oj = oi;
            while oj > 0 {
                oj -= 1;
                let left = snapshot[order[oj]];
                let l = left.position;
                if r - l >= 500 {
                    break;
                }
                if !(0 < r - l && r - l < 500) {
                    continue;
                }
                let l_other = *orfs.other_end.get(&l).unwrap_or(&l);
                let o1 = pstop_of(l);
                let o2 = pstop_of(r);
                let pstop = (o1 + o2) / 2.0;

                if left.frame * right.frame > 0 {
                    // same directions
                    if left.typ == NType::Stop && right.typ == NType::Start {
                        if left.frame > 0 {
                            let s = score_gap(r - l - 3, false, pgap);
                            g.add_edge(left, right, s);
                        } else if left.frame != right.frame && r < l_other && r_other < l {
                            let s = score_overlap(r - l + 3, false, pstop);
                            g.add_edge(right, left, s);
                        }
                    }
                    if left.typ == NType::Start && right.typ == NType::Stop {
                        if left.frame > 0 {
                            if left.frame != right.frame && r < l_other && r_other < l {
                                let s = score_overlap(r - l + 3, false, pstop);
                                g.add_edge(right, left, s);
                            }
                        } else {
                            let s = score_gap(r - l - 3, false, pgap);
                            g.add_edge(left, right, s);
                        }
                    }
                } else {
                    // different directions
                    if left.typ == NType::Stop && right.typ == NType::Stop {
                        if right.frame > 0 {
                            if r_other + 3 < l && r < l_other {
                                let s = score_overlap(r - l + 3, true, pstop);
                                g.add_edge(right, left, s);
                            }
                        } else {
                            let s = score_gap(r - l - 3, true, pgap);
                            g.add_edge(left, right, s);
                        }
                    }
                    if left.typ == NType::Start && right.typ == NType::Start {
                        if right.frame > 0 && r - l > 2 {
                            let s = score_gap(r - l - 3, true, pgap);
                            g.add_edge(left, right, s);
                        } else if right.frame < 0 && r_other < l && r < l_other {
                            let s = score_overlap(r - l + 3, true, pstop);
                            g.add_edge(right, left, s);
                        }
                    }
                }
            }
        }
    }

    // Source / target. The reference links them to ORF ends within 2 kb of the
    // contig ends. Those distances are measured from the outermost gap-free
    // segments that hold an ORF — identical for a contig without terminal
    // gaps, while a leading / trailing assembly gap (no ORF can reach into it)
    // no longer leaves the path without an entry or exit, which returned no
    // genes at all. Should no ORF end lie within the window, the window is
    // anchored on the outermost ORF end instead.
    let source = PNode { typ: NType::Source, frame: 0, position: 0 };
    let target = PNode { typ: NType::Target, frame: 0, position: n + 1 };
    let src_idx = g.node(source);
    let tgt_idx = g.node(target);
    {
        let snapshot: Vec<PNode> = g.nodes.clone();
        let left_end = |nd: &PNode| (nd.typ == NType::Start && nd.frame > 0) || (nd.typ == NType::Stop && nd.frame < 0);
        let right_end = |nd: &PNode| (nd.typ == NType::Start && nd.frame < 0) || (nd.typ == NType::Stop && nd.frame > 0);
        // (0-based offset, 1-based last base) of the segment holding position p
        let seg_bounds = |p: i64| -> (i64, i64) {
            let k = segments.partition_point(|&(o, _)| (o as i64) < p).saturating_sub(1);
            match segments.get(k) {
                Some(&(o, l)) => (o as i64, (o + l) as i64),
                None => (0, n),
            }
        };
        let (lo_ref, lo_win) = match snapshot.iter().filter(|nd| left_end(nd)).map(|nd| nd.position).min() {
            Some(p) => {
                let lo = seg_bounds(p).0;
                (lo, if p - lo <= 2000 { lo + 2000 } else { p + 2000 })
            }
            None => (0, 2000),
        };
        let (hi_ref, hi_win) = match snapshot.iter().filter(|nd| right_end(nd)).map(|nd| nd.position).max() {
            Some(p) => {
                let hi = seg_bounds(p).1;
                (hi, if hi - p <= 2000 { hi - 2000 } else { p - 2000 })
            }
            None => (n, n - 2000),
        };
        for node in &snapshot {
            if node.typ == NType::Source || node.typ == NType::Target {
                continue;
            }
            if node.position <= lo_win && left_end(node) {
                let s = score_gap(node.position - lo_ref, false, pgap);
                g.add_edge(source, *node, s);
            }
            if node.position >= hi_win && right_end(node) {
                let s = score_gap(hi_ref - node.position, false, pgap);
                g.add_edge(*node, target, s);
            }
        }
    }

    // ----------------------------------------------------------------------
    // Shortest path: SPFA (queue-based Bellman-Ford-Moore)
    // ----------------------------------------------------------------------
    let nv = g.nodes.len();
    let mut dist = vec![f64::INFINITY; nv];
    let mut pred = vec![usize::MAX; nv];
    let mut in_queue = vec![false; nv];
    let mut count = vec![0u32; nv];
    let mut queue: VecDeque<usize> = VecDeque::new();
    dist[src_idx] = 0.0;
    in_queue[src_idx] = true;
    queue.push_back(src_idx);
    let guard = (nv as u32).saturating_add(1);
    let mut aborted = false;

    while let Some(u) = queue.pop_front() {
        in_queue[u] = false;
        let du = dist[u];
        // clone out-edges to avoid borrow issues
        let edges = g.adj[u].clone();
        for (v, w) in edges {
            let nd = du + w;
            if nd < dist[v] {
                dist[v] = nd;
                pred[v] = u;
                if !in_queue[v] {
                    in_queue[v] = true;
                    count[v] += 1;
                    if count[v] > guard {
                        aborted = true;
                        break;
                    }
                    queue.push_back(v);
                }
            }
        }
        if aborted {
            break;
        }
    }

    if dist[tgt_idx].is_infinite() {
        return (Vec::new(), 0.0);
    }

    // reconstruct node path source -> target
    let mut path_idx: Vec<usize> = Vec::new();
    let mut cur = tgt_idx;
    let mut steps = 0usize;
    let mut seen = vec![false; nv];
    while cur != usize::MAX {
        if seen[cur] {
            break; // cycle guard
        }
        seen[cur] = true;
        path_idx.push(cur);
        if cur == src_idx {
            break;
        }
        cur = pred[cur];
        steps += 1;
        if steps > nv + 1 {
            break;
        }
    }
    path_idx.reverse();
    if path_idx.is_empty() || path_idx[0] != src_idx {
        return (Vec::new(), 0.0);
    }

    // drop source, then take genes in non-overlapping pairs
    let body: Vec<PNode> = path_idx[1..].iter().map(|&k| g.nodes[k]).collect();

    // build a Seq once for sequence filling
    let encoded = Seq::from_ascii(&dna);

    let mut genes: Vec<Gene> = Vec::new();
    let mut total_weight = 0.0f64;
    let mut gi = 0usize;
    let mut k = 0usize;
    while k + 1 < body.len() {
        let left = body[k];
        let right = body[k + 1];
        k += 2;
        if left.typ == NType::Target || right.typ == NType::Target {
            break;
        }

        // identify orf & strand
        let (orf_opt, strand): (Option<&Orf>, i8) = if left.frame > 0 {
            (orfs.get_orf(left.position, right.position), 1)
        } else {
            (orfs.get_orf(right.position, left.position), -1)
        };
        let orf = match orf_opt {
            Some(o) => o,
            None => continue,
        };

        let has_start = is_start(orf.start_codon());
        let has_stop = is_stop(orf.stop_codon());

        // partial flags (see module docs / write_output mapping)
        let (mut partial_left, mut partial_right);
        if strand > 0 {
            partial_left = !has_start;
            partial_right = !has_stop;
        } else {
            partial_left = !has_stop;
            partial_right = !has_start;
        }

        // coordinates
        let (mut begin, mut end) = if strand > 0 {
            (orf.start, orf.stop + 2)
        } else {
            (orf.stop, orf.start + 2)
        };
        if begin < 1 {
            begin = 1;
            partial_left = true;
        }
        if end > n {
            end = n;
            partial_right = true;
        }
        if begin > end {
            std::mem::swap(&mut begin, &mut end);
        }

        // start type
        let start_type = if (strand > 0 && !has_start) || (strand < 0 && !has_start) {
            3
        } else {
            let sc = orf.start_codon();
            if sc == b"atg" {
                0
            } else if sc == b"gtg" {
                1
            } else if sc == b"ttg" {
                2
            } else {
                3
            }
        };

        total_weight += orf.weight;
        gi += 1;
        let mut gene = Gene {
            contig: contig_id.to_string(),
            id: gi,
            begin: begin as usize,
            end: end as usize,
            strand,
            start_type,
            rbs_motif: "None".to_string(),
            rbs_spacer: "None".to_string(),
            partial_left,
            partial_right,
            score: orf.weight,
            cscore: 0.0,
            sscore: 0.0,
            conf: 0.0,
            gc: 0.0,
            trans_table: 0,
            engine: Engine::Phanotate,
            aa: String::new(),
            nt: String::new(),
        };
        gene.fill_sequences(&encoded, code);
        // GC content of the gene
        gene.gc = gc_of(&dna, gene.begin, gene.end);
        gene.trans_table = code.id;
        genes.push(gene);
    }

    // Shared finalisation, per gap-free segment (whose ends are open edges):
    // gap-codon trimming, codon-based partial inference (Prodigal's genomic
    // partial=XY), Edge start type, X for ambiguous codons. PHANOTATE opens
    // genes at ATG/GTG/TTG under every table.
    let mut by_seg: Vec<Vec<Gene>> = vec![Vec::new(); segments.len()];
    for g in genes {
        let k = segments.partition_point(|&(o, _)| o < g.begin).saturating_sub(1);
        by_seg[k].push(g);
    }
    let mut genes: Vec<Gene> = Vec::new();
    for (k, mut local) in by_seg.into_iter().enumerate() {
        if local.is_empty() {
            continue;
        }
        let (off, len) = segments[k];
        let seg_seq = if off == 0 && len == dna_ascii.len() {
            encoded.clone()
        } else {
            Seq::from_ascii(&dna_ascii[off..off + len])
        };
        for g in local.iter_mut() {
            g.begin -= off;
            g.end -= off;
        }
        let mut done = crate::finalize::finalize_genes(
            local,
            &seg_seq,
            code,
            closed,
            crate::finalize::SeqSource::Genome,
            crate::finalize::StartRule::Bacterial,
        );
        crate::finalize::shift_genes(&mut done, off, contig_id);
        genes.extend(done);
    }
    for g in genes.iter_mut() {
        g.gc = gc_of(&dna, g.begin, g.end);
    }
    genes.sort_by_key(|g| (g.begin, g.end));
    for (j, gene) in genes.iter_mut().enumerate() {
        gene.id = j + 1;
    }

    let _ = pstop_genome;
    (genes, -total_weight)
}

fn gc_of(dna: &[u8], begin: usize, end: usize) -> f64 {
    if begin == 0 || end == 0 || begin > end || end > dna.len() {
        return 0.0;
    }
    let mut gc = 0usize;
    let mut total = 0usize;
    for &b in &dna[begin - 1..end] {
        match b {
            b'g' | b'c' => {
                gc += 1;
                total += 1;
            }
            b'a' | b't' => total += 1,
            _ => {}
        }
    }
    if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scores_basic() {
        // score_gap with zero-length gap is 1; with diff adds 1/0.05
        let g = score_gap(0, false, 0.05);
        assert!((g - 1.0).abs() < 1e-9);
        let gd = score_gap(0, true, 0.05);
        assert!((gd - (1.0 + 20.0)).abs() < 1e-9);
        // overlap grows with length
        let o1 = score_overlap(10, false, 0.05);
        let o2 = score_overlap(20, false, 0.05);
        assert!(o2 > o1);
    }

    #[test]
    fn runs_and_returns_valid_genes() {
        // A degenerate synthetic sequence won't have the gene layout of a real
        // genome (artificial frames are often stop-free), so we don't assert a
        // specific call here — only that the engine runs, terminates, and emits
        // internally consistent genes. Real-genome correctness is covered by the
        // phiX174 / lambda integration tests.
        let mut dna = b"taa".to_vec();
        dna.extend_from_slice(b"atg");
        for _ in 0..40 {
            dna.extend_from_slice(b"aaa");
        }
        dna.extend_from_slice(b"taa");
        dna.extend_from_slice(b"acgtacgtacgtacgtacgtacgt");
        let n = dna.len();

        let code = GeneticCode::new(11).unwrap();
        let (genes, _score) = predict("test", &dna, &code, 90, true);
        assert!(!genes.is_empty(), "expected at least one gene");
        for g in &genes {
            assert!(g.begin >= 1 && g.end <= n, "coords out of bounds: {}..{}", g.begin, g.end);
            assert!(g.begin <= g.end, "begin past end");
            assert!(g.strand == 1 || g.strand == -1);
            assert!(g.start_type <= 3);
        }
    }
}
