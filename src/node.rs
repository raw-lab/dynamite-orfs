//! Prodigal-style node model: start/stop enumeration, the self-training
//! pipeline (GC-frame bias, hexamer coding model, Shine-Dalgarno start
//! training, upstream base-composition) and node scoring.
//!
//! This is a faithful Rust port of the core of Prodigal / Prodigal-GV
//! (`node.c`, Doug Hyatt et al., GPL-3.0). The dynamic-programming traceback
//! itself lives in [`crate::prodigal`].
//!
//! ## Faithful, with one documented simplification
//!
//! Prodigal can run one of two start trainers: a Shine-Dalgarno (SD) trainer
//! and a de-novo upstream-motif trainer for organisms that do not use an SD
//! motif. DYNAMITE implements the SD trainer in full (RBS-bin weights, ATG/
//! GTG/TTG type weights, and -1/-2/-15..-44 upstream base composition) and
//! always scores starts through that path. Organisms that Prodigal would flag
//! as non-SD are still scored with the (strong) combination of coding score,
//! start-type weight, SD-RBS score and upstream base composition; the only
//! piece not reproduced is de-novo discovery of a *non-canonical* RBS motif.
//! See the project README for the rationale.

use crate::gcframe::{max_fr, most_gc_frame};
use crate::gencode::GeneticCode;
use crate::seq::{is_gc, mer_ndx, Seq, A, G, T};

/// Minimum gene length (bp) for an internal start->stop ORF.
pub const MIN_GENE: i64 = 90;
/// Minimum gene length (bp) for an ORF that runs off a contig edge.
pub const MIN_EDGE_GENE: i64 = 60;
/// Maximum same-strand overlap (bp) recorded between a stop and nearby starts.
pub const MAX_SAM_OVLP: i64 = 60;
/// Operon distance constant used by the intergenic modifier.
pub const OPER_DIST: f64 = 60.0;
/// Bonus applied to genes that run off a contig edge.
pub const EDGE_BONUS: f64 = 0.74;
/// Upstream penalty applied near contig edges.
pub const EDGE_UPS: f64 = -1.00;
/// Coding penalty for short metagenomic fragments.
pub const META_PEN: f64 = 7.5;
/// Start weight constant (Prodigal's `tinf.st_wt`).
pub const ST_WT: f64 = 4.35;

/// Node type: start codons and the stop sentinel.
pub const ATG: i32 = 0;
pub const GTG: i32 = 1;
pub const TTG: i32 = 2;
pub const STOP: i32 = 3;

/// The trained model parameters (Prodigal's `struct _training`).
#[derive(Clone)]
pub struct Tinf {
    /// Genome GC fraction.
    pub gc: f64,
    /// NCBI translation table.
    pub trans_table: u32,
    /// Start weight (constant 4.35).
    pub st_wt: f64,
    /// GC-frame bias for the three codon positions.
    pub bias: [f64; 3],
    /// Log weights for ATG / GTG / TTG.
    pub type_wt: [f64; 3],
    /// Whether the organism uses a Shine-Dalgarno motif (reported only).
    pub uses_sd: bool,
    /// Log weights for the 28 RBS bins.
    pub rbs_wt: [f64; 28],
    /// Upstream base-composition log weights (32 positions x 4 bases).
    pub ups_comp: [[f64; 4]; 32],
    /// Hexamer (6-mer) coding log-likelihoods.
    pub gene_dc: Vec<f64>,
}

impl Tinf {
    /// A fresh model for a given genetic code and GC content.
    pub fn new(trans_table: u32, gc: f64) -> Tinf {
        Tinf {
            gc,
            trans_table,
            st_wt: ST_WT,
            bias: [0.0; 3],
            type_wt: [0.0; 3],
            uses_sd: true,
            rbs_wt: [0.0; 28],
            ups_comp: [[0.0; 4]; 32],
            gene_dc: vec![0.0; 4096],
        }
    }

    /// The genetic code object for this model.
    pub fn code(&self) -> GeneticCode {
        GeneticCode::new(self.trans_table).expect("validated translation table")
    }
}

/// One start or stop node (Prodigal's `struct _node`).
#[derive(Clone)]
pub struct Node {
    /// 0=ATG, 1=GTG, 2=TTG, 3=Stop.
    pub typ: i32,
    /// Runs off a contig edge.
    pub edge: bool,
    /// Position in the (forward) sequence of the node.
    pub ndx: i64,
    /// 1 = forward, -1 = reverse.
    pub strand: i32,
    /// For a start, the index of the stop that terminates its ORF; for a stop,
    /// the position of the previous in-frame stop.
    pub stop_val: i64,
    /// Indices of nearby in-frame starts within `MAX_SAM_OVLP` of a stop.
    pub star_ptr: [i64; 3],
    /// Frame of highest GC content within this node.
    pub gc_bias: usize,
    /// Per-codon-position GC scores.
    pub gc_score: [f64; 3],
    /// Coding score from the hexamer model.
    pub cscore: f64,
    /// GC content of the node's ORF.
    pub gc_cont: f64,
    /// Best SD RBS bin: `[0]` exact, `[1]` single-mismatch.
    pub rbs: [usize; 2],
    /// Upstream composition score.
    pub uscore: f64,
    /// Start-type score.
    pub tscore: f64,
    /// RBS-motif score.
    pub rscore: f64,
    /// Total start score (`tscore + rscore + uscore`).
    pub sscore: f64,
    /// Traceback pointer (dynamic programming).
    pub traceb: i64,
    /// Forward trace pointer.
    pub tracef: i64,
    /// Overlap-untangling marker.
    pub ov_mark: i64,
    /// Total path score to this node.
    pub score: f64,
    /// Eliminate this gene from the final set.
    pub elim: bool,
}

impl Default for Node {
    fn default() -> Node {
        Node {
            typ: STOP,
            edge: false,
            ndx: 0,
            strand: 1,
            stop_val: 0,
            star_ptr: [-1, -1, -1],
            gc_bias: 0,
            gc_score: [0.0; 3],
            cscore: 0.0,
            gc_cont: 0.0,
            rbs: [0, 0],
            uscore: 0.0,
            tscore: 0.0,
            rscore: 0.0,
            sscore: 0.0,
            traceb: -1,
            tracef: -1,
            ov_mark: -1,
            score: 0.0,
            elim: false,
        }
    }
}

#[inline]
fn base_at(seq: &[u8], i: i64) -> Option<u8> {
    if i < 0 {
        return None;
    }
    seq.get(i as usize).copied()
}

/// `true` if any base of the codon at `i` is ambiguous (N / IUPAC) or the
/// codon runs off the sequence. Such a codon is never a start or a stop.
#[inline]
fn codon_amb(amb: &[bool], i: i64) -> bool {
    if i < 0 {
        return true;
    }
    let u = i as usize;
    match (amb.get(u), amb.get(u + 1), amb.get(u + 2)) {
        (Some(&a), Some(&b), Some(&c)) => a || b || c,
        _ => true,
    }
}

#[inline]
fn is_stop_at(seq: &[u8], amb: &[bool], i: i64, code: &GeneticCode) -> bool {
    if codon_amb(amb, i) {
        return false;
    }
    match (base_at(seq, i), base_at(seq, i + 1), base_at(seq, i + 2)) {
        (Some(a), Some(b), Some(c)) => code.is_stop(a, b, c),
        _ => false,
    }
}

#[inline]
fn is_start_at(seq: &[u8], amb: &[bool], i: i64, code: &GeneticCode) -> bool {
    if codon_amb(amb, i) {
        return false;
    }
    match (base_at(seq, i), base_at(seq, i + 1), base_at(seq, i + 2)) {
        (Some(a), Some(b), Some(c)) => code.is_start(a, b, c),
        _ => false,
    }
}

#[inline]
fn codon_is(seq: &[u8], i: i64, x: u8, y: u8, z: u8) -> bool {
    base_at(seq, i) == Some(x) && base_at(seq, i + 1) == Some(y) && base_at(seq, i + 2) == Some(z)
}

#[inline]
fn is_atg_at(seq: &[u8], i: i64) -> bool {
    codon_is(seq, i, A, T, G)
}
#[inline]
fn is_gtg_at(seq: &[u8], i: i64) -> bool {
    codon_is(seq, i, G, T, G)
}
#[inline]
fn is_ttg_at(seq: &[u8], i: i64) -> bool {
    codon_is(seq, i, T, T, G)
}

#[inline]
fn dmin(a: f64, b: f64) -> f64 {
    a.min(b)
}
#[inline]
fn dmax(a: f64, b: f64) -> f64 {
    a.max(b)
}

/// Enumerate every start and stop node on both strands (Prodigal `add_nodes`).
pub fn add_nodes(seq: &Seq, code: &GeneticCode, closed: bool) -> Vec<Node> {
    add_nodes_masked(seq, code, closed, &[])
}

/// `true` if the half-open interval `[x, y)` overlaps a masked region —
/// pyrodigal's `Mask.intersects`, Prodigal's `cross_mask`. `masks` are sorted,
/// non-overlapping half-open `[begin, end)` intervals.
#[inline]
fn crosses_mask(masks: &[(usize, usize)], x: i64, y: i64) -> bool {
    if masks.is_empty() || y <= 0 {
        return false;
    }
    let x = x.max(0) as usize;
    let k = masks.partition_point(|&(_, e)| e <= x);
    k < masks.len() && (masks[k].0 as i64) < y
}

/// [`add_nodes`] with Prodigal 2.x's `-m` masking: no gene (start node) is
/// created whose span from start to stop crosses a masked region, so no gene
/// can be built across — or run into — a masked run of ambiguous bases.
/// `masks` are sorted half-open `[begin, end)` intervals
/// (see [`crate::seq::ambiguous_runs`]); with none this is exactly
/// [`add_nodes`].
pub fn add_nodes_masked(seq: &Seq, code: &GeneticCode, closed: bool, masks: &[(usize, usize)]) -> Vec<Node> {
    add_nodes_ends(seq, code, closed, closed, masks)
}

/// [`add_nodes_masked`] with each end of the sequence open or closed on its
/// own: `closed_left` / `closed_right` forbid genes running off the low /
/// high-coordinate end (Prodigal's `closed` closes both). Used for gap
/// mode 1, where a segment's ends at an assembly gap are closed while the
/// contig's own ends follow `--closed`.
pub fn add_nodes_ends(
    seq: &Seq,
    code: &GeneticCode,
    closed_left: bool,
    closed_right: bool,
    masks: &[(usize, usize)],
) -> Vec<Node> {
    let slen = seq.len() as i64;
    let fwd = &seq.fwd;
    let rev = &seq.rev;
    // Ambiguity masks: a codon touching an N / IUPAC base is never a start or
    // a stop (previously N was read as A, so e.g. `TNN` became a `TAA` stop
    // and genes could end inside assembly gaps).
    let famb: &[bool] = &seq.n;
    let ramb_v = seq.n_rev();
    let ramb: &[bool] = &ramb_v;
    let mut nodes: Vec<Node> = Vec::new();
    if slen < 3 {
        return nodes;
    }

    let mut last = [0i64; 3];
    let mut saw_start = [false; 3];
    let mut min_dist = [0i64; 3];
    let slmod = (slen % 3) as usize;

    // ---- Forward strand ----
    // (3' ends run off the right end; edge starts sit at the left end.)
    for i in 0..3i64 {
        let f = ((i + slmod as i64) % 3) as usize;
        last[f] = slen + i;
        saw_start[(i % 3) as usize] = false;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if !closed_right {
            while last[f] + 2 > slen - 1 {
                last[f] -= 3;
            }
        }
    }
    let mut i = slen - 3;
    while i >= 0 {
        let fr = (i % 3) as usize;
        if is_stop_at(fwd, famb, i, code) {
            if saw_start[fr] {
                let mut n = Node::default();
                if !is_stop_at(fwd, famb, last[fr], code) {
                    n.edge = true;
                }
                n.ndx = last[fr];
                n.typ = STOP;
                n.strand = 1;
                n.stop_val = i;
                nodes.push(n);
            }
            min_dist[fr] = MIN_GENE;
            last[fr] = i;
            saw_start[fr] = false;
            i -= 1;
            continue;
        }
        if last[fr] >= slen {
            i -= 1;
            continue;
        }
        if crosses_mask(masks, i, last[fr]) {
            i -= 1;
            continue;
        }
        let long_enough = (last[fr] - i + 3) >= min_dist[fr];
        if is_start_at(fwd, famb, i, code) && is_atg_at(fwd, i) && long_enough {
            push_start(&mut nodes, i, ATG, last[fr], 1, &mut saw_start, fr);
        } else if is_start_at(fwd, famb, i, code) && is_gtg_at(fwd, i) && long_enough {
            push_start(&mut nodes, i, GTG, last[fr], 1, &mut saw_start, fr);
        } else if is_start_at(fwd, famb, i, code) && is_ttg_at(fwd, i) && long_enough {
            push_start(&mut nodes, i, TTG, last[fr], 1, &mut saw_start, fr);
        } else if i <= 2 && !closed_left && (last[fr] - i) > MIN_EDGE_GENE {
            let mut n = Node::default();
            n.ndx = i;
            n.typ = ATG;
            n.edge = true;
            n.stop_val = last[fr];
            n.strand = 1;
            saw_start[fr] = true;
            nodes.push(n);
        }
        i -= 1;
    }
    for k in 0..3usize {
        if saw_start[k] {
            let mut n = Node::default();
            if !is_stop_at(fwd, famb, last[k], code) {
                n.edge = true;
            }
            n.ndx = last[k];
            n.typ = STOP;
            n.strand = 1;
            n.stop_val = k as i64 - 6;
            nodes.push(n);
        }
    }

    // ---- Reverse strand ----
    // (3' ends run off the left end; edge starts sit at the right end.)
    for i in 0..3i64 {
        let f = ((i + slmod as i64) % 3) as usize;
        last[f] = slen + i;
        saw_start[(i % 3) as usize] = false;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if !closed_left {
            while last[f] + 2 > slen - 1 {
                last[f] -= 3;
            }
        }
    }
    let mut i = slen - 3;
    while i >= 0 {
        let fr = (i % 3) as usize;
        if is_stop_at(rev, ramb, i, code) {
            if saw_start[fr] {
                let mut n = Node::default();
                if !is_stop_at(rev, ramb, last[fr], code) {
                    n.edge = true;
                }
                n.ndx = slen - last[fr] - 1;
                n.typ = STOP;
                n.strand = -1;
                n.stop_val = slen - i - 1;
                nodes.push(n);
            }
            min_dist[fr] = MIN_GENE;
            last[fr] = i;
            saw_start[fr] = false;
            i -= 1;
            continue;
        }
        if last[fr] >= slen {
            i -= 1;
            continue;
        }
        if crosses_mask(masks, slen - last[fr] - 1, slen - i - 1) {
            i -= 1;
            continue;
        }
        let long_enough = (last[fr] - i + 3) >= min_dist[fr];
        let nd = slen - i - 1;
        let sv = slen - last[fr] - 1;
        if is_start_at(rev, ramb, i, code) && is_atg_at(rev, i) && long_enough {
            push_start(&mut nodes, nd, ATG, sv, -1, &mut saw_start, fr);
        } else if is_start_at(rev, ramb, i, code) && is_gtg_at(rev, i) && long_enough {
            push_start(&mut nodes, nd, GTG, sv, -1, &mut saw_start, fr);
        } else if is_start_at(rev, ramb, i, code) && is_ttg_at(rev, i) && long_enough {
            push_start(&mut nodes, nd, TTG, sv, -1, &mut saw_start, fr);
        } else if i <= 2 && !closed_right && (last[fr] - i) > MIN_EDGE_GENE {
            let mut n = Node::default();
            n.ndx = nd;
            n.typ = ATG;
            n.edge = true;
            n.stop_val = sv;
            n.strand = -1;
            saw_start[fr] = true;
            nodes.push(n);
        }
        i -= 1;
    }
    for k in 0..3usize {
        if saw_start[k] {
            let mut n = Node::default();
            if !is_stop_at(rev, ramb, last[k], code) {
                n.edge = true;
            }
            n.ndx = slen - last[k] - 1;
            n.typ = STOP;
            n.strand = -1;
            n.stop_val = slen - k as i64 + 5;
            nodes.push(n);
        }
    }

    nodes
}

#[allow(clippy::too_many_arguments)]
fn push_start(
    nodes: &mut Vec<Node>,
    ndx: i64,
    typ: i32,
    stop_val: i64,
    strand: i32,
    saw_start: &mut [bool; 3],
    fr: usize,
) {
    let mut n = Node::default();
    n.ndx = ndx;
    n.typ = typ;
    n.stop_val = stop_val;
    n.strand = strand;
    nodes.push(n);
    saw_start[fr] = true;
}

/// Sort nodes by position (asc) then strand (fwd before rev) (`compare_nodes`).
pub fn sort_nodes(nodes: &mut [Node]) {
    nodes.sort_by(|a, b| {
        a.ndx
            .cmp(&b.ndx)
            .then(b.strand.cmp(&a.strand))
    });
}

/// Record GC-frame bias for each ORF and set `tinf.bias` (`record_gc_bias`).
pub fn record_gc_bias(gc: &[u8], nodes: &mut [Node], tinf: &mut Tinf) {
    let nn = nodes.len();
    if nn == 0 {
        return;
    }
    let mut ctr = [[0i64; 3]; 3];
    let mut last = [0i64; 3];
    // Prodigal stores -1 for the trailing bases of the GC frame plot.
    let gcv = |p: i64| -> i64 {
        match gc.get(p as usize) {
            Some(&v) if v != crate::gcframe::GC_FRAME_UNSET => v as i64,
            _ => -1,
        }
    };

    // Forward pass (high index to low).
    for idx in (0..nn).rev() {
        let fr = (nodes[idx].ndx % 3) as usize;
        let frmod = 3 - fr as i64;
        if nodes[idx].strand == 1 && nodes[idx].typ == STOP {
            for j in 0..3 {
                ctr[fr][j] = 0;
            }
            last[fr] = nodes[idx].ndx;
            let g = gcv(nodes[idx].ndx);
            ctr[fr][((g + frmod).rem_euclid(3)) as usize] = 1;
        } else if nodes[idx].strand == 1 {
            let mut j = last[fr] - 3;
            while j >= nodes[idx].ndx {
                let g = gcv(j);
                ctr[fr][((g + frmod).rem_euclid(3)) as usize] += 1;
                j -= 3;
            }
            let mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]) as usize;
            nodes[idx].gc_bias = mfr;
            let denom = (nodes[idx].stop_val - nodes[idx].ndx + 3) as f64;
            for j in 0..3 {
                nodes[idx].gc_score[j] = (3.0 * ctr[fr][j] as f64) / denom;
            }
            last[fr] = nodes[idx].ndx;
        }
    }

    // Reverse pass (low index to high).
    for idx in 0..nn {
        let fr = (nodes[idx].ndx % 3) as usize;
        let frmod = fr;
        if nodes[idx].strand == -1 && nodes[idx].typ == STOP {
            for j in 0..3 {
                ctr[fr][j] = 0;
            }
            last[fr] = nodes[idx].ndx;
            let g = gcv(nodes[idx].ndx);
            ctr[fr][(((3 - g) + frmod as i64).rem_euclid(3)) as usize] = 1;
        } else if nodes[idx].strand == -1 {
            let mut j = last[fr] + 3;
            while j <= nodes[idx].ndx {
                let g = gcv(j);
                ctr[fr][(((3 - g) + frmod as i64).rem_euclid(3)) as usize] += 1;
                j += 3;
            }
            let mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]) as usize;
            nodes[idx].gc_bias = mfr;
            let denom = (nodes[idx].ndx - nodes[idx].stop_val + 3) as f64;
            for j in 0..3 {
                nodes[idx].gc_score[j] = (3.0 * ctr[fr][j] as f64) / denom;
            }
            last[fr] = nodes[idx].ndx;
        }
    }

    tinf.bias = [0.0; 3];
    for idx in 0..nn {
        if nodes[idx].typ != STOP {
            let len = (nodes[idx].stop_val - nodes[idx].ndx).abs() + 1;
            let gb = nodes[idx].gc_bias;
            tinf.bias[gb] += (nodes[idx].gc_score[gb] * len as f64) / 1000.0;
        }
    }
    let tot = tinf.bias[0] + tinf.bias[1] + tinf.bias[2];
    if tot != 0.0 {
        for j in 0..3 {
            tinf.bias[j] *= 3.0 / tot;
        }
    }
}

/// Record nearby in-frame starts for each stop (`record_overlapping_starts`).
///
/// `flag = false` is used during training (first observed start), `flag =
/// true` during the final pass (best-scoring start).
pub fn record_overlapping_starts(nodes: &mut [Node], tinf: &Tinf, flag: bool) {
    let nn = nodes.len() as i64;
    for i in 0..nn as usize {
        for j in 0..3 {
            nodes[i].star_ptr[j] = -1;
        }
        if nodes[i].typ != STOP || nodes[i].edge {
            continue;
        }
        let i_ndx = nodes[i].ndx;
        if nodes[i].strand == 1 {
            let mut max_sc = -100.0f64;
            let mut j = i as i64 + 3;
            while j >= 0 {
                if j >= nn {
                    j -= 1;
                    continue;
                }
                let ju = j as usize;
                if nodes[ju].ndx > i_ndx + 2 {
                    j -= 1;
                    continue;
                }
                if nodes[ju].ndx + MAX_SAM_OVLP < i_ndx {
                    break;
                }
                if nodes[ju].strand == 1 && nodes[ju].typ != STOP {
                    if nodes[ju].stop_val <= i_ndx {
                        j -= 1;
                        continue;
                    }
                    let slot = (nodes[ju].ndx % 3) as usize;
                    if !flag {
                        if nodes[i].star_ptr[slot] == -1 {
                            nodes[i].star_ptr[slot] = j;
                        }
                    } else {
                        let v = nodes[ju].cscore
                            + nodes[ju].sscore
                            + intergenic_mod_idx(nodes, i, ju, tinf);
                        if v > max_sc {
                            nodes[i].star_ptr[slot] = j;
                            max_sc = v;
                        }
                    }
                }
                j -= 1;
            }
        } else {
            let mut max_sc = -100.0f64;
            let mut j = i as i64 - 3;
            while j < nn {
                if j < 0 {
                    j += 1;
                    continue;
                }
                let ju = j as usize;
                if nodes[ju].ndx < i_ndx - 2 {
                    j += 1;
                    continue;
                }
                if nodes[ju].ndx - MAX_SAM_OVLP > i_ndx {
                    break;
                }
                if nodes[ju].strand == -1 && nodes[ju].typ != STOP {
                    if nodes[ju].stop_val >= i_ndx {
                        j += 1;
                        continue;
                    }
                    let slot = (nodes[ju].ndx % 3) as usize;
                    if !flag {
                        if nodes[i].star_ptr[slot] == -1 {
                            nodes[i].star_ptr[slot] = j;
                        }
                    } else {
                        let v = nodes[ju].cscore
                            + nodes[ju].sscore
                            + intergenic_mod_idx(nodes, ju, i, tinf);
                        if v > max_sc {
                            nodes[i].star_ptr[slot] = j;
                            max_sc = v;
                        }
                    }
                }
                j += 1;
            }
        }
    }
}

/// Build the hexamer coding model from the training gene set selected by the
/// initial dynamic-programming pass (`calc_dicodon_gene`). `dbeg` is the index
/// returned by `dprog(flag=false)`.
pub fn calc_dicodon_gene(tinf: &mut Tinf, seq: &Seq, nodes: &[Node], dbeg: i64) {
    let slen = seq.len() as i64;
    let fwd = &seq.fwd;
    let rev = &seq.rev;
    let mut counts = vec![0i64; 4096];
    let bg = crate::seq::calc_mer_bg(6, fwd, rev);
    let mut glob = 0i64;
    let (mut left, mut right);
    left = -1i64;
    right = -1i64;
    let mut in_gene = 0i32;
    let mut path = dbeg;
    while path != -1 {
        let p = path as usize;
        if nodes[p].strand == -1 && nodes[p].typ != STOP {
            in_gene = -1;
            left = slen - nodes[p].ndx - 1;
        }
        if nodes[p].strand == 1 && nodes[p].typ == STOP {
            in_gene = 1;
            right = nodes[p].ndx + 2;
        }
        if in_gene == -1 && nodes[p].strand == -1 && nodes[p].typ == STOP {
            right = slen - nodes[p].ndx + 1;
            let mut i = left;
            while i < right - 5 {
                counts[mer_ndx(6, rev, i as usize)] += 1;
                glob += 1;
                i += 3;
            }
            in_gene = 0;
        }
        if in_gene == 1 && nodes[p].strand == 1 && nodes[p].typ != STOP {
            left = nodes[p].ndx;
            let mut i = left;
            while i < right - 5 {
                counts[mer_ndx(6, fwd, i as usize)] += 1;
                glob += 1;
                i += 3;
            }
            in_gene = 0;
        }
        path = nodes[p].traceb;
    }
    for i in 0..4096usize {
        let prob = if glob > 0 {
            counts[i] as f64 / glob as f64
        } else {
            0.0
        };
        let mut v = if prob == 0.0 && bg[i] != 0.0 {
            -5.0
        } else if bg[i] == 0.0 {
            0.0
        } else {
            (prob / bg[i]).ln()
        };
        if v > 5.0 {
            v = 5.0;
        }
        if v < -5.0 {
            v = -5.0;
        }
        tinf.gene_dc[i] = v;
    }
}

/// GC content for each start->stop pair (`calc_orf_gc`).
fn calc_orf_gc(seq: &Seq, nodes: &mut [Node]) {
    let fwd = &seq.fwd;
    let nn = nodes.len();
    let mut last = [0i64; 3];
    let mut gc = [0.0f64; 3];
    for idx in (0..nn).rev() {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == 1 && nodes[idx].typ == STOP {
            last[fr] = nodes[idx].ndx;
            let n = nodes[idx].ndx;
            gc[fr] = gcv(fwd, n) + gcv(fwd, n + 1) + gcv(fwd, n + 2);
        } else if nodes[idx].strand == 1 {
            let mut j = last[fr] - 3;
            while j >= nodes[idx].ndx {
                gc[fr] += gcv(fwd, j) + gcv(fwd, j + 1) + gcv(fwd, j + 2);
                j -= 3;
            }
            let gsize = (nodes[idx].stop_val - nodes[idx].ndx).abs() as f64 + 3.0;
            nodes[idx].gc_cont = gc[fr] / gsize;
            last[fr] = nodes[idx].ndx;
        }
    }
    gc = [0.0; 3];
    for idx in 0..nn {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == -1 && nodes[idx].typ == STOP {
            last[fr] = nodes[idx].ndx;
            let n = nodes[idx].ndx;
            gc[fr] = gcv(fwd, n) + gcv(fwd, n - 1) + gcv(fwd, n - 2);
        } else if nodes[idx].strand == -1 {
            let mut j = last[fr] + 1;
            while j <= nodes[idx].ndx {
                gc[fr] += gcv(fwd, j) + gcv(fwd, j + 1) + gcv(fwd, j + 2);
                j += 3;
            }
            let gsize = (nodes[idx].stop_val - nodes[idx].ndx).abs() as f64 + 3.0;
            nodes[idx].gc_cont = gc[fr] / gsize;
            last[fr] = nodes[idx].ndx;
        }
    }
}

#[inline]
fn gcv(seq: &[u8], i: i64) -> f64 {
    if i < 0 || i as usize >= seq.len() {
        0.0
    } else if is_gc(seq, i as usize) {
        1.0
    } else {
        0.0
    }
}

/// Raw coding score for every start node (`raw_coding_score`).
pub fn raw_coding_score(seq: &Seq, nodes: &mut [Node], tinf: &Tinf) {
    let slen = seq.len() as i64;
    let fwd = &seq.fwd;
    let rev = &seq.rev;
    let nn = nodes.len();
    let no_stop = if tinf.trans_table != 11 {
        let g = tinf.gc;
        let mut ns = ((1.0 - g) * (1.0 - g) * g) / 8.0;
        ns += ((1.0 - g) * (1.0 - g) * (1.0 - g)) / 8.0;
        1.0 - ns
    } else {
        let g = tinf.gc;
        let mut ns = ((1.0 - g) * (1.0 - g) * g) / 4.0;
        ns += ((1.0 - g) * (1.0 - g) * (1.0 - g)) / 8.0;
        1.0 - ns
    };

    let mut score = [0.0f64; 3];
    let mut last = [0i64; 3];
    for idx in (0..nn).rev() {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == 1 && nodes[idx].typ == STOP {
            last[fr] = nodes[idx].ndx;
            score[fr] = 0.0;
        } else if nodes[idx].strand == 1 {
            let mut j = last[fr] - 3;
            while j >= nodes[idx].ndx {
                score[fr] += tinf.gene_dc[mer_ndx(6, fwd, j as usize)];
                j -= 3;
            }
            nodes[idx].cscore = score[fr];
            last[fr] = nodes[idx].ndx;
        }
    }
    score = [0.0; 3];
    for idx in 0..nn {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == -1 && nodes[idx].typ == STOP {
            last[fr] = nodes[idx].ndx;
            score[fr] = 0.0;
        } else if nodes[idx].strand == -1 {
            let mut j = last[fr] + 3;
            while j <= nodes[idx].ndx {
                score[fr] += tinf.gene_dc[mer_ndx(6, rev, (slen - j - 1) as usize)];
                j += 3;
            }
            nodes[idx].cscore = score[fr];
            last[fr] = nodes[idx].ndx;
        }
    }

    // Second pass: penalize starts with ascending coding to their left.
    let mut score = [-10000.0f64; 3];
    for idx in 0..nn {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == 1 && nodes[idx].typ == STOP {
            score[fr] = -10000.0;
        } else if nodes[idx].strand == 1 {
            if nodes[idx].cscore > score[fr] {
                score[fr] = nodes[idx].cscore;
            } else {
                nodes[idx].cscore -= score[fr] - nodes[idx].cscore;
            }
        }
    }
    score = [-10000.0; 3];
    for idx in (0..nn).rev() {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == -1 && nodes[idx].typ == STOP {
            score[fr] = -10000.0;
        } else if nodes[idx].strand == -1 {
            if nodes[idx].cscore > score[fr] {
                score[fr] = nodes[idx].cscore;
            } else {
                nodes[idx].cscore -= score[fr] - nodes[idx].cscore;
            }
        }
    }

    // Third pass: length-based factor.
    let mut score = [-10000.0f64; 3];
    for idx in 0..nn {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == 1 && nodes[idx].typ == STOP {
            score[fr] = -10000.0;
        } else if nodes[idx].strand == 1 {
            length_factor(&mut nodes[idx], &mut score[fr], no_stop);
        }
    }
    let mut score = [-10000.0f64; 3];
    for idx in (0..nn).rev() {
        let fr = (nodes[idx].ndx % 3) as usize;
        if nodes[idx].strand == -1 && nodes[idx].typ == STOP {
            score[fr] = -10000.0;
        } else if nodes[idx].strand == -1 {
            length_factor(&mut nodes[idx], &mut score[fr], no_stop);
        }
    }
}

fn length_factor(node: &mut Node, score_fr: &mut f64, no_stop: f64) {
    let gsize = ((node.stop_val - node.ndx).abs() as f64 + 3.0) / 3.0;
    let mut lfac;
    if gsize > 1000.0 {
        lfac = ((1.0 - no_stop.powf(1000.0)) / no_stop.powf(1000.0)).ln();
        lfac -= ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
        lfac *= (gsize - 80.0) / 920.0;
    } else {
        lfac = ((1.0 - no_stop.powf(gsize)) / no_stop.powf(gsize)).ln();
        lfac -= ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
    }
    if lfac > *score_fr {
        *score_fr = lfac;
    } else {
        lfac -= dmax(dmin(*score_fr - lfac, lfac), 0.0);
    }
    if lfac > 3.0 && node.cscore < 0.5 * lfac {
        node.cscore = 0.5 * lfac;
    }
    node.cscore += lfac;
}

/// Score the SD RBS bins for every start node (`rbs_score`).
pub fn rbs_score(seq: &Seq, nodes: &mut [Node], tinf: &Tinf) {
    let slen = seq.len() as i64;
    let fwd = &seq.fwd;
    let rev = &seq.rev;
    for idx in 0..nodes.len() {
        if nodes[idx].typ == STOP || nodes[idx].edge {
            continue;
        }
        nodes[idx].rbs[0] = 0;
        nodes[idx].rbs[1] = 0;
        if nodes[idx].strand == 1 {
            let start = nodes[idx].ndx;
            let mut j = start - 20;
            while j <= start - 6 {
                let e = crate::rbs::sd_exact(fwd, j as isize, start as isize, &tinf.rbs_wt);
                let m = crate::rbs::sd_mm(fwd, j as isize, start as isize, &tinf.rbs_wt);
                if e > nodes[idx].rbs[0] {
                    nodes[idx].rbs[0] = e;
                }
                if m > nodes[idx].rbs[1] {
                    nodes[idx].rbs[1] = m;
                }
                j += 1;
            }
        } else {
            let start = slen - 1 - nodes[idx].ndx;
            let mut j = slen - nodes[idx].ndx - 21;
            while j <= slen - nodes[idx].ndx - 7 {
                let e = crate::rbs::sd_exact(rev, j as isize, start as isize, &tinf.rbs_wt);
                let m = crate::rbs::sd_mm(rev, j as isize, start as isize, &tinf.rbs_wt);
                if e > nodes[idx].rbs[0] {
                    nodes[idx].rbs[0] = e;
                }
                if m > nodes[idx].rbs[1] {
                    nodes[idx].rbs[1] = m;
                }
                j += 1;
            }
        }
    }
}

#[inline]
fn pick_max_rb(node: &Node, rbs_wt: &[f64; 28]) -> usize {
    let r0 = node.rbs[0];
    let r1 = node.rbs[1];
    if rbs_wt[r0] > rbs_wt[r1] + 1.0 || r1 == 0 {
        r0
    } else if rbs_wt[r0] < rbs_wt[r1] - 1.0 || r0 == 0 {
        r1
    } else {
        r0.max(r1)
    }
}

/// Train the Shine-Dalgarno start model: RBS-bin weights, ATG/GTG/TTG type
/// weights, and upstream base composition (`train_starts_sd`).
pub fn train_starts_sd(seq: &Seq, nodes: &mut [Node], tinf: &mut Tinf) {
    let fwd = &seq.fwd;
    let rev = &seq.rev;
    let nn = nodes.len();
    let wt = tinf.st_wt;
    tinf.type_wt = [0.0; 3];
    tinf.rbs_wt = [0.0; 28];
    tinf.ups_comp = [[0.0; 4]; 32];

    let mut tbg = [0.0f64; 3];
    for idx in 0..nn {
        if nodes[idx].typ == STOP {
            continue;
        }
        tbg[nodes[idx].typ as usize] += 1.0;
    }
    let sum: f64 = tbg.iter().sum();
    if sum > 0.0 {
        for v in tbg.iter_mut() {
            *v /= sum;
        }
    }

    let mut sthresh = 35.0f64;
    for iter in 0..10 {
        // RBS motif background from current weights.
        let mut rbg = [0.0f64; 28];
        for idx in 0..nn {
            if nodes[idx].typ == STOP || nodes[idx].edge {
                continue;
            }
            let mr = pick_max_rb(&nodes[idx], &tinf.rbs_wt);
            rbg[mr] += 1.0;
        }
        let s: f64 = rbg.iter().sum();
        if s > 0.0 {
            for v in rbg.iter_mut() {
                *v /= s;
            }
        }

        let mut rreal = [0.0f64; 28];
        let mut treal = [0.0f64; 3];

        // Forward strand pass.
        let mut best = [0.0f64; 3];
        let mut bndx = [-1i64; 3];
        let mut rbs = [0usize; 3];
        let mut typ = [0usize; 3];
        for idx in 0..nn {
            if nodes[idx].typ != STOP && nodes[idx].edge {
                continue;
            }
            let fr = (nodes[idx].ndx % 3) as usize;
            if nodes[idx].typ == STOP && nodes[idx].strand == 1 {
                if best[fr] >= sthresh
                    && bndx[fr] >= 0
                    && (nodes[bndx[fr] as usize].ndx % 3) as usize == fr
                {
                    rreal[rbs[fr]] += 1.0;
                    treal[typ[fr]] += 1.0;
                    if iter == 9 {
                        count_upstream_composition(fwd, 1, nodes[bndx[fr] as usize].ndx, tinf);
                    }
                }
                best[fr] = 0.0;
                bndx[fr] = -1;
                rbs[fr] = 0;
                typ[fr] = 0;
            } else if nodes[idx].strand == 1 {
                let mr = pick_max_rb(&nodes[idx], &tinf.rbs_wt);
                let cand =
                    nodes[idx].cscore + wt * tinf.rbs_wt[mr] + wt * tinf.type_wt[nodes[idx].typ as usize];
                if cand >= best[fr] {
                    best[fr] = nodes[idx].cscore + wt * tinf.rbs_wt[mr];
                    best[fr] += wt * tinf.type_wt[nodes[idx].typ as usize];
                    bndx[fr] = idx as i64;
                    typ[fr] = nodes[idx].typ as usize;
                    rbs[fr] = mr;
                }
            }
        }

        // Reverse strand pass.
        best = [0.0; 3];
        bndx = [-1; 3];
        rbs = [0; 3];
        typ = [0; 3];
        for idx in (0..nn).rev() {
            if nodes[idx].typ != STOP && nodes[idx].edge {
                continue;
            }
            let fr = (nodes[idx].ndx % 3) as usize;
            if nodes[idx].typ == STOP && nodes[idx].strand == -1 {
                if best[fr] >= sthresh
                    && bndx[fr] >= 0
                    && (nodes[bndx[fr] as usize].ndx % 3) as usize == fr
                {
                    rreal[rbs[fr]] += 1.0;
                    treal[typ[fr]] += 1.0;
                    if iter == 9 {
                        count_upstream_composition(rev, -1, nodes[bndx[fr] as usize].ndx, tinf);
                    }
                }
                best[fr] = 0.0;
                bndx[fr] = -1;
                rbs[fr] = 0;
                typ[fr] = 0;
            } else if nodes[idx].strand == -1 {
                let mr = pick_max_rb(&nodes[idx], &tinf.rbs_wt);
                let cand =
                    nodes[idx].cscore + wt * tinf.rbs_wt[mr] + wt * tinf.type_wt[nodes[idx].typ as usize];
                if cand >= best[fr] {
                    best[fr] = nodes[idx].cscore + wt * tinf.rbs_wt[mr];
                    best[fr] += wt * tinf.type_wt[nodes[idx].typ as usize];
                    bndx[fr] = idx as i64;
                    typ[fr] = nodes[idx].typ as usize;
                    rbs[fr] = mr;
                }
            }
        }

        // Convert RBS counts to log weights.
        let s: f64 = rreal.iter().sum();
        if s == 0.0 {
            tinf.rbs_wt = [0.0; 28];
        } else {
            for j in 0..28 {
                rreal[j] /= s;
                tinf.rbs_wt[j] = if rbg[j] != 0.0 {
                    (rreal[j] / rbg[j]).ln()
                } else {
                    -4.0
                };
                if tinf.rbs_wt[j] > 4.0 {
                    tinf.rbs_wt[j] = 4.0;
                }
                if tinf.rbs_wt[j] < -4.0 {
                    tinf.rbs_wt[j] = -4.0;
                }
            }
        }
        let s: f64 = treal.iter().sum();
        if s == 0.0 {
            tinf.type_wt = [0.0; 3];
        } else {
            for j in 0..3 {
                treal[j] /= s;
                tinf.type_wt[j] = if tbg[j] != 0.0 {
                    (treal[j] / tbg[j]).ln()
                } else {
                    -4.0
                };
                if tinf.type_wt[j] > 4.0 {
                    tinf.type_wt[j] = 4.0;
                }
                if tinf.type_wt[j] < -4.0 {
                    tinf.type_wt[j] = -4.0;
                }
            }
        }
        if s <= nn as f64 / 2000.0 {
            sthresh /= 2.0;
        }
    }

    // Convert upstream base composition to log scores.
    let gc = tinf.gc;
    for i in 0..32 {
        let s: f64 = tinf.ups_comp[i].iter().sum();
        if s == 0.0 {
            tinf.ups_comp[i] = [0.0; 4];
            continue;
        }
        for j in 0..4 {
            tinf.ups_comp[i][j] /= s;
            let p = tinf.ups_comp[i][j];
            tinf.ups_comp[i][j] = if gc > 0.1 && gc < 0.9 {
                if j == 0 || j == 3 {
                    (p * 2.0 / (1.0 - gc)).ln()
                } else {
                    (p * 2.0 / gc).ln()
                }
            } else if gc <= 0.1 {
                if j == 0 || j == 3 {
                    (p * 2.0 / 0.90).ln()
                } else {
                    (p * 2.0 / 0.10).ln()
                }
            } else if j == 0 || j == 3 {
                (p * 2.0 / 0.10).ln()
            } else {
                (p * 2.0 / 0.90).ln()
            };
            if tinf.ups_comp[i][j] > 4.0 {
                tinf.ups_comp[i][j] = 4.0;
            }
            if tinf.ups_comp[i][j] < -4.0 {
                tinf.ups_comp[i][j] = -4.0;
            }
        }
    }
}

/// Count upstream base composition for a putative real start
/// (`count_upstream_composition`).
fn count_upstream_composition(seq: &[u8], strand: i32, pos: i64, tinf: &mut Tinf) {
    let start = if strand == 1 {
        pos
    } else {
        seq.len() as i64 - 1 - pos
    };
    let mut count = 0usize;
    let mut i = 1i64;
    while i < 45 {
        if i > 2 && i < 15 {
            i += 1;
            continue;
        }
        if start - i >= 0 {
            let b = mer_ndx(1, seq, (start - i) as usize);
            tinf.ups_comp[count][b] += 1.0;
        }
        count += 1;
        i += 1;
    }
}

/// Score the upstream base composition of a start (`score_upstream_composition`).
fn score_upstream_composition(seq: &[u8], node: &mut Node, tinf: &Tinf) {
    let start = if node.strand == 1 {
        node.ndx
    } else {
        seq.len() as i64 - 1 - node.ndx
    };
    node.uscore = 0.0;
    let mut count = 0usize;
    let mut i = 1i64;
    while i < 45 {
        if i > 2 && i < 15 {
            i += 1;
            continue;
        }
        if start - i < 0 {
            i += 1;
            continue;
        }
        let b = mer_ndx(1, seq, (start - i) as usize);
        node.uscore += 0.4 * tinf.st_wt * tinf.ups_comp[count][b];
        count += 1;
        i += 1;
    }
}

/// Score every start node: coding + start signal (`score_nodes`).
pub fn score_nodes(seq: &Seq, nodes: &mut [Node], tinf: &Tinf, closed: bool, is_meta: bool) {
    score_nodes_ends(seq, nodes, tinf, closed, closed, is_meta)
}

/// [`score_nodes`] with each sequence end open or closed on its own (see
/// [`add_nodes_ends`]).
pub fn score_nodes_ends(
    seq: &Seq,
    nodes: &mut [Node],
    tinf: &Tinf,
    closed_left: bool,
    closed_right: bool,
    is_meta: bool,
) {
    let slen = seq.len() as i64;
    let fwd = &seq.fwd;
    let rev = &seq.rev;
    let nn = nodes.len();
    let code = tinf.code();
    let famb: &[bool] = &seq.n;
    let ramb_v = seq.n_rev();
    let ramb: &[bool] = &ramb_v;

    calc_orf_gc(seq, nodes);
    raw_coding_score(seq, nodes, tinf);
    rbs_score(seq, nodes, tinf);

    for i in 0..nn {
        if nodes[i].typ == STOP {
            continue;
        }
        let mut edge_gene = 0i32;
        if nodes[i].edge {
            edge_gene += 1;
        }
        let runs_off = (nodes[i].strand == 1 && !is_stop_at(fwd, famb, nodes[i].stop_val, &code))
            || (nodes[i].strand == -1 && !is_stop_at(rev, ramb, slen - 1 - nodes[i].stop_val, &code));
        if runs_off {
            edge_gene += 1;
        }

        if nodes[i].edge {
            nodes[i].tscore = EDGE_BONUS * tinf.st_wt / edge_gene as f64;
            nodes[i].uscore = 0.0;
            nodes[i].rscore = 0.0;
        } else {
            nodes[i].tscore = tinf.type_wt[nodes[i].typ as usize] * tinf.st_wt;
            let rbs1 = tinf.rbs_wt[nodes[i].rbs[0]];
            let rbs2 = tinf.rbs_wt[nodes[i].rbs[1]];
            nodes[i].rscore = dmax(rbs1, rbs2) * tinf.st_wt;

            if nodes[i].strand == 1 {
                score_upstream_composition(fwd, &mut nodes[i], tinf);
            } else {
                score_upstream_composition(rev, &mut nodes[i], tinf);
            }

            if !closed_left && nodes[i].ndx <= 2 && nodes[i].strand == 1 {
                nodes[i].uscore += EDGE_UPS * tinf.st_wt;
            } else if !closed_right && nodes[i].ndx >= slen - 3 && nodes[i].strand == -1 {
                nodes[i].uscore += EDGE_UPS * tinf.st_wt;
            } else if i < 500 && nodes[i].strand == 1 {
                for j in (0..i).rev() {
                    if nodes[j].edge && nodes[i].stop_val == nodes[j].stop_val {
                        nodes[i].uscore += EDGE_UPS * tinf.st_wt;
                        break;
                    }
                }
            } else if i + 500 >= nn && nodes[i].strand == -1 {
                for j in (i + 1)..nn {
                    if nodes[j].edge && nodes[i].stop_val == nodes[j].stop_val {
                        nodes[i].uscore += EDGE_UPS * tinf.st_wt;
                        break;
                    }
                }
            }
        }

        // Convert near-edge starts to edge genes when ends are open.
        if ((nodes[i].ndx <= 2 && nodes[i].strand == 1 && !closed_left)
            || (nodes[i].ndx >= slen - 3 && nodes[i].strand == -1 && !closed_right))
            && !nodes[i].edge
        {
            edge_gene += 1;
            nodes[i].edge = true;
            nodes[i].tscore = 0.0;
            nodes[i].uscore = EDGE_BONUS * tinf.st_wt / edge_gene as f64;
            nodes[i].rscore = 0.0;
        }

        if !nodes[i].edge && edge_gene == 1 {
            nodes[i].uscore -= 0.5 * EDGE_BONUS * tinf.st_wt;
        }

        // Penalize non-edge genes < 250 bp.
        let glen = (nodes[i].ndx - nodes[i].stop_val).abs();
        if edge_gene == 0 && glen < 250 {
            let negf = 250.0 / glen as f64;
            let posf = glen as f64 / 250.0;
            scale_scores(&mut nodes[i], negf, posf);
        }

        // Metagenomic short-fragment coding penalty.
        if is_meta && slen < 3000 && edge_gene == 0 && (nodes[i].cscore < 5.0 || glen < 120) {
            nodes[i].cscore -= META_PEN * dmax(0.0, (3000.0 - slen as f64) / 2700.0);
        }

        nodes[i].sscore = nodes[i].tscore + nodes[i].rscore + nodes[i].uscore;

        if nodes[i].cscore < 0.0 {
            if edge_gene > 0 && !nodes[i].edge {
                if !is_meta || slen > 1500 {
                    nodes[i].sscore -= tinf.st_wt;
                } else {
                    nodes[i].sscore -= 10.31 - 0.004 * slen as f64;
                }
            } else if is_meta && slen < 3000 && nodes[i].edge {
                let min_meta_len = (slen as f64).sqrt() * 5.0;
                if glen as f64 >= min_meta_len {
                    if nodes[i].cscore >= 0.0 {
                        nodes[i].cscore = -1.0;
                    }
                    nodes[i].sscore = 0.0;
                    nodes[i].uscore = 0.0;
                }
            } else {
                nodes[i].sscore -= 0.5;
            }
        } else if nodes[i].cscore < 5.0 && is_meta && glen < 120 && nodes[i].sscore < 0.0 {
            nodes[i].sscore -= tinf.st_wt;
        }
    }
}

fn scale_scores(node: &mut Node, negf: f64, posf: f64) {
    if node.rscore < 0.0 {
        node.rscore *= negf;
    }
    if node.uscore < 0.0 {
        node.uscore *= negf;
    }
    if node.tscore < 0.0 {
        node.tscore *= negf;
    }
    if node.rscore > 0.0 {
        node.rscore *= posf;
    }
    if node.uscore > 0.0 {
        node.uscore *= posf;
    }
    if node.tscore > 0.0 {
        node.tscore *= posf;
    }
}

/// Intergenic modifier between two nodes (`intergenic_mod`), index form.
pub fn intergenic_mod_idx(nodes: &[Node], i1: usize, i2: usize, tinf: &Tinf) -> f64 {
    intergenic_mod(&nodes[i1], &nodes[i2], tinf)
}

/// Intergenic modifier between two nodes (`intergenic_mod`).
pub fn intergenic_mod(n1: &Node, n2: &Node, tinf: &Tinf) -> f64 {
    let mut rval = 0.0f64;
    let mut ovlp = 0.0f64;
    let same_close = (n1.strand == 1
        && n2.strand == 1
        && (n1.ndx + 2 == n2.ndx || n1.ndx - 1 == n2.ndx))
        || (n1.strand == -1 && n2.strand == -1 && (n1.ndx + 2 == n2.ndx || n1.ndx - 1 == n2.ndx));
    if same_close {
        if n1.strand == 1 && n2.rscore < 0.0 {
            rval -= n2.rscore;
        }
        if n1.strand == -1 && n1.rscore < 0.0 {
            rval -= n1.rscore;
        }
        if n1.strand == 1 && n2.uscore < 0.0 {
            rval -= n2.uscore;
        }
        if n1.strand == -1 && n1.uscore < 0.0 {
            rval -= n1.uscore;
        }
    }
    let dist = (n1.ndx - n2.ndx).abs();
    if n1.strand == 1 && n2.strand == 1 && n1.ndx + 2 >= n2.ndx {
        ovlp = 1.0;
    } else if n1.strand == -1 && n2.strand == -1 && n1.ndx >= n2.ndx + 2 {
        ovlp = 1.0;
    }
    if dist as f64 > 3.0 * OPER_DIST || n1.strand != n2.strand {
        rval -= 0.15 * tinf.st_wt;
    } else if (dist as f64 <= OPER_DIST && ovlp == 0.0) || (dist as f64) < 0.25 * OPER_DIST {
        rval += (2.0 - dist as f64 / OPER_DIST) * 0.15 * tinf.st_wt;
    }
    rval
}

/// Whether the trained model uses an SD motif (Prodigal `determine_sd_usage`).
/// Reported for the user; DYNAMITE always scores via the SD path.
pub fn determine_sd_usage(tinf: &mut Tinf) {
    tinf.uses_sd = true;
    if tinf.rbs_wt[0] >= 0.0 {
        tinf.uses_sd = false;
    }
    if tinf.rbs_wt[16] < 1.0
        && tinf.rbs_wt[13] < 1.0
        && tinf.rbs_wt[15] < 1.0
        && (tinf.rbs_wt[0] >= -0.5
            || (tinf.rbs_wt[22] < 2.0 && tinf.rbs_wt[24] < 2.0 && tinf.rbs_wt[27] < 2.0))
    {
        tinf.uses_sd = false;
    }
}

/// Convenience: the per-base most-GC frame plot (trailing bases reported as 0;
/// the Prodigal trainer itself uses [`crate::gcframe::most_gc_frame_prodigal`]).
pub fn gc_frame_plot(seq: &Seq) -> Vec<u8> {
    most_gc_frame(&seq.fwd)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mask_crossing_is_half_open() {
        let masks = [(10usize, 20usize), (40, 50)];
        assert!(!crosses_mask(&masks, 0, 10)); // ends right before the mask
        assert!(crosses_mask(&masks, 0, 11));
        assert!(crosses_mask(&masks, 19, 30));
        assert!(!crosses_mask(&masks, 20, 40)); // between the masks
        assert!(crosses_mask(&masks, 20, 41));
        assert!(crosses_mask(&masks, 0, 100)); // spans both
        assert!(!crosses_mask(&[], 0, 100));
    }

    #[test]
    fn masked_and_one_sided_nodes_reduce_to_add_nodes() {
        let s = Seq::from_ascii(b"ATGAAACCCGGGTTTAAATGACCCATGGGTTTCCCAAATAGGCATGCATGCATGTAAATGCCC");
        let code = GeneticCode::new(11).unwrap();
        for closed in [false, true] {
            let a = add_nodes(&s, &code, closed);
            let b = add_nodes_ends(&s, &code, closed, closed, &[]);
            assert_eq!(a.len(), b.len());
            assert!(a.iter().zip(&b).all(|(x, y)| x.ndx == y.ndx && x.stop_val == y.stop_val && x.edge == y.edge));
        }
        // a mask over the whole sequence leaves no start node
        let all = add_nodes_masked(&s, &code, false, &[(0, s.len())]);
        assert!(all.iter().all(|n| n.typ == STOP));
    }
}
