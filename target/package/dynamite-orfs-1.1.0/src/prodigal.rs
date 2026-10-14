//! The Prodigal-style dynamic-programming engine: the DP traceback
//! (`dprog` / `score_connection`), gene extraction (`add_genes`,
//! `tweak_final_starts`), and the high-level self-training / prediction
//! orchestration including DYNAMITE's genetic-code auto-detection.
//!
//! Faithful Rust port of `dprog.c` and parts of `gene.c` (Prodigal /
//! Prodigal-GV, Doug Hyatt et al., GPL-3.0).

use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::node::{
    add_nodes_ends, add_nodes_masked, calc_dicodon_gene, determine_sd_usage, intergenic_mod,
    raw_coding_score, rbs_score, record_gc_bias, record_overlapping_starts, score_nodes_ends,
    sort_nodes, train_starts_sd, Node, Tinf, STOP,
};
use crate::seq::Seq;
use rayon::prelude::*;

const MAX_OPP_OVLP: i64 = 200;
const MAX_NODE_DIST: isize = 500;

/// The training concatenation separator. `TTAATTAATTAA` forces stop codons in
/// all six reading frames, preventing ORFs from spanning contig junctions
/// (exactly as Prodigal's `read_seq_training` does).
pub const TRAIN_GAP: &[u8] = b"TTAATTAATTAA";

/// Prodigal's Shine-Dalgarno motif / spacer labels for the 28 RBS bins
/// (`record_gene_data`, `sd_string` / `sd_spacer`).
const SD_STRING: [&str; 28] = [
    "None", "GGA/GAG/AGG", "3Base/5BMM", "4Base/6BMM", "AGxAG", "AGxAG", "GGA/GAG/AGG",
    "GGxGG", "GGxGG", "AGxAG", "AGGAG(G)/GGAGG", "AGGA/GGAG/GAGG", "AGGA/GGAG/GAGG",
    "GGA/GAG/AGG", "GGxGG", "AGGA", "GGAG/GAGG", "AGxAGG/AGGxGG", "AGxAGG/AGGxGG",
    "AGxAGG/AGGxGG", "AGGAG/GGAGG", "AGGAG", "AGGAG", "GGAGG", "GGAGG", "AGGAGG", "AGGAGG",
    "AGGAGG",
];
const SD_SPACER: [&str; 28] = [
    "None", "3-4bp", "13-15bp", "13-15bp", "11-12bp", "3-4bp", "11-12bp", "11-12bp", "3-4bp",
    "5-10bp", "13-15bp", "3-4bp", "11-12bp", "5-10bp", "5-10bp", "5-10bp", "5-10bp", "11-12bp",
    "3-4bp", "5-10bp", "11-12bp", "3-4bp", "5-10bp", "3-4bp", "5-10bp", "11-12bp", "3-4bp",
    "5-10bp",
];

/// The `rbs_motif` / `rbs_spacer` pair Prodigal reports for a start node.
/// DYNAMITE implements only the Shine-Dalgarno trainer, so organisms that
/// Prodigal classes as non-SD (where it would report a de-novo motif) get
/// `None` here.
fn rbs_strings(node: &Node, tinf: &Tinf) -> (&'static str, &'static str) {
    if !tinf.uses_sd {
        return ("None", "None");
    }
    let rbs1 = tinf.rbs_wt[node.rbs[0]] * tinf.st_wt;
    let rbs2 = tinf.rbs_wt[node.rbs[1]] * tinf.st_wt;
    let bin = if rbs1 > rbs2 { node.rbs[0] } else { node.rbs[1] };
    (SD_STRING[bin], SD_SPACER[bin])
}

/// Confidence from a model score (`calculate_confidence`).
pub fn calculate_confidence(score: f64, start_weight: f64) -> f64 {
    let mut conf;
    if score / start_weight < 41.0 {
        let c = (score / start_weight).exp();
        conf = (c / (c + 1.0)) * 100.0;
    } else {
        conf = 99.99;
    }
    if conf <= 50.0 {
        conf = 50.0;
    }
    conf
}

/// The core dynamic-programming routine (`dprog`). Returns the index of the
/// maximum-scoring terminal node, or `-1`.
pub fn dprog(nodes: &mut [Node], tinf: &Tinf, flag: bool) -> i64 {
    let nn = nodes.len();
    if nn == 0 {
        return -1;
    }
    for n in nodes.iter_mut() {
        n.score = 0.0;
        n.traceb = -1;
        n.tracef = -1;
    }

    for i in 0..nn {
        let i_i = i as isize;
        let mut min: isize = if i_i < MAX_NODE_DIST {
            0
        } else {
            i_i - MAX_NODE_DIST
        };
        // Giant-ORF exception (Prodigal `dprog`): walk `min` back to the node
        // that sits at this ORF's other end, then allow MAX_NODE_DIST more.
        let guard = (nodes[i].strand == -1 && nodes[i].typ != STOP)
            || (nodes[i].strand == 1 && nodes[i].typ == STOP);
        if guard && nodes[min as usize].ndx >= nodes[i].stop_val {
            while min >= 0 && nodes[min as usize].ndx != nodes[i].stop_val {
                min -= 1;
            }
        }
        if min < MAX_NODE_DIST {
            min = 0;
        } else {
            min -= MAX_NODE_DIST;
        }
        let mut j = min as usize;
        while j < i {
            score_connection(nodes, j, i, tinf, flag);
            j += 1;
        }
    }

    let mut max_sc = -1.0f64;
    let mut max_ndx: i64 = -1;
    for i in (0..nn).rev() {
        if nodes[i].strand == 1 && nodes[i].typ != STOP {
            continue;
        }
        if nodes[i].strand == -1 && nodes[i].typ == STOP {
            continue;
        }
        if nodes[i].score > max_sc {
            max_sc = nodes[i].score;
            max_ndx = i as i64;
        }
    }
    if max_ndx == -1 {
        return -1;
    }

    // First pass: untangle triple overlaps.
    let mut path = max_ndx;
    while nodes[path as usize].traceb != -1 {
        let nxt = nodes[path as usize].traceb;
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ == STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].typ == STOP
            && nodes[path as usize].ov_mark != -1
            && nodes[path as usize].ndx > nodes[nxt as usize].ndx
        {
            let ov = nodes[path as usize].ov_mark as usize;
            let tmp = nodes[path as usize].star_ptr[ov];
            let mut i = tmp;
            while nodes[i as usize].ndx != nodes[tmp as usize].stop_val {
                i -= 1;
            }
            nodes[path as usize].traceb = tmp;
            nodes[tmp as usize].traceb = i;
            nodes[i as usize].ov_mark = -1;
            nodes[i as usize].traceb = nxt;
        }
        path = nodes[path as usize].traceb;
    }

    // Second pass: untangle simple overlaps.
    let mut path = max_ndx;
    while nodes[path as usize].traceb != -1 {
        let nxt = nodes[path as usize].traceb;
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ != STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].typ == STOP
        {
            let mut i = path;
            while nodes[i as usize].ndx != nodes[path as usize].stop_val {
                i -= 1;
            }
            nodes[path as usize].traceb = i;
            nodes[i as usize].traceb = nxt;
        }
        if nodes[path as usize].strand == 1
            && nodes[path as usize].typ == STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].typ == STOP
        {
            let slot = (nodes[path as usize].ndx % 3) as usize;
            let tb = nodes[nxt as usize].star_ptr[slot];
            nodes[path as usize].traceb = tb;
            if tb >= 0 {
                nodes[tb as usize].traceb = nxt;
            }
        }
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ == STOP
            && nodes[nxt as usize].strand == -1
            && nodes[nxt as usize].typ == STOP
        {
            let slot = (nodes[nxt as usize].ndx % 3) as usize;
            let tb = nodes[path as usize].star_ptr[slot];
            nodes[path as usize].traceb = tb;
            if tb >= 0 {
                nodes[tb as usize].traceb = nxt;
            }
        }
        path = nodes[path as usize].traceb;
        if path < 0 {
            break;
        }
    }

    // Mark forward pointers.
    let mut path = max_ndx;
    while nodes[path as usize].traceb != -1 {
        let tb = nodes[path as usize].traceb;
        nodes[tb as usize].tracef = path;
        path = tb;
    }

    if nodes[max_ndx as usize].traceb == -1 {
        -1
    } else {
        max_ndx
    }
}

/// Score the connection between two nodes (`score_connection`), index form.
fn score_connection(nodes: &mut [Node], p1: usize, p2: usize, tinf: &Tinf, flag: bool) {
    // Copy the scalar fields we need (all Copy).
    let n1_type = nodes[p1].typ;
    let n1_strand = nodes[p1].strand;
    let n1_ndx = nodes[p1].ndx;
    let n1_traceb = nodes[p1].traceb;
    let n1_score = nodes[p1].score;
    let n1_cscore = nodes[p1].cscore;
    let n1_sscore = nodes[p1].sscore;
    let n1_gc = nodes[p1].gc_score;

    let n2_type = nodes[p2].typ;
    let n2_strand = nodes[p2].strand;
    let n2_ndx = nodes[p2].ndx;
    let n2_stop_val = nodes[p2].stop_val;
    let n2_score = nodes[p2].score;
    let n2_cscore = nodes[p2].cscore;
    let n2_sscore = nodes[p2].sscore;
    let n2_gc = nodes[p2].gc_score;

    let mut left = n1_ndx;
    let mut right = n2_ndx;
    let mut ovlp = 0i64;
    let mut maxfr: i64 = -1;
    let mut score = 0.0f64;
    let mut scr_mod = 0.0f64;

    // ---- Invalid connections ----
    if n1_type != STOP && n2_type != STOP && n1_strand == n2_strand {
        return;
    } else if n1_strand == 1 && n1_type != STOP && n2_strand == -1 {
        return;
    } else if n1_strand == -1 && n1_type == STOP && n2_strand == 1 {
        return;
    } else if n1_strand == -1 && n1_type != STOP && n2_strand == 1 && n2_type == STOP {
        return;
    }

    // ---- Edge artifacts ----
    if n1_traceb == -1 && n1_strand == 1 && n1_type == STOP {
        return;
    }
    if n1_traceb == -1 && n1_strand == -1 && n1_type != STOP {
        return;
    }

    let bias = tinf.bias;
    let gc_dot = |g: &[f64; 3]| bias[0] * g[0] + bias[1] * g[1] + bias[2] * g[2];

    // ---- Genes ----
    if n1_strand == n2_strand && n1_strand == 1 && n1_type != STOP && n2_type == STOP {
        if n2_stop_val >= n1_ndx {
            return;
        }
        if n1_ndx % 3 != n2_ndx % 3 {
            return;
        }
        right += 2;
        if !flag {
            scr_mod = gc_dot(&n1_gc);
        } else {
            score = n1_cscore + n1_sscore;
        }
    } else if n1_strand == n2_strand && n1_strand == -1 && n1_type == STOP && n2_type != STOP {
        if nodes[p1].stop_val <= n2_ndx {
            return;
        }
        if n1_ndx % 3 != n2_ndx % 3 {
            return;
        }
        left -= 2;
        if !flag {
            scr_mod = gc_dot(&n2_gc);
        } else {
            score = n2_cscore + n2_sscore;
        }
    }
    // ---- Intergenic ----
    else if n1_strand == 1 && n1_type == STOP && n2_strand == 1 && n2_type != STOP {
        left += 2;
        if left >= right {
            return;
        }
        if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    } else if n1_strand == 1 && n1_type == STOP && n2_strand == -1 && n2_type == STOP {
        left += 2;
        right -= 2;
        if left >= right {
            return;
        }
        maxfr = -1;
        let mut maxval = 0.0f64;
        for i in 0..3 {
            let sp = nodes[p2].star_ptr[i];
            if sp == -1 {
                continue;
            }
            let k = sp as usize;
            let n3_stop_val = nodes[k].stop_val;
            let n3_ndx = nodes[k].ndx;
            // Prodigal reuses the function-level `ovlp` here, so the value from
            // the last examined start survives the loop and feeds the flag=0
            // (GC-frame training) score below. Reproduce that exactly.
            ovlp = left - n3_stop_val + 3;
            let cur_ovlp = ovlp;
            if cur_ovlp <= 0 || cur_ovlp >= MAX_OPP_OVLP {
                continue;
            }
            if cur_ovlp >= n3_ndx - left {
                continue;
            }
            if n1_traceb == -1 {
                continue;
            }
            if cur_ovlp >= n3_stop_val - nodes[n1_traceb as usize].ndx - 2 {
                continue;
            }
            let val = if flag {
                nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf)
            } else {
                gc_dot(&nodes[k].gc_score)
            };
            if val > maxval {
                maxfr = i as i64;
                maxval = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf);
            }
        }
        if maxfr != -1 {
            let k = nodes[p2].star_ptr[maxfr as usize] as usize;
            if !flag {
                scr_mod = gc_dot(&nodes[k].gc_score);
            } else {
                score = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf);
            }
        } else if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    } else if n1_strand == -1 && n1_type != STOP && n2_strand == -1 && n2_type == STOP {
        right -= 2;
        if left >= right {
            return;
        }
        if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    } else if n1_strand == -1 && n1_type != STOP && n2_strand == 1 && n2_type != STOP {
        if left >= right {
            return;
        }
        if flag {
            score = intergenic_mod(&nodes[p1], &nodes[p2], tinf);
        }
    }
    // ---- Possible operons ----
    else if n1_strand == 1 && n2_strand == 1 && n1_type == STOP && n2_type == STOP {
        if n2_stop_val >= n1_ndx {
            return;
        }
        let slot = (n2_ndx % 3) as usize;
        if nodes[p1].star_ptr[slot] == -1 {
            return;
        }
        let k = nodes[p1].star_ptr[slot] as usize;
        left = nodes[k].ndx;
        right += 2;
        if !flag {
            scr_mod = gc_dot(&nodes[k].gc_score);
        } else {
            score = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[p1], &nodes[k], tinf);
        }
    } else if n1_strand == -1 && n1_type == STOP && n2_strand == -1 && n2_type == STOP {
        if nodes[p1].stop_val <= n2_ndx {
            return;
        }
        let slot = (n1_ndx % 3) as usize;
        if nodes[p2].star_ptr[slot] == -1 {
            return;
        }
        let k = nodes[p2].star_ptr[slot] as usize;
        left -= 2;
        right = nodes[k].ndx;
        if !flag {
            scr_mod = gc_dot(&nodes[k].gc_score);
        } else {
            score = nodes[k].cscore + nodes[k].sscore + intergenic_mod(&nodes[k], &nodes[p2], tinf);
        }
    }
    // ---- Overlapping opposite-strand 3' ends ----
    else if n1_strand == 1 && n1_type == STOP && n2_strand == -1 && n2_type != STOP {
        if n2_stop_val - 2 >= n1_ndx + 2 {
            return;
        }
        ovlp = (n1_ndx + 2) - (n2_stop_val - 2) + 1;
        if ovlp >= MAX_OPP_OVLP {
            return;
        }
        if (n1_ndx + 2 - n2_stop_val - 2 + 1) >= (n2_ndx - n1_ndx + 3 + 1) {
            return;
        }
        let bnd = if n1_traceb == -1 {
            0
        } else {
            nodes[n1_traceb as usize].ndx
        };
        if (n1_ndx + 2 - n2_stop_val - 2 + 1) >= (n2_stop_val - 3 - bnd + 1) {
            return;
        }
        left = n2_stop_val - 2;
        if !flag {
            scr_mod = gc_dot(&n2_gc);
        } else {
            score = n2_cscore + n2_sscore - 0.15 * tinf.st_wt;
        }
    }

    if !flag {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }

    if n1_score + score >= n2_score {
        nodes[p2].score = n1_score + score;
        nodes[p2].traceb = p1 as i64;
        nodes[p2].ov_mark = maxfr;
    }
}

/// Eliminate genes with negative scores (`eliminate_bad_genes`).
fn eliminate_bad_genes(nodes: &mut [Node], dbeg: i64, tinf: &Tinf) {
    if dbeg == -1 {
        return;
    }
    // Walk to the head.
    let mut head = dbeg;
    while nodes[head as usize].traceb != -1 {
        head = nodes[head as usize].traceb;
    }

    // Apply intergenic modifiers along the chain.
    let mut path = head;
    while nodes[path as usize].tracef != -1 {
        let tf = nodes[path as usize].tracef;
        if nodes[path as usize].strand == 1 && nodes[path as usize].typ == STOP {
            let m = intergenic_mod(&nodes[path as usize], &nodes[tf as usize], tinf);
            nodes[tf as usize].sscore += m;
        }
        if nodes[path as usize].strand == -1 && nodes[path as usize].typ != STOP {
            let m = intergenic_mod(&nodes[path as usize], &nodes[tf as usize], tinf);
            nodes[path as usize].sscore += m;
        }
        path = tf;
    }

    let mut path = head;
    while nodes[path as usize].tracef != -1 {
        let tf = nodes[path as usize].tracef;
        if nodes[path as usize].strand == 1
            && nodes[path as usize].typ != STOP
            && nodes[path as usize].cscore + nodes[path as usize].sscore < 0.0
        {
            nodes[path as usize].elim = true;
            nodes[tf as usize].elim = true;
        }
        if nodes[path as usize].strand == -1
            && nodes[path as usize].typ == STOP
            && nodes[tf as usize].cscore + nodes[tf as usize].sscore < 0.0
        {
            nodes[path as usize].elim = true;
            nodes[tf as usize].elim = true;
        }
        path = tf;
    }
}

/// Minimal gene record used between `add_genes` and `tweak_final_starts`.
#[derive(Clone, Copy)]
struct RawGene {
    begin: i64,
    end: i64,
    start_ndx: i64,
    stop_ndx: i64,
}

/// Extract the gene list from the traceback (`add_genes`).
fn add_genes(nodes: &[Node], dbeg: i64) -> Vec<RawGene> {
    let mut glist = Vec::new();
    if dbeg == -1 {
        return glist;
    }
    let mut path = dbeg;
    while nodes[path as usize].traceb != -1 {
        path = nodes[path as usize].traceb;
    }
    let mut cur = RawGene {
        begin: 0,
        end: 0,
        start_ndx: -1,
        stop_ndx: -1,
    };
    let mut have = false;
    let mut p = path;
    while p != -1 {
        let node = &nodes[p as usize];
        if node.elim {
            p = node.tracef;
            continue;
        }
        if node.strand == 1 && node.typ != STOP {
            cur.begin = node.ndx + 1;
            cur.start_ndx = p;
            have = true;
        }
        if node.strand == -1 && node.typ == STOP {
            cur.begin = node.ndx - 1;
            cur.stop_ndx = p;
            have = true;
        }
        if node.strand == 1 && node.typ == STOP {
            cur.end = node.ndx + 3;
            cur.stop_ndx = p;
            glist.push(cur);
            cur = RawGene {
                begin: 0,
                end: 0,
                start_ndx: -1,
                stop_ndx: -1,
            };
            have = false;
        }
        if node.strand == -1 && node.typ != STOP {
            cur.end = node.ndx + 1;
            cur.start_ndx = p;
            glist.push(cur);
            cur = RawGene {
                begin: 0,
                end: 0,
                start_ndx: -1,
                stop_ndx: -1,
            };
            have = false;
        }
        p = node.tracef;
    }
    let _ = have;
    glist
}

/// Refine close / rare-codon starts (`tweak_final_starts`).
fn tweak_final_starts(genes: &mut [RawGene], nodes: &[Node], tinf: &Tinf) {
    let ng = genes.len();
    let nn = nodes.len() as i64;
    for i in 0..ng {
        let ndx = genes[i].start_ndx;
        if ndx < 0 {
            continue;
        }
        let ndxu = ndx as usize;
        let sc = nodes[ndxu].sscore + nodes[ndxu].cscore;
        let mut igm = 0.0f64;
        if i > 0 {
            let prev_start = genes[i - 1].start_ndx;
            let prev_stop = genes[i - 1].stop_ndx;
            if prev_start >= 0 {
                let ps = nodes[prev_start as usize].strand;
                if nodes[ndxu].strand == 1 && ps == 1 && prev_stop >= 0 {
                    igm = intergenic_mod(&nodes[prev_stop as usize], &nodes[ndxu], tinf);
                } else if nodes[ndxu].strand == 1 && ps == -1 {
                    igm = intergenic_mod(&nodes[prev_start as usize], &nodes[ndxu], tinf);
                }
            }
        }
        if i + 1 < ng {
            let next_start = genes[i + 1].start_ndx;
            let next_stop = genes[i + 1].stop_ndx;
            if next_start >= 0 {
                let ns = nodes[next_start as usize].strand;
                if nodes[ndxu].strand == -1 && ns == 1 {
                    igm = intergenic_mod(&nodes[ndxu], &nodes[next_start as usize], tinf);
                } else if nodes[ndxu].strand == -1 && ns == -1 && next_stop >= 0 {
                    igm = intergenic_mod(&nodes[ndxu], &nodes[next_stop as usize], tinf);
                }
            }
        }

        let mut maxndx = [-1i64; 2];
        let mut maxsc = [0.0f64; 2];
        let mut maxigm = [0.0f64; 2];

        let lo = (ndx - 100).max(0);
        let hi = (ndx + 100).min(nn);
        let mut j = lo;
        while j < hi {
            if j == ndx {
                j += 1;
                continue;
            }
            let ju = j as usize;
            if nodes[ju].typ == STOP || nodes[ju].stop_val != nodes[ndxu].stop_val {
                j += 1;
                continue;
            }
            let mut tigm = 0.0f64;
            let mut skip = false;
            if i > 0 {
                let prev_start = genes[i - 1].start_ndx;
                let prev_stop = genes[i - 1].stop_ndx;
                if prev_start >= 0 {
                    let ps = nodes[prev_start as usize].strand;
                    if nodes[ju].strand == 1 && ps == 1 && prev_stop >= 0 {
                        if nodes[prev_stop as usize].ndx - nodes[ju].ndx > crate::node::MAX_SAM_OVLP {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[prev_stop as usize], &nodes[ju], tinf);
                        }
                    } else if nodes[ju].strand == 1 && ps == -1 {
                        if nodes[prev_start as usize].ndx - nodes[ju].ndx >= 0 {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[prev_start as usize], &nodes[ju], tinf);
                        }
                    }
                }
            }
            if !skip && i + 1 < ng {
                let next_start = genes[i + 1].start_ndx;
                let next_stop = genes[i + 1].stop_ndx;
                if next_start >= 0 {
                    let ns = nodes[next_start as usize].strand;
                    if nodes[ju].strand == -1 && ns == 1 {
                        if nodes[ju].ndx - nodes[next_start as usize].ndx >= 0 {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[ju], &nodes[next_start as usize], tinf);
                        }
                    } else if nodes[ju].strand == -1 && ns == -1 && next_stop >= 0 {
                        if nodes[ju].ndx - nodes[next_stop as usize].ndx > crate::node::MAX_SAM_OVLP {
                            skip = true;
                        } else {
                            tigm = intergenic_mod(&nodes[ju], &nodes[next_stop as usize], tinf);
                        }
                    }
                }
            }
            if skip {
                j += 1;
                continue;
            }

            let cs = nodes[ju].cscore + nodes[ju].sscore;
            if maxndx[0] == -1 {
                maxndx[0] = j;
                maxsc[0] = cs;
                maxigm[0] = tigm;
            } else if cs + tigm > maxsc[0] {
                maxndx[1] = maxndx[0];
                maxsc[1] = maxsc[0];
                maxigm[1] = maxigm[0];
                maxndx[0] = j;
                maxsc[0] = cs;
                maxigm[0] = tigm;
            } else if maxndx[1] == -1 || cs + tigm > maxsc[1] {
                maxndx[1] = j;
                maxsc[1] = cs;
                maxigm[1] = tigm;
            }
            j += 1;
        }

        for k in 0..2 {
            let mndx = maxndx[k];
            if mndx == -1 {
                continue;
            }
            let m = mndx as usize;
            let far = (nodes[m].ndx - nodes[ndxu].ndx).abs();
            if nodes[m].tscore < nodes[ndxu].tscore
                && maxsc[k] - nodes[m].tscore >= sc - nodes[ndxu].tscore + tinf.st_wt
                && nodes[m].rscore > nodes[ndxu].rscore
                && nodes[m].uscore > nodes[ndxu].uscore
                && nodes[m].cscore > nodes[ndxu].cscore
                && far > 15
            {
                maxsc[k] += nodes[ndxu].tscore - nodes[m].tscore;
            } else if far <= 15
                && nodes[m].rscore + nodes[m].tscore > nodes[ndxu].rscore + nodes[ndxu].tscore
                && !nodes[ndxu].edge
                && !nodes[m].edge
            {
                if nodes[ndxu].cscore > nodes[m].cscore {
                    maxsc[k] += nodes[ndxu].cscore - nodes[m].cscore;
                }
                if nodes[ndxu].uscore > nodes[m].uscore {
                    maxsc[k] += nodes[ndxu].uscore - nodes[m].uscore;
                }
                if igm > maxigm[k] {
                    maxsc[k] += igm - maxigm[k];
                }
            } else {
                maxsc[k] = -1000.0;
            }
        }

        let mut mndx: i64 = -1;
        for k in 0..2 {
            if maxndx[k] == -1 {
                continue;
            }
            if mndx == -1 && maxsc[k] + maxigm[k] > sc + igm {
                mndx = k as i64;
            } else if mndx >= 0
                && maxsc[k] + maxigm[k] > maxsc[mndx as usize] + maxigm[mndx as usize]
            {
                mndx = k as i64;
            }
        }
        if mndx != -1 {
            let chosen = maxndx[mndx as usize];
            let c = chosen as usize;
            if nodes[c].strand == 1 {
                genes[i].start_ndx = chosen;
                genes[i].begin = nodes[c].ndx + 1;
            } else {
                genes[i].start_ndx = chosen;
                genes[i].end = nodes[c].ndx + 1;
            }
        }
    }
}

/// Train a model on a (possibly concatenated) sequence under one genetic code.
pub fn train(seq: &Seq, trans_table: u32, closed: bool) -> Tinf {
    train_masked(seq, trans_table, closed, &[])
}

/// [`train`] with Prodigal 2.x `-m` masked regions (half-open `[begin, end)`;
/// see [`crate::seq::ambiguous_runs`]): genes crossing a mask are not used for
/// training, exactly as in Prodigal / pyrodigal `mask=True`.
pub fn train_masked(seq: &Seq, trans_table: u32, closed: bool, masks: &[(usize, usize)]) -> Tinf {
    let code = GeneticCode::new(trans_table).expect("validated table");
    // Prodigal: unambiguous G+C over the total length (N counts only in the
    // denominator). Identical to `seq.gc` for sequences without N.
    let mut tinf = Tinf::new(trans_table, seq.gc_prodigal());

    let mut nodes = add_nodes_masked(seq, &code, closed, masks);
    sort_nodes(&mut nodes);
    if nodes.is_empty() {
        return tinf;
    }
    let gc = crate::gcframe::most_gc_frame_prodigal(&seq.fwd);
    record_gc_bias(&gc, &mut nodes, &mut tinf);
    record_overlapping_starts(&mut nodes, &tinf, false);
    let ipath = dprog(&mut nodes, &tinf, false);
    calc_dicodon_gene(&mut tinf, seq, &nodes, ipath);
    raw_coding_score(seq, &mut nodes, &tinf);
    rbs_score(seq, &mut nodes, &tinf);
    train_starts_sd(seq, &mut nodes, &mut tinf);
    determine_sd_usage(&mut tinf);
    tinf
}

/// Predict genes on one contig with a trained model. Returns the genes and the
/// summed coding score (used by code auto-detection).
pub fn predict(seq: &Seq, contig: &str, tinf: &Tinf, closed: bool, is_meta: bool) -> (Vec<Gene>, f64) {
    predict_masked(seq, contig, tinf, closed, is_meta, &[])
}

/// [`predict`] with Prodigal 2.x `-m` masked regions: no gene is built across
/// (or into) a masked region.
pub fn predict_masked(
    seq: &Seq,
    contig: &str,
    tinf: &Tinf,
    closed: bool,
    is_meta: bool,
    masks: &[(usize, usize)],
) -> (Vec<Gene>, f64) {
    predict_ends(seq, contig, tinf, closed, closed, is_meta, masks)
}

/// [`predict_masked`] with each end of the sequence open or closed on its own
/// (`closed_left` / `closed_right`; see [`crate::node::add_nodes_ends`]).
pub fn predict_ends(
    seq: &Seq,
    contig: &str,
    tinf: &Tinf,
    closed_left: bool,
    closed_right: bool,
    is_meta: bool,
    masks: &[(usize, usize)],
) -> (Vec<Gene>, f64) {
    let code = tinf.code();
    let mut nodes = add_nodes_ends(seq, &code, closed_left, closed_right, masks);
    sort_nodes(&mut nodes);
    if nodes.is_empty() {
        return (Vec::new(), 0.0);
    }
    score_nodes_ends(seq, &mut nodes, tinf, closed_left, closed_right, is_meta);
    record_overlapping_starts(&mut nodes, tinf, true);
    let ipath = dprog(&mut nodes, tinf, true);
    eliminate_bad_genes(&mut nodes, ipath, tinf);
    let mut raw = add_genes(&nodes, ipath);
    tweak_final_starts(&mut raw, &nodes, tinf);

    let mut coding_total = 0.0f64;
    let mut genes = Vec::with_capacity(raw.len());
    for (k, g) in raw.iter().enumerate() {
        if g.start_ndx < 0 || g.stop_ndx < 0 {
            continue;
        }
        let sidx = g.start_ndx as usize;
        let stopidx = g.stop_ndx as usize;
        let strand = nodes[sidx].strand;
        let (begin, end) = if g.begin <= g.end {
            (g.begin, g.end)
        } else {
            (g.end, g.begin)
        };
        let begin = begin.max(1);
        let end = end.min(seq.len() as i64);
        if end <= begin {
            continue;
        }

        let cscore = nodes[sidx].cscore;
        let sscore = nodes[sidx].sscore;
        coding_total += cscore;
        let total = cscore + sscore;
        let conf = calculate_confidence(total, tinf.st_wt);

        // Prodigal `record_gene_data`: an edge start is reported as "Edge".
        let start_type = if nodes[sidx].edge {
            3u8
        } else {
            match nodes[sidx].typ {
                0 => 0u8,
                1 => 1,
                2 => 2,
                _ => 3,
            }
        };
        // Partial flags in Prodigal's genomic orientation (left = lower
        // coordinate, right = higher), taken straight from the edge nodes.
        // The old code additionally required begin == 1 / end == len, so edge
        // genes padded to start at base 2 or 3 (or end at len-1 / len-2 —
        // two of every three frames) were reported as complete.
        let (left_node, right_node) = if strand == 1 { (sidx, stopidx) } else { (stopidx, sidx) };
        let partial_left = nodes[left_node].edge;
        let partial_right = nodes[right_node].edge;
        let (rbs_motif, rbs_spacer) = rbs_strings(&nodes[sidx], tinf);

        let mut gene = Gene {
            contig: contig.to_string(),
            id: k + 1,
            begin: begin as usize,
            end: end as usize,
            strand: if strand == 1 { 1 } else { -1 },
            start_type,
            rbs_motif: rbs_motif.to_string(),
            rbs_spacer: rbs_spacer.to_string(),
            partial_left,
            partial_right,
            score: total,
            cscore,
            sscore,
            conf,
            gc: nodes[sidx].gc_cont,
            trans_table: tinf.trans_table,
            engine: Engine::Dynamic,
            aa: String::new(),
            nt: String::new(),
        };
        gene.fill_sequences(seq, &code);
        genes.push(gene);
    }
    // Renumber sequentially.
    for (k, g) in genes.iter_mut().enumerate() {
        g.id = k + 1;
    }
    (genes, coding_total)
}

/// Build a single training sequence from many contigs the way Prodigal does:
/// each contig followed by a stop-forcing linker (`TTAATTAATTAA`) — between
/// contigs *and* after the last one (a single contig is used as is).
pub fn build_training_seq(contigs: &[(String, Seq)]) -> Seq {
    if contigs.len() == 1 {
        return contigs[0].1.clone();
    }
    let mut bytes: Vec<u8> = Vec::new();
    for (_, s) in contigs {
        // Keep ambiguous bases ambiguous (`N`); decoding the 2-bit codes would
        // silently turn every N into a concrete base in the training set.
        bytes.extend_from_slice(&s.to_ascii_upper());
        bytes.extend_from_slice(TRAIN_GAP);
    }
    Seq::from_ascii(&bytes)
}

/// [`build_training_seq`] straight from raw FASTA sequences — exactly the
/// training input Prodigal builds from a multi-record file (assembly gaps
/// included; Prodigal stores N as a non-stop placeholder).
pub fn build_training_seq_ascii(records: &[&[u8]]) -> Seq {
    if records.len() == 1 {
        return Seq::from_ascii(records[0]);
    }
    // Prodigal (`read_seq_training`) and pyrodigal (`GeneFinder.train(*seqs)`)
    // put the linker after the last record too.
    let mut bytes: Vec<u8> = Vec::new();
    for r in records {
        bytes.extend_from_slice(r);
        bytes.extend_from_slice(TRAIN_GAP);
    }
    Seq::from_ascii(&bytes)
}

/// An alternative genetic code is adopted only if its coding density beats
/// the baseline and **every** evaluated candidate that keeps a superset of its
/// stop codons by more than this many density points (0.05 = 5 percentage
/// points; the same margin CheckM / GTDB-Tk use for table 4 vs 11).
pub const ALT_CODE_MIN_GAIN: f64 = 0.05;
/// ...and only if its own coding density exceeds this floor (CheckM / GTDB-Tk).
pub const ALT_CODE_MIN_DENSITY: f64 = 0.70;
/// No alternative code is adopted from less unambiguous sequence than this
/// (Prodigal's own minimum for self-training): the density estimate is too
/// noisy on a handful of genes. Force the table with `-g` instead.
pub const ALT_CODE_MIN_TRAIN_BP: usize = 20_000;
/// A genuine reassignment uses the former stop codon(s) as a *sense* codon
/// throughout genes. Inside the genes predicted under the alternative table,
/// the reassigned codons must occur at least this often per 1,000 codons
/// (real table-4 *Mycoplasma*: ~10; real table-15 phages Agate/Topaz: 11–14).
/// Fusion artefacts on a standard-code genome sit near 4.
pub const READTHROUGH_MIN_PER_KCODON: f64 = 5.0;
/// ...and at most this often: higher rates mean non-coding frames are being
/// read through, not recoded (table 6 on AT-rich *Mycoplasma*: ~85).
pub const READTHROUGH_MAX_PER_KCODON: f64 = 50.0;
/// In a recoded genome the alternative-code genes are the baseline-code gene
/// fragments **fused in the same frame**. At least this fraction of the bases
/// in the alternative-code genes must lie inside a baseline-code gene on the
/// same strand and frame (real recoded genomes: 0.53–0.69; heavy synthetic
/// recoding ≥ 0.30). Long ORFs that the baseline code calls non-coding —
/// random or non-coding sequence read under fewer stop codons — score ≤ 0.08.
pub const ALT_CODE_MIN_SUPPORT: f64 = 0.20;

/// In-gene usage (per 1,000 codons) of the codons that are stops under
/// `baseline` but sense codons under `alt` — the read-through rate. The last
/// codon of each gene is excluded; codons with ambiguous bases are skipped.
pub fn readthrough_per_kcodon(genes: &[Gene], baseline: &GeneticCode, alt: &GeneticCode) -> f64 {
    let (rt, codons) = readthrough_counts(genes.iter(), baseline, alt);
    if codons == 0 {
        0.0
    } else {
        1000.0 * rt as f64 / codons as f64
    }
}

fn readthrough_counts<'a, I: Iterator<Item = &'a Gene>>(genes: I, baseline: &GeneticCode, alt: &GeneticCode) -> (u64, u64) {
    let mut codons = 0u64;
    let mut rt = 0u64;
    for g in genes {
        let b = g.nt.as_bytes();
        let m = b.len() / 3;
        if m < 2 {
            continue;
        }
        for k in 0..m - 1 {
            let (c1, a1) = crate::seq::encode_base(b[3 * k]);
            let (c2, a2) = crate::seq::encode_base(b[3 * k + 1]);
            let (c3, a3) = crate::seq::encode_base(b[3 * k + 2]);
            if a1 || a2 || a3 {
                continue;
            }
            codons += 1;
            if baseline.is_stop(c1, c2, c3) && !alt.is_stop(c1, c2, c3) {
                rt += 1;
            }
        }
    }
    (rt, codons)
}

/// Number of bases of a `len`-long sequence covered by at least one gene.
fn covered_bases(genes: &[Gene], len: usize) -> usize {
    let mut cov = vec![false; len];
    for g in genes {
        let b = g.begin.max(1) - 1;
        let e = g.end.min(len);
        for c in cov.iter_mut().take(e).skip(b) {
            *c = true;
        }
    }
    cov.iter().filter(|&&c| c).count()
}

/// Fraction of `len` bases covered by at least one gene (coding density).
pub fn coding_density(genes: &[Gene], len: usize) -> f64 {
    if len == 0 {
        return 0.0;
    }
    covered_bases(genes, len) as f64 / len as f64
}

/// Reading frame key of a gene on a `len`-long sequence: 0–2 forward, 3–5 reverse.
#[inline]
fn frame_key(g: &Gene, len: usize) -> usize {
    if g.strand >= 0 {
        (g.begin - 1) % 3
    } else {
        3 + (len - g.end) % 3
    }
}

/// Bases of `alt` genes that lie inside a `base` gene on the same strand and
/// frame, and the total bases of `alt` genes (one sequence of length `len`).
fn same_frame_support(alt: &[Gene], base: &[Gene], len: usize) -> (u64, u64) {
    // Baseline genes in one strand/frame never overlap each other.
    let mut by_frame: [Vec<(usize, usize)>; 6] = Default::default();
    for g in base {
        by_frame[frame_key(g, len)].push((g.begin, g.end));
    }
    for v in by_frame.iter_mut() {
        v.sort_unstable();
    }
    let (mut sup, mut tot) = (0u64, 0u64);
    for g in alt {
        tot += (g.end - g.begin + 1) as u64;
        let v = &by_frame[frame_key(g, len)];
        // first interval that could overlap: begin > g.end are excluded
        let hi = v.partition_point(|&(b, _)| b <= g.end);
        for &(b, e) in v[..hi].iter().rev() {
            if e < g.begin {
                break;
            }
            let lo = b.max(g.begin);
            let up = e.min(g.end);
            if up >= lo {
                sup += (up - lo + 1) as u64;
            }
        }
    }
    (sup, tot)
}

/// Collapse a candidate list to one representative per *behaviourally
/// distinct* genetic code, in input order. Tables with the same stop codons and
/// the same start codons (e.g. 4 and 25, which differ only in what TGA encodes)
/// produce identical gene calls, so only the first one listed is evaluated.
pub fn distinct_candidates(candidates: &[u32]) -> Vec<u32> {
    let mut reps: Vec<(u32, (u64, u64))> = Vec::new();
    for &tt in candidates {
        if let Some(c) = GeneticCode::new(tt) {
            let sig = (c.stop_mask(), c.start_mask());
            if !reps.iter().any(|(t, s)| *t == tt || *s == sig) {
                reps.push((tt, sig));
            }
        }
    }
    reps.into_iter().map(|(t, _)| t).collect()
}

/// The baseline among `reps`: the first-listed candidate among those that
/// keep the most stop codons (the most conservative reading). For the standard
/// lists this is 11; for a list led by 1 (eukaryotes) it is 1.
fn baseline_code(reps: &[u32]) -> u32 {
    let stops = |tt: u32| GeneticCode::new(tt).map(|c| c.stop_mask().count_ones()).unwrap_or(0);
    let most = reps.iter().map(|&t| stops(t)).max().unwrap_or(0);
    reps.iter().copied().find(|&t| stops(t) == most).unwrap_or(reps[0])
}

/// Choose the genetic code by self-training under each candidate.
///
/// The old criterion — highest *summed coding score*, with a 10% margin over
/// the standard code — is not comparable across tables: Prodigal's length
/// factor uses a different stop-probability for every table other than 11, and
/// the start-codon set changes the score too. Even the real Prodigal scores
/// phage λ 28% higher under table 1 than table 11 (and 22% higher under 15),
/// although λ uses the standard code. DYNAMITE therefore picked table 1 for λ.
///
/// The decision now follows the established coding-density test (CheckM /
/// GTDB-Tk), generalised to several candidate codes:
///
/// 1. Candidates that make identical gene calls are collapsed
///    ([`distinct_candidates`]).
/// 2. The baseline is the first-listed candidate among those keeping the most
///    stop codons (11 for the standard lists).
/// 3. Each table is self-trained on `train_seq` and its **coding density**
///    (fraction of bases inside predicted genes) measured on the evaluation
///    sequences.
/// 4. A non-baseline table `T` is eligible only if its density exceeds
///    [`ALT_CODE_MIN_DENSITY`] *and* beats the baseline and **every** other
///    candidate whose stop-codon set contains `T`'s by more than
///    [`ALT_CODE_MIN_GAIN`]. Removing stop codons always lets ORFs read
///    further, so each removed stop has to earn its own gain: table 1 (same
///    stops as 11) can never displace 11, and table 6 (only TGA) must beat
///    both 11 and 15 — which rejects the read-through inflation 6 shows on
///    genuine table-15 phages.
/// 5. The codons `T` reassigns must be used inside its genes like sense
///    codons: between [`READTHROUGH_MIN_PER_KCODON`] and
///    [`READTHROUGH_MAX_PER_KCODON`] per 1,000 codons
///    ([`readthrough_per_kcodon`]).
/// 6. `T`'s genes must be baseline gene fragments fused in frame: at least
///    [`ALT_CODE_MIN_SUPPORT`] of their bases inside same-strand, same-frame
///    baseline genes. This rejects non-coding or random sequence, where fewer
///    stop codons simply create long ORFs the baseline code never called.
/// 7. The eligible table with the highest density wins (ties: list order);
///    if none is eligible — or there are fewer than
///    [`ALT_CODE_MIN_TRAIN_BP`] unambiguous bases — the baseline is kept.
///
/// Validated on: phage λ and 14 further table-11 phages, *E. coli* K-12
/// (→ 11); four *Mycoplasma mycoides* genomes (→ 4); the table-15 phages
/// Agate and Topaz (→ 15); random and reversed sequence (→ 11).
pub fn autodetect_code(train_seq: &Seq, candidates: &[u32], closed: bool, is_meta: bool) -> (u32, Tinf) {
    autodetect_code_on(train_seq, &[train_seq], candidates, closed, is_meta)
}

/// [`autodetect_code`] that trains on `train_seq` but measures the statistics
/// on `eval` (DYNAMITE passes the gap-free segments, so an assembly gap can
/// never count as coding sequence).
pub fn autodetect_code_on(
    train_seq: &Seq,
    eval: &[&Seq],
    candidates: &[u32],
    closed: bool,
    is_meta: bool,
) -> (u32, Tinf) {
    let reps = distinct_candidates(candidates);
    if reps.is_empty() {
        let tinf = train(train_seq, 11, closed);
        return (11, tinf);
    }
    let baseline = baseline_code(&reps);
    if reps.len() == 1 {
        return (reps[0], train(train_seq, reps[0], closed));
    }
    let unambiguous: usize = eval.iter().map(|s| s.n.iter().filter(|&&amb| !amb).count()).sum();
    if unambiguous < ALT_CODE_MIN_TRAIN_BP {
        return (baseline, train(train_seq, baseline, closed));
    }
    let base_code = GeneticCode::new(baseline).expect("validated table");
    let total_len: usize = eval.iter().map(|s| s.len()).sum();

    struct Eval {
        tt: u32,
        tinf: Tinf,
        density: f64,
        stops: u64,
        rt: f64,
        genes: Vec<Vec<Gene>>,
    }
    let mut results: Vec<Eval> = reps
        .par_iter()
        .map(|&tt| {
            let tinf = train(train_seq, tt, closed);
            let genes: Vec<Vec<Gene>> = eval
                .iter()
                .map(|s| predict(s, "__eval__", &tinf, closed, is_meta).0)
                .collect();
            let covered: usize = eval.iter().zip(&genes).map(|(s, g)| covered_bases(g, s.len())).sum();
            let code = GeneticCode::new(tt).expect("validated table");
            let (rt, codons) = readthrough_counts(genes.iter().flatten(), &base_code, &code);
            Eval {
                tt,
                tinf,
                density: covered as f64 / total_len.max(1) as f64,
                stops: code.stop_mask(),
                rt: if codons == 0 { 0.0 } else { 1000.0 * rt as f64 / codons as f64 },
                genes,
            }
        })
        .collect();

    let base_idx = results.iter().position(|r| r.tt == baseline).unwrap_or(0);
    let mut best = base_idx;
    for i in 0..results.len() {
        if i == base_idx {
            continue;
        }
        let r = &results[i];
        if r.density <= ALT_CODE_MIN_DENSITY
            || r.density - results[base_idx].density <= ALT_CODE_MIN_GAIN
            || !(READTHROUGH_MIN_PER_KCODON..=READTHROUGH_MAX_PER_KCODON).contains(&r.rt)
        {
            continue;
        }
        let beats_all_supersets = results.iter().enumerate().all(|(j, u)| {
            // `u` keeps every stop codon `T` keeps (u ⊇ T): T must clearly beat it.
            j == i || (u.stops & r.stops) != r.stops || r.density - u.density > ALT_CODE_MIN_GAIN
        });
        if !beats_all_supersets {
            continue;
        }
        let (mut sup, mut tot) = (0u64, 0u64);
        for (k, s) in eval.iter().enumerate() {
            let (a, b) = same_frame_support(&r.genes[k], &results[base_idx].genes[k], s.len());
            sup += a;
            tot += b;
        }
        if tot == 0 || (sup as f64 / tot as f64) < ALT_CODE_MIN_SUPPORT {
            continue;
        }
        if best == base_idx || r.density > results[best].density {
            best = i;
        }
    }
    let chosen = results.swap_remove(best);
    (chosen.tt, chosen.tinf)
}

// ---------------------------------------------------------------------------
// Preset entry points
// ---------------------------------------------------------------------------

/// Run the DP engine over a set of contigs with an explicit genetic-code policy
/// and engine tag.
///
/// * `fixed = Some(table)` forces a genetic code; `None` auto-detects among
///   `candidates` (see [`autodetect_code`]).
/// * Per-contig prediction is parallelised with rayon.
///
/// Returns `(genes, resolved_table)`. This is the shared core behind the
/// `prodigal`, `prodigalgv`, and `prodigal-meta` presets.
pub fn dp_call(
    contigs: &[(String, Seq)],
    fixed: Option<u32>,
    candidates: &[u32],
    is_meta: bool,
    closed: bool,
    engine: Engine,
) -> (Vec<Gene>, u32) {
    if contigs.is_empty() {
        return (Vec::new(), fixed.unwrap_or(11));
    }
    let train_seq = build_training_seq(contigs);
    dp_call_trained_on(&train_seq, contigs, fixed, candidates, is_meta, closed, engine)
}

/// [`dp_call`] with a separate training sequence: the model (and any genetic-
/// code search) is trained on `train_seq` exactly as Prodigal would — DYNAMITE
/// passes the input contigs as given, assembly gaps included — while genes are
/// predicted on `contigs` (the gap-free segments). Inserting a gap therefore
/// changes nothing away from the gap, just as in Prodigal.
///
/// An unsupported `fixed` table yields no genes (it used to panic).
#[allow(clippy::too_many_arguments)]
pub fn dp_call_trained_on(
    train_seq: &Seq,
    contigs: &[(String, Seq)],
    fixed: Option<u32>,
    candidates: &[u32],
    is_meta: bool,
    closed: bool,
    engine: Engine,
) -> (Vec<Gene>, u32) {
    dp_call_trained_on_ends(train_seq, contigs, &[], fixed, candidates, is_meta, closed, engine)
}

/// [`dp_call_trained_on`] with the two ends of each contig open or closed on
/// their own: `ends[k] = (closed_left, closed_right)` for `contigs[k]`
/// (missing entries use `closed`). Gap mode 1 closes the ends of the gap-free
/// segments that border an assembly gap. Training and code detection use
/// `closed`.
#[allow(clippy::too_many_arguments)]
pub fn dp_call_trained_on_ends(
    train_seq: &Seq,
    contigs: &[(String, Seq)],
    ends: &[(bool, bool)],
    fixed: Option<u32>,
    candidates: &[u32],
    is_meta: bool,
    closed: bool,
    engine: Engine,
) -> (Vec<Gene>, u32) {
    if contigs.is_empty() {
        return (Vec::new(), fixed.unwrap_or(11));
    }
    if let Some(id) = fixed {
        if GeneticCode::new(id).is_none() {
            return (Vec::new(), id);
        }
    }
    let (code_id, tinf) = match fixed {
        Some(id) => (id, train(train_seq, id, closed)),
        None => {
            let eval: Vec<&Seq> = contigs.iter().map(|(_, s)| s).collect();
            autodetect_code_on(train_seq, &eval, candidates, closed, is_meta)
        }
    };
    let code = tinf.code();
    let genes: Vec<Gene> = contigs
        .par_iter()
        .enumerate()
        .flat_map_iter(|(k, (id, s))| {
            let (closed_left, closed_right) = ends.get(k).copied().unwrap_or((closed, closed));
            let (mut g, _score) = predict_ends(s, id, &tinf, closed_left, closed_right, is_meta, &[]);
            for x in &mut g {
                x.engine = engine;
            }
            // Trim terminal gap codons / infer partial ends (no-op for clean
            // sequence: Prodigal's edge flags already agree with the codons).
            crate::finalize::finalize_genes_ends(
                g,
                s,
                &code,
                closed_left,
                closed_right,
                crate::finalize::SeqSource::Genome,
                crate::finalize::StartRule::Table,
            )
            .into_iter()
        })
        .collect();
    (genes, code_id)
}

/// The Prodigal engines on whole contigs as Prodigal 2.x calls them — the
/// core of gap mode 2 and `-m` (gap modes 0 and 1 split contigs at gaps
/// instead, see [`dp_call_trained_on_ends`]). The calls are Prodigal's own,
/// unchanged.
///
/// * `table` — the genetic code (already resolved; an unsupported one yields
///   no genes);
/// * `train_masks` — masked regions of `train_seq` (only `-m` masks training);
/// * `masks[k]` — masked regions of `contigs[k]` (half-open `[begin, end)`,
///   see [`crate::seq::ambiguous_runs`]); no gene crosses or runs into one.
#[allow(clippy::too_many_arguments)]
pub fn dp_call_masked(
    train_seq: &Seq,
    train_masks: &[(usize, usize)],
    contigs: &[(String, Seq)],
    masks: &[Vec<(usize, usize)>],
    table: u32,
    is_meta: bool,
    closed: bool,
    engine: Engine,
) -> Vec<Gene> {
    if GeneticCode::new(table).is_none() || contigs.is_empty() {
        return Vec::new();
    }
    let tinf = train_masked(train_seq, table, closed, train_masks);
    contigs
        .par_iter()
        .enumerate()
        .flat_map_iter(|(k, (id, s))| {
            let m: &[(usize, usize)] = masks.get(k).map(|v| v.as_slice()).unwrap_or(&[]);
            let (mut g, _score) = predict_masked(s, id, &tinf, closed, is_meta, m);
            for x in &mut g {
                x.engine = engine;
            }
            g.into_iter()
        })
        .collect()
}

/// **Prodigal** preset — bacteria and archaea: a single self-trained genome
/// under the standard bacterial/archaeal code 11 (override with `fixed`).
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    dp_call(contigs, Some(fixed.unwrap_or(11)), &[11], false, closed, Engine::Prodigal)
}

/// **Prodigal** preset with genetic-code auto-detection among `candidates`
/// (used for `-e prodigal -g auto`).
pub fn call_auto(contigs: &[(String, Seq)], candidates: &[u32], closed: bool) -> (Vec<Gene>, u32) {
    dp_call(contigs, None, candidates, false, closed, Engine::Prodigal)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_record_training_set_ends_with_a_linker() {
        // Prodigal / pyrodigal: rec1 + linker + rec2 + linker.
        let t = build_training_seq_ascii(&[b"ACGTAC".as_slice(), b"GGG".as_slice()]);
        assert_eq!(t.to_ascii_upper(), b"ACGTACTTAATTAATTAAGGGTTAATTAATTAA".to_vec());
        let one = build_training_seq_ascii(&[b"ACGTAC".as_slice()]);
        assert_eq!(one.to_ascii_upper(), b"ACGTAC".to_vec());
    }
}
