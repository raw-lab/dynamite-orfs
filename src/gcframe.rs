//! GC-by-frame plots.
//!
//! Two different but related constructions are needed:
//!
//! * [`most_gc_frame`] reproduces Prodigal's `calc_most_gc_frame`: for every
//!   base it returns which of the three reading frames (0/1/2) has the highest
//!   GC content in a 120 bp window. This seeds the self-training step of the
//!   dynamic-programming engine.
//!
//! * [`GcFrame`] reproduces PHANOTATE's sliding `GCframe`: it yields, for every
//!   genome position, the windowed GC counts of the three frames, which the
//!   PHANOTATE-style ORF scorer uses to up- or down-weight candidate ORFs.

use crate::seq::is_gc;

const WINDOW: usize = 120;

/// Value of [`most_gc_frame_prodigal`] for the trailing 1–2 bases that do not
/// form a full codon triple. Prodigal leaves these at `-1`; callers that need
/// Prodigal's exact arithmetic must treat this sentinel as `-1`.
pub const GC_FRAME_UNSET: u8 = u8::MAX;

/// Prodigal's GC frame plot: `out[i] ∈ {0,1,2}` is the codon frame with the
/// highest GC content in a window centred on base `i`. The last `len % 3`
/// bases (no full triple) are reported as `0`; see [`most_gc_frame_prodigal`]
/// for Prodigal's exact `-1` sentinel.
pub fn most_gc_frame(seq: &[u8]) -> Vec<u8> {
    most_gc_frame_prodigal(seq)
        .into_iter()
        .map(|v| if v == GC_FRAME_UNSET { 0 } else { v })
        .collect()
}

/// [`most_gc_frame`] with Prodigal's `-1` for the trailing bases, encoded as
/// [`GC_FRAME_UNSET`] (Prodigal's reverse-strand GC-bias pass reads them).
pub fn most_gc_frame_prodigal(seq: &[u8]) -> Vec<u8> {
    let slen = seq.len();
    let mut gp = vec![GC_FRAME_UNSET; slen.max(1)];
    if slen == 0 {
        return gp;
    }
    let mut fwd = vec![0i64; slen];
    let mut bwd = vec![0i64; slen];
    let mut tot = vec![0i64; slen];

    let gcv = |i: usize| -> i64 { if is_gc(seq, i) { 1 } else { 0 } };

    for i in 0..3 {
        let mut j = i;
        while j < slen {
            if j < 3 {
                fwd[j] = gcv(j);
            } else {
                fwd[j] = fwd[j - 3] + gcv(j);
            }
            let bj = slen - j - 1;
            if j < 3 {
                bwd[bj] = gcv(bj);
            } else {
                bwd[bj] = bwd[slen - j + 2] + gcv(bj);
            }
            j += 3;
        }
    }
    for i in 0..slen {
        tot[i] = fwd[i] + bwd[i] - gcv(i);
        if i >= WINDOW / 2 {
            tot[i] -= fwd[i - WINDOW / 2];
        }
        if i + WINDOW / 2 < slen {
            tot[i] -= bwd[i + WINDOW / 2];
        }
    }
    let mut i = 0;
    while i + 2 < slen {
        let win = max_fr(tot[i], tot[i + 1], tot[i + 2]);
        for j in 0..3 {
            gp[i + j] = win;
        }
        i += 3;
    }
    gp
}

/// Three-way arg-max returning `0`, `1` or `2` (Prodigal's `max_fr`).
#[inline]
pub fn max_fr(n1: i64, n2: i64, n3: i64) -> u8 {
    if n1 > n2 {
        if n1 > n3 {
            0
        } else {
            2
        }
    } else if n2 > n3 {
        1
    } else {
        2
    }
}

/// PHANOTATE-style 1-based arg-max over a codon's three frame GC counts.
#[inline]
pub fn max_idx(a: i64, b: i64, c: i64) -> usize {
    if a > b {
        if a > c {
            1
        } else {
            3
        }
    } else if b > c {
        2
    } else {
        3
    }
}

/// PHANOTATE-style 1-based arg-min over a codon's three frame GC counts.
#[inline]
pub fn min_idx(a: i64, b: i64, c: i64) -> usize {
    if a > b {
        if b > c {
            3
        } else {
            2
        }
    } else if a > c {
        3
    } else {
        1
    }
}

/// Faithful port of PHANOTATE's sliding `GCframe`.
pub struct GcFrame {
    window: usize,
    state: usize,
    // Per frame: window contents (true = GC base, false = AT, treated as a ring
    // via a VecDeque), the running GC count, and the per-step total series.
    bases: [std::collections::VecDeque<bool>; 3],
    gc: [i64; 3],
    total: [std::collections::VecDeque<i64>; 3],
}

impl Default for GcFrame {
    fn default() -> Self {
        GcFrame::new(WINDOW)
    }
}

impl GcFrame {
    /// New plot with the given nucleotide window (PHANOTATE uses 120).
    pub fn new(window: usize) -> Self {
        let w = window / 3;
        let mk = || std::collections::VecDeque::from(vec![false; w]);
        GcFrame {
            window: w,
            state: 0,
            bases: [mk(), mk(), mk()],
            gc: [0, 0, 0],
            total: [
                std::collections::VecDeque::new(),
                std::collections::VecDeque::new(),
                std::collections::VecDeque::new(),
            ],
        }
    }

    /// Feed one base (its raw 2-bit code is enough; only GC-ness matters).
    pub fn add_base(&mut self, code: u8) {
        let frame = self.state;
        self.state = (self.state + 1) % 3;
        let g = code == crate::seq::C || code == crate::seq::G;
        self.bases[frame].push_back(g);
        if g {
            self.gc[frame] += 1;
        }
        if let Some(old) = self.bases[frame].pop_front() {
            if old {
                self.gc[frame] -= 1;
            }
        }
        self.total[frame].push_back(self.gc[frame]);
    }

    fn close(&mut self) {
        for _ in 0..(self.window / 2) {
            for frame in 0..3 {
                self.total[frame].pop_front();
                if let Some(old) = self.bases[frame].pop_front() {
                    if old {
                        self.gc[frame] -= 1;
                    }
                }
                self.total[frame].push_back(self.gc[frame]);
            }
        }
    }

    /// Finalise and return the per-position frame-GC triples. The result is
    /// indexed by 1-based genome coordinate (`out[0]` is a sentinel).
    pub fn get(mut self) -> Vec<[i64; 3]> {
        self.close();
        let t: [Vec<i64>; 3] = [
            self.total[0].iter().copied().collect(),
            self.total[1].iter().copied().collect(),
            self.total[2].iter().copied().collect(),
        ];
        let mut freq: Vec<[i64; 3]> = Vec::new();
        freq.push([20, 20, 20]);
        let get = |v: &Vec<i64>, i: usize| -> i64 { v.get(i).copied().unwrap_or(0) };

        if t[2].len() >= 1 {
            let lim = t[2].len() - 1;
            let mut i = 0;
            while i < lim {
                freq.push([get(&t[0], i), get(&t[1], i), get(&t[2], i)]);
                freq.push([get(&t[1], i), get(&t[2], i), get(&t[0], i + 1)]);
                freq.push([get(&t[2], i), get(&t[0], i + 1), get(&t[1], i + 1)]);
                i += 1;
            }
            // tail (mirrors PHANOTATE's post-loop appends)
            freq.push([get(&t[0], i), get(&t[1], i), get(&t[2], i)]);
            if i < t[0].len().saturating_sub(1) {
                freq.push([get(&t[1], i), get(&t[2], i), get(&t[0], i + 1)]);
            }
            if i < t[1].len().saturating_sub(1) {
                freq.push([get(&t[2], i), get(&t[0], i + 1), get(&t[1], i + 1)]);
            }
        }
        freq
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Seq;

    #[test]
    fn gc_frame_lengths() {
        let s = Seq::from_ascii(&b"ATGAAACCCGGGTTTACGTACGTACGT".repeat(20));
        let gp = most_gc_frame(&s.fwd);
        assert_eq!(gp.len(), s.len());

        let mut f = GcFrame::default();
        for &c in &s.fwd {
            f.add_base(c);
        }
        let freq = f.get();
        // at least one entry per base region
        assert!(freq.len() >= s.len() / 2);
    }

    #[test]
    fn idx_helpers() {
        assert_eq!(max_idx(5, 1, 1), 1);
        assert_eq!(max_idx(1, 5, 1), 2);
        assert_eq!(max_idx(1, 1, 5), 3);
        assert_eq!(min_idx(5, 1, 9), 2);
        assert_eq!(max_fr(1, 5, 1), 1);
    }
}
