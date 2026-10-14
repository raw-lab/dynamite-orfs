//! Ribosome-binding-site (Shine-Dalgarno) scoring.
//!
//! Two scorers are provided, one per engine, both returning a motif *bin* in
//! `0..28`:
//!
//! * [`sd_exact`] / [`sd_mm`] are direct ports of Prodigal's
//!   `shine_dalgarno_exact` / `shine_dalgarno_mm`. They search a window
//!   upstream of a start for the best exact / single-mismatch match to the
//!   canonical `AGGAGG` motif and map `(score, spacer)` to a bin.
//!
//! * [`phanotate_score_rbs`] is a direct port of PHANOTATE's `score_rbs`,
//!   a hand-tuned lookup over a 21 bp reversed window.

use crate::seq::{A, G};

#[inline]
fn is_a(seq: &[u8], i: usize) -> bool {
    seq[i] == A
}
#[inline]
fn is_g(seq: &[u8], i: usize) -> bool {
    seq[i] == G
}

/// Map an exact-match SD `(score, spacer-flag)` pair to a bin (Prodigal).
fn exact_bin(cc: i32, dis: i32) -> usize {
    match (cc, dis) {
        (6, 2) => 1,
        (6, 3) => 2,
        (8, 3) => 3,
        (9, 3) => 3,
        (6, 1) => 6,
        (11, 3) => 10,
        (12, 3) => 10,
        (14, 3) => 10,
        (8, 2) => 11,
        (9, 2) => 11,
        (8, 1) => 12,
        (9, 1) => 12,
        (6, 0) => 13,
        (8, 0) => 15,
        (9, 0) => 16,
        (11, 2) => 20,
        (11, 1) => 21,
        (11, 0) => 22,
        (12, 2) => 20,
        (12, 1) => 23,
        (12, 0) => 24,
        (14, 2) => 25,
        (14, 1) => 26,
        (14, 0) => 27,
        _ => 0,
    }
}

/// Map a single-mismatch SD `(score, spacer-flag)` pair to a bin (Prodigal).
fn mm_bin(cc: i32, dis: i32) -> usize {
    match (cc, dis) {
        (6, 3) => 2,
        (7, 3) => 2,
        (9, 3) => 3,
        (6, 2) => 4,
        (6, 1) => 5,
        (6, 0) => 9,
        (7, 2) => 7,
        (7, 1) => 8,
        (7, 0) => 14,
        (9, 2) => 17,
        (9, 1) => 18,
        (9, 0) => 19,
        _ => 0,
    }
}

/// Prodigal's exact SD search. `pos` and `start` are signed because the search
/// window can begin upstream of the contig start.
pub fn sd_exact(seq: &[u8], pos: isize, start: isize, rwt: &[f64; 28]) -> usize {
    let limit = (start - 4 - pos).min(6);
    let mut m = [-10.0f64; 6];
    let lim_u = limit.max(0).min(6) as usize;
    for i in 0..lim_u {
        let p = pos + i as isize;
        if p >= 0 {
            let pu = p as usize;
            if pu < seq.len() {
                if i % 3 == 0 && is_a(seq, pu) {
                    m[i] = 2.0;
                } else if i % 3 != 0 && is_g(seq, pu) {
                    m[i] = 3.0;
                }
            }
        }
    }
    let mut max_val = 0usize;
    let mut i = limit;
    while i >= 3 {
        let span = i as usize;
        let mut j = 0isize;
        while j <= limit - i {
            let jj = j as usize;
            let mut cur = -2.0f64;
            let mut mism = 0;
            for k in jj..(jj + span) {
                cur += m[k];
                if m[k] < 0.0 {
                    mism += 1;
                }
            }
            if mism == 0 {
                let rdis = start - (pos + j + i);
                let dis = spacer_flag_exact(rdis, span);
                if !(rdis > 15 || cur < 6.0) {
                    let cc = cur.round() as i32;
                    let cv = exact_bin(cc, dis);
                    if !(rwt[cv] < rwt[max_val]
                        || (rwt[cv] == rwt[max_val] && cv < max_val))
                    {
                        max_val = cv;
                    }
                }
            }
            j += 1;
        }
        i -= 1;
    }
    max_val
}

fn spacer_flag_exact(rdis: isize, span: usize) -> i32 {
    if rdis < 5 && span < 5 {
        2
    } else if rdis < 5 && span >= 5 {
        1
    } else if rdis > 10 && rdis <= 12 && span < 5 {
        1
    } else if rdis > 10 && rdis <= 12 && span >= 5 {
        2
    } else if rdis >= 13 {
        3
    } else {
        0
    }
}

/// Prodigal's single-mismatch SD search.
pub fn sd_mm(seq: &[u8], pos: isize, start: isize, rwt: &[f64; 28]) -> usize {
    let limit = (start - 4 - pos).min(6);
    let mut m = [-10.0f64; 6];
    let lim_u = limit.max(0).min(6) as usize;
    for i in 0..lim_u {
        let p = pos + i as isize;
        if p >= 0 {
            let pu = p as usize;
            if pu < seq.len() {
                if i % 3 == 0 {
                    m[i] = if is_a(seq, pu) { 2.0 } else { -3.0 };
                } else {
                    m[i] = if is_g(seq, pu) { 3.0 } else { -2.0 };
                }
            }
        }
    }
    let mut max_val = 0usize;
    let mut i = limit;
    while i >= 5 {
        let span = i as usize;
        let mut j = 0isize;
        while j <= limit - i {
            let jj = j as usize;
            let mut cur = -2.0f64;
            let mut mism = 0;
            for k in jj..(jj + span) {
                cur += m[k];
                if m[k] < 0.0 {
                    mism += 1;
                    if k <= jj + 1 || k + 2 >= jj + span {
                        cur -= 10.0;
                    }
                }
            }
            if mism == 1 {
                let rdis = start - (pos + j + i);
                let dis = if rdis < 5 {
                    1
                } else if rdis > 10 && rdis <= 12 {
                    2
                } else if rdis >= 13 {
                    3
                } else {
                    0
                };
                if !(rdis > 15 || cur < 6.0) {
                    let cc = cur.round() as i32;
                    let cv = mm_bin(cc, dis);
                    if !(rwt[cv] < rwt[max_val]
                        || (rwt[cv] == rwt[max_val] && cv < max_val))
                    {
                        max_val = cv;
                    }
                }
            }
            j += 1;
        }
        i -= 1;
    }
    max_val
}

// ----------------------------------------------------------------------------
// PHANOTATE score_rbs
// ----------------------------------------------------------------------------

#[inline]
fn eq_any(rev: &[u8], motif: &[u8], starts: &[usize]) -> bool {
    let l = motif.len();
    for &s in starts {
        if s + l <= rev.len() && &rev[s..s + l] == motif {
            return true;
        }
    }
    false
}

/// PHANOTATE's `score_rbs`: classify the 21 bp window immediately upstream of a
/// start (passed 5'->3') into a bin `0..28`. The window is reversed internally,
/// matching the reference implementation.
pub fn phanotate_score_rbs(window: &[u8]) -> usize {
    // Reverse the (≤21 bp) window, lower-cased.
    let mut rev: Vec<u8> = window.iter().rev().map(|b| b.to_ascii_lowercase()).collect();
    // Pad to at least 21 so all slice offsets are valid (PHANOTATE relies on
    // Python slicing simply yielding shorter strings off the end).
    if rev.len() < 24 {
        rev.resize(24, b'-');
    }
    let s = &rev[..];

    let r = |a: usize, b: usize| -> &[u8] { &s[a..b.min(s.len())] };
    let _ = r; // (kept for parity/readability; eq_any used below)

    if eq_any(s, b"ggagga", &[5, 6, 7, 8, 9, 10]) {
        27
    } else if eq_any(s, b"ggagga", &[3, 4]) {
        26
    } else if eq_any(s, b"ggagga", &[11, 12]) {
        25
    } else if eq_any(s, b"ggagg", &[5, 6, 7, 8, 9, 10]) {
        24
    } else if eq_any(s, b"ggagg", &[3, 4]) {
        23
    } else if eq_any(s, b"gagga", &[5, 6, 7, 8, 9, 10]) {
        22
    } else if eq_any(s, b"gagga", &[3, 4]) {
        21
    } else if eq_any(s, b"gagga", &[11, 12]) || eq_any(s, b"ggagg", &[11, 12]) {
        20
    } else if eq_any(s, b"ggacga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggatga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggaaga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggcgga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggggga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggtgga", &[5, 6, 7, 8, 9, 10])
    {
        19
    } else if eq_any(s, b"ggaaga", &[3, 4])
        || eq_any(s, b"ggatga", &[3, 4])
        || eq_any(s, b"ggacga", &[3, 4])
        || eq_any(s, b"ggtgga", &[3, 4])
        || eq_any(s, b"ggggga", &[3, 4])
        || eq_any(s, b"ggcgga", &[3, 4])
    {
        18
    } else if eq_any(s, b"ggaaga", &[11, 12])
        || eq_any(s, b"ggatga", &[11, 12])
        || eq_any(s, b"ggacga", &[11, 12])
        || eq_any(s, b"ggtgga", &[11, 12])
        || eq_any(s, b"ggggga", &[11, 12])
        || eq_any(s, b"ggcgga", &[11, 12])
    {
        17
    } else if eq_any(s, b"ggag", &[5, 6, 7, 8, 9, 10]) || eq_any(s, b"gagg", &[5, 6, 7, 8, 9, 10])
    {
        16
    } else if eq_any(s, b"agga", &[5, 6, 7, 8, 9, 10]) {
        15
    } else if eq_any(s, b"ggtgg", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggggg", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"ggcgg", &[5, 6, 7, 8, 9, 10])
    {
        14
    } else if eq_any(s, b"agg", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"gag", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"gga", &[5, 6, 7, 8, 9, 10])
    {
        13
    } else if eq_any(s, b"agga", &[11, 12]) || eq_any(s, b"gagg", &[11, 12]) || eq_any(s, b"ggag", &[11, 12])
    {
        12
    } else if eq_any(s, b"agga", &[3, 4]) || eq_any(s, b"gagg", &[3, 4]) || eq_any(s, b"ggag", &[3, 4])
    {
        11
    } else if eq_any(s, b"gagga", &[13, 14, 15])
        || eq_any(s, b"ggagg", &[13, 14, 15])
        || eq_any(s, b"ggagga", &[13, 14, 15])
    {
        10
    } else if eq_any(s, b"gaaga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"gatga", &[5, 6, 7, 8, 9, 10])
        || eq_any(s, b"gacga", &[5, 6, 7, 8, 9, 10])
    {
        9
    } else if eq_any(s, b"ggtgg", &[3, 4]) || eq_any(s, b"ggggg", &[3, 4]) || eq_any(s, b"ggcgg", &[3, 4])
    {
        8
    } else if eq_any(s, b"ggtgg", &[11, 12]) || eq_any(s, b"ggggg", &[11, 12]) || eq_any(s, b"ggcgg", &[11, 12])
    {
        7
    } else if eq_any(s, b"agg", &[11, 12]) || eq_any(s, b"gag", &[11, 12]) || eq_any(s, b"gga", &[11, 12])
    {
        6
    } else if eq_any(s, b"gaaga", &[3, 4]) || eq_any(s, b"gatga", &[3, 4]) || eq_any(s, b"gacga", &[3, 4])
    {
        5
    } else if eq_any(s, b"gaaga", &[11, 12]) || eq_any(s, b"gatga", &[11, 12]) || eq_any(s, b"gacga", &[11, 12])
    {
        4
    } else if eq_any(s, b"agga", &[13, 14, 15]) || eq_any(s, b"gagg", &[13, 14, 15]) || eq_any(s, b"ggag", &[13, 14, 15])
    {
        3
    } else if eq_any(s, b"agg", &[13, 14, 15])
        || eq_any(s, b"gag", &[13, 14, 15])
        || eq_any(s, b"gga", &[13, 14, 15])
        || eq_any(s, b"ggaaga", &[13, 14, 15])
        || eq_any(s, b"ggatga", &[13, 14, 15])
        || eq_any(s, b"ggacga", &[13, 14, 15])
        || eq_any(s, b"ggtgg", &[13, 14, 15])
        || eq_any(s, b"ggggg", &[13, 14, 15])
        || eq_any(s, b"ggcgg", &[13, 14, 15])
    {
        2
    } else if eq_any(s, b"agg", &[3, 4]) || eq_any(s, b"gag", &[3, 4]) || eq_any(s, b"gga", &[3, 4]) {
        1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Seq;

    #[test]
    fn perfect_sd_exact() {
        // upstream ...AGGAGG.....ATG : an ideal SD should land in a high bin.
        // Build an encoded sequence containing AGGAGG then a spacer then ATG.
        let s = Seq::from_ascii(b"AAAAAAAGGAGGAAAAAAATG");
        // start (ATG) is at index 18; search window upstream.
        let rwt = [1.0f64; 28];
        let mut best = 0;
        for j in (18isize - 20)..=(18 - 6) {
            best = best.max(sd_exact(&s.fwd, j, 18, &rwt));
        }
        assert!(best > 0, "expected a non-zero SD bin, got {best}");
    }

    #[test]
    fn phanotate_rbs_detects_motif() {
        // window upstream of a start containing ggagg
        let w = b"aaaaaggaggaaaaaaaaaaa"; // 21 bp
        let bin = phanotate_score_rbs(w);
        assert!(bin > 0, "expected non-zero PHANOTATE RBS bin, got {bin}");
    }
}
