//! Codon usage and codon-bias statistics computed natively (no external deps).
//!
//! Given a set of coding sequences (in-frame nucleotides), this computes:
//! per-codon counts and frequencies, RSCU (relative synonymous codon usage),
//! GC / GC1 / GC2 / GC3, the effective number of codons (ENC, Wright 1990),
//! and codon density (codons per kb). Output is plain TSV.
//!
//! The same quantities are produced by the RAW-lab `rustyomestats` crate; this
//! module reimplements the codon-bias subset so DYNAMITE stays dependency-free
//! on the default build. Enable the `rustyomestats` feature to additionally
//! cross-check against that crate in environments with Rust 1.85+.

use crate::gencode::GeneticCode;
use crate::seq::encode_base;
use std::fmt::Write as _;

const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];

/// Result of a codon-usage analysis.
pub struct CodonStats {
    pub n_cds: usize,
    pub total_codons: u64,
    pub total_nt: u64,
    pub counts: [u64; 64],
    pub gc: f64,
    pub gc1: f64,
    pub gc2: f64,
    pub gc3: f64,
    pub enc: f64,
    pub trans_table: u32,
}

#[inline]
fn codon_index(b1: u8, b2: u8, b3: u8) -> Option<usize> {
    let (c1, a1) = encode_base(b1);
    let (c2, a2) = encode_base(b2);
    let (c3, a3) = encode_base(b3);
    if a1 || a2 || a3 {
        return None; // ambiguous base — skip
    }
    // index in TCAG order to match conventional codon tables
    let map = |c: u8| match c {
        3 => 0, // T
        1 => 1, // C
        0 => 2, // A
        2 => 3, // G
        _ => 0,
    };
    Some(map(c1) * 16 + map(c2) * 4 + map(c3))
}

fn codon_text(idx: usize) -> [u8; 3] {
    [BASES[idx / 16], BASES[(idx / 4) % 4], BASES[idx % 4]]
}

/// Amino acid for an ASCII codon. `GeneticCode::translate` takes 2-bit base
/// codes, so convert first.
fn aa_of(code: &GeneticCode, c: [u8; 3]) -> u8 {
    let (a, _) = encode_base(c[0]);
    let (b, _) = encode_base(c[1]);
    let (d, _) = encode_base(c[2]);
    code.translate(a, b, d)
}

/// Analyse a set of in-frame coding sequences.
pub fn analyze(cds: &[String], code: &GeneticCode) -> CodonStats {
    let mut counts = [0u64; 64];
    let mut gc = 0u64;
    let (mut gc1, mut gc2, mut gc3) = (0u64, 0u64, 0u64);
    let (mut n1, mut n2, mut n3) = (0u64, 0u64, 0u64);
    let mut total_nt = 0u64;

    for seq in cds {
        let b = seq.as_bytes();
        total_nt += b.len() as u64;
        let mut i = 0;
        while i + 3 <= b.len() {
            if let Some(ci) = codon_index(b[i], b[i + 1], b[i + 2]) {
                counts[ci] += 1;
            }
            for (k, (cnt, ncnt)) in [(&mut gc1, &mut n1), (&mut gc2, &mut n2), (&mut gc3, &mut n3)]
                .into_iter()
                .enumerate()
            {
                let (c, amb) = encode_base(b[i + k]);
                if !amb {
                    *ncnt += 1;
                    if c == 1 || c == 2 {
                        *cnt += 1;
                        gc += 1;
                    }
                }
            }
            i += 3;
        }
    }
    let total_codons: u64 = counts.iter().sum();
    let denom = (n1 + n2 + n3).max(1) as f64;
    let stats = CodonStats {
        n_cds: cds.len(),
        total_codons,
        total_nt,
        counts,
        gc: gc as f64 / denom,
        gc1: gc1 as f64 / n1.max(1) as f64,
        gc2: gc2 as f64 / n2.max(1) as f64,
        gc3: gc3 as f64 / n3.max(1) as f64,
        enc: 0.0,
        trans_table: code.id,
    };
    let enc = effective_number_of_codons(&stats, code);
    CodonStats { enc, ..stats }
}

/// RSCU for one codon: observed / expected-under-uniform-synonymous-usage.
fn rscu_table(stats: &CodonStats, code: &GeneticCode) -> [f64; 64] {
    // Group codon indices by amino acid (using the genetic code).
    let mut by_aa: std::collections::HashMap<u8, Vec<usize>> = std::collections::HashMap::new();
    for idx in 0..64 {
        let c = codon_text(idx);
        let aa = aa_of(code, c);
        by_aa.entry(aa).or_default().push(idx);
    }
    let mut rscu = [0.0f64; 64];
    for (_aa, group) in by_aa {
        let total: u64 = group.iter().map(|&i| stats.counts[i]).sum();
        let k = group.len() as f64;
        if total == 0 {
            continue;
        }
        for &i in &group {
            rscu[i] = stats.counts[i] as f64 * k / total as f64;
        }
    }
    rscu
}

/// Wright's effective number of codons (ENC), 1990. Averages within the
/// 2-, 3-, 4-, and 6-fold synonymous degeneracy classes.
fn effective_number_of_codons(stats: &CodonStats, code: &GeneticCode) -> f64 {
    let mut by_aa: std::collections::HashMap<u8, Vec<usize>> = std::collections::HashMap::new();
    for idx in 0..64 {
        let c = codon_text(idx);
        let aa = aa_of(code, c);
        if aa == b'*' {
            continue; // stops are not synonymous families
        }
        by_aa.entry(aa).or_default().push(idx);
    }
    // F-value per amino acid: F = (N*sum p_i^2 - 1)/(N-1), grouped by degeneracy.
    let mut fsum: std::collections::HashMap<usize, (f64, usize)> = std::collections::HashMap::new();
    for (_aa, group) in &by_aa {
        let deg = group.len();
        if deg < 2 {
            continue;
        }
        let n: u64 = group.iter().map(|&i| stats.counts[i]).sum();
        if n < 2 {
            continue;
        }
        let nf = n as f64;
        let sum_p2: f64 = group
            .iter()
            .map(|&i| {
                let p = stats.counts[i] as f64 / nf;
                p * p
            })
            .sum();
        let f = (nf * sum_p2 - 1.0) / (nf - 1.0);
        if f > 0.0 {
            let e = fsum.entry(deg).or_insert((0.0, 0));
            e.0 += f;
            e.1 += 1;
        }
    }
    let avg = |deg: usize| -> Option<f64> {
        fsum.get(&deg).map(|(s, c)| s / *c as f64)
    };
    // Standard Wright formula with class fallbacks.
    let f2 = avg(2).unwrap_or(0.5);
    let f3 = avg(3).or_else(|| {
        // interpolate if missing
        match (avg(2), avg(4)) {
            (Some(a), Some(b)) => Some((a + b) / 2.0),
            _ => Some(0.5),
        }
    }).unwrap_or(0.5);
    let f4 = avg(4).unwrap_or(0.5);
    let f6 = avg(6).unwrap_or(f4);
    let mut enc = 2.0 + 9.0 / f2 + 1.0 / f3 + 5.0 / f4 + 3.0 / f6;
    if enc > 61.0 {
        enc = 61.0;
    }
    enc
}

impl CodonStats {
    /// Render the analysis as TSV (summary header + per-codon table).
    pub fn to_tsv(&self, code: &GeneticCode) -> String {
        let rscu = rscu_table(self, code);
        let mut out = String::new();
        let _ = writeln!(out, "# DYNAMITE codon usage");
        let _ = writeln!(out, "# trans_table\t{}", self.trans_table);
        let _ = writeln!(out, "# CDS\t{}", self.n_cds);
        let _ = writeln!(out, "# total_codons\t{}", self.total_codons);
        let kb = (self.total_nt as f64 / 1000.0).max(1e-9);
        let _ = writeln!(out, "# codon_density_per_kb\t{:.3}", self.total_codons as f64 / kb);
        let _ = writeln!(out, "# GC\t{:.4}", self.gc);
        let _ = writeln!(out, "# GC1\t{:.4}", self.gc1);
        let _ = writeln!(out, "# GC2\t{:.4}", self.gc2);
        let _ = writeln!(out, "# GC3\t{:.4}", self.gc3);
        let _ = writeln!(out, "# ENC\t{:.2}", self.enc);
        let _ = writeln!(out, "codon\taa\tcount\tfraction\trscu");
        let total = self.total_codons.max(1) as f64;
        for idx in 0..64 {
            let c = codon_text(idx);
            let aa = aa_of(code, c) as char;
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{:.5}\t{:.3}",
                std::str::from_utf8(&c).unwrap(),
                aa,
                self.counts[idx],
                self.counts[idx] as f64 / total,
                rscu[idx],
            );
        }
        out
    }
}
