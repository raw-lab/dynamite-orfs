//! A 3-periodic (codon-position-specific) inhomogeneous Markov model of coding
//! potential, plus a homogeneous non-coding background — the statistical core
//! shared by the GeneMark- and FragGeneScan-style engines.
//!
//! This is an order-2 model: the probability of each base depends on the two
//! preceding bases and, for the coding model, on the codon position (0/1/2).
//! Scoring an ORF returns a coding/non-coding log-likelihood ratio (LLR);
//! positive means "looks coding".

use crate::seq::Seq;

const NCTX: usize = 16; // 4^2 contexts (previous two bases)

/// One contig and the coding ORFs (1-based `begin`, `end`, `strand`) found on it.
pub type TrainingSet<'a> = (&'a Seq, &'a [(usize, usize, i8)]);

/// A trained coding + background model.
#[derive(Clone)]
pub struct MarkovModel {
    /// log P(base | prev2, codon_pos) for the coding model.
    cod: [[[f64; 4]; NCTX]; 3],
    /// log P(base | prev2) for the non-coding background.
    nc: [[f64; 4]; NCTX],
    pub trained_on: usize,
}

impl MarkovModel {
    /// Train from a set of presumed-coding ORFs (global 1-based coords + strand)
    /// against a whole-sequence background.
    pub fn train(seq: &Seq, coding: &[(usize, usize, i8)]) -> MarkovModel {
        MarkovModel::train_pooled(&[(seq, coding)])
    }

    /// Train one model from several contigs, each paired with the coding ORFs
    /// that were found **on that contig** (1-based coords + strand). Background
    /// and coding counts are pooled across all of them. ORFs that do not fit
    /// their contig are skipped (they used to be applied to the wrong contig,
    /// which corrupted the model and could index out of bounds).
    pub fn train_pooled(sets: &[TrainingSet<'_>]) -> MarkovModel {
        // Pseudocounts of 1 everywhere.
        let mut cod = [[[1.0f64; 4]; NCTX]; 3];
        let mut nc = [[1.0f64; 4]; NCTX];
        let mut trained_on = 0usize;

        for &(seq, coding) in sets {
            // Background: both strands, position-independent.
            accumulate_background(&seq.fwd, &mut nc);
            accumulate_background(&seq.rev, &mut nc);

            // Coding: accumulate codon-position-specific counts.
            let n = seq.len();
            for &(begin, end, strand) in coding {
                if begin == 0 || end > n || begin > end {
                    continue;
                }
                accumulate_coding(seq, n, begin, end, strand, &mut cod);
                trained_on += 1;
            }
        }

        // Normalise rows to log-probabilities.
        let mut m = MarkovModel {
            cod: [[[0.0; 4]; NCTX]; 3],
            nc: [[0.0; 4]; NCTX],
            trained_on,
        };
        for ctx in 0..NCTX {
            let s: f64 = nc[ctx].iter().sum();
            for b in 0..4 {
                m.nc[ctx][b] = (nc[ctx][b] / s).ln();
            }
            for p in 0..3 {
                let cs: f64 = cod[p][ctx].iter().sum();
                for b in 0..4 {
                    m.cod[p][ctx][b] = (cod[p][ctx][b] / cs).ln();
                }
            }
        }
        m
    }

    /// Coding/non-coding log-likelihood ratio for an ORF (global 1-based,
    /// inclusive). Higher = more coding.
    pub fn llr(&self, seq: &Seq, begin: usize, end: usize, strand: i8) -> f64 {
        let n = seq.len();
        if begin == 0 || end > n || begin > end {
            return f64::NEG_INFINITY;
        }
        let (arr, lo, hi) = strand_window(seq, n, begin, end, strand);
        // lo..=hi are 0-based indices into `arr`, read 5'->3'.
        let mut score = 0.0;
        let mut j = 0usize;
        let mut idx = lo;
        while idx <= hi {
            if j >= 2 {
                let ctx = (arr[idx - 2] as usize) * 4 + (arr[idx - 1] as usize);
                let b = arr[idx] as usize;
                let pos = j % 3;
                score += self.cod[pos][ctx][b] - self.nc[ctx][b];
            }
            j += 1;
            idx += 1;
        }
        score
    }
}

/// Return (strand array, lo, hi) 0-based inclusive read 5'->3' for a gene.
fn strand_window<'a>(seq: &'a Seq, n: usize, begin: usize, end: usize, strand: i8) -> (&'a [u8], usize, usize) {
    if strand == 1 {
        (&seq.fwd, begin - 1, end - 1)
    } else {
        // global [begin,end] → reverse-array indices [n-end, n-begin], 5' at n-end
        (&seq.rev, n - end, n - begin)
    }
}

fn accumulate_background(arr: &[u8], nc: &mut [[f64; 4]; NCTX]) {
    if arr.len() < 3 {
        return;
    }
    for i in 2..arr.len() {
        let ctx = (arr[i - 2] as usize) * 4 + (arr[i - 1] as usize);
        nc[ctx][arr[i] as usize] += 1.0;
    }
}

fn accumulate_coding(seq: &Seq, n: usize, begin: usize, end: usize, strand: i8, cod: &mut [[[f64; 4]; NCTX]; 3]) {
    let (arr, lo, hi) = strand_window(seq, n, begin, end, strand);
    let mut j = 0usize;
    let mut idx = lo;
    while idx <= hi {
        if j >= 2 {
            let ctx = (arr[idx - 2] as usize) * 4 + (arr[idx - 1] as usize);
            let b = arr[idx] as usize;
            cod[j % 3][ctx][b] += 1.0;
        }
        j += 1;
        idx += 1;
    }
}
