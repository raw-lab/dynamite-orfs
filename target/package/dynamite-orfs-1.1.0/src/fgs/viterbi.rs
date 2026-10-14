use strum::EnumCount;
use strum::IntoEnumIterator;

use crate::fgs::dna::Nuc::{A, C, G, T};
use crate::fgs::dna::{count_cg_content, trinucleotide, Nuc};
use crate::fgs::{gene, hmm};

pub fn viterbi(
    global: &hmm::Global,
    locals: &Vec<hmm::Local>,
    head: Vec<u8>,
    seq: Vec<Nuc>,
    whole_genome: bool,
) -> gene::ReadPrediction {
    let local = &locals[count_cg_content(&seq)];
    let (alpha, path) = forward(global, local, &seq, whole_genome);
    let vpath = backtrack(&alpha, path);
    build_genes(local, head, seq, whole_genome, vpath, alpha)
}

pub fn forward(
    global: &hmm::Global,
    local: &hmm::Local,
    seq: &Vec<Nuc>,
    whole_genome: bool,
) -> (
    Vec<[f64; hmm::State::COUNT]>,
    Vec<[Option<hmm::State>; hmm::State::COUNT]>,
) {
    let mut alpha: Vec<[f64; hmm::State::COUNT]> = Vec::with_capacity(seq.len());
    let mut path: Vec<[Option<hmm::State>; hmm::State::COUNT]> = Vec::with_capacity(seq.len());
    let mut temp_i: [usize; hmm::PERIOD] = [0; hmm::PERIOD];
    let mut temp_i_1: [usize; hmm::PERIOD] = [0; hmm::PERIOD];

    for _ in 0..seq.len() {
        alpha.push([0.0; hmm::State::COUNT]);
        path.push([Some(hmm::State::S); hmm::State::COUNT]);
    }
    alpha[0].copy_from_slice(&global.pi);
    for i in &mut alpha[0] {
        *i *= -1.0
    }

    // If the sequence starts with a stop codon
    if seq[0] == T
        && ((seq[1] == A && seq[2] == A)
            || (seq[1] == A && seq[2] == G)
            || (seq[1] == G && seq[2] == A))
    {
        alpha[0][hmm::State::E] = f64::INFINITY;
        alpha[1][hmm::State::E] = f64::INFINITY;
        path[1][hmm::State::E] = Some(hmm::State::E);
        path[2][hmm::State::E] = Some(hmm::State::E);

        alpha[2][hmm::State::M6] = f64::INFINITY;
        alpha[1][hmm::State::M5] = f64::INFINITY;
        alpha[0][hmm::State::M4] = f64::INFINITY;
        alpha[2][hmm::State::M3] = f64::INFINITY;
        alpha[1][hmm::State::M2] = f64::INFINITY;
        alpha[0][hmm::State::M1] = f64::INFINITY;

        alpha[2][hmm::State::E] -= if seq[1] == A && seq[2] == A {
            0.53_f64.ln()
        } else if seq[1] == A && seq[2] == G {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    // If the sequence starts with a reverse stop codon
    if seq[2] == A
        && ((seq[1] == T && seq[0] == T)
            || (seq[1] == T && seq[0] == C)
            || (seq[1] == C && seq[0] == T))
    {
        alpha[0][hmm::State::Sr] = f64::INFINITY;
        alpha[1][hmm::State::Sr] = f64::INFINITY;
        alpha[2][hmm::State::Sr] = alpha[0][hmm::State::S];
        path[1][hmm::State::Sr] = Some(hmm::State::Sr);
        path[2][hmm::State::Sr] = Some(hmm::State::Sr);

        alpha[2][hmm::State::M3r] = f64::INFINITY;
        alpha[2][hmm::State::M6r] = f64::INFINITY;

        alpha[2][hmm::State::Sr] = if seq[1] == T && seq[0] == T {
            0.53_f64.ln()
        } else if seq[1] == T && seq[0] == C {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    let mut num_noncoding = 0; // number of invalid nts in sequence
    for t in 1..seq.len() {
        let from = (seq[t - 1]).to_int().unwrap_or(2);
        let from0 = if t > 1 {
            (seq[t - 2]).to_int().unwrap_or(2)
        } else {
            2
        };
        let to = (seq[t]).to_int().unwrap_or_else(|| {
            num_noncoding += 1;
            2
        });
        let from2 = from0 * 4 + from;

        // M state
        if alpha[t][hmm::State::M1].is_finite() {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M6, hmm::State::M1, local.e_m[0][from2][to], global.tr.gg);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5, hmm::State::M1, 2.0, local.e_m[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4, hmm::State::M1, 3.0, local.e_m[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3, hmm::State::M1, 4.0, local.e_m[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2, hmm::State::M1, 5.0, local.e_m[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1, hmm::State::M1, 6.0, local.e_m[0][from2][to]);
            }
            from_s_to_m(&mut alpha, &mut path, local, t, from2, to);
            #[rustfmt::skip] from_i_to_m(&mut alpha, &mut path, &seq, temp_i[5], global, t, hmm::State::I6, hmm::State::M1);
        }

        if alpha[t][hmm::State::M2].is_finite() {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M1, hmm::State::M2, local.e_m[1][from2][to], 0.0);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6, hmm::State::M2, 2.0, local.e_m[1][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5, hmm::State::M2, 3.0, local.e_m[1][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4, hmm::State::M2, 4.0, local.e_m[1][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3, hmm::State::M2, 5.0, local.e_m[1][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2, hmm::State::M2, 6.0, local.e_m[1][from2][to]);
            }

            #[rustfmt::skip] from_i_to_m(&mut alpha, &mut path, &seq, temp_i[0], global, t, hmm::State::I1, hmm::State::M2);
        }

        if alpha[t][hmm::State::M3].is_finite() {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M2, hmm::State::M3, local.e_m[2][from2][to], 0.0);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6, hmm::State::M3, 3.0, local.e_m[2][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5, hmm::State::M3, 4.0, local.e_m[2][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4, hmm::State::M3, 5.0, local.e_m[2][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3, hmm::State::M3, 6.0, local.e_m[2][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1, hmm::State::M3, 2.0, local.e_m[2][from2][to]);
            }

            #[rustfmt::skip] from_i_to_m(&mut alpha, &mut path, &seq, temp_i[1], global, t, hmm::State::I2, hmm::State::M3);
        }

        if alpha[t][hmm::State::M4].is_finite() {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M3, hmm::State::M4, local.e_m[3][from2][to], 0.0);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6, hmm::State::M4, 4.0, local.e_m[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5, hmm::State::M4, 5.0, local.e_m[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4, hmm::State::M4, 6.0, local.e_m[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2, hmm::State::M4, 2.0, local.e_m[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1, hmm::State::M4, 3.0, local.e_m[3][from2][to]);
            }

            #[rustfmt::skip] from_i_to_m(&mut alpha, &mut path, &seq, temp_i[2], global, t, hmm::State::I3, hmm::State::M4);
        }

        if alpha[t][hmm::State::M5].is_finite() {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M4, hmm::State::M5, local.e_m[4][from2][to], 0.0);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6, hmm::State::M5, 5.0, local.e_m[4][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5, hmm::State::M5, 6.0, local.e_m[4][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3, hmm::State::M5, 2.0, local.e_m[4][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2, hmm::State::M5, 3.0, local.e_m[4][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1, hmm::State::M5, 4.0, local.e_m[4][from2][to]);
            }

            #[rustfmt::skip] from_i_to_m(&mut alpha, &mut path, &seq, temp_i[3], global, t, hmm::State::I4, hmm::State::M5);
        }

        if alpha[t][hmm::State::M6].is_finite() {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M5, hmm::State::M6, local.e_m[5][from2][to], 0.0);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6, hmm::State::M6, 6.0, local.e_m[5][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4, hmm::State::M6, 2.0, local.e_m[5][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3, hmm::State::M6, 3.0, local.e_m[5][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2, hmm::State::M6, 4.0, local.e_m[5][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1, hmm::State::M6, 5.0, local.e_m[5][from2][to]);
            }

            #[rustfmt::skip] from_i_to_m(&mut alpha, &mut path, &seq, temp_i[4], global, t, hmm::State::I5, hmm::State::M6);
        }

        // I state
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I1);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I2);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I3);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I4);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I5);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I6);
        #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i[0], global, t, from, to, hmm::State::M1, hmm::State::I1, 0.0);
        #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i[1], global, t, from, to, hmm::State::M2, hmm::State::I2, 0.0);
        #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i[2], global, t, from, to, hmm::State::M3, hmm::State::I3, 0.0);
        #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i[3], global, t, from, to, hmm::State::M4, hmm::State::I4, 0.0);
        #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i[4], global, t, from, to, hmm::State::M5, hmm::State::I5, 0.0);
        #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i[5], global, t, from, to, hmm::State::M6, hmm::State::I6, global.tr.gg);

        // M' state
        if t >= 3
            && seq[t - 1] == A
            && ((seq[t - 2] == T && seq[t - 3] == T)
                || (seq[t - 2] == T && seq[t - 3] == C)
                || (seq[t - 2] == C && seq[t - 3] == T))
        {
            #[rustfmt::skip] from_s_to_m1(&mut alpha, &mut path, t, hmm::State::M1r, local.e_m1[0][from2][to]);
        } else {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M6r, hmm::State::M1r, local.e_m1[0][from2][to], global.tr.gg);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5r, hmm::State::M1r, 2.0, local.e_m1[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4r, hmm::State::M1r, 3.0, local.e_m1[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3r, hmm::State::M1r, 4.0, local.e_m1[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2r, hmm::State::M1r, 5.0, local.e_m1[0][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1r, hmm::State::M1r, 6.0, local.e_m1[0][from2][to]);
            }
            #[rustfmt::skip] from_i1_to_m1(&mut alpha, &mut path, &seq, temp_i_1[5], global, t, hmm::State::I6r, hmm::State::M1r);
        }

        #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M1r, hmm::State::M2r, local.e_m1[1][from2][to], 0.0);
        if !whole_genome {
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6r, hmm::State::M2r, 2.0, local.e_m1[1][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5r, hmm::State::M2r, 3.0, local.e_m1[1][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4r, hmm::State::M2r, 4.0, local.e_m1[1][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3r, hmm::State::M2r, 5.0, local.e_m1[1][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2r, hmm::State::M2r, 6.0, local.e_m1[1][from2][to]);
        }
        #[rustfmt::skip] from_i1_to_m1(&mut alpha, &mut path, &seq, temp_i_1[0], global, t, hmm::State::I1r, hmm::State::M2r);

        #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M2r, hmm::State::M3r, local.e_m1[2][from2][to], 0.0);
        if !whole_genome {
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6r, hmm::State::M3r, 3.0, local.e_m1[2][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5r, hmm::State::M3r, 4.0, local.e_m1[2][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4r, hmm::State::M3r, 5.0, local.e_m1[2][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3r, hmm::State::M3r, 6.0, local.e_m1[2][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1r, hmm::State::M3r, 2.0, local.e_m1[2][from2][to]);
        }
        #[rustfmt::skip] from_i1_to_m1(&mut alpha, &mut path, &seq, temp_i_1[1], global, t, hmm::State::I2r, hmm::State::M3r);

        if t >= 3
            && seq[t - 1] == A
            && ((seq[t - 2] == T && seq[t - 3] == T)
                || (seq[t - 2] == T && seq[t - 3] == C)
                || (seq[t - 2] == C && seq[t - 3] == T))
        {
            #[rustfmt::skip] from_s_to_m1(&mut alpha, &mut path, t, hmm::State::M4r, local.e_m1[3][from2][to]);
        } else {
            #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M3r, hmm::State::M4r, local.e_m1[3][from2][to], 0.0);
            if !whole_genome {
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6r, hmm::State::M4r, 4.0, local.e_m1[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5r, hmm::State::M4r, 5.0, local.e_m1[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4r, hmm::State::M4r, 6.0, local.e_m1[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2r, hmm::State::M4r, 2.0, local.e_m1[3][from2][to]);
                #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1r, hmm::State::M4r, 3.0, local.e_m1[3][from2][to]);
            }
            #[rustfmt::skip] from_i1_to_m1(&mut alpha, &mut path, &seq, temp_i_1[2], global, t, hmm::State::I3r, hmm::State::M4r);
        }

        #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M4r, hmm::State::M5r, local.e_m1[4][from2][to], 0.0);
        if !whole_genome {
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6r, hmm::State::M5r, 5.0, local.e_m1[4][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M5r, hmm::State::M5r, 6.0, local.e_m1[4][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3r, hmm::State::M5r, 2.0, local.e_m1[4][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2r, hmm::State::M5r, 3.0, local.e_m1[4][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1r, hmm::State::M5r, 4.0, local.e_m1[4][from2][to]);
        }
        #[rustfmt::skip] from_i1_to_m1(&mut alpha, &mut path, &seq, temp_i_1[3], global, t, hmm::State::I4r, hmm::State::M5r);

        #[rustfmt::skip] from_m_to_m(&mut alpha, &mut path, global, t, hmm::State::M5r, hmm::State::M6r, local.e_m1[5][from2][to], 0.0);
        if !whole_genome {
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M6r, hmm::State::M6r, 6.0, local.e_m1[5][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M4r, hmm::State::M6r, 2.0, local.e_m1[5][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M3r, hmm::State::M6r, 3.0, local.e_m1[5][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M2r, hmm::State::M6r, 4.0, local.e_m1[5][from2][to]);
            #[rustfmt::skip] from_d_to_m(&mut alpha, &mut path, global, t, hmm::State::M1r, hmm::State::M6r, 5.0, local.e_m1[5][from2][to]);
        }
        #[rustfmt::skip] from_i1_to_m1(&mut alpha, &mut path, &seq, temp_i_1[4], global, t, hmm::State::I5r, hmm::State::M6r);

        // I' state
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I1r);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I2r);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I3r);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I4r);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I5r);
        from_i_to_i(&mut alpha, &mut path, global, t, from, to, hmm::State::I6r);

        if (t >= 3 && path[t - 3][hmm::State::Sr] != Some(hmm::State::R))
            && (t >= 4 && path[t - 4][hmm::State::Sr] != Some(hmm::State::R))
            && (t >= 5 && path[t - 5][hmm::State::Sr] != Some(hmm::State::R))
        {
            #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i_1[0], global, t, from, to, hmm::State::M1r, hmm::State::I1r, 0.0);
            #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i_1[1], global, t, from, to, hmm::State::M2r, hmm::State::I2r, 0.0);
            #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i_1[2], global, t, from, to, hmm::State::M3r, hmm::State::I3r, 0.0);
            #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i_1[3], global, t, from, to, hmm::State::M4r, hmm::State::I4r, 0.0);
            #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i_1[4], global, t, from, to, hmm::State::M5r, hmm::State::I5r, 0.0);
            #[rustfmt::skip] from_m_to_i(&mut alpha, &mut path, &mut temp_i_1[5], global, t, from, to, hmm::State::M6r, hmm::State::I6r, global.tr.gg);
        }

        // non_coding state
        from_r_to_r(&mut alpha, &mut path, global, local, t, from, to);
        from_e_to_r(&mut alpha, &mut path, global, t, hmm::State::E);
        from_e_to_r(&mut alpha, &mut path, global, t, hmm::State::Er);

        // end state
        if alpha[t][hmm::State::E] == 0.0 {
            alpha[t][hmm::State::E] = f64::INFINITY;
            path[t][hmm::State::E] = None;

            if t < seq.len() - 2
                && seq[t] == T
                && ((seq[t + 1] == A && seq[t + 2] == A)
                    || (seq[t + 1] == A && seq[t + 2] == G)
                    || (seq[t + 1] == G && seq[t + 2] == A))
            {
                alpha[t + 2][hmm::State::E] = f64::INFINITY;

                // transition from frame4, frame5 and frame6
                let temp_alpha = alpha[t - 1][hmm::State::M6] - global.tr.ge;
                if temp_alpha < alpha[t + 2][hmm::State::E] {
                    alpha[t + 2][hmm::State::E] = temp_alpha;
                    path[t][hmm::State::E] = Some(hmm::State::M6);
                }

                // transition from frame1, frame2 and frame3
                let temp_alpha = alpha[t - 1][hmm::State::M3] - global.tr.ge;
                if temp_alpha < alpha[t + 2][hmm::State::E] {
                    alpha[t + 2][hmm::State::E] = temp_alpha;
                    path[t][hmm::State::E] = Some(hmm::State::M3);
                }

                alpha[t][hmm::State::E] = f64::INFINITY;
                alpha[t + 1][hmm::State::E] = f64::INFINITY;
                path[t + 1][hmm::State::E] = Some(hmm::State::E);
                path[t + 2][hmm::State::E] = Some(hmm::State::E);

                alpha[t + 2][hmm::State::M6] = f64::INFINITY;
                alpha[t + 1][hmm::State::M5] = f64::INFINITY;
                alpha[t][hmm::State::M4] = f64::INFINITY;
                alpha[t + 2][hmm::State::M3] = f64::INFINITY;
                alpha[t + 1][hmm::State::M2] = f64::INFINITY;
                alpha[t][hmm::State::M1] = f64::INFINITY;

                alpha[t + 2][hmm::State::E] -= if seq[t + 1] == A && seq[t + 2] == A {
                    0.54_f64.ln()
                } else if seq[t + 1] == A && seq[t + 2] == G {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                for i in (t.max(60) - 60)..=(t.max(3) - 3) {
                    start_freq -=
                        local.tr_e[i + 60 - t][trinucleotide(seq.get(i..).unwrap()).unwrap_or(0)];
                }
                if t < 60 {
                    start_freq *= 58.0 / (t.saturating_sub(3) + 1) as f64
                }
                modify_border_dist(&mut alpha[t + 2][hmm::State::E], &local.dist_e, start_freq);
            }
        }

        // start' state
        // originally stop codon of genes in - strand
        if alpha[t][hmm::State::Sr] == 0.0 {
            alpha[t][hmm::State::Sr] = f64::INFINITY;
            path[t][hmm::State::Sr] = None;

            if t < seq.len() - 2
                && seq[t + 2] == A
                && ((seq[t + 1] == T && seq[t] == T)
                    || (seq[t + 1] == T && seq[t] == C)
                    || (seq[t + 1] == C && seq[t] == T))
            {
                alpha[t][hmm::State::Sr] = f64::INFINITY;
                alpha[t + 1][hmm::State::Sr] = f64::INFINITY;
                alpha[t + 2][hmm::State::Sr] = alpha[t - 1][hmm::State::R] - global.tr.rs;
                path[t][hmm::State::Sr] = Some(hmm::State::R);
                path[t + 1][hmm::State::Sr] = Some(hmm::State::Sr);
                path[t + 2][hmm::State::Sr] = Some(hmm::State::Sr);

                let temp_alpha = alpha[t - 1][hmm::State::Er] - global.tr.es;
                if temp_alpha < alpha[t + 2][hmm::State::Sr] {
                    alpha[t + 2][hmm::State::Sr] = temp_alpha;
                    path[t][hmm::State::Sr] = Some(hmm::State::Er);
                }

                let temp_alpha = alpha[t - 1][hmm::State::E] - global.tr.es1;
                if temp_alpha < alpha[t + 2][hmm::State::Sr] {
                    alpha[t + 2][hmm::State::Sr] = temp_alpha;
                    path[t][hmm::State::Sr] = Some(hmm::State::E);
                }

                alpha[t + 2][hmm::State::M3r] = f64::INFINITY;
                alpha[t + 2][hmm::State::M6r] = f64::INFINITY;

                alpha[t + 2][hmm::State::Sr] -= if seq[t + 1] == T && seq[t] == T {
                    0.54_f64.ln()
                } else if seq[t + 1] == T && seq[t] == C {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                if t + 5 < seq.len() {
                    for i in (t + 3)..=(t + 60).min(seq.len() - 3) {
                        start_freq -= local.tr_s1[i - 3 - t]
                            [trinucleotide(seq.get(i..).unwrap()).unwrap_or(0)];
                    }
                }
                // TODO add similar limit to other 3 ends? proposal:
                //if t + 63 > seq.len() {
                //    start_freq *= 58.0 / (seq.len() - t - 5) as f64
                //}

                modify_border_dist(
                    &mut alpha[t + 2][hmm::State::Sr],
                    &local.dist_s1,
                    start_freq,
                );
            }
        }

        // start state
        if alpha[t][hmm::State::S] == 0.0 {
            alpha[t][hmm::State::S] = f64::INFINITY;
            path[t][hmm::State::S] = None;

            if t < seq.len() - 2
                && (seq[t] == A || seq[t] == G || seq[t] == T)
                && seq[t + 1] == T
                && seq[t + 2] == G
            {
                alpha[t][hmm::State::S] = f64::INFINITY;
                alpha[t + 1][hmm::State::S] = f64::INFINITY;
                alpha[t + 2][hmm::State::S] = alpha[t - 1][hmm::State::R] - global.tr.rs;
                path[t][hmm::State::S] = Some(hmm::State::R);
                path[t + 1][hmm::State::S] = Some(hmm::State::S);
                path[t + 2][hmm::State::S] = Some(hmm::State::S);

                let temp_alpha = alpha[t - 1][hmm::State::E] - global.tr.es;
                if temp_alpha < alpha[t + 2][hmm::State::S] {
                    alpha[t + 2][hmm::State::S] = temp_alpha;
                    path[t][hmm::State::S] = Some(hmm::State::E);
                }

                let temp_alpha = alpha[t - 1][hmm::State::Er] - global.tr.es1;
                if temp_alpha < alpha[t + 2][hmm::State::S] {
                    alpha[t + 2][hmm::State::S] = temp_alpha;
                    path[t][hmm::State::S] = Some(hmm::State::Er);
                }

                alpha[t + 2][hmm::State::S] -= if seq[t] == A {
                    0.83_f64.ln()
                } else if seq[t] == G {
                    0.10_f64.ln()
                } else {
                    0.07_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                for i in (t.max(30) - 30)..=(t + 30).min(seq.len() - 3) {
                    start_freq -=
                        local.tr_s[i + 30 - t][trinucleotide(seq.get(i..).unwrap()).unwrap_or(0)];
                }
                if t < 30 {
                    start_freq *= 61.0 / (t + 30 + 1) as f64;
                }
                modify_border_dist(&mut alpha[t + 2][hmm::State::S], &local.dist_s, start_freq);
            }
        }

        // end' state
        // originally start codon of genes in - strand
        if alpha[t][hmm::State::Er] == 0.0 {
            alpha[t][hmm::State::Er] = f64::INFINITY;
            path[t][hmm::State::Er] = None;

            if t < seq.len() - 2
                && seq[t] == C
                && seq[t + 1] == A
                && (seq[t + 2] == T || seq[t + 2] == C || seq[t + 2] == A)
            {
                // transition from frame6
                alpha[t][hmm::State::Er] = f64::INFINITY;
                alpha[t + 1][hmm::State::Er] = f64::INFINITY;
                alpha[t + 2][hmm::State::Er] = alpha[t - 1][hmm::State::M6r] - global.tr.ge;
                path[t][hmm::State::Er] = Some(hmm::State::M6r);
                path[t + 1][hmm::State::Er] = Some(hmm::State::Er);
                path[t + 2][hmm::State::Er] = Some(hmm::State::Er);

                alpha[t + 2][hmm::State::Er] -= if seq[t + 2] == T {
                    0.83_f64.ln()
                } else if seq[t + 2] == C {
                    0.10_f64.ln()
                } else {
                    0.07_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                for i in (t.max(30) - 30)..=(t + 30).min(seq.len() - 3) {
                    start_freq -=
                        local.tr_e1[i + 30 - t][trinucleotide(seq.get(i..).unwrap()).unwrap_or(0)];
                }
                if t < 30 {
                    start_freq *= 61.0 / (t + 30 + 1) as f64;
                }

                modify_border_dist(
                    &mut alpha[t + 2][hmm::State::Er],
                    &local.dist_e1,
                    start_freq,
                );
            }
        }

        if num_noncoding > 9 {
            for i in hmm::State::iter() {
                if i != hmm::State::R {
                    alpha[t][i] = f64::INFINITY;
                    path[t][i] = Some(hmm::State::R);
                }
            }
        }
    }

    (alpha, path)
}

fn backtrack(
    alpha: &Vec<[f64; hmm::State::COUNT]>,
    path: Vec<[Option<hmm::State>; hmm::State::COUNT]>,
) -> Vec<hmm::State> {
    // backtrack array to find the optimal path
    let mut vpath: Vec<hmm::State> = Vec::with_capacity(path.len());
    vpath.push(hmm::State::S); // or null
    let mut prob = f64::INFINITY;
    for (&prob_, i) in alpha
        .last()
        .expect("empty seq")
        .iter()
        .zip(hmm::State::iter())
    {
        if prob_ < prob {
            vpath[0] = i;
            prob = prob_;
        }
    }

    // backtrack the optimal path
    for t in (0..=path.len() - 2).rev() {
        vpath.push(path[t + 1][*vpath.last().unwrap()].unwrap_or_else(|| {
            eprintln!(
                "Warning: encountered None-state in chosen path, replacing with non-coding state"
            );
            hmm::State::R
        }));
    }
    vpath.reverse();
    vpath
}

fn build_genes(
    local: &hmm::Local,
    head: Vec<u8>,
    seq: Vec<Nuc>,
    whole_genome: bool,
    vpath: Vec<hmm::State>,
    alpha: Vec<[f64; hmm::State::COUNT]>,
) -> gene::ReadPrediction {
    let gene_len = if whole_genome { 120 } else { 60 }; // minimum length to be output
    let mut read_prediction = gene::ReadPrediction::new(head);
    let mut codon_start = 0; // ternaire boolean?
    let mut start_t: isize = -1;
    let mut dna_start_t_withstop: usize = 0;
    let mut dna_start_t: usize = 0;

    let mut dna: Vec<Nuc> = Vec::with_capacity(seq.len());
    let mut insert = vec![];
    let mut delete = vec![];

    let mut prev_match = hmm::State::S; // TODO or no state

    let mut start_orf = 0; // initialize?

    for t in 0..seq.len() {
        if codon_start == 0
            && start_t < 0
            && ((vpath[t] >= hmm::State::M1 && vpath[t] <= hmm::State::M6)
                || (vpath[t] >= hmm::State::M1r && vpath[t] <= hmm::State::M6r)
                || vpath[t] == hmm::State::S
                || vpath[t] == hmm::State::Sr)
        {
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            start_t = (t + 1) as isize;
            // introduce dna_start_t_withstop YY July 2018
        }

        if codon_start == 0
            && (vpath[t] == hmm::State::M1
                || vpath[t] == hmm::State::M4
                || vpath[t] == hmm::State::M1r
                || vpath[t] == hmm::State::M4r)
        {
            dna.clear();
            insert.clear();
            delete.clear();

            dna.push(seq[t]);
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            if vpath[t] == hmm::State::M1r || vpath[t] == hmm::State::M4r {
                if t > 2 {
                    dna_start_t_withstop = t - 2;
                    dna.splice(0..0, seq[t - 3..t].iter().cloned());
                }
            }

            start_orf = t + 1;
            prev_match = vpath[t];

            codon_start = if vpath[t] < hmm::State::M6 { 1 } else { -1 }
        } else if codon_start != 0
            && (vpath[t] == hmm::State::E || vpath[t] == hmm::State::Er || t == seq.len() - 1)
        {
            let mut end_t;
            if vpath[t] == hmm::State::E || vpath[t] == hmm::State::Er {
                end_t = t + 3
            } else {
                let mut temp_t = t;
                while vpath[temp_t] != hmm::State::M1
                    && vpath[temp_t] != hmm::State::M4
                    && vpath[temp_t] != hmm::State::M1r
                    && vpath[temp_t] != hmm::State::M4r
                {
                    dna.pop();
                    temp_t -= 1;
                }
                end_t = temp_t;
            }

            if dna.len()
                > gene_len
                    + if whole_genome && codon_start == -1 {
                        3
                    } else {
                        0
                    }
            {
                let final_score = (alpha[end_t - 4][vpath[end_t - 4]]
                    - alpha[start_t as usize + 2][vpath[start_t as usize + 2]])
                    / (end_t - start_t as usize - 5) as f64;
                let mut frame = start_orf % 3;
                if frame == 0 {
                    frame = 3
                }

                if codon_start == 1 {
                    if start_t == dna_start_t as isize - 3 {
                        dna_start_t -= 3;
                        dna.splice(0..0, seq[dna_start_t - 1..dna_start_t + 2].iter().cloned());
                    }

                    if whole_genome {
                        // add refinement of the start codons here
                        let start_old = start_t as usize;
                        let mut codon = &seq[start_old - 1..start_old + 2];
                        let mut s = 0;
                        // find the optimal start codon within 30bp up- and downstream of start codon
                        let mut e_save = 0.0;
                        let mut s_save = 0;
                        while !(codon == [T, A, A] || codon == [T, A, G] || codon == [T, G, A])
                            && start_old >= 1 + s + 35
                        {
                            if codon == [A, T, G] || codon == [G, T, G] || codon == [T, T, G] {
                                let utr = &seq[start_old - 1 - s - 30..start_old - 1 - s - 30 + 63];
                                let mut freq_sum = 0.0;
                                for j in 0..utr.len() - 2 {
                                    freq_sum -= local.tr_s[j]
                                        [trinucleotide(utr.get(j..).unwrap()).unwrap_or(0)];
                                }
                                if s == 0 {
                                    e_save = freq_sum;
                                    s_save = 0;
                                } else if freq_sum < e_save {
                                    e_save = freq_sum;
                                    s_save = s; // posivite chain, upstream s_save = -1 * 3
                                }
                            }
                            s += 3;
                            codon = &seq[start_old - 1 - s..start_old - 1 - s + 3];
                        }
                        dna_start_t -= s_save;
                        dna.splice(
                            0..0,
                            seq[dna_start_t - 1..dna_start_t + s_save - 1]
                                .iter()
                                .cloned(),
                        );

                        // add final codon
                        dna.push(seq[end_t - 3]);
                        dna.push(seq[end_t - 2]);
                        dna.push(seq[end_t - 1]);
                    }

                    read_prediction.genes.push(gene::Gene {
                        start: dna_start_t,
                        end: end_t,
                        frame: frame,
                        score: final_score,
                        dna: dna.clone(),
                        forward_strand: true,
                        inserted: insert.clone(),
                        deleted: delete.clone(),
                    });
                } else if codon_start == -1 {
                    if whole_genome {
                        // add refinement of the start codons here
                        // l 946
                        let end_old = end_t;
                        let mut codon = &seq[end_t - 1 - 2..end_t];
                        let mut s = 0;
                        // find the optimal start codon within 30bp up- and downstream of start codon
                        let mut e_save = 0.0;
                        let mut s_save = 0;
                        while !(codon == [T, T, A] || codon == [C, T, A] || codon == [T, C, A])
                            && end_old - 2 + s + 35 < seq.len()
                        {
                            if codon == [C, A, T] || codon == [C, A, C] || codon == [C, A, A] {
                                let utr = &seq[end_old - 1 - 2 + s - 30..end_old + s + 30];
                                let mut freq_sum = 0.0;
                                for j in 0..utr.len() - 2 {
                                    // TODO stop1? (their note)
                                    freq_sum -= local.tr_e1[j]
                                        [trinucleotide(utr.get(j..).unwrap()).unwrap_or(0)];
                                }
                                if s == 0 || freq_sum < e_save {
                                    e_save = freq_sum;
                                    s_save = s; // negative chain, s_save = s
                                }
                            }
                            s += 3;
                            codon = &seq[end_old - 1 - 2 + s..end_old + s];
                        }

                        end_t = end_old + s_save;

                        // add final codon
                        dna.push(seq[end_t - 3]);
                        dna.push(seq[end_t - 2]);
                        dna.push(seq[end_t - 1]);
                    }

                    read_prediction.genes.push(gene::Gene {
                        start: dna_start_t_withstop,
                        end: end_t,
                        frame: frame,
                        score: final_score,
                        dna: dna.clone(),
                        forward_strand: false,
                        inserted: insert.clone(),
                        deleted: delete.clone(),
                    });
                }
            }

            codon_start = 0;
            start_t = -1;
            dna.clear();
        } else if codon_start != 0
            && ((vpath[t] >= hmm::State::M1
                && vpath[t] <= hmm::State::M6
                && prev_match >= hmm::State::M1)
                || (vpath[t] >= hmm::State::M1r
                    && vpath[t] <= hmm::State::M6r
                    && prev_match >= hmm::State::M1r))
        {
            let out_nt = if vpath[t] < prev_match {
                vpath[t] as usize + 6 - prev_match as usize
            } else {
                vpath[t] as usize - prev_match as usize
            };
            for kk in 0..out_nt {
                // for deleted nt in reads
                dna.push(Nuc::Xi);
                if kk > 0 {
                    delete.push(t + 1);
                }
            }
            dna.pop();
            dna.push(seq[t]);
            prev_match = vpath[t];
        } else if codon_start != 0
            && ((vpath[t] >= hmm::State::I1 && vpath[t] <= hmm::State::I6)
                || (vpath[t] >= hmm::State::I1r && vpath[t] <= hmm::State::I6r))
        {
            dna.push(seq[t].to_lower());
            insert.push(t + 1);
        } else if codon_start != 0 && vpath[t] == hmm::State::R {
            // for long NNNNNNNN, pretend R state
            codon_start = 0;
            start_t = -1;
        }
    }

    read_prediction
}

#[inline]
fn from_m_to_m(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    global: &hmm::Global,
    t: usize,
    from_m: hmm::State,
    to_m: hmm::State,
    emission: f64,
    last_m: f64,
) {
    alpha[t][to_m] = alpha[t - 1][from_m] - last_m - global.tr.mm - emission;
    path[t][to_m] = Some(from_m);
}

#[inline]
fn from_d_to_m(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    global: &hmm::Global,
    t: usize,
    from_m: hmm::State,
    to_m: hmm::State,
    num_d: f64,
    emission: f64,
) {
    if num_d > 0.0 {
        let temp_alpha = alpha[t - 1][from_m]
            - global.tr.md
            - emission
            - 0.25_f64.ln() * (num_d - 1.0)
            - global.tr.dd * (num_d - 2.0)
            - global.tr.dm;
        if temp_alpha < alpha[t][to_m] {
            alpha[t][to_m] = temp_alpha;
            path[t][to_m] = Some(from_m);
        }
    }
}

#[inline]
fn from_s_to_m(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    local: &hmm::Local,
    t: usize,
    from2: usize,
    to: usize,
) {
    let temp_alpha = alpha[t - 1][hmm::State::S] - local.e_m[0][from2][to];
    if temp_alpha <= alpha[t][hmm::State::M1] {
        alpha[t][hmm::State::M1] = temp_alpha;
        path[t][hmm::State::M1] = Some(hmm::State::S);
    }
}

#[inline]
fn from_s_to_m1(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    t: usize,
    to_m: hmm::State,
    emission: f64,
) {
    // from Start state since this is actually a stop codon in minus strand
    alpha[t][to_m] = alpha[t - 1][hmm::State::Sr] - emission;
    path[t][to_m] = Some(hmm::State::Sr);
}

#[inline]
fn from_i_to_m(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    seq: &Vec<Nuc>,
    temp_i: usize,
    global: &hmm::Global,
    t: usize,
    from_i: hmm::State,
    to_m: hmm::State,
) {
    // to avoid stop codon
    if t < 2 {
    } else if (to_m == hmm::State::M2 || to_m == hmm::State::M5)
        && t + 1 < seq.len()
        && seq[temp_i] == T
        && ((seq[t] == A && seq[t + 1] == A)
            || (seq[t] == A && seq[t + 1] == G)
            || (seq[t] == G && seq[t + 1] == A))
    {
    } else if (to_m == hmm::State::M3 || to_m == hmm::State::M6)
        && temp_i > 0
        && (seq[temp_i - 1] == T)
        && ((seq[temp_i] == A && seq[t] == A)
            || (seq[temp_i] == A && seq[t] == G)
            || (seq[temp_i] == G && seq[t] == A))
    {
    } else {
        let temp_alpha = alpha[t - 1][from_i] - global.tr.im - 0.25_f64.ln();
        if temp_alpha < alpha[t][to_m] {
            alpha[t][to_m] = temp_alpha;
            path[t][to_m] = Some(from_i);
        }
    }
}

#[inline]
fn from_i_to_i(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    global: &hmm::Global,
    t: usize,
    from: usize,
    to: usize,
    i: hmm::State,
) {
    alpha[t][i] = alpha[t - 1][i] - global.tr.ii - global.tr_ii[from][to];
    path[t][i] = Some(i);
}

#[inline]
fn from_m_to_i(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    temp_i: &mut usize,
    global: &hmm::Global,
    t: usize,
    from: usize,
    to: usize,
    from_m: hmm::State,
    to_i: hmm::State,
    last_i: f64,
) {
    let temp_alpha = alpha[t - 1][from_m] - global.tr.mi - global.tr_mi[from][to] - last_i;
    if temp_alpha < alpha[t][to_i] {
        alpha[t][to_i] = temp_alpha;
        path[t][to_i] = Some(from_m);
        *temp_i = t - 1;
    }
}

#[inline]
fn from_i1_to_m1(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    seq: &Vec<Nuc>,
    temp_i_1: usize,
    global: &hmm::Global,
    t: usize,
    from_i: hmm::State,
    to_m: hmm::State,
) {
    // to avoid stop codon
    if t < 2 {
    } else if (to_m == hmm::State::M2r || to_m == hmm::State::M5r)
        && t + 1 < seq.len()
        && seq[t + 1] == A
        && ((seq[t] == T && seq[temp_i_1] == T)
            || (seq[t] == T && seq[temp_i_1] == C)
            || (seq[t] == A && seq[temp_i_1] == T))
    {
    } else if (to_m == hmm::State::M3r || to_m == hmm::State::M6r)
        && seq[t] == A
        && temp_i_1 > 1
        && ((seq[temp_i_1] == T && seq[temp_i_1 - 1] == T)
            || (seq[temp_i_1] == T && seq[temp_i_1 - 1] == C)
            || (seq[temp_i_1] == C && seq[temp_i_1 - 1] == T))
    {
    } else {
        let temp_alpha = alpha[t - 1][from_i] - global.tr.im - 0.25_f64.ln();
        if temp_alpha < alpha[t][to_m] {
            alpha[t][to_m] = temp_alpha;
            path[t][to_m] = Some(from_i);
        }
    }
}

#[inline]
fn from_r_to_r(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    global: &hmm::Global,
    local: &hmm::Local,
    t: usize,
    from: usize,
    to: usize,
) {
    alpha[t][hmm::State::R] =
        alpha[t - 1][hmm::State::R] - local.tr_rr[from][to] - global.tr.rr - 0.95_f64.ln();
    path[t][hmm::State::R] = Some(hmm::State::R);
}

#[inline]
fn from_e_to_r(
    alpha: &mut Vec<[f64; hmm::State::COUNT]>,
    path: &mut Vec<[Option<hmm::State>; hmm::State::COUNT]>,
    global: &hmm::Global,
    t: usize,
    from_e: hmm::State,
) {
    let temp_alpha = alpha[t - 1][from_e] - global.tr.er - 0.95_f64.ln();
    if temp_alpha < alpha[t][hmm::State::R] {
        alpha[t][hmm::State::R] = temp_alpha;
        path[t][hmm::State::R] = Some(from_e);
    }
}

fn modify_border_dist(cell: &mut f64, values: &[f64], start_freq: f64) {
    let h_kd =
        values[2] * (-1.0 * (start_freq - values[1]).powi(2) / (values[0]).powi(2) / 2.0).exp();
    let r_kd =
        values[5] * (-1.0 * (start_freq - values[4]).powi(2) / (values[3]).powi(2) / 2.0).exp();
    *cell -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use crate::fgs::dna::dna;

    fn test_viterbi(file: &str, seq: Vec<Nuc>, genes: Vec<gene::Gene>) {
        let (global, locals) =
            hmm::get_train_from_file(PathBuf::from("train"), PathBuf::from(file)).unwrap();
        let received = viterbi(&global, &locals, vec![], seq, false);
        let expected = gene::ReadPrediction {
            head: vec![],
            genes: genes,
        };
        assert_eq!(expected, received);
        // assert_eq!(expected, received)
        // assert_matches!(received, pattern)
    }

    #[test]
    fn test_forward_read() {
        test_viterbi(
            "454_10",
            dna("TGTTCGCTGGCGGTGCTTTAGGGGGAGATGCGCAGAATGTCTCAGCCCTGCAAGGGGGTACTCTAGACATGGTGGTATTGAATTCCGGTATCTTAGCTTCGCAAGTGAAAGAGTTTGCCGTGTATGACTTTCCGTTCATGTTCAATAACTCTGAAGAAGCCGATAAAGTCGTCGACGGTGAATTTGGTAAGGCTTTGCACGCGAAATTAGAAGAGAAAGGTATTATCGGTTTAGCCTATTGGGAGTTGGGCTTCCGCGATATGACTAATAGCAAACACCCTATCACTAAGGTCGAGGATATTAAAGGTCTGAAACTGCGTGTTATCCC"),
            vec![
                gene::Gene {
                    start: 3,
                    end: 326,
                    frame: 3,
                    score: 1.310173669505397,
                    dna: dna("TTCGCTGGCGGTGCTTTAGGGGGAGATGCGCAGAATGTCTCAGCCCTGCAAGGGGGTACTCTAGACATGGTGGTATTGAATTCCGGTATCTTAGCTTCGCAAGTGAAAGAGTTTGCCGTGTATGACTTTCCGTTCATGTTCAATAACTCTGAAGAAGCCGATAAAGTCGTCGACGGTGAATTTGGTAAGGCTTTGCACGCGAAATTAGAAGAGAAAGGTATTATCGGTTTAGCCTATTGGGAGTTGGGCTTCCGCGATATGACTAATAGCAAACACCCTATCACTAAGGTCGAGGATATTAAAGGTCTGAAACTGCGTGTTATC"),
                    forward_strand: true,
                    inserted: vec![],
                    deleted: vec![],
                }
            ]
        );
    }

    #[test]
    fn test_reverse_read() {
        test_viterbi(
            "454_10",
            dna("GGGATAACACGCAGTTTCAGACCTTTAATATCCTCGACCTTAGTGATAGGGTGTTTGCTATTAGTCATATCGCGGAAGCCCAACTCCCAATAGGCTAAACCGATAATACCTTTCTCTTCTAATTTCGCGTGCAAAGCCTTACCAAATTCACCGTCGACGACTTTATCGGCTTCTTCAGAGTTATTGAACATGAACGGAAAGTCATACACGGCAAACTCTTTCACTTGCGAAGCTAAGATACCGGAATTCAATACCACCATGTCTAGAGTACCCCCTTGCAGGGCTGAGACATTCTGCGCATCTCCCCCTAAAGCACCGCCAGCGAACA"),
            vec![
                gene::Gene {
                    start: 3,
                    end: 326,
                    frame: 3,
                    score: 1.3039173404556075,
                    dna: dna("GATAACACGCAGTTTCAGACCTTTAATATCCTCGACCTTAGTGATAGGGTGTTTGCTATTAGTCATATCGCGGAAGCCCAACTCCCAATAGGCTAAACCGATAATACCTTTCTCTTCTAATTTCGCGTGCAAAGCCTTACCAAATTCACCGTCGACGACTTTATCGGCTTCTTCAGAGTTATTGAACATGAACGGAAAGTCATACACGGCAAACTCTTTCACTTGCGAAGCTAAGATACCGGAATTCAATACCACCATGTCTAGAGTACCCCCTTGCAGGGCTGAGACATTCTGCGCATCTCCCCCTAAAGCACCGCCAGCGAA"),
                    forward_strand: false,
                    inserted: vec![],
                    deleted: vec![],
                }
            ]
        );
    }

    #[test]
    fn test_deletion() {
        test_viterbi(
            "454_10",
            dna("GTCGACAGTGTAGTAACCAGTGCTCACGATACCATTGTGGGATCAGCGACCAGAGTTGCTGCAACATTTCACCGCTGGTAACAACGACCATCG"),
            vec![
                gene::Gene {
                    start: 1,
                    end: 91,
                    frame: 1,
                    score: 1.3668123444828921,
                    dna: dna("GTCGACAGTGTAGTAACCAGTGCTCACGATACCATTGtGGGATCAGCGACCAGAGTTGCTGCAACATTTCACCGCTGGTAACAACGACCAT"),
                    forward_strand: false,
                    inserted: vec![38],
                    deleted: vec![],
                }
            ]
        );
    }
}
