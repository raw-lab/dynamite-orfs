//! Vendored core of **FragGeneScanRs** (Van der Jeugt, Dawyndt & Mesuere, 2022),
//! the maintained Rust reimplementation of FragGeneScan (Rho, Tang & Ye, 2010).
//!
//! Copied verbatim from FragGeneScanRs so that DYNAMITE's `fraggenescan`
//! engine runs the real FragGeneScan HMM — full forward/Viterbi decode,
//! frameshift (insertion/deletion) handling, and the original trained models
//! (embedded via `include_bytes!`). Only the cross-module paths were
//! rewritten; the algorithm and model data are unchanged.
//! Upstream: https://github.com/unipept/FragGeneScanRs
//!
//! Licence: the FragGeneScanRs authors license this code and its models under
//! GPL-3.0-or-later, and this directory keeps that licence. The rest of
//! DYNAMITE is CC BY-NC 4.0 (see README, "License").
pub mod dna;
pub mod gene;
pub mod hmm;
pub mod viterbi;
