//! Vendored core of **FragGeneScanRs** (Van der Jeugt, Dawyndt & Mesuere, 2022),
//! the maintained Rust reimplementation of FragGeneScan (Rho, Tang & Ye, 2010).
//!
//! Copied verbatim from FragGeneScanRs (GPL-3.0-or-later, compatible with
//! DYNAMITE's licence) so that DYNAMITE's `fraggenescan` engine runs the real
//! FragGeneScan HMM — full forward/Viterbi decode, frameshift (insertion/
//! deletion) handling, and the original trained models (embedded via
//! `include_bytes!`). Only the cross-module paths were rewritten; the algorithm
//! and model data are unchanged. Upstream: https://github.com/unipept/FragGeneScanRs
pub mod dna;
pub mod gene;
pub mod hmm;
pub mod viterbi;
