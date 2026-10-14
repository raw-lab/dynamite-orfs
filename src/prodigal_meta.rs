//! **Prodigal-meta** preset — fragmented, mixed metagenomic input.
//!
//! Wraps the shared Prodigal dynamic-programming core in metagenomic mode:
//! genes are allowed to run off contig ends (partials), and the genetic code is
//! auto-detected among the common prokaryotic/viral tables. Use this for
//! assembled contigs of unknown or mixed origin.

use crate::gene::{Engine, Gene};
use crate::prodigal::dp_call;
use crate::seq::Seq;

/// Candidate genetic codes for general metagenomes.
pub const META_CANDIDATES: &[u32] = &[11, 4, 15, 25];

/// Run the Prodigal-meta preset. `fixed = Some(table)` forces a code; `None`
/// auto-detects among [`META_CANDIDATES`].
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    call_with_candidates(contigs, fixed, META_CANDIDATES, closed)
}

/// [`call`] with an explicit candidate list for auto-detection.
pub fn call_with_candidates(
    contigs: &[(String, Seq)],
    fixed: Option<u32>,
    candidates: &[u32],
    closed: bool,
) -> (Vec<Gene>, u32) {
    dp_call(contigs, fixed, candidates, true, closed, Engine::ProdigalMeta)
}
