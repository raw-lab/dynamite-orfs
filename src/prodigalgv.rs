//! **Prodigal-GV** preset — phages, giant viruses (NCLDV / *Nucleocytoviricota*)
//! and crassphages.
//!
//! Wraps the shared Prodigal dynamic-programming core but auto-detects the
//! genetic code among the tables seen across large/▷giant viruses and
//! crassphages (which frequently reassign stop codons — e.g. code 15 reading
//! `TAG`→Q in some crassphages). Detection uses DYNAMITE's standard-code bias
//! so ordinary code-11 genomes are not spuriously reassigned.

use crate::gene::{Engine, Gene};
use crate::prodigal::dp_call;
use crate::seq::Seq;

/// Candidate genetic codes for giant viruses and crassphages (used when the
/// engine is forced with `-e prodigalgv`; the `virus` / `giant-virus` /
/// `crassphage` modes supply their own lists). Tables that cannot change the
/// gene calls relative to another candidate (1 ≡ 11 stops, 25 ≡ 4) are never
/// adopted on their own — see [`crate::prodigal::autodetect_code`].
pub const GV_CANDIDATES: &[u32] = &[11, 1, 4, 15, 25, 6];

/// Run the Prodigal-GV preset. `fixed = Some(table)` forces a code; `None`
/// auto-detects among [`GV_CANDIDATES`]. Treated as metagenome-like so multiple
/// contigs with mixed coding are handled gracefully.
pub fn call(contigs: &[(String, Seq)], fixed: Option<u32>, closed: bool) -> (Vec<Gene>, u32) {
    call_with_candidates(contigs, fixed, GV_CANDIDATES, closed)
}

/// [`call`] with an explicit candidate list for auto-detection (the per-mode
/// lists in [`crate::model::config_for`]).
pub fn call_with_candidates(
    contigs: &[(String, Seq)],
    fixed: Option<u32>,
    candidates: &[u32],
    closed: bool,
) -> (Vec<Gene>, u32) {
    dp_call(contigs, fixed, candidates, true, closed, Engine::ProdigalGv)
}
