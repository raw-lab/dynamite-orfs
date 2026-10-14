//! # DYNAMITE core
//!
//! `dynamite_core` is the engine library behind the DYNAMITE gene caller — an
//! ORF / gene caller for **all domains of life**. It exposes a family of
//! callers that share a common sequence model, genetic-code table, and
//! [`gene::Gene`] output type:
//!
//! * [`phanotate`] — PHANOTATE-style graph / shortest-path caller (phages).
//! * [`prodigal`] — Prodigal dynamic-programming engine; the `prodigal`,
//!   [`prodigalgv`], and [`prodigal_meta`] presets configure it for bacteria /
//!   archaea, phages + giant viruses / crassphages, and metagenomes.
//! * [`genemark`] — self-training 3-periodic Markov caller (GeneMark-style).
//! * [`augustus`] — eukaryotic single-exon caller (AUGUSTUS-style).
//! * [`fraggenescan`] — read caller (the vendored FragGeneScanRs HMM).
//!
//! Plus helpers: [`sixframe`] (six-frame translation), [`codon`] (codon-bias /
//! density statistics), [`rna`] (heuristic tRNA / rRNA detection), and
//! [`doctor`] (self-test). [`run_records`] ties the callers together, resolving
//! the engine and genetic code from a [`model::Mode`] (or an explicit engine
//! override) and predicting every contig into one gene list.
//!
//! Every engine's calls pass through [`finalize`]: contigs are split at
//! assembly gaps (runs of N), terminal gap codons are trimmed, and partial
//! ends are inferred from the start / stop codons and reported in Prodigal's
//! genomic `partial=XY` (left/right) convention on both strands.

pub mod augustus;
pub mod banner;

/// DYNAMITE's licence, as an SPDX identifier: Creative Commons
/// Attribution-NonCommercial 4.0 International (CC BY-NC 4.0). See
/// `LICENSE.txt`; bundled third-party code keeps its own licence (README,
/// "License").
pub const LICENSE: &str = "CC-BY-NC-4.0";

const fn same_str(a: &str, b: &str) -> bool {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    if a.len() != b.len() {
        return false;
    }
    let mut i = 0;
    while i < a.len() {
        if a[i] != b[i] {
            return false;
        }
        i += 1;
    }
    true
}

// Licence guard: the crate does not compile (so it cannot be built, tested or
// published) unless Cargo.toml declares exactly `license = "CC-BY-NC-4.0"`.
const _: () = assert!(
    same_str(env!("CARGO_PKG_LICENSE"), LICENSE),
    "Cargo.toml must declare license = \"CC-BY-NC-4.0\": DYNAMITE is CC BY-NC 4.0 throughout"
);
pub mod codon;
pub mod dedup;
pub mod doctor;
pub mod fasta;
pub mod annotate;
pub mod fgs;
pub mod finalize;
pub mod fraggenescan;
pub mod gcframe;
pub mod gencode;
pub mod gene;
pub mod genemark;
pub mod markov;
pub mod model;
pub mod node;
pub mod orf;
pub mod phanotate;
pub mod prodigal;
pub mod prodigal_meta;
pub mod prodigalgv;
pub mod rbs;
pub mod rna;
pub mod seq;
pub mod sixframe;

use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::gencode::GeneticCode;
use crate::gene::{Engine, Gene};
use crate::model::{config_for, CodeChoice, EngineKind, Mode};
use crate::seq::Seq;

/// Everything needed to run a prediction over a set of records.
#[derive(Clone, Debug)]
pub struct RunConfig {
    /// Organism preset.
    pub mode: Mode,
    /// Explicit engine override (takes precedence over the mode's default).
    pub engine: Option<EngineKind>,
    /// Genetic-code choice (overrides the mode's default when not `Default`).
    pub code: CodeChoice,
    /// Minimum ORF length in nucleotides (PHANOTATE engine). `<= 0` uses the
    /// engine default.
    pub min_len: i64,
    /// Disallow genes that run off the ends of a contig.
    pub closed: bool,
    /// Worker threads (`0` = all available).
    pub threads: usize,
}

impl Default for RunConfig {
    fn default() -> Self {
        RunConfig {
            mode: Mode::Auto,
            engine: None,
            code: CodeChoice::Default,
            min_len: 0,
            closed: false,
            threads: 0,
        }
    }
}

/// How runs of ambiguous bases (assembly gaps) are handled by the Prodigal
/// engines (`prodigal`, `prodigalgv`, `prodigal-meta`): Prodigal's gap
/// options (`-e` of the Prodigal 3 development version, `-m` of Prodigal 2.x).
/// The other engines always split contigs at gaps ([`GapMode::Partial`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GapMode {
    /// `--gap-mode 0` — the default (Prodigal 3 `-e 0`): runs of at least
    /// [`seq::GAP_MIN_RUN`] ambiguous bases are gaps. Contigs are split
    /// there, so no gene spans one, and genes run into a gap as partial genes,
    /// exactly as into a contig end.
    #[default]
    Partial,
    /// `--gap-mode 1` (Prodigal 3 `-e 1`): the same gaps, but genes cannot run
    /// into them: every gene next to a gap has its own start and stop codon.
    Closed,
    /// `--gap-mode 2` (Prodigal 3 `-e 2`; what Prodigal 2.x does without `-m`):
    /// N's are not gaps. Genes may span runs of N (translated as `X`);
    /// ambiguous bases still never form start or stop codons.
    NoGaps,
    /// `-m` (Prodigal 2.x `-m`, pyrodigal `mask=True`): runs of at least
    /// [`seq::MASK_MIN_RUN`] ambiguous bases are masked and no gene is built
    /// across or into them, in training and prediction; shorter runs are
    /// treated as with [`GapMode::NoGaps`].
    Mask,
}

impl GapMode {
    /// Parse a Prodigal `-e` value: `0`, `1` or `2`.
    pub fn from_prodigal_e(value: &str) -> Option<GapMode> {
        match value.trim() {
            "0" => Some(GapMode::Partial),
            "1" => Some(GapMode::Closed),
            "2" => Some(GapMode::NoGaps),
            _ => None,
        }
    }

    /// Short description for logs.
    pub fn label(&self) -> &'static str {
        match self {
            GapMode::Partial => "0 (partial genes run into gaps)",
            GapMode::Closed => "1 (genes cannot run into gaps)",
            GapMode::NoGaps => "2 (N's are not gaps)",
            GapMode::Mask => "mask (-m: no gene crosses a run of >= 50 N)",
        }
    }
}

/// Options for [`run_records_with`] beyond [`RunConfig`] (kept apart so that
/// `RunConfig` stays source-compatible). Start from `RunOptions::default()`.
#[derive(Debug, Clone, Default)]
#[non_exhaustive]
pub struct RunOptions {
    /// Gap handling for the Prodigal engines (see [`GapMode`]).
    pub gap_mode: GapMode,
}

impl RunOptions {
    /// Set the gap handling for the Prodigal engines.
    pub fn with_gap_mode(mut self, gap_mode: GapMode) -> Self {
        self.gap_mode = gap_mode;
        self
    }
}

/// Whether an engine honours [`GapMode`] (the Prodigal engines do).
pub fn engine_supports_gap_mode(engine: EngineKind) -> bool {
    matches!(engine, EngineKind::Prodigal | EngineKind::ProdigalGv | EngineKind::ProdigalMeta)
}

/// The result of a run: every gene called, plus the `(id, length)` of every
/// input contig (needed for GFF/GenBank headers).
pub struct RunOutput {
    pub genes: Vec<Gene>,
    pub contigs: Vec<(String, usize)>,
    /// Which engine ran.
    pub engine: Engine,
    /// The genetic code that was used (resolved, after any auto-detection).
    pub trans_table: u32,
}

fn resolve_fixed(code: &CodeChoice, mode_default: Option<u32>) -> Option<u32> {
    match code {
        CodeChoice::Fixed(id) => Some(*id),
        CodeChoice::Default => mode_default,
        CodeChoice::Auto => None,
    }
}

fn kind_to_engine(k: EngineKind) -> Engine {
    match k {
        EngineKind::Phanotate => Engine::Phanotate,
        EngineKind::Prodigal => Engine::Prodigal,
        EngineKind::ProdigalGv => Engine::ProdigalGv,
        EngineKind::ProdigalMeta => Engine::ProdigalMeta,
        EngineKind::GeneMark => Engine::GeneMark,
        EngineKind::Augustus => Engine::Augustus,
        EngineKind::FragGeneScan => Engine::FragGeneScan,
    }
}

/// The genetic code an engine uses when the code is neither fixed nor
/// auto-detected.
fn engine_default_code(engine: EngineKind) -> u32 {
    match engine {
        EngineKind::Augustus => 1,
        _ => 11,
    }
}

/// Candidate codes for auto-detection: the mode's own list when the mode's
/// preset engine is the one running, otherwise the engine's built-in list.
/// An explicit `-g auto` with a single-table preset (bacteria, archaea,
/// phage, reads …) searches the standard alternative-code set instead of
/// silently returning that one table.
fn candidate_codes(engine: EngineKind, mode_cfg: &crate::model::ModeConfig, explicit_auto: bool) -> Vec<u32> {
    let list = if engine == mode_cfg.engine && !mode_cfg.auto_candidates.is_empty() {
        mode_cfg.auto_candidates.clone()
    } else {
        match engine {
            EngineKind::ProdigalGv => prodigalgv::GV_CANDIDATES.to_vec(),
            EngineKind::ProdigalMeta => prodigal_meta::META_CANDIDATES.to_vec(),
            _ => mode_cfg.auto_candidates.clone(),
        }
    };
    if explicit_auto && list.len() < 2 {
        let mut wide = list;
        for &t in crate::model::DEFAULT_CANDIDATES {
            if !wide.contains(&t) {
                wide.push(t);
            }
        }
        return wide;
    }
    list
}

/// Run a prediction over the given FASTA records.
///
/// Contigs are split at assembly gaps (runs of at least
/// [`seq::GAP_MIN_RUN`] N / ambiguous bases) and genes are called on each
/// gap-free segment as on a contig whose ends are open (unless `closed`): no
/// gene can span a gap, and genes running into one are reported partial at
/// that end. Self-training (and genetic-code detection) still uses the contigs
/// exactly as given, as Prodigal does, so genes away from a gap are the ones
/// Prodigal calls on the same sequence.
/// Coordinates are mapped back to the original contigs, genes are ordered by
/// input record then position, and ids are renumbered per record (records
/// that share an id are kept apart, as in 1.0.0).
///
/// Genetic code: a fixed table (`-g N` or the mode's default) is used as is.
/// Otherwise the Prodigal-family engines detect it among the candidate tables
/// (see [`prodigal::autodetect_code`]); with an explicit `-g auto`, GeneMark,
/// AUGUSTUS and PHANOTATE use the table detected the same way. FragGeneScan's
/// HMM only supports table 11.
///
/// This is [`run_records_with`] with the default [`RunOptions`]
/// ([`GapMode::Partial`]).
pub fn run_records(records: &[fasta::Record], cfg: &RunConfig) -> Result<RunOutput> {
    run_records_with(records, cfg, &RunOptions::default())
}

/// [`run_records`] with extra [`RunOptions`] — Prodigal's gap handling
/// ([`GapMode`]; honoured by the Prodigal engines, see
/// [`engine_supports_gap_mode`]). [`GapMode::Partial`] and [`GapMode::Closed`]
/// call the gap-free segments (with the ends at a gap open or closed);
/// [`GapMode::NoGaps`] and [`GapMode::Mask`] call whole contigs, exactly as
/// Prodigal 2.x does without / with `-m`. The genetic code is chosen the same
/// way in every gap mode, on the gap-free sequence.
pub fn run_records_with(records: &[fasta::Record], cfg: &RunConfig, opts: &RunOptions) -> Result<RunOutput> {
    let mode_cfg = config_for(cfg.mode);
    let engine = cfg.engine.unwrap_or(mode_cfg.engine);
    let closed = cfg.closed || mode_cfg.closed;
    let min_len = if cfg.min_len <= 0 {
        phanotate::MIN_ORF_LEN
    } else {
        cfg.min_len
    };
    let fixed = resolve_fixed(&cfg.code, mode_cfg.code);
    // Validate an explicitly requested table up front: an unsupported number
    // used to panic deep inside the Prodigal engines (`-g 7`).
    if let Some(id) = fixed {
        if GeneticCode::new(id).is_none() {
            return Err(anyhow!(
                "unsupported genetic code {id} (supported NCBI tables: {})",
                gencode::SUPPORTED.iter().map(|t| t.to_string()).collect::<Vec<_>>().join(", ")
            ));
        }
    }
    let explicit_auto = matches!(cfg.code, CodeChoice::Auto);
    let candidates = candidate_codes(engine, &mode_cfg, explicit_auto);
    // Prodigal's gap options apply to the Prodigal engines; the others always
    // split contigs at gaps.
    let gap_mode = if engine_supports_gap_mode(engine) { opts.gap_mode } else { GapMode::Partial };
    // Gap mode 2 and -m call whole contigs (Prodigal 2.x); modes 0 and 1 call
    // the gap-free segments.
    let whole_contigs = matches!(gap_mode, GapMode::NoGaps | GapMode::Mask);

    let contig_list: Vec<(String, usize)> =
        records.iter().map(|r| (r.id.clone(), r.seq.len())).collect();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.threads)
        .build()
        .map_err(|e| anyhow!("failed to build thread pool: {e}"))?;

    // The training input Prodigal itself would build: every record as given.
    let training_seq = || {
        let raw: Vec<&[u8]> = records.iter().map(|r| r.seq.as_slice()).collect();
        prodigal::build_training_seq_ascii(&raw)
    };
    // Gap-free segments: (record index, offset, length).
    let gap_free: Vec<(usize, usize, usize)> = records
        .iter()
        .enumerate()
        .flat_map(|(ri, r)| finalize::segments(&r.seq).into_iter().map(move |(o, l)| (ri, o, l)))
        .collect();
    // What the engine is run on: the gap-free segments, or (Prodigal gap mode
    // 2 / -m) each whole record as one segment.
    let segs: Vec<(usize, usize, usize)> = if whole_contigs {
        records.iter().enumerate().map(|(ri, r)| (ri, 0, r.seq.len())).collect()
    } else {
        gap_free.clone()
    };
    let seg_name = |k: usize| format!("dynamite-segment-{k}");
    let contigs_of = |parts: &[(usize, usize, usize)]| -> Vec<(String, Seq)> {
        parts
            .iter()
            .enumerate()
            .map(|(k, &(ri, off, len))| (seg_name(k), Seq::from_ascii(&records[ri].seq[off..off + len])))
            .collect()
    };
    let seg_contigs = || -> Vec<(String, Seq)> { contigs_of(&segs) };
    // Genetic code for engines that cannot search for it themselves: fixed, or
    // (explicit `-g auto`) detected with the Prodigal trainer.
    let detected_or = |default: u32, contigs: &[(String, Seq)]| -> u32 {
        match fixed {
            Some(id) => id,
            None if explicit_auto && !contigs.is_empty() => {
                let train = training_seq();
                let eval: Vec<&Seq> = contigs.iter().map(|(_, s)| s).collect();
                let tt = pool.install(|| prodigal::autodetect_code_on(&train, &eval, &candidates, closed, true).0);
                // The search only tells stop-codon sets apart (1 vs 11 look the
                // same), so keep the engine's own table unless the stops changed.
                match (GeneticCode::new(tt), GeneticCode::new(default)) {
                    (Some(a), Some(b)) if a.stop_mask() == b.stop_mask() => default,
                    _ => tt,
                }
            }
            None => default,
        }
    };

    // PHANOTATE works on whole contigs (it keeps its genome-wide statistics and
    // only restricts ORFs to the gap-free segments itself).
    if engine == EngineKind::Phanotate {
        let code_id = if fixed.is_none() && explicit_auto {
            detected_or(11, &seg_contigs())
        } else {
            fixed.unwrap_or(11)
        };
        let code = GeneticCode::new(code_id)
            .ok_or_else(|| anyhow!("unsupported genetic code {code_id}"))?;
        let tagged: Vec<(usize, Gene)> = pool.install(|| {
            records
                .par_iter()
                .enumerate()
                .flat_map_iter(|(ri, r)| {
                    let (g, _score) = phanotate::predict(&r.id, &r.seq, &code, min_len, closed);
                    g.into_iter().map(move |g| (ri, g))
                })
                .collect()
        });
        let genes = finalize::order_by_record(tagged);
        return Ok(RunOutput {
            genes,
            contigs: contig_list,
            engine: Engine::Phanotate,
            trans_table: code_id,
        });
    }

    let contigs = seg_contigs();
    let (genes, table) = match engine {
        EngineKind::Prodigal | EngineKind::ProdigalGv | EngineKind::ProdigalMeta => {
            let (fixed_for, is_meta, tag) = match engine {
                // The plain Prodigal preset keeps table 11 unless asked to search.
                EngineKind::Prodigal => (
                    if fixed.is_none() && !explicit_auto { Some(11) } else { fixed },
                    false,
                    Engine::Prodigal,
                ),
                EngineKind::ProdigalGv => (fixed, true, Engine::ProdigalGv),
                _ => (fixed, true, Engine::ProdigalMeta),
            };
            let train = training_seq();
            if !whole_contigs {
                // Gap mode 0: segment ends at a gap are open, like contig ends
                // (`--closed` closes all). Gap mode 1 closes them.
                let ends: Vec<(bool, bool)> = segs
                    .iter()
                    .map(|&(ri, off, len)| {
                        let gap_left = off > 0;
                        let gap_right = off + len < records[ri].seq.len();
                        let shut = gap_mode == GapMode::Closed;
                        (closed || (shut && gap_left), closed || (shut && gap_right))
                    })
                    .collect();
                pool.install(|| {
                    prodigal::dp_call_trained_on_ends(&train, &contigs, &ends, fixed_for, &candidates, is_meta, closed, tag)
                })
            } else {
                // The genetic code is chosen exactly as in the default gap
                // mode (statistics on the gap-free sequence); the whole
                // contigs are then called as Prodigal 2.x does, with `-m`'s
                // masks if asked.
                let table = match fixed_for {
                    Some(id) => id,
                    None => {
                        let parts = contigs_of(&gap_free);
                        let eval: Vec<&Seq> = parts.iter().map(|(_, s)| s).collect();
                        pool.install(|| prodigal::autodetect_code_on(&train, &eval, &candidates, closed, is_meta).0)
                    }
                };
                let mask = gap_mode == GapMode::Mask;
                let masks: Vec<Vec<(usize, usize)>> = contigs
                    .iter()
                    .map(|(_, s)| if mask { seq::ambiguous_runs(s, seq::MASK_MIN_RUN) } else { Vec::new() })
                    .collect();
                let train_masks = if mask { seq::ambiguous_runs(&train, seq::MASK_MIN_RUN) } else { Vec::new() };
                let genes = pool.install(|| {
                    prodigal::dp_call_masked(&train, &train_masks, &contigs, &masks, table, is_meta, closed, tag)
                });
                (genes, table)
            }
        }
        EngineKind::GeneMark => {
            let tt = detected_or(engine_default_code(engine), &contigs);
            pool.install(|| genemark::call(&contigs, Some(tt), closed))
        }
        EngineKind::Augustus => {
            let tt = detected_or(engine_default_code(engine), &contigs);
            pool.install(|| augustus::call(&contigs, Some(tt), closed))
        }
        EngineKind::FragGeneScan => {
            let whole = !matches!(cfg.mode, crate::model::Mode::Reads | crate::model::Mode::LongReads);
            pool.install(|| fraggenescan::call_closed(&contigs, whole, closed))
        }
        EngineKind::Phanotate => unreachable!(),
    };

    // Map segment-local genes back onto their records (by index, so records
    // that share an id stay apart).
    let index: std::collections::HashMap<String, usize> =
        (0..segs.len()).map(|k| (seg_name(k), k)).collect();
    let tagged: Vec<(usize, Gene)> = genes
        .into_iter()
        .map(|mut g| match index.get(&g.contig) {
            Some(&k) => {
                let (ri, off, _) = segs[k];
                g.begin += off;
                g.end += off;
                g.contig = records[ri].id.clone();
                (ri, g)
            }
            None => (usize::MAX, g),
        })
        .collect();
    let genes = finalize::order_by_record(tagged);

    Ok(RunOutput {
        genes,
        contigs: contig_list,
        engine: kind_to_engine(engine),
        trans_table: table,
    })
}
