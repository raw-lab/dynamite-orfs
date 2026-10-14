//! Organism *modes* and how they map onto an engine, a genetic code, and the
//! meta / self-training switch.
//!
//! A mode is a convenience preset that resolves to one of DYNAMITE's engines
//! (see [`EngineKind`]). The engine can also be selected directly on the command
//! line, overriding the mode's default.

/// The concrete gene-calling engines DYNAMITE ships.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EngineKind {
    /// PHANOTATE graph / shortest-path caller (phages).
    Phanotate,
    /// Prodigal preset — bacteria & archaea (single self-trained genome).
    Prodigal,
    /// Prodigal-GV preset — phages, giant viruses, crassphages.
    ProdigalGv,
    /// Prodigal metagenomic preset — fragmented / mixed input.
    ProdigalMeta,
    /// GeneMark-style self-training Markov caller.
    GeneMark,
    /// AUGUSTUS-style eukaryotic caller (single-exon).
    Augustus,
    /// FragGeneScan read caller (vendored FragGeneScanRs HMM, `complete` model).
    FragGeneScan,
}

impl EngineKind {
    /// Parse an engine name (the module names, plus a few aliases).
    pub fn parse(s: &str) -> Option<EngineKind> {
        match s.to_ascii_lowercase().replace(['_', ' '], "-").as_str() {
            "phanotate" => Some(EngineKind::Phanotate),
            "prodigal" => Some(EngineKind::Prodigal),
            "prodigalgv" | "prodigal-gv" | "gv" => Some(EngineKind::ProdigalGv),
            "prodigal-meta" | "prodigalmeta" | "meta" => Some(EngineKind::ProdigalMeta),
            "genemark" => Some(EngineKind::GeneMark),
            "augustus" => Some(EngineKind::Augustus),
            "fraggenescan" | "fgs" => Some(EngineKind::FragGeneScan),
            _ => None,
        }
    }

    /// Short label.
    pub fn label(self) -> &'static str {
        match self {
            EngineKind::Phanotate => "phanotate",
            EngineKind::Prodigal => "prodigal",
            EngineKind::ProdigalGv => "prodigalgv",
            EngineKind::ProdigalMeta => "prodigal-meta",
            EngineKind::GeneMark => "genemark",
            EngineKind::Augustus => "augustus",
            EngineKind::FragGeneScan => "fraggenescan",
        }
    }
}

/// Organism preset selected on the command line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Phage,
    Virus,
    GiantVirus,
    Crassphage,
    Bacteria,
    Archaea,
    Meta,
    Eukaryote,
    Reads,
    LongReads,
    Auto,
}

impl Mode {
    /// Parse a mode name (case-insensitive, a few aliases allowed).
    pub fn parse(s: &str) -> Option<Mode> {
        match s.to_ascii_lowercase().replace(['_', ' '], "-").as_str() {
            "phage" | "phages" => Some(Mode::Phage),
            "virus" | "viral" | "viruses" => Some(Mode::Virus),
            "giant-virus" | "giant" | "ncldv" | "nucleocytoviricota" => Some(Mode::GiantVirus),
            "crassphage" | "crass" | "crassvirales" => Some(Mode::Crassphage),
            "bacteria" | "bacterium" | "bact" => Some(Mode::Bacteria),
            "archaea" | "archaeon" | "arch" => Some(Mode::Archaea),
            "meta" | "metagenome" | "metagenomic" => Some(Mode::Meta),
            "eukaryote" | "eukaryotic" | "euk" => Some(Mode::Eukaryote),
            "reads" | "read" | "short-reads" => Some(Mode::Reads),
            "long-reads" | "long" | "longreads" | "nanopore" | "pacbio" | "ont" | "hifi" => Some(Mode::LongReads),
            "auto" | "default" => Some(Mode::Auto),
            _ => None,
        }
    }

    /// Human-readable name.
    pub fn label(self) -> &'static str {
        match self {
            Mode::Phage => "phage",
            Mode::Virus => "virus",
            Mode::GiantVirus => "giant-virus",
            Mode::Crassphage => "crassphage",
            Mode::Bacteria => "bacteria",
            Mode::Archaea => "archaea",
            Mode::Meta => "meta",
            Mode::Eukaryote => "eukaryote",
            Mode::Reads => "reads",
            Mode::LongReads => "long-reads",
            Mode::Auto => "auto",
        }
    }
}

/// Whether the genetic code is fixed, auto-detected, or left to the mode.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CodeChoice {
    /// Use whatever the mode prefers.
    Default,
    /// Auto-detect among the mode's candidate codes.
    Auto,
    /// Force a specific NCBI translation table.
    Fixed(u32),
}

/// Resolved configuration for a mode.
#[derive(Clone, Debug)]
pub struct ModeConfig {
    /// Which engine to run.
    pub engine: EngineKind,
    /// Fixed genetic code, or `None` to auto-detect from `auto_candidates`.
    pub code: Option<u32>,
    /// Candidate codes to try when auto-detecting.
    pub auto_candidates: Vec<u32>,
    /// Whether to treat input as metagenomic / fragmented.
    pub is_meta: bool,
    /// Default for `closed` (don't allow genes to run off contig ends).
    pub closed: bool,
}

/// The standard, broad auto-detection candidate set.
pub const DEFAULT_CANDIDATES: &[u32] = &[11, 4, 15, 25];
/// A wider set used for giant viruses (adds code 6).
pub const GIANT_CANDIDATES: &[u32] = &[11, 1, 4, 15, 25, 6];

/// Build the preset configuration for a mode.
pub fn config_for(mode: Mode) -> ModeConfig {
    use EngineKind::*;
    match mode {
        Mode::Phage => ModeConfig {
            engine: Phanotate,
            code: Some(11),
            auto_candidates: vec![11],
            is_meta: true,
            closed: false,
        },
        Mode::Virus => ModeConfig {
            engine: ProdigalGv,
            code: None,
            auto_candidates: DEFAULT_CANDIDATES.to_vec(),
            is_meta: true,
            closed: false,
        },
        Mode::GiantVirus => ModeConfig {
            engine: ProdigalGv,
            code: None,
            auto_candidates: GIANT_CANDIDATES.to_vec(),
            is_meta: true,
            closed: false,
        },
        Mode::Crassphage => ModeConfig {
            engine: ProdigalGv,
            code: None,
            auto_candidates: vec![11, 15, 4, 25],
            is_meta: true,
            closed: false,
        },
        Mode::Bacteria => ModeConfig {
            engine: Prodigal,
            code: Some(11),
            auto_candidates: vec![11],
            is_meta: false,
            closed: false,
        },
        Mode::Archaea => ModeConfig {
            engine: Prodigal,
            code: Some(11),
            auto_candidates: vec![11],
            is_meta: false,
            closed: false,
        },
        Mode::Meta => ModeConfig {
            engine: ProdigalMeta,
            code: None,
            auto_candidates: DEFAULT_CANDIDATES.to_vec(),
            is_meta: true,
            closed: false,
        },
        Mode::Eukaryote => ModeConfig {
            engine: Augustus,
            code: Some(1),
            auto_candidates: vec![1],
            is_meta: false,
            closed: false,
        },
        Mode::Reads => ModeConfig {
            engine: FragGeneScan,
            code: Some(11),
            auto_candidates: vec![11],
            is_meta: true,
            closed: false,
        },
        Mode::LongReads => ModeConfig {
            engine: FragGeneScan,
            code: Some(11),
            auto_candidates: vec![11],
            is_meta: true,
            closed: false,
        },
        Mode::Auto => ModeConfig {
            engine: ProdigalMeta,
            code: None,
            auto_candidates: DEFAULT_CANDIDATES.to_vec(),
            is_meta: true,
            closed: false,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_modes() {
        assert_eq!(Mode::parse("Phage"), Some(Mode::Phage));
        assert_eq!(Mode::parse("giant_virus"), Some(Mode::GiantVirus));
        assert_eq!(Mode::parse("BACTERIA"), Some(Mode::Bacteria));
        assert_eq!(Mode::parse("eukaryote"), Some(Mode::Eukaryote));
        assert_eq!(Mode::parse("nonsense"), None);
    }

    #[test]
    fn parse_engines() {
        assert_eq!(EngineKind::parse("prodigalgv"), Some(EngineKind::ProdigalGv));
        assert_eq!(EngineKind::parse("prodigal-meta"), Some(EngineKind::ProdigalMeta));
        assert_eq!(EngineKind::parse("fraggenescan"), Some(EngineKind::FragGeneScan));
    }

    #[test]
    fn presets() {
        assert!(matches!(config_for(Mode::Phage).engine, EngineKind::Phanotate));
        assert!(matches!(config_for(Mode::Bacteria).engine, EngineKind::Prodigal));
        assert_eq!(config_for(Mode::Bacteria).code, Some(11));
        assert!(config_for(Mode::Meta).code.is_none());
        assert!(matches!(config_for(Mode::Eukaryote).engine, EngineKind::Augustus));
    }
}
