//! DYNAMITE command-line interface.
//!
//! ```text
//! dynamite call -i genome.fna -p bacteria -f gff3 -o genes.gff
//! dynamite call -i phage.fna  -p phage    -f genbank
//! dynamite call -i ncldv.fna  -e prodigalgv -g auto -f faa
//! dynamite translate -i seqs.fna --orfs
//! dynamite codon -i genome.fna -p meta -o codons.tsv
//! dynamite rna   -i genome.fna -o ncrna.gff
//! dynamite doctor
//! ```

use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::{anyhow, Context, Result};
use clap::{Args, Parser, Subcommand};

use dynamite::fasta::{self, Record};
use dynamite::gencode::GeneticCode;
use dynamite::gene::Format;
use dynamite::model::{CodeChoice, EngineKind, Mode};
use dynamite::annotate::{annotate_with, AnnotateConfig};
use dynamite::{codon, doctor, rna, run_records, run_records_with, sixframe, GapMode, RunConfig, RunOptions};

const BANNER: &str = concat!(
    r#"                                                      _.-^^~~~~~^^-._
                                                  _.-~ ░░▒▒▓▓▓▓▒▒░░ ~-._
                                                .'  ░▒▓████████████▓▒░  '.
                                               ( ░▒▓██████████████████▓▒░ )
                                                '._ ░▒▓████████████▓▒░ _.'
                                                   '~-.░▒▓██████▓▒░.-~'
 ______   ___   _    _    __  __ ___ _____ _____      \ ░▒▓███▓▒░ /
|  _ \ \ / / \ | |  / \  |  \/  |_ _|_   _| ____|      \ ▒▓███▓▒ /
| | | \ V /|  \| | / _ \ | |\/| || |  | | |  _|   ╭──╮  \░▒▓█▓▒░/
| |_| || | | |\  |/ ___ \| |  | || |  | | | |___  ││▒▓██╼╾╼ ✸
|____/ |_| |_| \_/_/   \_\_|  |_|___| |_| |_____| ╰──╯"#,
    "\nv", env!("CARGO_PKG_VERSION"),
    "  \u{00b7}  ORF / gene caller for all domains of life\n",
    "bacteria \u{00b7} archaea \u{00b7} phage \u{00b7} giant virus \u{00b7} crassphage \u{00b7} metagenome \u{00b7} eukaryote \u{00b7} reads"
);

/// Shown after the option list for `dynamite call --help`.
const ENGINES_HELP: &str = "\
ENGINE MODULES  (pure-Rust re-writes — DYNAMITE never shells out to an external binary)

  engine name    DYNAMITE module       re-writes / faithful to
  -----------    ---------------       -----------------------
  prodigal       prodigal.rs           Prodigal (Hyatt 2010) single-genome dynamic-
                                       programming caller — bacteria & archaea
  prodigalgv     prodigalgv.rs         Prodigal-GV (Camargo) — viruses & giant viruses;
                                       self-trains + alt-code search [11,1,4,15,25,6]
  prodigal-meta  prodigal_meta.rs      Prodigal anonymous / metagenomic mode
  phanotate      phanotate.rs          PHANOTATE (McNair 2019) — phage gene caller
  genemark       genemark.rs           GeneMark-style 3-periodic Markov caller
  augustus       augustus.rs           AUGUSTUS-style single-exon eukaryote caller
  fraggenescan   fgs/  (vendored)      FragGeneScanRs (Van der Jeugt 2022) — the REAL
                                       FragGeneScan HMM (Viterbi) with its `complete`
                                       model (no sequencing-error frameshift
                                       correction); short & long reads
  (six-frame)    sixframe.rs           6-frame translation / ORF finder — see
                                       `dynamite translate [--orfs]`

MODE -> ENGINE
  bacteria,archaea            -> prodigal        phage                    -> phanotate
  virus,giant-virus,crassphage-> prodigalgv      meta,auto                -> prodigal-meta
  eukaryote                   -> augustus        reads (short)            -> fraggenescan
  long-reads (nanopore|pacbio|ont|hifi)          -> fraggenescan (same decode as reads)

GAPS (runs of N)  prodigal, prodigalgv, prodigal-meta — Prodigal's -e / -m
  --gap-mode 0   default: runs of >= 10 N are gaps; partial genes run into them
                 as into a contig end (as if the contig were cut there)
  --gap-mode 1   genes cannot run into gaps (each gene next to a gap is complete)
  --gap-mode 2   N's are not gaps: genes may span them (Prodigal 2.x default)
  -m, --mask     Prodigal 2.x -m: no gene is built across or into >= 50 N
  Prodigal's own spelling works too: -e 0|1|2 (e.g. -e prodigal -e 1).
  The other engines always use gap mode 0.

INPUT   nucleotide FASTA — contigs, short reads, or long reads (plain or .gz)
OUTPUT  -f  gff3 | gtf | genbank | faa | fna (alias: ffn) | tsv";

/// Shown after the option list for `dynamite all --help`.
const ALL_HELP: &str = "\
WHAT `all` WRITES  (one file per format — raw gene calls, no annotation layer)

  {prefix}.gff    gene calls in GFF3              {prefix}.faa    CDS proteins (amino acids)
  {prefix}.gtf    gene calls in GTF               {prefix}.ffn    CDS gene nucleotides
  {prefix}.gbk    GenBank: features + ORIGIN seq  {prefix}.fna    input contig sequences
  {prefix}.tsv    coordinate / feature table

EXAMPLE
  dynamite all -i genome.fna -p bacteria -o results/ --prefix mygenome
    -> results/mygenome.gff  .gtf  .gbk  .faa  .ffn  .fna  .tsv

ENGINE / MODE selection is shared with `dynamite call` (run `dynamite call --help` for the full
engine-module map). For locus_tags, products, tRNA/rRNA and NCBI .tbl/.fsa/.sqn inputs, use
`dynamite annotate` instead.";

/// DYNAMITE — an ORF / gene caller for all domains of life.
#[derive(Parser, Debug)]
#[command(
    name = "dynamite",
    version,
    about = "DYNAMITE: ORF / gene caller for phages, viruses, giant viruses, bacteria, archaea, eukaryotes and reads.",
    before_help = BANNER,
    before_long_help = BANNER,
    arg_required_else_help = true,
    propagate_version = true
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Call genes / ORFs (default workflow).
    Call(CallArgs),
    /// Run the caller once and write every output format (no annotation layer).
    All(AllArgs),
    /// Whole-genome annotation: write the full Prokka-style output bundle.
    Annotate(AnnotateArgs),
    /// Six-frame translation (and optional six-frame ORF listing).
    Translate(TranslateArgs),
    /// Codon usage, bias (RSCU / ENC) and codon density.
    Codon(CodonArgs),
    /// Heuristic tRNA / rRNA (SSU + LSU) detection.
    Rna(RnaArgs),
    /// Environment report and end-to-end self-test.
    Doctor,
}

#[derive(Args, Debug)]
#[command(after_help = ENGINES_HELP, after_long_help = ENGINES_HELP)]
struct CallArgs {
    /// Input nucleotide FASTA (optionally gzip-compressed).
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    input: Option<PathBuf>,
    /// Positional input FASTA (alternative to -i).
    #[arg(value_name = "FASTA")]
    input_pos: Option<PathBuf>,
    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,
    /// Output format: gff3 | gtf | genbank | faa | fna | ffn | tsv.
    #[arg(short = 'f', long = "format", default_value = "gff3")]
    format: String,
    /// Organism mode: phage | virus | giant-virus | crassphage | bacteria |
    /// archaea | meta | eukaryote | reads | long-reads | auto.
    #[arg(short = 'p', long = "mode", default_value = "auto")]
    mode: String,
    /// Engine override (a module from the list below): phanotate | prodigal |
    /// prodigalgv | prodigal-meta | genemark | augustus | fraggenescan. A
    /// number instead sets the gap mode, as Prodigal's own -e does (-e 0|1|2,
    /// see --gap-mode); give -e twice for both, e.g. -e prodigal -e 1.
    #[arg(short = 'e', long = "engine", value_name = "ENGINE|0|1|2")]
    engine: Vec<String>,
    /// Genetic code: an NCBI table number (e.g. 11, 4, 15, 25) or "auto".
    #[arg(short = 'g', long = "code", value_name = "CODE|auto")]
    code: Option<String>,
    /// Gap handling for runs of >= 10 N (prodigal, prodigalgv, prodigal-meta):
    /// 0 = partial genes run into gaps (default); 1 = genes cannot run into
    /// gaps; 2 = N's are not gaps (Prodigal 2.x without -m).
    #[arg(long = "gap-mode", value_name = "0|1|2")]
    gap_mode: Option<String>,
    /// Prodigal 2.x -m: treat runs of >= 50 N as masked sequence; no gene is
    /// built across or into them (prodigal, prodigalgv, prodigal-meta).
    #[arg(short = 'm', long = "mask", default_value_t = false)]
    mask: bool,
    /// Minimum ORF length in nucleotides (PHANOTATE engine; default 90).
    #[arg(short = 'l', long = "minlen", default_value_t = 0)]
    minlen: i64,
    /// Do not allow genes to run off the ends of a contig.
    #[arg(long = "closed", default_value_t = false)]
    closed: bool,
    /// Accepted for compatibility only: genes are never wrapped across the
    /// origin. Contig ends are open (partial genes allowed) unless --closed.
    #[arg(long = "circular", default_value_t = false)]
    circular: bool,
    /// Worker threads (0 = all available cores).
    #[arg(short = 't', long = "threads", default_value_t = 0)]
    threads: usize,
}

#[derive(Args, Debug)]
#[command(
    about = "Run the caller once and write every output format (no annotation layer).",
    long_about = "Run the chosen engine once and write the gene calls to EVERY output format \
                  in one pass: .gff .gtf .gbk .faa .ffn .fna .tsv. This is plain gene-calling \
                  output — no locus_tags, products, tRNA/rRNA, or NCBI submission files (use \
                  `dynamite annotate` for those).",
    after_help = ALL_HELP,
    after_long_help = ALL_HELP
)]
struct AllArgs {
    /// Input nucleotide FASTA: contigs, short reads, or long reads (.gz ok). Use this or the positional arg.
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    input: Option<PathBuf>,
    /// Input FASTA given positionally (alternative to -i/--input).
    #[arg(value_name = "FASTA")]
    input_pos: Option<PathBuf>,
    /// Directory to write the output files into (created if missing).
    #[arg(short = 'o', long = "outdir", value_name = "DIR", default_value = "dynamite_out")]
    outdir: PathBuf,
    /// Base name shared by every output file: {prefix}.gff, {prefix}.faa, ... .
    #[arg(long = "prefix", value_name = "NAME", default_value = "dynamite")]
    prefix: String,
    /// Organism mode: phage | virus | giant-virus | crassphage | bacteria | archaea | meta | eukaryote | reads | long-reads | auto.
    #[arg(short = 'p', long = "mode", value_name = "MODE", default_value = "auto")]
    mode: String,
    /// Force a specific engine: phanotate | prodigal | prodigalgv | prodigal-meta | genemark | augustus | fraggenescan (default: from --mode). A number instead sets the gap mode, as Prodigal's own -e does (-e 0|1|2; e.g. -e prodigal -e 1).
    #[arg(short = 'e', long = "engine", value_name = "ENGINE|0|1|2")]
    engine: Vec<String>,
    /// Genetic code: an NCBI translation-table number (e.g. 11) or "auto" to detect.
    #[arg(short = 'g', long = "code", value_name = "CODE|auto")]
    code: Option<String>,
    /// Gap handling for runs of >= 10 N (prodigal engines): 0 = partial genes run into gaps (default); 1 = genes cannot run into gaps; 2 = N's are not gaps.
    #[arg(long = "gap-mode", value_name = "0|1|2")]
    gap_mode: Option<String>,
    /// Prodigal 2.x -m: no gene is built across or into a run of >= 50 N (prodigal engines).
    #[arg(short = 'm', long = "mask", default_value_t = false)]
    mask: bool,
    /// Require complete genes only (do not allow genes to run off contig ends).
    #[arg(long = "closed", default_value_t = false)]
    closed: bool,
    /// Minimum ORF length in nucleotides (PHANOTATE engine only).
    #[arg(short = 'l', long = "minlen", value_name = "NT", default_value_t = 0)]
    minlen: i64,
    /// Worker threads; 0 = use all available CPU cores.
    #[arg(short = 't', long = "threads", value_name = "N", default_value_t = 0)]
    threads: usize,
}

#[derive(Args, Debug)]
#[command(after_help = ENGINES_HELP, after_long_help = ENGINES_HELP)]
struct AnnotateArgs {
    /// Input nucleotide FASTA (contigs/reads; optionally gzip-compressed).
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    input: Option<PathBuf>,
    /// Positional input FASTA (alternative to -i).
    #[arg(value_name = "FASTA")]
    input_pos: Option<PathBuf>,
    /// Output directory for the annotation bundle.
    #[arg(short = 'o', long = "outdir", value_name = "DIR", default_value = "dynamite_out")]
    outdir: PathBuf,
    /// Base filename for outputs ({prefix}.gff, {prefix}.faa, ...).
    #[arg(long = "prefix", default_value = "dynamite")]
    prefix: String,
    /// locus_tag prefix (default: uppercased --prefix).
    #[arg(long = "locus", value_name = "TAG")]
    locus: Option<String>,
    /// Organism name (GenBank / Sequin metadata).
    #[arg(long = "organism", default_value = "Unknown organism")]
    organism: String,
    /// Strain (GenBank / Sequin metadata).
    #[arg(long = "strain", default_value = "")]
    strain: String,
    /// Organism mode (see `dynamite call --help`).
    #[arg(short = 'p', long = "mode", default_value = "auto")]
    mode: String,
    /// Engine override (a module from the list below). A number instead sets
    /// the gap mode, as Prodigal's own -e does (-e 0|1|2; e.g. -e prodigal -e 1).
    #[arg(short = 'e', long = "engine", value_name = "ENGINE|0|1|2")]
    engine: Vec<String>,
    /// Genetic code: an NCBI table number or "auto".
    #[arg(short = 'g', long = "code", value_name = "CODE|auto")]
    code: Option<String>,
    /// Gap handling for runs of >= 10 N (prodigal engines): 0 = partial genes
    /// run into gaps (default); 1 = genes cannot run into gaps; 2 = N's are not gaps.
    #[arg(long = "gap-mode", value_name = "0|1|2")]
    gap_mode: Option<String>,
    /// Prodigal 2.x -m: no gene is built across or into a run of >= 50 N
    /// (prodigal engines).
    #[arg(short = 'm', long = "mask", default_value_t = false)]
    mask: bool,
    /// Skip tRNA / rRNA detection.
    #[arg(long = "no-rna", default_value_t = false)]
    no_rna: bool,
    /// Do not allow genes to run off contig ends.
    #[arg(long = "closed", default_value_t = false)]
    closed: bool,
    /// Worker threads (0 = all available cores).
    #[arg(short = 't', long = "threads", default_value_t = 0)]
    threads: usize,
}

#[derive(Args, Debug)]
struct TranslateArgs {
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    input: Option<PathBuf>,
    #[arg(value_name = "FASTA")]
    input_pos: Option<PathBuf>,
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,
    /// Genetic code (NCBI table number; default 1 = standard).
    #[arg(short = 'g', long = "code", default_value_t = 1)]
    code: u32,
    /// Instead of six-frame protein FASTA, list six-frame ORFs.
    #[arg(long = "orfs", default_value_t = false)]
    orfs: bool,
    /// Output format for --orfs: gff3 | gtf | genbank | faa | fna | tsv (default faa).
    #[arg(short = 'f', long = "format", default_value = "faa")]
    format: String,
    /// Minimum ORF length for --orfs (nt).
    #[arg(short = 'l', long = "minlen", default_value_t = 90)]
    minlen: usize,
    /// Require a start codon for --orfs (otherwise list all stop-to-stop ORFs,
    /// start-less ones flagged 5'-partial).
    #[arg(long = "require-start", default_value_t = false)]
    require_start: bool,
}

#[derive(Args, Debug)]
struct CodonArgs {
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    input: Option<PathBuf>,
    #[arg(value_name = "FASTA")]
    input_pos: Option<PathBuf>,
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,
    /// Genetic code (NCBI table number; default 11).
    #[arg(short = 'g', long = "code", default_value_t = 11)]
    code: u32,
    /// Treat the input FASTA as in-frame CDS (skip gene calling).
    #[arg(long = "cds", default_value_t = false)]
    cds: bool,
    /// When not --cds, the mode used to call genes first (default meta).
    #[arg(short = 'p', long = "mode", default_value = "meta")]
    mode: String,
}

#[derive(Args, Debug)]
struct RnaArgs {
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    input: Option<PathBuf>,
    #[arg(value_name = "FASTA")]
    input_pos: Option<PathBuf>,
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,
    /// Genetic code used to translate predicted anticodons (default 11).
    #[arg(short = 'g', long = "code", default_value_t = 11)]
    code: u32,
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => code,
        Err(e) => {
            eprintln!("dynamite: error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<ExitCode> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Call(a) => cmd_call(a).map(|_| ExitCode::SUCCESS),
        Cmd::All(a) => cmd_all(a).map(|_| ExitCode::SUCCESS),
        Cmd::Annotate(a) => cmd_annotate(a).map(|_| ExitCode::SUCCESS),
        Cmd::Translate(a) => cmd_translate(a).map(|_| ExitCode::SUCCESS),
        Cmd::Codon(a) => cmd_codon(a).map(|_| ExitCode::SUCCESS),
        Cmd::Rna(a) => cmd_rna(a).map(|_| ExitCode::SUCCESS),
        Cmd::Doctor => {
            let (report, ok) = doctor::run();
            println!("{report}");
            Ok(if ok { ExitCode::SUCCESS } else { ExitCode::FAILURE })
        }
    }
}

fn read_input(input: Option<PathBuf>, input_pos: Option<PathBuf>) -> Result<(PathBuf, Vec<Record>)> {
    let path = input
        .or(input_pos)
        .ok_or_else(|| anyhow!("no input FASTA given (use -i/--input or a positional path)"))?;
    let records = fasta::read_fasta(&path)
        .with_context(|| format!("failed to read FASTA from {}", path.display()))?;
    if records.is_empty() {
        return Err(anyhow!("no sequences found in {}", path.display()));
    }
    Ok((path, records))
}

fn write_out(output: &Option<PathBuf>, buf: &str) -> Result<()> {
    match output {
        Some(path) => fs::write(path, buf.as_bytes())
            .with_context(|| format!("failed to write output to {}", path.display())),
        None => {
            let stdout = io::stdout();
            let mut h = stdout.lock();
            h.write_all(buf.as_bytes()).context("failed to write to stdout")
        }
    }
}

/// Pre-flight warnings for a gene-calling run (stderr only).
fn warn_run(mode: Mode, engine: Option<EngineKind>, code: &CodeChoice, records: &[Record], gap_given: bool) {
    let eng = engine.unwrap_or_else(|| dynamite::model::config_for(mode).engine);
    if gap_given && !dynamite::engine_supports_gap_mode(eng) {
        eprintln!(
            "dynamite: warning: gap modes (-e 0|1|2, --gap-mode, -m) apply to the prodigal, prodigalgv and \
             prodigal-meta engines; {} always splits contigs at gaps (gap mode 0)",
            eng.label()
        );
    }
    if eng == EngineKind::FragGeneScan {
        match code {
            CodeChoice::Fixed(n) if *n != 11 => eprintln!(
                "dynamite: warning: the FragGeneScan HMM is trained for code 11; -g {n} is ignored by this engine"
            ),
            CodeChoice::Auto => eprintln!(
                "dynamite: warning: the FragGeneScan HMM is trained for code 11; -g auto is ignored by this engine"
            ),
            _ => {}
        }
    }
    // Sequence ids name the records in every output format; same-named records
    // are called separately but cannot be told apart downstream.
    {
        let mut seen = std::collections::HashSet::new();
        let mut dup_set = std::collections::HashSet::new();
        let mut dups: Vec<&str> = Vec::new();
        for r in records {
            if !seen.insert(r.id.as_str()) && dup_set.insert(r.id.as_str()) {
                dups.push(r.id.as_str());
            }
        }
        if !dups.is_empty() {
            let shown: Vec<&str> = dups.iter().take(5).copied().collect();
            eprintln!(
                "dynamite: warning: {} sequence id(s) occur more than once ({}{}); outputs identify sequences \
                 by id, so give each record a unique id",
                dups.len(),
                shown.join(", "),
                if dups.len() > 5 { ", …" } else { "" }
            );
        }
    }
    if matches!(eng, EngineKind::Prodigal | EngineKind::ProdigalGv | EngineKind::ProdigalMeta) {
        let bases: usize = records
            .iter()
            .map(|r| r.seq.iter().filter(|&&b| !dynamite::seq::is_ambiguous_byte(b)).count())
            .sum();
        if bases < 20_000 {
            eprintln!(
                "dynamite: warning: the {} engine self-trains on the input, but only {bases} bp were given \
                 (Prodigal needs >= 20,000 bp to train); gene calls and genetic-code detection may be unreliable",
                eng.label()
            );
        }
    }
}

fn parse_code(s: Option<&str>) -> Result<CodeChoice> {
    Ok(match s {
        None => CodeChoice::Default,
        Some(s) if s.eq_ignore_ascii_case("auto") => CodeChoice::Auto,
        Some(s) => {
            let n: u32 = s
                .parse()
                .with_context(|| format!("invalid genetic code '{s}' (expected a number or 'auto')"))?;
            CodeChoice::Fixed(n)
        }
    })
}

/// Resolve the `-e/--engine` values — an engine name, or Prodigal's `-e`
/// gap mode `0|1|2` — together with `--gap-mode` and `-m/--mask`. Returns the
/// engine override, the gap mode, and whether a gap option was given.
fn parse_engine_and_gaps(
    values: &[String],
    gap_mode: Option<&str>,
    mask: bool,
) -> Result<(Option<EngineKind>, GapMode, bool)> {
    let mut engine: Option<EngineKind> = None;
    let mut gap: Option<GapMode> = None;
    let mut first = String::new();
    let mut set_gap = |g: GapMode, what: String| -> Result<()> {
        match gap {
            Some(prev) if prev != g => Err(anyhow!("conflicting gap modes: {first} and {what}")),
            Some(_) => Ok(()),
            None => {
                gap = Some(g);
                first = what;
                Ok(())
            }
        }
    };
    for v in values {
        if let Some(g) = GapMode::from_prodigal_e(v) {
            set_gap(g, format!("-e {}", v.trim()))?;
        } else {
            let e = EngineKind::parse(v).ok_or_else(|| {
                anyhow!(
                    "unknown engine '{v}' (try: phanotate, prodigal, prodigalgv, prodigal-meta, genemark, augustus, \
                     fraggenescan — or a gap mode 0, 1, 2)"
                )
            })?;
            if engine.is_some_and(|prev| prev != e) {
                return Err(anyhow!("more than one engine given with -e/--engine"));
            }
            engine = Some(e);
        }
    }
    if let Some(v) = gap_mode {
        let g = GapMode::from_prodigal_e(v)
            .ok_or_else(|| anyhow!("invalid --gap-mode '{v}' (expected 0, 1 or 2)"))?;
        set_gap(g, format!("--gap-mode {}", v.trim()))?;
    }
    if mask {
        if gap.is_some() {
            return Err(anyhow!(
                "-m/--mask (Prodigal 2.x masking) and a gap mode (-e 0|1|2, --gap-mode) are alternatives; give one"
            ));
        }
        return Ok((engine, GapMode::Mask, true));
    }
    Ok((engine, gap.unwrap_or_default(), gap.is_some()))
}

/// ` gap-mode=…` for the run summary when it is not the default.
fn gap_note(engine: EngineKind, gap: GapMode) -> String {
    if gap != GapMode::Partial && dynamite::engine_supports_gap_mode(engine) {
        match gap {
            GapMode::Mask => " gap-mode=mask(-m)".to_string(),
            GapMode::Closed => " gap-mode=1".to_string(),
            GapMode::NoGaps => " gap-mode=2".to_string(),
            GapMode::Partial => String::new(),
        }
    } else {
        String::new()
    }
}

fn engine_kind_of(engine: dynamite::gene::Engine) -> Option<EngineKind> {
    use dynamite::gene::Engine as E;
    Some(match engine {
        E::Prodigal => EngineKind::Prodigal,
        E::ProdigalGv => EngineKind::ProdigalGv,
        E::ProdigalMeta => EngineKind::ProdigalMeta,
        _ => return None,
    })
}

fn cmd_call(a: CallArgs) -> Result<()> {
    let mode = Mode::parse(&a.mode).ok_or_else(|| {
        anyhow!(
            "unknown mode '{}' (try: phage, virus, giant-virus, crassphage, bacteria, archaea, meta, eukaryote, reads, auto)",
            a.mode
        )
    })?;
    let (engine, gap, gap_given) = parse_engine_and_gaps(&a.engine, a.gap_mode.as_deref(), a.mask)?;
    let format = Format::parse(&a.format)
        .ok_or_else(|| anyhow!("unknown format '{}' (try: gff3, gtf, genbank, faa, fna, ffn, tsv)", a.format))?;
    let code = parse_code(a.code.as_deref())?;
    if a.circular && a.closed {
        return Err(anyhow!("--circular and --closed are mutually exclusive"));
    }
    if a.circular {
        eprintln!("dynamite: note: --circular is accepted for compatibility; genes are not wrapped across the origin");
    }

    let (_path, records) = read_input(a.input, a.input_pos)?;
    warn_run(mode, engine, &code, &records, gap_given);
    let cfg = RunConfig {
        mode,
        engine,
        code,
        min_len: a.minlen,
        closed: a.closed,
        threads: a.threads,
    };
    let out = run_records_with(&records, &cfg, &RunOptions::default().with_gap_mode(gap))?;

    let mut buf = String::new();
    if format == Format::Genbank {
        // A GenBank record needs its sequence (ORIGIN); the coordinate-only
        // writer produced records that GenBank parsers reject.
        let seqs: Vec<(String, Vec<u8>)> = records
            .iter()
            .map(|r| (r.id.clone(), r.seq.iter().map(|b| b.to_ascii_uppercase()).collect()))
            .collect();
        dynamite::gene::write_genbank_full(&mut buf, &out.genes, &seqs);
    } else {
        dynamite::gene::write_all(&mut buf, &out.genes, format, &out.contigs);
    }

    eprintln!(
        "dynamite: mode={} engine={} code={} contigs={} genes={}{}",
        mode.label(),
        out.engine.label(),
        out.trans_table,
        out.contigs.len(),
        out.genes.len(),
        engine_kind_of(out.engine).map(|k| gap_note(k, gap)).unwrap_or_default(),
    );
    write_out(&a.output, &buf)
}

fn cmd_all(a: AllArgs) -> Result<()> {
    let mode = Mode::parse(&a.mode)
        .ok_or_else(|| anyhow!("unknown mode '{}' (see `dynamite call --help`)", a.mode))?;
    let (engine, gap, gap_given) = parse_engine_and_gaps(&a.engine, a.gap_mode.as_deref(), a.mask)?;
    let code = parse_code(a.code.as_deref())?;
    let (_path, records) = read_input(a.input, a.input_pos)?;
    warn_run(mode, engine, &code, &records, gap_given);
    let cfg = RunConfig {
        mode,
        engine,
        code,
        min_len: a.minlen,
        closed: a.closed,
        threads: a.threads,
    };
    let out = run_records_with(&records, &cfg, &RunOptions::default().with_gap_mode(gap))?;

    std::fs::create_dir_all(&a.outdir)
        .map_err(|e| anyhow!("cannot create outdir {:?}: {e}", a.outdir))?;

    // Raw input contigs (uppercased) for .fna and the GenBank ORIGIN.
    let contigs_seq: Vec<(String, Vec<u8>)> = records
        .iter()
        .map(|r| (r.id.clone(), r.seq.iter().map(|b| b.to_ascii_uppercase()).collect()))
        .collect();

    let write = |ext: &str, content: &str| -> Result<()> {
        let p = a.outdir.join(format!("{}.{}", a.prefix, ext));
        std::fs::write(&p, content).map_err(|e| anyhow!("write {:?}: {e}", p))
    };

    use dynamite::gene::{write_all, write_genbank_full};
    let mut s = String::new();
    write_all(&mut s, &out.genes, Format::Gff3, &out.contigs);
    write("gff", &s)?;
    s.clear();
    write_all(&mut s, &out.genes, Format::Gtf, &out.contigs);
    write("gtf", &s)?;
    s.clear();
    write_all(&mut s, &out.genes, Format::Faa, &out.contigs);
    write("faa", &s)?;
    s.clear();
    write_all(&mut s, &out.genes, Format::Fna, &out.contigs);
    write("ffn", &s)?;
    s.clear();
    write_all(&mut s, &out.genes, Format::Tsv, &out.contigs);
    write("tsv", &s)?;
    s.clear();
    write_genbank_full(&mut s, &out.genes, &contigs_seq);
    write("gbk", &s)?;

    // .fna = raw input contig sequences
    s.clear();
    for (id, seq) in &contigs_seq {
        s.push('>');
        s.push_str(id);
        s.push('\n');
        let mut i = 0;
        while i < seq.len() {
            let e = (i + 70).min(seq.len());
            s.push_str(std::str::from_utf8(&seq[i..e]).unwrap_or(""));
            s.push('\n');
            i = e;
        }
    }
    write("fna", &s)?;

    eprintln!(
        "dynamite all: engine={}{} contigs={} genes={} -> {}/{}.{{gff,gtf,gbk,faa,ffn,fna,tsv}}",
        out.engine.label(),
        engine_kind_of(out.engine).map(|k| gap_note(k, gap)).unwrap_or_default(),
        out.contigs.len(),
        out.genes.len(),
        a.outdir.display(),
        a.prefix,
    );
    Ok(())
}

fn cmd_annotate(a: AnnotateArgs) -> Result<()> {
    let mode = Mode::parse(&a.mode)
        .ok_or_else(|| anyhow!("unknown mode '{}' (see `dynamite call --help`)", a.mode))?;
    let (engine, gap, gap_given) = parse_engine_and_gaps(&a.engine, a.gap_mode.as_deref(), a.mask)?;
    let code = parse_code(a.code.as_deref())?;
    let (_path, records) = read_input(a.input, a.input_pos)?;
    warn_run(mode, engine, &code, &records, gap_given);

    let locus_prefix = a
        .locus
        .clone()
        .unwrap_or_else(|| a.prefix.to_ascii_uppercase());
    let rna_code = match &code {
        CodeChoice::Fixed(n) => *n,
        _ => 11,
    };
    let cfg = AnnotateConfig {
        run: RunConfig {
            mode,
            engine,
            code,
            min_len: 0,
            closed: a.closed,
            threads: a.threads,
        },
        locus_prefix,
        organism: a.organism.clone(),
        strain: a.strain.clone(),
        detect_rna: !a.no_rna,
        rna_code,
    };

    let bundle = annotate_with(&records, &cfg, &RunOptions::default().with_gap_mode(gap))?;

    std::fs::create_dir_all(&a.outdir)
        .map_err(|e| anyhow!("cannot create outdir {:?}: {e}", a.outdir))?;
    let write = |ext: &str, content: &str| -> Result<()> {
        let p = a.outdir.join(format!("{}.{}", a.prefix, ext));
        std::fs::write(&p, content).map_err(|e| anyhow!("write {:?}: {e}", p))
    };
    write("gff", &bundle.gff)?;
    write("gbk", &bundle.gbk)?;
    write("fna", &bundle.fna)?;
    write("faa", &bundle.faa)?;
    write("ffn", &bundle.ffn)?;
    write("tsv", &bundle.tsv)?;
    write("txt", &bundle.txt)?;
    write("tbl", &bundle.tbl)?;
    write("fsa", &bundle.fsa)?;
    write("err", &bundle.err)?;
    write("log", &bundle.log)?;

    eprintln!(
        "dynamite annotate: {} CDS, {} tRNA, {} rRNA -> {}/{}.{{gff,gbk,fna,faa,ffn,tsv,txt,tbl,fsa,err,log}}",
        bundle.n_cds,
        bundle.n_trna,
        bundle.n_rrna,
        a.outdir.display(),
        a.prefix,
    );
    eprintln!(
        "note: for the NCBI .sqn (+ official .err), run:  table2asn -i {0}/{1}.fsa -f {0}/{1}.tbl -o {0}/{1}.sqn",
        a.outdir.display(),
        a.prefix,
    );
    Ok(())
}

fn cmd_translate(a: TranslateArgs) -> Result<()> {
    let code = GeneticCode::new(a.code)
        .ok_or_else(|| anyhow!("unsupported genetic code {}", a.code))?;
    let (_path, records) = read_input(a.input, a.input_pos)?;

    let mut buf = String::new();
    if a.orfs {
        let format = Format::parse(&a.format)
            .ok_or_else(|| anyhow!("unknown format '{}' (try: gff3, faa, fna, tsv)", a.format))?;
        let mut genes = Vec::new();
        let mut contigs = Vec::new();
        for r in &records {
            let s = dynamite::seq::Seq::from_ascii(&r.seq);
            genes.extend(sixframe::orfs(&r.id, &s, &code, a.minlen, a.require_start));
            contigs.push((r.id.clone(), r.seq.len()));
        }
        if format == Format::Genbank {
            let seqs: Vec<(String, Vec<u8>)> = records
                .iter()
                .map(|r| (r.id.clone(), r.seq.iter().map(|b| b.to_ascii_uppercase()).collect()))
                .collect();
            dynamite::gene::write_genbank_full(&mut buf, &genes, &seqs);
        } else {
            dynamite::gene::write_all(&mut buf, &genes, format, &contigs);
        }
        eprintln!("dynamite: translate --orfs code={} contigs={} orfs={}", a.code, records.len(), genes.len());
    } else {
        for r in &records {
            let s = dynamite::seq::Seq::from_ascii(&r.seq);
            let frames = sixframe::translate_six(&s, &code);
            buf.push_str(&sixframe::to_fasta(&r.id, &frames));
        }
        eprintln!("dynamite: translate code={} contigs={} (six frames each)", a.code, records.len());
    }
    write_out(&a.output, &buf)
}

fn cmd_codon(a: CodonArgs) -> Result<()> {
    let code = GeneticCode::new(a.code)
        .ok_or_else(|| anyhow!("unsupported genetic code {}", a.code))?;
    let (_path, records) = read_input(a.input, a.input_pos)?;

    let cds: Vec<String> = if a.cds {
        records
            .iter()
            .map(|r| String::from_utf8_lossy(&r.seq).to_string())
            .collect()
    } else {
        let mode = Mode::parse(&a.mode)
            .ok_or_else(|| anyhow!("unknown mode '{}'", a.mode))?;
        let cfg = RunConfig {
            mode,
            engine: None,
            code: CodeChoice::Fixed(a.code),
            min_len: 0,
            closed: false,
            threads: 0,
        };
        let out = run_records(&records, &cfg)?;
        eprintln!(
            "dynamite: codon — called {} genes ({}), computing usage",
            out.genes.len(),
            out.engine.label()
        );
        out.genes.iter().map(|g| g.nt.clone()).collect()
    };

    let stats = codon::analyze(&cds, &code);
    let buf = stats.to_tsv(&code);
    write_out(&a.output, &buf)
}

fn cmd_rna(a: RnaArgs) -> Result<()> {
    let code = GeneticCode::new(a.code)
        .ok_or_else(|| anyhow!("unsupported genetic code {}", a.code))?;
    let (_path, records) = read_input(a.input, a.input_pos)?;

    // One GFF3 header for the whole file (it used to be repeated per contig).
    let mut buf = String::from("##gff-version 3\n");
    let mut total = 0usize;
    for r in &records {
        let s = dynamite::seq::Seq::from_ascii(&r.seq);
        let feats = rna::find_all(&s, &code);
        total += feats.len();
        buf.push_str(&rna::to_gff_records(&r.id, &feats));
    }
    eprintln!(
        "dynamite: rna — {} features across {} contigs (heuristic; verify with tRNAscan-SE / barrnap)",
        total,
        records.len()
    );
    write_out(&a.output, &buf)
}
