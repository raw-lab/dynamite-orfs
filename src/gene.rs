//! The engine-independent gene representation and all output formats.

use crate::gencode::GeneticCode;
use crate::seq::{decode_base, revcomp, Seq};
use std::fmt::Write as _;

/// Which engine produced a call.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Engine {
    /// PHANOTATE-style weighted-graph / shortest-path engine (phages).
    Phanotate,
    /// Prodigal-style self-training dynamic-programming engine (generic).
    Dynamic,
    /// Prodigal preset for bacteria / archaea (code 11, single genome).
    Prodigal,
    /// Prodigal-GV preset for phages, giant viruses and crassphages.
    ProdigalGv,
    /// Prodigal metagenomic preset (fragmented, mixed-code input).
    ProdigalMeta,
    /// GeneMark-style 3-periodic inhomogeneous Markov gene finder.
    GeneMark,
    /// AUGUSTUS-style eukaryotic (intron-aware) gene finder.
    Augustus,
    /// FragGeneScan codon HMM (vendored FragGeneScanRs; short and long reads).
    FragGeneScan,
    /// Plain six-frame translation (every ORF between stops).
    SixFrame,
}

impl Engine {
    /// Short label used in output attributes.
    pub fn label(self) -> &'static str {
        match self {
            Engine::Phanotate => "phanotate",
            Engine::Dynamic => "dynamic",
            Engine::Prodigal => "prodigal",
            Engine::ProdigalGv => "prodigalgv",
            Engine::ProdigalMeta => "prodigal-meta",
            Engine::GeneMark => "genemark",
            Engine::Augustus => "augustus",
            Engine::FragGeneScan => "fraggenescan",
            Engine::SixFrame => "sixframe",
        }
    }
}

/// A predicted gene / ORF on a contig. Coordinates are 1-based inclusive.
#[derive(Clone)]
pub struct Gene {
    /// Contig id this gene belongs to.
    pub contig: String,
    /// 1-based gene index within the contig.
    pub id: usize,
    /// Leftmost (smallest) coordinate, inclusive.
    pub begin: usize,
    /// Rightmost (largest) coordinate, inclusive (includes the stop codon for
    /// complete genes).
    pub end: usize,
    /// `+1` forward, `-1` reverse.
    pub strand: i8,
    /// Start-codon type: 0=ATG, 1=GTG, 2=TTG, 3=edge/none.
    pub start_type: u8,
    /// RBS motif text (e.g. `AGGAGG`) or `None`.
    pub rbs_motif: String,
    /// RBS spacer text (e.g. `5-10bp`) or `None`.
    pub rbs_spacer: String,
    /// The gene is incomplete at its **left** (lower-coordinate, `begin`) end —
    /// Prodigal's genomic convention, independent of strand. On the forward
    /// strand that is the 5' end (no start codon); on the reverse strand it is
    /// the 3' end (no stop codon). Use [`Gene::partial_5p`] /
    /// [`Gene::partial_3p`] for the strand-relative view.
    pub partial_left: bool,
    /// The gene is incomplete at its **right** (higher-coordinate, `end`) end
    /// (genomic convention, see `partial_left`).
    pub partial_right: bool,
    /// Total model score for the gene (engine-specific scale).
    pub score: f64,
    /// Coding-potential component (dynamic engine; 0 for PHANOTATE).
    pub cscore: f64,
    /// Start-signal component (dynamic engine; 0 for PHANOTATE).
    pub sscore: f64,
    /// Confidence percentage (dynamic engine).
    pub conf: f64,
    /// GC fraction of the gene.
    pub gc: f64,
    /// NCBI translation table used.
    pub trans_table: u32,
    /// Engine that produced the call.
    pub engine: Engine,
    /// Protein translation (filled by [`Gene::fill_sequences`]).
    pub aa: String,
    /// Nucleotide coding sequence, 5'->3' (filled by [`Gene::fill_sequences`]).
    pub nt: String,
}

impl Gene {
    /// Reading-frame partial flag in Prodigal's `XY` form: `X` = left (begin)
    /// end, `Y` = right (end) end, `1` = incomplete — for **both** strands.
    ///
    /// (Previously the two flags were swapped for reverse-strand genes, so
    /// `partial=` was written in 5'/3' order on the minus strand.)
    pub fn partial_str(&self) -> String {
        format!("{}{}", self.partial_left as u8, self.partial_right as u8)
    }

    /// `true` if the gene's 5' end (start-codon side) is incomplete.
    #[inline]
    pub fn partial_5p(&self) -> bool {
        if self.strand >= 0 {
            self.partial_left
        } else {
            self.partial_right
        }
    }

    /// `true` if the gene's 3' end (stop-codon side) is incomplete.
    #[inline]
    pub fn partial_3p(&self) -> bool {
        if self.strand >= 0 {
            self.partial_right
        } else {
            self.partial_left
        }
    }

    /// Set the strand-relative partial flags (stored in genomic orientation).
    pub fn set_partial_5p_3p(&mut self, p5: bool, p3: bool) {
        if self.strand >= 0 {
            self.partial_left = p5;
            self.partial_right = p3;
        } else {
            self.partial_left = p3;
            self.partial_right = p5;
        }
    }

    /// GenBank / INSDC location string, with `<` / `>` on the incomplete
    /// genomic ends (e.g. `complement(<1..>300)`).
    pub fn insdc_location(&self) -> String {
        let lb = if self.partial_left { "<" } else { "" };
        let rb = if self.partial_right { ">" } else { "" };
        if self.strand >= 0 {
            format!("{lb}{}..{rb}{}", self.begin, self.end)
        } else {
            format!("complement({lb}{}..{rb}{})", self.begin, self.end)
        }
    }

    /// Start-codon text.
    pub fn start_codon_str(&self) -> &'static str {
        match self.start_type {
            0 => "ATG",
            1 => "GTG",
            2 => "TTG",
            _ => "Edge",
        }
    }

    /// Extract the nucleotide coding sequence and translate it to protein,
    /// using the supplied genetic code. The leading start codon is rendered as
    /// `M` for genes with a complete 5' end (matching Prodigal), and the
    /// trailing stop codon is dropped from the protein. Ambiguous bases are
    /// written as `N` and any codon containing one translates to `X` — never
    /// as a concrete base (N used to be read as A, e.g. `NNN` → `K`).
    pub fn fill_sequences(&mut self, seq: &Seq, code: &GeneticCode) {
        let b = self.begin.saturating_sub(1);
        let e = self.end.min(seq.len());
        if b >= e {
            return;
        }
        let (coding, amb): (Vec<u8>, Vec<bool>) = if self.strand >= 0 {
            (seq.fwd[b..e].to_vec(), seq.n[b..e].to_vec())
        } else {
            (revcomp(&seq.fwd[b..e]), seq.n[b..e].iter().rev().copied().collect())
        };

        // Nucleotide string.
        let mut nt = String::with_capacity(coding.len());
        for (&c, &a) in coding.iter().zip(amb.iter()) {
            nt.push(if a { 'N' } else { decode_base(c) as char });
        }
        self.nt = nt;

        // Protein string.
        let ncodons = coding.len() / 3;
        let p5 = self.partial_5p();
        let p3 = self.partial_3p();
        let mut aa = String::with_capacity(ncodons);
        for k in 0..ncodons {
            let i = k * 3;
            let (c1, c2, c3) = (coding[i], coding[i + 1], coding[i + 2]);
            let ambiguous = amb[i] || amb[i + 1] || amb[i + 2];
            let residue = if ambiguous {
                b'X'
            } else if k == 0 && !p5 && (code.is_start(c1, c2, c3) || GeneticCode::start_type(c1, c2, c3) < 3) {
                // A complete 5' end begins with the gene's start codon, read
                // as Met (PHANOTATE opens genes at ATG/GTG/TTG under any table).
                b'M'
            } else {
                code.translate(c1, c2, c3)
            };
            // Drop the final stop codon from the protein for complete 3' ends.
            if k == ncodons - 1 && residue == b'*' && !p3 {
                break;
            }
            aa.push(residue as char);
        }
        self.aa = aa;
    }

    /// Protein length (residues).
    pub fn protein_len(&self) -> usize {
        self.aa.len()
    }
}

/// Output formats DYNAMITE can emit.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Format {
    /// GFF3 (default).
    Gff3,
    /// GTF / GFF2 (gene_id / transcript_id attributes).
    Gtf,
    /// GenBank flat file with CDS features.
    Genbank,
    /// Protein FASTA.
    Faa,
    /// Nucleotide (CDS) FASTA — also written for the `ffn` alias.
    Fna,
    /// Tab-separated PHANOTATE-style coordinate table.
    Tsv,
}

impl Format {
    /// Parse a format name (case-insensitive). Accepts common aliases.
    pub fn parse(s: &str) -> Option<Format> {
        match s.to_ascii_lowercase().as_str() {
            "gff" | "gff3" => Some(Format::Gff3),
            "gtf" | "gff2" => Some(Format::Gtf),
            "genbank" | "gbk" | "gb" => Some(Format::Genbank),
            "faa" | "protein" | "proteins" => Some(Format::Faa),
            "fna" | "ffn" | "nucleotide" | "genes" | "fasta" => Some(Format::Fna),
            "tsv" | "tabular" | "coords" | "tab" => Some(Format::Tsv),
            _ => None,
        }
    }
}

/// A GenBank `LOCUS` line. The name and length are always separated by
/// whitespace (a contig id of 16+ characters used to run straight into the
/// length, e.g. `NC_001416.1_lambda48502 bp`, which parsers reject).
pub fn locus_line(id: &str, len: usize) -> String {
    format!("LOCUS       {id:<16} {len:>11} bp    DNA     linear   UNK 01-JAN-2025")
}

/// FASTA-style descriptive header attributes shared by several formats.
fn attrs(g: &Gene) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        "ID={contig}_{id};partial={partial};start_type={st};rbs_motif={rm};rbs_spacer={rs};gc_cont={gc:.3};conf={conf:.2};score={score:.2};engine={eng};transl_table={tt}",
        contig = g.contig,
        id = g.id,
        partial = g.partial_str(),
        st = g.start_codon_str(),
        rm = g.rbs_motif,
        rs = g.rbs_spacer,
        gc = g.gc,
        conf = g.conf,
        score = g.score,
        eng = g.engine.label(),
        tt = g.trans_table,
    );
    s
}

/// Write all genes (already grouped/ordered by the caller) in the chosen
/// format to `out`.
pub fn write_all(out: &mut String, genes: &[Gene], format: Format, contigs: &[(String, usize)]) {
    match format {
        Format::Gff3 => write_gff3(out, genes, contigs),
        Format::Gtf => write_gtf(out, genes),
        Format::Genbank => write_genbank(out, genes, contigs),
        Format::Faa => write_faa(out, genes),
        Format::Fna => write_fna(out, genes),
        Format::Tsv => write_tsv(out, genes),
    }
}

fn write_gtf(out: &mut String, genes: &[Gene]) {
    for g in genes {
        let strand = if g.strand >= 0 { '+' } else { '-' };
        let _ = writeln!(
            out,
            "{seqid}\tDYNAMITE\tCDS\t{start}\t{end}\t{score:.2}\t{strand}\t0\tgene_id \"{contig}_{id}\"; transcript_id \"{contig}_{id}.t1\"; start_type \"{st}\"; partial \"{partial}\"; engine \"{eng}\"; transl_table \"{tt}\";",
            seqid = g.contig,
            start = g.begin,
            end = g.end,
            score = g.score,
            strand = strand,
            contig = g.contig,
            id = g.id,
            st = g.start_codon_str(),
            partial = g.partial_str(),
            eng = g.engine.label(),
            tt = g.trans_table,
        );
    }
}

fn write_gff3(out: &mut String, genes: &[Gene], contigs: &[(String, usize)]) {
    out.push_str("##gff-version 3\n");
    for (id, len) in contigs {
        let _ = writeln!(out, "##sequence-region {id} 1 {len}");
    }
    for g in genes {
        let strand = if g.strand >= 0 { '+' } else { '-' };
        let _ = writeln!(
            out,
            "{seqid}\tDYNAMITE\tCDS\t{start}\t{end}\t{score:.2}\t{strand}\t0\t{attrs}",
            seqid = g.contig,
            start = g.begin,
            end = g.end,
            score = g.score,
            strand = strand,
            attrs = attrs(g),
        );
    }
}

fn fasta_wrap(out: &mut String, seq: &str) {
    const W: usize = 60;
    let bytes = seq.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let end = (i + W).min(bytes.len());
        out.push_str(std::str::from_utf8(&bytes[i..end]).unwrap_or(""));
        out.push('\n');
        i = end;
    }
    if bytes.is_empty() {
        out.push('\n');
    }
}

fn write_faa(out: &mut String, genes: &[Gene]) {
    for g in genes {
        let strand = if g.strand >= 0 { 1 } else { -1 };
        let _ = writeln!(
            out,
            ">{contig}_{id} # {begin} # {end} # {strand} # {attrs}",
            contig = g.contig,
            id = g.id,
            begin = g.begin,
            end = g.end,
            strand = strand,
            attrs = attrs(g),
        );
        fasta_wrap(out, &g.aa);
    }
}

fn write_fna(out: &mut String, genes: &[Gene]) {
    for g in genes {
        let strand = if g.strand >= 0 { 1 } else { -1 };
        let _ = writeln!(
            out,
            ">{contig}_{id} # {begin} # {end} # {strand} # {attrs}",
            contig = g.contig,
            id = g.id,
            begin = g.begin,
            end = g.end,
            strand = strand,
            attrs = attrs(g),
        );
        fasta_wrap(out, &g.nt);
    }
}

fn write_tsv(out: &mut String, genes: &[Gene]) {
    out.push_str("#START\tSTOP\tFRAME\tCONTIG\tSCORE\tSTART_TYPE\tRBS_MOTIF\tPARTIAL\tENGINE\n");
    for g in genes {
        // PHANOTATE convention: START is the 5' coordinate, STOP the 3'.
        let (start, stop, frame) = if g.strand >= 0 {
            (g.begin, g.end, "+")
        } else {
            (g.end, g.begin, "-")
        };
        let _ = writeln!(
            out,
            "{start}\t{stop}\t{frame}\t{contig}\t{score:E}\t{st}\t{rm}\t{partial}\t{eng}",
            start = start,
            stop = stop,
            frame = frame,
            contig = g.contig,
            score = g.score,
            st = g.start_codon_str(),
            rm = g.rbs_motif,
            partial = g.partial_str(),
            eng = g.engine.label(),
        );
    }
}

fn write_genbank(out: &mut String, genes: &[Gene], contigs: &[(String, usize)]) {
    use std::collections::HashMap;
    let mut by_contig: HashMap<&str, Vec<&Gene>> = HashMap::new();
    for g in genes {
        by_contig.entry(g.contig.as_str()).or_default().push(g);
    }
    for (id, len) in contigs {
        let _ = writeln!(out, "{}", locus_line(id, *len));
        let _ = writeln!(out, "DEFINITION  {id} DYNAMITE gene predictions.");
        out.push_str("FEATURES             Location/Qualifiers\n");
        let _ = writeln!(out, "     source          1..{len}");
        if let Some(gs) = by_contig.get(id.as_str()) {
            for g in gs {
                let loc = g.insdc_location();
                let _ = writeln!(out, "     CDS             {loc}");
                let _ = writeln!(out, "                     /note=\"score:{:E};engine:{}\"", g.score, g.engine.label());
                let _ = writeln!(out, "                     /transl_table={}", g.trans_table);
                let _ = writeln!(out, "                     /translation=\"{}\"", g.aa);
            }
        }
        out.push_str("//\n");
    }
}

/// Write a complete GenBank flat file (CDS features + wrapped translation +
/// the `ORIGIN` sequence) for raw gene predictions. Unlike [`write_genbank`],
/// this takes the contig **sequences** so the records are self-contained.
pub fn write_genbank_full(out: &mut String, genes: &[Gene], contigs: &[(String, Vec<u8>)]) {
    use std::collections::HashMap;
    let mut by_contig: HashMap<&str, Vec<&Gene>> = HashMap::new();
    for g in genes {
        by_contig.entry(g.contig.as_str()).or_default().push(g);
    }
    for (id, seq) in contigs {
        let len = seq.len();
        let _ = writeln!(out, "{}", locus_line(id, len));
        let _ = writeln!(out, "DEFINITION  {id} DYNAMITE gene predictions.");
        let _ = writeln!(out, "ACCESSION   {id}");
        out.push_str("FEATURES             Location/Qualifiers\n");
        let _ = writeln!(out, "     source          1..{len}");
        if let Some(gs) = by_contig.get(id.as_str()) {
            for g in gs {
                let loc = g.insdc_location();
                let _ = writeln!(out, "     CDS             {loc}");
                let _ = writeln!(out, "                     /note=\"score:{:E};engine:{}\"", g.score, g.engine.label());
                let _ = writeln!(out, "                     /transl_table={}", g.trans_table);
                // wrapped translation
                let aa = g.aa.as_bytes();
                if aa.is_empty() {
                    let _ = writeln!(out, "                     /translation=\"\"");
                } else {
                    let mut i = 0;
                    let mut first = true;
                    while i < aa.len() {
                        let e = (i + 58).min(aa.len());
                        let chunk = std::str::from_utf8(&aa[i..e]).unwrap_or("");
                        if first {
                            let _ = writeln!(out, "                     /translation=\"{chunk}");
                            first = false;
                        } else {
                            let _ = writeln!(out, "                     {chunk}");
                        }
                        i = e;
                    }
                    let trimmed = out.trim_end_matches('\n').len();
                    out.truncate(trimmed);
                    out.push_str("\"\n");
                }
            }
        }
        out.push_str("ORIGIN\n");
        let lower = seq.to_ascii_lowercase();
        let mut pos = 0;
        while pos < lower.len() {
            let line_end = (pos + 60).min(lower.len());
            let mut line = String::new();
            let _ = write!(line, "{:>9}", pos + 1);
            let mut j = pos;
            while j < line_end {
                let g_end = (j + 10).min(line_end);
                line.push(' ');
                line.push_str(std::str::from_utf8(&lower[j..g_end]).unwrap_or(""));
                j = g_end;
            }
            let _ = writeln!(out, "{line}");
            pos = line_end;
        }
        out.push_str("//\n");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_gene() -> Gene {
        Gene {
            contig: "c1".into(),
            id: 1,
            begin: 1,
            end: 9,
            strand: 1,
            start_type: 0,
            rbs_motif: "None".into(),
            rbs_spacer: "None".into(),
            partial_left: false,
            partial_right: false,
            score: -12.3,
            cscore: 0.0,
            sscore: 0.0,
            conf: 99.0,
            gc: 0.5,
            trans_table: 11,
            engine: Engine::Dynamic,
            aa: String::new(),
            nt: String::new(),
        }
    }

    #[test]
    fn translate_simple_gene() {
        // ATG AAA TAA  -> M K (stop dropped)
        let s = Seq::from_ascii(b"ATGAAATAA");
        let code = GeneticCode::new(11).unwrap();
        let mut g = dummy_gene();
        g.fill_sequences(&s, &code);
        assert_eq!(g.aa, "MK");
        assert_eq!(g.nt, "ATGAAATAA");
    }

    #[test]
    fn genbank_full_emits_origin() {
        let s = Seq::from_ascii(b"ATGAAATAA");
        let code = GeneticCode::new(11).unwrap();
        let mut g = dummy_gene();
        g.fill_sequences(&s, &code);
        let contigs = vec![("c1".to_string(), b"ATGAAATAA".to_vec())];
        let mut out = String::new();
        write_genbank_full(&mut out, &[g], &contigs);
        assert!(out.contains("LOCUS       c1"));
        assert!(out.contains("/translation=\"MK\""));
        assert!(out.contains("ORIGIN"));
        assert!(out.contains("atgaaataa"));
        assert!(out.trim_end().ends_with("//"));
    }

    #[test]
    fn translate_reverse_gene() {
        // forward TTA TTT CAT ; reverse-complement = ATG AAA TAA -> M K
        let s = Seq::from_ascii(b"TTATTTCAT");
        let code = GeneticCode::new(11).unwrap();
        let mut g = dummy_gene();
        g.strand = -1;
        g.fill_sequences(&s, &code);
        assert_eq!(g.aa, "MK");
        assert_eq!(g.nt, "ATGAAATAA");
    }

    #[test]
    fn format_parsing() {
        assert_eq!(Format::parse("GFF3"), Some(Format::Gff3));
        assert_eq!(Format::parse("faa"), Some(Format::Faa));
        assert_eq!(Format::parse("nonsense"), None);
    }
}
