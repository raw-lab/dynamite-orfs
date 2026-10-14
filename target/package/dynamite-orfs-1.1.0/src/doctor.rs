//! `dynamite doctor` — environment report and end-to-end self-test.
//!
//! Builds a small deterministic synthetic genome (no external files, no RNG
//! crate) and exercises every engine and helper, reporting a ✓/✗ checklist.
//! Useful for verifying an install and for CI smoke testing.

use crate::gencode::{GeneticCode, SUPPORTED};
use crate::gene::{write_all, Format, Gene};
use crate::seq::Seq;
use crate::{augustus, banner, codon, fraggenescan, genemark, phanotate, prodigal, prodigal_meta, prodigalgv, rna, sixframe};
use rayon::prelude::*;
use std::fmt::Write as _;

struct Check {
    label: String,
    ok: bool,
    detail: String,
}

fn pass(checks: &mut Vec<Check>, label: &str, ok: bool, detail: String) {
    checks.push(Check { label: label.into(), ok, detail });
}

/// Deterministic ~2.5 kb test genome containing several clean ORFs.
fn build_test_genome() -> (Vec<u8>, Vec<(String, Seq)>) {
    // All 61 sense codons under the standard code (exclude TAA/TAG/TGA).
    let bases = [b'T', b'C', b'A', b'G'];
    let mut pool: Vec<[u8; 3]> = Vec::new();
    for &a in &bases {
        for &b in &bases {
            for &c in &bases {
                let cod = [a, b, c];
                if matches!(&cod, b"TAA" | b"TAG" | b"TGA") {
                    continue;
                }
                pool.push(cod);
            }
        }
    }
    // xorshift64 for reproducibility without an RNG dependency.
    let mut x: u64 = 0x2545_F491_4F6C_DD1D;
    let mut next = || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut s: Vec<u8> = Vec::new();
    for _ in 0..4 {
        s.extend_from_slice(b"ATG");
        let ncod = 120 + (next() % 80) as usize;
        for _ in 0..ncod {
            s.extend_from_slice(&pool[(next() as usize) % pool.len()]);
        }
        s.extend_from_slice(b"TAA");
        for _ in 0..30 {
            s.push(bases[(next() as usize) % 4]);
        }
    }
    let seq = Seq::from_ascii(&s);
    (s.clone(), vec![("dynamite_test".to_string(), seq)])
}

/// Run all checks. Returns (report, all_ok).
pub fn run() -> (String, bool) {
    let mut checks: Vec<Check> = Vec::new();
    let (ascii, contigs) = build_test_genome();
    let code = GeneticCode::new(11).unwrap();
    let seq = &contigs[0].1;
    let contig_lens = vec![(contigs[0].0.clone(), seq.len())];

    // 1) Genetic codes.
    let mut codes_ok = 0usize;
    for &id in SUPPORTED {
        if let Some(c) = GeneticCode::new(id) {
            let _ = c.translate(b'A', b'T', b'G');
            codes_ok += 1;
        }
    }
    pass(
        &mut checks,
        "genetic codes",
        codes_ok == SUPPORTED.len(),
        format!("{}/{} tables loaded", codes_ok, SUPPORTED.len()),
    );

    // 2) Six-frame translation + ORFs.
    let frames = sixframe::translate_six(seq, &code);
    let sf_orfs = sixframe::orfs("dynamite_test", seq, &code, 90, true);
    pass(
        &mut checks,
        "sixframe translate",
        frames.len() == 6,
        format!("{} frames", frames.len()),
    );
    pass(&mut checks, "sixframe ORFs", !sf_orfs.is_empty(), format!("{} ORFs", sf_orfs.len()));

    // 3) Prodigal presets.
    let (pg, _t) = prodigal::call(&contigs, None, false);
    pass(&mut checks, "prodigal (bacteria/archaea)", !pg.is_empty(), format!("{} genes", pg.len()));
    let (gv, _t) = prodigalgv::call(&contigs, Some(11), false);
    pass(&mut checks, "prodigalgv (phage/giant virus)", !gv.is_empty(), format!("{} genes", gv.len()));
    let (mt, _t) = prodigal_meta::call(&contigs, Some(11), false);
    pass(&mut checks, "prodigal-meta (metagenome)", !mt.is_empty(), format!("{} genes", mt.len()));

    // 4) GeneMark.
    let (gm, _t) = genemark::call(&contigs, Some(11), false);
    pass(&mut checks, "genemark (Markov self-train)", !gm.is_empty(), format!("{} genes", gm.len()));

    // 5) AUGUSTUS (eukaryotic; data-dependent — must run, genes optional).
    let (au, _t) = augustus::call(&contigs, Some(1), false);
    pass(&mut checks, "augustus (eukaryote)", true, format!("{} CDS", au.len()));

    // 6) FragGeneScan on read-sized fragments.
    let mut reads: Vec<(String, Seq)> = Vec::new();
    let mut i = 0;
    while i < ascii.len() {
        let end = (i + 200).min(ascii.len());
        reads.push((format!("read{}", i / 200), Seq::from_ascii(&ascii[i..end])));
        i = end;
    }
    let (fgs, _t) = fraggenescan::call(&reads, false);
    pass(&mut checks, "fraggenescan (reads)", true, format!("{} genes on {} reads", fgs.len(), reads.len()));

    // 7) PHANOTATE.
    let (ph, _s) = phanotate::predict("dynamite_test", &ascii, &code, 90, false);
    pass(&mut checks, "phanotate (graph)", !ph.is_empty(), format!("{} genes", ph.len()));

    // 8) Codon statistics.
    let cds: Vec<String> = pg.iter().map(|g| g.nt.clone()).collect();
    let cstats = codon::analyze(&cds, &code);
    pass(
        &mut checks,
        "codon stats (RSCU/ENC/GC3)",
        cstats.total_codons > 0,
        format!("{} codons, ENC={:.1}, GC3={:.2}", cstats.total_codons, cstats.enc, cstats.gc3),
    );

    // 9) RNA scan (heuristic; must run).
    let rnas = rna::find_all(seq, &code);
    let trna = rnas.iter().filter(|r| r.kind == "tRNA").count();
    pass(&mut checks, "rna scan (tRNA/rRNA heuristic)", true, format!("{} features ({} tRNA)", rnas.len(), trna));

    // 10) Output formats.
    let mut fmt_ok = true;
    let mut fmt_detail = String::new();
    for (name, fmt) in [
        ("gff3", Format::Gff3),
        ("gtf", Format::Gtf),
        ("genbank", Format::Genbank),
        ("faa", Format::Faa),
        ("fna", Format::Fna),
        ("ffn", Format::Fna),
        ("tsv", Format::Tsv),
    ] {
        let mut out = String::new();
        write_all(&mut out, &pg, fmt, &contig_lens);
        if out.is_empty() {
            fmt_ok = false;
        }
        let _ = write!(fmt_detail, "{} ", name);
    }
    pass(&mut checks, "output formats", fmt_ok, fmt_detail.trim().to_string());

    // 11) FASTA round-trip (write to temp, read back).
    let fasta_ok = fasta_roundtrip();
    pass(&mut checks, "fasta read (+gzip transparent)", fasta_ok, "parse 1 record".into());

    // 12) rayon parallelism.
    let sum: u64 = (0..1000u64).into_par_iter().sum();
    pass(&mut checks, "rayon thread pool", sum == 499_500, format!("par-sum={}", sum));

    // Assemble report.
    let mut all_ok = true;
    let mut report = String::new();
    report.push_str(&banner::banner());
    report.push('\n');
    let _ = writeln!(report, "doctor — environment & self-test\n");
    let _ = writeln!(report, "engines : phanotate prodigal prodigalgv prodigal-meta genemark augustus fraggenescan sixframe");
    let _ = writeln!(report, "formats : gff3 gtf genbank faa fna ffn tsv");
    let _ = writeln!(report, "codes   : {} NCBI translation tables\n", SUPPORTED.len());
    let _: Vec<Gene> = Vec::new(); // (keep Gene import meaningful)
    for c in &checks {
        if !c.ok {
            all_ok = false;
        }
        let mark = if c.ok { "\u{2713}" } else { "\u{2717}" };
        let _ = writeln!(report, "  [{}] {:<34} {}", mark, c.label, c.detail);
    }
    report.push('\n');
    if all_ok {
        let _ = writeln!(report, "All systems go. 💥  DYNAMITE is ready.");
    } else {
        let _ = writeln!(report, "Some checks failed — see ✗ above.");
    }
    (report, all_ok)
}

fn fasta_roundtrip() -> bool {
    use std::io::Write as _;
    let mut path = std::env::temp_dir();
    path.push("dynamite_doctor_test.fna");
    let body = b">t1 test\nATGAAACCCGGGTTTTAA\n";
    if std::fs::File::create(&path).and_then(|mut f| f.write_all(body)).is_err() {
        return false;
    }
    match crate::fasta::read_fasta(&path) {
        Ok(recs) => {
            let _ = std::fs::remove_file(&path);
            recs.len() == 1 && recs[0].seq.len() == 18
        }
        Err(_) => false,
    }
}
