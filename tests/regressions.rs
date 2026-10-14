//! Regression tests for the bugs fixed in 1.1.0. Expected values for the
//! Prodigal engine were produced with pyrodigal 3.7.1 (the reference Prodigal
//! implementation) on the same sequences.

use dynamite::fasta::{read_fasta, Record};
use dynamite::gene::Gene;
use dynamite::model::{CodeChoice, EngineKind, Mode};
use dynamite::{run_records, RunConfig};

fn data(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn lambda() -> Vec<u8> {
    read_fasta(data("lambda.fa")).unwrap().remove(0).seq
}

fn record(id: &str, seq: Vec<u8>) -> Vec<Record> {
    vec![Record {
        id: id.to_string(),
        desc: id.to_string(),
        seq,
    }]
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn run(records: &[Record], mode: Mode, engine: Option<EngineKind>, code: CodeChoice) -> (Vec<Gene>, u32) {
    let cfg = RunConfig {
        mode,
        engine,
        code,
        min_len: 0,
        closed: false,
        threads: 1,
    };
    let out = run_records(records, &cfg).unwrap();
    (out.genes, out.trans_table)
}

fn find(genes: &[Gene], begin: usize, end: usize, strand: i8) -> &Gene {
    genes
        .iter()
        .find(|g| g.begin == begin && g.end == end && g.strand == strand)
        .unwrap_or_else(|| panic!("gene {begin}..{end} ({strand}) not called"))
}

/// Bug 1: `partial=` was written in 5'/3' order for reverse-strand genes.
/// Prodigal's convention is left/right on both strands.
#[test]
fn reverse_strand_partial_flags_follow_prodigal() {
    let l = lambda();
    let rc = revcomp(&l[1001..40001]);
    let (genes, _) = run(&record("cut_rc", rc), Mode::Bacteria, Some(EngineKind::Prodigal), CodeChoice::Fixed(11));
    // pyrodigal: 1..420 (-) partial=10 ATG ; 37366..39000 (-) partial=01 Edge
    let g = find(&genes, 1, 420, -1);
    assert_eq!(g.partial_str(), "10");
    assert_eq!(g.start_codon_str(), "ATG");
    assert!(g.partial_3p() && !g.partial_5p());
    let g = find(&genes, 37366, 39000, -1);
    assert_eq!(g.partial_str(), "01");
    assert_eq!(g.start_codon_str(), "Edge");
    assert!(!g.aa.ends_with('*'));
    // GenBank locations carry the marker on the incomplete genomic end
    assert_eq!(find(&genes, 1, 420, -1).insdc_location(), "complement(<1..420)");
    assert_eq!(find(&genes, 37366, 39000, -1).insdc_location(), "complement(37366..>39000)");
}

/// Bug 2: edge genes padded to start at base 2/3 (or to end at len-1/len-2)
/// were reported complete (and with a real start type).
#[test]
fn frame_padded_edge_genes_are_partial() {
    let l = lambda();
    let cut = l[1003..40003].to_vec(); // 39,000 bp
    let (genes, _) = run(&record("cut", cut), Mode::Bacteria, Some(EngineKind::Prodigal), CodeChoice::Fixed(11));
    // pyrodigal: 2..1633 (+) partial=10 Edge ; 38579..38998 (+) partial=01 ATG
    let g = find(&genes, 2, 1633, 1);
    assert_eq!((g.partial_str().as_str(), g.start_codon_str()), ("10", "Edge"));
    let g = find(&genes, 38579, 38998, 1);
    assert_eq!((g.partial_str().as_str(), g.start_codon_str()), ("01", "ATG"));
}

/// Bug 3/4: N was read as A, so `T` + `NN` became a TAA stop and genes ended
/// inside assembly gaps; partial ends were not inferred and gap codons were
/// not trimmed.
#[test]
fn genes_never_end_in_or_span_assembly_gaps() {
    let l = lambda();
    // In-frame T of the lambda gene at 711..2636 (+), then a 100-N gap.
    let mut gapped = l[..1320].to_vec();
    gapped.extend(std::iter::repeat(b'N').take(100));
    gapped.extend_from_slice(&l[1322..]);
    let gap = 1321..=1420usize;
    for engine in [
        EngineKind::Prodigal,
        EngineKind::ProdigalGv,
        EngineKind::Phanotate,
        EngineKind::GeneMark,
        EngineKind::Augustus,
        EngineKind::FragGeneScan,
    ] {
        let (genes, _) = run(&record("gapped", gapped.clone()), Mode::Auto, Some(engine), CodeChoice::Fixed(11));
        assert!(!genes.is_empty(), "{engine:?}");
        for g in &genes {
            assert!(!gap.contains(&g.begin) && !gap.contains(&g.end), "{engine:?}: {}..{} ends in the gap", g.begin, g.end);
            assert!(!(g.begin < *gap.start() && g.end > *gap.end()), "{engine:?}: {}..{} spans the gap", g.begin, g.end);
            assert!(!g.aa.contains('X') && !g.nt.contains('N'), "{engine:?}: gap codons kept");
            if g.end + 3 >= *gap.start() && g.end < *gap.start() {
                assert!(g.partial_right, "{engine:?}: {}..{} abuts the gap but is not right-partial", g.begin, g.end);
            }
            if g.begin <= *gap.end() + 3 && g.begin > *gap.end() {
                assert!(g.partial_left, "{engine:?}: {}..{} abuts the gap but is not left-partial", g.begin, g.end);
            }
        }
    }
    // The gene that ran into the gap stops at its last full codon, 3'-partial.
    let (genes, _) = run(&record("gapped", gapped), Mode::Bacteria, Some(EngineKind::Prodigal), CodeChoice::Fixed(11));
    let g = genes.iter().find(|g| g.strand == 1 && g.end == 1319).expect("gene truncated at the gap");
    assert_eq!(g.partial_str(), "01");
}

/// Partial ends inferred from codons for every engine: a complete 5' end
/// starts with a start codon (protein begins with M); a complete 3' end ends
/// with a stop codon.
#[test]
fn partial_flags_agree_with_codons_for_all_engines() {
    let l = lambda();
    let reads: Vec<Record> = (0..40)
        .map(|k| Record {
            id: format!("r{k}"),
            desc: String::new(),
            seq: l[k * 1200..k * 1200 + 300].to_vec(),
        })
        .collect();
    let code = dynamite::gencode::GeneticCode::new(11).unwrap();
    let starts = [b"ATG", b"GTG", b"TTG"];
    let stops = [b"TAA", b"TAG", b"TGA"];
    let mut cases: Vec<(Vec<Record>, Mode, EngineKind)> = Vec::new();
    for e in [EngineKind::Prodigal, EngineKind::Phanotate, EngineKind::GeneMark, EngineKind::Augustus] {
        cases.push((record("cut", l[1003..40003].to_vec()), Mode::Auto, e));
    }
    cases.push((reads, Mode::Reads, EngineKind::FragGeneScan));
    for (recs, mode, engine) in cases {
        let (genes, _) = run(&recs, mode, Some(engine), CodeChoice::Fixed(11));
        assert!(!genes.is_empty());
        for g in &genes {
            let nt = g.nt.as_bytes();
            if !g.partial_5p() {
                assert!(starts.iter().any(|s| &nt[..3] == *s), "{engine:?} {}..{} complete 5' without start", g.begin, g.end);
                assert!(g.aa.starts_with('M'));
            }
            if !g.partial_3p() {
                assert!(stops.iter().any(|s| &nt[nt.len() - 3..] == *s), "{engine:?} {}..{} complete 3' without stop", g.begin, g.end);
            }
            assert!(!g.aa.contains('*'), "{engine:?} {}..{}: stop in protein", g.begin, g.end);
        }
        let _ = &code;
    }
}

/// Bug 5: the code search returned table 1 for phage lambda (pyrodigal-gv: 11).
#[test]
fn lambda_is_standard_code_in_every_auto_mode() {
    let recs = read_fasta(data("lambda.fa")).unwrap();
    for mode in [Mode::Virus, Mode::GiantVirus, Mode::Crassphage, Mode::Meta, Mode::Auto] {
        let (_, tt) = run(&recs, mode, None, CodeChoice::Default);
        assert_eq!(tt, 11, "{mode:?}");
    }
    let (_, tt) = run(&recs, Mode::Auto, Some(EngineKind::ProdigalGv), CodeChoice::Auto);
    assert_eq!(tt, 11);
}

/// ...while a genuine TAG→Gln phage is still detected as table 15.
#[test]
fn topaz_is_detected_as_code_15() {
    let recs = read_fasta(data("topaz.fa")).unwrap();
    for mode in [Mode::Virus, Mode::GiantVirus, Mode::Crassphage, Mode::Meta] {
        let (genes, tt) = run(&recs, mode, None, CodeChoice::Default);
        assert_eq!(tt, 15, "{mode:?}");
        assert!(genes.iter().all(|g| g.trans_table == 15));
    }
}

/// An explicit `-g auto` was silently ignored by single-table presets (e.g.
/// `-e prodigal -g auto` always used 11). It now searches with every engine
/// that can, and keeps the engine's own table when no alternative stop set
/// is found (AUGUSTUS: 1, not 11).
#[test]
fn explicit_auto_searches_for_every_engine() {
    let topaz = read_fasta(data("topaz.fa")).unwrap();
    for engine in [EngineKind::Prodigal, EngineKind::GeneMark, EngineKind::Augustus, EngineKind::Phanotate] {
        let (genes, tt) = run(&topaz, Mode::Bacteria, Some(engine), CodeChoice::Auto);
        assert_eq!(tt, 15, "{engine:?}");
        assert!(genes.iter().all(|g| g.trans_table == 15), "{engine:?}");
    }
    let recs = read_fasta(data("lambda.fa")).unwrap();
    assert_eq!(run(&recs, Mode::Bacteria, Some(EngineKind::Prodigal), CodeChoice::Auto).1, 11);
    assert_eq!(run(&recs, Mode::Auto, Some(EngineKind::Augustus), CodeChoice::Auto).1, 1);
}

/// A leading or trailing assembly gap must not wipe a contig's PHANOTATE
/// calls (the path's entry / exit windows are measured from the outermost
/// gap-free segments), and calls away from the gap are unchanged.
#[test]
fn phanotate_survives_terminal_gaps() {
    let l = lambda();
    let (plain, _) = run(&record("lam", l.clone()), Mode::Phage, None, CodeChoice::Default);
    assert!(plain.len() > 50);
    let mut lead = vec![b'N'; 2100];
    lead.extend_from_slice(&l);
    let mut trail = l.clone();
    trail.extend(std::iter::repeat(b'N').take(3000));
    for (name, seq, shift) in [("lead", lead, 2100usize), ("trail", trail, 0usize)] {
        for closed in [false, true] {
            let cfg = RunConfig { mode: Mode::Phage, engine: None, code: CodeChoice::Default, min_len: 0, closed, threads: 1 };
            let genes = run_records(&record(name, seq.clone()), &cfg).unwrap().genes;
            assert!(!genes.is_empty(), "{name} closed={closed}: no genes");
            if !closed {
                let same = genes
                    .iter()
                    .filter(|g| plain.iter().any(|p| p.begin + shift == g.begin && p.end + shift == g.end && p.strand == g.strand))
                    .count();
                assert!(same * 10 >= plain.len() * 9, "{name}: only {same}/{} genes unchanged", plain.len());
            }
        }
    }
}

/// PHANOTATE reads U as T, like every other engine.
#[test]
fn phanotate_reads_u_as_t() {
    let l = lambda();
    let u: Vec<u8> = l.iter().map(|&b| if b == b'T' { b'U' } else { b }).collect();
    let (dna, _) = run(&record("lam", l), Mode::Phage, None, CodeChoice::Default);
    let (rna, _) = run(&record("lam", u), Mode::Phage, None, CodeChoice::Default);
    let key = |g: &Gene| (g.begin, g.end, g.strand, g.aa.clone());
    assert_eq!(dna.iter().map(key).collect::<Vec<_>>(), rna.iter().map(key).collect::<Vec<_>>());
}

/// Records that share an id are called separately and kept apart, in input
/// order, with ids numbered per record (as 1.0.0 did).
#[test]
fn duplicate_record_ids_stay_apart() {
    let l = lambda();
    let recs = vec![
        Record { id: "a".into(), desc: String::new(), seq: l[..12000].to_vec() },
        Record { id: "b".into(), desc: String::new(), seq: l[12000..24000].to_vec() },
        Record { id: "a".into(), desc: String::new(), seq: l[30000..42000].to_vec() },
    ];
    for engine in [EngineKind::Prodigal, EngineKind::Phanotate, EngineKind::GeneMark] {
        let (genes, _) = run(&recs, Mode::Auto, Some(engine), CodeChoice::Fixed(11));
        // contig order a, b, a; positions ascend and ids restart within each block
        let blocks: Vec<&str> = genes.iter().map(|g| g.contig.as_str()).fold(Vec::new(), |mut v, c| {
            if v.last() != Some(&c) {
                v.push(c);
            }
            v
        });
        assert_eq!(blocks, ["a", "b", "a"], "{engine:?}");
        for w in genes.windows(2) {
            if w[0].contig == w[1].contig && w[1].id != 1 {
                assert_eq!(w[1].id, w[0].id + 1, "{engine:?}");
                assert!(w[0].begin <= w[1].begin, "{engine:?}");
            }
        }
        // PHANOTATE models each record on its own: the first block is exactly
        // the first record's calls (the other engines train on all records)
        if engine == EngineKind::Phanotate {
            let first_a: Vec<&Gene> = genes.iter().take_while(|g| g.contig == "a").collect();
            let single = run(&recs[..1], Mode::Auto, Some(engine), CodeChoice::Fixed(11)).0;
            assert_eq!(first_a.len(), single.len());
            assert!(first_a.iter().zip(&single).all(|(a, b)| a.begin == b.begin && a.aa == b.aa));
        }
    }
}

/// Prodigal's gap options for the Prodigal engines (expected calls from
/// pyrodigal 3.7.1, mask=False / mask=True, on λ with a 100-N gap at 1321).
#[test]
fn gap_modes_follow_prodigal() {
    use dynamite::{run_records_with, GapMode, RunOptions};
    let l = lambda();
    let mut gapped = l[..1320].to_vec();
    gapped.extend(std::iter::repeat(b'N').take(100));
    gapped.extend_from_slice(&l[1322..]);
    let recs = record("gapped", gapped);
    let cfg = RunConfig { mode: Mode::Bacteria, engine: Some(EngineKind::Prodigal), code: CodeChoice::Default, min_len: 0, closed: false, threads: 1 };
    let call = |gm: GapMode| run_records_with(&recs, &cfg, &RunOptions::default().with_gap_mode(gm)).unwrap().genes;
    let first = |genes: &[Gene]| -> Vec<(usize, usize, i8, String)> {
        genes.iter().filter(|g| g.begin < 2800).map(|g| (g.begin, g.end, g.strand, g.partial_str())).collect()
    };
    let t = |b: usize, e: usize, s: i8, p: &str| (b, e, s, p.to_string());

    // 0 (default): partial genes run into the gap from both sides
    let g0 = call(GapMode::Partial);
    assert_eq!(first(&g0), vec![t(191, 736, 1, "00"), t(711, 1319, 1, "01"), t(1421, 2734, 1, "10"), t(2731, 2937, 1, "00")]);
    assert_eq!(g0.len(), run(&recs, Mode::Bacteria, Some(EngineKind::Prodigal), CodeChoice::Default).0.len());
    // 1: genes cannot run into the gap
    let g1 = call(GapMode::Closed);
    assert_eq!(first(&g1), vec![t(191, 736, 1, "00"), t(771, 998, -1, "00"), t(1457, 2734, 1, "00"), t(2731, 2937, 1, "00")]);
    // 2: N's are not gaps — Prodigal 2.x / pyrodigal default: 711..1460 spans it
    let g2 = call(GapMode::NoGaps);
    assert_eq!(first(&g2), vec![t(191, 736, 1, "00"), t(711, 1460, 1, "00"), t(1457, 2734, 1, "00"), t(2731, 2937, 1, "00")]);
    assert_eq!(g2.len(), 61);
    assert_eq!(g2.iter().find(|g| g.begin == 711).unwrap().aa.matches('X').count(), 34);
    // -m: pyrodigal mask=True — nothing is built across the masked run
    let gm = call(GapMode::Mask);
    assert_eq!(first(&gm), vec![t(191, 736, 1, "00"), t(771, 998, -1, "00"), t(1457, 2734, 1, "00"), t(2731, 2937, 1, "00")]);
    assert_eq!(gm.len(), 61);
    // other engines ignore the option (they always split at gaps)
    let cfg_gm = RunConfig { engine: Some(EngineKind::GeneMark), ..cfg.clone() };
    let a = run_records_with(&recs, &cfg_gm, &RunOptions::default().with_gap_mode(GapMode::NoGaps)).unwrap().genes;
    let b = run_records(&recs, &cfg_gm).unwrap().genes;
    assert_eq!(a.iter().map(|g| (g.begin, g.end)).collect::<Vec<_>>(), b.iter().map(|g| (g.begin, g.end)).collect::<Vec<_>>());
}

/// The GenBank LOCUS line must separate a long sequence id from the length
/// (1.0.0 wrote `..._cov_34761361 bp` for a 761,361 bp contig).
#[test]
fn genbank_locus_line_separates_long_ids() {
    let id = "KLEBSIELLA_00_KL11214_0_1_length_761361_cov_34";
    let l = lambda();
    let recs = record(id, l.clone());
    let (genes, _) = run(&recs, Mode::Bacteria, Some(EngineKind::Prodigal), CodeChoice::Default);
    let mut out = String::new();
    dynamite::gene::write_genbank_full(&mut out, &genes, &[(id.to_string(), l.clone())]);
    let locus = out.lines().next().unwrap();
    let tok: Vec<&str> = locus.split_whitespace().collect();
    assert_eq!(&tok[..4], &["LOCUS", id, &l.len().to_string(), "bp"], "{locus}");
    for id16 in ["exactly16chars__", "short"] {
        let tok: Vec<String> = dynamite::gene::locus_line(id16, 123).split_whitespace().map(String::from).collect();
        assert_eq!(tok[1..3], [id16.to_string(), "123".to_string()]);
    }
}

/// AUGUSTUS trained on contig 0 with every contig's coordinates and aborted
/// (index out of bounds) when a later contig was longer than the first.
#[test]
fn augustus_multi_contig_does_not_panic() {
    let l = lambda();
    let recs = vec![
        Record { id: "short".into(), desc: String::new(), seq: l[..5000].to_vec() },
        Record { id: "long".into(), desc: String::new(), seq: l.clone() },
    ];
    let (genes, tt) = run(&recs, Mode::Eukaryote, None, CodeChoice::Default);
    assert_eq!(tt, 1);
    assert!(genes.iter().any(|g| g.contig == "long"));
}

/// FragGeneScan aborted on 1–2 bp records (and gap splitting can create tiny
/// pieces); every engine must survive degenerate input.
#[test]
fn degenerate_records_do_not_panic() {
    let l = lambda();
    let mut tiny_pieces = l[..3000].to_vec();
    tiny_pieces.extend(std::iter::repeat(b'N').take(20));
    tiny_pieces.extend_from_slice(b"AC");
    tiny_pieces.extend(std::iter::repeat(b'N').take(20));
    tiny_pieces.extend_from_slice(&l[3000..6000]);
    let recs = vec![
        Record { id: "one".into(), desc: String::new(), seq: b"A".to_vec() },
        Record { id: "two".into(), desc: String::new(), seq: b"AC".to_vec() },
        Record { id: "allN".into(), desc: String::new(), seq: vec![b'N'; 500] },
        Record { id: "pieces".into(), desc: String::new(), seq: tiny_pieces },
    ];
    for engine in [
        EngineKind::Prodigal,
        EngineKind::ProdigalGv,
        EngineKind::ProdigalMeta,
        EngineKind::Phanotate,
        EngineKind::GeneMark,
        EngineKind::Augustus,
        EngineKind::FragGeneScan,
    ] {
        for mode in [Mode::Auto, Mode::Reads] {
            let _ = run(&recs, mode, Some(engine), CodeChoice::Default);
        }
    }
}

/// An unsupported table used to panic inside the Prodigal engine.
#[test]
fn unsupported_code_is_an_error_not_a_panic() {
    let cfg = RunConfig {
        mode: Mode::Bacteria,
        engine: Some(EngineKind::Prodigal),
        code: CodeChoice::Fixed(7),
        min_len: 0,
        closed: false,
        threads: 1,
    };
    let err = run_records(&record("x", lambda()), &cfg).err().expect("error");
    assert!(err.to_string().contains("unsupported genetic code 7"));
}
