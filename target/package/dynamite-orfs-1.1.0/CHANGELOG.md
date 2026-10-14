# Changelog

## 1.1.0 — Prodigal gap options + bug fixes

Drafted as 1.0.2; the new gap options make it a minor release under semver. Source audited against
GitHub `main` and the crates.io `dynamite-orfs` 1.0.0 / 1.0.1 packages (the `src/` trees of all three
are byte-identical, so every fix below applies to both). Reference behaviour was measured with
pyrodigal 3.7.1, pyrodigal-gv 0.3.2 and Biopython.

### Added (requested)

- [x] **Prodigal's gap handling for the Prodigal engines** (`prodigal`, `prodigalgv`,
  `prodigal-meta`): `--gap-mode 0|1|2`, the values of Prodigal's `-e`. DYNAMITE's `-e` already picks
  the engine, so it also takes Prodigal's spelling: a number is the gap mode (`-e 1`, or
  `-e prodigal -e 1`). Other engines always use mode 0 (and warn if given a gap option).
  - `0` (default, the 1.0.2-draft behaviour): runs of ≥ 10 N are gaps; partial genes run into them as
    into a contig end.
  - `1`: genes cannot run into gaps — each gap-free segment is called with its gap ends closed.
  - `2`: N's are not gaps (Prodigal 2.6.3 without `-m`): genes may span runs of N (`X` in the protein).
- [x] **`-m` / `--mask`** — Prodigal 2.x `-m`: runs of ≥ 50 N are masked and no gene is built across
  or into them, in training and prediction (pyrodigal `mask=True`; trailing runs included, which
  Prodigal 2.6.3 itself missed). `-m` and a gap mode are alternatives.
- [x] Library: `run_records_with` / `annotate_with` take `RunOptions` (`GapMode`); `run_records` /
  `annotate` and `RunConfig` are unchanged. `annotate`'s `.log` records the genetic code and gap mode.
- [x] **Checked:** `--gap-mode 2` and `-m` are identical to pyrodigal (`mask=False` / `True`, open and
  `--closed`) on 10 test sets — 40 runs, 25,608 genes: coordinates, `partial=`, start type, RBS,
  protein. *Kutzneria* KK037166 differs in every mode because Prodigal classes it as
  non-Shine–Dalgarno (known caveat). `--gap-mode 1`: no gene runs into a gap (all three Prodigal
  engines, open and `--closed`, 9 inputs, Biopython audit clean). Gap mode 0 output is byte-identical
  to the 1.0.2 draft on single-record inputs. On λ with one 100-N gap: `711..1460` spans the gap
  under `2` (as in Prodigal 2.6.3); `-m` replaces it with `771..998` (as pyrodigal `mask=True`), and
  so does `1`; `0` lets genes run into the gap from both sides (`711..1319`, `1421..2734`).

### Reported bugs

- [x] **`partial=` written in 5′/3′ order for reverse-strand genes.** `Gene::partial_str` swapped the
  flags on the minus strand, and the engines disagreed on what the flags meant (Prodigal / PHANOTATE:
  genomic left/right; the shared ORF scanner used by GeneMark / AUGUSTUS / six-frame: 5′/3′). One
  convention now everywhere — Prodigal's genomic `partial=XY` (left, right) — with `partial_5p()` /
  `partial_3p()` for strand-relative uses. GFF3 / GTF / FASTA headers / TSV, both GenBank writers
  (`<` / `>`), the annotate `.tbl` (5′/3′ columns), `.gbk` and `.err`, and the M-start / stop-strip
  translation rules were all corrected.
- [x] **Edge genes padded to start at base 2/3 reported complete.** The Prodigal engine required
  `begin == 1` / `end == len` in addition to the edge flag, so two of every three frames' edge genes
  were `partial=00`; edge starts were also labelled `ATG` instead of `Edge`.
- [x] **N read as A, so genes could end inside assembly gaps.** `T`+`NN` became a `TAA` stop, `NTG`
  an `ATG` start, `NNN` translated as `K` and `.ffn` printed `A` for `N`. Ambiguous bases are now
  flagged (stored as Prodigal's placeholder `C` on both strands) and never form start/stop codons in
  any engine; they translate to `X` and are written as `N`. Runs of ≥ 10 N are assembly gaps: contigs
  are split there for prediction, so no gene spans, starts in or ends in a gap. Self-training still
  uses the contigs exactly as given (as Prodigal does), so genes away from a gap are the ones
  Prodigal calls on the same sequence.
- [x] **Partial ends not inferred from start/stop codons; gap codons not trimmed.** New `finalize`
  step for every engine: terminal codons touching an ambiguous base are trimmed in frame, a gene is
  5′-partial unless it begins with a start codon (or starts in the first codon of an open edge, as in
  Prodigal) and 3′-partial unless it ends with a stop codon; `start_type=Edge` for 5′-partial genes.
  FragGeneScan previously reported every gene `partial=00`.
- [x] **Code search returned table 1 for phage λ (pyrodigal-gv: 11).** Summed coding scores are not
  comparable across tables (even real Prodigal scores λ 28% higher under table 1). Replaced by a
  coding-density test (CheckM / GTDB-Tk margin) applied hierarchically over stop-codon sets: an
  alternative table must reach > 70% density and beat the standard code and every table with a
  superset of its stops by > 5 points, the reassigned codons must be read through at a sense-codon-like
  rate (5–50 per 1,000 codons), at least 20% of its coding sequence must be in genes the standard code
  also finds in the same frame, and the input must be ≥ 20 kb. Tables that cannot change the calls
  (1 vs 11, 4 vs 25) are collapsed. The per-mode candidate lists are now honoured (they were ignored).
  λ → 11 in every mode; four *Mycoplasma mycoides* → 4; phages Agate / Topaz → 15; *E. coli* and
  14 further phages → 11 (all as pyrodigal-gv); random sequence and reversed λ → 11.

### Other bugs found

- [x] **Prodigal port deviated from Prodigal.** `dprog`'s giant-ORF window compared the wrong node
  (the DP considered every earlier node), Prodigal's `ovlp` carry in `score_connection` and the `-1`
  GC-frame sentinel were not reproduced. λ matched Prodigal on only 56/61 genes; the trained model is
  now bit-identical and calls match Prodigal on all 8,248 genes tested.
- [x] Prodigal `rbs_motif` / `rbs_spacer` were always `None`; now reported as Prodigal does (SD).
- [x] ORF scanner (GeneMark / AUGUSTUS / six-frame) emitted interior start-less stop-to-stop regions
  as partial genes, and ignored `require_start` in the trailing region. (`translate --orfs` keeps its
  documented stop-to-stop listing; `--require-start` now applies to every region.)
- [x] AUGUSTUS trained on contig 0 with every contig's coordinates → wrong model and an
  index-out-of-bounds **abort** when a later contig was longer; GeneMark mixed coordinates across
  contigs. Training is now pooled per contig.
- [x] An unsupported table (`-g 7`) **panicked** inside the Prodigal engines; now a clean error.
- [x] `-g auto` was silently ignored by every single-table preset (`-e prodigal`, `-p bacteria`,
  GeneMark, AUGUSTUS, PHANOTATE all used a fixed table) — it now searches, keeping the engine's own
  table unless the stop codons change (AUGUSTUS stays on 1). FragGeneScan's HMM is code-11 only:
  `-g` used to be dropped silently, now with a warning.
- [x] GC fraction counted S/B/V as GC in the numerator only (could exceed 1.0); multi-contig training
  decoded N as a concrete base.
- [x] Greedy overlap resolution re-ordered contigs alphabetically (`contig10` before `contig2`) and was
  quadratic.
- [x] PHANOTATE ignored `--closed`; FragGeneScan ignored `--closed` and `--threads`.
- [x] PHANOTATE read `U` as an unknown base (an RNA-alphabet FASTA of λ gave a single "gene" instead
  of 68); every other engine already read `U` as `T`.
- [x] Records sharing an id: GeneMark / AUGUSTUS resolved overlaps *across* them (a gene of one record
  could be dropped for overlapping another record's coordinates) and numbered them as one contig.
  Every engine now calls, orders and numbers them separately, and `call` / `all` / `annotate` warn,
  because every output format identifies sequences by id (GenBank / `annotate` still list a shared
  id's genes under each same-named record).
- [x] FragGeneScan was documented as a "frameshift-aware decode" for error-prone short and long
  reads, but DYNAMITE always loads FragGeneScan's error-free `complete` model (insertion / deletion
  transitions 1e-18): 0 of 4,580 genes were frameshift-corrected on λ reads with 0.5% insertions +
  0.5% deletions, and `reads` / `long-reads` give identical calls. Help text, README and module docs
  now say so (switching models is listed below).
- [x] FragGeneScan **aborted** on records shorter than 3 bp (vendored decoder indexes the first
  codon); such records cannot hold a gene and are now skipped.
- [x] FragGeneScan (vendored, upstream behaviour): minus-strand start refinement moves the coordinate
  but drops the codons in between from its own sequence strings — proteins disagreed with
  coordinates. Frameshift-free genes now take their sequences from the genome at FragGeneScan's
  (unchanged) coordinates.
- [x] **GenBank: a contig id of 16+ characters ran into the length on the `LOCUS` line**
  (`LOCUS KLEBSIELLA_00_KL11214_0_1_length_761361_cov_34761361 bp`, reported for `all -p meta`); now
  always whitespace-separated, in `call`, `all`, `annotate` and `translate` output (Biopython reads
  the id and length back). `call -f genbank` / `translate --orfs -f genbank` also wrote records
  without `ORIGIN` — both rejected by GenBank parsers.
- [x] Multi-contig training omitted Prodigal's last linker: Prodigal and pyrodigal join contigs as
  `c1 + TTAATTAATTAA + c2 + TTAATTAATTAA`; DYNAMITE dropped the final one, so training on many short
  records could differ slightly (e.g. 300 λ reads, 333 → 335 genes).
- [x] `annotate`: tRNA product was the raw scanner note (`anticodon=…;aa=…`) instead of `tRNA-Xxx`;
  anticodons were read with table 11 even when another code was detected (the *Mycoplasma*
  TGA-reading tRNA-Trp was dropped); `.gbk` CDS had no partial markers.
- [x] `rna` repeated the `##gff-version 3` header for every contig.
- [x] FASTA reader rejected files with non-UTF-8 header bytes.
- [x] `--circular` was documented as "currently the default" but does nothing (help text corrected).
- [x] README told users to `cargo install dynamite`, which installs an unrelated crate; the crate is
  `dynamite-orfs`.
- [x] `cargo clippy` aborted on a deny-level lint (`almost_swapped`), hiding the MSRV lint.

### Performance

Same or faster everywhere: *E. coli* K-12 `-p bacteria` 4.5 → 4.4 s, `-p auto` 20.6 → 9.2 s (equivalent
tables are no longer re-run); λ `-p virus` 0.19 → 0.07 s; FragGeneScan on 8,000 × 300 bp reads with
`-t 2` 0.86 → 0.58 s.

### Tests

`cargo test` now runs 49 unit tests and 15 regression tests (`tests/regressions.rs`, using λ and the
table-15 phage Topaz in `tests/data/`). 12 of the 14 regression tests that run on 1.0.0 / 1.0.1 fail
there, including the long-id `LOCUS` line; the gap-mode test needs the new API. The two that pass
guard against regressions: the old search did detect genuine code-15 genomes (Topaz), and PHANOTATE
must keep its calls on contigs with a leading or trailing gap (a fault in an unreleased draft of
this fix).
