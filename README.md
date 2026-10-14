<div align="center">

```
                                                      _.-^^~~~~~^^-._
                                                  _.-~ ░░▒▒▓▓▓▓▒▒░░ ~-._
                                                .'  ░▒▓████████████▓▒░  '.
                                               ( ░▒▓██████████████████▓▒░ )
                                                '._ ░▒▓████████████▓▒░ _.'
                                                   '~-.░▒▓██████▓▒░.-~'
 ______   ___   _    _    __  __ ___ _____ _____      \ ░▒▓███▓▒░ /
|  _ \ \ / / \ | |  / \  |  \/  |_ _|_   _| ____|      \ ▒▓███▓▒ /
| | | \ V /|  \| | / _ \ | |\/| || |  | | |  _|   ╭──╮  \░▒▓█▓▒░/
| |_| || | | |\  |/ ___ \| |  | || |  | | | |___  ││▒▓██╼╾╼ ✸
|____/ |_| |_| \_/_/   \_\_|  |_|___| |_| |_____| ╰──╯
```

### an ORF / gene caller for **all domains of life**

[![Rust](https://img.shields.io/badge/Rust-1.85%2B-orange?logo=rust)](https://www.rust-lang.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.1.0-success)](https://github.com/raw-lab/dynamite)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS%20%7C%20windows-lightgrey)](https://github.com/raw-lab/dynamite)
[![Single binary](https://img.shields.io/badge/single%20binary-~8%20MB-brightgreen)](https://github.com/raw-lab/dynamite)

*One pure-Rust binary that finds genes in phages, viruses, giant viruses, crassphages,*
*bacteria, archaea, eukaryotes — and raw sequencing reads.*

</div>

---

## 💡 What is DYNAMITE?

DYNAMITE detonates a genome into its genes. Point it at a FASTA, tell it *what kind of organism*
it is (or let it guess), and it picks the right algorithm and the right genetic code, then writes
genes in the format you want. No Python, no BLAST, no shell-outs — just one small binary.

It bundles seven gene-calling engines — six pure-Rust re-implementations plus the **vendored real
FragGeneScanRs HMM** — behind **one command and one output schema**, plus helpers for translation,
codon bias, and RNA genes.

## ⚡ Quick start

```sh
# install from crates.io (the crate is `dynamite-orfs`; it installs the `dynamite` binary —
# note `cargo install dynamite` would fetch an unrelated crate)
cargo install dynamite-orfs --locked

# …or from a local checkout (use --locked for the known-good, reproducible dependency set)
cargo install --path . --locked

# 1) a bacterial genome → GFF3
dynamite call -i genome.fna -p bacteria -o genes.gff

# 2) a phage → GenBank (PHANOTATE engine)
dynamite call -i phage.fna -p phage -f genbank

# 3) a giant virus, auto-detecting the genetic code → proteins
dynamite call -i ncldv.fna -p giant-virus -g auto -f faa

# kick the tires anytime:
dynamite doctor
```

## 🧭 Which mode should I use?

Pick a **mode** with `-p`; DYNAMITE chooses the engine and a sensible genetic code for you.

| Your input | `-p mode` | Engine used | Code |
|---|---|---|---|
| Phage genome | `phage` | PHANOTATE | 11 |
| Bacterium / archaeon | `bacteria` / `archaea` | Prodigal | 11 |
| Giant virus (NCLDV) | `giant-virus` | Prodigal-GV | auto |
| Crassphage | `crassphage` | Prodigal-GV | auto |
| Any virus | `virus` | Prodigal-GV | auto |
| Metagenome contigs | `meta` | Prodigal-meta | auto |
| Eukaryote | `eukaryote` | AUGUSTUS | 1 |
| Short reads (Illumina) | `reads` | FragGeneScanRs* | 11 |
| Long reads (Nanopore/PacBio) | `long-reads` | FragGeneScanRs* | 11 |
| Not sure | `auto` | Prodigal-meta | auto |

<sub>*the `reads` and `long-reads` engines run the **vendored, verbatim FragGeneScanRs HMM** — see below.</sub>

Want a specific engine instead of a mode? Use `-e/--engine`. Force a genetic code with
`-g/--code <N|auto>`.

## 🧬 The engines

| `-e engine` | Modelled on | Sweet spot |
|---|---|---|
| `phanotate` | PHANOTATE (graph / shortest-path) | dense, overlapping phage genomes |
| `prodigal` | Prodigal (self-training DP) | single bacterial / archaeal genomes |
| `prodigalgv` | Prodigal-GV | phages, giant viruses, crassphages (alt codes) |
| `prodigal-meta` | Prodigal `-p meta` | fragmented / mixed metagenomes |
| `genemark` | GeneMark (3-periodic Markov) | a fast Markov second opinion |
| `augustus` | AUGUSTUS (eukaryotic) | single-exon eukaryotic CDS |
| `fraggenescan` | **FragGeneScanRs (vendored HMM)** | short & long reads, gene fragments |

The three `prodigal*` presets share one validated Prodigal dynamic-programming core; `genemark`
reuses its 3-periodic Markov coding model. **`fraggenescan` is different in kind:** rather than a
re-implementation, DYNAMITE *vendors* the real [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs)
HMM (forward/Viterbi decode, original trained models embedded via
`include_bytes!`) — so its output is **bit-for-bit identical** to FragGeneScanRs. Every engine is
pure Rust; DYNAMITE never shells out to an external binary. Run `dynamite call --help` for the full
engine→module map.

## 🛠️ Commands

```sh
dynamite call       -i in.fna -p bacteria -f gff3 -o out.gff   # gene calling (one format)
dynamite all        -i in.fna -p bacteria -o outdir            # gene calling → every format
dynamite annotate   -i in.fna -p bacteria -o outdir            # full Prokka-style bundle
dynamite translate  -i in.fna -g 11                            # six-frame protein FASTA
dynamite translate  -i in.fna --orfs -f tsv                    # …or six-frame ORF table
dynamite codon      -i in.fna -p meta -o codons.tsv            # codon usage / RSCU / ENC / density
dynamite codon      -i cds.fna --cds                           # …when input is already CDS
dynamite rna        -i in.fna -o ncrna.gff                     # heuristic tRNA / rRNA
dynamite doctor                                                # environment + self-test
```

**Output formats** (`-f`): `gff3` · `gtf` · `genbank` · `faa` (proteins) · `fna` (gene nucleotides;
alias `ffn`) · `tsv`. Input is nucleotide FASTA — **contigs, short reads, or long reads** — plain or
gzip-compressed (detected by content). `-t/--threads` parallelises by contig. Give every record a
unique id: outputs identify sequences by id (DYNAMITE warns about duplicates).

**Partial genes & assembly gaps.** Every engine reports partial genes the Prodigal way:
`partial=XY` where `X` is the **left** (lower-coordinate) end and `Y` the **right** end, on *both*
strands (`1` = incomplete), with `<`/`>` on the same ends in GenBank and 5′/3′ columns in the NCBI
`.tbl`. A gene is 5′-partial unless it begins with a start codon of the active genetic code (a start
in the first codon of an open contig end counts as partial, as in Prodigal) and 3′-partial unless it
ends with a stop codon; 5′-partial genes get `start_type=Edge`. N / IUPAC bases never form start or
stop codons and translate to `X`. Runs of ≥ 10 N are treated as **assembly gaps**: contigs are split
there, so no gene spans, starts in or ends in a gap — genes that run into one are truncated at their
last full codon and flagged partial on that side. Self-training still sees the contigs exactly as
given (as Prodigal does), so genes away from a gap are the ones Prodigal calls on the same sequence.

**Gap handling for the Prodigal engines** (`prodigal`, `prodigalgv`, `prodigal-meta`) follows
Prodigal's own options. The other engines always use mode `0`.

| Option | Prodigal | What happens at a run of N |
|---|---|---|
| `--gap-mode 0` *(default)* | `-e 0` | ≥ 10 N is a gap; partial genes run into it as into a contig end (as if the contig were cut there) |
| `--gap-mode 1` | `-e 1` | ≥ 10 N is a gap; genes cannot run into it — every gene next to a gap has its own start and stop codon |
| `--gap-mode 2` | `-e 2`, or 2.x without `-m` | N's are not gaps: genes may span them (`X` in the protein) — identical to Prodigal 2.6.3 / pyrodigal |
| `-m`, `--mask` | 2.x `-m` | ≥ 50 N is masked: no gene is built across or into it, in training and prediction — identical to pyrodigal `mask=True` |

DYNAMITE's `-e` already picks the engine, so it also accepts Prodigal's spelling: a number given to
`-e` is the gap mode (`-e 1`, or `-e prodigal -e 1` for both).

## 🗂️ Every format at once — the `all` command

Use **`all`** when you want the raw gene-calling output in **every format from a single run** — no
annotation layer (no `locus_tag`s, no `product=`, no tRNA/rRNA, no NCBI submission files). For the
richer Prokka-style bundle, use [`annotate`](#-whole-genome-annotation-annotate) instead.

```sh
# minimal — auto-detect the organism/code, write to ./dynamite_out/dynamite.*
dynamite all -i genome.fna

# typical — bacterial genome, custom directory + filename prefix
dynamite all -i genome.fna -p bacteria -o results/ --prefix mygenome

# force an engine + genetic code, complete genes only, 8 threads
dynamite all -i contigs.fna -e prodigalgv -g 11 --closed -t 8 -o out/ --prefix sample1
```

**The seven files written** (all share `--prefix`, default `dynamite`):

| File | Contents |
|---|---|
| `{prefix}.gff` | gene calls in **GFF3** (DYNAMITE attributes) |
| `{prefix}.gtf` | gene calls in **GTF** |
| `{prefix}.gbk` | **GenBank** — CDS features **and** the `ORIGIN` sequence |
| `{prefix}.faa` | translated **CDS proteins** |
| `{prefix}.ffn` | **CDS gene nucleotides** |
| `{prefix}.fna` | the **input contig** sequences |
| `{prefix}.tsv` | coordinate / feature table |

**All `all` options** (also shown by `dynamite all --help`):

| Option | Default | Meaning |
|---|---|---|
| `-i, --input <FASTA>` | – | Input nucleotide FASTA — contigs, short reads, or long reads (`.gz` ok). A positional `FASTA` works too. |
| `-o, --outdir <DIR>` | `dynamite_out` | Directory to write the files into (created if missing). |
| `--prefix <NAME>` | `dynamite` | Base name shared by every output file (`{prefix}.gff`, `{prefix}.faa`, …). |
| `-p, --mode <MODE>` | `auto` | `phage`·`virus`·`giant-virus`·`crassphage`·`bacteria`·`archaea`·`meta`·`eukaryote`·`reads`·`long-reads`·`auto` |
| `-e, --engine <ENGINE>` | from `--mode` | Force an engine: `phanotate`·`prodigal`·`prodigalgv`·`prodigal-meta`·`genemark`·`augustus`·`fraggenescan` (a number `0`/`1`/`2` sets the gap mode instead, as Prodigal's `-e`) |
| `-g, --code <CODE\|auto>` | from `--mode` | NCBI translation-table number (e.g. `11`) or `auto` to detect. |
| `--gap-mode <0\|1\|2>` | `0` | Prodigal engines: how runs of N are handled (see *Gap handling* above). |
| `-m, --mask` | off | Prodigal engines: Prodigal 2.x `-m` — no gene across or into ≥ 50 N. |
| `--closed` | off | Require complete genes only (no genes running off contig ends). |
| `-l, --minlen <NT>` | `0` | Minimum ORF length in nt (PHANOTATE engine only). |
| `-t, --threads <N>` | `0` | Worker threads; `0` = all CPU cores. |

## 🧬 Whole-genome annotation (`annotate`)

One command runs the chosen engine **plus** heuristic tRNA/rRNA detection, assigns sequential
`locus_tag`s along each contig, and writes the full Prokka-style file set in one pass — for any engine
(`prodigal`, `prodigalgv`, …):

```sh
dynamite annotate -i genome.fna -p bacteria \
  --prefix mygenome --locus MYG --organism "Escherichia coli" --strain K12 -o out/
```

| File | Contents |
|---|---|
| `.gff` | Master **GFF3** (annotations + `##FASTA` block) — opens in Artemis / IGV |
| `.gbk` | **GenBank** with every feature (CDS + tRNA/rRNA), `/translation`, ORIGIN |
| `.fna` | Nucleotide FASTA of the input contigs |
| `.faa` | Protein FASTA of translated CDS |
| `.ffn` | Nucleotide FASTA of **all** transcripts (CDS + tRNA + rRNA) |
| `.tsv` | Feature table: `locus_tag · ftype · length_bp · gene · EC_number · COG · product` |
| `.txt` | Feature-count statistics |
| `.tbl` | NCBI **Feature Table** (input to `table2asn`) |
| `.fsa` | Contigs with Sequin tags (input to `table2asn`) |
| `.err` | DYNAMITE QA report (short/partial/internal-stop checks) |
| `.log` | Run record (engine, settings, counts) |

The `.sqn` (ASN.1 Sequin) and the *official* NCBI discrepancy `.err` are produced by NCBI's
`table2asn`, which is a separate tool DYNAMITE does not bundle (keeping it shell-out-free). DYNAMITE
writes the `.tbl` + `.fsa` it needs and prints the exact command:
`table2asn -i out/mygenome.fsa -f out/mygenome.tbl -o out/mygenome.sqn`. (tmRNA/misc_RNA are not
detected by the heuristic RNA scanner — use **tRNAscan-SE** / **barrnap** / **Aragorn** for production
ncRNA annotation.)

## ✅ Validation

Measured with **pyrodigal 3.7.1** (the reference Prodigal implementation), **pyrodigal-gv 0.3.2** and
Biopython. Every called protein is re-translated from the genome and every partial flag is checked
against the actual start / stop codons — **0 mismatches across all seven engines** on λ (both strands,
edge-cut contigs, 15 inserted N-gaps, scattered IUPAC codes, 300 bp reads with N runs), *Kutzneria*
KK037166 and *E. coli* MIIJ01000039 (real assembly gaps), *E. coli* K-12 with 20 × 100-N gaps, phage
Topaz (code 15) and *Yersinia* phage vB_YenM_201.16. No gene starts in, ends in or spans a gap. Gaps
change nothing away from them: genes more than 1 kb from a gap are identical to Prodigal's on the same
gapped sequence (4,263 / 4,263 on *E. coli*).

- **`prodigal` ≡ Prodigal.** On *E. coli* K-12 MG1655, *C. diphtheriae* NCTC11397, *M. mycoides* GM12
  (table 4), phage Agate (table 15), metagenomic contig SRR492066, *Acinetobacter* phage Apostate,
  *Yersinia* phage vB_YenM_201.16 and λ (plus two edge-cut / reverse-complemented λ contigs), DYNAMITE's
  calls are **identical to Prodigal for all 8,248 genes** — coordinates, strand, `partial=`,
  `start_type`, `rbs_motif`, `rbs_spacer` and protein sequence (the trained model is bit-identical).
  Draft assemblies train exactly as Prodigal does (contigs joined by `TTAATTAATTAA` linkers, one after
  the last too): a 50-contig *E. coli* and a 12-contig *C. diphtheriae* draft give identical calls.
- **Gap modes ≡ Prodigal.** `--gap-mode 2` is identical to Prodigal 2.6.3 (pyrodigal) and `-m` to
  pyrodigal `mask=True` — coordinates, `partial=`, start type, RBS and protein, with and without
  `--closed` — on λ with 7 or 15 gaps and IUPAC runs, a multi-record set with leading / trailing /
  all-N records, λ with 2–5 kb terminal gaps, *E. coli* MIIJ01000039 (one 433 kb N-run), an 8-contig
  *Mycoplasma* (table 4) and *E. coli* K-12 with 20 gaps (4,319 and 4,322 genes). `--gap-mode 1`
  lets no gene run into a gap on any of them. (*Kutzneria* KK037166 differs in every mode because
  Prodigal classes it as non-Shine–Dalgarno; see the caveats.)
- **Genetic-code detection** (`virus` / `giant-virus` / `crassphage` / `meta` / `auto`, and any engine
  with `-g auto`): λ, *E. coli* and 14 further phages → **11**; four *Mycoplasma mycoides* genomes →
  **4**; phages Agate and Topaz → **15** — the same tables pyrodigal-gv reports for all 21 phage /
  *Mycoplasma* genomes. Synthetic TGG→TGA / CAG→TAG recodings of λ and *C. diphtheriae* are detected as
  4 / 15. Random sequence (30–70% GC) and reversed λ → 11: no false alternative code.
- **`prodigalgv` vs Prodigal-GV.** `prodigalgv` self-trains (Prodigal single-genome mode) instead of
  using Prodigal-GV's pre-trained metagenomic models. On λ it reproduces **85%** of pyrodigal-gv's gene
  3′-ends and the **same start 91%** of the time where they agree — exactly the agreement real Prodigal
  single mode has with pyrodigal-gv.
- **`fraggenescan` == FragGeneScanRs.** Because DYNAMITE vendors the real HMM, gene intervals and strands
  are **identical** to FragGeneScanRs on gap-free input (183/183 genes on 161 × 300 bp λ reads).

*(Reproduce: `pip install pyrodigal pyrodigal-gv`, run them and DYNAMITE on the same FASTA, compare
GFF coordinates / `partial=` / proteins. `cargo test` runs the λ / Topaz regression suite.)*

## 🔬 Codon stats & `rustyomestats`

`dynamite codon` computes codon usage natively — counts, fractions, **RSCU**, **ENC** (Wright 1990),
**GC / GC1 / GC2 / GC3**, and **codon density** — with zero extra dependencies. For the full RAW-lab
statistics suite (Castro *U50/UG50*, etc.) enable the optional crate (Rust 1.85+):

```sh
cargo install --path . --locked --features rustyomestats
```

## 📦 Build notes

- **Use `--locked`.** It installs the pinned, tested dependency set. (Plain `cargo install` re-resolves
  and may pull newer transitive crates than your toolchain supports.)
- Release profile: fat LTO, one codegen unit, `panic = "abort"`, stripped. The binary is ~8 MB because
  it **embeds FragGeneScanRs's trained HMM models** (~6.4 MB) via `include_bytes!`, so no external data
  files are ever needed; the non-FGS code is ~1.3 MB.
- Practical MSRV **Rust 1.85+** (driven by transitive crates; the algorithms are edition-2021). The
  pinned `Cargo.lock` (via `--locked`) additionally builds on Rust 1.75.

<details>
<summary>⚠️ <b>Honest scope &amp; caveats</b> (please read before publishing annotations)</summary>

DYNAMITE is faithful **in spirit** to its namesakes; some engines are deliberate simplifications.

- **`phanotate`** reproduces PHANOTATE's phiX174 calls *and* path scores to ~7 significant figures.
- **`prodigal` / `prodigalgv` / `prodigal-meta`** are a Rust port of the Prodigal / Prodigal-GV DP core
  (self-trained coding model, Shine–Dalgarno RBS scoring, start tweaking). For Shine–Dalgarno organisms
  the output is identical to Prodigal (see above). Prodigal's *non-SD* start trainer is not
  implemented: for organisms Prodigal classes as non-SD, start choices can differ and `rbs_motif` is
  reported as `None`. Self-training needs sequence — below 20 kb (Prodigal's own minimum) DYNAMITE warns,
  and genetic-code auto-detection keeps the standard code.
  <br><b>Genetic-code auto-detection</b> self-trains under each candidate table and adopts an
  alternative only on clear evidence: its coding density must exceed 70% and beat the standard code,
  and every candidate that keeps a superset of its stop codons, by > 5 points (the CheckM / GTDB-Tk
  test, generalised); the reassigned codons must be used inside genes like sense codons (5–50 per 1,000
  codons); and at least 20% of its coding sequence must lie in genes the standard code also finds in
  the same frame (random or non-coding sequence otherwise "gains" density under tables with fewer
  stops). Tables that differ only in start codons (1 vs 11) or only in what a reassigned codon encodes
  (4 vs 25) cannot be told apart from the sequence and are never auto-selected over the first-listed one
  (or, with `-g auto`, over the engine's own table) — force them with `-g`.
  <br><b>On pyrodigal / pyrodigal-gv / pyrodigal-rv:</b> these are <i>not</i> a better algorithm —
  pyrodigal is a Cython binding to the <i>same</i> Prodigal C dynamic-programming method (bug fixes +
  a Python API). pyrodigal-<b>gv</b> and pyrodigal-<b>rv</b> only ship extra <i>pre-trained metagenomic
  models</i> for giant viruses (62 models) and RNA viruses (64 models), distributed as 24 MB / 34 MB
  JSON. Embedding those is out of scope (size, and they assume Prodigal's exact Tinf scoring); DYNAMITE's
  `prodigalgv` instead self-trains and searches alternative genetic codes, which is lighter but does not
  reproduce Prodigal-GV's pre-trained anonymous model library.
- **`genemark`** captures GeneMark's statistical core (codon-position Markov chains, self-training,
  coding/non-coding LLR) but is **not** the full GeneMark.hmm (no whole-sequence Viterbi with
  duration/RBS submodels, no GeneMarkS motif refinement). Overlaps are resolved greedily.
- **`augustus`** is a **single-exon** eukaryotic caller (coding Markov model + Kozak start preference).
  It does **not** implement AUGUSTUS's intron/exon GHMM or trained splice models — for production
  eukaryotic annotation use **AUGUSTUS**, **BRAKER**, or **Helixer**.
- **`fraggenescan`** is **not** a simplification — DYNAMITE vendors the actual
  [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs) core (`src/fgs/`,
  GPL-3.0, unmodified algorithm) and runs its full forward/Viterbi HMM with the original embedded
  models. Output is identical to FragGeneScanRs run with `-t complete`. **That is FragGeneScan's
  error-free model** (insertion / deletion transitions of 1e-18), so sequencing-error frameshifts are
  *not* corrected — on λ reads with 0.5% insertions + 0.5% deletions, 0 of 4,580 genes were
  frameshift-corrected. FragGeneScan's error models (`illumina_1/5/10`, `454_5/10/30`, `sanger_5/10`)
  are embedded but not selectable yet. `reads` and `long-reads` both use the short-read/fragment
  decode (identical calls); forcing `-e fraggenescan` on assembled contigs uses the complete-genome
  decode. DYNAMITE adds what FragGeneScan does not report — partial flags inferred from the start /
  stop codons — and, for genes without frameshifts, takes the nucleotide / protein strings from the
  genome at FragGeneScan's coordinates (upstream's minus-strand start refinement moves the coordinate
  but drops the codons in between from its own strings). The HMM is trained for code 11 only, so
  `-g` is ignored by this engine (with a warning).
- **`rna`** is **heuristic, not a covariance-model search.** The tRNA scanner requires a near-canonical
  cloverleaf (7-bp acceptor stem, 5/5 anticodon &amp; T-arm stems, conserved TΨC loop): high precision,
  modest sensitivity (0 false positives on phiX174 / λ; finds canonical tRNAs with the right anticodon).
  rRNA detection just anchors on conserved SSU/LSU motifs. For annotation use **tRNAscan-SE** and
  **barrnap** / **Infernal + Rfam**.

</details>

<details>
<summary>📚 <b>Citations</b></summary>

```bibtex
@article{hyatt2010prodigal,
  title={Prodigal: prokaryotic gene recognition and translation initiation site identification},
  author={Hyatt, Doug and Chen, Gwo-Liang and LoCascio, Philip F and Land, Miriam L and Larimer, Frank W and Hauser, Loren J},
  journal={BMC Bioinformatics}, volume={11}, pages={119}, year={2010}}

@article{camargo2023genomad,
  title={Identification of mobile genetic elements with geNomad},
  author={Camargo, Antonio Pedro and others}, journal={Nature Biotechnology}, year={2023},
  note={Prodigal-GV: alternative genetic codes for giant viruses and crassphages}}

@article{mcnair2019phanotate,
  title={PHANOTATE: a novel approach to gene identification in phage genomes},
  author={McNair, Katelyn and Zhou, Carol and Dinsdale, Elizabeth A and Souza, Brian and Edwards, Robert A},
  journal={Bioinformatics}, volume={35}, number={22}, pages={4537--4542}, year={2019}}

@article{besemer2001genemarks,
  title={GeneMarkS: a self-training method for prediction of gene starts in microbial genomes},
  author={Besemer, John and Lomsadze, Alexandre and Borodovsky, Mark},
  journal={Nucleic Acids Research}, volume={29}, number={12}, pages={2607--2618}, year={2001}}

@article{stanke2003augustus,
  title={Gene prediction with a hidden Markov model and a new intron submodel},
  author={Stanke, Mario and Waack, Stephan},
  journal={Bioinformatics}, volume={19}, number={suppl_2}, pages={ii215--ii225}, year={2003}}

@article{rho2010fraggenescan,
  title={FragGeneScan: predicting genes in short and error-prone reads},
  author={Rho, Mina and Tang, Haixu and Ye, Yuzhen},
  journal={Nucleic Acids Research}, volume={38}, number={20}, pages={e191}, year={2010}}

@article{vanderjeugt2022fraggenescanrs,
  title={FragGeneScanRs: faster gene prediction for short reads},
  author={Van der Jeugt, Felix and Dawyndt, Peter and Mesuere, Bart},
  journal={BMC Bioinformatics}, volume={23}, pages={198}, year={2022}}

@article{chan2019trnascanse2,
  title={tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes},
  author={Chan, Patricia P and Lowe, Todd M},
  journal={Nucleic Acids Research}, volume={47}, number={16}, pages={e76}, year={2019}}

@article{wright1990enc,
  title={The 'effective number of codons' used in a gene},
  author={Wright, Frank}, journal={Gene}, volume={87}, number={1}, pages={23--29}, year={1990}}

@misc{seemann2013barrnap,
  title={barrnap: BAsic Rapid Ribosomal RNA Predictor},
  author={Seemann, Torsten}, year={2013}, howpublished={\url{https://github.com/tseemann/barrnap}}}
```

</details>

---

<div align="center">
<sub><b>DYNAMITE</b> · RAW Lab · Richard Allen White III · GPL-3.0-or-later · github.com/raw-lab/dynamite</sub>
</div>
