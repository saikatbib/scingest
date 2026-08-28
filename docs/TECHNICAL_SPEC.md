# scIngest — Technical Specification

**A cross-platform desktop application for building Seurat and Scanpy objects from GEO accessions and local 10x directories.**

Version 0.2 (design draft) · 28 August 2026 · Author: drafted with Claude for Saikat

> **Name under review.** The earlier working name *CellForge* is unusable: it is a claimed trademark (CellForge™, iOrganBio, Oct 2025) and the name of a 2025 arXiv single-cell method. **scIngest** is the current placeholder; *CellIntake* and *scGather* are the alternatives. See `LICENSING_AND_LEGAL_AUDIT.md` §2. Nothing should be published under any name until clearance completes.

---

## 1. Purpose and scope

### 1.1 Problem statement

A wet-lab scientist finds a published scRNA-seq dataset, wants the counts matrix as a Seurat or AnnData object, and hits a wall. The barrier is not conceptual — it is that GEO supplementary files are unstandardised, and turning them into an object requires knowing which of a dozen layouts a given submitter used. The existing options are all wrong for this user:

- `GEOquery::getGEOSuppFiles()` downloads bytes and stops there.
- `Read10X()` works only when the files are already named and arranged exactly as CellRanger emits them.
- Web platforms (CELLxGENE, Broad SCP, UCSC Cell Browser) only host curated re-processed data; most GSEs never appear there.
- Every lab has a bespoke, undocumented `load_data.R` that breaks on the next dataset.

CellForge closes exactly this gap: **accession or folder in, validated object out, with a reproducible record of what happened.**

### 1.2 In scope (v1)

- Ingest from (a) a GEO Series or Sample accession, (b) a local directory, (c) drag-and-dropped files.
- Automatic detection of the layout and format of what was found.
- A mandatory human confirmation step before any object is built.
- Automatic extraction of sample-level metadata from GEO MINiML records.
- Export to `.h5ad` (AnnData) and `.rds` (Seurat), guaranteed equivalent.
- A provenance bundle: manifest, checksums, software versions, and generated `build.py` / `build.R` reproducing the run outside the app.

### 1.3 Explicitly out of scope (v1)

Per the scope decision, CellForge is an **ingestion tool, not an analysis tool**. It does not filter, normalise, scale, integrate, cluster, or embed. The object it emits contains raw counts and metadata, and nothing else.

Rationale worth stating in the README, because users will ask: QC thresholds are a scientific judgement that depends on tissue, chemistry, and study design. A GUI that ships a default `percent.mt < 5` will have that default silently propagated into hundreds of papers by people who did not know it was a choice. Refusing to filter is a feature.

**One boundary case to decide (see §12, D-4):** a *post-build report* — dimensions, per-sample cell/feature counts, sparsity, and distributions of counts-per-cell — is reporting, not processing. It changes no data. Recommendation: include it as a read-only summary so users can sanity-check what they built. If you want a stricter line, it drops out cleanly.

### 1.4 Target user

Someone who can install an application and identify their samples, but cannot debug an R error. Every failure mode must produce a sentence they can act on, not a stack trace. A secondary user — the bioinformatician they hand the object to — is served by the provenance bundle.

---

## 2. Locked decisions

| # | Decision | Choice | Rationale |
|---|---|---|---|
| L-1 | R runtime | Optional, micromamba-managed | `.h5ad` works with zero extra download; Seurat support installs on demand behind one button. Confirmed viable: `r-seurat` 5.5.1 is on conda-forge for win-64, osx-64, osx-arm64, and linux-64. |
| L-2 | Processing scope | Raw object only | No filtering, normalisation, or downstream steps. See §1.3. |
| L-3 | Shell framework | Electron + React + TypeScript | Mature cross-platform CI, mature component ecosystem, straightforward sidecar supervision. |
| L-4 | First deliverable | This document | Architecture reviewed before code. |
| L-5 | Python runtime (was D-1) | Pre-solved conda environment | One mechanism shared with the R install; avoids PyInstaller hidden-import work in scipy/h5py/numba across three platforms. Costs ~500 MB in the installer. |
| L-6 | Windows R (was D-2) | Native conda-forge `win-64`, no WSL | `r-seurat` 5.5.1 publishes a `win-64` build. **Must be smoke-tested in week 1 of Phase 1** — see §13. |
| L-7 | Extra modalities in h5ad (was D-3) | Companion `<name>_adt.h5ad`, referenced from the manifest | Avoids a MuData dependency and an unfamiliar format. Seurat needs no equivalent workaround — extra assays are native. Revisit before Phase 3 if MuData adoption grows. |
| L-8 | Post-build report (was D-4) | Yes, read-only | Dimensions, per-sample counts, sparsity, counts-per-cell distribution. Reporting only; transforms nothing. Does not breach L-2. |
| L-9 | Licensing (was D-5) | Open source: MIT root, GPL-3 for `core/r/` | Full audit in `LICENSING_AND_LEGAL_AUDIT.md`. **Release is blocked on institutional clearance** — see §16. |
| L-10 | Mitochondrial gene lists (was D-6) | Not in v1 | Annotating MT features is metadata, but `percent.mt` is the first step onto the analysis path and invites a filtering slider next to it. |

---

## 3. Architecture

### 3.1 Process model

```
┌──────────────────────────────────────────────────────┐
│ Electron main process (Node)                          │
│  · window + menu management                           │
│  · sidecar lifecycle: spawn, health-check, restart    │
│  · environment manager (micromamba orchestration)     │
│  · native file dialogs, auto-update                   │
└────────────┬─────────────────────────────────────────┘
             │ contextBridge IPC
┌────────────┴─────────────────────────────────────────┐
│ Renderer (React + TypeScript)                         │
│  · all UI; no direct filesystem or network access     │
│  · talks to the core over localhost HTTP + SSE        │
└────────────┬─────────────────────────────────────────┘
             │ 127.0.0.1:<ephemeral port>, bearer token
┌────────────┴─────────────────────────────────────────┐
│ cellforge-core — FastAPI sidecar (Python 3.11)        │
│  · GEO client, format sniffer, readers, builders      │
│  · job queue, progress events, cache management       │
│  scanpy · anndata · scipy · pandas · h5py · lxml      │
└────────────┬─────────────────────────────────────────┘
             │ subprocess, JSON manifest on stdin
┌────────────┴─────────────────────────────────────────┐
│ Rscript build_seurat.R — thin, ~60 lines              │
│  requires only: Seurat, Matrix                        │
└──────────────────────────────────────────────────────┘
```

### 3.2 Why this shape

**Python owns the logic.** Everything hard — sniffing, downloading, parsing, metadata joining, validation — happens once, in one language, in a package that is testable without launching a GUI. The Electron layer is a client.

**R is a printer, not a translator.** The obvious design is to write `.h5ad` and convert to Seurat via `SeuratDisk`, `sceasy`, or `zellkonverter`. Do not do this. Those bridges break across Seurat/anndata releases, silently drop layers, and mangle `obs` column types. Instead, the Python side writes an unambiguous plain-text intermediate — MatrixMarket + three CSVs — and the R script does nothing but read it and call `CreateSeuratObject()`. The R dependency surface shrinks to `Seurat` + `Matrix`, and the R code is short enough to read in full and verify by eye.

**Sidecar over embedding.** `rpy2` and `PyO3`-style embedding put R and Python in one address space, where a segfault in either takes down the app. Subprocesses fail independently and are killable from the UI.

### 3.3 Repository layout

```
cellforge/
├── core/                       # pip-installable, GUI-independent
│   ├── cellforge_core/
│   │   ├── geo/                # accession parsing, FTP listing, MINiML
│   │   ├── sniff/              # format detection cascade
│   │   ├── readers/            # mtx, h5, dense, h5ad, loom
│   │   ├── model/              # IngestPlan, SampleSpec, pydantic schemas
│   │   ├── build/              # AnnData assembly, R bridge, exporters
│   │   ├── provenance/         # manifest, checksums, script generation
│   │   └── server/             # FastAPI app, job queue, SSE
│   ├── r/build_seurat.R
│   └── tests/
│       ├── unit/
│       └── golden/             # real GSE fixtures, see §10
├── app/                        # Electron
│   ├── main/                   # sidecar supervisor, env manager, updater
│   ├── preload/
│   └── renderer/               # React screens
├── resources/
│   ├── micromamba/{win,mac,linux}
│   └── envs/                   # pre-solved environment specs (lockfiles)
└── build/                      # electron-builder config, CI, signing
```

The `core/` package must be independently usable: `pip install cellforge-core && cellforge build GSE123456 --out ./objects`. This is not a nice-to-have — it is what makes the sniffer testable at scale and gives you a CLI for power users at near-zero extra cost.

---

## 4. The format sniffer

This is where the product lives or dies. Everything else is plumbing.

### 4.1 Detection cascade

Applied per file, cheapest test first. Never trust the extension alone.

| Priority | Signature | Verdict |
|---|---|---|
| 1 | Magic bytes `\x89HDF\r\n\x1a\n` + `/matrix/data` dataset | CellRanger HDF5 v3 |
| 2 | Magic bytes HDF5 + single top-level genome group containing `genes`, `barcodes`, `data` | CellRanger HDF5 v2 |
| 3 | Magic bytes HDF5 + `/X` + `/obs` + `/var` | `.h5ad` (AnnData) |
| 4 | Magic bytes HDF5 + `/matrix` + `/row_attrs` + `/col_attrs` | `.loom` |
| 5 | Magic bytes HDF5 + `/assays` + `/meta.data` | `.h5Seurat` |
| 6 | First line matches `%%MatrixMarket matrix coordinate` | MatrixMarket |
| 7 | RDS magic (`\x1f\x8b` → decompress → `X\n\x00\x00\x00`) | R serialised object |
| 8 | tar magic at offset 257 (`ustar`) | Archive — recurse into members |
| 9 | Text, ≤4 columns, no header, ≥90% of col 1 matching a barcode or gene pattern | Sidecar TSV (barcodes / features) |
| 10 | Text, wide, first row and column look like labels | Dense matrix |

Compressed files (`.gz`, `.bz2`, `.xz`, `.zst`) are transparently decompressed for sniffing; the first 64 KB is enough for every test above.

### 4.2 Grouping MTX triplets into samples

The single most common GEO layout. Files carry a submitter-chosen prefix:

```
GSM4008658_ACC1_barcodes.tsv.gz
GSM4008658_ACC1_features.tsv.gz
GSM4008658_ACC1_matrix.mtx.gz
```

Algorithm:

1. For each candidate file, strip a trailing token matching
   `[._-]?(matrix\.mtx|barcodes\.tsv|features\.tsv|genes\.tsv|counts\.mtx)(\.(gz|bz2|xz|zst))?$`
   (case-insensitive). The remainder is the **group key**.
2. Group by key. A complete group has exactly one matrix, one barcode file, and one feature file.
3. Incomplete groups (2 of 3) are surfaced as warnings, not silently dropped — a missing features file is often present under a name the regex did not anticipate, and the UI lets the user assign it manually.
4. Derive a display name: prefer an embedded GSM ID matched against MINiML titles; else the group key with the GSM prefix stripped; else the folder name.

Known variants the regex must tolerate: `genes.tsv` (v2), no separator (`GSM123barcodes.tsv.gz`), capitalised (`Barcodes.tsv`), `.mtx` without `matrix` in the name, per-sample subdirectories with unprefixed canonical names, and the CellRanger convention of a `filtered_feature_bc_matrix/` or `raw_feature_bc_matrix/` directory (prefer *filtered* by default, flag the choice in the UI).

### 4.3 Orientation and dimension resolution

MatrixMarket is a bare `rows cols nnz` header with no semantics. 10x writes features × barcodes; other tools write the transpose.

Resolution, in order:
1. If `len(features) == n_rows` and `len(barcodes) == n_cols` → features × barcodes. Standard.
2. If `len(features) == n_cols` and `len(barcodes) == n_rows` → transposed; transpose on read.
3. If both sidecar lengths equal both dimensions (a square matrix) → unresolvable automatically; ask.
4. If neither matches → hard error naming the three numbers. This almost always means the wrong files were grouped together, and saying so is more useful than any guess.

For **dense text matrices** there are no sidecar files, so the heuristic is content-based:
- Barcode-like: `^[ACGT]{12,20}(-\d+)?$` — count matches in the first row and first column.
- Gene-like: matches `^ENS[A-Z]*G\d{6,}` (Ensembl), or is in a bundled HGNC/MGI symbol set.
- Whichever axis scores higher on gene-likeness is the feature axis.
- Tie-break: the axis whose length falls in 10,000–60,000 is more likely features.
- Fallback default: genes × cells (the prevailing GEO convention for dense text).
- **The result is always shown and always flippable.** The heuristic is a suggestion.

### 4.4 Is this actually raw counts?

Submitters frequently upload normalised data labelled "counts". Detect and warn — never silently correct.

Sample up to 10⁵ non-zero values and evaluate:

| Observation | Conclusion | Action |
|---|---|---|
| All integral, min ≥ 0 | Raw counts | Proceed |
| Non-integral, min ≥ 0, max ≲ 15 | Probably log1p-transformed | Warn prominently; offer to store in a `logcounts` layer rather than `X` |
| Non-integral, min ≥ 0, max large, row sums ≈ constant | Probably CPM/TPM | Warn |
| Any negative value | Scaled/centred (z-scored) | Warn strongly; this cannot be used as counts |
| Integral but max < 10 and ≥99% zeros | Possibly binarised or ATAC | Warn |

The verdict is written to the manifest as `matrix_value_class` so a downstream analyst can see what CellForge thought.

### 4.5 Feature table handling

- **v2 (`genes.tsv`, 2 columns):** `[ensembl_id, symbol]`. Synthesise `feature_type = "Gene Expression"`.
- **v3 (`features.tsv`, 3 columns):** `[id, name, feature_type]`.
- **4+ columns:** CellRanger ARC / ATAC peak annotations — read the first three, retain the rest as extra `var` columns.
- **Mixed feature types** (Gene Expression + Antibody Capture + CRISPR Guide Capture + Multiplexing Capture): split by type. Gene Expression becomes the primary matrix; each other type becomes an additional modality — a second AnnData in `.h5ad` terms (written as a companion file or `.obsm` block, see D-3), and an additional `Assay` (`ADT`, `HTO`, `CRISPR`) in the Seurat object. This is common enough in CITE-seq and hashing datasets that ignoring it would be a visible gap.
- **Species inference:** `ENSG`/`ENSMUSG`/`ENSDARG` prefixes are decisive. Failing that, the ratio of all-uppercase to Title-case symbols separates human from mouse. Recorded in `uns['species_guess']`; never used to alter data.

### 4.6 Name sanitisation — the cross-language correctness trap

**This is the subtlest bug in the whole design and it must be handled once, in Python, before either exporter runs.**

Three transformations differ between the ecosystems:

1. **Duplicate feature names.** Scanpy offers `var_names_make_unique()` (appends `-1`, `-2`); R's `CreateSeuratObject` applies `make.unique()` (appends `.1`, `.2`). Left to each side, the same gene ends up as `TBCE-1` in the `.h5ad` and `TBCE.1` in the `.rds`.
2. **Underscores.** Seurat rejects underscores in feature names and silently rewrites them to dashes, emitting only a warning: *"Feature names cannot have underscores ('_'), replacing with dashes ('-')"*. AnnData does not care. Any hashtag antibody named `HTO_1` diverges immediately.
3. **Barcode collisions on merge.** Seurat's `merge()` prefixes via `add.cell.ids`; scanpy's `concat` suffixes via `index_unique`. Different separators, different order.

**Rule:** Python performs all three normalisations exactly once, records every rename in `manifest.renames`, and hands the final strings to both exporters. The R script must be given already-unique, already-sanitised names and must not be allowed to alter them. Cell IDs are formed in Python as `{sample_id}_{barcode}` and `add.cell.ids` is left `NULL`.

The verification step (§8) exists largely to catch regressions here.

---

## 5. GEO ingestion

### 5.1 Resolving an accession

Accepted input: `GSE123456`, `GSM4008658`, a full GEO URL, or a whitespace/comma-separated list.

Metadata is fetched via NCBI E-utilities (`esearch` + `esummary` against the `gds` database). Files are fetched over HTTPS from the FTP tree, using the documented bucketing convention — **the last three digits of the accession are replaced by `nnn`**:

```
Series supplementary : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/
Series MINiML        : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/miniml/GSE123456_family.xml.tgz
Series SOFT          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/soft/GSE123456_family.soft.gz
Sample supplementary : https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4008nnn/GSM4008658/suppl/
```

Note the asymmetry: series files are frequently bundled as one `GSE123456_RAW.tar`, while the per-GSM `suppl/` directories hold the same content unbundled. **Prefer per-GSM downloads when the user wants a subset of samples** — it avoids pulling a 40 GB tar to extract 2 GB. Offer the tar only when the user wants everything, or when per-GSM files are absent.

### 5.2 Rate limits and etiquette

E-utilities permits 3 requests/second without an API key and 10/second with one. Expose an optional NCBI API key field in Settings. Set a descriptive `User-Agent` including the app name and version. Serialise FTP listing requests. Downloads are resumable via HTTP `Range` and are checksummed on completion.

### 5.3 Metadata extraction — the second high-value feature

The MINiML family XML contains, per sample: title, `source_name`, organism, library strategy, platform, and an arbitrary number of `Characteristics` tags with submitter-defined `tag` attributes (`tissue`, `treatment`, `patient_id`, `timepoint`, `condition`, …).

CellForge parses these into a tidy per-sample table and offers each column for inclusion in `obs` / `meta.data`. Presented as a preview grid where the user can rename columns, drop them, and correct types (categorical vs numeric). Values are broadcast to every cell of the corresponding sample.

This is worth building carefully. It is the step users are most likely to do by hand, most likely to get wrong, and least likely to document. Characteristic tags are inconsistently cased and often contain `key: value` inside a single string — normalise keys to snake_case and split on the first colon when the tag is generic.

### 5.4 Caching

Downloads are content-addressed under the app data directory, keyed by URL plus ETag/size. Re-running the same accession never re-downloads. The Dashboard shows cache size with a per-accession eviction control — scRNA-seq downloads are large and users will fill their disk otherwise.

---

## 6. The build path

### 6.1 Intermediate representation

Everything converges on an in-memory `AnnData` per sample:

- `X` — `scipy.sparse.csr_matrix`, cells × genes, raw counts, dtype preserved (`int32` where it fits, else `float32`).
- `obs` — barcode index plus `sample_id` plus joined MINiML metadata.
- `var` — feature index (sanitised symbols) plus `gene_id`, `feature_type`, and any extra columns.
- `uns` — provenance block: source URLs, sniffer verdicts, species guess, value class, app version.
- Additional modalities held alongside as separate `AnnData` objects keyed by feature type.

Multi-sample builds concatenate with `anndata.concat(join="outer", label="sample_id")`. Outer join matters: samples from different CellRanger reference versions have non-identical feature sets, and an inner join silently discards genes. Missing entries are filled with zero and the count of unshared features is reported in the UI.

### 6.2 Writing `.h5ad`

Direct: `adata.write_h5ad(path, compression="gzip")`. For very large objects, `compression="lzf"` trades size for speed — expose as a setting.

### 6.3 Writing `.rds` via the R bridge

Python writes a staging directory:

```
_staging/
├── matrix.mtx          # MatrixMarket, features × cells (R/Seurat orientation)
├── features.csv        # final sanitised names + gene_id + feature_type
├── barcodes.csv        # final cell IDs
├── metadata.csv        # obs, one row per cell, in barcode order
├── modalities/         # optional: adt/matrix.mtx, adt/features.csv, ...
└── manifest.json       # dims, checksums, assay name, object version target
```

`build_seurat.R` then:

```r
suppressPackageStartupMessages({ library(Matrix); library(Seurat) })
args <- commandArgs(trailingOnly = TRUE)
staging <- args[1]; out <- args[2]
mf <- jsonlite::fromJSON(file.path(staging, "manifest.json"))

if (identical(mf$object_version, "v4")) {
  options(Seurat.object.assay.version = "v3")   # emit legacy Assay, not Assay5
}

counts <- readMM(file.path(staging, "matrix.mtx"))
features <- read.csv(file.path(staging, "features.csv"), stringsAsFactors = FALSE)
barcodes <- read.csv(file.path(staging, "barcodes.csv"), stringsAsFactors = FALSE)
rownames(counts) <- features$name          # already unique and sanitised
colnames(counts) <- barcodes$cell_id
counts <- as(counts, "CsparseMatrix")

meta <- read.csv(file.path(staging, "metadata.csv"), row.names = 1,
                 stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(rownames(meta), colnames(counts)))

obj <- CreateSeuratObject(counts = counts, meta.data = meta,
                          project = mf$project,
                          min.cells = 0, min.features = 0)   # no filtering — L-2

# optional additional modalities
for (m in mf$modalities) { ... CreateAssayObject(...) ... }

obj@misc$cellforge <- mf$provenance
saveRDS(obj, file = out, compress = TRUE)
cat(sprintf("OK\t%d\t%d\n", nrow(obj), ncol(obj)))
```

`min.cells = 0, min.features = 0` is load-bearing. Seurat's defaults are already 0, but stating them explicitly documents the no-filtering guarantee and protects against a future default change.

**Seurat v4 vs v5.** Seurat 5 stores counts in `Assay5` objects with layers; Seurat 4 pipelines cannot read them. Expose a radio button — *Seurat v5 (default)* / *Seurat v4-compatible* — implemented via `options(Seurat.object.assay.version = "v3")`. Many labs are still on v4 and this one toggle prevents a whole class of support requests.

### 6.4 Other export targets

- `.qs` / `.qs2` — optional, dramatically faster load than `.rds` for large objects. One extra R dependency.
- `.h5Seurat` — **not recommended.** SeuratDisk is effectively unmaintained. Offer only if users demand it, and label it as such.
- `.loom` — cheap to add via `loompy`; low demand. Defer.

---

## 7. Environment management

### 7.1 Python

Ship a pre-solved conda environment (built with `conda-lock` + `constructor`, or extracted from a `conda-pack` tarball at install time) containing Python 3.11, scanpy, anndata, scipy, h5py, pandas, lxml, and FastAPI. Roughly 400–600 MB installed. This is present after install; the app is usable immediately.

An alternative is PyInstaller-freezing the core, which produces a smaller bundle but requires fighting hidden-import issues in scipy, h5py and numba on three platforms. See D-1.

### 7.2 R, on demand

A per-platform `micromamba` binary ships in `resources/`. When the user first requests a Seurat export, or clicks **Enable Seurat support** in Settings:

```
micromamba create -p <appdata>/envs/r -c conda-forge \
    r-base=4.4 r-seurat=5.5.1 r-matrix r-jsonlite --yes
```

Verified available: `r-seurat` 5.5.1 is published on conda-forge for `win-64`, `osx-64`, `osx-arm64`, and `linux-64`, so a single recipe covers every target. Roughly 1.2–1.8 GB, several minutes on a normal connection.

Requirements for this flow:
- Streamed progress with a real percentage, not a spinner. Users abandon silent multi-minute downloads.
- Fully cancellable and resumable; partial environments are removed on failure.
- Pinned versions in a lockfile, so two users on the same app version get identical R environments.
- A **Detect existing R** escape hatch in Settings for users who already have a suitable Seurat installation and object to the download. Validate it by running a probe script that checks `packageVersion("Seurat")`.
- Offline machines: document a manual environment path. Institutional clusters and secure environments often block conda channels.

---

## 8. Validation and verification

Before an export is reported as successful:

1. **Structural** — `n_obs`, `n_vars`, and `nnz` identical between the AnnData and the Seurat object (the R script returns its dimensions on stdout; Python asserts).
2. **Identity** — sorted `var_names` and sorted `cell_ids` hash-match across both outputs.
3. **Numeric spot-check** — 1,000 random non-zero coordinates compared value-for-value by re-reading both objects. Cheap, and it catches transposition and dtype bugs that dimension checks miss.
4. **Metadata** — every `obs` column present in both, with matching value counts per column.

Failures block the export and produce a plain-language message. This suite is also the core regression test.

---

## 9. User interface

### 9.1 Design language

Clean, dense, scientific. White or near-white background; a restrained accent (a single blue) with colourblind-safe categorical colours where any encoding is needed; Inter or IBM Plex Sans; generous table density rather than card-based dashboards. Any plot follows publication conventions — no chartjunk, no gradients, embedded fonts on export. Dark mode as a preference, not a default.

### 9.2 Screens

**1 — Dashboard**
Recent jobs with status and duration; the object library (path, samples, cells, features, formats present, build date, one-click "reveal in folder"); cached downloads with sizes and eviction; environment status card (Python ✓, R ✓/Install).

**2 — New build**
Three entry points: accession field with inline validation; folder picker; drag-and-drop zone. Recent accessions are remembered.

**3 — Discovery** *(GEO path only)*
The supplementary file listing: name, size, sniffer verdict, and per-sample grouping preview. Checkboxes to include or exclude samples *before* downloading — critical for large series. A running total of bytes to download.

**4 — Sample map** *(the critical screen)*
An editable table, one row per detected sample:

| Sample name | Matrix file | Barcodes | Features | Orientation | Features detected | Cells detected | Value class |
|---|---|---|---|---|---|---|---|

Every cell is editable. Unassigned or ambiguous files sit in a side panel and can be dragged onto a row. Warnings appear inline and are dismissible only by acknowledgement. **The Build button is disabled until every row resolves cleanly or the user explicitly overrides a warning.** This screen is the product's trust boundary; give it the most design attention.

**5 — Metadata** *(GEO path only)*
The MINiML characteristics grid: sample rows × characteristic columns, with per-column include/rename/type controls and a live preview of the resulting `obs` columns.

**6 — Build**
Format checkboxes (`.h5ad`, `.rds`, `.qs`); Seurat v5/v4 radio; output directory; project name. Then a live log with a real progress bar, per-stage, cancellable.

**7 — Result**
Dimensions, per-sample breakdown, validation results (all green, or what failed), buttons to open the output folder, view the manifest, and copy the generated `build.py` / `build.R`.

### 9.3 Error philosophy

Every user-facing error carries three things: what happened in plain words, why CellForge thinks it happened, and what the user can do. The raw traceback lives behind a "Technical details" disclosure and a "Copy report" button that bundles logs, manifest, and versions for a GitHub issue.

Example — do this:

> **Couldn't match the features file to the matrix for sample ACC1.**
> The matrix has 32,738 rows, but the features file lists 33,694 genes. These usually come from different CellRanger runs. Check whether `GSM4008658_ACC1_features.tsv.gz` belongs with this matrix, or assign a different file on the Sample map screen.

Not this:

> `ValueError: shape mismatch (32738,) vs (33694,)`

---

## 10. Testing

### 10.1 Golden corpus

The regression suite is a curated set of real GEO accessions, one or more per layout, with expected dimensions committed as fixtures. Layouts to cover:

| Case | What it exercises |
|---|---|
| `GSE_RAW.tar` with per-sample prefixed triplets | Archive recursion + grouping |
| Per-GSM `suppl/` triplets, no tar | GSM-level download path |
| CellRanger v2 (`genes.tsv`, 2-column) | Legacy feature table |
| CellRanger v3 (`features.tsv`, 3-column) | Standard path |
| `filtered_feature_bc_matrix.h5` | HDF5 v3 reader |
| CITE-seq with Antibody Capture rows | Modality splitting |
| Cell hashing with `HTO_` feature names | Underscore sanitisation |
| Dense `.csv.gz`, genes × cells | Orientation heuristic |
| Dense `.txt.gz`, cells × genes | Orientation heuristic, transposed |
| Log-normalised values in a file named "counts" | Value-class warning |
| Pre-built `.h5ad` supplementary | Passthrough path |
| Duplicated gene symbols | Uniquification parity |
| Samples with differing reference versions | Outer-join concatenation |
| >1M cells | Memory and chunking |

**Curate the actual accession list in Phase 0 by browsing GEO** rather than trusting any list from memory — accessions get updated and re-uploaded, and a fixture that does not match reality is worse than none. Store only checksums and expected dimensions in the repo; download fixtures in CI with caching.

### 10.2 Levels

- **Unit** — sniffer verdicts against synthetic files for every signature in §4.1; regex grouping against a corpus of ~200 real GEO filenames.
- **Integration** — full build from small fixtures, asserting the §8 validation suite.
- **Parity** — the h5ad/rds equivalence checks, run on every golden case that has R available.
- **UI** — Playwright over the packaged Electron app for the happy paths on all three platforms.

---

## 11. Packaging and distribution

`electron-builder` targeting NSIS (Windows), DMG + universal binary (macOS arm64 + x64), and AppImage + `.deb` (Linux). Auto-update via `electron-updater` against GitHub Releases.

**Code signing is not optional if the audience is non-computational.** An unsigned Windows installer produces a SmartScreen block that most biologists will not click through; an unsigned macOS app is refused outright by Gatekeeper. Budget an Apple Developer membership (~$99/yr) and a Windows OV or EV code-signing certificate (~$200–500/yr; EV avoids SmartScreen reputation build-up but requires hardware token handling in CI). Notarisation of the macOS build must be part of the release pipeline, not a manual step.

Installer size expectation: ~150 MB Electron shell + ~500 MB Python environment ≈ **650 MB**, with the R environment adding ~1.5 GB on demand. State this on the download page; scientists on institutional networks plan around it.

---

## 12. Licensing and legal constraints on the design

D-1 through D-6 are resolved and now appear as L-5 through L-10 in §2. The full analysis lives in the companion document `LICENSING_AND_LEGAL_AUDIT.md`; what follows is only what constrains the engineering.

### 12.1 Seurat is MIT, R is not

Seurat itself is MIT licensed (Copyright © 2021 Seurat authors), which removes the copyleft concern from the package the app actually calls. But **R itself is GPL-2 | GPL-3 and the `Matrix` package is GPL (>= 2)**, so the R environment as a whole is GPL.

### 12.2 Three design rules that exist for licensing reasons

These are not stylistic preferences. Violating any of them changes the project's legal position.

1. **The R interface stays at arm's length.** Files in, exit code and a status line out. Separate programs communicating via pipes, arguments, and temp files are separate works; components sharing an address space are not.
2. **Never use `rpy2` or any in-process R embedding.** This would put GPL code in the application's address space and destroy the separation in rule 1. The sidecar architecture is now load-bearing for two reasons, not one.
3. **`build_seurat.R` stays thin and free of application logic**, both so it is plainly an independent script and so that GPL-3 licensing that one directory costs nothing.

### 12.3 Repository license structure

```
LICENSE            MIT      — everything except core/r/
core/r/LICENSE     GPL-3    — the Seurat bridge script (it calls Matrix, GPL >= 2)
```

The split is explained in three sentences in the README. Dual-licensing by directory is a normal pattern and is cheaper than arguing about whether 60 lines of glue are derivative.

### 12.4 L-1 now has a legal dimension

Because R is installed on demand from conda-forge by the user's own machine, **the installer distributes no GPL code and you are not a distributor of R.** If R were ever bundled instead, GPL-3 §6 obligations would attach — conveying corresponding source, or a written offer valid for three years, for every GPL component shipped.

L-1 was chosen for installer size. It should not be revisited on size grounds alone.

### 12.5 Constraints carried into other sections

- **No GEO data in the repository** (§10.1). Already the plan for size reasons; also the right call legally. Fixtures download in CI; only checksums and expected dimensions are committed.
- **Honest `User-Agent` and rate-limit compliance** (§5.2). Abusive traffic from a distributed application gets an IP range blocked, affecting every user simultaneously.
- **Nominative use only for Seurat and Scanpy** — descriptive references are fine; no logos, no use in the product or repository name, no implied endorsement, and a disclaimer in the README.
- **Fonts** (§9.1) — Inter and IBM Plex are SIL OFL 1.1. Bundle freely; retain notices.
- **CI license gates** — `pip-licenses` and `license-checker` run on every build against an allowlist, and `THIRD_PARTY_NOTICES.md` is regenerated on every dependency change rather than maintained by hand.

### 12.6 Remaining open items

| # | Question | Status |
|---|---|---|
| O-1 | Final project name | **Blocked.** *CellForge* is unusable. *scIngest* is the placeholder pending USPTO / GitHub / PyPI / npm / domain / Scholar clearance. |
| O-2 | Copyright holder on the LICENSE file | **Blocked** on §16. Likely the institution rather than an individual. |
| O-3 | Revisit MuData for modalities | Deferred to Phase 3. |

---

## 13. Risks

| Risk | Severity | Mitigation |
|---|---|---|
| GEO format heterogeneity exceeds the sniffer | High | Never claim full coverage. The Sample map screen's manual assignment is the universal fallback; make it good enough that a failed sniff is an inconvenience, not a dead end. |
| Windows R compilation failures | High | Validate the conda-forge `win-64` Seurat install in the first week of Phase 1, before any UI work. |
| Memory blow-up on large series | Medium | Chunked MatrixMarket reading; `backed="r"` for h5ad; stream to disk rather than concatenating in RAM above a size threshold; warn before starting a build the machine cannot hold. |
| Seurat / anndata API drift | Medium | Pin everything in lockfiles; the parity suite in CI catches breakage before users do. |
| NCBI blocks the app's traffic | Medium | Respect rate limits, identify the client honestly in `User-Agent`, support an API key, back off on 429. |
| Silent scientific error (wrong orientation, normalised data treated as counts) | **Critical** | The confirmation screen, the value-class warnings, and the provenance manifest exist for this. A wrong object that looks right is far worse than a crash. |

---

## 14. Roadmap

| Phase | Deliverable | Rough effort |
|---|---|---|
| **0** | `cellforge-core` as a CLI: sniffer, readers, AnnData builder, R bridge, validation suite, golden corpus | 3–4 weeks |
| **1** | Electron shell, sidecar supervision, local-folder path end to end, environment manager | 3–4 weeks |
| **2** | GEO ingestion: listing, download, cache, MINiML metadata, Discovery + Sample map screens | 3–4 weeks |
| **3** | Multi-sample merge, modalities, export options, result screen, generated scripts | 2–3 weeks |
| **4** | Packaging, signing, notarisation, auto-update, cross-platform UI tests, docs | 2–3 weeks |

Phase 0 first is deliberate. It produces a genuinely useful tool on its own, it is where all the difficult logic lives, and it gives you a test suite that makes every later phase safe to change.

---

## 15. Immediate next steps

1. **Institutional clearance** (§16) — blocks everything public.
2. **Name clearance** — settle O-1 and reserve the GitHub, PyPI, and npm names.
3. Curate the golden accession list — 15–20 real GSEs spanning the §10.1 table.
4. Smoke-test the riskiest technical assumption: `micromamba create -c conda-forge r-base r-seurat` on Windows, then round-trip a small matrix through `.h5ad` and `.rds` and diff them.
5. Begin Phase 0 in a private repository while 1 and 2 resolve.

---

## 16. Publishing to GitHub

### 16.1 Clearance gate

MD Anderson's Intellectual Property Policy states that IP created within the course of employment is owned by the Board, that it must be disclosed to the Office of Technology Commercialization via an invention disclosure report, and that **no public disclosure should be made until authorized**. A public repository is a public disclosure, and an MIT license is a grant of rights that only the copyright holder can make.

This is very likely a routine approval — institutions release research software openly all the time, and an ingestion utility with no patentable claim is the easy case. But it must be the OTC's decision, in writing, and it determines the copyright line (O-2). Development in a **private** repository in the meantime is the sensible path.

### 16.2 Mechanics, once cleared

Verified by direct test in both shells:

| Capability | Cloud container | Linked Windows device (`device_bash`) |
|---|---|---|
| `git` | ✓ | ✓ |
| `curl`, `wget`, `python3` | ✓ | ✓ |
| `gh` CLI | ✗ not installed; `cli.github.com` blocked | ✗ not installed (not needed) |
| `git clone` / `push` over HTTPS | ✓ | ✓ |
| `api.github.com` | ✗ 403 from the egress proxy | **✓ 200** — `POST /user/repos` returns 401 (auth required), not a proxy block |

**The device shell has full GitHub API access.** Repository creation, topics, releases, pull requests, branch protection, and Actions configuration are all automatable from there via `curl` against `api.github.com`; `gh` is a convenience, not a requirement.

**Therefore the repository lives on the device, not in the cloud container** — under `D:\Saikat\Claude\<name>`. This is the better arrangement anyway: the working tree sits on Saikat's own disk where he can open it in an editor, and all git and API operations run from the same place.

The cloud container remains useful for anything needing PyPI or tooling the device VM lacks; files move across with the staging tools when a step genuinely requires it.

Working arrangement:

1. Saikat issues a **fine-grained** personal access token. For repository creation it needs `Administration: read and write` on the target account plus `Contents: read and write`; if he creates the empty repo himself first, `Contents: read and write` scoped to that one repository is sufficient and is the safer option.
2. The token is stored on the device in a file outside the mounted folders (so it is never committed and never staged), and consumed via `git`'s credential helper and `curl -H "Authorization: Bearer …"`.
3. Revoked when the work is done.

**Token hygiene, stated plainly:** any token pasted into the conversation is recorded in the transcript. Use a fine-grained token, never a classic token with full account scopes; prefer the single-repository scope by pre-creating the repo; set the shortest workable expiry; revoke on completion. A token that can only write to one repository has a worst case of a revert.

---

### Sources consulted

- [Programmatic access to GEO — NCBI](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html) — FTP URL patterns and the `nnn` bucketing convention
- [GEO FTP site](https://www.ncbi.nlm.nih.gov/geo/info/download.html) — download layouts
- [r-seurat on conda-forge](https://anaconda.org/conda-forge/r-seurat) — version 5.5.1, platform coverage
- [satijalab/seurat issue #1648](https://github.com/satijalab/seurat/issues/1648) — underscore-to-dash feature name rewriting
- [Seurat v5 essential commands](https://satijalab.org/seurat/articles/seurat5_essential_commands) — assay version option and v5 layer model
- [Seurat license](https://satijalab.org/seurat/license) — MIT
- [Matrix DESCRIPTION](https://github.com/cran/Matrix/blob/master/DESCRIPTION) — GPL (>= 2)
- [MD Anderson Intellectual Property Policy](https://www.mdanderson.org/about-md-anderson/business-legal/legal-and-policy/legal-statements/intellectual-property-policy.html) — ownership and pre-disclosure authorization

Full licensing analysis: `LICENSING_AND_LEGAL_AUDIT.md`
