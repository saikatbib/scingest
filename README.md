# scIngest

**Build Seurat and Scanpy objects from GEO accessions or local 10x directories — without writing code.**


---

## What this is

A wet-lab scientist finds a published scRNA-seq dataset, wants the counts matrix as a
Seurat or AnnData object, and hits a wall. The barrier is not conceptual — it is that GEO
supplementary files are unstandardised, and turning them into an object requires knowing
which of a dozen layouts a given submitter used.

scIngest closes exactly that gap: **accession or folder in, validated object out, with a
reproducible record of what happened.**

## What this is not

scIngest is an **ingestion tool, not an analysis tool**. It does not filter, normalise,
scale, integrate, cluster, or embed. The object it emits contains raw counts and metadata,
and nothing else.

This is deliberate. QC thresholds are a scientific judgement that depends on tissue,
chemistry, and study design. A GUI that ships a default `percent.mt < 5` would have that
default silently propagated into hundreds of papers by people who did not know it was a
choice. Refusing to filter is a feature.

## Design principles

1. **Never guess silently.** Every inference the format sniffer makes — sample grouping,
   matrix orientation, whether the values are actually raw counts — is shown to the user
   and is overridable before anything is built.
2. **A wrong object that looks right is worse than a crash.** Validation runs before any
   export is reported as successful.
3. **Reproducible by construction.** Every build emits a manifest and a runnable
   `build.py` / `build.R` that reproduces it outside the app.
4. **Errors are sentences, not stack traces.**

## Status

Phase 0 — `scingest-core`, the GUI-independent Python package. Nothing works yet.

## Repository layout

```
core/                  Python package — sniffer, readers, builders, provenance
  scingest_core/
  r/                   Seurat bridge script  (GPL-3, see below)
  tests/
docs/                  Technical specification and legal audit
app/                   Electron shell (Phase 1)
```

## Licensing

This repository is **dual-licensed by directory**:

| Path | License |
|---|---|
| Everything except `core/r/` | MIT — see `LICENSE` |
| `core/r/` | GPL-3.0 — see `core/r/LICENSE` |

`core/r/build_seurat.R` calls the R `Matrix` package, which is GPL (>= 2). Rather than
argue about whether ~60 lines of glue constitute a derivative work, that directory is
simply licensed GPL-3.

The application itself distributes no GPL code. R and Seurat are installed on demand from
conda-forge by the user's own machine. Seurat itself is MIT licensed.

Full analysis: `docs/LICENSING_AND_LEGAL_AUDIT.md`.

## Data, trademarks, and scope

- scIngest neither hosts, curates, nor vouches for GEO content. Users are responsible for
  the terms attaching to any dataset they download.
- Seurat and Scanpy are the work of their respective authors. This project is independent
  and is **not affiliated with or endorsed by** them.
- **Not validated for PHI or patient-identifiable data.** Do not point it at such data.

## Citation

If scIngest is useful to you, please cite it (see `CITATION.cff`) — and please cite Seurat,
Scanpy, AnnData, and the original study behind whatever accession you processed. A tool
that makes data easy to obtain should push people toward crediting those who generated it.
