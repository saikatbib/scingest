# Licensing and Legal Audit

**Project:** scRNA-seq object builder (working name under review — see §2)
**Date:** 28 August 2026
**Status:** Pre-release review. Two blocking items.

> I am not a lawyer and this is not legal advice. It is a technical review of license texts and published policies, intended to tell you which questions are settled, which need a human decision, and which need your institution's counsel. Items marked **BLOCKER** should be cleared by a person with authority before any public release.

---

## 1. Headline findings

| # | Finding | Severity |
|---|---|---|
| F-1 | MD Anderson's IP policy claims ownership of employee-created IP and requires disclosure to the Office of Technology Commercialization **before** any public disclosure. Publishing to a personal GitHub account without clearance would conflict with this. | **BLOCKER** |
| F-2 | "CellForge" is already in use twice in this exact field, once as a claimed trademark (CellForge™, iOrganBio, launched Oct 2025) and once as a 2025 arXiv single-cell method. The name must change. | **BLOCKER** |
| F-3 | Seurat is **MIT licensed**, not GPL. This removes the largest copyleft risk in the design. | Positive |
| F-4 | GPL exposure is confined to R itself and the `Matrix` package. The on-demand install architecture (decision L-1) means you never distribute either one. | Manageable |
| F-5 | Electron bundles ffmpeg under LGPL-2.1+; notice obligations apply. Standard and routine, but must not be skipped. | Minor |

---

## 2. F-2 — the name

Two prior uses, both in single-cell biology:

- **CellForge™** — iOrganBio's AI cell-manufacturing platform, launched from stealth in October 2025 with $2M funding, using the ™ symbol.
- **CellForge** — *CellForge: Agentic Design of Virtual Cell Models*, arXiv:2508.02276 (2025).

The trademarked commercial use is the real problem. Trademark protection turns on likelihood of confusion within a field, and "single-cell biology software" is unambiguously the same field. A ™ claim without registration is weaker than a registered mark, but not worth testing, and rebranding after a paper and a user base exist is far more expensive than picking a different name now.

**Candidates that returned no conflicting software in a search of the single-cell ecosystem:**

| Name | Notes |
|---|---|
| **scIngest** | Descriptive, uses the standard `sc` prefix convention, no hits. My preference. |
| **CellIntake** | Plain-language, communicates the scope limit (ingestion, not analysis), no hits. |
| **scGather** | No hits, slightly softer. |
| **MatrixMaker** | Memorable but generic; higher trademark risk in software generally. |

Names to avoid: anything containing *Seurat* or *Scanpy* (see §6), and *Anvil* (collides with NHGRI's AnVIL platform).

**Before committing to any name, run:** a USPTO TESS search, a GitHub org/repo search, a PyPI and npm name check, a domain check, and a Google Scholar search. Cheap now, expensive later.

---

## 3. F-1 — institutional ownership

MD Anderson's published Intellectual Property Policy states that intellectual property created by an individual within the course of their employment responsibilities is owned by the Board, that employees must disclose IP via an invention disclosure report, that disclosure should occur prior to any public disclosure or publication, and that no public disclosure should be made until authorized by MD Anderson. The Office of Technology Commercialization administers this.

Pushing this repository public under your personal account is a public disclosure. Two consequences:

1. **You may not have the right to grant an MIT license** on institutionally-owned code. An open-source license is a grant of rights by the copyright holder; if the Board holds copyright, only the Board can grant it.
2. Doing it first and asking later is the failure mode that damages relationships with a tech transfer office, and it is not reversible — you cannot un-publish source that has been cloned.

**This is very likely a routine approval.** Institutions release research software openly all the time, and a data-ingestion utility with no patentable claim and no commercial pathway is the easy case for an OTC. But it needs to be the OTC's decision, in writing, not an assumption.

**What to do:** submit an invention/software disclosure to the OTC describing the tool, state that you intend an open-source release under MIT, and ask for written authorization. Ask specifically whether they want a copyright line naming the Board (typical form: `Copyright (c) 2026 The University of Texas M.D. Anderson Cancer Center`) rather than your name.

Two situations change this analysis, and only you know which applies:

- If the tool is genuinely developed outside the scope of your employment, on personal time and personal equipment, and does not use institutional resources or unpublished institutional data, the ownership claim weakens considerably. Note that your prior side projects establish that you do build software independently, which helps this argument — but "related to my job function as a computational cancer researcher" is exactly the language such policies turn on, and a scRNA-seq tool is squarely within it. Do not rely on this without confirming.
- If you build it after leaving MD Anderson, or under the India lab, the question moves to that institution's rules.

**I will not push anything to a public repository until you confirm this is cleared.** A private repository under your account is lower risk and is a reasonable place to develop in the meantime — but confirm even that, since some institutions treat private repos outside institutional control as an issue for data-handling reasons.

---

## 4. Dependency license inventory

### 4.1 Python core — all permissive

| Package | License | Notes |
|---|---|---|
| CPython | PSF-2.0 | Permissive, GPL-compatible |
| numpy, scipy, pandas | BSD-3-Clause | |
| scanpy | BSD-3-Clause | |
| anndata | BSD-3-Clause | |
| h5py | BSD-3-Clause | Wraps HDF5 (HDF5 License — BSD-style, permissive) |
| lxml | BSD-3-Clause | Wraps libxml2 / libxslt (MIT) |
| FastAPI, Starlette, pydantic | MIT | |
| uvicorn | BSD-3-Clause | |
| httpx / requests | BSD-3 / Apache-2.0 | |
| GEOparse | BSD-3-Clause | Verify at pin time |
| tqdm, natsort | MIT / MIT | |

**Verdict:** nothing here constrains the project license. MIT or BSD-3 for the core is clean.

### 4.2 R environment — GPL, but not distributed

| Package | License | |
|---|---|---|
| R (r-base) | GPL-2 \| GPL-3 | ⚠ Copyleft |
| Matrix | GPL (>= 2) | ⚠ Copyleft |
| **Seurat** | **MIT** | ✓ Copyright (c) 2021 Seurat authors |
| SeuratObject | MIT | Verify at pin time |
| jsonlite | MIT | |
| Seurat's transitive deps (Rcpp, RcppEigen, igraph, sctransform, …) | Mixed, largely GPL-2/GPL-3 | ⚠ |

The R environment taken as a whole is effectively GPL. Seurat being MIT is a genuine relief — it means the package you actually call is permissive — but it does not sanitise the environment, because `Matrix` and R itself are not.

### 4.3 Desktop shell

| Component | License | Notes |
|---|---|---|
| Electron | MIT | Bundles Chromium (BSD-3 + a long notice file) and Node (MIT) |
| — ffmpeg inside Electron | LGPL-2.1+ | ⚠ Shipped as a separate shared library; notice + relinking rights required |
| React, React DOM | MIT | |
| electron-builder | MIT | |
| micromamba / mamba | BSD-3-Clause | |
| Inter, IBM Plex Sans | SIL OFL 1.1 | Free to bundle; keep notices; do not sell the fonts alone; do not reuse the Reserved Font Name on a modified font |

---

## 5. The GPL analysis

Three distinct questions, often conflated.

### Q1 — Does invoking R as a subprocess make the app a derivative work of R?

**No,** on the standard and widely-accepted reading. The FSF's own position is that separate programs communicating at arm's length — pipes, command-line arguments, temporary files — are separate works, not a single combined program. The relationship here is the same one `git` has with `$EDITOR`.

The app spawns `Rscript`, hands it a directory path, and reads back a status line. There is no linking, no shared address space, no shared data structures, no plugin API.

**Design constraints that must be preserved to keep this true:**

- The interface stays text files plus exit codes. No shared memory, no sockets carrying application objects.
- **Never use `rpy2`** or any in-process R embedding. That puts GPL code in your address space and the arm's-length argument evaporates.
- Keep `build_seurat.R` free of application business logic, so it is plainly a small independent script rather than a component of the app.

This is an additional and substantial reason for the sidecar architecture beyond crash isolation.

### Q2 — Is `build_seurat.R` itself derivative of GPL code?

It calls `library(Matrix)`, which is GPL (>= 2). A script that only functions in combination with a GPL library is, on a conservative reading, a derivative work of it.

**Recommendation: dual-license by directory.** Root of the repository MIT; `core/r/` GPL-3. The cost is essentially nil — that directory is ~60 lines of glue with no reuse value — and it removes the argument entirely rather than inviting it.

```
LICENSE                 MIT — everything except core/r/
core/r/LICENSE          GPL-3 — the Seurat bridge script
README.md               explains the split in three sentences
```

Do not contort the design to avoid `Matrix`; sparse matrix construction in R needs it, and dual-licensing a directory is a normal, well-understood pattern.

### Q3 — Does the installer distribute GPL code?

**With decision L-1 as taken: no.** The installer contains Electron, the Python conda environment, and the micromamba binary — all permissive. R, Matrix, and Seurat are fetched from conda-forge by the user's own machine, at the user's explicit instruction, after installation. You are not a distributor of R.

⚠️ **This flips if you ever bundle R in the installer.** Then GPL-3 §6 obligations attach: you must convey the corresponding source, or a written offer valid for three years, for every GPL component you ship. It is achievable (conda-forge sources are public; ship the lockfile plus a source offer) but it is a real, ongoing obligation with real work attached.

**Add this to the spec as a standing note:** the size argument for bundling R is not the only argument, and revisiting L-1 has legal consequences, not just installer-size consequences.

---

## 6. Other legal considerations

### 6.1 GEO data

GEO records are largely US Government works and NCBI states its data may be used without restriction, while disclaiming warranty and noting that submitters may retain rights they assert. Practical rules for this project:

- **Do not commit GEO data to the repository.** Test fixtures are downloaded in CI and cached; only checksums and expected dimensions live in git. This was already the design for size reasons; it is also the right call legally.
- Respect E-utilities rate limits (3 requests/second, 10 with an API key) and send an honest `User-Agent` identifying the application and version. Abusive traffic from a distributed application gets an IP range blocked and would affect every user at once.
- State in the README that the project neither hosts, curates, nor vouches for GEO content, and that users are responsible for the terms attaching to any dataset they download.

### 6.2 Trademarks and naming

Referring to Seurat and Scanpy descriptively — "builds Seurat and Scanpy objects" — is nominative fair use and is fine. To stay clearly inside it:

- Do not use their logos or visual identity.
- Do not incorporate "Seurat" or "Scanpy" into the product name, the repository name, or a domain.
- Do not imply endorsement or affiliation.
- Add a short notice to the README: *"Seurat and Scanpy are the work of their respective authors. This project is independent and is not affiliated with or endorsed by them."*

### 6.3 Fonts

Inter and IBM Plex are SIL OFL 1.1. Bundling and embedding is permitted. Retain the copyright notice and the OFL text in `THIRD_PARTY_NOTICES`, do not sell the fonts on their own, and do not ship a modified font under the Reserved Font Name.

### 6.4 Code signing certificates

Not a licensing matter, but connected to F-1: certificate authorities issue code-signing certificates to a verified legal entity or individual. If the software is institutionally owned, the certificate arguably should be institutional too, and the OTC or IT may have an existing certificate and a policy about its use. Settle ownership before spending money on a certificate.

### 6.5 Data handling disclaimer

The tool handles public data only and performs no encryption beyond TLS, so there is no export-control concern. Given the environments you work in, add an explicit README warning that the application is not validated for PHI or patient-identifiable data and should not be pointed at it. Someone will eventually try.

### 6.6 Academic citation

Add a `CITATION.cff` and a README section directing users to cite Seurat, Scanpy, AnnData, and — importantly — the original study behind whatever GEO accession they processed. A tool that makes data trivially easy to obtain should push people toward crediting the people who generated it.

---

## 7. Recommended license structure

```
Root:            MIT
core/r/:         GPL-3
Copyright line:  to be confirmed with the OTC — likely
                 "The University of Texas M.D. Anderson Cancer Center"
                 rather than an individual
```

MIT over BSD-3 for the root: functionally equivalent for your purposes, shorter, and the more common default in the JavaScript half of the stack, which reduces friction for contributors.

---

## 8. Compliance checklist

**Blocking**

- [ ] OTC software disclosure submitted and written authorization received (F-1)
- [ ] Name changed and cleared — USPTO, GitHub, PyPI, npm, domain, Scholar (F-2)
- [ ] Copyright holder confirmed with the OTC

**Before first public commit**

- [ ] `LICENSE` (MIT) at repository root
- [ ] `core/r/LICENSE` (GPL-3) with the split explained in the README
- [ ] `THIRD_PARTY_NOTICES.md`, generated not hand-written
- [ ] Electron/Chromium/ffmpeg notices included, LGPL relinking statement present
- [ ] Font OFL notices included
- [ ] NCBI/GEO disclaimer in README
- [ ] Trademark disclaimer in README
- [ ] `CITATION.cff`
- [ ] `.gitignore` excludes all data files; CI verifies no fixture data is committed

**In CI, ongoing**

- [ ] `pip-licenses` over the Python environment, failing the build on any license outside an allowlist
- [ ] `license-checker` over the npm tree, same policy
- [ ] `reuse lint` for SPDX header coverage (optional, but it makes the next audit trivial)
- [ ] Regenerate `THIRD_PARTY_NOTICES.md` on every dependency change

---

## 9. Bottom line

The dependency situation is better than expected. Seurat being MIT was the open question and it resolved favourably; the only copyleft in the stack is R and `Matrix`, and the on-demand install architecture you already chose means you never distribute either. A directory-scoped GPL-3 license on the 60-line R bridge closes the remaining gap at no practical cost.

The real obstacles are not technical. They are the institutional clearance and the name — and both are much cheaper to resolve now than after there is a repository, a user base, and a preprint.

---

### Sources

- [MD Anderson Intellectual Property Policy](https://www.mdanderson.org/about-md-anderson/business-legal/legal-and-policy/legal-statements/intellectual-property-policy.html)
- [Seurat license (MIT)](https://satijalab.org/seurat/license) · [seurat/LICENSE on GitHub](https://github.com/satijalab/seurat/blob/master/LICENSE)
- [Matrix DESCRIPTION — GPL (>= 2)](https://github.com/cran/Matrix/blob/master/DESCRIPTION)
- [r-seurat on conda-forge](https://anaconda.org/conda-forge/r-seurat)
- [iOrganBio launches CellForge™](https://www.businesswire.com/news/home/20251029644221/en/iOrganBio-Emerges-from-Stealth-with-$2M-and-Launches-CellForge-the-First-AI-Powered-Cell-Manufacturing-Platform)
- [CellForge: Agentic Design of Virtual Cell Models, arXiv:2508.02276](https://arxiv.org/abs/2508.02276)
- [Programmatic access to GEO — NCBI](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html)
