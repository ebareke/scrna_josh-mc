# E9.5 Mouse Single-Cell Atlas — Spliceosome Mutants (*Eftud2* · *Sf3b4* · Control)

Interactive single-cell RNA-seq dashboard and analysis pipeline comparing two
craniofacial **spliceosomopathy mutants** against a wild-type baseline at
embryonic day **E9.5** in the mouse.

> **🔗 Live dashboard: https://ebareke.github.io/scrna_josh-mc/**

> ⚠️ **Restricted, unpublished research data.** Access is limited to the named
> project collaborators. Do not redistribute. See [`LICENSE`](LICENSE).

---

## What this is

Three conditions are compared, each pooling **two biological replicates**:

| Condition | Role | Description |
|-----------|------|-------------|
| **Control** | baseline / reference | wild-type embryo |
| ***Eftud2*** | mutant | spliceosome (U5 snRNP) mutant |
| ***Sf3b4*** | mutant | spliceosome (SF3b complex) mutant |

Two complementary atlases are analysed and presented:

- **All cell types** — the whole-embryo atlas across all annotated cell types.
- **Neural crest cells (NCC)** — a re-clustered neural-crest compartment with
  NCC-specific sub-types (Cranial, Cardiac, Trunk, Premigratory).

## What the dashboard shows

For **both** atlases, the dashboard visualises the differences between the
three conditions across:

1. **Cell-type composition** — per-condition annotated UMAP/t-SNE, stacked /
   grouped / delta composition, and compositional statistics (pooled χ²,
   per-cell-type Fisher odds ratios, replicate-aware CLR linear-mixed-model
   omnibus + pairwise contrasts).
2. **Differential expression (per cell type)** — interactive volcano plots,
   ranked top-gene tables and matched GO over-representation enrichment, for
   each pairwise contrast (`de/` and `enrichment/de/`).
3. **Differential expression (pseudo-bulk)** — replicate-level DESeq2 contrasts
   and their GO enrichment (`pseudobulk/` and `enrichment/pseudo_bulk/`).
4. **Splicing ratios** — per-cell `unspliced / (spliced + 1)`, summarised
   globally and per cell type per condition.
5. **RNA velocity** — per-condition velocity speed, coherence and latent-time
   summaries.

The three pairwise contrasts throughout are:
`Eftud2 vs Control`, `Sf3b4 vs Control`, and `Eftud2 vs Sf3b4`.

## Repository layout

```
.
├── index.html                 # the interactive dashboard (single page)
├── assets/
│   ├── css/style.css
│   ├── js/app.js              # dashboard logic (Plotly-based)
│   └── js/data.js             # generated data layer (summaries embedded as JSON)
│   └── figures/               # curated, web-optimised figures from the pipeline
├── data/                      # small derived summary tables (CSV) for provenance
├── process_ALL_CELLS.R        # pipeline → whole-embryo atlas
├── process_NEURAL_CREST_CELLS.R  # pipeline → neural-crest atlas
├── build_dashboard.py         # builds data.js and curates figures from pipeline output
├── CITATION.cff
└── LICENSE                    # Restricted Research Use License
```

## Methods (summary)

- **Platform / design.** 10x Genomics scRNA-seq of E9.5 mouse embryos; spliced
  and unspliced layers quantified for RNA velocity. Six samples (3 conditions ×
  2 replicates).
- **Processing.** Per-sample loading, QC, normalisation, integration,
  PCA → UMAP/t-SNE, graph clustering and cell-type annotation (Seurat v3 assay
  model; SingleCellExperiment / scran / scater).
- **Composition statistics.** Pooled χ²; per-cell-type Fisher exact tests with
  BH-adjusted odds ratios; centred-log-ratio (CLR) linear mixed model
  (replicate-aware omnibus LRT + pairwise contrasts).
- **Differential expression.** Cell-type-resolved Wilcoxon tests
  (`Seurat::FindMarkers`); significance at **FDR < 0.05** and **|log₂FC| ≥ 0.25**.
  Replicate-level **pseudo-bulk** with **DESeq2** (FDR < 0.05, |log₂FC| ≥ 1).
- **Functional enrichment.** GO BP/CC/MF over-representation analysis
  (`clusterProfiler`, `org.Mm.eg.db`), split by direction of change.
- **Splicing.** Per-cell ratio = `unspliced / (spliced + 1)`.
- **RNA velocity.** scVelo via `velociraptor`/`reticulate`; per-condition speed,
  coherence and latent-time summaries plus differential-velocity genes.

Full, runnable pipelines are in
[`process_ALL_CELLS.R`](process_ALL_CELLS.R) and
[`process_NEURAL_CREST_CELLS.R`](process_NEURAL_CREST_CELLS.R).

### Key R / Bioconductor packages

`Seurat`, `SingleCellExperiment`, `velociraptor`, `scran`, `scater`, `DESeq2`,
`limma`, `clusterProfiler`, `enrichplot`, `org.Mm.eg.db`, `EnhancedVolcano`,
`ggplot2`, `patchwork`, `pheatmap`, `viridis`.

## Rebuilding the dashboard data

The data layer and curated figures are regenerated from the pipeline output
folders (`ALL_CELLS/`, `NEURAL_CREST_CELLS/`) with:

```bash
python3 build_dashboard.py
```

This parses the DE / pseudo-bulk / enrichment / composition / splicing outputs,
writes `assets/js/data.js`, and copies (down-scaling large PNGs) the figures the
dashboard references into `assets/figures/`.

## Citation

Please cite as described in [`CITATION.cff`](CITATION.cff). Results are
unpublished — see the license terms before any use or disclosure.

## License & contact

**Restricted Research Use** — permitted users only (Prof. Loydie A. Majewska,
Dr. Joshua Moore, Dr. Marie-Claude Beauchamp, Dr. Eric Bareke). See
[`LICENSE`](LICENSE).

Analysis & dashboard by **Eric Bareke** (`eric.bareke@mcgill.ca`),
Department of Human Genetics, McGill University.
