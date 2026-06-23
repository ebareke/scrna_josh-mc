# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdf5r)
  library(SingleCellExperiment)
  library(velociraptor)
  library(scran)
  library(scater)
  library(BiocParallel)
  library(RhpcBLASctl)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(viridis)
  library(pheatmap)
  library(data.table)
  library(ggridges)
  library(matrixTests)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(scales)
  library(DESeq2)
  library(limma)       # for lmFit-based pseudobulk velocity DE
  library(ggsignif)    # significance brackets on ggplots
})

options(Seurat.object.assay.version = "v3")

# =============================================================================
# 0. Configuration
# =============================================================================
IN_BASE        <- "/Users/ebareke/Documents/Projects/scRNA_Eftud2-Sf3b4-Ctrl"
LOOM_DIR       <- file.path(IN_BASE, "LOOM_FILES")
OUT_BASE       <- "/Users/ebareke/Documents/Projects/scRNA_Eftud2-Sf3b4-Ctrl"
CELLRANGER_DIR <- file.path(IN_BASE, "SAMPLES")
OUTDIR         <- file.path(OUT_BASE, "NEURAL_CREST_CELLS")
N_WORKERS      <- 6
BPPARAM        <- MulticoreParam(workers = N_WORKERS)

SAMPLES <- list(
  Control_1_C1_802932 = list(condition = "Control", replicate = "rep1", label = "Ctrl_1"),
  Control_2_C2_802934 = list(condition = "Control", replicate = "rep2", label = "Ctrl_2"),
  Eftud2_1_E1_802933  = list(condition = "Eftud2",  replicate = "rep1", label = "Eftud2_1"),
  Eftud2_2_E2_802935  = list(condition = "Eftud2",  replicate = "rep2", label = "Eftud2_2"),
  Sf3b4_2_S2_802936   = list(condition = "Sf3b4",   replicate = "rep1", label = "Sf3b4_1"),
  Sf3b4_3_S3_802937   = list(condition = "Sf3b4",   replicate = "rep2", label = "Sf3b4_2")
)
CONDITIONS  <- c("Control", "Eftud2", "Sf3b4")
COMPARISONS <- list(c("Eftud2", "Control"), c("Sf3b4", "Control"), c("Eftud2", "Sf3b4"))

FORCE_RERUN           <- FALSE
DIFF_KINETICS_ENABLED <- TRUE
NC_MIN_CELLS          <- 10L   # minimum Neural_crest cells required to run Level-2 analysis

LOOM_FILES <- list(
  per_replicate = setNames(file.path(LOOM_DIR, paste0(names(SAMPLES), ".loom")), names(SAMPLES)),
  per_condition = setNames(file.path(LOOM_DIR, paste0(CONDITIONS, "_condition.loom")), CONDITIONS),
  qc_metrics  = file.path(LOOM_DIR, "AllSamples_combined.loom")
)

# ── AAAS colour palette ───────────────────────────────────────────────────────
AAAS_COLS <- c(
  "#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
  "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"
)
COND_PALETTE <- c(Control = AAAS_COLS[1], Eftud2 = AAAS_COLS[2], Sf3b4 = AAAS_COLS[3])
REP_PALETTE  <- c(rep1    = AAAS_COLS[4], rep2   = AAAS_COLS[5])
# Safe single-key lookup from a named character vector.
.pal_get <- function(pal, key) { v <- pal[key]; if (is.na(v)) NULL else unname(v) }

PLT <- list(
  dpi        = 300,
  base_size  = 14,
  w_single   = 14, h_single = 12,
  w_double   = 20, h_double = 11,
  w_quad     = 24, h_quad   = 20,
  w_triple   = 30, h_triple = 11,
  pt_size    = 0.6,
  label_size = 6
)

AAAS_SPEED_COLS     <- c("grey92", AAAS_COLS[2])
AAAS_COHER_COLS     <- c("grey92", AAAS_COLS[3])
AAAS_LATENT_COLS    <- c(AAAS_COLS[1], "grey95", AAAS_COLS[2])
AAAS_PSEUDO_COLS    <- colorRampPalette(c(AAAS_COLS[5], "grey95", AAAS_COLS[6]))(50)
AAAS_SELFTRANS_COLS <- c("grey92", AAAS_COLS[7])
AAAS_RATIO_COLS     <- colorRampPalette(c(AAAS_COLS[1], "grey95", AAAS_COLS[2]))(50)

# ── Mouse embryo marker genes (E8.25–E8.5 optimized) ─────────────────────────

MOUSE_MARKERS <- list(
  Gut                         = c("Apela","Krt8","Sox17"),
  Caudal_lateral_epiblast     = c("Nkx1-2","Cdx2","Gbx2","T"),
  Neuroectoderm               = c("Hesx1","Hes3"),
  Surface_ectoderm            = c("Foxg1","Grhl2","Grhl3"),
  Extraembryonic_mesoderm     = c("Bmp4","Twist1"),
  Paraxial_mesoderm           = c("Uncx","Hes7","Tcf21","Mesp2"),
  Intermediate_mesoderm       = c("Osr1","Lhx1","Wt1","Foxc2"),
  First_heart_field           = c("Tbx5","Hcn4","Gata4","Sema3c"),
  Allantois                   = c("Hoxa11","Tfap2c"),
  Endothelium                 = c("Kdr","Pecam1","Cdh5"),
  Somatic_mesoderm            = c("Lix1","Msx1","Tbx18"),
  Head_mesoderm               = c("Mesp1","Alx1","Alx3","Pitx1"),
  Forebrain_midbrain          = c("Igfbp2","Otx2","Pcsk1n"),
  Hindbrain                   = c("Crabp1","Fst","Hoxa2"),
  Amniochorionic_mesoderm     = c("Hlx","Postn","Bmp2"),
  Second_heart_field          = c("Isl1","Fgf10","Hand2","Smarcd3"),
  Pre_epidermal_keratinocytes = c("Tfap2b","Trp63","Egfr"),
  Optic_placode               = c("Pou2f1","Hes1"),
  Floor_plate                 = c("Shh","Ntn1","Bmp7"),
  Sclerotome                  = c("Pax1","Sox5","Nkx3-2","Sfrp2"),
  Osteoblast_progenitors      = c("Twist2","Meis2"),
  Skeletal_muscle_progenitors = c("Fap","Sim1","Il1a","Il1b"),
  Limb_mesenchyme_progenitors = c("Cpa2","Wnt5a","Lmx1b"),
  Pharyngeal_arch_epithelium  = c("Eya1","Six1"),
  Olfactory_epithelium        = c("Hes5","Ngn1","Ebf2","Ascl1","Dlx5"),
  Otic_epithelium             = c("Lmx1a","Gata3","Otx1"),
  Gut_and_lung_epithelium     = c("Prox1","Foxp2"),
  Endoderm                    = c("Hex","Foxa2","Afp","Lefty1"),
  Mid_hind_brain              = c("Pax5","En1","Dmbx1"),
  Neuron_progenitor_cells     = c("Mybl1","Prmt8","Gadd45g","Cdkn1c","Btg2"),
  Di_telencephalon            = c("Emx2","Sox5","Sox6","Wnt8b"),
  Neural_crest                = c("GFP", "CRE", "Sox10","Tfap2a","Foxd3","Sox9","Dlx2"),
  Spinal_cord_ventral         = c("Foxb1","Crabp2"),
  Spinal_cord_dorsal          = c("Pth2r","Fabp7","Fzd10","Atoh1","Ptf1a","Wnt1"),
  Apical_ectodermal_ridge     = c("Msx2","Rspo2"),
  Motor_neurons               = c("Mnx1","Isl2","Lhx3","Lhx4")
)

Neural_crest_L2 <- list(
  Premigratory_NC = c(
    "Foxd3","Snai2","Tfap2a","Sox9"
  ),
  Migratory_NC = c(
    "Sox10","Plp1","Ednrb","Foxd3"
  ),
  Cranial_NC = c(
    "Dlx5","Dlx6","Ets1","Tfap2a"
  ),
  Cardiac_NC = c(
    "Hand2","Pax3","Foxc1","Msx1"
  ),
  Trunk_NC = c(
    "Pax3","Sox10","Erbb3","Hes6"
  )
)

.build_ct_palette <- function(cell_types) {
  cell_types <- unique(cell_types)
  n <- length(cell_types)
  base_cols <- c(AAAS_COLS,
                 RColorBrewer::brewer.pal(12, "Paired"),
                 RColorBrewer::brewer.pal(8,  "Dark2"),
                 RColorBrewer::brewer.pal(8,  "Set2"))
  setNames(rep_len(base_cols, n), cell_types)
}

save_plot <- function(plot, path_base, width, height,
                       dpi = PLT$dpi, limitsize = FALSE) {
  ggsave(paste0(path_base, ".png"), plot,
         width = width, height = height, dpi = dpi, limitsize = limitsize)
  suppressMessages(
    ggsave(paste0(path_base, ".pdf"), plot,
           width = width, height = height, limitsize = limitsize)
  )
}

save_pheatmap <- function(ph, path_base, width_px, height_px, dpi = PLT$dpi) {
  png(paste0(path_base, ".png"), width = width_px, height = height_px, res = dpi)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  pdf(paste0(path_base, ".pdf"), width = width_px/dpi, height = height_px/dpi)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Seurat v3 / v5 compatible layer accessor.
# Uses GetAssayData(..., slot = layer) which works reliably under both Seurat v3
# and v5 (the slot= argument is deprecated in v5 but still fully functional).
# Deprecation warnings are suppressed to keep the log clean.
.get_layer <- function(obj, assay = "RNA", layer = "counts") {
  suppressWarnings(GetAssayData(obj, assay = assay, slot = layer))
}

.theme_velo <- function(base_size = PLT$base_size) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(size = base_size * 1.1, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size * 0.9, hjust = 0.5),
      axis.title       = element_text(size = base_size,       face = "bold"),
      axis.text        = element_text(size = base_size * 0.9, face = "bold"),
      legend.title     = element_text(size = base_size,       face = "bold"),
      legend.text      = element_text(size = base_size * 0.85),
      strip.text       = element_text(size = base_size * 0.95, face = "bold"),
      legend.key.size  = unit(0.45, "cm")
    )
}

.auto_dims_v <- function(n_items,
                          base_w = PLT$w_single, base_h = PLT$h_single,
                          orientation = "h",
                          per_item_w  = 0.55, per_item_h = 0.38,
                          max_w = 60, max_h = 50) {
  if (orientation == "h") {
    w <- base_w + n_items * per_item_w; h <- base_h
  } else {
    w <- base_w; h <- base_h + n_items * per_item_h
  }
  list(w = min(w, max_w), h = min(h, max_h))
}

# ── Output directories ────────────────────────────────────────────────────────
for (.d in c("per_replicate", "per_condition", "qc_metrics",
             "differential",
             "de", "enrichment", "enrichment/pseudo_bulk", "enrichment/de",
             "pseudobulk",
             "replicate_concordance", "splicing_ratios",
             "velocity_programs", "composition_stats",
             "quality_control",
             "splicing_efficiency",
             "velocity_pseudotime_corr",
             "celltype_transitions",
             "gene_velocity_coherence",
             "latent_time_ks_factorial",
             "differentially_variable",
             "custom_summary_plots",
             "nc_annotation")) {
  dir.create(file.path(OUTDIR, .d), recursive = TRUE, showWarnings = FALSE)
}
message("Output directory: ", OUTDIR)

# =============================================================================
# 1. Helper functions (loom loading, preprocessing, annotation)
# =============================================================================

.parse_cell_metadata <- function(cell_ids, samples_meta) {
  sample_prefix <- sub(":.*$", "", cell_ids)
  clean_ids <- paste0(
    sample_prefix, ":",
    sub("x$", "", sub("^.*:", "", cell_ids)), "-1"
  )
  condition <- vapply(sample_prefix, function(p) {
    m <- samples_meta[[p]]; if (!is.null(m)) m$condition else "Unknown"
  }, character(1))
  replicate <- vapply(sample_prefix, function(p) {
    m <- samples_meta[[p]]; if (!is.null(m)) m$replicate else "Unknown"
  }, character(1))
  label <- vapply(sample_prefix, function(p) {
    m <- samples_meta[[p]]; if (!is.null(m)) m$label else p
  }, character(1))
  data.frame(orig_id = cell_ids, clean_id = clean_ids,
             sample = sample_prefix, condition = condition,
             replicate = replicate, label = label,
             stringsAsFactors = FALSE)
}

.h5_read_vec <- function(h5, path) as.character(unlist(h5[[path]]$read()))

.read_gene_ids <- function(h5, loom_path) {
  for (key in c("row_attrs/Gene","row_attrs/gene","row_attrs/GeneID","row_attrs/Accession")) {
    ids <- tryCatch(.h5_read_vec(h5, key), error = function(e) NULL)
    if (length(ids) > 0) return(make.unique(ids))
  }
  stop("No gene identifiers found in row_attrs of: ", loom_path)
}

.read_layer <- function(h5, h5path, gene_ids, cell_ids, chunk_size = 2000L) {
  n_g <- length(gene_ids); n_c <- length(cell_ids)
  h5d <- h5[[h5path]]; h5dim <- h5d$dims
  message("      [", h5path, "] HDF5 dims: ", h5dim[1], " x ", h5dim[2],
          "  |  expected genes=", n_g, " cells=", n_c)
  genes_first <- (h5dim[1] == n_g && h5dim[2] == n_c)
  cells_first <- (h5dim[1] == n_c && h5dim[2] == n_g)
  if (!genes_first && !cells_first) stop("Layer '", h5path, "' dims mismatch")
  dim1   <- h5dim[1]
  chunks <- split(seq_len(dim1), ceiling(seq_len(dim1) / chunk_size))
  sparse_chunks <- lapply(chunks, function(idx) {
    i0 <- idx[1L]; i1 <- idx[length(idx)]
    slab <- h5d[i0:i1, ]
    if (cells_first) slab <- t(slab)
    as(slab, "dgCMatrix")
  })
  mat <- if (genes_first) do.call(rbind, sparse_chunks) else do.call(cbind, sparse_chunks)
  rownames(mat) <- gene_ids; colnames(mat) <- cell_ids
  blank <- which(rownames(mat) == "")
  if (length(blank)) rownames(mat)[blank] <- paste0("gene_", blank)
  mat
}

# ---------------------------------------------------------------------------
# Rescue intronless reporter genes (GFP, RFP, …) for velocity loom loading.
#
# Problem 1 (intronless): velocyto classifies ALL reads from single-exon
#           (intronless) genes as "ambiguous".  Spliced count ≈ 0, so Seurat's
#           min.cells filter silently drops them.
# Problem 2 (repeat-masked): if the repeat-masker GTF (-m) covers the
#           transgene contig, velocyto masks every read → ALL layers are zero.
#
# Strategy (two-tier):
#   A. Fold ambiguous counts into spliced for reporters with negligible
#      spliced counts.  Biologically correct for intronless transcripts.
#   B. If all loom layers are zero (repeat-masker wipe-out), fall back to
#      CellRanger's filtered_feature_bc_matrix for the reporter counts and
#      inject them into the spliced matrix.  This preserves gene counts
#      for downstream analysis while the remaining genes come from the loom.
# ---------------------------------------------------------------------------
load_replicate_loom <- function(sample_name) {
  loom_path <- LOOM_FILES$per_replicate[[sample_name]]
  message("  Loading replicate loom: ", basename(loom_path))
  h5 <- hdf5r::H5File$new(loom_path, mode = "r"); on.exit(h5$close_all())
  gene_ids <- .read_gene_ids(h5, loom_path)
  raw_bc   <- .h5_read_vec(h5, "col_attrs/CellID")
  barcodes <- paste0(sub("x$", "", sub("^.*:", "", raw_bc)), "-1")
  spliced   <- .read_layer(h5, "layers/spliced",   gene_ids, barcodes)
  unspliced <- .read_layer(h5, "layers/unspliced", gene_ids, barcodes)
  ambiguous <- .read_layer(h5, "layers/ambiguous", gene_ids, barcodes)
  cr_matrix_dir <- file.path(CELLRANGER_DIR, sample_name,
                              "outs/filtered_feature_bc_matrix")
  bc_path <- file.path(cr_matrix_dir, "barcodes.tsv.gz")
  if (file.exists(bc_path)) {
    cr_bc <- as.character(read.table(gzfile(bc_path), header = FALSE)[[1]])
    keep  <- intersect(barcodes, cr_bc)
    if (length(keep) == 0) stop("No barcode overlap for: ", sample_name)
    spliced   <- spliced[,   keep, drop = FALSE]
    unspliced <- unspliced[, keep, drop = FALSE]
    ambiguous <- ambiguous[, keep, drop = FALSE]
  }
  seu <- CreateSeuratObject(counts = spliced, project = sample_name,
                             min.cells = 3, min.features = 200)
  kg <- rownames(seu); kc <- colnames(seu)
  seu[["unspliced"]] <- CreateAssayObject(counts = unspliced[kg[kg %in% rownames(unspliced)],
                                                               kc, drop = FALSE])
  seu[["ambiguous"]] <- CreateAssayObject(counts = ambiguous[kg[kg %in% rownames(ambiguous)],
                                                               kc, drop = FALSE])
  meta <- SAMPLES[[sample_name]]
  seu$sample    <- sample_name
  seu$condition <- meta$condition
  seu$replicate <- meta$replicate
  seu$label     <- meta$label
  message("    ", ncol(seu), " cells  |  ", nrow(seu), " genes")
  seu
}

load_merged_loom <- function(loom_path, project_name = "merged") {
  message("  Loading merged loom: ", basename(loom_path))
  h5 <- hdf5r::H5File$new(loom_path, mode = "r"); on.exit(h5$close_all())
  gene_ids <- .read_gene_ids(h5, loom_path)
  cell_ids <- .h5_read_vec(h5, "col_attrs/CellID")
  meta_df  <- .parse_cell_metadata(cell_ids, SAMPLES)
  spliced   <- .read_layer(h5, "layers/spliced",   gene_ids, meta_df$clean_id)
  unspliced <- .read_layer(h5, "layers/unspliced", gene_ids, meta_df$clean_id)
  ambiguous <- .read_layer(h5, "layers/ambiguous", gene_ids, meta_df$clean_id)

  seu <- CreateSeuratObject(counts = spliced, project = project_name,
                             min.cells = 1, min.features = 200)
  kg <- rownames(seu); kc <- colnames(seu)
  seu[["unspliced"]] <- CreateAssayObject(counts = unspliced[kg[kg %in% rownames(unspliced)],
                                                               kc, drop = FALSE])
  meta_flt  <- meta_df[match(kc, meta_df$clean_id), ]
  seu$sample    <- meta_flt$sample
  seu$condition <- meta_flt$condition
  seu$replicate <- meta_flt$replicate
  seu$label     <- meta_flt$label
  message("    ", ncol(seu), " cells  |  ", nrow(seu), " genes")
  seu
}

preprocess_seurat <- function(seu, label, use_harmony = FALSE,
                               harmony_vars = "replicate",
                               n_dims = 30, resolution = 0.5) {
  message("  Preprocessing: ", label)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = n_dims, verbose = FALSE)
  if (use_harmony) {
    if (requireNamespace("harmony", quietly = TRUE)) {
      message("    Harmony batch correction on: ", harmony_vars)
      seu <- harmony::RunHarmony(seu, group.by.vars = harmony_vars,
                                  plot_convergence = FALSE, verbose = FALSE)
      reduction <- "harmony"
    } else {
      message("    NOTE: harmony not available; skipping batch correction.")
      reduction <- "pca"
    }
  } else {
    reduction <- "pca"
  }
  seu <- RunUMAP(seu, reduction = reduction, dims = 1:n_dims, verbose = FALSE)
  seu <- RunTSNE(seu, reduction = reduction, dims = 1:n_dims, verbose = FALSE,
                 check_duplicates = FALSE)
  seu <- FindNeighbors(seu, reduction = reduction, dims = 1:n_dims, verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  seu
}

annotate_by_markers <- function(seu, marker_list) {
  message("  Annotating cells by markers...")

  DefaultAssay(seu) <- "RNA"
  mat <- .get_layer(seu, "RNA", "counts")

  if (is.null(mat) || !length(marker_list)) {
    warning("Cannot annotate: empty matrix or marker list")
    seu$all_cell_types <- "Unknown"
    seu$cell_type      <- "Unassigned"
    seu$cluster_majority_type <- "Unassigned"
    seu$cluster_majority_id   <- paste0("C", seu$seurat_clusters, ": Unassigned")
    return(seu)
  }

  feats     <- rownames(mat)
  sets      <- names(marker_list)
  sets_clean <- lapply(marker_list, function(v) {
    unique(intersect(gsub(" \\(.*\\)$", "", v), feats))
  })

  score <- matrix(0L, nrow = ncol(mat), ncol = length(sets),
                  dimnames = list(colnames(mat), sets))

  for (ct in sets) {
    present <- sets_clean[[ct]]
    if (length(present) > 0)
      score[, ct] <- as.integer(Matrix::colSums(mat[present, , drop = FALSE] > 0))
  }

  best <- apply(score, 1, function(x) {
    if (all(x == 0)) return("Unknown")
    w <- which(x == max(x))
    if (length(w) == 1) colnames(score)[w] else "Ambiguous"
  })

  seu$cell_type_marker <- unname(best)
  seu$all_cell_types   <- seu$cell_type_marker

  if ("seurat_clusters" %in% names(seu@meta.data)) {
    cluster_labels <- sapply(levels(seu$seurat_clusters), function(cl) {
      cells <- WhichCells(seu, idents = cl)
      tab   <- table(seu$all_cell_types[cells])
      tab   <- tab[!names(tab) %in% c("Unknown", "Ambiguous")]
      if (length(tab) == 0) return("Unknown")
      names(tab)[which.max(tab)]
    })

    final_labels <- cluster_labels[as.character(seu$seurat_clusters)]
    names(final_labels) <- colnames(seu)
    seu$all_cell_types <- factor(final_labels)
  }

  Idents(seu) <- "all_cell_types"

  message(sprintf("    Annotation complete. %d cell types identified.",
                  length(unique(seu$all_cell_types))))

  # Bridge to downstream column names expected by the rest of the pipeline
  seu$cell_type <- as.character(seu$all_cell_types)
  seu$cell_type[seu$cell_type %in% c("Unknown", "Ambiguous")] <- "Unassigned"

  cluster_ids <- as.character(seu$seurat_clusters)
  seu$cluster_majority_type <- seu$cell_type
  seu$cluster_majority_id   <- paste0("C", cluster_ids, ": ", seu$cell_type)
  message("    Majority-rule types assigned: ", length(unique(seu$cluster_majority_type)), " distinct types")

  seu
}

# =============================================================================
# 1f. Neural crest subset + Level-2 re-annotation
# =============================================================================

.subset_to_nc_and_reannotate <- function(seu, label,
                                          use_harmony  = FALSE,
                                          harmony_vars = "replicate") {
  nc_cells <- colnames(seu)[seu$cell_type == "Neural_crest"]
  n_total  <- ncol(seu)
  n_nc     <- length(nc_cells)

  if (n_nc < NC_MIN_CELLS) {
    message("  WARNING: only ", n_nc, "/", n_total, " Neural_crest cells in ", label,
            " (threshold = ", NC_MIN_CELLS, ") -- skipping Level-2 analysis.")
    return(seu)
  }

  message("  Neural_crest subset: ", n_nc, "/", n_total,
          sprintf(" (%.1f%%)", 100 * n_nc / n_total), " -- ", label)

  seu_nc <- seu[, nc_cells]
  n_dims <- max(5L, min(20L, n_nc - 1L))

  seu_nc <- preprocess_seurat(seu_nc,
                               label        = paste0(label, "_NC"),
                               use_harmony  = use_harmony,
                               harmony_vars = harmony_vars,
                               n_dims       = n_dims,
                               resolution   = 0.4)

  seu_nc <- annotate_by_markers(seu_nc, Neural_crest_L2)

  message("  Level-2 annotation complete for ", label, " (", n_nc, " NC cells).")
  seu_nc
}

# =============================================================================
# 1g. RNA velocity (scvelo via reticulate)
# =============================================================================
.norm_log1p <- function(mat) {
  sf <- Matrix::colSums(mat); sf[sf == 0] <- 1
  target <- median(sf[sf > 0])
  log1p(Matrix::t(Matrix::t(mat) * (target / sf)))
}

run_velociraptor <- function(seu, label) {
  message("  velociraptor (reticulate): ", label)
  scv <- reticulate::import("scvelo",       convert = FALSE)
  np  <- reticulate::import("numpy",        convert = FALSE)
  sp  <- reticulate::import("scipy.sparse", convert = FALSE)
  ad  <- reticulate::import("anndata",      convert = FALSE)
  sc  <- reticulate::import("scanpy",       convert = FALSE)

  n_cells <- ncol(seu)
  message("    Cells available: ", n_cells)
  if (n_cells < 10L) {
    message("    Too few cells (", n_cells, ") for velocity  -  skipping ", label)
    return(list(seu = seu, velo_out = NULL, sce_full = NULL))
  }

  spliced_mat   <- as(.get_layer(seu, "RNA", "counts"), "dgCMatrix")
  unspliced_mat <- as(.get_layer(seu, "unspliced", "counts"), "dgCMatrix")
  common_genes  <- intersect(rownames(spliced_mat), rownames(unspliced_mat))
  if (length(common_genes) == 0) stop("No shared genes for: ", label)
  spliced_mat   <- spliced_mat[common_genes,   colnames(seu), drop = FALSE]
  unspliced_mat <- unspliced_mat[common_genes, colnames(seu), drop = FALSE]

  message("    Normalising spliced / unspliced ...")
  spliced_norm   <- .norm_log1p(spliced_mat)
  unspliced_norm <- .norm_log1p(unspliced_mat)

  hvgs <- intersect(VariableFeatures(seu), common_genes)
  if (length(hvgs) < 100) hvgs <- common_genes
  message("    HVGs for velocity: ", length(hvgs))

  # Adapt PCA dims and neighbor count to the actual data size
  n_pca_available <- if ("pca" %in% names(seu@reductions)) ncol(Embeddings(seu, "pca")) else 30L
  n_pcs_use       <- as.integer(min(30L, n_pca_available, n_cells - 1L))
  n_neighbors_use <- as.integer(min(30L, max(5L, floor(sqrt(n_cells)))))
  message("    Adaptive params: n_pcs=", n_pcs_use,
          ", n_neighbors=", n_neighbors_use, " (", n_cells, " cells)")

  to_scipy_csr <- function(mat) {
    sp$csc_matrix(
      reticulate::tuple(
        np$array(mat@x, dtype = np$float32),
        np$array(mat@i, dtype = np$int32),
        np$array(mat@p, dtype = np$int32)
      ),
      shape = reticulate::tuple(nrow(mat), ncol(mat))
    )$T
  }
  spliced_py   <- to_scipy_csr(spliced_norm)
  unspliced_py <- to_scipy_csr(unspliced_norm)
  adata <- ad$AnnData(
    X      = spliced_py$copy(),
    layers = reticulate::dict(spliced   = spliced_py$copy(),
                               unspliced = unspliced_py$copy())
  )
  adata$obs_names <- reticulate::r_to_py(colnames(spliced_mat))
  adata$var_names <- reticulate::r_to_py(rownames(spliced_mat))

  if ("pca"  %in% names(seu@reductions))
    reticulate::py_set_item(adata$obsm, "X_pca",
      np$array(Embeddings(seu, "pca")[colnames(seu), , drop = FALSE], dtype = np$float32))
  if ("umap" %in% names(seu@reductions))
    reticulate::py_set_item(adata$obsm, "X_umap",
      np$array(Embeddings(seu, "umap")[colnames(seu), , drop = FALSE], dtype = np$float32))
  reticulate::py_set_item(adata$var, "highly_variable",
    reticulate::r_to_py(rownames(spliced_mat) %in% hvgs))

  .step <- function(name, expr) {
    message("    [scvelo] ", name, " ...")
    tryCatch(expr, error = function(e) {
      message("    [scvelo] FAILED at step '", name, "' for ", label)
      stop(name, call. = FALSE)
    })
  }

  velo_out <- tryCatch({
    .step("neighbors",           sc$pp$neighbors(adata, n_pcs = n_pcs_use, n_neighbors = n_neighbors_use, use_rep = "X_pca"))
    .step("moments",             scv$pp$moments(adata, n_neighbors = n_neighbors_use))
    .step("recover_dynamics",    scv$tl$recover_dynamics(adata, n_jobs = N_WORKERS))
    .step("velocity",            scv$tl$velocity(adata, mode = "dynamical"))
    .step("velocity_graph",      scv$tl$velocity_graph(adata))
    .step("velocity_confidence", scv$tl$velocity_confidence(adata))
    .step("velocity_pseudotime", scv$tl$velocity_pseudotime(adata))
    .step("latent_time",         scv$tl$latent_time(adata))
    adata
  }, error = function(e) {
    message("    [scvelo] pipeline aborted at: ", conditionMessage(e)); NULL
  })

  if (!is.null(velo_out)) {
    vel_mat <- tryCatch({
      sp_mod <- reticulate::import("scipy.sparse", convert = FALSE)
      np_mod <- reticulate::import("numpy",        convert = FALSE)
      vel_py  <- velo_out$layers["velocity"]
      vel_csc <- sp_mod$csc_matrix(vel_py)$T$tocsc()
      x_r <- as.numeric(reticulate::py_to_r(np_mod$array(vel_csc$data,    dtype = np_mod$float64)))
      i_r <- as.integer(reticulate::py_to_r(np_mod$array(vel_csc$indices, dtype = np_mod$int32)))
      p_r <- as.integer(reticulate::py_to_r(np_mod$array(vel_csc$indptr,  dtype = np_mod$int32)))
      shp <- reticulate::py_to_r(vel_csc$shape)
      m <- Matrix::sparseMatrix(i = i_r, p = p_r, x = x_r,
                                 dims = c(shp[[1L]], shp[[2L]]), index1 = FALSE)
      rownames(m) <- reticulate::py_to_r(velo_out$var_names$tolist())
      colnames(m) <- reticulate::py_to_r(velo_out$obs_names$tolist())
      m
    }, error = function(e) { message("    velocity extraction failed: ", e$message); NULL })

    if (!is.null(vel_mat)) {
      kc <- intersect(colnames(seu), colnames(vel_mat))
      seu[["velocity"]] <- CreateAssayObject(counts = vel_mat[, kc, drop = FALSE])
    }

    cell_ids <- reticulate::py_to_r(velo_out$obs_names$tolist())
    obs_cols <- tryCatch(reticulate::py_to_r(velo_out$obs$columns$tolist()), error = function(e) character(0))
    for (metric in c("velocity_length","velocity_confidence","velocity_pseudotime",
                     "latent_time","root_cells","end_points","velocity_self_transition")) {
      if (!metric %in% obs_cols) next
      tryCatch({
        vals        <- as.numeric(reticulate::py_to_r(velo_out$obs[[metric]]$to_numpy()))
        names(vals) <- cell_ids
        kc          <- intersect(colnames(seu), names(vals))
        seu[[metric]]             <- NA_real_
        seu@meta.data[kc, metric] <- vals[kc]
      }, error = function(e) message("    metric '", metric, "' skipped: ", e$message))
    }
    n_vg <- tryCatch(sum(reticulate::py_to_r(velo_out$var[["velocity_genes"]]), na.rm = TRUE), error = function(e) NA)
    message("    velocity_genes: ", n_vg)

    tryCatch({
      var_df   <- reticulate::py_to_r(velo_out$var)
      kin_cols <- intersect(c("fit_alpha","fit_beta","fit_gamma",
                               "fit_likelihood","fit_r2"), colnames(var_df))
      if (length(kin_cols) > 0) {
        kin_df <- var_df[, kin_cols, drop = FALSE]
        kin_df$gene <- rownames(kin_df)
        seu@misc$velocity_kinetics <- kin_df
        message("    kinetics params stored: ", paste(kin_cols, collapse = ", "))
      }
    }, error = function(e) message("    kinetics extraction skipped: ", e$message))

    tryCatch({
      if ("Ms" %in% names(reticulate::py_to_r(velo_out$layers))) {
        ms_py   <- velo_out$layers["Ms"]
        sp_mod2 <- reticulate::import("scipy.sparse", convert = FALSE)
        ms_csc  <- sp_mod2$csc_matrix(ms_py)$T$tocsc()
        np_mod2 <- reticulate::import("numpy", convert = FALSE)
        x_r2 <- as.numeric(reticulate::py_to_r(np_mod2$array(ms_csc$data,    dtype = np_mod2$float64)))
        i_r2 <- as.integer(reticulate::py_to_r(np_mod2$array(ms_csc$indices, dtype = np_mod2$int32)))
        p_r2 <- as.integer(reticulate::py_to_r(np_mod2$array(ms_csc$indptr,  dtype = np_mod2$int32)))
        shp2 <- reticulate::py_to_r(ms_csc$shape)
        ms_r <- Matrix::sparseMatrix(i = i_r2, p = p_r2, x = x_r2,
                                      dims = c(shp2[[1L]], shp2[[2L]]), index1 = FALSE)
        rownames(ms_r) <- reticulate::py_to_r(velo_out$var_names$tolist())
        colnames(ms_r) <- reticulate::py_to_r(velo_out$obs_names$tolist())
        seu@misc$smoothed_spliced <- ms_r
        message("    Smoothed spliced (Ms) stored for residual analysis.")
      }
    }, error = function(e) message("    Ms extraction skipped: ", e$message))
  }
  list(seu = seu, velo_out = velo_out, sce_full = NULL)
}

compute_per_cell_velocity_metrics <- function(seu) {
  if (!"velocity" %in% names(seu@assays)) return(seu)
  vel_mat <- .get_layer(seu, "velocity", "counts")

  seu$velocity_magnitude <- as.numeric(Matrix::colMeans(abs(vel_mat), na.rm = TRUE))

  .entropy <- function(v) {
    v <- abs(v[!is.na(v) & v != 0])
    if (length(v) == 0) return(NA_real_)
    p <- v / sum(v)
    -sum(p * log2(p + 1e-12))
  }
  seu$velocity_entropy <- vapply(seq_len(ncol(vel_mat)),
                                  function(j) .entropy(as.numeric(vel_mat[, j])),
                                  numeric(1))

  # Spliced / unspliced fractions per cell
  if ("RNA" %in% names(seu@assays) && "unspliced" %in% names(seu@assays)) {
    spl_counts   <- Matrix::colSums(.get_layer(seu, "RNA", "counts"))
    unspl_counts <- Matrix::colSums(.get_layer(seu, "unspliced", "counts"))
    total        <- spl_counts + unspl_counts + 1e-6
    seu$spliced_fraction   <- spl_counts   / total
    seu$unspliced_fraction <- unspl_counts / total
    seu$splicing_ratio     <- unspl_counts / (spl_counts + 1)
  }

  message("    Per-cell velocity metrics computed (magnitude, entropy, fractions, splicing_ratio).")
  seu
}

compute_replicate_qc <- function(rep_objects, outdir) {
  message("  Computing per-replicate QC metrics ...")
  rows <- lapply(names(rep_objects), function(sn) {
    seu  <- rep_objects[[sn]]
    meta <- SAMPLES[[sn]]
    n_vg <- if ("velocity" %in% names(seu@assays))
               nrow(.get_layer(seu, "velocity", "counts"))
             else NA_integer_
    mean_spl_ratio <- if ("splicing_ratio" %in% colnames(seu@meta.data))
                        mean(seu$splicing_ratio, na.rm = TRUE)
                      else NA_real_
    data.frame(
      sample            = meta$label,
      condition         = meta$condition,
      replicate         = meta$replicate,
      n_cells           = ncol(seu),
      n_genes           = nrow(seu),
      median_lib_size   = median(Matrix::colSums(.get_layer(seu, "RNA", "counts")), na.rm = TRUE),
      median_genes_det  = median(Matrix::colSums(.get_layer(seu, "RNA", "counts") > 0), na.rm = TRUE),
      n_velocity_genes  = n_vg,
      mean_splicing_ratio = mean_spl_ratio,
      stringsAsFactors  = FALSE
    )
  })
  qc_df <- do.call(rbind, rows)
  if (is.null(qc_df) || nrow(qc_df) == 0) {
    message("    No replicate QC data available."); return(invisible(NULL))
  }
  write.csv(qc_df, file.path(outdir, "per_replicate_QC_metrics.csv"), row.names = FALSE)

  # Plot: library size and genes detected
  qc_long <- tidyr::pivot_longer(qc_df,
    cols = c(median_lib_size, median_genes_det, mean_splicing_ratio),
    names_to = "metric", values_to = "value")
  p_qc <- ggplot(qc_long, aes(x = sample, y = value, fill = condition)) +
    geom_col() +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = COND_PALETTE) +
    labs(x = NULL, y = NULL, title = "Per-replicate QC metrics") +
    .theme_velo(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_plot(p_qc, file.path(outdir, "per_replicate_QC_metrics"), 14, 12)
  message("    QC metrics saved.")
  qc_df
}

run_splicing_ratio_analysis <- function(seu_list, group_labels, outdir) {
  message("  Splicing ratio analysis ...")
  ratio_rows <- lapply(seq_along(seu_list), function(i) {
    seu <- seu_list[[i]]; lbl <- group_labels[i]
    if (!"splicing_ratio" %in% colnames(seu@meta.data)) return(NULL)
    df <- data.frame(
      group         = lbl,
      cell_type     = seu$cell_type,
      condition     = if ("condition" %in% colnames(seu@meta.data)) seu$condition else lbl,
      replicate     = if ("replicate" %in% colnames(seu@meta.data)) seu$replicate else "rep1",
      splicing_ratio = seu$splicing_ratio,
      stringsAsFactors = FALSE
    )
    df
  })
  ratio_df <- do.call(rbind, Filter(Negate(is.null), ratio_rows))
  if (is.null(ratio_df) || nrow(ratio_df) == 0) return(invisible(NULL))

  write.csv(ratio_df, file.path(outdir, "splicing_ratios_per_cell.csv"), row.names = FALSE)

  # Violin by condition with pairwise Wilcoxon significance brackets
  conds_present <- sort(unique(ratio_df$condition))
  pairs_sr <- combn(conds_present, 2, simplify = FALSE)

  # Run pairwise Wilcoxon tests
  wtest_sr <- lapply(pairs_sr, function(pr) {
    v1 <- ratio_df$splicing_ratio[ratio_df$condition == pr[1]]
    v2 <- ratio_df$splicing_ratio[ratio_df$condition == pr[2]]
    wt <- wilcox.test(v1, v2, exact = FALSE)
    pv <- wt$p.value
    stars <- dplyr::case_when(
      pv < 0.001 ~ "***",
      pv < 0.01  ~ "**",
      pv < 0.05  ~ "*",
      TRUE       ~ "ns"
    )
    list(comparisons = pr, p.value = pv, label = sprintf("%s\np=%.2e", stars, pv))
  })

  # Stagger bracket heights to avoid overlap
  y_max_sr   <- max(ratio_df$splicing_ratio, na.rm = TRUE)
  y_step_sr  <- (y_max_sr - min(ratio_df$splicing_ratio, na.rm = TRUE)) * 0.12
  y_starts   <- y_max_sr + y_step_sr * seq_along(pairs_sr)

  sig_layers <- mapply(function(wt, y_pos) {
    ggsignif::geom_signif(
      comparisons  = list(wt$comparisons),
      annotations  = wt$label,
      y_position   = y_pos,
      tip_length   = 0.015,
      textsize     = 3.2,
      vjust        = 0.3,
      parse        = FALSE
    )
  }, wtest_sr, y_starts, SIMPLIFY = FALSE)

  y_top_sr <- max(y_starts) + y_step_sr * 0.8

  p_cond <- ggplot(ratio_df, aes(x = condition, y = splicing_ratio, fill = condition)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.07, outlier.size = 0.4, fill = "white") +
    sig_layers +
    scale_fill_manual(values = COND_PALETTE) +
    coord_cartesian(ylim = c(NA, y_top_sr)) +
    labs(y = "Splicing ratio (unspliced / (spliced+1))", x = NULL,
         title = "Splicing ratio by condition") +
    .theme_velo() + theme(legend.position = "none")
  save_plot(p_cond, file.path(outdir, "splicing_ratio_by_condition"), 12, 10)

  # By cell type + condition
  ct_ratio <- ratio_df %>%
    dplyr::filter(cell_type != "Unassigned") %>%
    dplyr::group_by(cell_type, condition) %>%
    dplyr::summarise(mean_ratio = mean(splicing_ratio, na.rm = TRUE),
                     sd_ratio   = sd(splicing_ratio,   na.rm = TRUE),
                     n          = dplyr::n(), .groups = "drop")
  write.csv(ct_ratio, file.path(outdir, "splicing_ratio_by_celltype_condition.csv"), row.names = FALSE)

  n_ct <- length(unique(ct_ratio$cell_type))
  p_ct <- ggplot(ct_ratio, aes(x = reorder(cell_type, mean_ratio), y = mean_ratio, fill = condition)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio),
                  position = position_dodge(0.9), width = 0.3) +
    scale_fill_manual(values = COND_PALETTE) +
    coord_flip() +
    labs(x = NULL, y = "Mean splicing ratio", title = "Splicing ratio by cell type & condition") +
    .theme_velo(base_size = 11)
  dims <- .auto_dims_v(n_ct, base_h = 8, orientation = "v")
  save_plot(p_ct, file.path(outdir, "splicing_ratio_by_celltype"), 16, dims$h, limitsize = FALSE)
  message("    Splicing ratio analysis complete.")
}

run_replicate_concordance <- function(rep_objects, outdir) {
  message("\n--- Replicate concordance analysis ---")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sample_names <- names(rep_objects)
  pb_list <- lapply(sample_names, function(sn) {
    seu <- rep_objects[[sn]]
    Matrix::rowMeans(.get_layer(seu, "RNA", "counts"))
  })
  names(pb_list) <- vapply(sample_names, function(s) SAMPLES[[s]]$label, character(1))

  common_genes <- Reduce(intersect, lapply(pb_list, names))
  if (length(common_genes) > 100) {
    pb_mat <- do.call(cbind, lapply(pb_list, function(v) log1p(v[common_genes])))
    cor_mat <- cor(pb_mat, method = "pearson")

    # Heatmap
    ph_cor <- pheatmap(cor_mat,
                       color        = colorRampPalette(c("white", AAAS_COLS[1]))(100),
                       breaks       = seq(0.5, 1, length.out = 101),
                       cluster_rows = TRUE, cluster_cols = TRUE,
                       main         = "Pseudobulk expression correlation (Pearson)",
                       silent       = TRUE)
    save_pheatmap(ph_cor, file.path(outdir, "pseudobulk_expression_correlation"), 1800, 1600)

    # Spearman
    cor_spear <- cor(pb_mat, method = "spearman")
    write.csv(as.data.frame(cor_spear),
              file.path(outdir, "pseudobulk_spearman_correlation.csv"))
    write.csv(as.data.frame(cor_mat),
              file.path(outdir, "pseudobulk_pearson_correlation.csv"))
    message("    Pseudobulk expression correlation saved.")

    cond_rep_pairs <- list(
      c("Ctrl_1",   "Ctrl_2"),
      c("Eftud2_1", "Eftud2_2"),
      c("Sf3b4_1",  "Sf3b4_2")
    )
    ba_plots <- lapply(cond_rep_pairs, function(pair) {
      r1 <- pair[1]; r2 <- pair[2]
      if (!r1 %in% colnames(pb_mat) || !r2 %in% colnames(pb_mat)) return(NULL)
      df_ba <- data.frame(
        mean = (pb_mat[, r1] + pb_mat[, r2]) / 2,
        diff =  pb_mat[, r1] - pb_mat[, r2]
      )
      mn   <- mean(df_ba$diff, na.rm = TRUE)
      sd15 <- 1.96 * sd(df_ba$diff, na.rm = TRUE)
      ggplot(df_ba, aes(mean, diff)) +
        geom_point(alpha = 0.25, size = 0.6, colour = "grey40") +
        geom_hline(yintercept = mn,          colour = AAAS_COLS[2], linetype = "solid") +
        geom_hline(yintercept = mn + sd15,   colour = AAAS_COLS[1], linetype = "dashed") +
        geom_hline(yintercept = mn - sd15,   colour = AAAS_COLS[1], linetype = "dashed") +
        labs(x = "Mean log1p expression", y = "Difference (rep1 - rep2)",
             title = paste("Bland-Altman:", r1, "vs", r2)) +
        .theme_velo(base_size = 12)
    })
    ba_plots <- Filter(Negate(is.null), ba_plots)
    if (length(ba_plots) > 0) {
      p_ba <- wrap_plots(ba_plots, ncol = min(3, length(ba_plots)))
      save_plot(p_ba, file.path(outdir, "bland_altman_replicates"),
                PLT$w_triple, PLT$h_single, limitsize = FALSE)
      message("    Bland-Altman plots saved.")
    }
  }

  all_meta <- do.call(rbind, lapply(sample_names, function(sn) {
    seu <- rep_objects[[sn]]
    data.frame(
      label     = SAMPLES[[sn]]$label,
      condition = seu$condition,
      replicate = seu$replicate,
      stringsAsFactors = FALSE
    )
  }))

  pca_rows <- lapply(sample_names, function(sn) {
    seu <- rep_objects[[sn]]
    if (!"pca" %in% names(seu@reductions)) return(NULL)
    emb <- Embeddings(seu, "pca")[, 1:2, drop = FALSE]
    data.frame(
      PC1 = mean(emb[,1]), PC2 = mean(emb[,2]),
      label     = SAMPLES[[sn]]$label,
      condition = SAMPLES[[sn]]$condition,
      replicate = SAMPLES[[sn]]$replicate,
      stringsAsFactors = FALSE
    )
  })
  pca_df <- do.call(rbind, Filter(Negate(is.null), pca_rows))
  if (!is.null(pca_df) && nrow(pca_df) > 1) {
    p_pca_rep <- ggplot(pca_df, aes(PC1, PC2, colour = condition, shape = replicate, label = label)) +
      geom_point(size = 6) +
      ggrepel::geom_text_repel(size = 4, fontface = "bold") +
      scale_colour_manual(values = COND_PALETTE) +
      labs(title = "PCA centroids coloured by replicate/condition (bias detection)") +
      .theme_velo()
    save_plot(p_pca_rep, file.path(outdir, "PCA_replicate_bias_detection"), 14, 10)
    message("    PCA replicate bias plot saved.")
  }

  vel_pairs <- list(
    list(s1 = "Control_1_C1_802932", s2 = "Control_2_C2_802934", lbl = "Control_rep1_vs_rep2"),
    list(s1 = "Eftud2_1_E1_802933",  s2 = "Eftud2_2_E2_802935",  lbl = "Eftud2_rep1_vs_rep2"),
    list(s1 = "Sf3b4_2_S2_802936",   s2 = "Sf3b4_3_S3_802937",   lbl = "Sf3b4_rep1_vs_rep2")
  )
  vel_cors <- lapply(vel_pairs, function(p) {
    s1 <- rep_objects[[p$s1]]; s2 <- rep_objects[[p$s2]]
    if (is.null(s1) || is.null(s2)) return(NULL)
    if (!"velocity" %in% names(s1@assays) || !"velocity" %in% names(s2@assays)) return(NULL)
    v1 <- Matrix::rowMeans(.get_layer(s1, "velocity", "counts"), na.rm = TRUE)
    v2 <- Matrix::rowMeans(.get_layer(s2, "velocity", "counts"), na.rm = TRUE)
    cg <- intersect(names(v1), names(v2))
    if (length(cg) < 50) return(NULL)
    r_pearson <- cor(v1[cg], v2[cg], method = "pearson",  use = "complete.obs")
    r_spear   <- cor(v1[cg], v2[cg], method = "spearman", use = "complete.obs")
    data.frame(comparison = p$lbl, pearson = r_pearson, spearman = r_spear,
               n_genes = length(cg))
  })
  vel_cor_df <- do.call(rbind, Filter(Negate(is.null), vel_cors))
  if (!is.null(vel_cor_df) && nrow(vel_cor_df) > 0) {
    write.csv(vel_cor_df, file.path(outdir, "velocity_replicate_correlation.csv"), row.names = FALSE)
    message("    Velocity replicate correlation saved.")
  }

  if (requireNamespace("variancePartition", quietly = TRUE)) {
    message("    variancePartition available  -  run separately on pseudobulk matrix if needed.")
  } else {
    message("    NOTE: Install variancePartition for formal variance decomposition.")
  }

  message("  Replicate concordance analysis complete.")
}

run_latent_time_ks_tests <- function(seu_q, seu_r, query_label, ref_label, outdir) {
  message("  Latent time KS tests: ", query_label, " vs ", ref_label)
  if (!"latent_time" %in% colnames(seu_q@meta.data) ||
      !"latent_time" %in% colnames(seu_r@meta.data)) {
    message("    latent_time missing in one/both objects  -  skipping"); return(invisible(NULL))
  }
  if (!"cell_type" %in% colnames(seu_q@meta.data)) return(invisible(NULL))

  cts_q <- setdiff(unique(seu_q$cell_type), "Unassigned")
  cts_r <- setdiff(unique(seu_r$cell_type), "Unassigned")
  cts   <- sort(intersect(cts_q, cts_r))

  ks_rows <- lapply(cts, function(ct) {
    lt_q <- seu_q$latent_time[seu_q$cell_type == ct]
    lt_r <- seu_r$latent_time[seu_r$cell_type == ct]
    lt_q <- lt_q[!is.na(lt_q)]; lt_r <- lt_r[!is.na(lt_r)]
    if (length(lt_q) < 5 || length(lt_r) < 5) return(NULL)
    tst <- ks.test(lt_q, lt_r)
    data.frame(cell_type = ct,
               n_query   = length(lt_q), n_ref = length(lt_r),
               D_stat    = tst$statistic, pvalue = tst$p.value,
               mean_lt_query = mean(lt_q), mean_lt_ref = mean(lt_r),
               delta_mean_lt = mean(lt_q) - mean(lt_r),
               stringsAsFactors = FALSE)
  })
  ks_df <- do.call(rbind, Filter(Negate(is.null), ks_rows))
  if (is.null(ks_df) || nrow(ks_df) == 0) return(invisible(NULL))
  ks_df$padj <- p.adjust(ks_df$pvalue, method = "BH")
  ks_df <- ks_df[order(ks_df$pvalue), ]
  lbl <- paste0(query_label, "_vs_", ref_label)
  write.csv(ks_df, file.path(outdir, paste0(lbl, "_latent_time_KS_tests.csv")), row.names = FALSE)

  p_lt <- ggplot(ks_df, aes(x = reorder(cell_type, delta_mean_lt), y = delta_mean_lt,
                              fill = padj < 0.05)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = AAAS_COLS[2], `FALSE` = "grey70"),
                      labels = c(`TRUE` = "FDR<0.05", `FALSE` = "n.s.")) +
    coord_flip() +
    labs(x = NULL, y = paste0("Delta  Mean latent time (", query_label, " - ", ref_label, ")"),
         title = paste("Latent time shift by cell type:", lbl), fill = NULL) +
    .theme_velo(base_size = 12)
  dims <- .auto_dims_v(nrow(ks_df), base_h = 7, orientation = "v")
  save_plot(p_lt, file.path(outdir, paste0(lbl, "_latent_time_shift")),
            14, dims$h, limitsize = FALSE)
  message("    Latent time KS tests saved: ", nrow(ks_df), " cell types")
  ks_df
}

run_velocity_gene_programs <- function(seu, label, outdir) {
  message("  Velocity gene program analysis: ", label)
  if (!"velocity" %in% names(seu@assays)) {
    message("    No velocity assay  -  skipping"); return(invisible(NULL))
  }
  vel_mat <- .get_layer(seu, "velocity", "counts")
  # Keep only genes with non-trivial velocity in ≥10 cells
  nz <- Matrix::rowSums(vel_mat != 0, na.rm = TRUE)
  vel_active <- vel_mat[nz >= 10, , drop = FALSE]
  if (nrow(vel_active) < 20) {
    message("    Too few active velocity genes  -  skipping"); return(invisible(NULL))
  }
  MAX_CELLS_KM <- 5000L
  use_cells <- if (ncol(vel_active) > MAX_CELLS_KM)
                 sample(seq_len(ncol(vel_active)), MAX_CELLS_KM) else seq_len(ncol(vel_active))
  vel_dense <- as.matrix(vel_active[, use_cells, drop = FALSE])

  # Cluster velocity genes by pattern (k-means on gene  x  cell velocity profiles)
  set.seed(42)
  n_clust <- min(6L, floor(nrow(vel_dense) / 20))
  if (n_clust < 2) { message("    Too few genes for clustering"); return(invisible(NULL)) }
  km <- tryCatch(kmeans(vel_dense, centers = n_clust, nstart = 10, iter.max = 50),
                 error = function(e) NULL)
  if (is.null(km)) return(invisible(NULL))

  gene_clusters <- data.frame(gene    = rownames(vel_active),
                               cluster = km$cluster,
                               stringsAsFactors = FALSE)
  write.csv(gene_clusters, file.path(outdir, paste0(label, "_velocity_gene_clusters.csv")),
            row.names = FALSE)

  # Mean velocity per cluster per cell-type  -  heatmap
  if ("cell_type" %in% colnames(seu@meta.data)) {
    ct_vals_vp <- seu$cell_type[match(colnames(vel_dense), colnames(seu))]
    ct_mean <- lapply(seq_len(n_clust), function(k) {
      g_k  <- gene_clusters$gene[gene_clusters$cluster == k]
      vm_k <- colMeans(vel_dense[g_k, , drop = FALSE], na.rm = TRUE)
      tapply(vm_k, ct_vals_vp, mean, na.rm = TRUE)
    })
    ct_mean_mat <- do.call(rbind, ct_mean)
    rownames(ct_mean_mat) <- paste0("Cluster_", seq_len(n_clust))
    # Scale pixel dims based on matrix size to avoid oversized text
    n_r_vp <- nrow(ct_mean_mat); n_c_vp <- ncol(ct_mean_mat)
    h_vp   <- max(900,  n_r_vp * 80);  w_vp <- max(1800, n_c_vp * 55)
    fs_r_vp <- max(5, min(10, floor(h_vp / n_r_vp / PLT$dpi * 72 * 0.55)))
    fs_c_vp <- max(4, min(8,  floor(w_vp / n_c_vp / PLT$dpi * 72 * 0.45)))
    ph_vp <- tryCatch(
      pheatmap(ct_mean_mat,
               color        = colorRampPalette(c(AAAS_COLS[1],"grey95",AAAS_COLS[2]))(100),
               cluster_rows = nrow(ct_mean_mat) >= 2, cluster_cols = ncol(ct_mean_mat) >= 2,
               fontsize_row = fs_r_vp, fontsize_col = fs_c_vp,
               main         = paste0("Velocity gene programs  x  cell types  -  ", label),
               silent       = TRUE),
      error = function(e) { message("    pheatmap failed: ", e$message); NULL })
    if (!is.null(ph_vp))
      save_pheatmap(ph_vp, file.path(outdir, paste0(label, "_velocity_programs_heatmap")), w_vp, h_vp)

    # ── Per-cluster summary bar plots ─────────────────────────────────────────
    # For each cluster: n_genes, mean velocity per cell type, top cell types
    cluster_summaries <- lapply(seq_len(n_clust), function(k) {
      g_k   <- gene_clusters$gene[gene_clusters$cluster == k]
      vm_k  <- colMeans(vel_dense[g_k, , drop = FALSE], na.rm = TRUE)
      ct_vm <- tapply(vm_k, ct_vals_vp, mean, na.rm = TRUE)
      data.frame(cell_type  = names(ct_vm),
                 mean_vel   = as.numeric(ct_vm),
                 cluster    = paste0("Cluster_", k),
                 n_genes    = length(g_k),
                 stringsAsFactors = FALSE)
    })
    clust_df <- do.call(rbind, cluster_summaries)
    clust_df  <- clust_df[!is.na(clust_df$mean_vel), ]

    if (nrow(clust_df) > 0) {
      # Combined faceted bar plot: mean velocity per cell type, facet by cluster
      n_cls <- n_clust
      p_cbar <- ggplot(clust_df, aes(y = reorder(cell_type, mean_vel),
                                      x = mean_vel, fill = mean_vel > 0)) +
        geom_col(alpha = 0.85) +
        scale_fill_manual(values = c(`TRUE` = AAAS_COLS[2], `FALSE` = AAAS_COLS[1]),
                          guide = "none") +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
        facet_wrap(~ paste0(cluster, " (", n_genes, " genes)"),
                   scales = "free_y", ncol = min(3L, n_cls)) +
        labs(x = "Mean velocity", y = NULL,
             title    = paste0(label, "  -  Velocity cluster metrics by cell type"),
             subtitle = "Positive = upregulated (induction); negative = downregulated (repression)") +
        .theme_velo(base_size = 10) +
        theme(strip.text = element_text(size = 9, face = "bold"),
              axis.text.y = element_text(size = 7))
      h_cbar <- max(8, length(unique(clust_df$cell_type)) * 0.35 + 2) *
                  ceiling(n_cls / 3)
      save_plot(p_cbar,
                file.path(outdir, paste0(label, "_velocity_cluster_metrics_by_celltype")),
                PLT$w_triple * 0.85, h_cbar, limitsize = FALSE)

      # Per-cluster individual ranked bar (top 10 cell types)
      for (k in seq_len(n_clust)) {
        df_k <- clust_df[clust_df$cluster == paste0("Cluster_", k), ]
        df_k <- df_k[order(abs(df_k$mean_vel), decreasing = TRUE), ]
        df_k <- head(df_k, 15)
        if (nrow(df_k) < 2) next
        df_k$cell_type <- factor(df_k$cell_type,
                                  levels = df_k$cell_type[order(df_k$mean_vel)])
        p_k <- ggplot(df_k, aes(y = cell_type, x = mean_vel, fill = mean_vel > 0)) +
          geom_col(alpha = 0.88) +
          scale_fill_manual(values = c(`TRUE` = AAAS_COLS[2], `FALSE` = AAAS_COLS[1]),
                            guide = "none") +
          geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
          geom_text(aes(label = cell_type,
                        x = ifelse(mean_vel >= 0, 0.001, -0.001),
                        hjust = ifelse(mean_vel >= 0, 0, 1)),
                    size = 3, colour = "white", fontface = "bold") +
          scale_y_discrete(labels = NULL) +
          labs(x = "Mean velocity", y = NULL,
               title    = paste0(label, "  -  Velocity Cluster_", k,
                                  " (", df_k$n_genes[1], " genes)"),
               subtitle = "Top cell types by mean velocity magnitude") +
          .theme_velo(base_size = 12)
        save_plot(p_k,
                  file.path(outdir, paste0(label, "_velocity_cluster", k, "_celltype_bar")),
                  12, max(5, nrow(df_k) * 0.4), limitsize = FALSE)
      }
    }
  }

  # Pathway enrichment per velocity gene cluster
  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    all_genes_ent <- tryCatch(
      suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
        keys = rownames(vel_dense), column = "ENTREZID",
        keytype = "SYMBOL", multiVals = "first")),
      error = function(e) NULL)
    if (!is.null(all_genes_ent)) {
      all_genes_ent <- na.omit(all_genes_ent)
      for (k in seq_len(n_clust)) {
        g_k   <- gene_clusters$gene[gene_clusters$cluster == k]
        g_ent <- tryCatch(
          suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
            keys = g_k, column = "ENTREZID",
            keytype = "SYMBOL", multiVals = "first")),
          error = function(e) NULL)
        if (is.null(g_ent) || length(na.omit(g_ent)) < 5) next
        tryCatch({
          ego <- clusterProfiler::enrichGO(
            gene = na.omit(g_ent), universe = all_genes_ent,
            OrgDb = org.Mm.eg.db::org.Mm.eg.db,
            keyType = "ENTREZID", ont = "BP",
            pAdjustMethod = "BH", pvalueCutoff = 0.1, readable = TRUE)
          ego_df <- as.data.frame(ego)
          if (nrow(ego_df) > 0) {
            write.csv(ego_df,
                      file.path(outdir, paste0(label, "_velocity_cluster", k, "_GO_BP.csv")),
                      row.names = FALSE)
            # Save barplot of top GO terms
            tryCatch({
              p_go <- enrichplot::barplot(ego, showCategory = 15) +
                ggtitle(paste0(label, "  -  Velocity Cluster_", k, " GO-BP")) +
                .theme_velo(base_size = 10)
              save_plot(p_go,
                        file.path(outdir, paste0(label, "_velocity_cluster", k, "_GO_BP_barplot")),
                        14, 8)
            }, error = function(e) NULL)
          }
        }, error = function(e) NULL)
      }
    }
  }
  message("    Velocity gene program analysis complete.")
  gene_clusters
}

run_velocity_residuals <- function(seu, label, outdir) {
  message("  Velocity residuals: ", label)
  if (is.null(seu@misc$smoothed_spliced)) {
    message("    Smoothed spliced (Ms) not stored  -  skipping"); return(invisible(NULL))
  }
  if (!"velocity" %in% names(seu@assays)) return(invisible(NULL))

  Ms  <- seu@misc$smoothed_spliced
  vel <- .get_layer(seu, "velocity", "counts")
  cg  <- intersect(rownames(Ms), rownames(vel))
  cc  <- intersect(colnames(Ms), colnames(vel))
  if (length(cg) < 10 || length(cc) < 10) return(invisible(NULL))

  vel_sp   <- vel[cg, cc, drop = FALSE]
  mean_vel <- as.numeric(Matrix::rowMeans(vel_sp, na.rm = TRUE))
  mean_vsq <- as.numeric(Matrix::rowMeans(vel_sp * vel_sp, na.rm = TRUE))
  sd_vel   <- sqrt(pmax(0, mean_vsq - mean_vel ^ 2))

  resid_stats <- data.frame(
    gene        = cg,
    mean_resid  = -mean_vel,
    sd_resid    = sd_vel,
    stringsAsFactors = FALSE
  )
  resid_stats <- resid_stats[order(abs(resid_stats$mean_resid), decreasing = TRUE), ]
  write.csv(resid_stats,
            file.path(outdir, paste0(label, "_velocity_residuals.csv")), row.names = FALSE)

  # Plot top 20 highest-residual genes
  top20 <- head(resid_stats, 20)
  p_res <- ggplot(top20, aes(x = reorder(gene, abs(mean_resid)), y = mean_resid,
                               fill = mean_resid > 0)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_resid - sd_resid, ymax = mean_resid + sd_resid),
                  width = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = AAAS_COLS[2], `FALSE` = AAAS_COLS[1]),
                      guide = "none") +
    labs(x = NULL, y = "Mean residual (observed - predicted)",
         title = paste("Top velocity residual genes  - ", label)) +
    .theme_velo(base_size = 12)
  save_plot(p_res, file.path(outdir, paste0(label, "_velocity_residuals_top20")), 12, 9)
  message("    Velocity residuals saved.")
}

compute_splicing_efficiency_index <- function(seu, label, outdir) {
  message("  Splicing efficiency index (gene  x  cell type): ", label)
  if (!"RNA"       %in% names(seu@assays) ||
      !"unspliced" %in% names(seu@assays)) {
    message("    Missing RNA or unspliced assay  -  skipping"); return(invisible(NULL))
  }
  if (!"cell_type" %in% colnames(seu@meta.data)) {
    message("    No cell_type column  -  skipping"); return(invisible(NULL))
  }
  spl   <- .get_layer(seu, "RNA", "counts")
  unspl <- .get_layer(seu, "unspliced", "counts")
  cg    <- intersect(rownames(spl), rownames(unspl))
  cc    <- intersect(colnames(spl), colnames(unspl))
  spl   <- spl[cg, cc, drop = FALSE]
  unspl <- unspl[cg, cc, drop = FALSE]
  ct_vals    <- seu@meta.data[cc, "cell_type"]
  cell_types <- setdiff(unique(ct_vals), "Unassigned")
  if (length(cell_types) == 0) return(invisible(NULL))

  sei_mat <- do.call(cbind, lapply(cell_types, function(ct) {
    mask <- ct_vals == ct
    if (sum(mask) < 5) return(rep(NA_real_, length(cg)))
    spl_m   <- Matrix::rowMeans(spl[,   mask, drop = FALSE], na.rm = TRUE)
    unspl_m <- Matrix::rowMeans(unspl[, mask, drop = FALSE], na.rm = TRUE)
    unspl_m / (unspl_m + spl_m + 1e-6)
  }))
  rownames(sei_mat) <- cg
  colnames(sei_mat) <- cell_types

  write.csv(as.data.frame(sei_mat),
            file.path(outdir, paste0(label, "_splicing_efficiency_index.csv")))

  # Heatmap: top 80 most variable genes across cell types
  gv    <- apply(sei_mat, 1, var, na.rm = TRUE)
  top_g <- names(sort(gv, decreasing = TRUE))[seq_len(min(80, sum(!is.na(gv))))]
  sei_top <- sei_mat[top_g, , drop = FALSE]
  sei_top[is.na(sei_top)] <- 0.5   # impute NA as mid-value for display

  if (nrow(sei_top) >= 3 && ncol(sei_top) >= 2) {
    h_px <- max(1200, nrow(sei_top) * 16); w_px <- max(800, ncol(sei_top) * 60)
    fs_r <- max(4, min(8, floor(h_px / nrow(sei_top) / PLT$dpi * 72 * 0.50)))
    fs_c <- max(5, min(9, floor(w_px / ncol(sei_top) / PLT$dpi * 72 * 0.50)))
    ph_sei <- tryCatch(
      pheatmap(sei_top,
               color        = colorRampPalette(c(AAAS_COLS[1], "grey95", AAAS_COLS[2]))(100),
               breaks       = seq(0, 1, length.out = 101),
               cluster_rows = nrow(sei_top) >= 2, cluster_cols = ncol(sei_top) >= 2,
               fontsize      = 7,
               fontsize_row = fs_r, fontsize_col = fs_c,
               main         = paste0("Splicing efficiency index ",
                                     "(unspliced / (unspliced+spliced))  -  top variable genes\n", label),
               silent       = TRUE),
      error = function(e) { message("    SEI pheatmap failed: ", e$message); NULL })
    if (!is.null(ph_sei))
      save_pheatmap(ph_sei, file.path(outdir, paste0(label, "_splicing_efficiency_heatmap")), w_px, h_px * 2)
    message("    SEI heatmap saved: ", nrow(sei_top), " genes  x  ", ncol(sei_top), " cell types.")
  }
  invisible(sei_mat)
}

compute_velocity_pseudotime_correlation <- function(seu, label, outdir) {
  message("  Velocity-pseudotime correlation: ", label)
  pt_col <- if ("latent_time"         %in% colnames(seu@meta.data)) "latent_time" else
            if ("velocity_pseudotime" %in% colnames(seu@meta.data)) "velocity_pseudotime" else NULL
  if (is.null(pt_col)) {
    message("    No pseudotime column  -  skipping"); return(invisible(NULL))
  }
  if (!"velocity" %in% names(seu@assays)) {
    message("    No velocity assay  -  skipping"); return(invisible(NULL))
  }
  pt_vals <- seu@meta.data[[pt_col]]
  valid   <- !is.na(pt_vals)
  if (sum(valid) < 20) { message("    Too few cells with pseudotime  -  skipping"); return(invisible(NULL)) }

  vel_mat <- .get_layer(seu, "velocity", "counts")
  vel_mat <- vel_mat[, valid, drop = FALSE]
  pt_sub  <- pt_vals[valid]
  nz      <- Matrix::rowSums(vel_mat != 0, na.rm = TRUE)
  genes_u <- names(nz)[nz >= 10]
  if (length(genes_u) < 5) { message("    Too few velocity genes  -  skipping"); return(invisible(NULL)) }

  vel_sp_u <- vel_mat[genes_u, , drop = FALSE]
  cor_v <- vapply(genes_u, function(g) {
    v  <- as.numeric(vel_sp_u[g, ])
    ok <- !is.na(v)
    if (sum(ok) < 10) return(NA_real_)
    cor(v[ok], pt_sub[ok], method = "pearson", use = "complete.obs")
  }, numeric(1))

  cor_df <- data.frame(gene = genes_u, pearson_r = cor_v, stringsAsFactors = FALSE)
  cor_df <- cor_df[!is.na(cor_df$pearson_r), ]
  cor_df <- cor_df[order(abs(cor_df$pearson_r), decreasing = TRUE), ]
  write.csv(cor_df, file.path(outdir, paste0(label, "_velocity_pseudotime_correlation.csv")),
            row.names = FALSE)

  # Bar chart: top +/-15 genes
  top_pos <- head(cor_df[cor_df$pearson_r > 0, ], 15)
  top_neg <- head(cor_df[cor_df$pearson_r < 0, ], 15)
  top_df  <- rbind(top_pos, top_neg); top_df <- top_df[!is.na(top_df$pearson_r), ]
  if (nrow(top_df) > 1) {
    top_df$gene <- factor(top_df$gene, levels = top_df$gene[order(top_df$pearson_r)])
    p_bar <- ggplot(top_df, aes(y = gene, x = pearson_r, fill = pearson_r > 0)) +
      geom_col(alpha = 0.88) +
      scale_fill_manual(values = c(`TRUE` = AAAS_COLS[2], `FALSE` = AAAS_COLS[1]), guide = "none") +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      labs(x = paste0("Pearson r  (velocity ~ ", pt_col, ")"), y = NULL,
           title = paste0(label, "  -  Velocity-pseudotime correlation (top genes)")) +
      .theme_velo(base_size = 12)
    save_plot(p_bar, file.path(outdir, paste0(label, "_velocity_pseudotime_cor_top_genes")),
              12, max(6, nrow(top_df) * 0.40), limitsize = FALSE)
  }

  # Distribution histogram across all velocity genes
  p_hist <- ggplot(cor_df, aes(x = pearson_r)) +
    geom_histogram(bins = 40, fill = .pal_get(COND_PALETTE, label) %||% AAAS_COLS[1],
                   alpha = 0.75, colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = paste0("Pearson r  (velocity ~ ", pt_col, ")"), y = "Gene count",
         title    = paste0(label, "  -  Velocity-pseudotime correlation distribution"),
         subtitle = paste0("n = ", nrow(cor_df), "  |  median r = ",
                           round(median(cor_df$pearson_r, na.rm = TRUE), 3))) +
    .theme_velo(base_size = 12)
  save_plot(p_hist, file.path(outdir, paste0(label, "_velocity_pseudotime_cor_hist")), 10, 7)
  message("    Velocity-pseudotime correlation saved: ", nrow(cor_df), " genes.")
  invisible(cor_df)
}

extract_transition_matrix <- function(seu, velo_out, label, outdir) {
  message("  Cell-type transition probability heatmap: ", label)
  if (is.null(velo_out))                              { message("    No velo_out  -  skipping"); return(invisible(NULL)) }
  if (!"cell_type" %in% colnames(seu@meta.data))     { message("    No cell_type  -  skipping"); return(invisible(NULL)) }

  trans_mat <- tryCatch({
    np <- reticulate::import("numpy",        convert = FALSE)
    sp <- reticulate::import("scipy.sparse", convert = FALSE)
    uns_keys <- tryCatch(names(reticulate::py_to_r(velo_out$uns)), error = function(e) character(0))
    if (!"velocity_graph" %in% uns_keys) {
      message("    velocity_graph not in adata.uns  -  skipping"); return(NULL)
    }
    vg_py  <- velo_out$uns[["velocity_graph"]]
    vg_csr <- sp$csr_matrix(vg_py)
    x_r <- as.numeric(reticulate::py_to_r(np$array(vg_csr$data,    dtype = np$float64)))
    i_r <- as.integer(reticulate::py_to_r(np$array(vg_csr$indices, dtype = np$int32)))
    p_r <- as.integer(reticulate::py_to_r(np$array(vg_csr$indptr,  dtype = np$int32)))
    shp <- reticulate::py_to_r(vg_csr$shape)
    m   <- Matrix::sparseMatrix(i = i_r, p = p_r, x = x_r,
                                 dims = c(shp[[1L]], shp[[2L]]), index1 = FALSE)
    cell_ids    <- reticulate::py_to_r(velo_out$obs_names$tolist())
    rownames(m) <- cell_ids; colnames(m) <- cell_ids
    m
  }, error = function(e) { message("    Transition matrix extraction failed: ", e$message); NULL })
  if (is.null(trans_mat)) return(invisible(NULL))

  common_cells <- intersect(rownames(trans_mat), colnames(seu))
  if (length(common_cells) < 20) { message("    Too few common cells  -  skipping"); return(invisible(NULL)) }
  trans_sub  <- trans_mat[common_cells, common_cells, drop = FALSE]
  ct_vals    <- seu@meta.data[common_cells, "cell_type"]
  cell_types <- setdiff(unique(ct_vals), "Unassigned")
  if (length(cell_types) < 2) { message("    Too few cell types  -  skipping"); return(invisible(NULL)) }

  # Mean outgoing transition probability from each source cell type to each target
  ct_trans_mat <- do.call(rbind, lapply(cell_types, function(from) {
    from_c <- common_cells[ct_vals == from]
    vapply(cell_types, function(to) {
      to_c <- common_cells[ct_vals == to]
      if (length(from_c) == 0 || length(to_c) == 0) return(NA_real_)
      mean(Matrix::rowMeans(trans_sub[from_c, to_c, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    }, numeric(1))
  }))
  rownames(ct_trans_mat) <- cell_types; colnames(ct_trans_mat) <- cell_types
  write.csv(as.data.frame(ct_trans_mat),
            file.path(outdir, paste0(label, "_celltype_transition_probabilities.csv")))

  # Row-normalise so each source sums to 1 (relative preference)
  rm_mx <- apply(ct_trans_mat, 1, max, na.rm = TRUE)
  rm_mx[rm_mx == 0 | is.na(rm_mx)] <- 1
  ct_norm <- sweep(ct_trans_mat, 1, rm_mx, "/"); ct_norm[is.nan(ct_norm)] <- 0

  n_ct <- length(cell_types)

  # Compute canvas size that comfortably accommodates the cell matrix PLUS all
  # pheatmap chrome: row/col labels, 45deg angled column text, both dendrograms,
  # multi-line title, and the colour legend.
  #
  # Rule of thumb verified empirically with pheatmap at 300 DPI:
  #   • each cell   ≈ 55 px  (auto-scaled  -  no fixed cellwidth/cellheight)
  #   • col chrome  ≈ 620 px (angled labels + top dendrogram + title lines)
  #   • row chrome  ≈ 520 px (row labels + left dendrogram + legend)
  max_lbl_chars <- max(nchar(cell_types), na.rm = TRUE)
  lbl_px        <- max(300L, max_lbl_chars * 9L)   # ~9 px per character at 300 DPI
  cell_px       <- max(55L,  floor(300L / n_ct))    # fill ~ 1 in of height per n_ct
  col_chrome    <- 620L + lbl_px                    # angled labels can be long
  row_chrome    <- 520L + lbl_px

  h_px <- as.integer(n_ct * cell_px + col_chrome)
  w_px <- as.integer(n_ct * cell_px + row_chrome)

  # Ensure minimum canvas so title + legend are never clipped
  h_px <- max(h_px, 1800L)
  w_px <- max(w_px, 2000L)

  fs_lbl <- max(9, min(13, 11))   # fixed legible label size

  ph_tr <- pheatmap(ct_norm,
                    color        = colorRampPalette(c("grey95", AAAS_COLS[7]))(100),
                    cluster_rows = TRUE, cluster_cols = TRUE,
                    fontsize     = 9,       # title / legend font (reduced)
                    fontsize_row = fs_lbl,
                    fontsize_col = fs_lbl,
                    angle_col    = 45,      # angled col labels prevent overlap
                    # No cellwidth/cellheight: let pheatmap auto-scale cells
                    # to fill the canvas so nothing is ever clipped.
                    main         = paste0("Cell-type transition probabilities (row-normalised)\n",
                                          label),
                    silent       = TRUE)
  save_pheatmap(ph_tr,
                file.path(outdir, paste0(label, "_celltype_transition_heatmap")),
                w_px, h_px)
  message("    Transition heatmap saved: ", n_ct, "  x  ", n_ct, " cell types.")
  invisible(ct_trans_mat)
}

compute_gene_velocity_coherence <- function(seu_list, cond_labels, outdir) {
  message("  Gene-level velocity coherence per condition ...")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  coh_rows <- lapply(seq_along(seu_list), function(i) {
    seu  <- seu_list[[i]]; cond <- cond_labels[i]
    if (!"velocity" %in% names(seu@assays)) return(NULL)
    vm_sp <- .get_layer(seu, "velocity", "counts")
    nz    <- Matrix::rowSums(vm_sp != 0, na.rm = TRUE)
    vm_sp <- vm_sp[nz >= 10, , drop = FALSE]
    if (nrow(vm_sp) < 5) return(NULL)
    gm  <- as.numeric(Matrix::rowMeans(vm_sp, na.rm = TRUE))
    names(gm) <- rownames(vm_sp)
    coh <- vapply(rownames(vm_sp), function(g) {
      v <- as.numeric(vm_sp[g, ]); v <- v[!is.na(v)]
      if (length(v) == 0) return(NA_real_)
      mean(sign(v) == sign(gm[g]), na.rm = TRUE)
    }, numeric(1))
    data.frame(gene = rownames(vm_sp), coherence = coh, mean_velocity = gm,
               condition = cond, stringsAsFactors = FALSE)
  })
  coh_df <- do.call(rbind, Filter(Negate(is.null), coh_rows))
  if (is.null(coh_df) || nrow(coh_df) == 0) return(invisible(NULL))

  write.csv(coh_df, file.path(outdir, "gene_velocity_coherence_per_condition.csv"),
            row.names = FALSE)

  if (length(unique(coh_df$condition)) >= 2) {
    # Wide matrix for heatmap
    coh_wide <- tidyr::pivot_wider(coh_df[, c("gene", "coherence", "condition")],
                                    names_from = "condition", values_from = "coherence")
    coh_wide <- as.data.frame(coh_wide)
    rownames(coh_wide) <- coh_wide$gene; coh_wide$gene <- NULL
    coh_wide <- coh_wide[complete.cases(coh_wide), , drop = FALSE]

    if (nrow(coh_wide) >= 5) {
      gv_coh  <- apply(coh_wide, 1, var, na.rm = TRUE)
      top_coh <- names(sort(gv_coh, decreasing = TRUE))[seq_len(min(80, length(gv_coh)))]
      mat_coh <- as.matrix(coh_wide[top_coh, , drop = FALSE])
      h_px <- max(1200, nrow(mat_coh) * 16); w_px <- max(600, ncol(mat_coh) * 90)
      fs_r <- max(4, min(8, floor(h_px / nrow(mat_coh) / PLT$dpi * 72 * 0.5)))
      ph_coh <- tryCatch(
        pheatmap(mat_coh,
                 color        = colorRampPalette(c(AAAS_COLS[1], "grey95", AAAS_COLS[3]))(100),
                 breaks       = seq(0, 1, length.out = 101),
                 cluster_rows = nrow(mat_coh) >= 2, cluster_cols = FALSE,
                 fontsize      = 7,
                 fontsize_row = fs_r, fontsize_col = 10,
                 main = "Gene-level velocity coherence across conditions (top variable genes)",
                 silent = TRUE),
        error = function(e) { message("    Coherence pheatmap failed: ", e$message); NULL })
      if (!is.null(ph_coh))
        save_pheatmap(ph_coh, file.path(outdir, "gene_velocity_coherence_heatmap"), w_px * 2, h_px)
    }

    # Violin per condition with pairwise significance brackets
    viol_data  <- coh_df[!is.na(coh_df$coherence), ]
    viol_data$condition <- factor(viol_data$condition, levels = CONDITIONS)
    comp_pairs <- lapply(COMPARISONS, as.character)
    y_max      <- max(viol_data$coherence, na.rm = TRUE)
    step       <- (y_max - min(viol_data$coherence, na.rm = TRUE)) * 0.12
    sig_map    <- function(p) {
      if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else
        sprintf("p=%.2f", p)
    }
    p_viol_coh <- ggplot(viol_data, aes(x = condition, y = coherence, fill = condition)) +
      geom_violin(trim = FALSE, scale = "width", alpha = 0.85, colour = NA) +
      geom_boxplot(width = 0.08, outlier.size = 0.5, fill = "white",
                   colour = "grey30", linewidth = 0.45) +
      ggsignif::geom_signif(
        comparisons      = comp_pairs,
        test             = "wilcox.test",
        map_signif_level = sig_map,
        step_increase    = 0.12,
        tip_length       = 0.015,
        textsize         = 3.5,
        vjust            = 0.3,
        colour           = "grey20",
        fontface         = "bold") +
      scale_fill_manual(values = COND_PALETTE) +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.22))) +
      labs(x = NULL, y = "Velocity coherence (gene level)",
           title    = "Gene-level velocity coherence by condition",
           subtitle = "Wilcoxon test  |  * p<0.05  ** p<0.01  *** p<0.001") +
      .theme_velo(base_size = 12) +
      theme(legend.position  = "none",
            axis.text.x      = element_text(size = 12, face = "bold"),
            axis.text.y      = element_text(size = 11),
            axis.title.y     = element_text(size = 12),
            plot.title       = element_text(size = 13, face = "bold"),
            plot.subtitle    = element_text(size = 10, colour = "grey40"))
    save_plot(p_viol_coh, file.path(outdir, "gene_velocity_coherence_violin"), 12, 11)
  }
  message("    Gene-level coherence scores saved.")
  invisible(coh_df)
}

run_latent_time_ks_factorial <- function(cond_objects, outdir) {
  message("  Latent time KS  -  factorial summary (cell type  x  condition) ...")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  all_ks <- lapply(COMPARISONS, function(comp) {
    ql <- comp[1]; rl <- comp[2]
    seu_q <- cond_objects[[ql]]; seu_r <- cond_objects[[rl]]
    if (is.null(seu_q) || is.null(seu_r)) return(NULL)
    ks <- run_latent_time_ks_tests(seu_q, seu_r, ql, rl, outdir)
    if (!is.null(ks) && nrow(ks) > 0) { ks$comparison <- paste0(ql, "_vs_", rl); ks } else NULL
  })
  ks_all <- do.call(rbind, Filter(Negate(is.null), all_ks))
  if (is.null(ks_all) || nrow(ks_all) == 0) return(invisible(NULL))

  write.csv(ks_all, file.path(outdir, "latent_time_KS_factorial_all.csv"), row.names = FALSE)

  cts  <- sort(unique(ks_all$cell_type))
  cmps <- unique(ks_all$comparison)
  if (length(cts) < 2 || length(cmps) < 1) return(invisible(ks_all))

  d_mat <- do.call(cbind, lapply(cmps, function(cmp) {
    sub <- ks_all[ks_all$comparison == cmp, ]
    v   <- setNames(sub$delta_mean_lt, sub$cell_type)
    v[cts]
  }))
  colnames(d_mat) <- cmps; rownames(d_mat) <- cts

  lim  <- max(abs(d_mat), na.rm = TRUE) * 0.95; lim <- max(lim, 0.001)
  h_px <- max(900, length(cts) * 50); w_px <- max(600, length(cmps) * 120)
  fs_r <- max(4, min(9, floor(h_px / length(cts) / PLT$dpi * 72 * 0.5)))
  ph_ks <- tryCatch(
    pheatmap(d_mat,
             color        = colorRampPalette(c(AAAS_COLS[1], "grey95", AAAS_COLS[2]))(100),
             breaks       = seq(-lim, lim, length.out = 101),
             cluster_rows = nrow(d_mat) >= 2, cluster_cols = FALSE,
             fontsize_row = fs_r, fontsize_col = 9,
             na_col       = "grey90",
             main         = "Delta  Mean latent time  -  cell type  x  condition comparison",
             silent       = TRUE),
    error = function(e) { message("    KS factorial pheatmap failed: ", e$message); NULL })
  if (!is.null(ph_ks))
    save_pheatmap(ph_ks, file.path(outdir, "latent_time_KS_factorial_heatmap"), w_px, h_px)
  message("    Factorial KS heatmap saved: ",
          length(cts), " cell types  x  ", length(cmps), " comparisons.")
  invisible(ks_all)
}

run_differentially_variable_genes <- function(cond_objects, outdir) {
  message("  Differentially variable genes (scran::modelGeneVar) ...")
  if (!requireNamespace("scran", quietly = TRUE)) {
    message("    scran not available  -  skipping"); return(invisible(NULL))
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  dvg_results <- list()
  for (comp in COMPARISONS) {
    ql <- comp[1]; rl <- comp[2]; lbl <- paste0(ql, "_vs_", rl)
    seu_q <- cond_objects[[ql]]; seu_r <- cond_objects[[rl]]
    if (is.null(seu_q) || is.null(seu_r)) next

    tryCatch({
      mat_q <- .get_layer(seu_q, "RNA", "counts")
      mat_r <- .get_layer(seu_r, "RNA", "counts")
      cg    <- intersect(rownames(mat_q), rownames(mat_r))
      if (length(cg) < 100) { message("    Too few common genes  -  skipping ", lbl); next }

      # Log-normalise each condition separately
      norm_q <- .norm_log1p(mat_q[cg, , drop = FALSE])
      norm_r <- .norm_log1p(mat_r[cg, , drop = FALSE])
      sce_q  <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = norm_q))
      sce_r  <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = norm_r))

      dv_q   <- scran::modelGeneVar(sce_q)
      dv_r   <- scran::modelGeneVar(sce_r)
      cg2    <- intersect(rownames(dv_q), rownames(dv_r))

      var_df <- data.frame(
        gene         = cg2,
        total_var_q  = as.numeric(dv_q[cg2, "total"]),
        total_var_r  = as.numeric(dv_r[cg2, "total"]),
        bio_var_q    = as.numeric(dv_q[cg2, "bio"]),
        bio_var_r    = as.numeric(dv_r[cg2, "bio"]),
        tech_var_q   = as.numeric(dv_q[cg2, "tech"]),
        tech_var_r   = as.numeric(dv_r[cg2, "tech"]),
        p_val_q      = as.numeric(dv_q[cg2, "p.value"]),
        p_val_r      = as.numeric(dv_r[cg2, "p.value"]),
        stringsAsFactors = FALSE)
      var_df$delta_total_var <- var_df$total_var_q - var_df$total_var_r
      var_df$delta_bio_var   <- var_df$bio_var_q   - var_df$bio_var_r
      var_df$padj_q          <- p.adjust(var_df$p_val_q, "BH")
      var_df$padj_r          <- p.adjust(var_df$p_val_r, "BH")
      var_df$avg_bio_var     <- (var_df$bio_var_q + var_df$bio_var_r) / 2

      # Category: significantly higher biological variance in one condition only
      var_df$dv_category <- "Other"
      var_df$dv_category[var_df$padj_q < 0.05 & var_df$padj_r >= 0.05 &
                         var_df$delta_bio_var > 0] <- paste0("HighVar_", ql)
      var_df$dv_category[var_df$padj_r < 0.05 & var_df$padj_q >= 0.05 &
                         var_df$delta_bio_var < 0] <- paste0("HighVar_", rl)

      var_df <- var_df[order(abs(var_df$delta_bio_var), decreasing = TRUE), ]
      write.csv(var_df, file.path(outdir, paste0(lbl, "_differentially_variable_genes.csv")),
                row.names = FALSE)

      n_dv <- sum(var_df$dv_category != "Other")
      message("    ", lbl, ": ", n_dv, " DVGs  (HighVar_", ql, ": ",
              sum(var_df$dv_category == paste0("HighVar_", ql)), "  HighVar_", rl, ": ",
              sum(var_df$dv_category == paste0("HighVar_", rl)), ")")

      # Scatter: Δbio variance vs mean bio variance
      sig_dv  <- var_df[var_df$dv_category != "Other", ]
      top_dvg <- rbind(head(sig_dv[order(-sig_dv$delta_bio_var), ], 8),
                       head(sig_dv[order( sig_dv$delta_bio_var), ], 8))
      cat_cols_dv <- c(Other = "grey80")
      if (paste0("HighVar_", ql) %in% var_df$dv_category)
        cat_cols_dv[paste0("HighVar_", ql)] <- .pal_get(COND_PALETTE, ql) %||% AAAS_COLS[2]
      if (paste0("HighVar_", rl) %in% var_df$dv_category)
        cat_cols_dv[paste0("HighVar_", rl)] <- .pal_get(COND_PALETTE, rl) %||% AAAS_COLS[1]

      p_dv <- ggplot(var_df, aes(avg_bio_var, delta_bio_var, colour = dv_category)) +
        geom_point(alpha = 0.40, size = 0.9) +
        ggrepel::geom_label_repel(data = top_dvg, aes(label = gene),
                                   size = 2.8, show.legend = FALSE, max.overlaps = 20) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_colour_manual(values = cat_cols_dv) +
        labs(x = "Mean biological variance",
             y = paste0("Delta  Bio variance (", ql, " - ", rl, ")"),
             title    = paste("Differentially variable genes:", lbl),
             subtitle = paste0(n_dv, " DVGs  |  HighVar_", ql, ": ",
                               sum(var_df$dv_category == paste0("HighVar_", ql)),
                               "  HighVar_", rl, ": ",
                               sum(var_df$dv_category == paste0("HighVar_", rl))),
             colour = NULL) +
        .theme_velo(base_size = 13)
      save_plot(p_dv, file.path(outdir, paste0(lbl, "_differential_variance_scatter")),
                PLT$w_single, PLT$h_single)

      dvg_results[[lbl]] <- var_df
    }, error = function(e) message("    DVG failed for ", lbl, ": ", e$message))
  }
  invisible(dvg_results)
}

plot_custom_summary_visualizations <- function(base_outdir) {

  pdir <- file.path(base_outdir, "custom_summary_plots")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  ## ── internal colour scheme (mirrors PPTX theme) ──────────────────────────
  .DARK_NAV <- "#1B2A4A"
  .TEAL     <- "#0D9E8F"
  .GOLD     <- "#E8A838"
  .RED      <- "#C05050"
  .LGRAY    <- "#D1D5DB"
  .COND_COL <- setNames(c(AAAS_COLS[1], AAAS_COLS[2], AAAS_COLS[3]),
                        c("Control", "Eftud2", "Sf3b4"))

  ## ── internal theme ────────────────────────────────────────────────────────
  .tv <- function(base_size = 13) {
    theme_minimal(base_size = base_size) +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border       = element_blank(),
        axis.line          = element_line(colour = "#666666", linewidth = 0.5),
        plot.title         = element_text(face = "bold", colour = .DARK_NAV, size = 16),
        plot.subtitle      = element_text(colour = "#444444", size = 13),
        plot.caption       = element_text(colour = "#666666", size = 10, face = "italic"),
        legend.position    = "bottom",
        legend.title       = element_text(face = "bold", size = 12),
        strip.text         = element_text(face = "bold", size = 13)
      )
  }

  ## ── safe CSV reader (returns NULL if file missing) ────────────────────────
  .rc <- function(path) {
    if (!file.exists(path)) {
      message("    [custom_plots] missing: ", basename(path), "  -  skipping.")
      return(NULL)
    }
    tryCatch(read.csv(path, stringsAsFactors = FALSE),
             error = function(e) { message("    [custom_plots] read error: ", e$message); NULL })
  }

  ## ── data paths ────────────────────────────────────────────────────────────
  diff_dir   <- file.path(base_outdir, "differential")
  combo_dir  <- file.path(base_outdir, "qc_metrics")

  ## ── Load shared datasets ─────────────────────────────────────────────────
  comp <- .rc(file.path(combo_dir, "conditions_celltype_composition_metrics.csv"))
  if (!is.null(comp)) {
    comp <- comp[!tolower(comp$cell_type) %in% c("unassigned", "all", "total"), ]
    # ensure expected condition columns exist
    missing_conds <- setdiff(CONDITIONS, colnames(comp))
    if (length(missing_conds) > 0) {
      message("    [custom_plots] composition CSV missing columns: ",
              paste(missing_conds, collapse = ", "), "  -  composition plots skipped.")
      comp <- NULL
    }
  }

  ctrl_vel <- .rc(file.path(diff_dir, "Control_celltype_velocity_summary.csv"))
  eft_vel  <- .rc(file.path(diff_dir, "Eftud2_celltype_velocity_summary.csv"))
  sf3_vel  <- .rc(file.path(diff_dir, "Sf3b4_celltype_velocity_summary.csv"))
  eft_ct   <- .rc(file.path(diff_dir, "Eftud2_vs_Control_by_celltype_summary.csv"))
  sf3_ct   <- .rc(file.path(diff_dir, "Sf3b4_vs_Control_by_celltype_summary.csv"))
  shared   <- .rc(file.path(diff_dir, "shared_significant_genes.csv"))
  eft_sig  <- .rc(file.path(diff_dir, "Eftud2_vs_Control_significant_diff_velocity.csv"))
  sf3_sig  <- .rc(file.path(diff_dir, "Sf3b4_vs_Control_significant_diff_velocity.csv"))

  message("  [custom_plots] Starting custom summary visualisations …")

  if (!is.null(comp) && all(c("cell_type", CONDITIONS) %in% colnames(comp))) {
    message("    Plot 1: Composition dodged bars")
    top30 <- comp
    top30$max_pct <- do.call(pmax, as.list(comp[, CONDITIONS, drop = FALSE]))
    top30 <- top30[order(-top30$max_pct), ][seq_len(min(30L, nrow(top30))), ]
    top30 <- top30[order(top30$Control), ]
    top30$cell_type <- factor(top30$cell_type, levels = top30$cell_type)

    df_long1 <- tidyr::pivot_longer(top30[, c("cell_type", CONDITIONS)],
                                    cols = all_of(CONDITIONS),
                                    names_to  = "Condition",
                                    values_to = "pct")
    df_long1$Condition <- factor(df_long1$Condition, levels = CONDITIONS)

    p1 <- ggplot(df_long1, aes(x = pct, y = cell_type, fill = Condition)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7) +
      geom_text(
        aes(label = ifelse(pct >= 2, sprintf("%.1f%%", pct), "")),
        position = position_dodge(width = 0.75), hjust = -0.15,
        size = 3.2, fontface = "bold"
      ) +
      geom_vline(xintercept = c(5, 10), linetype = "dashed",
                 colour = "#CCCCCC", linewidth = 0.5) +
      scale_fill_manual(values = .COND_COL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
      scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
      labs(
        title   = "Cell Type Composition Across Conditions (Top Cell Types)",
        x       = "Cell Proportion (%)",
        y       = NULL,
        fill    = "Condition",
        caption = "Dashed lines at 5% and 10% thresholds"
      ) +
      .tv() +
      theme(panel.grid.major.y = element_line(colour = "#EEEEEE"),
            legend.position    = "top")
    save_plot(p1, file.path(pdir, "01_composition_dodged_bars"), 16, 10)
  }

  if (!is.null(comp) && all(c("cell_type", "Control", "Eftud2") %in% colnames(comp))) {
    message("    Plot 2: Lollipop Eftud2 composition")
    df_e <- comp
    df_e$log2fc    <- log2((df_e$Eftud2 + 0.1) / (df_e$Control + 0.1))
    df_e$direction <- ifelse(df_e$log2fc > 0, "Increased", "Decreased")
    df_e$cell_type <- gsub("_", " ", df_e$cell_type)
    df_e <- df_e[order(df_e$log2fc), ]
    df_e$cell_type <- factor(df_e$cell_type, levels = df_e$cell_type)
    top5e <- df_e[rank(-df_e$log2fc) <= 5 | rank(df_e$log2fc) <= 5, ]

    p2 <- ggplot(df_e, aes(x = log2fc, y = cell_type, colour = direction)) +
      geom_segment(aes(x = 0, xend = log2fc, yend = cell_type),
                   linewidth = 1.4, alpha = 0.7) +
      geom_point(size = 3.5) +
      geom_vline(xintercept = 0, linewidth = 1.0, colour = "#444444") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                 colour = "#CCCCCC", linewidth = 0.7) +
      ggrepel::geom_text_repel(
        data = top5e,
        aes(label = sprintf("%.2f x ", log2fc)),
        size = 3.5, fontface = "bold", show.legend = FALSE,
        direction = "y", box.padding = 0.3
      ) +
      scale_colour_manual(values = c(Increased = .RED, Decreased = .TEAL)) +
      labs(
        title   = "Eftud2 vs Control  -  Cell Type Composition Changes",
        x       = expression(Log[2]~"Fold-Change vs Control"),
        y       = NULL,
        colour  = NULL,
        caption = "Dashed lines at +/-1 fold-change threshold"
      ) +
      .tv() +
      theme(legend.position    = "top",
      panel.grid.major.x = element_line(colour = "#EEEEEE"),
            panel.grid.major.y = element_blank())
    save_plot(p2, file.path(pdir, "02_lollipop_eftud2_composition"), 10, 14)
  }

  if (!is.null(comp) && all(c("cell_type", "Control", "Sf3b4") %in% colnames(comp))) {
    message("    Plot 3: Lollipop Sf3b4 composition")
    df_s <- comp
    df_s$log2fc    <- log2((df_s$Sf3b4 + 0.1) / (df_s$Control + 0.1))
    df_s$direction <- ifelse(df_s$log2fc > 0, "Increased", "Decreased")
    df_s$cell_type <- gsub("_", " ", df_s$cell_type)
    df_s <- df_s[order(df_s$log2fc), ]
    df_s$cell_type <- factor(df_s$cell_type, levels = df_s$cell_type)
    top5s <- df_s[rank(-df_s$log2fc) <= 5 | rank(df_s$log2fc) <= 5, ]

    p3 <- ggplot(df_s, aes(x = log2fc, y = cell_type, colour = direction)) +
      geom_segment(aes(x = 0, xend = log2fc, yend = cell_type),
                   linewidth = 1.4, alpha = 0.7) +
      geom_point(size = 3.5) +
      geom_vline(xintercept = 0, linewidth = 1.0, colour = "#444444") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                 colour = "#CCCCCC", linewidth = 0.7) +
      ggrepel::geom_text_repel(
        data = top5s,
        aes(label = sprintf("%.2f x ", log2fc)),
        size = 3.5, fontface = "bold", show.legend = FALSE,
        direction = "y", box.padding = 0.3
      ) +
      scale_colour_manual(values = c(Increased = .RED, Decreased = .TEAL)) +
      labs(
        title   = "Sf3b4 vs Control  -  Cell Type Composition Changes",
        x       = expression(Log[2]~"Fold-Change vs Control"),
        y       = NULL,
        colour  = NULL
      ) +
      .tv() +
      theme(legend.position    = "top",
            panel.grid.major.x = element_line(colour = "#EEEEEE"),
            panel.grid.major.y = element_blank())
    save_plot(p3, file.path(pdir, "03_lollipop_sf3b4_composition"), 10, 14)
  }

  if (!is.null(comp) && all(c("cell_type", "Control", "Eftud2", "Sf3b4") %in% colnames(comp))) {
    message("    Plot 4: Lollipop both mutants (faceted)")
    df_both <- comp
    df_both$eft_fc    <- log2((df_both$Eftud2 + 0.1) / (df_both$Control + 0.1))
    df_both$sf3_fc    <- log2((df_both$Sf3b4  + 0.1) / (df_both$Control + 0.1))
    df_both$cell_type <- gsub("_", " ", df_both$cell_type)
    df_both <- df_both[order(df_both$sf3_fc), ]
    df_both$cell_type <- factor(df_both$cell_type, levels = df_both$cell_type)

    df_long4 <- tidyr::pivot_longer(df_both[, c("cell_type", "eft_fc", "sf3_fc")],
                                    cols      = c("eft_fc", "sf3_fc"),
                                    names_to  = "mutant",
                                    values_to = "log2fc")
    df_long4$mutant    <- ifelse(df_long4$mutant == "eft_fc", "Eftud2", "Sf3b4")
    df_long4$direction <- ifelse(df_long4$log2fc > 0, "Increased", "Decreased")

    p4 <- ggplot(df_long4, aes(x = log2fc, y = cell_type, colour = direction)) +
      geom_segment(aes(x = 0, xend = log2fc, yend = cell_type),
                   linewidth = 1.2, alpha = 0.75) +
      geom_point(size = 3) +
      geom_vline(xintercept = 0, linewidth = 0.8, colour = "#444444") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                 colour = "#CCCCCC", linewidth = 0.6) +
      scale_colour_manual(values = c(Increased = .RED, Decreased = .TEAL)) +
      facet_wrap(~mutant, nrow = 1, scales = "free_x") +
      labs(
        title    = "Cell Type Composition Fold-Changes vs Control (Log2 Scale)",
        subtitle = "Cell types sorted by Sf3b4 fold-change",
        x        = expression(Log[2]~"Fold-Change"),
        y        = NULL,
        colour   = NULL
      ) +
      .tv() +
      theme(
        legend.position    = "top",
        panel.grid.major.x = element_line(colour = "#EEEEEE"),
        panel.grid.major.y = element_blank(),
        strip.background   = element_rect(fill = .DARK_NAV, colour = NA),
        strip.text         = element_text(colour = "white", face = "bold", size = 14)
      )
    save_plot(p4, file.path(pdir, "04_lollipop_both_mutants"), 18, 14)
  }

  if (!is.null(comp) && !is.null(eft_ct) && !is.null(sf3_ct) &&
      all(c("cell_type", "n_sig_genes") %in% colnames(eft_ct)) &&
      all(c("cell_type", "n_sig_genes") %in% colnames(sf3_ct))) {
    message("    Plot 5: Bubble composition  x  sig-gene overlay")

    ct_genes <- merge(
      eft_ct[, c("cell_type", "n_sig_genes")],
      sf3_ct[, c("cell_type", "n_sig_genes")],
      by = "cell_type", suffixes = c("_eft", "_sf3"), all = TRUE)
    ct_genes$n_sig_genes_eft[is.na(ct_genes$n_sig_genes_eft)] <- 0
    ct_genes$n_sig_genes_sf3[is.na(ct_genes$n_sig_genes_sf3)] <- 0
    ct_genes$avg_ngenes <- (ct_genes$n_sig_genes_eft + ct_genes$n_sig_genes_sf3) / 2

    df5 <- comp
    df5$eft_fc    <- log2((df5$Eftud2 + 0.1) / (df5$Control + 0.1))
    df5$sf3_fc    <- log2((df5$Sf3b4  + 0.1) / (df5$Control + 0.1))
    df5$sf3_delta <- df5$Sf3b4 - df5$Control
    df5$label     <- gsub("_", " ", df5$cell_type)
    df5 <- merge(df5, ct_genes[, c("cell_type", "avg_ngenes")], by = "cell_type", all.x = TRUE)
    df5$avg_ngenes[is.na(df5$avg_ngenes)] <- median(ct_genes$avg_ngenes, na.rm = TRUE)

    p5 <- ggplot(df5, aes(x = eft_fc, y = sf3_fc)) +
      annotate("rect", xmin = 0,   xmax = Inf,  ymin = 0,   ymax = Inf,
               fill = .RED,  alpha = 0.04) +
      annotate("rect", xmin = -Inf, xmax = 0,  ymin = -Inf, ymax = 0,
               fill = .TEAL, alpha = 0.04) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888") +
      geom_hline(yintercept = c(-1, 1), linetype = "dotted",
                 colour = "#CCCCCC", linewidth = 0.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dotted",
                 colour = "#CCCCCC", linewidth = 0.5) +
      geom_point(aes(size = avg_ngenes, fill = sf3_delta),
                 shape = 21, colour = "#555555", stroke = 0.5, alpha = 0.88) +
      ggrepel::geom_text_repel(
        aes(label = label),
        size = 3.0, max.overlaps = 30,
        segment.colour = "#AAAAAA", segment.size = 0.4,
        box.padding = 0.4, point.padding = 0.2
      ) +
      scale_fill_gradient2(
        low = .TEAL, mid = "#F5F5F5", high = .RED, midpoint = 0,
        name = "Sf3b4 Delta  proportion (pp)"
      ) +
      scale_size_continuous(name = "Avg # sig. velocity genes", range = c(3, 15)) +
      annotate("text", x =  Inf, y =  Inf, label = "Both Expanded",
               hjust = 1.1, vjust = 1.5, colour = .RED,  size = 4, fontface = "italic") +
      annotate("text", x = -Inf, y =  Inf, label = "Only Sf3b4 Expanded",
               hjust = -0.1, vjust = 1.5, colour = "#999999", size = 4, fontface = "italic") +
      annotate("text", x =  Inf, y = -Inf, label = "Only Eftud2 Expanded",
               hjust = 1.1, vjust = -0.5, colour = "#999999", size = 4, fontface = "italic") +
      annotate("text", x = -Inf, y = -Inf, label = "Both Contracted",
               hjust = -0.1, vjust = -0.5, colour = .TEAL, size = 4, fontface = "italic") +
      labs(
        title    = "Cell Type Composition Shifts: Eftud2 vs Sf3b4",
        subtitle = "Bubble size = average # significant velocity genes",
        x        = "Eftud2 vs Control  (Log2 FC)",
        y        = "Sf3b4 vs Control  (Log2 FC)"
      ) +
      .tv() +
      theme(legend.box = "vertical")
    save_plot(p5, file.path(pdir, "05_bubble_composition_signgene"), 13, 10)
  }

  vel_csvs_ok <- !is.null(ctrl_vel) && !is.null(eft_vel) && !is.null(sf3_vel) &&
    all(c("cell_type", "mean_speed", "mean_coherence", "mean_latent_time") %in% colnames(ctrl_vel))
  if (vel_csvs_ok) {
    message("    Plot 6: Velocity metrics tile heatmap")

    vel_all <- merge(
      ctrl_vel[, c("cell_type", "mean_speed", "mean_coherence", "mean_latent_time")],
      eft_vel[,  c("cell_type", "mean_speed", "mean_coherence", "mean_latent_time")],
      by = "cell_type", suffixes = c("_ctrl", "_eft"))
    vel_all <- merge(vel_all,
      sf3_vel[,  c("cell_type", "mean_speed", "mean_coherence", "mean_latent_time")],
      by = "cell_type", suffixes = c("", "_sf3"))
    names(vel_all)[names(vel_all) == "mean_speed"]        <- "mean_speed_sf3"
    names(vel_all)[names(vel_all) == "mean_coherence"]    <- "mean_coherence_sf3"
    names(vel_all)[names(vel_all) == "mean_latent_time"]  <- "mean_latent_time_sf3"
    vel_all <- vel_all[order(vel_all$mean_speed_ctrl, na.last = TRUE), ]

    .make_tile <- function(df, ctrl_col, eft_col, sf3_col,
                           title, palette = "YlOrRd", direction = 1) {
      df_long <- data.frame(
        cell_type = rep(df$cell_type, 3),
        Condition = rep(c("Control", "Eftud2", "Sf3b4"), each = nrow(df)),
        value     = c(df[[ctrl_col]], df[[eft_col]], df[[sf3_col]])
      )
      df_long$cell_type <- factor(df_long$cell_type, levels = df$cell_type)
      df_long$Condition <- factor(df_long$Condition, levels = CONDITIONS)
      ggplot(df_long, aes(x = Condition, y = cell_type, fill = value)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(aes(label = ifelse(is.na(value), "NA",
                                     sprintf("%.3f", value))),
                  size = 2.4, colour = "black") +
        scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
        scale_fill_distiller(palette = palette, direction = direction,
                             name = title, na.value = "#EEEEEE") +
        labs(title = title, x = NULL, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.x     = element_text(face = "bold", size = 12),
          axis.text.y     = element_text(size = 9),
          plot.title      = element_text(face = "bold", colour = .DARK_NAV, hjust = 0.5),
          legend.position = "right"
        )
    }

    p_speed <- .make_tile(vel_all, "mean_speed_ctrl", "mean_speed_eft", "mean_speed_sf3",
                          "Mean Speed",     "YlOrRd",  1)
    p_coh   <- .make_tile(vel_all, "mean_coherence_ctrl", "mean_coherence_eft", "mean_coherence_sf3",
                          "Mean Coherence", "Blues",   1)
    p_lat   <- .make_tile(vel_all, "mean_latent_time_ctrl", "mean_latent_time_eft", "mean_latent_time_sf3",
                          "Latent Time",    "RdYlGn", -1)

    p6 <- (p_speed | p_coh | p_lat) +
      patchwork::plot_annotation(
        title = "Velocity Metrics per Cell Type Across Conditions",
        theme = theme(plot.title = element_text(face = "bold", colour = .DARK_NAV, size = 16))
      )
    save_plot(p6, file.path(pdir, "06_velocity_metrics_tile_heatmap"), 22, 14)
  }

  if (!is.null(eft_ct) && !is.null(sf3_ct) &&
      all(c("cell_type", "mean_vel_query", "mean_vel_ref") %in% colnames(eft_ct)) &&
      all(c("cell_type", "mean_vel_query", "mean_vel_ref") %in% colnames(sf3_ct))) {
    message("    Plot 7: Delta-velocity scatter")

    df7 <- merge(
      within(eft_ct, delta_eft <- mean_vel_query - mean_vel_ref)[, c("cell_type", "delta_eft")],
      within(sf3_ct, delta_sf3 <- mean_vel_query - mean_vel_ref)[, c("cell_type", "delta_sf3")],
      by = "cell_type")
    df7$quad <- with(df7, dplyr::case_when(
      delta_eft > 0 & delta_sf3 > 0 ~ "Both increased",
      delta_eft < 0 & delta_sf3 < 0 ~ "Both decreased",
      delta_eft > 0                  ~ "Only Eftud2 increased",
      TRUE                           ~ "Only Sf3b4 increased"
    ))
    df7$label <- gsub("_", " ", df7$cell_type)

    quad_col <- c(
      "Both increased"        = "#8B5E3C",
      "Both decreased"        = .TEAL,
      "Only Eftud2 increased" = .GOLD,
      "Only Sf3b4 increased"  = .RED
    )
    lim <- range(c(df7$delta_eft, df7$delta_sf3), na.rm = TRUE)

    p7 <- ggplot(df7, aes(x = delta_eft, y = delta_sf3, colour = quad)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                  colour = "#CCCCCC", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888") +
      geom_point(size = 4, alpha = 0.85) +
      ggrepel::geom_text_repel(
        aes(label = label),
        size = 3.2, max.overlaps = 35,
        segment.colour = "#AAAAAA", segment.size = 0.4,
        box.padding = 0.35, point.padding = 0.2
      ) +
      scale_colour_manual(values = quad_col, name = NULL) +
      annotate("text", x = max(lim) * 0.95, y = max(lim) * 0.6,
               label = "y = x", colour = "#AAAAAA", size = 3.5, fontface = "italic") +
      labs(
        title    = "Per-Cell-Type Velocity Change: Eftud2 vs Sf3b4",
        subtitle = "Cell types above diagonal = stronger Sf3b4 effect",
        x        = expression(Delta*"Velocity (Eftud2 - Control)"),
        y        = expression(Delta*"Velocity (Sf3b4 - Control)")
      ) +
      .tv() +
      theme(legend.position    = "top",
            panel.grid.major   = element_line(colour = "#EEEEEE"))
    save_plot(p7, file.path(pdir, "07_delta_velocity_scatter"), 12, 10)
  }

  if (vel_csvs_ok) {
    message("    Plot 8: Condition-level velocity summary bars (data-driven)")

    .weighted_mean <- function(df, col) {
      has_ncells <- "n_cells" %in% colnames(df)
      w <- if (has_ncells) df$n_cells else rep(1, nrow(df))
      wm <- weighted.mean(df[[col]], w = w, na.rm = TRUE)
      round(wm, 4)
    }

    vel_sum <- data.frame(
      Condition   = factor(CONDITIONS, levels = CONDITIONS),
      Speed       = c(.weighted_mean(ctrl_vel, "mean_speed"),
                      .weighted_mean(eft_vel,  "mean_speed"),
                      .weighted_mean(sf3_vel,  "mean_speed")),
      Coherence   = c(.weighted_mean(ctrl_vel, "mean_coherence"),
                      .weighted_mean(eft_vel,  "mean_coherence"),
                      .weighted_mean(sf3_vel,  "mean_coherence")),
      Latent_Time = c(.weighted_mean(ctrl_vel, "mean_latent_time"),
                      .weighted_mean(eft_vel,  "mean_latent_time"),
                      .weighted_mean(sf3_vel,  "mean_latent_time")),
      stringsAsFactors = FALSE
    )
    vel_long8 <- tidyr::pivot_longer(vel_sum, cols = -Condition,
                                     names_to = "Metric", values_to = "Value")
    vel_long8$Metric <- dplyr::recode(vel_long8$Metric,
                                      Speed = "Mean Speed",
                                      Coherence = "Coherence",
                                      Latent_Time = "Latent Time")

    p8 <- ggplot(vel_long8, aes(x = Condition, y = Value, fill = Condition)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.3f", Value)),
                vjust = -0.4, size = 4.5, fontface = "bold") +
      scale_fill_manual(values = .COND_COL) +
      facet_wrap(~Metric, scales = "free_y", nrow = 1) +
      labs(
        title   = "RNA Velocity Summary Metrics Across Conditions",
        x       = NULL,
        y       = "Metric Value",
        caption = "n_cells-weighted means across all cell types"
      ) +
      .tv() +
      theme(
        legend.position    = "none",
        panel.grid.major.y = element_line(colour = "#EEEEEE"),
        strip.background   = element_rect(fill = .DARK_NAV, colour = NA),
        strip.text         = element_text(colour = "white", face = "bold", size = 14)
      )
    save_plot(p8, file.path(pdir, "08_velocity_summary_bars"), 14, 7)

    # Also write the aggregated values as a CSV for reproducibility
    write.csv(vel_sum, file.path(pdir, "08_velocity_condition_summary.csv"),
              row.names = FALSE)
  }

  if (!is.null(eft_ct) && !is.null(sf3_ct) &&
      all(c("cell_type", "n_sig_genes") %in% colnames(eft_ct)) &&
      all(c("cell_type", "n_sig_genes") %in% colnames(sf3_ct))) {
    message("    Plot 9: Sig genes per cell type grouped bar")

    ct_merge <- merge(
      eft_ct[, c("cell_type", "n_sig_genes")],
      sf3_ct[, c("cell_type", "n_sig_genes")],
      by = "cell_type", suffixes = c("_eft", "_sf3"), all = TRUE)
    ct_merge$n_sig_genes_eft[is.na(ct_merge$n_sig_genes_eft)] <- 0
    ct_merge$n_sig_genes_sf3[is.na(ct_merge$n_sig_genes_sf3)] <- 0
    ct_merge <- ct_merge[order(-ct_merge$n_sig_genes_sf3), ]
    ct_merge$cell_type_label <- factor(gsub("_", " ", ct_merge$cell_type),
                                       levels = rev(gsub("_", " ", ct_merge$cell_type)))


    ct_long9 <- rbind(
      data.frame(cell_type = ct_merge$cell_type_label,
                 Mutant = "Eftud2", n_genes = ct_merge$n_sig_genes_eft),
      data.frame(cell_type = ct_merge$cell_type_label,
                 Mutant = "Sf3b4",  n_genes = ct_merge$n_sig_genes_sf3)
    )
    ct_long9$Mutant <- factor(ct_long9$Mutant, levels = c("Eftud2", "Sf3b4"))

    mean_eft9 <- mean(ct_merge$n_sig_genes_eft, na.rm = TRUE)
    mean_sf39 <- mean(ct_merge$n_sig_genes_sf3, na.rm = TRUE)

    p9 <- ggplot(ct_long9, aes(x = cell_type, y = n_genes, fill = Mutant)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7) +
      geom_text(
        aes(label = ifelse((Mutant == "Eftud2" & n_genes > mean_eft9 + 5) |
                             (Mutant == "Sf3b4"  & n_genes > mean_sf39 + 5),
                           as.character(n_genes), "")),
        position = position_dodge(width = 0.75),
        vjust = -0.4, size = 3.2, fontface = "bold"
      ) +
      geom_hline(yintercept = mean_eft9, colour = .GOLD, linewidth = 1.2,
                 linetype = "dashed", alpha = 0.8) +
      geom_hline(yintercept = mean_sf39, colour = .RED,  linewidth = 1.2,
                 linetype = "dashed", alpha = 0.8) +
      annotate("text", x = 0.6, y = mean_eft9 + 1.5,
               label = paste("Eftud2 mean:", round(mean_eft9, 1)),
               hjust = 0, size = 3.2, colour = .GOLD, fontface = "italic") +
      annotate("text", x = 0.6, y = mean_sf39 + 1.5,
               label = paste("Sf3b4 mean:", round(mean_sf39, 1)),
               hjust = 0, size = 3.2, colour = .RED, fontface = "italic") +
      scale_fill_manual(values = c(Eftud2 = .GOLD, Sf3b4 = .RED)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(
        title   = "Differentially Altered Velocity Genes per Cell Type",
        x       = NULL,
        y       = "# Significant Velocity Genes",
        fill    = "Mutant",
        caption = "FDR < 0.05, Wilcoxon rank-sum test"
      ) +
      .tv() +
      theme(
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.major.y = element_line(colour = "#EEEEEE"),
        legend.position    = "top"
      )
    save_plot(p9, file.path(pdir, "09_sig_genes_per_celltype"), 20, 9)
  }

  if (!is.null(shared) && !is.null(eft_sig) && !is.null(sf3_sig) &&
      all(c("gene") %in% colnames(shared)) &&
      all(c("gene", "delta_velocity") %in% colnames(eft_sig)) &&
      all(c("gene", "delta_velocity") %in% colnames(sf3_sig))) {
    message("    Plot 10: Shared genes pheatmap (top 30 by mean |ΔVelocity|)")

    shared_genes <- shared$gene
    eft_sub <- eft_sig[eft_sig$gene %in% shared_genes, c("gene", "delta_velocity")]
    sf3_sub <- sf3_sig[sf3_sig$gene %in% shared_genes, c("gene", "delta_velocity")]
    names(eft_sub)[2] <- "Eftud2"
    names(sf3_sub)[2] <- "Sf3b4"

    combined10 <- merge(eft_sub, sf3_sub, by = "gene", all = FALSE)
    combined10$abs_mean <- (abs(combined10$Eftud2) + abs(combined10$Sf3b4)) / 2
    combined10 <- combined10[order(-combined10$abs_mean), ]
    combined10 <- combined10[seq_len(min(30L, nrow(combined10))), ]
    combined10 <- combined10[, c("gene", "Eftud2", "Sf3b4")]

    if (nrow(combined10) >= 2) {
      mat10 <- as.matrix(combined10[, c("Eftud2", "Sf3b4")])
      rownames(mat10) <- combined10$gene
      vmax10 <- max(abs(mat10), na.rm = TRUE)

      ph10 <- pheatmap::pheatmap(
        mat10,
        color           = colorRampPalette(c(.TEAL, "#F5F5F5", .RED))(101),
        breaks          = seq(-vmax10, vmax10, length.out = 102),
        cluster_rows    = FALSE,
        cluster_cols    = FALSE,
        display_numbers = TRUE,
        number_format   = "%.4f",
        number_color    = "black",
        fontsize        = 11,
        fontsize_row    = 11,
        fontsize_col    = 13,
        main            = "Top 30 Shared Differential Velocity Genes\n(Delta Velocity relative to Control)",
        border_color    = "white",
        angle_col       = 0,
        silent          = TRUE
      )
      save_pheatmap(ph10,
                    file.path(pdir, "10_shared_genes_deltavelocity_heatmap"), 2100, 3600)
    }
  }

  message("  [custom_plots] Done. Outputs in: ", pdir)
  invisible(NULL)
}


compute_replicate_composition <- function(rep_objects) {
  rows <- lapply(names(rep_objects), function(sname) {
    seu  <- rep_objects[[sname]]
    meta <- SAMPLES[[sname]]
    if (is.null(meta) || !"cell_type" %in% colnames(seu@meta.data)) return(NULL)
    ct_vec <- seu$cell_type[seu$cell_type != "Unassigned"]
    if (length(ct_vec) == 0) return(NULL)
    tab <- table(ct_vec)
    data.frame(sample      = sname,
               label       = meta$label,
               condition   = meta$condition,
               replicate   = meta$replicate,
               cell_type   = names(tab),
               n_cells     = as.integer(tab),
               total_cells = sum(tab),
               stringsAsFactors = FALSE)
  })
  comp_df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(comp_df) || nrow(comp_df) == 0) return(NULL)
  comp_df$prop <- comp_df$n_cells / comp_df$total_cells
  comp_df <- do.call(rbind, lapply(split(comp_df, comp_df$sample), function(df) {
    pseudo <- 0.5 / df$total_cells[1]
    p_adj  <- df$prop + pseudo
    df$clr <- log(p_adj) - mean(log(p_adj))
    df
  }))
  comp_df$condition <- factor(comp_df$condition, levels = CONDITIONS)
  comp_df
}

run_composition_stats <- function(rep_objects, cond_objects, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("  Compositional statistics (CLR-LME + Fisher reference) ...")
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for replicate-aware composition tests.")

  comp_df <- compute_replicate_composition(rep_objects)
  if (is.null(comp_df)) { message("    No composition data  -  skipping."); return(invisible(NULL)) }
  all_cts <- sort(unique(comp_df$cell_type))

  message("    CLR-LME omnibus (", length(all_cts), " cell types) ...")
  omnibus_rows <- lapply(all_cts, function(ct) {
    ct_sub <- comp_df[comp_df$cell_type == ct, ]
    if (length(unique(ct_sub$condition)) < 2) return(NULL)
    suppressWarnings(tryCatch({
      f1  <- lme4::lmer(clr ~ condition + (1 | replicate), data = ct_sub, REML = FALSE)
      f0  <- lme4::lmer(clr ~ 1         + (1 | replicate), data = ct_sub, REML = FALSE)
      lrt <- anova(f0, f1)
      data.frame(cell_type = ct, LRT_chi2 = lrt$Chisq[2], df = lrt$Df[2],
                 pvalue = lrt[["Pr(>Chisq)"]][2], stringsAsFactors = FALSE)
    }, error = function(e) NULL))
  })
  omnibus_df <- do.call(rbind, Filter(Negate(is.null), omnibus_rows))
  if (!is.null(omnibus_df) && nrow(omnibus_df) > 0) {
    omnibus_df$padj <- p.adjust(omnibus_df$pvalue, method = "BH")
    omnibus_df      <- omnibus_df[order(omnibus_df$pvalue), ]
    write.csv(omnibus_df, file.path(outdir, "composition_clr_lme_omnibus.csv"), row.names = FALSE)
    message("    CLR-LME omnibus: ", sum(!is.na(omnibus_df$padj) & omnibus_df$padj < 0.05),
            " cell types significant (FDR<0.05)")
  }

  message("    CLR-LME pairwise ...")
  pair_rows_all <- lapply(COMPARISONS, function(comp) {
    trt <- comp[1]; ref <- comp[2]
    sub <- comp_df[comp_df$condition %in% c(trt, ref), ]
    lapply(all_cts, function(ct) {
      ct_sub <- sub[sub$cell_type == ct, ]
      if (nrow(ct_sub) < 4 || length(unique(ct_sub$condition)) < 2) return(NULL)
      ct_sub$cond_bin <- ct_sub$condition == trt
      suppressWarnings(tryCatch({
        f1  <- lme4::lmer(clr ~ cond_bin + (1 | replicate), data = ct_sub, REML = FALSE)
        f0  <- lme4::lmer(clr ~ 1        + (1 | replicate), data = ct_sub, REML = FALSE)
        lrt <- anova(f0, f1)
        sc  <- summary(f1)$coefficients
        idx <- grep("cond_bin", rownames(sc))[1]
        data.frame(
          comparison   = paste0(trt, "_vs_", ref),
          cell_type    = ct,
          clr_effect   = sc[idx, "Estimate"],
          clr_se       = sc[idx, "Std. Error"],
          delta_clr    = mean(ct_sub$clr[ct_sub$condition == trt]) -
                         mean(ct_sub$clr[ct_sub$condition == ref]),
          LRT_chi2     = lrt$Chisq[2],
          pvalue       = lrt[["Pr(>Chisq)"]][2],
          mean_pct_trt = mean(ct_sub$prop[ct_sub$condition == trt]) * 100,
          mean_pct_ref = mean(ct_sub$prop[ct_sub$condition == ref]) * 100,
          stringsAsFactors = FALSE)
      }, error = function(e) NULL))
    })
  })
  lme_pair_df <- do.call(rbind, Filter(Negate(is.null),
                                       unlist(pair_rows_all, recursive = FALSE)))
  if (!is.null(lme_pair_df) && nrow(lme_pair_df) > 0) {
    lme_pair_df <- do.call(rbind, lapply(split(lme_pair_df, lme_pair_df$comparison), function(df) {
      df$padj <- p.adjust(df$pvalue, method = "BH"); df
    }))
    lme_pair_df <- lme_pair_df[order(lme_pair_df$comparison, lme_pair_df$pvalue), ]
    write.csv(lme_pair_df, file.path(outdir, "composition_clr_lme_pairwise.csv"), row.names = FALSE)
    message("    CLR-LME pairwise: ", sum(!is.na(lme_pair_df$padj) & lme_pair_df$padj < 0.05),
            " significant cell-type shifts (FDR<0.05)")

    lme_pair_df$sig_lbl    <- ifelse(!is.na(lme_pair_df$padj) & lme_pair_df$padj < 0.01, "**",
                               ifelse(!is.na(lme_pair_df$padj) & lme_pair_df$padj < 0.05, "*", ""))
    comp_levels            <- vapply(COMPARISONS, paste, character(1), collapse = "_vs_")
    lme_pair_df$comparison <- factor(lme_pair_df$comparison, levels = comp_levels)
    ct_ord <- if (!is.null(omnibus_df)) omnibus_df$cell_type else sort(all_cts)
    lme_pair_df$cell_type  <- factor(lme_pair_df$cell_type,
                                     levels = intersect(ct_ord, unique(lme_pair_df$cell_type)))

    p_forest <- ggplot(lme_pair_df,
                       aes(x = clr_effect, y = cell_type,
                           xmin = clr_effect - 1.96 * clr_se,
                           xmax = clr_effect + 1.96 * clr_se,
                           colour = comparison)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(height = 0.25, linewidth = 0.7, position = position_dodge(0.6)) +
      geom_point(size = 2.2, position = position_dodge(0.6)) +
      geom_text(aes(x = clr_effect + 1.96 * clr_se + 0.02, label = sig_lbl),
                position = position_dodge(0.6), size = 3.5, hjust = 0, colour = "black") +
      scale_colour_manual(values = setNames(AAAS_COLS[c(2, 3, 5)], levels(lme_pair_df$comparison))) +
      labs(x = "CLR effect (log-ratio units, 95% CI)", y = NULL,
           title    = "Cell-type composition  -  CLR-LME pairwise effects",
           subtitle = "Points = LME estimate +/- 95% CI  |  * FDR<0.05  ** FDR<0.01",
           colour   = "Comparison") +
      .theme_velo(base_size = 12) + theme(legend.position = "bottom")
    save_plot(p_forest, file.path(outdir, "composition_clr_lme_forest"),
              16, max(9, length(all_cts) * 0.38), limitsize = FALSE)

    for (comp_lbl in levels(lme_pair_df$comparison)) {
      fd <- lme_pair_df[lme_pair_df$comparison == comp_lbl, ]
      if (nrow(fd) < 2) next
      fd$neg_logp <- -log10(fd$pvalue + 1e-300)
      fd$n_mean   <- (fd$mean_pct_trt + fd$mean_pct_ref) / 2
      col_parts   <- strsplit(as.character(comp_lbl), "_vs_")[[1]]
      col_q <- .pal_get(COND_PALETTE, col_parts[1]) %||% AAAS_COLS[2]
      col_r <- .pal_get(COND_PALETTE, col_parts[2]) %||% AAAS_COLS[1]
      p_bub <- ggplot(fd, aes(y = reorder(cell_type, neg_logp), x = neg_logp,
                               size = n_mean, fill = delta_clr)) +
        geom_point(shape = 21, colour = "black", stroke = 0.45, alpha = 0.88) +
        geom_text(aes(label = sig_lbl), size = 4, hjust = -0.5,
                  colour = "black", fontface = "bold") +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                   colour = "red", linewidth = 0.5) +
        scale_size_continuous(range = c(3, 14), name = "Mean % cells") +
        scale_fill_gradient2(low = col_r, mid = "white", high = col_q,
                             midpoint = 0, name = "Delta CLR") +
        labs(x = expression(-log[10] * "(LRT p-value)"), y = NULL,
             title    = paste0("Cell-type composition  -  ", gsub("_", " ", comp_lbl)),
             subtitle = "Bubble size = mean % cells; fill = Delta CLR; * FDR<0.05; ** FDR<0.01") +
        .theme_velo(base_size = 12) + theme(legend.position = "right")
      save_plot(p_bub, file.path(outdir, paste0(comp_lbl, "_by_celltype_bubble")),
                PLT$w_single, max(8, nrow(fd) * 0.4), limitsize = FALSE)
    }
  }

  message("    Chi-squared / Fisher tests (pooled-replicate reference) ...")
  meta_all <- do.call(rbind, lapply(seq_along(cond_objects), function(i) {
    seu <- cond_objects[[i]]
    if (!"cell_type" %in% colnames(seu@meta.data)) return(NULL)
    data.frame(group = names(cond_objects)[i], cell_type = seu$cell_type,
               stringsAsFactors = FALSE)
  }))
  if (!is.null(meta_all) && nrow(meta_all) > 0) {
    meta_all <- meta_all[!is.na(meta_all$cell_type) & meta_all$cell_type != "Unassigned", ]
  }
  if (!is.null(meta_all) && nrow(meta_all) > 0) {
    ct_table <- table(meta_all$cell_type, meta_all$group)
    chi_test <- tryCatch(chisq.test(ct_table), error = function(e) NULL)
    if (!is.null(chi_test)) {
      write.csv(
        data.frame(test = "Chi-squared (pooled replicates)",
                   X2 = chi_test$statistic, df = chi_test$parameter,
                   pvalue = chi_test$p.value, stringsAsFactors = FALSE),
        file.path(outdir, "composition_chisq_test.csv"), row.names = FALSE)
      message("    Chi-squared (pooled): X2=", round(chi_test$statistic, 2),
              "  p=", formatC(chi_test$p.value, format = "e", digits = 3))
    }
    pair_rows <- lapply(COMPARISONS, function(comp) {
      g1 <- comp[1]; g2 <- comp[2]
      if (!g1 %in% colnames(ct_table) || !g2 %in% colnames(ct_table)) return(NULL)
      lapply(rownames(ct_table), function(ct) {
        n1 <- ct_table[ct, g1]; n2 <- ct_table[ct, g2]
        o1 <- sum(ct_table[, g1]) - n1; o2 <- sum(ct_table[, g2]) - n2
        if (any(c(n1, n2, o1, o2) < 1)) return(NULL)
        ft <- tryCatch(fisher.test(matrix(c(n1, o1, n2, o2), nrow = 2)),
                       error = function(e) NULL)
        if (is.null(ft)) return(NULL)
        data.frame(comparison = paste0(g1, "_vs_", g2), cell_type = ct,
                   OR = ft$estimate, pvalue = ft$p.value,
                   pct_g1 = n1 / sum(ct_table[, g1]) * 100,
                   pct_g2 = n2 / sum(ct_table[, g2]) * 100,
                   n_g1 = n1, n_g2 = n2, stringsAsFactors = FALSE)
      })
    })
    fisher_df <- do.call(rbind, Filter(Negate(is.null), unlist(pair_rows, recursive = FALSE)))
    if (!is.null(fisher_df) && nrow(fisher_df) > 0) {
      fisher_df$padj   <- p.adjust(fisher_df$pvalue, method = "BH")
      fisher_df$log2OR <- log2(fisher_df$OR + 1e-6)
      fisher_df        <- fisher_df[order(fisher_df$pvalue), ]
      write.csv(fisher_df, file.path(outdir, "composition_fisher_tests.csv"), row.names = FALSE)
      message("    Fisher (pooled): ", sum(!is.na(fisher_df$padj) & fisher_df$padj < 0.05),
              " significant shifts (FDR<0.05)  -  exploratory only")

      all_cts2 <- sort(unique(fisher_df$cell_type))
      all_cmp  <- unique(fisher_df$comparison)
      if (length(all_cts2) >= 2 && length(all_cmp) >= 2) {
        or_mat <- do.call(cbind, lapply(all_cmp, function(comp) {
          fd <- fisher_df[fisher_df$comparison == comp, ]
          v  <- setNames(fd$log2OR, fd$cell_type); v[all_cts2]
        }))
        colnames(or_mat) <- all_cmp; rownames(or_mat) <- all_cts2

        # ── Sanitise or_mat before clustering ─────────────────────────────────
        # fisher.test() can return OR = Inf (complete separation), making
        # log2OR = Inf.  v[all_cts2] produces NA for cell types absent from a
        # given comparison.  Both cause dist() -> hclust() to crash with
        # "NA/NaN/Inf in foreign function call (arg 10)".
        #
        # Fix: (1) replace NaN/NA with 0 (log2OR of 1 = no enrichment),
        #       (2) cap +/-Inf to the largest observed finite |log2OR|,
        #       (3) disable row clustering if the matrix is still degenerate
        #           (e.g. single unique row or all-zero).
        or_mat[is.nan(or_mat)] <- 0
        or_mat[is.na(or_mat)]  <- 0
        finite_vals <- or_mat[is.finite(or_mat)]
        cap_val     <- if (length(finite_vals) && max(abs(finite_vals)) > 0)
                         max(abs(finite_vals)) else 1
        or_mat[or_mat ==  Inf] <-  cap_val
        or_mat[or_mat == -Inf] <- -cap_val
        # Safe to cluster only if ≥2 rows, all finite, and rows are not identical
        can_cluster_rows <- nrow(or_mat) >= 2 &&
                            all(is.finite(or_mat)) &&
                            (nrow(unique(or_mat)) > 1)
        if (!can_cluster_rows)
          message("    NOTE: OR heatmap row clustering disabled (degenerate matrix).")

        lim_or <- max(max(abs(or_mat), na.rm = TRUE) * 0.95, 0.1)
        h_or <- max(1200, length(all_cts2) * 60); w_or <- max(900, length(all_cmp) * 120)
        ph_or <- pheatmap(or_mat,
                          color        = colorRampPalette(c(AAAS_COLS[1], "white", AAAS_COLS[2]))(100),
                          breaks       = seq(-lim_or, lim_or, length.out = 101),
                          cluster_rows = can_cluster_rows, cluster_cols = FALSE,
                          fontsize      = 7,
                          fontsize_row = max(5, min(9, floor(h_or / length(all_cts2) / PLT$dpi * 72 * 0.5))),
                          fontsize_col = 8, na_col = "grey90",
                          main         = "Cell-type composition log2(OR)  -  Fisher (pooled replicates)",
                          silent       = TRUE)
        save_pheatmap(ph_or, file.path(outdir, "composition_OR_heatmap"), w_or * 2, h_or)
      }
    }
  }
  message("  Compositional statistics complete.")
}

plot_composition_ranked_bars <- function(rep_objects, cond_objects, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("  Composition ranked bar plots ...")

  comp_df <- compute_replicate_composition(rep_objects)
  if (is.null(comp_df)) { message("    No composition data  -  skipping."); return(invisible(NULL)) }

  sum_df <- do.call(rbind, by(comp_df, list(comp_df$condition, comp_df$cell_type), function(df) {
    data.frame(condition = df$condition[1], cell_type = df$cell_type[1],
               mean_pct  = mean(df$prop) * 100,
               se_pct    = sd(df$prop) / sqrt(nrow(df)) * 100,
               stringsAsFactors = FALSE)
  }))
  sum_df$condition <- factor(sum_df$condition, levels = CONDITIONS)

  ct_order_all <- sum_df %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(m = mean(mean_pct), .groups = "drop") %>%
    dplyr::arrange(m) %>% dplyr::pull(cell_type)

  cond_rep_map <- lapply(setNames(CONDITIONS, CONDITIONS), function(cond) {
    names(Filter(function(s) s$condition == cond, SAMPLES))
  })
  for (cond in CONDITIONS) {
    avail <- intersect(cond_rep_map[[cond]], names(rep_objects))
    if (length(avail) < 2) next
    rep_labels <- vapply(avail, function(s) SAMPLES[[s]]$label, character(1))
    df_rep     <- comp_df[comp_df$sample %in% avail, ]
    ct_ord_rep <- sum_df$cell_type[sum_df$condition == cond][
      order(sum_df$mean_pct[sum_df$condition == cond])]
    df_rep$cell_type <- factor(df_rep$cell_type, levels = ct_ord_rep)
    df_rep$label     <- factor(df_rep$label,     levels = rev(rep_labels))
    pal_rep          <- setNames(c(AAAS_COLS[4], AAAS_COLS[5]), rep_labels)
    p_rep <- ggplot(df_rep, aes(y = cell_type, x = prop * 100, fill = label)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.88) +
      geom_label(aes(label = sprintf("%.1f%%", prop * 100)),
                 position = position_dodge(width = 0.75),
                 hjust = -0.05, vjust = 0.5, size = 3,
                 fontface = "bold", colour = "black",
                 fill = scales::alpha("white", 0.75), label.size = 0) +
      scale_fill_manual(values = pal_rep) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
      labs(x = "% cells (per replicate)", y = NULL,
           title    = paste0("Cell-type composition  -  ", cond, " replicates"),
           subtitle = paste(rep_labels, collapse = " vs "),
           fill = "Replicate") +
      .theme_velo(base_size = 11) + theme(legend.position = "bottom")
    save_plot(p_rep,
              file.path(outdir, paste0(cond, "_replicate_composition_ranked_bar")),
              16, max(8, length(ct_ord_rep) * 0.42), limitsize = FALSE)
  }

  sum_df$cell_type <- factor(sum_df$cell_type, levels = ct_order_all)
  h_cond <- max(9, length(ct_order_all) * 0.42)

  p_cond <- ggplot(sum_df, aes(y = cell_type, x = mean_pct, fill = condition)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72, alpha = 0.88) +
    geom_errorbarh(aes(xmin = pmax(0, mean_pct - se_pct), xmax = mean_pct + se_pct),
                   position = position_dodge(width = 0.78),
                   height = 0.3, linewidth = 0.5, colour = "grey30") +
    geom_label(aes(x = mean_pct + se_pct, label = sprintf("%.1f%%", mean_pct)),
               position = position_dodge(width = 0.78),
               hjust = -0.05, vjust = 0.5, size = 2.8,
               fontface = "bold", colour = "black",
               fill = scales::alpha("white", 0.75), label.size = 0) +
    scale_fill_manual(values = COND_PALETTE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.26))) +
    labs(x = "Mean % cells +/- SE (across replicates)", y = NULL,
         title    = "Cell-type composition across conditions",
         subtitle = paste0(paste(CONDITIONS, collapse = " vs "),
                           "  |  bars = mean, whiskers = SE"),
         fill = "Condition") +
    .theme_velo(base_size = 11) + theme(legend.position = "bottom")
  save_plot(p_cond, file.path(outdir, "conditions_composition_ranked_bar"),
            18, h_cond, limitsize = FALSE)

  comp_df$cell_type <- factor(comp_df$cell_type, levels = ct_order_all)
  comp_df$condition <- factor(comp_df$condition, levels = CONDITIONS)

  p_dot <- ggplot(sum_df, aes(y = cell_type, x = mean_pct, colour = condition)) +
    geom_errorbarh(aes(xmin = pmax(0, mean_pct - se_pct), xmax = mean_pct + se_pct),
                   position = position_dodge(width = 0.65), height = 0.2, linewidth = 0.6) +
    geom_point(size = 4, position = position_dodge(width = 0.65), shape = 18) +
    geom_point(data = comp_df, aes(y = cell_type, x = prop * 100, colour = condition),
               position = position_dodge(width = 0.65),
               size = 1.8, alpha = 0.55, shape = 16) +
    scale_colour_manual(values = COND_PALETTE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "% cells", y = NULL,
         title    = "Cell-type composition  -  replicate scatter with condition mean +/- SE",
         subtitle = "Diamonds = condition mean; dots = individual replicates",
         colour   = "Condition") +
    .theme_velo(base_size = 11) + theme(legend.position = "bottom")
  save_plot(p_dot, file.path(outdir, "conditions_composition_dotplot"),
            16, h_cond, limitsize = FALSE)

  for (comp in COMPARISONS) {
    g1 <- comp[1]; g2 <- comp[2]
    df2 <- sum_df[sum_df$condition %in% c(g1, g2), ]
    if (nrow(df2) < 2) next
    ct_ord_p <- df2 %>% dplyr::group_by(cell_type) %>%
      dplyr::summarise(m = mean(mean_pct), .groups = "drop") %>%
      dplyr::arrange(m) %>% dplyr::pull(cell_type)
    df2$cell_type <- factor(df2$cell_type, levels = ct_ord_p)
    df2$condition <- factor(df2$condition, levels = c(g1, g2))
    lbl <- paste0(g1, "_vs_", g2)
    p_pair <- ggplot(df2, aes(y = cell_type, x = mean_pct, fill = condition)) +
      geom_col(position = position_dodge(width = 0.76), width = 0.70, alpha = 0.88) +
      geom_errorbarh(aes(xmin = pmax(0, mean_pct - se_pct), xmax = mean_pct + se_pct),
                     position = position_dodge(width = 0.76),
                     height = 0.3, linewidth = 0.5, colour = "grey30") +
      geom_label(aes(x = mean_pct + se_pct, label = sprintf("%.1f%%", mean_pct)),
                 position = position_dodge(width = 0.76),
                 hjust = -0.05, vjust = 0.5, size = 3,
                 fontface = "bold", colour = "black",
                 fill = scales::alpha("white", 0.75), label.size = 0) +
      scale_fill_manual(values = COND_PALETTE[c(g1, g2)]) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
      labs(x = "Mean % cells +/- SE (across replicates)", y = NULL,
           title = paste0("Cell-type composition  -  ", g1, " vs ", g2),
           fill  = "Condition") +
      .theme_velo(base_size = 11) + theme(legend.position = "bottom")
    save_plot(p_pair, file.path(outdir, paste0(lbl, "_composition_ranked_bar")),
              16, max(8, length(ct_ord_p) * 0.42), limitsize = FALSE)
  }
  message("    Composition ranked bar plots saved.")
}

run_velocity_expression_integration <- function(seu_q, seu_r,
                                                 vel_de_df, de_df,
                                                 query_label, ref_label, outdir) {
  message("  Velocity-expression integration (discordant genes): ",
          query_label, " vs ", ref_label)
  if (is.null(vel_de_df) || is.null(de_df)) return(invisible(NULL))
  if (nrow(vel_de_df) == 0 || nrow(de_df) == 0) return(invisible(NULL))

  # Ensure gene column in both
  if (!"gene" %in% colnames(vel_de_df)) vel_de_df$gene <- rownames(vel_de_df)
  if (!"gene" %in% colnames(de_df))     de_df$gene     <- rownames(de_df)

  merged <- merge(vel_de_df[, c("gene","delta_velocity","padj")],
                  de_df[, c("gene","avg_log2FC","p_val_adj")],
                  by = "gene", suffixes = c("_vel","_expr"))
  if (nrow(merged) < 5) return(invisible(NULL))

  merged$sig_vel  <- !is.na(merged$padj)       & merged$padj       < 0.05
  merged$sig_expr <- !is.na(merged$p_val_adj)  & merged$p_val_adj  < 0.1

  # Discordant: velocity and expression in opposite directions (both significant)
  merged$direction_vel  <- sign(merged$delta_velocity)
  merged$direction_expr <- sign(merged$avg_log2FC)
  merged$discordant <- merged$sig_vel & merged$sig_expr &
                       (merged$direction_vel != merged$direction_expr)
  merged$concordant <- merged$sig_vel & merged$sig_expr &
                       (merged$direction_vel == merged$direction_expr)

  lbl <- paste0(query_label, "_vs_", ref_label)
  write.csv(merged, file.path(outdir, paste0(lbl, "_velocity_expression_integration.csv")),
            row.names = FALSE)

  merged$category <- "Other"
  merged$category[merged$discordant] <- "Discordant"
  merged$category[merged$concordant] <- "Concordant"
  cat_cols <- c(Discordant = AAAS_COLS[2], Concordant = AAAS_COLS[3], Other = "grey80")

  disc_genes <- merged[merged$discordant, ]
  conc_genes <- merged[merged$concordant, ]
  message("    Discordant: ", nrow(disc_genes), "  Concordant: ", nrow(conc_genes))

  top_disc <- head(disc_genes[order(-abs(disc_genes$delta_velocity)), ], 10)
  p_integ <- ggplot(merged, aes(avg_log2FC, delta_velocity, colour = category)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_label_repel(data = top_disc, aes(label = gene),
                              size = 3, show.legend = FALSE) +
    scale_colour_manual(values = cat_cols) +
    labs(x = "log2FC (expression)", y = "Delta Velocity",
         title = paste("Velocity-Expression integration:", lbl),
         subtitle = paste0("Discordant: ", nrow(disc_genes),
                           "  Concordant: ", nrow(conc_genes)),
         colour = NULL) +
    .theme_velo(base_size = 13)
  save_plot(p_integ, file.path(outdir, paste0(lbl, "_velocity_expression_scatter")),
            PLT$w_single, PLT$h_single)
  message("    Integration plot saved.")

  if (nrow(disc_genes) >= 5 &&
      requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace("org.Mm.eg.db",    quietly = TRUE)) {
    message("    Running GO enrichment on ", nrow(disc_genes), " discordant genes ...")
    all_ent  <- .to_entrez(merged$gene)
    disc_ent <- .to_entrez(disc_genes$gene)
    if (length(na.omit(disc_ent)) >= 5) {
      disc_enr_dir <- file.path(outdir, paste0(lbl, "_discordant_GO"))
      dir.create(disc_enr_dir, recursive = TRUE, showWarnings = FALSE)
      for (ont in c("BP", "MF", "CC")) {
        tryCatch({
          ego_disc <- clusterProfiler::enrichGO(
            gene          = na.omit(disc_ent),
            universe      = na.omit(all_ent),
            OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
            keyType       = "ENTREZID", ont = ont,
            pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,
            pvalueCutoff  = 0.1, readable = TRUE)
          df_disc_go <- as.data.frame(ego_disc)
          if (nrow(df_disc_go) > 0) {
            write.csv(df_disc_go,
                      file.path(disc_enr_dir, paste0(lbl, "_discordant_GO_", ont, ".csv")),
                      row.names = FALSE)
            .save_enrich_plots_v(ego_disc,
                                  file.path(disc_enr_dir, paste0(lbl, "_discordant_GO_", ont)),
                                  disc_enr_dir,
                                  paste0("GO-", ont, " | velocity-discordant genes | ", lbl))
            message("      GO-", ont, ": ", nrow(df_disc_go), " enriched terms")
          }
        }, error = function(e) message("      GO-", ont, " failed: ", e$message))
      }
    } else {
      message("    Insufficient Entrez IDs for discordant GO  -  skipping")
    }
  }

  merged
}

.pseudobulk_velocity_lm <- function(seu_q, seu_r, query_label, ref_label,
                                     rep_objects, outdir) {
  message("  Pseudobulk-LM velocity DE: ", query_label, " vs ", ref_label)
  lbl <- paste0(query_label, "_vs_", ref_label)

  .pb_vel <- function(seu, cond) {
    if (!"velocity" %in% names(seu@assays)) return(NULL)
    vm   <- .get_layer(seu, "velocity", "counts")
    reps <- if ("replicate" %in% colnames(seu@meta.data)) seu$replicate else rep("rep1", ncol(seu))
    do.call(cbind, lapply(unique(reps), function(r) {
      cells_r <- colnames(seu)[reps == r]
      if (length(cells_r) == 0) return(NULL)
      col_nm <- paste0(cond, "_", r)
      vm_r   <- Matrix::rowMeans(vm[, cells_r, drop = FALSE], na.rm = TRUE)
      setNames(data.frame(vm_r), col_nm)
    }))
  }
  pb_q <- .pb_vel(seu_q, query_label)
  pb_r <- .pb_vel(seu_r, ref_label)
  if (is.null(pb_q) || is.null(pb_r)) {
    message("    Could not build pseudobulk velocity matrices  -  skipping"); return(NULL)
  }

  cg    <- intersect(rownames(seu_q@assays$velocity), rownames(seu_r@assays$velocity))
  pb_q  <- pb_q[cg, , drop = FALSE]
  pb_r  <- pb_r[cg, , drop = FALSE]
  pb    <- as.matrix(cbind(pb_q, pb_r))
  group <- factor(c(rep(query_label, ncol(pb_q)), rep(ref_label, ncol(pb_r))),
                  levels = c(ref_label, query_label))
  design <- model.matrix(~ group)

  if (ncol(pb) < 3) {
    message("    Too few pseudobulk samples for limma  -  skipping"); return(NULL)
  }

  tryCatch({
    fit  <- limma::lmFit(pb, design)
    fit  <- limma::eBayes(fit)
    res  <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
    res$gene <- rownames(res)
    colnames(res)[colnames(res) == "logFC"]       <- "delta_velocity"
    colnames(res)[colnames(res) == "adj.P.Val"]   <- "padj"
    colnames(res)[colnames(res) == "P.Value"]     <- "pvalue"
    write.csv(res, file.path(outdir, paste0(lbl, "_pseudobulk_velocity_LM.csv")),
              row.names = FALSE)
    n_sig <- sum(!is.na(res$padj) & res$padj < 0.05)
    message("    Pseudobulk-LM: ", n_sig, " sig genes (FDR<0.05)")
    res
  }, error = function(e) {
    message("    Pseudobulk-LM failed: ", e$message); NULL
  })
}


plot_overview <- function(seu, label, outdir, extra_groups = character(0)) {
  ct_palette <- .build_ct_palette(unique(seu$cell_type))
  p_cluster <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
    pt.size = PLT$pt_size, label.size = PLT$label_size) +
    ggtitle(paste(label, " -  Clusters")) + NoLegend()
  p_ct <- DimPlot(seu, group.by = "cell_type", cols = ct_palette, pt.size = PLT$pt_size,
    label = TRUE, repel = TRUE, label.size = PLT$label_size) +
    ggtitle(paste(label, " -  Cell type")) +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  has_speed <- "velocity_length"     %in% colnames(seu@meta.data)
  has_coh   <- "velocity_confidence" %in% colnames(seu@meta.data)
  has_lt    <- "latent_time"         %in% colnames(seu@meta.data)
  has_pt    <- "velocity_pseudotime" %in% colnames(seu@meta.data)
  p_speed <- if (has_speed) FeaturePlot(seu, features = "velocity_length",
    cols = AAAS_SPEED_COLS, pt.size = PLT$pt_size) + ggtitle("Speed") else NULL
  p_coh   <- if (has_coh)   FeaturePlot(seu, features = "velocity_confidence",
    cols = AAAS_COHER_COLS, pt.size = PLT$pt_size) + ggtitle("Coherence") else NULL
  p_lt    <- if (has_lt)    FeaturePlot(seu, features = "latent_time",
    cols = colorRampPalette(rev(AAAS_LATENT_COLS))(50), pt.size = PLT$pt_size) + ggtitle("Latent time") else NULL
  p_pt    <- if (has_pt)    FeaturePlot(seu, features = "velocity_pseudotime",
    cols = AAAS_PSEUDO_COLS, pt.size = PLT$pt_size) + ggtitle("Velocity pseudotime") else NULL
  row1_plots <- list(p_cluster, p_ct)
  if (!is.null(p_speed)) row1_plots <- c(row1_plots, list(p_speed))
  row2_real  <- Filter(Negate(is.null), list(p_coh, p_lt, p_pt))
  row2_plots <- c(row2_real, rep(list(plot_spacer()), max(0, 3L - length(row2_real))))
  layout_main <- wrap_plots(row1_plots, ncol = 3L) / wrap_plots(row2_plots, ncol = 3L)
  save_plot(layout_main, file.path(outdir, paste0(label, "_overview_grid")),
            PLT$w_triple * 1.5, PLT$h_double * 2.0)
  save_plot(p_ct, file.path(outdir, paste0(label, "_umap_cell_type")), PLT$w_single, PLT$h_single)
  save_plot(p_cluster + NoLegend(),
            file.path(outdir, paste0(label, "_umap_clusters")), PLT$w_single, PLT$h_single)
  if (!is.null(p_lt))
    save_plot(p_lt, file.path(outdir, paste0(label, "_latent_time")), PLT$w_single, PLT$h_single)
  for (grp in extra_groups) {
    if (!grp %in% colnames(seu@meta.data)) next
    p_grp <- DimPlot(seu, group.by = grp, pt.size = PLT$pt_size,
      cols = if (grp == "condition") COND_PALETTE else if (grp == "replicate") REP_PALETTE else NULL) +
      ggtitle(paste(label, " - ", grp))
    save_plot(p_grp, file.path(outdir, paste0(label, "_umap_", grp)), PLT$w_single, PLT$h_single)
  }
  message("  Plots saved: ", outdir, "/", label, "_*")
}

plot_ct_composition <- function(seu, label, outdir) {
  if (!"cell_type" %in% colnames(seu@meta.data)) return(invisible(NULL))
  ct_df <- as.data.frame(table(cell_type = seu$cell_type)) %>%
    dplyr::mutate(pct = Freq / sum(Freq) * 100) %>% dplyr::arrange(dplyr::desc(pct))
  ct_palette <- .build_ct_palette(ct_df$cell_type)
  dims <- .auto_dims_v(nrow(ct_df), base_w = 14, base_h = 9, per_item_w = 0.45)
  p <- ggplot(ct_df, aes(x = reorder(cell_type, -pct), y = pct, fill = cell_type)) +
    geom_col(show.legend = FALSE) + scale_fill_manual(values = ct_palette) +
    labs(x = NULL, y = "% cells", title = paste("Cell-type composition  - ", label)) +
    .theme_velo() + theme(axis.text.x = element_text(angle = 55, hjust = 1))
  save_plot(p, file.path(outdir, paste0(label, "_celltype_composition")), dims$w, dims$h)
}

plot_velocity_arrows <- function(seu, velo_out, sce_full, label, outdir) {
  if (is.null(velo_out)) { message("  Arrow plots skipped (no velocity)"); return(invisible(NULL)) }
  if (!"umap" %in% names(seu@reductions)) { message("  Arrow plots skipped (no UMAP)"); return(invisible(NULL)) }
  grid.df <- tryCatch({
    scv <- reticulate::import("scvelo", convert = FALSE)
    scv$tl$velocity_embedding(velo_out, basis = "umap")
    vel_umap  <- reticulate::py_to_r(velo_out$obsm[["velocity_umap"]])
    adata_obs <- reticulate::py_to_r(velo_out$obs_names$tolist())
    umap_mat  <- Embeddings(seu, "umap")
    common_cells <- intersect(rownames(umap_mat), adata_obs)
    if (length(common_cells) < 10L) stop("< 10 common cells")
    idx      <- match(common_cells, adata_obs)
    umap_sub <- umap_mat[common_cells, , drop = FALSE]
    vel_sub  <- vel_umap[idx, , drop = FALSE]
    n_grid <- 20L
    xb <- .bincode(umap_sub[,1], seq(min(umap_sub[,1]), max(umap_sub[,1]), length.out = n_grid+1L), include.lowest = TRUE)
    yb <- .bincode(umap_sub[,2], seq(min(umap_sub[,2]), max(umap_sub[,2]), length.out = n_grid+1L), include.lowest = TRUE)
    gid  <- paste(xb, yb, sep = "_")
    rows <- lapply(unique(gid), function(g) {
      i <- which(gid == g); if (length(i) < 3L) return(NULL)
      data.frame(start.1 = mean(umap_sub[i,1]), start.2 = mean(umap_sub[i,2]),
                 end.1   = mean(umap_sub[i,1]) + mean(vel_sub[i,1]),
                 end.2   = mean(umap_sub[i,2]) + mean(vel_sub[i,2]))
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }, error = function(e) { message("  Arrow grid failed: ", e$message); NULL })
  maybe_arrows <- if (!is.null(grid.df))
    geom_segment(data = grid.df,
                 mapping = aes(x = start.1, y = start.2, xend = end.1, yend = end.2, colour = NULL),
                 arrow = arrow(length = unit(0.12,"inches"), type = "closed", angle = 20),
                 colour = "black", alpha = 0.80, linewidth = 0.8, lineend = "round", linejoin = "mitre")
  else NULL
  ct_palette <- .build_ct_palette(unique(seu$cell_type))
  p_ct <- DimPlot(seu, group.by = "cell_type", cols = ct_palette, pt.size = PLT$pt_size,
    label = TRUE, repel = TRUE, label.size = PLT$label_size) +
    ggtitle(paste(label, " -  velocity (cell type)")) +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  if (!is.null(maybe_arrows)) p_ct <- p_ct + maybe_arrows
  save_plot(p_ct, file.path(outdir, paste0(label, "_velocity_arrows_celltype")), PLT$w_single, PLT$h_single)
  if ("condition" %in% colnames(seu@meta.data) && length(unique(seu$condition)) > 1L) {
    p_cond <- DimPlot(seu, group.by = "condition", cols = COND_PALETTE, pt.size = PLT$pt_size) +
      ggtitle(paste(label, " -  velocity (condition)"))
    if (!is.null(maybe_arrows)) p_cond <- p_cond + maybe_arrows
    save_plot(p_cond, file.path(outdir, paste0(label, "_velocity_arrows_condition")), PLT$w_single, PLT$h_single)
  }
  if ("velocity_pseudotime" %in% colnames(seu@meta.data)) {
    p_pt <- FeaturePlot(seu, features = "velocity_pseudotime", cols = AAAS_PSEUDO_COLS, pt.size = PLT$pt_size) +
      ggtitle(paste(label, " -  velocity pseudotime"))
    if (!is.null(maybe_arrows)) p_pt <- p_pt + maybe_arrows
    save_plot(p_pt, file.path(outdir, paste0(label, "_velocity_pseudotime")), PLT$w_single, PLT$h_single)
  }
  message("  Arrow plots saved.")
}

plot_pca_panels <- function(seu, label, outdir) {
  if (!"pca" %in% names(seu@reductions)) return(invisible(NULL))
  message("  PCA panels: ", label)
  p_elbow <- ElbowPlot(seu, ndims = 40) +
    ggtitle(paste(label, " -  Variance explained per PC")) + theme_classic(base_size = 14)
  save_plot(p_elbow, file.path(outdir, paste0(label, "_pca_elbow")), 10, 7)
  ct_palette <- .build_ct_palette(unique(seu$cell_type))
  p_pca_clust <- DimPlot(seu, reduction = "pca", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, pt.size = PLT$pt_size, label.size = PLT$label_size) +
    ggtitle(paste(label, " -  PCA (clusters)")) + NoLegend()
  p_pca_ct <- DimPlot(seu, reduction = "pca", group.by = "cell_type",
    cols = ct_palette, pt.size = PLT$pt_size, label = TRUE, repel = TRUE,
    label.size = PLT$label_size) + ggtitle(paste(label, " -  PCA (cell type)")) +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  panels <- list(p_pca_clust, p_pca_ct)
  if ("condition" %in% colnames(seu@meta.data) && length(unique(seu$condition)) > 1L) {
    panels <- c(panels, list(DimPlot(seu, reduction = "pca", group.by = "condition",
      cols = COND_PALETTE, pt.size = PLT$pt_size) + ggtitle(paste(label, " -  PCA (condition)"))))
  }
  save_plot(wrap_plots(panels, ncol = min(length(panels), 3L)),
            file.path(outdir, paste0(label, "_pca_dimplots")), PLT$w_triple, PLT$h_single)
  save_plot(VizDimLoadings(seu, dims = 1:2, reduction = "pca", nfeatures = 15, balanced = TRUE) &
              theme_classic(base_size = 12),
            file.path(outdir, paste0(label, "_pca_loadings")), PLT$w_double, PLT$h_single)
}

plot_tsne_panels <- function(seu, label, outdir) {
  if (!"tsne" %in% names(seu@reductions)) return(invisible(NULL))
  message("  t-SNE panels: ", label)
  ct_palette <- .build_ct_palette(unique(seu$cell_type))
  p_clust <- DimPlot(seu, reduction = "tsne", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, pt.size = PLT$pt_size, label.size = PLT$label_size) +
    ggtitle(paste(label, " -  t-SNE (clusters)")) + NoLegend()
  p_ct <- DimPlot(seu, reduction = "tsne", group.by = "cell_type",
    cols = ct_palette, pt.size = PLT$pt_size, label = TRUE, repel = TRUE,
    label.size = PLT$label_size) + ggtitle(paste(label, " -  t-SNE (cell type)")) +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  save_plot(p_ct, file.path(outdir, paste0(label, "_tsne_cell_type")), PLT$w_single, PLT$h_single)
  panels <- list(p_clust, p_ct)
  if ("cluster_majority_type" %in% colnames(seu@meta.data)) {
    maj_pal <- .build_ct_palette(unique(seu$cluster_majority_type))
    p_maj <- DimPlot(seu, reduction = "tsne", group.by = "cluster_majority_type",
      cols = maj_pal, pt.size = PLT$pt_size, label = TRUE, repel = TRUE,
      label.size = PLT$label_size) +
      ggtitle(paste(label, " -  t-SNE (majority cluster type)")) +
      theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
    panels <- c(panels, list(p_maj))
    save_plot(p_maj, file.path(outdir, paste0(label, "_tsne_majority_clusters")), PLT$w_single, PLT$h_single)
  }
  if ("condition" %in% colnames(seu@meta.data) && length(unique(seu$condition)) > 1L)
    panels <- c(panels, list(DimPlot(seu, reduction = "tsne", group.by = "condition",
      cols = COND_PALETTE, pt.size = PLT$pt_size) + ggtitle(paste(label, " -  t-SNE (condition)"))))
  save_plot(wrap_plots(panels, ncol = min(length(panels), 3L)),
            file.path(outdir, paste0(label, "_tsne_grid")), PLT$w_triple, PLT$h_single)
  pt_feats <- intersect(c("latent_time","velocity_pseudotime"), colnames(seu@meta.data))
  if (length(pt_feats) > 0)
    save_plot(FeaturePlot(seu, reduction = "tsne", features = pt_feats,
      cols = AAAS_PSEUDO_COLS, pt.size = PLT$pt_size, ncol = length(pt_feats)),
      file.path(outdir, paste0(label, "_tsne_pseudotime")),
      PLT$w_single * length(pt_feats), PLT$h_single)
}

plot_majority_dimplots <- function(seu, label, outdir) {
  if (!"cluster_majority_type" %in% colnames(seu@meta.data)) return(invisible(NULL))
  message("  Majority-rule UMAP DimPlots: ", label)
  maj_pal <- .build_ct_palette(unique(seu$cluster_majority_type))
  shared_theme <- theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  if ("umap" %in% names(seu@reductions)) {
    p_umap_maj <- DimPlot(seu, reduction = "umap", group.by = "cluster_majority_type",
      cols = maj_pal, pt.size = PLT$pt_size, label = TRUE, repel = TRUE,
      label.size = PLT$label_size) +
      ggtitle(paste(label, " -  UMAP (majority cluster type)")) + shared_theme
    save_plot(p_umap_maj, file.path(outdir, paste0(label, "_umap_majority_clusters")), PLT$w_single, PLT$h_single)
  }
  if ("umap" %in% names(seu@reductions) && "tsne" %in% names(seu@reductions)) {
    p_u <- DimPlot(seu, reduction = "umap", group.by = "cluster_majority_type",
      cols = maj_pal, pt.size = PLT$pt_size, label = TRUE, repel = TRUE,
      label.size = PLT$label_size) + ggtitle("UMAP  -  majority type") + shared_theme + NoLegend()
    p_t <- DimPlot(seu, reduction = "tsne", group.by = "cluster_majority_type",
      cols = maj_pal, pt.size = PLT$pt_size, label = TRUE, repel = TRUE,
      label.size = PLT$label_size) + ggtitle("t-SNE  -  majority type") + shared_theme
    save_plot(p_u | p_t, file.path(outdir, paste0(label, "_umap_tsne_majority_side_by_side")),
              PLT$w_double, PLT$h_single)
  }
}

plot_feature_plots <- function(seu, label, outdir, n_types_max = 8L, n_genes_per_type = 2L,
                               marker_list = MOUSE_MARKERS) {
  if (!"cell_type" %in% colnames(seu@meta.data)) return(invisible(NULL))
  message("  Feature plots: ", label)
  type_counts <- sort(table(seu$cell_type[seu$cell_type != "Unassigned"]), decreasing = TRUE)
  top_types   <- names(type_counts)[seq_len(min(n_types_max, length(type_counts)))]
  feat_genes  <- character(0)
  for (ct in top_types) {
    found <- intersect(marker_list[[ct]], rownames(seu))
    if (length(found) == 0) next
    feat_genes <- unique(c(feat_genes, head(found, n_genes_per_type)))
  }
  if (length(feat_genes) == 0) return(invisible(NULL))
  feat_genes <- head(feat_genes, 24L); ncols <- min(4L, length(feat_genes))
  p_feat <- FeaturePlot(seu, features = feat_genes, cols = AAAS_SPEED_COLS,
    pt.size = PLT$pt_size * 0.7, ncol = ncols, order = TRUE) &
    theme(plot.title = element_text(size = 9), axis.title = element_blank())
  h_calc <- ceiling(length(feat_genes) / ncols) * 3.5
  save_plot(p_feat, file.path(outdir, paste0(label, "_feature_plots")),
            ncols * 4, h_calc, limitsize = FALSE)
}

plot_dot_plot <- function(seu, label, outdir, n_genes_per_type = 2L,
                          marker_list = MOUSE_MARKERS) {
  if (!"cluster_majority_type" %in% colnames(seu@meta.data)) return(invisible(NULL))
  message("  Dot plot: ", label)
  maj_types <- unique(seu$cluster_majority_type); maj_types <- maj_types[maj_types != "Unassigned"]
  dot_genes <- character(0)
  for (ct in maj_types) {
    found <- intersect(marker_list[[ct]], rownames(seu)); if (length(found) == 0) next
    dot_genes <- unique(c(dot_genes, head(found, n_genes_per_type)))
  }
  if (length(dot_genes) < 3) return(invisible(NULL))
  dot_genes <- head(dot_genes, 40L)
  Idents(seu) <- "cluster_majority_type"
  n_ct_dot <- length(unique(seu$cluster_majority_type))
  dims_dot <- .auto_dims_v(n_ct_dot, base_w = 12, base_h = max(10, length(dot_genes) * 0.45), per_item_w = 0.8)
  p_dot <- DotPlot(seu, features = dot_genes, cols = AAAS_SPEED_COLS, dot.scale = 5) + coord_flip() +
    labs(title = paste(label, " -  Dot Plot"), x = "Gene", y = "Cell type") +
    .theme_velo(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_dot, file.path(outdir, paste0(label, "_dot_plot")), dims_dot$w, dims_dot$h, limitsize = FALSE)
  Idents(seu) <- "seurat_clusters"
}

plot_violin_plots <- function(seu, label, outdir, n_genes = 9L,
                              marker_list = MOUSE_MARKERS) {
  message("  Violin plots: ", label)
  type_counts <- sort(table(seu$cell_type[seu$cell_type != "Unassigned"]), decreasing = TRUE)
  top_types   <- names(type_counts)[seq_len(min(6L, length(type_counts)))]
  vln_genes   <- character(0)
  for (ct in top_types) {
    found <- intersect(marker_list[[ct]], rownames(seu)); if (length(found) == 0) next
    vln_genes <- unique(c(vln_genes, head(found, 2L)))
  }
  if (length(vln_genes) == 0) return(invisible(NULL))
  vln_genes <- head(vln_genes, n_genes); n_rows <- ceiling(length(vln_genes) / 3L)
  Idents(seu) <- "seurat_clusters"
  p_vln_clust <- VlnPlot(seu, features = vln_genes, pt.size = 0, ncol = 3L) &
    .theme_velo(base_size = 10) & theme(axis.title.x = element_blank())
  dims_vln <- .auto_dims_v(length(unique(seu$seurat_clusters)), base_w = 18, base_h = n_rows * 4)
  save_plot(p_vln_clust, file.path(outdir, paste0(label, "_violin_per_cluster")),
            dims_vln$w, dims_vln$h, limitsize = FALSE)
  if ("condition" %in% colnames(seu@meta.data) && length(unique(seu$condition)) > 1L) {
    Idents(seu) <- "condition"
    p_vln_cond <- VlnPlot(seu, features = vln_genes, cols = COND_PALETTE, pt.size = 0, ncol = 3L) &
      .theme_velo(base_size = 10) & theme(axis.title.x = element_blank())
    save_plot(p_vln_cond, file.path(outdir, paste0(label, "_violin_per_condition")), 18, n_rows * 4, limitsize = FALSE)
  }
  vel_metrics <- intersect(c("velocity_length","velocity_confidence","velocity_magnitude","velocity_entropy"),
                            colnames(seu@meta.data))
  if (length(vel_metrics) > 0 && "condition" %in% colnames(seu@meta.data)) {
    Idents(seu) <- "condition"
    p_vmet <- VlnPlot(seu, features = vel_metrics, cols = COND_PALETTE, pt.size = 0,
                       ncol = min(4L, length(vel_metrics))) &
      .theme_velo(base_size = 11) & theme(axis.title.x = element_blank())
    save_plot(p_vmet, file.path(outdir, paste0(label, "_violin_velocity_metrics")),
              14 * ceiling(length(vel_metrics)/2), 7 * ceiling(length(vel_metrics)/2))
  }
  Idents(seu) <- "seurat_clusters"
}

plot_stacked_composition <- function(seu, label, outdir) {
  if (!"cell_type" %in% colnames(seu@meta.data)) return(invisible(NULL))
  message("  Stacked composition: ", label)
  ct_palette <- .build_ct_palette(sort(unique(seu$cell_type)))
  .stacked_bar <- function(df, x_col, title, pal) {
    df2 <- df %>%
      dplyr::group_by(.data[[x_col]], cell_type) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(.data[[x_col]]) %>%
      dplyr::mutate(pct = n / sum(n) * 100) %>% dplyr::ungroup()
    ggplot(df2, aes(x = .data[[x_col]], y = pct, fill = cell_type)) +
      geom_col(position = "stack", width = 0.8) + scale_fill_manual(values = pal) +
      labs(x = NULL, y = "% cells", title = title, fill = "Cell type") +
      .theme_velo() + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      guides(fill = guide_legend(ncol = 1))
  }
  meta <- seu@meta.data
  if ("condition" %in% colnames(meta) && length(unique(meta$condition)) > 1L) {
    p_cond <- .stacked_bar(meta, "condition", paste(label, " -  Composition by condition"), ct_palette)
    save_plot(p_cond, file.path(outdir, paste0(label, "_stacked_composition_condition")),
              max(8, length(unique(meta$condition)) * 2.5), 9)
  }
  for (grp in c("label","sample")) {
    if (!grp %in% colnames(meta) || length(unique(meta[[grp]])) < 2L) next
    p_grp <- .stacked_bar(meta, grp, paste(label, paste0(" -  Composition by ", grp)), ct_palette)
    save_plot(p_grp, file.path(outdir, paste0(label, "_stacked_composition_", grp)),
              max(10, length(unique(meta[[grp]])) * 2.5), 9)
    break
  }
}

plot_velocity_density <- function(seu, label, outdir) {
  metrics <- intersect(c("velocity_length","velocity_confidence","velocity_magnitude",
                          "velocity_entropy","latent_time","splicing_ratio"),
                        colnames(seu@meta.data))
  if (length(metrics) == 0) return(invisible(NULL))
  message("  Velocity density/CDF plots: ", label)
  group_col <- if ("condition" %in% colnames(seu@meta.data) &&
                    length(unique(seu$condition)) > 1L) "condition" else
               if ("cell_type" %in% colnames(seu@meta.data)) "cell_type" else NULL
  if (is.null(group_col)) return(invisible(NULL))
  pal <- if (group_col == "condition") COND_PALETTE else
           .build_ct_palette(unique(seu@meta.data[[group_col]]))

  meta_df <- seu@meta.data

  density_plots <- lapply(metrics, function(m) {
    df <- meta_df[!is.na(meta_df[[m]]), c(group_col, m)]
    if (nrow(df) < 10) return(NULL)
    ggplot(df, aes(x = .data[[m]], colour = .data[[group_col]], fill = .data[[group_col]])) +
      geom_density(alpha = 0.25) +
      scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
      labs(x = gsub("_"," ", m), y = "Density", title = gsub("_"," ", m),
           colour = NULL, fill = NULL) +
      .theme_velo(base_size = 11)
  })
  density_plots <- Filter(Negate(is.null), density_plots)
  if (length(density_plots) > 0) {
    p_dens <- wrap_plots(density_plots, ncol = min(3, length(density_plots)))
    save_plot(p_dens, file.path(outdir, paste0(label, "_velocity_density_plots")),
              PLT$w_triple, PLT$h_single * ceiling(length(density_plots)/3), limitsize = FALSE)
  }

  if ("latent_time" %in% metrics && "condition" %in% colnames(meta_df)) {
    lt_df <- meta_df[!is.na(meta_df$latent_time), c("condition","latent_time")]
    p_cdf <- ggplot(lt_df, aes(latent_time, colour = condition)) +
      stat_ecdf(linewidth = 0.9) +
      scale_colour_manual(values = COND_PALETTE) +
      labs(x = "Latent time", y = "Cumulative fraction",
           title = paste("Latent time CDF  - ", label), colour = NULL) +
      .theme_velo()
    save_plot(p_cdf, file.path(outdir, paste0(label, "_latent_time_CDF")), 12, 8)
  }
}

plot_trajectory <- function(seu, label, outdir) {
  has_lt  <- "latent_time"         %in% colnames(seu@meta.data)
  has_vpt <- "velocity_pseudotime" %in% colnames(seu@meta.data)
  if (!has_lt && !has_vpt) return(invisible(NULL))
  if (!"umap" %in% names(seu@reductions)) return(invisible(NULL))
  message("  Trajectory / pseudotime: ", label)
  umap_df <- as.data.frame(Embeddings(seu, "umap")); colnames(umap_df) <- c("UMAP_1","UMAP_2")
  umap_df$cell_type <- seu$cell_type
  ct_palette <- .build_ct_palette(unique(umap_df$cell_type))
  if (has_lt) {
    umap_df$latent_time <- seu$latent_time
    lt_ordered <- umap_df[order(umap_df$latent_time, na.last = TRUE), ]
    p_lt_traj <- ggplot(lt_ordered, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(colour = latent_time), size = PLT$pt_size * 0.9, alpha = 0.85) +
      scale_colour_gradientn(colours = colorRampPalette(AAAS_LATENT_COLS)(100),
                              na.value = "grey80", name = "Latent time") +
      labs(title = paste(label, " -  Trajectory (latent time)")) +
      theme_void(base_size = PLT$base_size) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_plot(p_lt_traj, file.path(outdir, paste0(label, "_trajectory_latent_time")), PLT$w_single, PLT$h_single)
    if (has_lt && requireNamespace("ggridges", quietly = TRUE)) {
      pd <- umap_df[!is.na(umap_df$latent_time) & umap_df$cell_type != "Unassigned", ]
      if (nrow(pd) > 100) {
        p_ridge <- ggplot(pd, aes(x = latent_time, y = cell_type, fill = cell_type)) +
          ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9, colour = "white") +
          scale_fill_manual(values = ct_palette) +
          labs(x = "Latent time", y = NULL,
               title = paste(label, " -  Pseudotime distribution by cell type")) +
          .theme_velo(base_size = PLT$base_size * 1.8) + theme(legend.position = "none")
        h_calc <- max(7, length(unique(pd$cell_type)) * 0.5) * 2.0
        save_plot(p_ridge, file.path(outdir, paste0(label, "_pseudotime_ridges")),
                  18, h_calc, limitsize = FALSE)
      }
    }
  }
  if (has_vpt) {
    umap_df$vpt <- seu$velocity_pseudotime
    vpt_ordered <- umap_df[order(umap_df$vpt, na.last = TRUE), ]
    p_vpt_traj <- ggplot(vpt_ordered, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(colour = vpt), size = PLT$pt_size * 0.9, alpha = 0.85) +
      scale_colour_gradientn(colours = colorRampPalette(AAAS_PSEUDO_COLS)(100),
                              na.value = "grey80", name = "Pseudotime") +
      labs(title = paste(label, " -  Trajectory (velocity pseudotime)")) +
      theme_void(base_size = PLT$base_size) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_plot(p_vpt_traj, file.path(outdir, paste0(label, "_trajectory_velocity_pseudotime")),
              PLT$w_single, PLT$h_single)
  }
}

plot_expression_heatmap <- function(seu, label, outdir,
                                      n_genes_per_cluster = 5L,
                                      max_cells_per_cluster = 100L) {
  message("  Expression heatmap: ", label)
  Idents(seu) <- "seurat_clusters"
  markers <- tryCatch(
    FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3,
                   test.use = "wilcox", verbose = FALSE),
    error = function(e) { message("    FindAllMarkers failed: ", e$message); NULL })
  if (is.null(markers) || nrow(markers) == 0) return(invisible(NULL))
  top_genes <- markers %>% dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = n_genes_per_cluster) %>%
    pull(gene) %>% unique()
  if (length(top_genes) < 3) return(invisible(NULL))
  cells_use <- Seurat::WhichCells(seu, downsample = max_cells_per_cluster)
  base_label_size <- max(5, 8 - length(top_genes) %/% 10)
  p_doheat <- DoHeatmap(seu, features = top_genes, cells = cells_use, size = 6,
    angle = 90, hjust = 0.5, draw.lines = TRUE) +
    ggtitle(paste(label, " -  Top marker gene heatmap")) +
    theme(axis.text.y = element_text(size = base_label_size * 2, face = "bold"),
          axis.text.x = element_text(size = 8),
          plot.title  = element_text(size = 16, hjust = 0.5, face = "bold"))
  save_plot(p_doheat, file.path(outdir, paste0(label, "_heatmap_top_markers")),
            PLT$w_single * 1.5, max(12, length(top_genes) * 0.32) * 2.0, limitsize = FALSE)
  expr_mat <- AverageExpression(seu, features = top_genes, group.by = "seurat_clusters",
                                 assays = "RNA", slot = "data")$RNA
  if (!is.null(expr_mat) && nrow(expr_mat) > 1) {
    expr_z <- t(scale(t(as.matrix(expr_mat)))); expr_z[is.nan(expr_z)] <- 0
    lim  <- min(3, max(abs(expr_z), na.rm = TRUE))
    h_px <- max(600, nrow(expr_z) * 14) * 5L; w_px <- max(600, ncol(expr_z) * 55) * 3L
    ph <- pheatmap(expr_z,
                   color = colorRampPalette(c("#2166AC","#F7F7F7","#D6604D"))(100),
                   breaks = seq(-lim, lim, length.out = 101),
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   fontsize_row = max(6, 10 - nrow(expr_z) %/% 10), fontsize_col = 9,
                   main = paste("Avg. expression (z-score)\n", label),
                   border_color = NA, silent = TRUE)
    save_pheatmap(ph, file.path(outdir, paste0(label, "_heatmap_avg_expr")), w_px, h_px)
  }
  write.csv(markers, file.path(outdir, paste0(label, "_cluster_markers.csv")), row.names = FALSE)
  message("    Heatmaps + markers CSV saved.")
  return(invisible(markers))
}

run_per_replicate <- function() {
  message("\n", strrep("=", 70))
  message("PER-REPLICATE NEURAL CREST RNA VELOCITY")
  message(strrep("=", 70))
  outdir      <- file.path(OUTDIR, "per_replicate")
  rep_objects <- list()

  for (sname in names(SAMPLES)) {
    lbl      <- SAMPLES[[sname]]$label
    rds_path <- file.path(outdir, paste0(lbl, "_velocity.rds"))
    if (!FORCE_RERUN && file.exists(rds_path)) {
      message("  [CHECKPOINT] ", lbl, " -- loading cached RDS")
      rep_objects[[sname]] <- readRDS(rds_path)
      next
    }
    message("\n--- Sample: ", lbl, " ---")

    seu <- load_replicate_loom(sname)
    seu <- preprocess_seurat(seu, lbl, use_harmony = FALSE)

    # Level-1: classify all cells using the whole-embryo marker panel
    seu <- annotate_by_markers(seu, MOUSE_MARKERS)

    # Subset to Neural_crest cells and apply Level-2 NC sub-state annotation
    seu <- .subset_to_nc_and_reannotate(seu, lbl, use_harmony = FALSE)
    if (ncol(seu) < NC_MIN_CELLS) {
      message("  Skipping ", lbl, ": insufficient Neural_crest cells after subsetting.")
      next
    }

    vr_res <- run_velociraptor(seu, lbl)
    seu    <- vr_res$seu
    seu    <- compute_per_cell_velocity_metrics(seu)

    plot_overview(seu, lbl, outdir)
    plot_velocity_arrows(seu, vr_res$velo_out, vr_res$sce_full, lbl, outdir)
    extract_transition_matrix(seu, vr_res$velo_out, lbl,
                              file.path(OUTDIR, "celltype_transitions"))
    plot_ct_composition(seu, lbl, outdir)
    plot_pca_panels(seu, lbl, outdir)
    plot_tsne_panels(seu, lbl, outdir)
    plot_majority_dimplots(seu, lbl, outdir)
    plot_feature_plots(seu, lbl, outdir, marker_list = Neural_crest_L2)
    plot_dot_plot(seu, lbl, outdir, marker_list = Neural_crest_L2)
    plot_violin_plots(seu, lbl, outdir, marker_list = Neural_crest_L2)
    plot_stacked_composition(seu, lbl, outdir)
    plot_trajectory(seu, lbl, outdir)
    plot_velocity_density(seu, lbl, outdir)
    plot_expression_heatmap(seu, lbl, outdir)
    run_velocity_residuals(seu, lbl, outdir)
    compute_splicing_efficiency_index(seu, lbl,
                                      file.path(OUTDIR, "splicing_efficiency"))
    compute_velocity_pseudotime_correlation(seu, lbl,
                                            file.path(OUTDIR, "velocity_pseudotime_corr"))

    saveRDS(seu, rds_path)
    rep_objects[[sname]] <- seu
    message("  Done: ", lbl)
  }

  rep_objects
}

run_per_condition <- function() {
  message("\n", strrep("=", 70))
  message("PER-CONDITION NEURAL CREST RNA VELOCITY")
  message(strrep("=", 70))
  outdir       <- file.path(OUTDIR, "per_condition")
  cond_objects <- list()

  for (cond in CONDITIONS) {
    rds_path <- file.path(outdir, paste0(cond, "_velocity.rds"))
    if (!FORCE_RERUN && file.exists(rds_path)) {
      message("  [CHECKPOINT] ", cond, " -- loading cached RDS")
      cond_objects[[cond]] <- readRDS(rds_path)
      next
    }
    message("\n--- Condition: ", cond, " ---")
    seu <- load_merged_loom(LOOM_FILES$per_condition[[cond]], project_name = cond)
    seu <- preprocess_seurat(seu, cond, use_harmony = TRUE, harmony_vars = "replicate")

    # Level-1: classify all cells using the whole-embryo marker panel
    seu <- annotate_by_markers(seu, MOUSE_MARKERS)

    # Subset to Neural_crest cells and apply Level-2 NC sub-state annotation
    seu <- .subset_to_nc_and_reannotate(seu, cond,
                                         use_harmony  = TRUE,
                                         harmony_vars = "replicate")
    if (ncol(seu) < NC_MIN_CELLS) {
      message("  Skipping ", cond, ": insufficient Neural_crest cells after subsetting.")
      next
    }

    vr_res <- run_velociraptor(seu, cond)
    seu    <- vr_res$seu
    seu    <- compute_per_cell_velocity_metrics(seu)

    plot_overview(seu, cond, outdir, extra_groups = c("replicate","condition"))
    plot_velocity_arrows(seu, vr_res$velo_out, vr_res$sce_full, cond, outdir)
    extract_transition_matrix(seu, vr_res$velo_out, cond,
                              file.path(OUTDIR, "celltype_transitions"))
    plot_ct_composition(seu, cond, outdir)
    plot_pca_panels(seu, cond, outdir)
    plot_tsne_panels(seu, cond, outdir)
    plot_majority_dimplots(seu, cond, outdir)
    plot_feature_plots(seu, cond, outdir, marker_list = Neural_crest_L2)
    plot_dot_plot(seu, cond, outdir, marker_list = Neural_crest_L2)
    plot_violin_plots(seu, cond, outdir, marker_list = Neural_crest_L2)
    plot_stacked_composition(seu, cond, outdir)
    plot_trajectory(seu, cond, outdir)
    plot_velocity_density(seu, cond, outdir)
    plot_expression_heatmap(seu, cond, outdir)
    run_velocity_residuals(seu, cond, outdir)
    compute_splicing_efficiency_index(seu, cond,
                                      file.path(OUTDIR, "splicing_efficiency"))
    compute_velocity_pseudotime_correlation(seu, cond,
                                            file.path(OUTDIR, "velocity_pseudotime_corr"))

    saveRDS(seu, rds_path)
    cond_objects[[cond]] <- seu
    message("  Done: ", cond)
  }
  cond_objects
}

.fast_mwu_velocity <- function(vel_q, vel_r, common_genes, query_label, ref_label) {
  vq_sub <- vel_q[common_genes, , drop = FALSE]
  vr_sub <- vel_r[common_genes, , drop = FALSE]
  nz_q <- Matrix::rowSums(vq_sub != 0, na.rm = TRUE)
  nz_r <- Matrix::rowSums(vr_sub != 0, na.rm = TRUE)
  testable <- common_genes[nz_q >= 5L & nz_r >= 5L]
  if (length(testable) == 0) return(data.frame())
  message("    Testable velocity genes (>=5 non-zero): ", length(testable), " / ", length(common_genes))
  mat_q <- as.matrix(vq_sub[testable, , drop = FALSE])
  mat_r <- as.matrix(vr_sub[testable, , drop = FALSE])
  if (requireNamespace("matrixTests", quietly = TRUE)) {
    mat_q[is.nan(mat_q)] <- NA; mat_r[is.nan(mat_r)] <- NA
    res <- tryCatch(matrixTests::row_wilcoxon_twosample(mat_q, mat_r, alternative = "two.sided"),
                    error = function(e) NULL)
    if (!is.null(res)) {
      rows <- data.frame(
        gene           = testable,
        delta_velocity = rowMeans(mat_q, na.rm = TRUE) - rowMeans(mat_r, na.rm = TRUE),
        mean_vel_query = rowMeans(mat_q, na.rm = TRUE),
        mean_vel_ref   = rowMeans(mat_r, na.rm = TRUE),
        pct_pos_query  = rowMeans(mat_q > 0, na.rm = TRUE) * 100,
        pct_pos_ref    = rowMeans(mat_r > 0, na.rm = TRUE) * 100,
        n_cells_query  = rowSums(!is.na(mat_q)),
        n_cells_ref    = rowSums(!is.na(mat_r)),
        W_stat         = res$statistic,
        pvalue         = res$pvalue,
        stringsAsFactors = FALSE)
      rows <- rows[!is.na(rows$pvalue), ]
      rows$padj <- p.adjust(rows$pvalue, method = "BH")
      return(rows[order(rows$padj), ])
    }
  }
  rows <- lapply(seq_len(length(testable)), function(i) {
    vq <- mat_q[i, ]; vq <- vq[!is.na(vq)]
    vr <- mat_r[i, ]; vr <- vr[!is.na(vr)]
    if (length(vq) < 5 || length(vr) < 5) return(NULL)
    test <- wilcox.test(vq, vr, exact = FALSE)
    data.frame(gene = testable[i],
               delta_velocity = mean(vq) - mean(vr),
               mean_vel_query = mean(vq), mean_vel_ref = mean(vr),
               pct_pos_query = mean(vq > 0) * 100, pct_pos_ref = mean(vr > 0) * 100,
               n_cells_query = length(vq), n_cells_ref = length(vr),
               W_stat = unname(test$statistic), pvalue = test$p.value, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$padj <- p.adjust(df$pvalue, method = "BH"); df[order(df$padj), ]
}

.compute_speed_from_matrix <- function(vel_mat_sub) colMeans(abs(vel_mat_sub), na.rm = TRUE)
.compute_coherence_from_matrix <- function(vel_mat_sub) {
  cell_mean <- colMeans(vel_mat_sub, na.rm = TRUE)
  agree     <- sweep(sign(vel_mat_sub), 2, sign(cell_mean), `==`)
  colMeans(agree, na.rm = TRUE)
}

.celltype_velocity_summary <- function(seu, label, outdir) {
  if (!"cell_type" %in% colnames(seu@meta.data) || !"velocity" %in% names(seu@assays))
    return(data.frame())
  vel_mat <- .get_layer(seu, "velocity", "counts")
  ct_vals <- seu$cell_type; ct_list <- setdiff(unique(ct_vals), "Unassigned")
  if (length(ct_list) == 0) return(data.frame())
  has_speed <- "velocity_length"     %in% colnames(seu@meta.data) && !all(is.na(seu$velocity_length))
  has_coher <- "velocity_confidence" %in% colnames(seu@meta.data) && !all(is.na(seu$velocity_confidence))
  has_lt    <- "latent_time"         %in% colnames(seu@meta.data)
  has_mag   <- "velocity_magnitude"  %in% colnames(seu@meta.data)
  rows <- lapply(ct_list, function(ct) {
    mask <- ct_vals == ct; n <- sum(mask); if (n < 5) return(NULL)
    vm_sp <- vel_mat[, mask, drop = FALSE]
    spd <- if (has_speed) mean(seu$velocity_length[mask], na.rm = TRUE) else
              mean(.compute_speed_from_matrix(as.matrix(vm_sp)), na.rm = TRUE)
    coh <- if (has_coher) mean(seu$velocity_confidence[mask], na.rm = TRUE) else
              mean(.compute_coherence_from_matrix(as.matrix(vm_sp)), na.rm = TRUE)
    mag <- if (has_mag) mean(seu$velocity_magnitude[mask], na.rm = TRUE) else NA_real_
    cell_mean_vel <- as.numeric(Matrix::colMeans(abs(vm_sp), na.rm = TRUE))
    pct_dynamic   <- mean(cell_mean_vel > quantile(cell_mean_vel, 0.25, na.rm = TRUE), na.rm = TRUE) * 100
    list(cell_type = ct, n_cells = n, mean_speed = spd, mean_coherence = coh,
         mean_magnitude = mag,
         mean_latent_time = if (has_lt) mean(seu$latent_time[mask], na.rm = TRUE) else NA_real_,
         pct_pos_velocity = mean(rowMeans(vm_sp, na.rm = TRUE) > 0, na.rm = TRUE) * 100,
         pct_dynamic_cells = pct_dynamic,
         overall_mean_vel = mean(vm_sp, na.rm = TRUE))
  })
  df <- do.call(rbind, lapply(Filter(Negate(is.null), rows), as.data.frame))
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  df$mean_speed[is.nan(df$mean_speed)] <- NA_real_
  df$mean_coherence[is.nan(df$mean_coherence)] <- NA_real_
  if (nrow(df) >= 2) {
    df_plot <- df[!is.na(df$mean_speed) | !is.na(df$mean_coherence), ]
    if (nrow(df_plot) >= 2) {
      p_metrics <- df_plot %>%
        tidyr::pivot_longer(cols = c(mean_speed, mean_coherence, pct_dynamic_cells),
                            names_to = "metric", values_to = "value") %>%
        dplyr::filter(!is.na(value)) %>%
        ggplot(aes(y = reorder(cell_type, value), x = value, fill = metric)) +
        geom_col(position = "dodge", alpha = 0.8) +
        scale_fill_manual(values = c(mean_speed = AAAS_COLS[2],
                                      mean_coherence = AAAS_COLS[1],
                                      pct_dynamic_cells = AAAS_COLS[3])) +
        labs(x = "Mean value per cell type", y = NULL,
             title = paste(label, " -  Velocity metrics by cell type"), fill = NULL) +
        .theme_velo(base_size = 12) + theme(legend.position = "bottom") +
        facet_wrap(~metric, scales = "free_x")
      save_plot(p_metrics, file.path(outdir, paste0(label, "_celltype_velocity_metrics")),
                PLT$w_double, max(6, nrow(df_plot) * 0.45), limitsize = FALSE)
    }
  }
  df
}

.plot_diff_results <- function(df, sig, query_label, ref_label, label, outdir) {
  col_up <- .pal_get(COND_PALETTE, query_label) %||% AAAS_COLS[2]
  col_dn <- .pal_get(COND_PALETTE, ref_label)   %||% AAAS_COLS[1]
  df$category <- "n.s."
  df$category[!is.na(df$padj) & df$padj < 0.05 & df$delta_velocity > 0] <- paste0("Up in ", query_label)
  df$category[!is.na(df$padj) & df$padj < 0.05 & df$delta_velocity < 0] <- paste0("Up in ", ref_label)
  df$neg_log10_p <- -log10(df$pvalue + 1e-300)
  df$mean_vel    <- (df$mean_vel_query + df$mean_vel_ref) / 2
  sig       <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$delta_velocity) > 0.01, ]
  top_labs  <- rbind(head(sig[sig$delta_velocity > 0, ], 8), head(sig[sig$delta_velocity < 0, ], 7))
  n_up      <- sum(sig$delta_velocity > 0); n_dn <- sum(sig$delta_velocity < 0)
  cat_cols  <- setNames(c("grey80", col_up, col_dn),
                         c("n.s.", paste0("Up in ", query_label), paste0("Up in ", ref_label)))
  p_vol <- ggplot(df, aes(delta_velocity, neg_log10_p, colour = category)) +
    geom_point(alpha = 0.5, size = 1.2) + geom_point(data = sig, size = 2.5) +
    geom_label_repel(data = top_labs, aes(label = gene), size = 3,
                     max.overlaps = 25, box.padding = 0.35, show.legend = FALSE) +
    scale_colour_manual(values = cat_cols) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
    annotate("text", x = Inf,  y = Inf, label = paste0("up  ", n_up),
             hjust = 1.15, vjust = 1.8, size = 5, fontface = "bold", colour = col_up) +
    annotate("text", x = -Inf, y = Inf, label = paste0(n_dn, " dn "),
             hjust = -0.15, vjust = 1.8, size = 5, fontface = "bold", colour = col_dn) +
    labs(x = bquote(Delta * "Velocity (" * .(query_label) * " - " * .(ref_label) * ")"),
         y = expression(-log[10] * "(p-value)"),
         title = paste("Differential RNA Velocity:", label),
         subtitle = sprintf("Up in %s: %d  |  Up in %s: %d", query_label, n_up, ref_label, n_dn),
         colour = "") + .theme_velo(base_size = 13) + theme(legend.position = "bottom")
  save_plot(p_vol, file.path(outdir, paste0(label, "_volcano")), PLT$w_single, PLT$h_single)
  p_ma <- ggplot(df, aes(mean_vel, delta_velocity, colour = category)) +
    geom_point(alpha = 0.45, size = 1.2) + geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = cat_cols) +
    labs(x = "Mean velocity",
         y = bquote(Delta * "Velocity (" * .(query_label) * " - " * .(ref_label) * ")"),
         title = paste("MA Plot  - ", label), colour = "") + .theme_velo(base_size = 13)
  save_plot(p_ma, file.path(outdir, paste0(label, "_MAplot")), PLT$w_single, PLT$h_single)
  top50 <- head(sig, 50)
  if (nrow(top50) >= 5) {
    heat_mat <- as.matrix(top50[, c("mean_vel_query","mean_vel_ref")])
    rownames(heat_mat) <- top50$gene; colnames(heat_mat) <- c(query_label, ref_label)
    lim <- max(abs(heat_mat), na.rm = TRUE)
    ph  <- pheatmap(heat_mat,
                    color = colorRampPalette(c(col_dn,"#F7F7F7",col_up))(100),
                    breaks = seq(-lim, lim, length.out = 101),
                    cluster_cols = FALSE, fontsize_row = 9,
                    main = paste("Top differential velocity genes\n", label), silent = TRUE)
    save_pheatmap(ph, file.path(outdir, paste0(label, "_heatmap")), 2700, max(3000, nrow(heat_mat)*110))
    top50$marker_for_celltypes <- vapply(top50$gene, function(g) {
      cts <- which(vapply(MOUSE_MARKERS, function(mk) g %in% mk, logical(1)))
      if (length(cts) == 0) " - " else paste(names(cts), collapse = ", ")
    }, character(1))
    write.csv(top50, file.path(outdir, paste0(label, "_top_diff_genes_annotated.csv")), row.names = FALSE)
  }
}

.jaccard_velocity_genes <- function(seu_q, seu_r, query_label, ref_label, outdir) {
  label <- paste0(query_label, "_vs_", ref_label)
  get_vgenes <- function(seu) {
    if (!"velocity" %in% names(seu@assays)) return(character(0))
    rownames(.get_layer(seu, "velocity", "counts"))
  }
  vg_q <- get_vgenes(seu_q); vg_r <- get_vgenes(seu_r)
  jacc <- length(intersect(vg_q, vg_r)) / max(length(union(vg_q, vg_r)), 1)
  overlap_df <- data.frame(
    metric = c("n_velocity_genes_query","n_velocity_genes_ref",
               "n_shared","n_unique_query","n_unique_ref","jaccard"),
    value  = c(length(vg_q), length(vg_r),
               length(intersect(vg_q, vg_r)), length(setdiff(vg_q, vg_r)),
               length(setdiff(vg_r, vg_q)), jacc))
  write.csv(overlap_df, file.path(outdir, paste0(label, "_velocity_gene_overlap.csv")), row.names = FALSE)
  message("  Jaccard overlap: ", round(jacc, 3)); overlap_df
}

.diff_velocity_by_celltype <- function(seu_q, seu_r, query_label, ref_label, outdir) {
  label     <- paste0(query_label, "_vs_", ref_label)
  ct_outdir <- file.path(outdir, "by_celltype")
  dir.create(ct_outdir, recursive = TRUE, showWarnings = FALSE)
  message("  Per-cell-type differential: ", label)
  if (!"velocity" %in% names(seu_q@assays) || !"velocity" %in% names(seu_r@assays)) {
    message("    velocity assay missing  -  skipping per-CT analysis"); return(list())
  }
  vel_q  <- .get_layer(seu_q, "velocity", "counts")
  vel_r  <- .get_layer(seu_r, "velocity", "counts")
  common_genes <- intersect(rownames(vel_q), rownames(vel_r))
  cts_q <- setdiff(unique(seu_q$cell_type), "Unassigned")
  cts_r <- setdiff(unique(seu_r$cell_type), "Unassigned")
  cell_types <- sort(intersect(cts_q, cts_r))
  message("    Cell types in both conditions: ", length(cell_types))
  ct_results <- list(); ct_summary <- list(); ct_top_genes <- list()
  for (ct in cell_types) {
    mask_q <- seu_q$cell_type == ct; mask_r <- seu_r$cell_type == ct
    n_q <- sum(mask_q); n_r <- sum(mask_r)
    if (n_q < 10 || n_r < 10) next
    df_ct <- .fast_mwu_velocity(vel_q[, mask_q, drop = FALSE],
                                  vel_r[, mask_r, drop = FALSE],
                                  common_genes, query_label, ref_label)
    if (nrow(df_ct) == 0) next
    ct_safe <- gsub("[/ ]", "_", ct)
    write.csv(df_ct, file.path(ct_outdir, paste0(label, "_", ct_safe, "_differential.csv")), row.names = FALSE)
    col_up <- .pal_get(COND_PALETTE, query_label) %||% AAAS_COLS[2]
    col_dn <- .pal_get(COND_PALETTE, ref_label)   %||% AAAS_COLS[1]
    df_ct$neg_log10_p <- -log10(df_ct$pvalue + 1e-300)
    df_ct$cat <- ifelse(!is.na(df_ct$padj) & df_ct$padj < 0.05 & df_ct$delta_velocity > 0, paste0("Up ", query_label),
                  ifelse(!is.na(df_ct$padj) & df_ct$padj < 0.05 & df_ct$delta_velocity < 0, paste0("Up ", ref_label), "n.s."))
    sig_ct <- df_ct[!is.na(df_ct$padj) & df_ct$padj < 0.05 & abs(df_ct$delta_velocity) > 0.01, ]
    n_sig  <- nrow(sig_ct)
    if (n_sig >= 3) {
      n_arm   <- min(15L, ceiling(n_sig / 2))
      top_up  <- head(sig_ct[order(-sig_ct$delta_velocity), ], n_arm)
      top_dn  <- head(sig_ct[order( sig_ct$delta_velocity), ], n_arm)
      ll_df   <- unique(rbind(top_up, top_dn))
      ll_df   <- ll_df[order(ll_df$delta_velocity), ]
      ll_df$gene      <- factor(ll_df$gene, levels = ll_df$gene)
      ll_df$dot_col   <- ifelse(ll_df$delta_velocity > 0, col_up, col_dn)
      ll_df$neglogfdr <- -log10(ll_df$padj + 1e-300)
      p_lolly <- ggplot(ll_df, aes(x = delta_velocity, y = gene)) +
        geom_segment(aes(x = 0, xend = delta_velocity, yend = gene),
                     colour = "grey75", linewidth = 0.7) +
        geom_point(aes(size = neglogfdr, colour = dot_col, fill = dot_col),
                   shape = 21, stroke = 0.4, alpha = 0.88) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   colour = "grey35", linewidth = 0.5) +
        scale_colour_identity(aesthetics = c("colour","fill")) +
        scale_size_continuous(range = c(2.5, 7),
                              name  = expression(-log[10](FDR))) +
        labs(x        = paste0("Delta Velocity  (", query_label, " - ", ref_label, ")"),
             y        = NULL,
             title    = paste0("Differential velocity  -  ", gsub("_"," ", ct)),
             subtitle = paste0(query_label, " vs ", ref_label,
                               "   up ", sum(sig_ct$delta_velocity > 0), " up",
                               "   dn ", sum(sig_ct$delta_velocity < 0), " down",
                               "   FDR < 0.05")) +
        .theme_velo(base_size = 11) +
        theme(legend.position      = "right",
              legend.title         = element_text(size = 9),
              legend.text          = element_text(size = 8),
              axis.text.y          = element_text(size = 9, face = "italic"),
              axis.text.x          = element_text(size = 10),
              axis.title.x         = element_text(size = 10),
              plot.title           = element_text(size = 11, face = "bold"),
              plot.subtitle        = element_text(size = 9, colour = "grey40"),
              panel.grid.major.y   = element_line(colour = "grey93", linewidth = 0.35))
      save_plot(p_lolly,
                file.path(ct_outdir, paste0(label, "_", ct_safe, "_lollipop")),
                PLT$w_single, max(5, nrow(ll_df) * 0.30 + 2.5),
                limitsize = FALSE)
    }
    message("    ", ct, ": n_q=", n_q, ", n_r=", n_r, ", sig_genes=", n_sig)
    ct_results[[ct]] <- df_ct
    ct_summary[[ct]] <- data.frame(cell_type = ct, n_cells_query = n_q, n_cells_ref = n_r,
      n_sig_genes = n_sig, mean_delta_vel_sig = if (n_sig > 0) mean(sig_ct$delta_velocity, na.rm = TRUE) else 0,
      frac_sig_up_query = if (n_sig > 0) mean(sig_ct$delta_velocity > 0, na.rm = TRUE) else 0,
      mean_vel_query = mean(df_ct$mean_vel_query, na.rm = TRUE),
      mean_vel_ref   = mean(df_ct$mean_vel_ref,   na.rm = TRUE))
    ct_top_genes[[ct]] <- c(head(sig_ct$gene[order(-sig_ct$delta_velocity)], 5),
                              head(sig_ct$gene[order( sig_ct$delta_velocity)], 5))
  }
  if (length(ct_results) == 0) return(list())
  summary_df <- do.call(rbind, ct_summary)
  summary_df <- summary_df[order(-summary_df$n_sig_genes), ]
  write.csv(summary_df, file.path(outdir, paste0(label, "_by_celltype_summary.csv")), row.names = FALSE)

  n_ct   <- length(ct_results); n_cols <- min(4L, n_ct); n_rows <- ceiling(n_ct / n_cols)
  panel_plots <- lapply(names(ct_results), function(ct) {
    df_ct  <- ct_results[[ct]]
    col_up <- .pal_get(COND_PALETTE, query_label) %||% AAAS_COLS[2]
    col_dn <- .pal_get(COND_PALETTE, ref_label)   %||% AAAS_COLS[1]
    sig_ct <- df_ct[!is.na(df_ct$padj) & df_ct$padj < 0.05 & abs(df_ct$delta_velocity) > 0.01, ]
    n_q    <- ct_summary[[ct]]$n_cells_query; n_r <- ct_summary[[ct]]$n_cells_ref
    ggplot(df_ct, aes(delta_velocity, neg_log10_p, colour = cat)) +
      geom_point(size = 0.8, alpha = 0.5) + geom_point(data = sig_ct, size = 1.5) +
      scale_colour_manual(values = setNames(c("grey80", col_up, col_dn),
        c("n.s.", paste0("Up ", query_label), paste0("Up ", ref_label)))) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "red", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      labs(title = paste0(gsub("_"," ",ct), "\nn=", n_q, "/", n_r,
                           "  up ", sum(sig_ct$delta_velocity > 0),
                           "  dn ", sum(sig_ct$delta_velocity < 0)),
           x = "Delta Velocity", y = "-log10(p)") +
      .theme_velo(base_size = 9) + theme(legend.position = "none")
  })
  save_plot(wrap_plots(panel_plots, ncol = n_cols),
            file.path(outdir, paste0(label, "_by_celltype_multipanel_volcano")),
            max(6, PLT$w_single/2) * n_cols, max(5, PLT$h_single/2) * n_rows, limitsize = FALSE)

  all_top <- unique(unlist(ct_top_genes))[seq_len(min(60, length(unique(unlist(ct_top_genes)))))]
  if (length(all_top) >= 3 && length(ct_results) >= 2) {
    heat_mat <- do.call(rbind, lapply(names(ct_results), function(ct) {
      df_idx <- ct_results[[ct]]
      row    <- setNames(rep(NA_real_, length(all_top)), all_top)
      for (g in all_top) { idx <- match(g, df_idx$gene); if (!is.na(idx)) row[g] <- df_idx$delta_velocity[idx] }
      row
    }))
    rownames(heat_mat) <- names(ct_results)
    lim <- max(max(abs(heat_mat), na.rm = TRUE) * 0.95, 0.001)
    h_px_ht <- max(1400, nrow(heat_mat) * 70)
    w_px_ht <- max(2600, ncol(heat_mat) * 28)
    fs_row  <- max(5, min(9,  floor(h_px_ht / nrow(heat_mat) / PLT$dpi * 72 * 0.55)))
    fs_col  <- max(4, min(7,  floor(w_px_ht / ncol(heat_mat) / PLT$dpi * 72 * 0.45)))
    ph_ht <- pheatmap(heat_mat,
                      color = colorRampPalette(c(.pal_get(COND_PALETTE, ref_label) %||% AAAS_COLS[1],
                                                  "#F7F7F7", .pal_get(COND_PALETTE, query_label) %||% AAAS_COLS[2]))(100),
                      breaks = seq(-lim, lim, length.out = 101), na_col = "grey90",
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      fontsize_row = fs_row, fontsize_col = fs_col,
                      main = paste0("Per-CT Delta Velocity  -  ", label), silent = TRUE)
    save_pheatmap(ph_ht, file.path(outdir, paste0(label, "_by_celltype_heatmap")),
                  w_px_ht, h_px_ht)
  }
  ct_results
}

.compare_velocity_metrics <- function(cond_objects, outdir) {
  message("  Comparing velocity metrics across conditions ...")
  metric_df <- do.call(rbind, lapply(names(cond_objects), function(cond) {
    seu <- cond_objects[[cond]]
    data.frame(
      condition           = cond,
      velocity_length     = if ("velocity_length"     %in% colnames(seu@meta.data)) seu$velocity_length     else NA,
      velocity_confidence = if ("velocity_confidence" %in% colnames(seu@meta.data)) seu$velocity_confidence else NA,
      velocity_magnitude  = if ("velocity_magnitude"  %in% colnames(seu@meta.data)) seu$velocity_magnitude  else NA,
      stringsAsFactors    = FALSE)
  }))
  ks_rows <- lapply(COMPARISONS, function(comp) {
    q <- metric_df$velocity_length[metric_df$condition == comp[1]]
    r <- metric_df$velocity_length[metric_df$condition == comp[2]]
    q <- q[!is.na(q)]; r <- r[!is.na(r)]
    if (length(q) < 2 || length(r) < 2) return(NULL)
    t <- ks.test(q, r)
    data.frame(comparison = paste(comp, collapse = "_vs_"), D_stat = t$statistic, pvalue = t$p.value)
  })
  ks_df <- do.call(rbind, Filter(Negate(is.null), ks_rows))
  if (!is.null(ks_df)) write.csv(ks_df, file.path(outdir, "velocity_speed_KS_test.csv"), row.names = FALSE)
  for (metric in c("velocity_length","velocity_confidence","velocity_magnitude")) {
    md <- metric_df[!is.na(metric_df[[metric]]), ]
    if (nrow(md) < 10) next
    p <- ggplot(md, aes(x = condition, y = .data[[metric]], fill = condition)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_boxplot(width = 0.08, outlier.size = 0.5, fill = "white") +
      scale_fill_manual(values = COND_PALETTE) +
      labs(y = gsub("_"," ", metric), x = NULL, title = paste(gsub("_"," ", metric), "by condition")) +
      theme_classic(base_size = 15) + theme(legend.position = "none")
    save_plot(p, file.path(outdir, paste0(metric, "_violin")), 14, 10)
  }
}

.plot_kinetics_scatter <- function(df_mwu, query_label, ref_label, label, outdir) {
  params <- c("fit_alpha","fit_beta","fit_gamma","fit_likelihood","fit_r2")
  avail  <- intersect(paste0(params, "_query"), colnames(df_mwu))
  if (length(avail) < 2) return(invisible(NULL))
  param_names <- sub("_query$", "", avail)
  plot_list <- lapply(param_names, function(param) {
    cx <- paste0(param, "_query"); cy <- paste0(param, "_ref")
    if (!cx %in% colnames(df_mwu) || !cy %in% colnames(df_mwu)) return(NULL)
    d <- df_mwu[!is.na(df_mwu[[cx]]) & !is.na(df_mwu[[cy]]), ]
    if (nrow(d) < 5) return(NULL)
    is_sig <- !is.na(d$padj) & d$padj < 0.05
    ggplot(d, aes(.data[[cx]], .data[[cy]])) +
      geom_point(aes(colour = is_sig), size = 1.2, alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = AAAS_COLS[2]),
                          labels = c(`FALSE` = "n.s.", `TRUE` = "FDR<0.05"), name = NULL) +
      labs(x = paste(query_label, param), y = paste(ref_label, param),
           title = gsub("fit_", "", param)) +
      .theme_velo(base_size = 11) + theme(legend.position = "bottom")
  })
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0) return(invisible(NULL))
  n_cols <- min(length(plot_list), 3L)
  p_kin  <- wrap_plots(plot_list, ncol = n_cols) +
    plot_annotation(title = paste("Kinetics parameters:", label),
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  save_plot(p_kin, file.path(outdir, paste0(label, "_kinetics_scatter")),
            PLT$w_double + n_cols * 2, max(PLT$h_single * ceiling(length(plot_list)/n_cols) * 0.6, 8),
            limitsize = FALSE)
}

run_differential <- function(cond_objects) {
  message("\n", strrep("=", 70))
  message("DIFFERENTIAL RNA VELOCITY  (MWU + pseudobulk-LM + Jaccard + per-cell-type)")
  message(strrep("=", 70))
  outdir       <- file.path(OUTDIR, "differential")
  diff_results <- list(); diff_by_ct <- list(); all_de_for_integ <- list()

  for (comp in COMPARISONS) {
    ql <- comp[1]; rl <- comp[2]
    lbl <- paste0(ql, "_vs_", rl)
    message("\n--- Comparison: ", lbl, " ---")
    diff_cache_path <- file.path(outdir, paste0(lbl, "_diff_cache.rds"))
    if (!FORCE_RERUN && file.exists(diff_cache_path)) {
      message("  [CHECKPOINT] ", lbl, "  -  loading cached diff results")
      cache <- readRDS(diff_cache_path)
      diff_results[[lbl]] <- cache$df
      diff_by_ct[[lbl]]   <- cache$ct_dfs
      next
    }
    seu_q <- cond_objects[[ql]]; seu_r <- cond_objects[[rl]]
    if (is.null(seu_q) || is.null(seu_r)) {
      message("  Missing condition object  -  skipping ", lbl); next
    }
    if (!"velocity" %in% names(seu_q@assays) || !"velocity" %in% names(seu_r@assays)) {
      message("  velocity assay missing  -  skipping ", lbl); next
    }
    vel_q <- .get_layer(seu_q, "velocity", "counts")
    vel_r <- .get_layer(seu_r, "velocity", "counts")
    common_genes <- intersect(rownames(vel_q), rownames(vel_r))
    message("  Common velocity genes: ", length(common_genes))
    if (length(common_genes) < 3) next

    sv_q <- .celltype_velocity_summary(seu_q, ql, outdir)
    sv_r <- .celltype_velocity_summary(seu_r, rl, outdir)
    if (nrow(sv_q) > 0) write.csv(sv_q, file.path(outdir, paste0(ql, "_celltype_velocity_summary.csv")), row.names = FALSE)
    if (nrow(sv_r) > 0) write.csv(sv_r, file.path(outdir, paste0(rl, "_celltype_velocity_summary.csv")), row.names = FALSE)

    df  <- .fast_mwu_velocity(vel_q, vel_r, common_genes, ql, rl)
    if (nrow(df) == 0) { message("  No results  -  skipping plots"); next }
    sig <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$delta_velocity) > 0.01, ]
    message("  Significant (FDR<0.05): ", nrow(sig), " genes")
    write.csv(df,  file.path(outdir, paste0(lbl, "_differential_velocity.csv")),  row.names = FALSE)
    write.csv(sig, file.path(outdir, paste0(lbl, "_significant_diff_velocity.csv")), row.names = FALSE)

    kin_q <- seu_q@misc$velocity_kinetics; kin_r <- seu_r@misc$velocity_kinetics
    if (!is.null(kin_q) && !is.null(kin_r)) {
      shared_kin <- intersect(setdiff(colnames(kin_q),"gene"), setdiff(colnames(kin_r),"gene"))
      if (length(shared_kin) > 0) {
        kq <- kin_q[, c("gene", shared_kin), drop = FALSE]; kr <- kin_r[, c("gene", shared_kin), drop = FALSE]
        colnames(kq)[-1] <- paste0(shared_kin,"_query"); colnames(kr)[-1] <- paste0(shared_kin,"_ref")
        df <- merge(df, kq, by = "gene", all.x = TRUE)
        df <- merge(df, kr, by = "gene", all.x = TRUE)
        for (kp in shared_kin) df[[paste0("delta_",kp)]] <- df[[paste0(kp,"_query")]] - df[[paste0(kp,"_ref")]]
        sig <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$delta_velocity) > 0.01, ]
        message("  Kinetics params joined: ", paste(shared_kin, collapse = ", "))
      }
    }

    .plot_kinetics_scatter(df, ql, rl, lbl, outdir)
    .pseudobulk_velocity_lm(seu_q, seu_r, ql, rl, cond_objects, outdir)
    run_latent_time_ks_tests(seu_q, seu_r, ql, rl, outdir)

    .jaccard_velocity_genes(seu_q, seu_r, ql, rl, outdir)

    if (DIFF_KINETICS_ENABLED) {
      tryCatch({
        scv <- reticulate::import("scvelo",  convert = FALSE)
        sc  <- reticulate::import("scanpy",  convert = FALSE)
        .seu_to_adata <- function(seu_x, cond_lbl) {
          vel_mat   <- .get_layer(seu_x, "velocity", "counts")
          spl_mat   <- .get_layer(seu_x, "RNA", "counts")
          unspl_mat <- .get_layer(seu_x, "unspliced", "counts")
          if (is.null(vel_mat) || is.null(spl_mat) || is.null(unspl_mat))
            stop("Missing required assay layer for AnnData conversion")
          cg  <- intersect(rownames(vel_mat), intersect(rownames(spl_mat), rownames(unspl_mat)))
          if (length(cg) < 3) stop("Fewer than 3 common genes for differential kinetics")
          cc  <- colnames(spl_mat)
          spl <- t(as.matrix(spl_mat[cg, cc])); unsp <- t(as.matrix(unspl_mat[cg, cc]))
          vel <- t(as.matrix(vel_mat[cg, cc]))
          adata <- sc$AnnData(X = reticulate::r_to_py(spl),
            layers = reticulate::dict(spliced = reticulate::r_to_py(spl),
                                      unspliced = reticulate::r_to_py(unsp),
                                      velocity  = reticulate::r_to_py(vel)))
          adata$obs_names <- reticulate::r_to_py(cc); adata$var_names <- reticulate::r_to_py(cg)
          obs_df <- seu_x@meta.data[cc, , drop = FALSE]; obs_df$condition_label <- cond_lbl
          adata$obs <- reticulate::r_to_py(obs_df); adata
        }
        adata_joint <- sc$concat(list(.seu_to_adata(seu_q, ql), .seu_to_adata(seu_r, rl)),
          join = "inner", keys = reticulate::r_to_py(list(ql, rl)), label = "batch")
        n_joint    <- as.integer(reticulate::py_to_r(adata_joint$n_obs))
        n_pcs_dk   <- as.integer(min(30L, n_joint - 1L))
        n_neigh_dk <- as.integer(min(30L, max(5L, floor(sqrt(n_joint)))))
        scv$pp$moments(adata_joint, n_pcs = n_pcs_dk, n_neighbors = n_neigh_dk)
        scv$tl$recover_dynamics(adata_joint, n_jobs = N_WORKERS)
        scv$tl$velocity(adata_joint, mode = "dynamical")
        scv$tl$differential_kinetic_test(adata_joint, groupby = "batch")
        dk_var <- reticulate::py_to_r(adata_joint$var)
        if (all(c("fit_diff_kinetics","fit_pval_kinetics") %in% colnames(dk_var))) {
          dk_out <- dk_var[, c("fit_diff_kinetics","fit_pval_kinetics"), drop = FALSE]
          write.csv(dk_out[order(dk_out$fit_pval_kinetics), ],
                    file.path(outdir, paste0(lbl, "_differential_kinetics.csv")))
          message("  [diff kinetics] Saved.")
        }
      }, error = function(e) message("  [diff kinetics] Skipped: ", conditionMessage(e)))
    }

    .plot_diff_results(df, sig, ql, rl, lbl, outdir)
    diff_results[[lbl]] <- df
    ct_dfs <- .diff_velocity_by_celltype(seu_q, seu_r, ql, rl, outdir)
    diff_by_ct[[lbl]] <- ct_dfs
    tryCatch(saveRDS(list(df = df, ct_dfs = ct_dfs), diff_cache_path),
             error = function(e) message("  [CHECKPOINT] Could not save diff cache: ", e$message))
    message("  [CHECKPOINT] Saved: ", diff_cache_path)
  }

  if (length(diff_results) >= 2) {
    k1 <- names(diff_results)[1]; k2 <- names(diff_results)[2]
    s1 <- diff_results[[k1]]$gene[!is.na(diff_results[[k1]]$padj) & diff_results[[k1]]$padj < 0.05]
    s2 <- diff_results[[k2]]$gene[!is.na(diff_results[[k2]]$padj) & diff_results[[k2]]$padj < 0.05]
    shared <- intersect(s1, s2)
    message("  Shared significant genes (first two comparisons): ", length(shared))
    if (length(shared) > 0)
      write.csv(data.frame(gene = shared), file.path(outdir, "shared_significant_genes.csv"), row.names = FALSE)
  }
  has_velocity <- any(vapply(cond_objects, function(seu)
    "velocity_length" %in% colnames(seu@meta.data), logical(1)))
  if (has_velocity) .compare_velocity_metrics(cond_objects, outdir)
  list(diff_results = diff_results, diff_by_ct = diff_by_ct)
}

.filter_ribo_mito <- function(de_res) {
  if (is.null(de_res) || nrow(de_res) == 0) return(de_res)
  rm_idx <- grepl("^(Rpl|Rps|RPL|RPS|mt-|MT-)", rownames(de_res), ignore.case = TRUE)
  if (any(rm_idx)) message("    Filtered ", sum(rm_idx), " ribo/mito genes")
  de_res[!rm_idx, , drop = FALSE]
}

.to_entrez <- function(syms) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db",  quietly = TRUE)) return(character(0))
  syms <- unique(as.character(syms))
  out  <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                           keys = syms, column = "ENTREZID",
                                           keytype = "SYMBOL", multiVals = "first")),
    error = function(e) setNames(rep(NA_character_, length(syms)), syms))
  stats::na.omit(out)
}

.save_enrich_plots_v <- function(enr, stub, outdir, title_pfx = "") {
  df <- tryCatch(as.data.frame(enr), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  tryCatch({
    p_bar <- enrichplot::barplot(enr, showCategory = 10) +
      ggtitle(paste0(title_pfx, " - Barplot")) + .theme_velo(base_size = 11)
    save_plot(p_bar, file.path(outdir, paste0(stub, "_barplot")), 14, 9)
  }, error = function(e) NULL)
  tryCatch({
    p_dot <- enrichplot::dotplot(enr, showCategory = 10) +
      ggtitle(paste0(title_pfx, " - Dotplot")) + .theme_velo(base_size = 11)
    save_plot(p_dot, file.path(outdir, paste0(stub, "_dotplot")), 14, 9)
  }, error = function(e) NULL)
  if (nrow(df) >= 3) tryCatch({
    enr_ts <- suppressWarnings(enrichplot::pairwise_termsim(enr))
    p_emap <- enrichplot::emapplot(enr_ts, showCategory = min(10, nrow(df))) +
      ggtitle(paste0(title_pfx, " - Enrichment Map")) + .theme_velo(base_size = 10)
    save_plot(p_emap, file.path(outdir, paste0(stub, "_emapplot")), 14, 11)
  }, error = function(e) NULL)
  invisible(NULL)
}

.run_go_enrichment_v <- function(de_res, all_genes, outdir, file_stem) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(invisible(NULL))
  if (is.null(de_res) || nrow(de_res) == 0) return(invisible(NULL))
  de_res <- de_res[!is.na(de_res$p_val_adj) & !is.na(de_res$avg_log2FC), , drop = FALSE]
  uni_ent  <- .to_entrez(all_genes)
  up_syms  <- rownames(de_res)[de_res$avg_log2FC >  0.25 & de_res$p_val_adj < 0.1]
  dn_syms  <- rownames(de_res)[de_res$avg_log2FC < -0.25 & de_res$p_val_adj < 0.1]
  up_ent   <- .to_entrez(up_syms); dn_ent <- .to_entrez(dn_syms)
  rnk_sym  <- setNames(de_res$avg_log2FC, rownames(de_res))
  rnk_sym  <- sort(tapply(rnk_sym, names(rnk_sym), function(v) v[which.max(abs(v))]), decreasing = TRUE)
  sym2ent  <- .to_entrez(names(rnk_sym))
  rnk_ent  <- sort(tapply(setNames(de_res$avg_log2FC, sym2ent[rownames(de_res)]),
                           sym2ent[rownames(de_res)],
                           function(v) v[which.max(abs(v))]), decreasing = TRUE)
  rnk_ent  <- rnk_ent[!is.na(names(rnk_ent))]
  ora_dir  <- file.path(outdir, "ORA");  dir.create(ora_dir,  showWarnings = FALSE, recursive = TRUE)
  for (ont in c("BP", "MF", "CC")) {
    if (length(up_ent) >= 5) tryCatch({
      ego_up <- clusterProfiler::enrichGO(gene = up_ent, universe = uni_ent,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.1, readable = TRUE)
      if (nrow(as.data.frame(ego_up)) > 0) {
        write.csv(as.data.frame(ego_up), file.path(ora_dir, sprintf("%s_ORA_GO_%s_UP.csv", file_stem, ont)), row.names = FALSE)
        .save_enrich_plots_v(ego_up, sprintf("%s_ORA_GO_%s_UP", file_stem, ont), ora_dir, paste0("ORA GO-", ont, " UP"))
      }
    }, error = function(e) NULL)
    if (length(dn_ent) >= 5) tryCatch({
      ego_dn <- clusterProfiler::enrichGO(gene = dn_ent, universe = uni_ent,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.1, readable = TRUE)
      if (nrow(as.data.frame(ego_dn)) > 0) {
        write.csv(as.data.frame(ego_dn), file.path(ora_dir, sprintf("%s_ORA_GO_%s_DOWN.csv", file_stem, ont)), row.names = FALSE)
        .save_enrich_plots_v(ego_dn, sprintf("%s_ORA_GO_%s_DOWN", file_stem, ont), ora_dir, paste0("ORA GO-", ont, " DOWN"))
      }
    }, error = function(e) NULL)
  }
  invisible(NULL)
}

.run_pseudobulk_v <- function(seu_q, seu_r, group1, group2, pb_dir, enrich_dir) {
  tryCatch({
    n_q <- ncol(seu_q); n_r <- ncol(seu_r)
    seu_joint <- merge(seu_q, seu_r, add.cell.ids = c(group1, group2))
    seu_joint$group_pb <- c(rep(group1, n_q), rep(group2, n_r))
    rep_col <- if ("replicate" %in% colnames(seu_joint@meta.data)) "replicate" else "orig.ident"

    agg  <- AggregateExpression(seu_joint, assays = "RNA",
                                 group.by = c(rep_col, "group_pb"),
                                 slot = "counts", return.seurat = FALSE)
    bulk <- if (is.list(agg)) as.matrix(agg[[1]]) else as.matrix(agg)

    parsed     <- strsplit(colnames(bulk), "_")
    grp_parsed <- vapply(parsed, function(x) tail(x, 1), character(1))
    rep_parsed <- vapply(parsed, function(x) paste(x[-length(x)], collapse = "_"), character(1))
    colnames(bulk) <- make.names(paste0(rep_parsed, "_", grp_parsed), unique = TRUE)
    coldata <- data.frame(replicate = rep_parsed,
                          group     = factor(grp_parsed),
                          row.names = colnames(bulk))
    write.csv(as.data.frame(bulk), file.path(pb_dir, "Pseudobulk_counts_per_sample.csv"), row.names = TRUE)
    rep_counts <- table(coldata$group)
    file_stem  <- paste0("Pseudobulk_", group1, "_vs_", group2)

    if (all(rep_counts >= 2) && requireNamespace("DESeq2", quietly = TRUE)) {
      has_rep_variation <- tryCatch(length(unique(coldata$replicate)) > 1, error = function(e) FALSE)
      design_formula <- if (has_rep_variation) ~ replicate + group else ~ group
      message("    Pseudo-bulk: DESeq2 (design: ", deparse(design_formula), ")")
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(bulk),
                                             colData = coldata,
                                             design  = design_formula)
      dds <- dds[rowSums(DESeq2::counts(dds) >= 3) >= 2, ]
      dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE),
                      error = function(e) { message("DESeq2 failed: ", e$message); NULL })
      if (is.null(dds)) return(invisible(NULL))
      res    <- DESeq2::results(dds, contrast = c("group", group1, group2))
      res_df <- as.data.frame(res)
      res_df$avg_log2FC <- res_df$log2FoldChange
      res_df$p_val_adj  <- res_df$padj
      res_df$gene       <- rownames(res_df)
      res_df <- .filter_ribo_mito(res_df)
      write.csv(res_df, file.path(pb_dir, "Pseudobulk_DESeq2_results.csv"), row.names = FALSE)
      if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
        res_plot <- res_df[!is.na(res_df$p_val_adj), ]
        up_n <- sum(res_plot$avg_log2FC >  0.25 & res_plot$p_val_adj < 0.1, na.rm = TRUE)
        dn_n <- sum(res_plot$avg_log2FC < -0.25 & res_plot$p_val_adj < 0.1, na.rm = TRUE)
        p_vol <- tryCatch(EnhancedVolcano::EnhancedVolcano(res_plot, lab = rownames(res_plot),
          x = "avg_log2FC", y = "p_val_adj",
          title    = paste("Pseudo-bulk DESeq2:", group1, "vs", group2),
          subtitle = sprintf("Up: %d  |  Down: %d", up_n, dn_n),
          pCutoff = 0.1, FCcutoff = 0.25, labSize = 3, drawConnectors = TRUE,
          col = c("grey70", AAAS_COLS[5], AAAS_COLS[1], AAAS_COLS[2])
        ) + .theme_velo(), error = function(e) NULL)
        if (!is.null(p_vol))
          save_plot(p_vol, file.path(pb_dir, "Pseudobulk_Volcano"), PLT$w_single, PLT$h_single)
      }
      .run_go_enrichment_v(res_df, rownames(bulk), enrich_dir, file_stem)
    } else {
      message("    Pseudo-bulk: CPM fold-change fallback  -  fewer than 2 replicates per group detected")
      lib  <- colSums(bulk); lib[lib == 0] <- 1
      cpm  <- sweep(bulk, 2, lib / 1e6, "/")
      keep <- rowSums(cpm >= 1) >= 1; cpm <- cpm[keep, ]
      idx1 <- which(coldata$group == group1); idx2 <- which(coldata$group == group2)
      lfc  <- log2((rowMeans(cpm[, idx1, drop = FALSE]) + 1) /
                   (rowMeans(cpm[, idx2, drop = FALSE]) + 1))
      res_df <- data.frame(avg_log2FC = lfc,
                           p_val_adj  = ifelse(abs(lfc) >= 1, 0.01, 0.5),
                           gene = names(lfc), row.names = names(lfc))
      res_df <- .filter_ribo_mito(res_df)
      write.csv(res_df, file.path(pb_dir, "Pseudobulk_FoldChange_results.csv"), row.names = FALSE)
      writeLines(
        c("NOTE: CPM fold-change fallback used  -  fewer than 2 replicates per group were detected.",
          "With 2 replicates per condition this should NOT occur. Check that 'replicate' metadata",
          "is present and that AggregateExpression column-parsing succeeds."),
        file.path(pb_dir, "Analysis_Note.txt"))
      .run_go_enrichment_v(res_df, rownames(bulk), enrich_dir, file_stem)
    }
  }, error = function(e) message("    Pseudo-bulk failed: ", e$message))
  invisible(NULL)
}

run_de_analysis <- function(cond_objects) {
  message("\n", strrep("=", 70))
  message("DIFFERENTIAL GENE EXPRESSION  (replicate-aware pseudobulk + per-CT Wilcoxon)")
  message(strrep("=", 70))
  de_outdir     <- file.path(OUTDIR, "de")
  enrich_outdir <- file.path(OUTDIR, "enrichment")
  pb_outdir     <- file.path(OUTDIR, "pseudobulk")
  for (.d in c(de_outdir, enrich_outdir, pb_outdir))
    dir.create(.d, showWarnings = FALSE, recursive = TRUE)

  all_de_results <- list()   # for velocity-expression integration

  for (comp in COMPARISONS) {
    group1 <- comp[1]; group2 <- comp[2]
    comp_lbl <- paste0(group1, "_vs_", group2)
    message("\n--- DE Comparison: ", comp_lbl, " ---")
    seu_q <- cond_objects[[group1]]; seu_r <- cond_objects[[group2]]
    if (is.null(seu_q) || is.null(seu_r)) { message("  Missing condition object  -  skipping"); next }

    comp_de_dir       <- file.path(de_outdir,              comp_lbl)
    comp_enr_pb_dir   <- file.path(enrich_outdir, "pseudo_bulk", comp_lbl)
    comp_enr_de_dir   <- file.path(enrich_outdir, "de",          comp_lbl)
    comp_pb_dir       <- file.path(pb_outdir,              comp_lbl)
    for (.d in c(comp_de_dir, comp_enr_pb_dir, comp_enr_de_dir, comp_pb_dir))
      dir.create(.d, showWarnings = FALSE, recursive = TRUE)

    de_cache <- file.path(comp_de_dir, paste0(comp_lbl, "_de_cache.rds"))
    if (!FORCE_RERUN && file.exists(de_cache)) {
      message("  [CHECKPOINT] Loading cached DE results  -  skipping recomputation")
      all_de_results[[comp_lbl]] <- readRDS(de_cache)
      next
    }

    # NOTE: Avoid merging the two full condition objects up-front  -  for large
    # datasets with multiple assays (RNA + spliced + unspliced + reductions),
    # Seurat::merge() can blow past size_t limits in C++ (`vector::reserve`)
    # and crash the R session even when overall RAM is sufficient. Instead,
    # determine cell types from the source objects, then subset and merge
    # one cell type at a time on RNA-only diet copies.
    ct_col <- if ("cell_type" %in% colnames(seu_q@meta.data) &&
                  "cell_type" %in% colnames(seu_r@meta.data)) "cell_type" else NULL
    all_de_cache <- list()

    if (!is.null(ct_col)) {
      cell_types <- setdiff(union(unique(seu_q@meta.data[[ct_col]]),
                                  unique(seu_r@meta.data[[ct_col]])),
                            c("Unknown", "Ambiguous", "Unassigned"))
      message("  Cell types to test: ", length(cell_types))
      for (ct in cell_types) {
        q_cells <- colnames(seu_q)[seu_q@meta.data[[ct_col]] == ct]
        r_cells <- colnames(seu_r)[seu_r@meta.data[[ct_col]] == ct]
        n1 <- length(q_cells); n2 <- length(r_cells)
        if (n1 < 3 || n2 < 3) next
        tryCatch({
          # Strip to RNA assay only, drop reductions/graphs, then merge the
          # small per-CT subsets  -  keeps memory footprint bounded.
          seu_q_ct <- DietSeurat(subset(seu_q, cells = q_cells),
                                 assays = "RNA", dimreducs = NULL, graphs = NULL)
          seu_r_ct <- DietSeurat(subset(seu_r, cells = r_cells),
                                 assays = "RNA", dimreducs = NULL, graphs = NULL)
          seu_ct <- merge(seu_q_ct, seu_r_ct, add.cell.ids = c(group1, group2))
          seu_ct$group_de <- c(rep(group1, n1), rep(group2, n2))
          Idents(seu_ct) <- "group_de"
          DefaultAssay(seu_ct) <- "RNA"
          # Seurat 5 layered counts must be joined before FindMarkers
          if (inherits(seu_ct[["RNA"]], "Assay5"))
            seu_ct <- JoinLayers(seu_ct, assay = "RNA")

          de_res <- FindMarkers(seu_ct, ident.1 = group1, ident.2 = group2,
                                test.use = "wilcox", logfc.threshold = 0.1,
                                min.pct = 0.1, verbose = FALSE,
                                max.cells.per.ident = 300)
          de_res <- .filter_ribo_mito(de_res)
          if (nrow(de_res) > 0) {
            de_res$gene <- rownames(de_res)
            ct_safe   <- gsub("[^[:alnum:]]", "_", ct)
            ct_de_dir <- file.path(comp_de_dir,     ct_safe)
            ct_en_dir <- file.path(comp_enr_de_dir, ct_safe)
            dir.create(ct_de_dir, showWarnings = FALSE, recursive = TRUE)
            dir.create(ct_en_dir, showWarnings = FALSE, recursive = TRUE)
            write.csv(de_res, file.path(ct_de_dir, "DE_results.csv"), row.names = FALSE)
            up_n <- sum(de_res$avg_log2FC >  0.25 & de_res$p_val_adj < 0.1, na.rm = TRUE)
            dn_n <- sum(de_res$avg_log2FC < -0.25 & de_res$p_val_adj < 0.1, na.rm = TRUE)
            message("    ", ct, ": n=", n1, "/", n2, "  up=", up_n, "  down=", dn_n)
            if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
              de_plot <- de_res[!is.na(de_res$p_val_adj), ]
              p_vol <- tryCatch(EnhancedVolcano::EnhancedVolcano(de_plot, lab = rownames(de_plot),
                x = "avg_log2FC", y = "p_val_adj",
                title = paste(ct, ":", group1, "vs", group2),
                subtitle = sprintf("Up: %d  |  Down: %d", up_n, dn_n),
                pCutoff = 0.1, FCcutoff = 0.25, labSize = 3, drawConnectors = TRUE,
                col = c("grey70", AAAS_COLS[5], AAAS_COLS[1], AAAS_COLS[2])
              ) + .theme_velo(), error = function(e) NULL)
              if (!is.null(p_vol))
                save_plot(p_vol, file.path(ct_de_dir, "Volcano"), PLT$w_single, PLT$h_single)
            }
            .run_go_enrichment_v(de_res, rownames(seu_ct), ct_en_dir,
                                 paste(ct_safe, comp_lbl, sep = "_"))
            all_de_cache[[ct]] <- de_res
          }
          rm(seu_q_ct, seu_r_ct, seu_ct); gc(verbose = FALSE)
        }, error = function(e) message("    DE failed for ", ct, ": ", e$message))
      }
    }

    message("  Running replicate-aware pseudobulk (", comp_lbl, ") ...")
    .run_pseudobulk_v(seu_q, seu_r, group1, group2, comp_pb_dir, comp_enr_pb_dir)

    tryCatch(saveRDS(all_de_cache, de_cache),
             error = function(e) message("  [CHECKPOINT] DE cache save failed: ", e$message))
    message("  [CHECKPOINT] DE saved: ", de_cache)
    all_de_results[[comp_lbl]] <- all_de_cache
  }
  message("\nDE analysis complete. Results in: de/  enrichment/  pseudobulk/")
  all_de_results
}

run_qc_metrics <- function(rep_objects, cond_objects) {
  message("\n", strrep("=", 70))
  message("ALL-COMBINED PLOTS")
  message(strrep("=", 70))
  outdir <- file.path(OUTDIR, "qc_metrics")
  .make_dim_panel <- function(seu, title, reduction, group_col, pal) {
    if (!reduction %in% names(seu@reductions)) return(plot_spacer())
    if (!group_col %in% colnames(seu@meta.data)) return(plot_spacer())
    DimPlot(seu, reduction = reduction, group.by = group_col,
            cols = pal, pt.size = PLT$pt_size * 0.5, label = TRUE, repel = TRUE, label.size = 3.5) +
      ggtitle(title) + NoLegend() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  }
  sample_order <- names(SAMPLES)
  all_ct <- unique(unlist(lapply(c(rep_objects, as.list(cond_objects)),
                                  function(s) unique(s$cluster_majority_type))))
  all_ct <- all_ct[!is.na(all_ct)]
  maj_pal_global <- .build_ct_palette(all_ct)
  all_celltypes  <- unique(unlist(lapply(c(rep_objects, as.list(cond_objects)),
                                          function(s) unique(s$cell_type))))
  ct_pal_global  <- .build_ct_palette(all_celltypes)

  message("  Combined sample UMAP+t-SNE majority grid ...")
  sample_panels <- lapply(sample_order, function(sn) {
    seu <- rep_objects[[sn]]; if (is.null(seu)) return(list(u = plot_spacer(), t = plot_spacer()))
    lbl <- SAMPLES[[sn]]$label
    grp <- if ("cluster_majority_type" %in% colnames(seu@meta.data)) "cluster_majority_type" else "cell_type"
    pal <- if (grp == "cluster_majority_type") maj_pal_global else .build_ct_palette(unique(seu$cell_type))
    list(u = .make_dim_panel(seu, paste0(lbl,"  -  UMAP"), "umap", grp, pal),
         t = .make_dim_panel(seu, paste0(lbl,"  -  t-SNE"), "tsne", grp, pal))
  })
  umap_row <- wrap_plots(lapply(sample_panels, `[[`, "u"), nrow = 1)
  tsne_row <- wrap_plots(lapply(sample_panels, `[[`, "t"), nrow = 1)
  save_plot((umap_row / tsne_row) + plot_annotation(title = "All samples  -  UMAP & t-SNE",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))),
    file.path(outdir, "samples_umap_tsne_majority_side_by_side"),
    PLT$w_triple, PLT$h_double, limitsize = FALSE)

  message("  Combined condition UMAP+t-SNE grid ...")
  cond_panels <- lapply(CONDITIONS, function(cond) {
    seu <- cond_objects[[cond]]; if (is.null(seu)) return(list(u = plot_spacer(), t = plot_spacer()))
    grp <- if ("cluster_majority_type" %in% colnames(seu@meta.data)) "cluster_majority_type" else "cell_type"
    pal <- if (grp == "cluster_majority_type") maj_pal_global else .build_ct_palette(unique(seu$cell_type))
    list(u = .make_dim_panel(seu, paste0(cond,"  -  UMAP"), "umap", grp, pal),
         t = .make_dim_panel(seu, paste0(cond,"  -  t-SNE"), "tsne", grp, pal))
  })
  cond_u_row <- wrap_plots(lapply(cond_panels, `[[`, "u"), nrow = 1)
  cond_t_row <- wrap_plots(lapply(cond_panels, `[[`, "t"), nrow = 1)
  save_plot((cond_u_row / cond_t_row) + plot_annotation(title = "All conditions  -  UMAP & t-SNE",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))),
    file.path(outdir, "conditions_umap_tsne_majority_side_by_side"),
    PLT$w_triple, PLT$h_single, limitsize = FALSE)

  message("  Combined stacked composition plots ...")
  sample_labels <- vapply(sample_order, function(s) SAMPLES[[s]]$label, character(1))
  .stacked_with_labels <- function(objects, names_vec, title, outpath, table_path) {
    meta_list <- lapply(seq_along(objects), function(i) {
      df <- objects[[i]]@meta.data
      if (!"cell_type" %in% colnames(df)) return(NULL)
      data.frame(grp = names_vec[i], cell_type = unname(df$cell_type), stringsAsFactors = FALSE)
    })
    meta_all <- do.call(rbind, Filter(Negate(is.null), meta_list))
    if (nrow(meta_all) == 0) return(invisible(NULL))
    meta_all$grp <- factor(meta_all$grp, levels = unique(names_vec))
    df2 <- meta_all %>%
      dplyr::group_by(grp, cell_type) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(grp) %>% dplyr::mutate(pct = n / sum(n) * 100) %>% dplyr::ungroup() %>%
      dplyr::arrange(grp, dplyr::desc(pct))
    df2$label_txt <- ifelse(df2$pct >= 2, sprintf("%.1f%%", df2$pct), "")
    p <- ggplot(df2, aes(x = grp, y = pct, fill = cell_type)) +
      geom_col(position = "stack", width = 0.82) +
      geom_text(aes(label = label_txt), position = position_stack(vjust = 0.5),
                size = 3.2, fontface = "bold", colour = "white") +
      scale_fill_manual(values = ct_pal_global) +
      labs(x = NULL, y = "% cells", title = title, fill = "Cell type") +
      .theme_velo(base_size = 14) +
      theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 13, face = "bold")) +
      guides(fill = guide_legend(ncol = 1))
    save_plot(p, outpath, max(10, length(unique(names_vec)) * 2.5),
              max(12, length(unique(df2$cell_type)) * 0.35 + 6), limitsize = FALSE)
    write.csv(df2 %>% dplyr::select(grp, cell_type, pct) %>%
                tidyr::pivot_wider(names_from = grp, values_from = pct, values_fill = 0),
              table_path, row.names = FALSE)
  }
  .stacked_with_labels(rep_objects[sample_order], sample_labels,
    "All samples  -  Cell type proportions",
    file.path(outdir, "samples_celltype_composition"),
    file.path(outdir, "samples_celltype_composition_metrics.csv"))
  .stacked_with_labels(as.list(cond_objects), CONDITIONS,
    "All conditions  -  Cell type proportions",
    file.path(outdir, "conditions_celltype_composition"),
    file.path(outdir, "conditions_celltype_composition_metrics.csv"))
  message("  All-combined plots saved to: ", outdir)
}

write_summary <- function(rep_objects, cond_objects, diff_list) {
  diff_results <- diff_list$diff_results; diff_by_ct <- diff_list$diff_by_ct
  lines <- c("RNA Velocity R Pipeline (Refactored)  -  Summary Report", strrep("=", 60), "", "Per-replicate:")
  for (s in names(rep_objects)) {
    seu <- rep_objects[[s]]; m <- SAMPLES[[s]]
    lines <- c(lines, sprintf("  %-12s  cond=%-8s  rep=%s  cells=%d  genes=%d",
                               m$label, m$condition, m$replicate, ncol(seu), nrow(seu)))
  }
  lines <- c(lines, "", "Per-condition:")
  for (cond in names(cond_objects)) {
    seu <- cond_objects[[cond]]
    lines <- c(lines, sprintf("  %-8s  cells=%d  genes=%d", cond, ncol(seu), nrow(seu)))
  }
  lines <- c(lines, "", "Differential velocity (FDR < 0.05):")
  for (comp in names(diff_results)) {
    df  <- diff_results[[comp]]
    sig <- df[!is.na(df$padj) & df$padj < 0.05, ]
    lines <- c(lines, sprintf("  %s: %d sig genes  (up query: %d, up ref: %d)",
                               comp, nrow(sig),
                               sum(sig$delta_velocity > 0, na.rm = TRUE),
                               sum(sig$delta_velocity < 0, na.rm = TRUE)))
  }
  lines <- c(lines, "", "DE analysis: see de/  enrichment/  pseudobulk/  (replicate-aware DESeq2)")
  lines <- c(lines, "", "Per-cell-type velocity differential:")
  for (comp in names(diff_by_ct)) {
    lines <- c(lines, paste0("  ", comp, ":"))
    ct_dfs <- diff_by_ct[[comp]]
    if (length(ct_dfs) == 0) { lines <- c(lines, "    (no results)"); next }
    for (ct in names(ct_dfs)) {
      df_ct  <- ct_dfs[[ct]]
      sig_ct <- df_ct[!is.na(df_ct$padj) & df_ct$padj < 0.05 & abs(df_ct$delta_velocity) > 0.01, ]
      lines  <- c(lines, sprintf("    %-40s  %4d sig genes", ct, nrow(sig_ct)))
    }
  }
  lines <- c(lines, "",
    "Annotation strategy:",
    "  Level-1: all cells classified with MOUSE_MARKERS (whole-embryo panel)",
    "  Level-2: Neural_crest subset re-embedded and annotated with Neural_crest_L2",
    "           sub-states: Premigratory_NC, Migratory_NC, Cranial_NC, Cardiac_NC, Trunk_NC",
    "",
    "Output directories:",
    "  per_replicate/          -- per-sample velocity, UMAP, dot/violin/feature plots",
    "  per_condition/          -- per-condition merged objects (Harmony batch correction)",
    "  splicing_efficiency/    -- unspliced/(unspliced+spliced) index per gene x cell type",
    "  velocity_pseudotime_corr/ -- velocity ~ pseudotime correlation",
    "  celltype_transitions/   -- Markov transition probability heatmaps",
    "  gene_velocity_coherence/ -- gene-level velocity coherence per condition",
    "  latent_time_ks_factorial/ -- latent time KS tests (cell type x comparison)",
    "  differentially_variable/ -- scran::modelGeneVar DVGs between conditions",
    "  differential/           -- velocity DE, expression DE, discordant GO enrichment",
    "  de/                     -- replicate-aware DESeq2 + per-cell-type Wilcoxon",
    "  composition_stats/      -- CLR-LME composition + ranked bar plots",
    "  replicate_concordance/  -- expression/velocity correlation, BA plots",
    "  splicing_ratios/        -- spliced/unspliced ratio analysis")
  report <- paste(lines, collapse = "\n")
  cat(report, "\n")
  writeLines(report, file.path(OUTDIR, "R_analysis_summary.txt"))
}

# =============================================================================
# 7. Main
# =============================================================================
message("Starting Neural Crest RNA Velocity Pipeline")
message("Strategy    : Level-1 whole-embryo annotation -> Neural_crest subset -> Level-2 NC sub-states")
message("CellRanger  : ", CELLRANGER_DIR)
message("Output      : ", OUTDIR)

# ── Main pipeline ─────────────────────────────────────────────────────────────
rep_objects  <- run_per_replicate()
cond_objects <- run_per_condition()

# Per-replicate QC metrics
compute_replicate_qc(rep_objects, file.path(OUTDIR, "per_replicate"))

# Replicate concordance analysis
run_replicate_concordance(rep_objects, file.path(OUTDIR, "replicate_concordance"))

# Splicing ratio analysis (per-condition objects)
run_splicing_ratio_analysis(as.list(cond_objects), names(cond_objects),
                             file.path(OUTDIR, "splicing_ratios"))

run_composition_stats(rep_objects, cond_objects,
                      file.path(OUTDIR, "composition_stats"))

plot_composition_ranked_bars(rep_objects, cond_objects,
                             file.path(OUTDIR, "composition_stats"))

# Differential velocity (MWU + pseudobulk-LM + per-CT + latent-time KS)
diff_list    <- run_differential(cond_objects)

# Replicate-aware DESeq2 DE + Wilcoxon per cell type + GO enrichment
de_results   <- run_de_analysis(cond_objects)

# Velocity gene programs (per condition)
for (cond in names(cond_objects)) {
  tryCatch(
    run_velocity_gene_programs(cond_objects[[cond]], cond,
                                file.path(OUTDIR, "velocity_programs")),
    error = function(e) message("  [velocity_programs] ERROR for ", cond, ": ", e$message))
}

# Gene-level velocity coherence per condition
tryCatch(
  compute_gene_velocity_coherence(as.list(cond_objects), names(cond_objects),
                                   file.path(OUTDIR, "gene_velocity_coherence")),
  error = function(e) message("  [gene_velocity_coherence] ERROR: ", e$message))

# Latent time KS test  -  factorial summary (cell type  x  condition)
tryCatch(
  run_latent_time_ks_factorial(cond_objects,
                                file.path(OUTDIR, "latent_time_ks_factorial")),
  error = function(e) message("  [latent_time_ks_factorial] ERROR: ", e$message))

# Differentially variable genes (scran::modelGeneVar) between conditions
tryCatch(
  run_differentially_variable_genes(cond_objects,
                                     file.path(OUTDIR, "differentially_variable")),
  error = function(e) message("  [differentially_variable] ERROR: ", e$message))

for (comp in COMPARISONS) {
  tryCatch({
    ql <- comp[1]; rl <- comp[2]; lbl <- paste0(ql, "_vs_", rl)
    vel_df <- diff_list$diff_results[[lbl]]
    de_df  <- tryCatch({
      cache_path <- file.path(OUTDIR, "de", lbl, paste0(lbl, "_de_cache.rds"))
      if (file.exists(cache_path)) {
        cache <- readRDS(cache_path)
        raw   <- do.call(rbind, Filter(is.data.frame, cache))
        if (!is.null(raw) && nrow(raw) > 0 && "gene" %in% colnames(raw)) {
          raw <- raw[order(raw$p_val_adj), ]
          raw[!duplicated(raw$gene), ]
        } else NULL
      } else NULL
    }, error = function(e) NULL)
    if (!is.null(vel_df) && !is.null(de_df)) {
      run_velocity_expression_integration(
        cond_objects[[ql]], cond_objects[[rl]],
        vel_df, de_df, ql, rl,
        file.path(OUTDIR, "differential"))
    }
  }, error = function(e) message("  [vel_expr_integration] ERROR for ", paste(comp, collapse="_vs_"), ": ", e$message))
}

# All-combined cross-sample plots
tryCatch(
  run_qc_metrics(rep_objects, cond_objects),
  error = function(e) message("  [qc_metrics] ERROR: ", e$message))

# Custom summary visualisations (Plots 1–10 from generate_plots.R, fully data-driven)
tryCatch(
  plot_custom_summary_visualizations(OUTDIR),
  error = function(e) message("  [custom_summary_plots] ERROR: ", e$message))

# Summary report
tryCatch(
  write_summary(rep_objects, cond_objects, diff_list),
  error = function(e) message("  [write_summary] ERROR: ", e$message))

message("\n=== Refactored R Pipeline complete!  Results: ", OUTDIR, " ===")
      
      
 


