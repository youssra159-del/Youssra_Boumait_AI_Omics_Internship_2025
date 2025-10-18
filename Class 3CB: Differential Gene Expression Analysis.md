# =====================================================================
#               AI and Biotechnology / Bioinformatics
# =====================================================================
#              AI and Omics Research Internship (2025)
# ---------------------------------------------------------------------
#             Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------------------
#                     Microarray Data Analysis
# =====================================================================

# =====================================================================
#  Probe→Gene mapping • limma DEG • Volcano • Heatmap • Summary
#  Input: Illumina_Preprocessed.rds  (from your Part I)
# =====================================================================

gc()
set.seed(42)

## -------------------- Setup & folders --------------------------------
RESULTS_DIR <- "Results2"
PLOTS_DIR   <- file.path(RESULTS_DIR, "Result_Plots2")
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(PLOTS_DIR,   showWarnings = FALSE)

## -------------------- Packages ---------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr","tibble","ggplot2","pheatmap","RColorBrewer")
bioc_pkgs <- c("limma","AnnotationDbi","Biobase","illuminaHumanv3.db","org.Hs.eg.db")

need_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need_cran)) install.packages(need_cran, quiet = TRUE)

need_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(need_bioc)) BiocManager::install(need_bioc, ask = FALSE, update = FALSE)

invisible(lapply(c(cran_pkgs, bioc_pkgs), require, character.only = TRUE))

## -------------------- Load your ExpressionSet -------------------------
# (Run from folder containing Illumina_Preprocessed.rds)
eset_path <- "Illumina_Preprocessed.rds"
stopifnot(file.exists(eset_path))
eset <- readRDS(eset_path)
stopifnot(methods::is(eset, "ExpressionSet"))

expr_mat <- Biobase::exprs(eset)     # features x samples
pheno    <- Biobase::pData(eset)

## --------- Probe→Gene mapping (Illumina-aware) -----------------------
# Detect whether rownames are ILMN probe IDs or ENTREZ IDs
ids <- rownames(expr_mat)
stopifnot(length(ids) > 0)

is_ilmn   <- grepl("^ILMN[_-]?[0-9]+$", ids, ignore.case = TRUE)
is_entrez <- grepl("^[0-9]+$", ids)

if (any(is_ilmn) && !any(!is_ilmn)) {
  message("Detected ILMN probe IDs; mapping via illuminaHumanv3.db (PROBEID → SYMBOL).")
  gene_symbols <- AnnotationDbi::mapIds(illuminaHumanv3.db,
                                        keys     = ids,
                                        keytype  = "PROBEID",
                                        column   = "SYMBOL",
                                        multiVals = "first")
} else if (any(is_entrez) && !any(!is_entrez)) {
  message("Detected ENTREZ IDs; mapping via org.Hs.eg.db (ENTREZID → SYMBOL).")
  gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                        keys     = ids,
                                        keytype  = "ENTREZID",
                                        column   = "SYMBOL",
                                        multiVals = "first")
} else {
  stop("Row names are neither uniform ILMN IDs nor uniform ENTREZ IDs. Inspect head(rownames(expr_mat)).")
}

# Build mapping table and merge with expression
gene_map_df <- tibble::tibble(PROBEID = ids, SYMBOL = unname(gene_symbols))

processed_data_df <- expr_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::left_join(gene_map_df, by = "PROBEID") %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

# Keep only numeric expression columns
expr_only <- processed_data_df %>% dplyr::select(-PROBEID, -SYMBOL)

## --------- Handle multiple probes per gene (collapse duplicates) -----
# Average probes per gene with limma::avereps → one row per gene (SYMBOL)
expr_gene_all <- limma::avereps(as.matrix(expr_only), ID = processed_data_df$SYMBOL)

# Name columns as GSM sample names
colnames(expr_gene_all) <- Biobase::sampleNames(eset)

## --------- Define groups & align samples ------------------------------
# Prefer curated column if present
lvl <- c("Control","Active_TB_BeforeTx","Active_TB_2m","Active_TB_12m")

if ("group_label" %in% colnames(pheno) && any(!is.na(pheno$group_label))) {
  groups <- factor(as.character(pheno$group_label), levels = lvl)
} else {
  # Fallback: derive from sample names (CON / PTB ... long_0/2/12)
  samps <- Biobase::sampleNames(eset)
  groups <- dplyr::case_when(
    grepl("^CON",     samps, ignore.case = TRUE) ~ "Control",
    grepl("long_0_",  samps, ignore.case = TRUE) ~ "Active_TB_BeforeTx",
    grepl("long_2_",  samps, ignore.case = TRUE) ~ "Active_TB_2m",
    grepl("long_12_", samps, ignore.case = TRUE) ~ "Active_TB_12m",
    TRUE ~ NA_character_
  )
  groups <- factor(groups, levels = lvl)
}

# Keep labeled samples and align matrix + factor
keep <- !is.na(groups)
expr_gene <- expr_gene_all[, keep, drop = FALSE]
groups    <- droplevels(groups[keep])

# Optional: order columns for nicer plots: Control → TB0 → TB2 → TB12
ord <- order(groups)
expr_gene <- expr_gene[, ord, drop = FALSE]
groups    <- droplevels(groups[ord])

# Sanity checks
print(table(groups))
stopifnot(ncol(expr_gene) == length(groups))

## ---------------- limma design & contrast ----------------------------
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

contr <- limma::makeContrasts(
  TB0_vs_CTRL  = Active_TB_BeforeTx - Control,
  TB2_vs_TB0   = Active_TB_2m       - Active_TB_BeforeTx,
  TB12_vs_CTRL = Active_TB_12m      - Control,
  levels = design
)

fit  <- limma::lmFit(expr_gene, design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr))

## ---------------- Extract DE for the primary contrast ----------------
coef_name  <- "TB0_vs_CTRL"   # change if you want a different main contrast
deg_results <- limma::topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")

# Classify for plots/files
deg_results$threshold <- with(deg_results,
  ifelse(adj.P.Val < 0.05 & logFC >  1, "Upregulated",
  ifelse(adj.P.Val < 0.05 & logFC < -1, "Downregulated", "No"))
)

up   <- subset(deg_results, threshold == "Upregulated")
down <- subset(deg_results, threshold == "Downregulated")
ud   <- rbind(up, down)

## ---------------- Save outputs (per course) --------------------------
write.csv(deg_results, file.path(RESULTS_DIR, "DEGs_Results.csv"),       row.names = TRUE)
write.csv(up,          file.path(RESULTS_DIR, "Upregulated_DEGs.csv"),   row.names = TRUE)
write.csv(down,        file.path(RESULTS_DIR, "Downregulated_DEGs.csv"), row.names = TRUE)
write.csv(ud,          file.path(RESULTS_DIR, "Updown_DEGs.csv"),        row.names = TRUE)

## ---------------- Volcano (keep your colors) -------------------------
volc_cols <- c(Upregulated="#0ea5b1", Downregulated="#1d4ed8", No="grey70")

png(file.path(PLOTS_DIR, "volcano_plot.png"), width=2000, height=1500, res=300)
ggplot2::ggplot(deg_results, ggplot2::aes(logFC, -log10(adj.P.Val), color = threshold)) +
  ggplot2::geom_point(alpha=.75, size=1.8) +
  ggplot2::scale_color_manual(values=volc_cols) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype="dashed", linewidth=.4) +
  ggplot2::geom_vline(xintercept = c(-1,1), linetype="dotted",  linewidth=.4) +
  ggplot2::theme_minimal(base_size=12) +
  ggplot2::labs(title=paste("Volcano:", coef_name), x="log2 Fold Change", y="-log10(FDR-adjusted P-value)", color="Regulation")
dev.off()

## ---------------- Heatmap Top 25 (keep your gradient) ----------------
top25 <- rownames(ud[order(ud$adj.P.Val), ])[seq_len(min(25, nrow(ud)))]
hm <- expr_gene[top25, , drop=FALSE]
lbl <- as.character(groups)
colnames(hm) <- ave(lbl, lbl, FUN=function(x) paste0(x, "_", seq_along(x)))

png(file.path(PLOTS_DIR, "heatmap_top25_DEGs.png"), width=2000, height=1500, res=300)
pheatmap::pheatmap(hm, scale="none", cluster_rows=FALSE, cluster_cols=TRUE,
                   show_rownames=TRUE, show_colnames=TRUE,
                   color=colorRampPalette(c("#1d4ed8","white","#0ea5b1"))(101),
                   fontsize_row=7, fontsize_col=9,
                   main=paste("Top 25 DEGs:", coef_name))
dev.off()

## ---------------- Short summary  --------------------------
n_up <- nrow(up); n_down <- nrow(down)
writeLines(c(
  "Multiple probes mapped to the same gene; used illuminaHumanv3.db when ILMN probe IDs were present, or ENTREZ→SYMBOL via org.Hs.eg.db otherwise. Duplicates collapsed with limma::avereps (mean).",
  paste0("Differential expression tested for ", gsub("_", " ", coef_name), " using limma (no-intercept design; empirical Bayes moderation)."),
  sprintf("Thresholds: FDR < 0.05 and |log2FC| > 1 → %d upregulated, %d downregulated genes.", n_up, n_down),
  "Exports: DEGs_Results.csv (+ up/down splits), volcano_plot.png, heatmap_top25_DEGs.png.",
  "Groups used: Control, Active_TB_BeforeTx, Active_TB_2m, Active_TB_12m."
), file.path(RESULTS_DIR, "RESULTS_SUMMARY.txt"))

