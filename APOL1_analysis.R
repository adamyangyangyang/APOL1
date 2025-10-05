suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# ================================
# Config
# ================================
counts_fp <- "/workdir/cfrna/alignment/takara_human_V3/output/APOL_urine/APOL_urine_feature_counts.tsv"
out_dir   <- "/workdir/qy236/apol_kidney/QC_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Complexity parameters
K_TOP          <- 1500    # evaluate fraction of reads in top-K genes (lower is better)
COMPLEXITY_MAX <- 0.75    # FAIL if frac_topK > this cutoff

# Output files
QC_TSV   <- file.path(out_dir, sprintf("complexity_qc_K%d_cut%.2f.tsv", K_TOP, COMPLEXITY_MAX))
PLOT_PNG <- file.path(out_dir, sprintf("library_complexity_K%d_cut%.2f.png", K_TOP, COMPLEXITY_MAX))
PASS_COUNTS_TSV <- file.path(out_dir, "APOL_urine_feature_counts_passComplexity.tsv")

# IRB / negative-control core IDs (no "_bp..." suffix)
irb_ids <- c(
  "cfrna_APOL_11","cfrna_APOL_123","cfrna_APOL_139","cfrna_APOL_150",
  "cfrna_APOL_23","cfrna_APOL_39","cfrna_APOL_51","cfrna_APOL_63",
  "cfrna_APOL_75","cfrna_APOL_87","cfrna_APOL_99","cfrna_APOL_111",
  "cfrna_APOL_12","cfrna_APOL_124","cfrna_APOL_140","cfrna_APOL_151",
  "cfrna_APOL_24","cfrna_APOL_40","cfrna_APOL_52","cfrna_APOL_64",
  "cfrna_APOL_76","cfrna_APOL_88","cfrna_APOL_100","cfrna_APOL_112"
)

# ================================
# Load counts (genes x samples)
# ================================
raw_count <- readr::read_tsv(counts_fp, show_col_types = FALSE)
stopifnot("geneID" %in% names(raw_count))
gene_ids <- raw_count$geneID
count_mat <- as.matrix(raw_count[ , -1, drop = FALSE])
rownames(count_mat) <- gene_ids

# Normalize sample IDs (optional but safer)
norm_id <- function(x) gsub("-", "_", gsub("\\.", "_", x))
colnames(count_mat) <- norm_id(colnames(count_mat))

# ================================
# Exclude IRB/negative controls
# ================================
# Column names look like "cfrna_APOL_67_bp150x150"; core id is before "_bp"
col_core <- sub("_bp.*$", "", colnames(count_mat))
keep_col <- !(col_core %in% irb_ids)
if (any(!keep_col)) {
  message("Excluding IRB/neg-control samples: ",
          paste(unique(col_core[!keep_col]), collapse = ", "))
}
count_mat <- count_mat[, keep_col, drop = FALSE]

# Keep genes with any counts to speed up
count_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]

# ================================
# Complexity curves & metric
# ================================
cum_curve <- function(x, k = NULL) {
  if (length(x) == 0 || sum(x, na.rm = TRUE) == 0) {
    return(data.frame(rank = integer(0), frac = numeric(0)))
  }
  x <- sort(x, decreasing = TRUE, na.last = TRUE)
  if (!is.null(k) && length(x) > k) x <- x[seq_len(k)]
  cf <- cumsum(x) / sum(x)
  data.frame(rank = seq_along(cf), frac = cf)
}

TOP_RANK <- max(K_TOP, 10000)  # display up to this many ranks (at least K_TOP)

curves_list <- lapply(colnames(count_mat), function(sid) {
  df <- cum_curve(count_mat[, sid], k = TOP_RANK)
  if (nrow(df) == 0) return(NULL)
  df$sample_id <- sid
  df
})
curves <- dplyr::bind_rows(curves_list)

# Compute frac_topK per sample (fraction of reads within top-K genes)
complexity_metric <- curves %>%
  group_by(sample_id) %>%
  summarize(
    frac_topK = {
      rmax <- max(rank, na.rm = TRUE)
      k    <- min(K_TOP, rmax)
      approx(x = rank, y = frac, xout = k, rule = 2)$y
    },
    .groups = "drop"
  ) %>%
  mutate(
    fail_complexity = frac_topK > COMPLEXITY_MAX,
    passComplexity  = !fail_complexity
  )

# Save complexity QC table
readr::write_tsv(complexity_metric, QC_TSV)
message("Complexity QC written to: ", QC_TSV)

# ================================
# Plot curves colored by complexity pass/fail
# ================================
curves_plot <- curves %>%
  dplyr::left_join(dplyr::select(complexity_metric, sample_id, passComplexity),
                   by = "sample_id") %>%
  dplyr::mutate(color_group = ifelse(passComplexity, "Pass (complexity)", "Fail (complexity)"))

idx_unique <- !duplicated(curves_plot$sample_id)
n_pass <- sum(curves_plot$passComplexity[idx_unique], na.rm = TRUE)
n_fail <- sum(!curves_plot$passComplexity[idx_unique], na.rm = TRUE)

p <- ggplot(curves_plot, aes(x = rank, y = frac, group = sample_id, color = color_group)) +
  geom_line(linewidth = 0.45) +
  scale_color_manual(values = c("Pass (complexity)" = "steelblue",
                                "Fail (complexity)" = "firebrick"),
                     name = NULL) +
  labs(
    title    = "RNA Library Complexity (complexity-only QC, IRB-excluded)",
    subtitle = paste0("Top-K = ", K_TOP,
                      ", Fail if fraction > ", COMPLEXITY_MAX,
                      " | Pass = ", n_pass, ", Fail = ", n_fail),
    x = "Ranked genes (descending by count)",
    y = "Fraction of all counts (cumulative)"
  ) +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "top")

# Save plot
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(PLOT_PNG, width = 1200, height = 700, res = 150)
  print(p); dev.off()
} else {
  ggsave(PLOT_PNG, p, width = 10, height = 6, dpi = 150)
}
message("Complexity plot saved to: ", PLOT_PNG)



# ================================
# Build counts_pass (complexity-only) and save
# ================================
pass_ids <- complexity_metric$sample_id[complexity_metric$passComplexity]
pass_ids <- intersect(pass_ids, colnames(count_mat))

counts_pass <- count_mat[, pass_ids, drop = FALSE]

# (Optional) filter extremely low-signal genes to reduce downstream noise
keep_genes <- rowSums(counts_pass) >= 10
counts_pass <- counts_pass[keep_genes, , drop = FALSE]

# Write counts_pass
out_tbl <- tibble(geneID = rownames(counts_pass)) %>%
  bind_cols(as_tibble(counts_pass, .name_repair = "minimal"))
readr::write_tsv(out_tbl, PASS_COUNTS_TSV)
message("Complexity-pass counts written to: ", PASS_COUNTS_TSV)

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(DESeq2)
  library(AnnotationDbi)
})

# t <- read_tsv("/workdir/qy236/apol_kidney/QC_output/DESeq2_significant_genes_ABC_vs_DEF_padj0.05_LFC0.5.tsv")
# View(t)
# ---------- inputs ----------
grp_csv <- "/workdir/qy236/apol_kidney/QC_output/group_sample_with_met_date.csv"
out_dir <- "/workdir/qy236/apol_kidney/QC_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 0) checks ----------
if (!exists("counts_pass")) stop("counts_pass not found. Create it first (genes x samples).")
count_matrix <- as.matrix(counts_pass)
if (is.null(rownames(count_matrix))) stop("counts_pass must have Ensembl IDs as rownames.")


# ---------- 1) map Ensembl -> SYMBOL (collapse duplicates by sum) ----------
# strip Ensembl version if present (ENSG... .xx)
ens_core <- sub("\\..*$", "", rownames(count_matrix))

# Ensure org.Hs.eg.db is available
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Package 'org.Hs.eg.db' is required. Please install it:\nBiocManager::install('org.Hs.eg.db')")
}
suppressPackageStartupMessages(library(org.Hs.eg.db))

# Map to gene symbols
sym_vec <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                 keys = ens_core,
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL",
                                 multiVals = "first")

# Keep rows with a mapped symbol
keep_idx <- !is.na(sym_vec) & nzchar(sym_vec)
if (!any(keep_idx)) stop("No Ensembl IDs mapped to gene symbols.")

# Collapse duplicates by SYMBOL using rowsum
mat_mapped <- count_matrix[keep_idx, , drop = FALSE]
sym_used   <- sym_vec[keep_idx]
counts_by_symbol <- rowsum(mat_mapped, group = sym_used)  # rows named by SYMBOL

# Save raw counts-by-symbol
out_counts_sym <- file.path(out_dir, "APOL_urine_feature_counts_passQC_bySymbol.tsv")
write_tsv(
  cbind.data.frame(gene=rownames(counts_by_symbol),
                   as.data.frame(counts_by_symbol),
                   row.names = NULL),
  out_counts_sym
)

# ---------- 2) read date mapping and align to columns ----------
grp <- readr::read_csv(grp_csv, show_col_types = FALSE)   # expects columns: group, sample_id, met_id, date

cols_full <- colnames(counts_by_symbol)                    # e.g., cfrna_APOL_89_bp150x150
cols_core <- sub("_bp.*$", "", cols_full)                  # -> cfrna_APOL_89

m_date <- match(cols_core, grp$met_id)
date_vec <- grp$date[m_date]
if (any(is.na(m_date))) {
  warning("Date missing for ", sum(is.na(m_date)), " samples: ",
          paste(cols_full[is.na(m_date)], collapse = ", "))
  date_vec[is.na(date_vec)] <- "Unknown"
}
date_fac <- factor(date_vec)  # 4 date categories expected

# ---------- 3) DESeq2 VST normalization on counts_by_symbol ----------
storage.mode(counts_by_symbol) <- "integer"

meta <- data.frame(sample_id = cols_full,
                   date = date_fac,
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

# Optional gene filter to reduce noise (keep genes with some counts)
keep_genes <- rowSums(counts_by_symbol) >= 10
counts_for_vst <- counts_by_symbol[keep_genes, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts_for_vst,
                              colData   = meta,
                              design    = ~ 1)
vsd <- vst(dds, blind = TRUE)    # size-factor normalization + variance-stabilizing transform
vst_mat <- assay(vsd)            # rows = SYMBOL, cols = samples

# Save VST matrix by symbol (optional)
out_vst_sym <- file.path(out_dir, "APOL_urine_vst_passQC_bySymbol.tsv")
write_tsv(
  cbind.data.frame(gene=rownames(vst_mat),
                   as.data.frame(vst_mat),
                   row.names = NULL),
  out_vst_sym
)

# ---------- 4) Pull APOL1 and build plot data ----------
if (!("APOL1" %in% rownames(vst_mat))) {
  stop("APOL1 not found in VST matrix. Check mapping or gene symbol availability.")
}
apol1_expr <- as.numeric(vst_mat["APOL1", ])
plot_df <- data.frame(
  sample_id = colnames(vst_mat),
  date      = meta$date,
  APOL1_vst = apol1_expr,
  stringsAsFactors = FALSE
)




# ---- enforce desired date order and re-plot ----
desired_levels <- c("one_month","three_month","six_month","twelve_month")

# coerce to plain character first, then set an ordered factor
plot_df$date <- as.character(plot_df$date)
plot_df$date <- factor(plot_df$date, levels = desired_levels, ordered = TRUE)

# (optional) quick sanity check on any unexpected labels
unexpected <- setdiff(unique(plot_df$date), desired_levels)
if (length(unexpected) > 0) {
  warning("Some samples have date labels not in desired_levels and will appear as NA: ",
          paste(setdiff(unique(as.character(plot_df$date)), desired_levels), collapse=", "))
}

p <- ggplot(plot_df, aes(x = date, y = APOL1_vst)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 1.8) +
  xlab("Date category") +
  ylab("APOL1 expression (VST)") +
  ggtitle("APOL1 VST expression across date categories") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

print(p)
ggsave(file.path(out_dir, "APOL1_VST_by_date.png"),
       p, width = 8, height = 6, dpi = 400)














suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(AnnotationDbi)
})

# ---- inputs ----
grp_csv <- "/workdir/qy236/apol_kidney/QC_output/group_sample_with_met_date.csv"
out_dir <- "/workdir/qy236/apol_kidney/QC_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# counts_pass must exist (genes x samples, rownames = Ensembl IDs, cols like cfrna_APOL_89_bp150x150)
if (!exists("counts_pass")) stop("counts_pass not found in environment.")
count_matrix <- as.matrix(counts_pass)
storage.mode(count_matrix) <- "integer"

# ---- Map Ensembl -> SYMBOL ----
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Please install org.Hs.eg.db: BiocManager::install('org.Hs.eg.db')")
}
suppressPackageStartupMessages(library(org.Hs.eg.db))

ens_core <- sub("\\..*$", "", rownames(count_matrix))
sym_vec <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                 keys     = ens_core,
                                 keytype  = "ENSEMBL",
                                 column   = "SYMBOL",
                                 multiVals = "first")

keep_idx <- !is.na(sym_vec) & nzchar(sym_vec)
mat_mapped <- count_matrix[keep_idx, , drop = FALSE]
sym_used   <- sym_vec[keep_idx]
counts_sym <- rowsum(mat_mapped, group = sym_used)

# ---- Compute log2 CPM ----
lib_sizes <- colSums(counts_sym)
lib_sizes[lib_sizes == 0] <- NA_real_

cpm_mat <- sweep(counts_sym, 2, lib_sizes, "/") * 1e6
log2cpm_mat <- log2(cpm_mat + 1)

# ---- Extract APOL1 ----
if (!("APOL1" %in% rownames(log2cpm_mat))) {
  stop("APOL1 not found in the log2CPM matrix.")
}
apol1_log2cpm <- as.numeric(log2cpm_mat["APOL1", ])
sample_ids <- colnames(log2cpm_mat)

# ---- Read date mapping ----
grp <- readr::read_csv(grp_csv, show_col_types = FALSE)   # columns: group, sample_id, met_id, date
cols_core <- sub("_bp.*$", "", sample_ids)
m <- match(cols_core, grp$met_id)
date_vec <- grp$date[m]
date_vec[is.na(date_vec)] <- "Unknown"

# ---- Build plot data ----
plot_df <- data.frame(
  sample_id   = sample_ids,
  date        = as.character(date_vec),
  APOL1_log2CPM = apol1_log2cpm,
  stringsAsFactors = FALSE
)

# enforce desired order
desired_levels <- c("one_month","three_month","six_month","twelve_month")
plot_df$date <- factor(plot_df$date, levels = desired_levels, ordered = TRUE)

# ---- Plot ----
p <- ggplot(plot_df, aes(x = date, y = APOL1_log2CPM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 1.8) +
  xlab("Date category") +
  ylab("APOL1 expression (log2 CPM)") +
  ggtitle("APOL1 log2 CPM across date categories") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

print(p)
ggsave(file.path(out_dir, "APOL1_log2CPM_by_date.png"),
       p, width = 8, height = 6, dpi = 400)






suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(ggplot2)
  library(AnnotationDbi)
})

## ========== INPUTS ==========
grp_csv <- "/workdir/qy236/apol_kidney/QC_output/group_sample_with_met.csv"  # cols: group, sample_id, met_id
out_dir <- "/workdir/qy236/apol_kidney/QC_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========== 1) COUNTS ==========
if (!exists("counts_pass")) stop("counts_pass not found (genes x samples).")
count_matrix <- as.matrix(counts_pass)
storage.mode(count_matrix) <- "integer"

## ========== 2) SAMPLE GROUPS: ABC vs DEF ==========
grp <- readr::read_csv(grp_csv, show_col_types = FALSE)  # has group + met_id
cols_full <- colnames(count_matrix)
cols_core <- sub("_bp.*$", "", cols_full)

m <- match(cols_core, grp$met_id)
group_orig <- grp$group[m]
grp_upper <- toupper(as.character(group_orig))
supergroup <- ifelse(grp_upper %in% c("A","B","C"), "ABC",
                     ifelse(grp_upper %in% c("D","E","F"), "DEF", NA))

keep <- !is.na(supergroup)
if (!all(keep)) {
  warning("Dropping samples not in A–F (cannot map to ABC/DEF): ",
          paste(cols_full[!keep], collapse = ", "))
}
count_matrix <- count_matrix[, keep, drop = FALSE]
supergroup <- supergroup[keep]
cols_full <- cols_full[keep]

## ========== 3) DESeq2 DATASET ==========
meta <- data.frame(sample_id = cols_full,
                   supergroup = factor(supergroup, levels = c("DEF","ABC")),
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData   = meta,
                              design    = ~ supergroup)

## ========== 4) RUN DESeq2 ==========
dds <- DESeq(dds)

## ========== 5) RESULTS (ABC vs DEF) + optional shrinkage ==========
res <- results(dds, contrast = c("supergroup","ABC","DEF"))
# Try apeglm shrinkage if available
if (requireNamespace("apeglm", quietly = TRUE)) {
  res <- lfcShrink(dds, coef = "supergroup_ABC_vs_DEF", type = "apeglm")
} else {
  message("Package 'apeglm' not installed; using unshrunk log2FC from results().")
}

res_df <- as.data.frame(res)
res_df$ensembl <- rownames(res_df)

## ========== 6) MAP ENSEMBL → SYMBOL ==========
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Please install org.Hs.eg.db:\n  BiocManager::install('org.Hs.eg.db')")
}
suppressPackageStartupMessages(library(org.Hs.eg.db))

ens_core <- sub("\\..*$", "", res_df$ensembl)  # strip version
sym_vec <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                 keys     = ens_core,
                                 keytype  = "ENSEMBL",
                                 column   = "SYMBOL",
                                 multiVals = "first")
res_df$symbol <- sym_vec

# For plotting labels, prefer symbol; fallback to Ensembl if missing
res_df$label <- ifelse(!is.na(res_df$symbol) & nzchar(res_df$symbol),
                       res_df$symbol, res_df$ensembl)

## ========== 7) FILTER: padj < 0.05 & |LFC| > 0.5 ==========
# thresholds
padj_cut <- 0.05
lfc_cut  <- 1

# add derived cols on the full results DF
res_df$negLog10Padj <- -log10(res_df$padj)
res_df$status <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cut & abs(res_df$log2FoldChange) > lfc_cut,
                        "Significant", "Not significant")

# pick up to top 20 significant genes to label
sig_idx  <- !is.na(res_df$padj) & res_df$padj < padj_cut & abs(res_df$log2FoldChange) > lfc_cut
label_df <- res_df[sig_idx, ]
label_df <- label_df[order(label_df$padj), ]
#if (nrow(label_df) > 20) label_df <- label_df[seq_len(20), ]

# make volcano
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  warning("Install ggrepel for nicer labels: install.packages('ggrepel')")
  use_repel <- FALSE
} else {
  use_repel <- TRUE
  library(ggrepel)
}

p_volc <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = status), alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "firebrick")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4) +
  xlab("log2 Fold Change (ABC vs DEF)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano: ABC vs DEF (padj < 0.05 & |LFC| > 0.5)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

if (nrow(label_df) > 0) {
  if (use_repel) {
    p_volc <- p_volc +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = log2FoldChange, y = -log10(padj), label = label),
        max.overlaps = 50,
        size = 3,
        box.padding = 0.35,
        point.padding = 0.25,
        min.segment.length = 0.1,
        inherit.aes = FALSE
      )
  } else {
    p_volc <- p_volc +
      geom_text(
        data = label_df,
        aes(x = log2FoldChange, y = -log10(padj), label = label),
        size = 3, vjust = -0.3,
        inherit.aes = FALSE
      )
  }
}

print(p_volc)
ggsave(file.path(out_dir, "Volcano_ABC_vs_DEF_padj0.05_LFC0.5.png"),
       p_volc, width = 10, height = 15, dpi = 400)


res <- read_tsv("/workdir/qy236/apol_kidney/QC_output/DESeq2_significant_genes_ABC_vs_DEF_padj0.05_LFC0.5.tsv")
View(res)





suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# thresholds
padj_cut <- 0.05
lfc_cut  <- 1

# add extra columns
res_df$negLog10Padj <- -log10(res_df$padj)

res_df$regulation <- "Not significant"
res_df$regulation[!is.na(res_df$padj) & res_df$padj < padj_cut & res_df$log2FoldChange > lfc_cut] <- "Up in ABC"
res_df$regulation[!is.na(res_df$padj) & res_df$padj < padj_cut & res_df$log2FoldChange < -lfc_cut] <- "Down in ABC"

# pick top genes to label
label_df <- res_df[res_df$regulation != "Not significant", ]
label_df <- label_df[order(label_df$padj), ]
if (nrow(label_df) > 20) label_df <- label_df[1:20, ]  # top 20

# volcano plot
p_volc <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = regulation), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up in ABC" = "firebrick",
                                "Down in ABC" = "steelblue",
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4) +
  xlab("log2 Fold Change (ABC vs DEF)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot: ABC vs DEF") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

if (nrow(label_df) > 0) {
  p_volc <- p_volc +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = log2FoldChange, y = -log10(padj), label = label, color = regulation),
      inherit.aes = FALSE,
      size = 3,
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = 50
    )
}

print(p_volc)

# save
volc_png <- file.path(out_dir, "Volcano_ABC_vs_DEF_redUp_blueDown.png")
ggsave(volc_png, p_volc, width = 10, height = 15, dpi = 400)
message("Volcano plot saved to: ", volc_png)










suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(ggplot2)
  library(AnnotationDbi)
})

## ---- thresholds ----
padj_cut <- 0.05
lfc_cut  <- 1.0   # you said you changed it to 1

## ---- 1) Identify significant genes by Ensembl ID ----
# If you used shrinkage --> res_df already contains final LFC/padj, else build from results(dds)
if (!exists("res_df")) {
  res <- results(dds, contrast = c("supergroup","ABC","DEF"))
  res_df <- as.data.frame(res)
  res_df$ensembl <- rownames(res_df)
}

sig_idx <- !is.na(res_df$padj) & (res_df$padj < padj_cut) & (abs(res_df$log2FoldChange) > lfc_cut)
sig_ens <- res_df$ensembl[sig_idx]

if (length(sig_ens) == 0) {
  stop("No significant genes at padj < 0.05 & |log2FC| > 1. Try lowering lfc_cut or checking your design.")
}

## ---- 2) Make sure we have a VST matrix for the *same* samples used in DE ----
# Reuse 'dds' you used for DE (has supergroup), then compute VST
vsd_plot <- vst(dds, blind = TRUE)
vst_mat  <- assay(vsd_plot)              # rows=genes (Ensembl IDs), cols=samples

# Keep only significant genes present in VST
sig_ens <- intersect(sig_ens, rownames(vst_mat))
if (length(sig_ens) == 0) stop("Significant genes not found in VST matrix.")

## ---- 3) Map Ensembl -> gene symbols for nicer facet labels (optional but recommended) ----
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Please install org.Hs.eg.db:\n  BiocManager::install('org.Hs.eg.db')")
}
suppressPackageStartupMessages(library(org.Hs.eg.db))

ens_core <- sub("\\..*$", "", sig_ens)
sym_vec  <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                  keys     = ens_core,
                                  keytype  = "ENSEMBL",
                                  column   = "SYMBOL",
                                  multiVals = "first")
# Build a lookup that preserves order/duplicates handling
sym_lookup <- setNames(sym_vec, sig_ens)
gene_labels <- ifelse(!is.na(sym_lookup[sig_ens]) & nzchar(sym_lookup[sig_ens]),
                      sym_lookup[sig_ens], sig_ens)

## ---- 4) Build long data.frame for plotting (expr + date + group) ----
# Expression for significant genes only
expr_sub <- vst_mat[sig_ens, , drop = FALSE]

# Sample metadata already has supergroup in dds/vsd colData
samps_full <- colnames(expr_sub)                  # e.g., cfrna_APOL_89_bp150x150
samps_core <- sub("_bp.*$", "", samps_full)

# Read date mapping (met_id + date)
grp_csv <- "/workdir/qy236/apol_kidney/QC_output/group_sample_with_met_date.csv"
grp_map <- readr::read_csv(grp_csv, show_col_types = FALSE)  # expects: group, sample_id, met_id, date

# Match to date via met_id
m_date   <- match(samps_core, grp_map$met_id)
date_vec <- grp_map$date[m_date]
if (any(is.na(date_vec))) {
  warning("Date missing for ", sum(is.na(date_vec)), " samples; labeling as 'Unknown'")
  date_vec[is.na(date_vec)] <- "Unknown"
}
# Use same group definition as DE (ABC vs DEF) using the dds colData
group_vec <- as.character(colData(vsd_plot)$supergroup)
if (is.null(group_vec)) {
  stop("supergroup not found in colData(vsd_plot). Make sure your 'dds' was built with 'supergroup' in colData.")
}

# Order timepoints
desired_levels <- c("one_month","three_month","six_month","twelve_month")
date_fac <- factor(as.character(date_vec), levels = desired_levels, ordered = TRUE)

# Melt to long format manually (no pipes)
# Rows = length(sig_ens) * nSamples
ng <- length(sig_ens)
ns <- ncol(expr_sub)
plot_df <- data.frame(
  gene    = rep(gene_labels, each = ns),
  gene_id = rep(sig_ens, each = ns),
  sample  = rep(samps_full, times = ng),
  date    = rep(date_fac, times = ng),
  group   = rep(group_vec, times = ng),
  expr    = as.numeric(t(expr_sub)),      # VST values
  stringsAsFactors = FALSE
)


## counts per gene × date × group
cnt <- aggregate(expr ~ gene + date + group, data = plot_df,
                 FUN = function(z) sum(!is.na(z)))
names(cnt)[names(cnt) == "expr"] <- "n"

# per-facet y for labels (headroom)
y_max <- tapply(plot_df$expr, plot_df$gene, function(v) max(v, na.rm = TRUE))
y_rng <- tapply(plot_df$expr, plot_df$gene, function(v) diff(range(v, na.rm = TRUE)))
pad   <- ifelse(is.na(y_rng) | y_rng == 0, 0.2, 0.08 * y_rng)
cnt$y <- y_max[as.character(cnt$gene)] + pad[as.character(cnt$gene)]

# same cnt / y label prep as you already computed above

dodge_w <- 0.6
pos_jit_dodge <- position_jitterdodge(jitter.width = 0.10, dodge.width = dodge_w)
pos_dodge     <- position_dodge(width = dodge_w)

p_sig <- ggplot(plot_df, aes(x = date, y = expr)) +
  # separate boxplots per group at each timepoint
  geom_boxplot(aes(fill = group),
               position = pos_dodge,
               width = 0.55,
               color = "grey40",
               outlier.shape = NA,
               alpha = 0.85) +
  # jittered points, dodged the same as boxplots
  geom_jitter(aes(color = group),
              size = 1.7, alpha = 0.85,
              position = pos_jit_dodge) +
  # "n=" labels per gene × date × group, dodged to sit above the correct box
  geom_text(data = cnt,
            aes(x = date, y = y, label = paste0("n=", n), color = group),
            position = pos_dodge,
            size = 3, vjust = -0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c(ABC = "mistyrose3",  # light version for boxes
                               DEF = "lightsteelblue3",
                               Other = "grey80"),
                    name = "Group") +
  scale_color_manual(values = c(ABC = "firebrick",
                                DEF = "steelblue",
                                Other = "grey60"),
                     name = "Group") +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  xlab("Timepoint") +
  ylab("Expression (VST)") +
  ggtitle("Significant genes (padj < 0.05 & |log2FC| > 1): VST expression by timepoint") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  ) +
  # extra headroom so the 'n=' labels don't clip
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))

print(p_sig)
ggsave(file.path(out_dir, "SigGenes_byTime_VST_groupSplitBoxes_withCounts.png"),
       p_sig, width = 21, height = 9, dpi = 400)



dodge_w <- 0.6
pos_jit_dodge <- position_jitterdodge(jitter.width = 0.10, dodge.width = dodge_w)
pos_dodge     <- position_dodge(width = dodge_w)

p_sig <- ggplot(plot_df, aes(x = date, y = expr)) +
  # separate violins per group at each timepoint
  geom_violin(aes(fill = group),
              position = pos_dodge,
              width = 0.75,
              alpha = 0.7,
              trim = FALSE) +
  # jittered points, dodged the same as violins
  geom_jitter(aes(color = group),
              size = 1.5, alpha = 0.8,
              position = pos_jit_dodge) +
  # add sample counts per gene × date × group
  geom_text(data = cnt,
            aes(x = date, y = y, label = paste0("n=", n), color = group),
            position = pos_dodge,
            size = 3, vjust = -0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c(ABC = "mistyrose3",
                               DEF = "lightsteelblue3",
                               Other = "grey80"),
                    name = "Group") +
  scale_color_manual(values = c(ABC = "firebrick",
                                DEF = "steelblue",
                                Other = "grey60"),
                     name = "Group") +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  xlab("Timepoint") +
  ylab("Expression (VST)") +
  ggtitle("Significant genes (padj < 0.05 & |log2FC| > 1): VST expression by timepoint") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))

print(p_sig)
ggsave(file.path(out_dir, "SigGenes_byTime_VST_violin_coloredByGroup.png"),
       p_sig, width = 21, height = 9, dpi = 400)


