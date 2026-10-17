#!/usr/bin/env Rscript
# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# =============================================================================
# ALS vs spaceflight manuscript | FIGURES 2, 3 AND SUPPLEMENTARY FIGURE S1
# Unified NASA + JAXA ACTUAL-SPACEFLIGHT results (not unloading analogues)
#
# Figure 2: nine prespecified metabolic gene sets across 16 primary comparisons
#            A: six exact SOD1G93A-reference KEGG sets
#            B: three broad Hallmark bioenergetic sets
# Figure 3: selected within- and between-mission genome-wide LFC correlations
#            A: RR-23 | B: BION-M1 | C: RR-5 | D: other comparisons
# Supplementary S1: MHU-1 3-group PCA, sample correlations, and
#            independent limma-voom/fgsea leave-one-GC-out sensitivity refits
#
# Input:  two CSV tables ALREADY produced by unified spaceflight pipeline.
# NO differential-expression or gene-set enrichment models are rerun.
# Primary comparisons only; secondary flight-vs-vivarium / MHU1 ground
# controls are not silently included in either main figure.
#
# RStudio: set working directory to project root, then
#   source("08_spaceflight_figures_2_3_and_s1.R")
# =============================================================================

# ---------------- USER SETTINGS ---------------------------------------------
PROJECT_DIR <- getwd()  # Change if script is run outside the main project folder
RESULTS_DIR <- file.path(PROJECT_DIR, "UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results")
FIGURE_DIR  <- file.path(PROJECT_DIR, "MANUSCRIPT_FIGURES_2_3_UNIFIED")
FIG2_WIDTH  <- 8.8    # compact PORTRAIT heatmap, readable at print width
FIG2_HEIGHT <- 8.8   # shorter y-axis: wider-than-tall heatmap tiles
FIG3_WIDTH  <- 13.5
FIG3_HEIGHT <- 8.8
DPI         <- 450
ALPHA       <- 0.05
# -----------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2: install.packages('ggplot2')")
}
if (!requireNamespace("svglite", quietly = TRUE)) {
  stop("Please install svglite: install.packages('svglite')")
}

dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# Locate source CSVs in the normal analysis folder, an extracted results/ folder,
# or results.zip. Unzip only the two required CSVs if necessary.
find_output <- function(rel_path) {
  probes <- c(
    file.path(RESULTS_DIR, rel_path),
    file.path(PROJECT_DIR, "results", rel_path),
    file.path(PROJECT_DIR, rel_path),
    file.path(PROJECT_DIR, basename(rel_path))
  )
  found <- probes[file.exists(probes)]
  if (length(found)) return(normalizePath(found[1], winslash = "/"))

  for (z in c(file.path(PROJECT_DIR, "results.zip"),
              file.path(RESULTS_DIR, "results.zip"))) {
    if (!file.exists(z)) next
    members <- utils::unzip(z, list = TRUE)$Name
    if (!rel_path %in% members) next
    dest <- file.path(FIGURE_DIR, "_inputs_from_results_zip")
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(z, files = rel_path, exdir = dest, overwrite = TRUE)
    unpacked <- file.path(dest, rel_path)
    if (file.exists(unpacked)) return(normalizePath(unpacked, winslash = "/"))
  }
  stop("Cannot find required unified output: ", rel_path,
       "\nExpected under: ", RESULTS_DIR,
       "\nAlternatively put results.zip in: ", PROJECT_DIR)
}

PATHWAYS_FILE <- find_output("GSEA/six_SOD1_KEGG_and_three_Hallmarks.csv")
CORRS_FILE    <- find_output("concordance/genomewide_cross_study_correlations.csv")
message("Pathway input: ", PATHWAYS_FILE)
message("Correlation input: ", CORRS_FILE)

p <- utils::read.csv(PATHWAYS_FILE, check.names = FALSE, stringsAsFactors = FALSE)
corr <- utils::read.csv(CORRS_FILE, check.names = FALSE, stringsAsFactors = FALSE)

require_columns <- function(dat, cols, source) {
  missing <- setdiff(cols, colnames(dat))
  if (length(missing)) stop("Missing columns in ", source, ": ",
                            paste(missing, collapse = ", "))
}
require_columns(p, c("dataset", "mission", "tissue", "cohort", "contrast",
                     "primary", "pathway", "NES", "padj"), "pathway file")
require_columns(corr, c("contrast1", "contrast2", "n_shared", "spearman_rho"),
                "correlation file")

# ============================================================================
# FIGURE 2: TWO PRESPECIFIED PATHWAY PANELS
# ============================================================================

is_primary <- tolower(trimws(as.character(p$primary))) %in%
  c("true", "t", "1", "yes")
p <- p[is_primary, , drop = FALSE]

pathway_ids <- c(
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "KEGG_PURINE_METABOLISM",
  "KEGG_STARCH_AND_SUCROSE_METABOLISM",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM"
)
# Short labels are intentional: nine pathways must fit across a PORTRAIT page.
# Full canonical names are retained above (pathway_ids) and in source-values CSV.
# Abbreviations are defined in the proposed figure legend / README.
pathway_labels <- c(
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS" = "Gly/GNG",
  "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM" = "Fru/Man",
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY" = "PPP",
  "KEGG_INSULIN_SIGNALING_PATHWAY" = "Insulin",
  "KEGG_PURINE_METABOLISM" = "Purine",
  "KEGG_STARCH_AND_SUCROSE_METABOLISM" = "Starch/Suc",
  "HALLMARK_GLYCOLYSIS" = "Glycolysis",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "OXPHOS",
  "HALLMARK_FATTY_ACID_METABOLISM" = "FA metab"
)

# Exact, manuscript-specified cohort registry; also acts as an input QC check.
# The two OSD-684 donor groups and both OSD-401 muscles are separate contrasts.
registry <- data.frame(
  dataset = c("OSD101", "OSD419", "OSD104",
              "OSD401", "OSD401",
              "OSD770", "OSD665", "OSD666",
              "OSD111", "OSD111", "OSD21",
              "OSD684", "OSD684", "MHU1", "MHU4", "MHU5"),
  tissue = c("GASTROCNEMIUS", "GASTROCNEMIUS", "SOLEUS",
             "GASTROCNEMIUS", "QUADRICEPS",
             "SOLEUS", "EDL", "QUADRICEPS",
             "SOLEUS", "EDL", "GASTROCNEMIUS",
             "MYOBUNDLE", "MYOBUNDLE", "SOLEUS", "SOLEUS", "SOLEUS"),
  cohort = c(rep("ALL", 11), "YOUNG", "OLD", "ALL", "ALL", "ALL"),
  contrast = c(rep("FLT_vs_HGC", 3), rep("FLT_vs_GC", 2),
               rep("FLT_vs_HGC", 5), "FLT_vs_GC", "FLT_vs_GC", "FLT_vs_GC",
               "MG_vs_AG", "PG_vs_GC", "PG_vs_GC"),
  plot_label = c(
    "RR-1 | Gastrocnemius L (OSD-101)",
    "RR-1 | Gastrocnemius R (OSD-419)",
    "RR-1 | Soleus (OSD-104)",
    "RR-5 | Gastrocnemius",
    "RR-5 | Quadriceps",
    "RR-23 | Soleus (OSD-770)",
    "RR-23 | EDL (OSD-665)",
    "RR-23 | Quadriceps (OSD-666)",
    "BION-M1 | Soleus",
    "BION-M1 | EDL",
    "STS-108 | Gastrocnemius",
    "CRS-21 | Young human myobundles",
    "CRS-21 | Older human myobundles",
    "MHU-1 | Soleus (MG vs AG)",
    "MHU-4 | Soleus (PG vs GC)",
    "MHU-5 | Soleus (PG vs GC)"
  ), stringsAsFactors = FALSE
)

make_key <- function(x) paste(x$dataset, x$tissue, x$cohort, x$contrast, sep="||")
registry$key <- make_key(registry)
p$key <- make_key(p)
observed <- unique(p$key)
missing <- setdiff(registry$key, observed)
extra <- setdiff(observed, registry$key)
if (length(missing) || length(extra)) {
  stop("The 16 primary comparisons differ from the expected unified analysis.\n",
       "Missing: ", paste(missing, collapse="; "), "\n",
       "Unexpected: ", paste(extra, collapse="; "),
       "\nInspect your primary contrast registry before plotting.")
}
if (length(observed) != 16L) stop("Expected exactly 16 primary contrasts.")
if (!setequal(unique(p$pathway), pathway_ids)) {
  stop("The pathway panel differs from the fixed six KEGG and three Hallmarks.")
}
if (anyDuplicated(paste(p$key,p$pathway))) {
  stop("Duplicate primary contrast/pathway combinations found.")
}
if (nrow(p) != 16L * 9L) {
  stop("Expected 144 pathway rows (16 contrasts x 9 pathways); got ", nrow(p))
}
if (any(!is.finite(p$NES)) || any(!is.finite(p$padj))) {
  stop("Non-finite NES or per-contrast pathway FDR; check GSEA output.")
}

p$plot_label <- registry$plot_label[match(p$key, registry$key)]
p$plot_label <- factor(p$plot_label, levels = rev(registry$plot_label))
p$pathway <- factor(p$pathway, levels = pathway_ids)
p$panel <- ifelse(grepl("^KEGG_",as.character(p$pathway)),
                  "A.  ALS-reference KEGG",
                  "B.  Hallmark")
p$panel <- factor(p$panel, levels = c(
  "A.  ALS-reference KEGG",
  "B.  Hallmark"))
p$asterisk <- ifelse(p$padj < ALPHA, "*", "")
p$asterisk_colour <- ifelse(abs(p$NES) >= 1.9, "white", "#20252B")

# Within-contrast FDR (padj) is from the full MSigDB gene-set testing universe;
# do NOT substitute the targeted 9-panel-wide FDR without also changing the text.
maximum <- max(3.5, ceiling(max(abs(p$NES)) * 2) / 2)
fig2 <- ggplot2::ggplot(p, ggplot2::aes(x = pathway, y = plot_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.65, height = 0.79) +
  ggplot2::geom_text(ggplot2::aes(label = asterisk,
                                 colour = asterisk_colour),
                     size = 5.1, fontface = "bold", show.legend = FALSE) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(low = "#C15B5B", mid = "#F5F5F5",
                                high = "#4A7EB3", midpoint = 0,
                                limits = c(-maximum, maximum), name = "NES") +
  ggplot2::facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  ggplot2::scale_x_discrete(labels = pathway_labels, drop = TRUE) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 15) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 13.5, face = "bold"),
    strip.background = ggplot2::element_rect(fill = "#EDEFF2", colour = NA),
    panel.spacing.x = grid::unit(0.45, "lines"),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 12.6, angle = 45,
                                       hjust = 1, vjust = 1, colour = "#282828"),
    axis.text.y = ggplot2::element_text(size = 12.5, colour = "#282828"),
    axis.ticks = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(size = 10.4, colour = "#4A4A4A", hjust = 0, lineheight = 1.03),
    legend.title = ggplot2::element_text(face = "bold", size = 11.5),
    legend.text = ggplot2::element_text(size = 11.5),
    legend.position = "right",
    plot.margin = ggplot2::margin(14, 12, 14, 12)
  )

# ============================================================================
# FIGURE 3: SELECTED GENOME-WIDE CORRELATIONS
# ============================================================================

# Selected comparisons follow EXACTLY the claims in the current Results text.
# Full set of 253 exploratory correlations is in the supplementary table.
fig3_specs <- data.frame(
  panel = c(rep("A  RR-23: muscle comparisons", 3),
            "B  BION-M1: muscle comparison",
            "C  RR-5: muscle comparison",
            rep("D  RR-1 and cross-mission comparisons", 3)),
  label = c("EDL vs quadriceps", "Soleus vs quadriceps", "Soleus vs EDL",
            "Soleus vs EDL", "Gastrocnemius vs quadriceps",
            "Left vs right gastrocnemius (RR-1)",
            "RR-1 vs STS-108 gastrocnemius",
            "MHU-4 vs MHU-5 soleus"),
  id1 = c(
    "OSD665__EDL__ALL__FLT_vs_HGC",
    "OSD770__SOLEUS__ALL__FLT_vs_HGC",
    "OSD770__SOLEUS__ALL__FLT_vs_HGC",
    "OSD111__SOLEUS__ALL__FLT_vs_HGC",
    "OSD401__GASTROCNEMIUS__ALL__FLT_vs_GC",
    "OSD101__GASTROCNEMIUS__ALL__FLT_vs_HGC",
    "OSD101__GASTROCNEMIUS__ALL__FLT_vs_HGC",
    "MHU4__SOLEUS__ALL__PG_vs_GC"),
  id2 = c(
    "OSD666__QUADRICEPS__ALL__FLT_vs_HGC",
    "OSD666__QUADRICEPS__ALL__FLT_vs_HGC",
    "OSD665__EDL__ALL__FLT_vs_HGC",
    "OSD111__EDL__ALL__FLT_vs_HGC",
    "OSD401__QUADRICEPS__ALL__FLT_vs_GC",
    "OSD419__GASTROCNEMIUS__ALL__FLT_vs_HGC",
    "OSD21__GASTROCNEMIUS__ALL__FLT_vs_GC",
    "MHU5__SOLEUS__ALL__PG_vs_GC"),
  stringsAsFactors = FALSE
)

get_pair <- function(id1, id2) {
  keep <- (corr$contrast1 == id1 & corr$contrast2 == id2) |
          (corr$contrast1 == id2 & corr$contrast2 == id1)
  z <- corr[keep,,drop = FALSE]
  if (nrow(z) != 1L) {
    stop("Expected exactly one cross-study pair:\n",id1,"\n",id2,
         "\nFound: ",nrow(z))
  }
  if (!is.finite(z$spearman_rho) || !is.finite(z$n_shared)) {
    stop("Missing rho or gene coverage for: ",id1," vs ",id2)
  }
  z
}
fig3_data <- do.call(rbind, lapply(seq_len(nrow(fig3_specs)), function(i) {
  x <- fig3_specs[i,,drop = FALSE]
  r <- get_pair(x$id1, x$id2)
  data.frame(panel=x$panel,label=x$label,
             spearman_rho=r$spearman_rho, n_shared=r$n_shared,
             stringsAsFactors=FALSE)
}))
fig3_data$panel <- factor(fig3_data$panel, levels=unique(fig3_specs$panel))
# Unique y factor levels for the repeated "Soleus vs EDL" names across panels.
fig3_data$plot_key <- paste0(seq_len(nrow(fig3_data))," ",fig3_data$label)
fig3_data$plot_key <- factor(fig3_data$plot_key,
                             levels=rev(fig3_data$plot_key))
fig3_data$summary <- sprintf("rho = %.3f    n = %s",
  fig3_data$spearman_rho,
  format(fig3_data$n_shared, big.mark = ",", scientific = FALSE, trim=TRUE))

fig3 <- ggplot2::ggplot(fig3_data, ggplot2::aes(y=plot_key)) +
  ggplot2::geom_vline(xintercept=0, linetype="dashed",
                      linewidth=0.4, colour="#8F969D") +
  ggplot2::geom_segment(ggplot2::aes(x=0,xend=spearman_rho,yend=plot_key),
                        colour="#B0BAC6", linewidth=1.2) +
  ggplot2::geom_point(ggplot2::aes(x=spearman_rho),
                      size=3.5,shape=21,fill="#357AA5", colour="#265878",
                      stroke=0.55) +
  ggplot2::geom_text(ggplot2::aes(x=spearman_rho,label=summary),
                     nudge_x=0.022,hjust=0,size=4.2,
                     colour="#283441") +
  ggplot2::facet_grid(panel ~ ., scales="free_y",space="free_y", switch="y") +
  ggplot2::scale_y_discrete(labels=function(x) sub("^[0-9]+ ","",x)) +
  ggplot2::scale_x_continuous(breaks=seq(0,0.8,0.2),
                             limits=c(-0.05,1.08), expand=c(0,0)) +
  ggplot2::labs(x="Spearman correlation of shared-gene log2 fold changes",
    y=NULL) +
  ggplot2::theme_minimal(base_size=13.5) +
  ggplot2::theme(
    strip.placement="outside",
    strip.background=ggplot2::element_rect(fill="#EDEFF2",colour=NA),
    strip.text.y.left=ggplot2::element_text(angle=0,face="bold",size=11.8),
    axis.text.y=ggplot2::element_text(size=12.0,colour="#283441"),
    axis.text.x=ggplot2::element_text(size=11.2),
    axis.title.x=ggplot2::element_text(size=13.2,face="bold",
                                       margin=ggplot2::margin(t=12)),
    panel.grid.major.y=ggplot2::element_blank(),
    panel.grid.minor=ggplot2::element_blank(),
    panel.spacing.y=grid::unit(1.0,"lines"),
    plot.caption=ggplot2::element_text(size=10.3,colour="#4A4A4A",hjust=0, lineheight = 1.03),
    plot.margin=ggplot2::margin(14,18,14,12)
  )

# ---------------- OUTPUT -----------------------------------------------------
write_fig <- function(plot, stem, width, height) {
  for (extension in c("png","svg")) {
    name <- file.path(FIGURE_DIR,paste0(stem,".",extension))
    if (extension == "png") {
      ggplot2::ggsave(name, plot=plot, width=width, height=height,
                      units="in",dpi=DPI,bg="white", limitsize=FALSE)
    } else {
      ggplot2::ggsave(name,plot=plot,width=width,height=height,
                      units="in", device=svglite::svglite,
                      bg="white", limitsize=FALSE)
    }
    message("Wrote: ",normalizePath(name,winslash="/"))
  }
}
write_fig(fig2,"Figure_2_UNIFIED_spaceflight_metabolic_pathways_PORTRAIT",
          FIG2_WIDTH,FIG2_HEIGHT)
write_fig(fig3,"Figure_3_UNIFIED_spaceflight_genomewide_concordance",
          FIG3_WIDTH,FIG3_HEIGHT)

utils::write.csv(p[,c("dataset","mission","tissue","cohort","contrast",
                      "pathway","NES","padj","plot_label")],
  file.path(FIGURE_DIR,"Figure_2_source_values.csv"),row.names=FALSE)
utils::write.csv(fig3_data[,c("panel","label","spearman_rho","n_shared")],
  file.path(FIGURE_DIR,"Figure_3_source_values.csv"),row.names=FALSE)
writeLines(c(
  "FIGURE 2: Compact portrait layout (9.5 x 10.2 inches), rectangular cells and enlarged manuscript fonts.",
  "Pathway abbreviations: Gly/GNG=glycolysis/gluconeogenesis; Fru/Man=fructose and mannose metabolism; PPP=pentose-phosphate pathway; Starch/Suc=starch and sucrose metabolism; OXPHOS=oxidative phosphorylation; FA metab=fatty-acid metabolism.",
  "FIGURE 2: ALS-reference and broader Hallmark GSEA; 16 primary actual-flight comparisons.",
  "Figure 2 pathway stars: FDR < 0.05 from full MSigDB testing universe within each contrast (CSV padj).",
  "Negative NES (warm red) and positive NES (blue) indicate down/up-enrichment.",
  "Figure 3: only manuscript-mentioned correlations. All 253 pairs remain in the original CSV.",
  "The Figure 3 correlations are computed between per-gene log2FC vectors, not expression profiles or individual mice.",
  "Figure 4 unloading pathway plot is UNCHANGED and comes from the separate unloading analysis.",
  paste0("Pathway data: ",PATHWAYS_FILE),
  paste0("Correlation data: ",CORRS_FILE)
),file.path(FIGURE_DIR,"FIGURE_PROVENANCE.txt"))
message("\nFigures 2 and 3 saved; now generating Supplementary Figure S1 from MHU-1 count data...")


# ==========================================================================
# SUPPLEMENTARY FIGURE S1: MHU-1 3-GROUP QC AND LEAVE-ONE-GC-OUT ROBUSTNESS
# ==========================================================================
# In contrast to the older DESeq2-based QC/sensitivity figure, all sensitivity
# comparisons here are independently refitted using the UNIFIED analysis model:
# fractional estimated counts -> edgeR filterByExpr -> TMM -> limma::voom ->
# robust limma eBayes -> mouse-symbol Wald-like moderated t-rank ->
# fgseaMultilevel across the FULL KEGG Legacy + Hallmark collection.
# S1A/B use all nine MHU-1 soleus animals and TMM-normalized logCPM.
# D000017 is NEVER excluded from the full primary fit, only in named sensitivity.
# MHU-1 MG vs AG is onboard-only and unaffected by GC leave-one-out.
# --------------------------------------------------------------------------

S1_FIG_WIDTH  <- 9.5
S1_FIG_HEIGHT <- 12.2
S1_SEED       <- 20261015L
set.seed(S1_SEED)

s1_pkg <- c("data.table", "edgeR", "limma", "fgsea", "msigdbr",
            "org.Mm.eg.db", "AnnotationDbi", "patchwork", "statmod",
            "BiocParallel")
s1_missing <- s1_pkg[!vapply(s1_pkg, requireNamespace, logical(1), quietly=TRUE)]
if (length(s1_missing)) stop(
  "Supplementary Figure S1 needs installed R packages: ",
  paste(s1_missing, collapse=", "),
  ".\nCRAN: install.packages(c('data.table','msigdbr','patchwork','statmod'))",
  "\nBioconductor: BiocManager::install(c('edgeR','limma','fgsea','org.Mm.eg.db','AnnotationDbi','BiocParallel'))"
)

locate_ibSLS <- function(stem) {
  # The deposited files were often named '*_for_MHU1(1).tsv' by Windows.
  nm <- c(paste0(stem,".tsv"),paste0(stem,"(1).tsv"))
  direct <- file.path(PROJECT_DIR,nm)
  existing <- direct[file.exists(direct)]
  if (length(existing)) return(normalizePath(existing[1],winslash="/"))
  possible <- list.files(PROJECT_DIR,pattern="\\.tsv$",recursive=TRUE,
                         full.names=TRUE,ignore.case=TRUE)
  possible <- possible[basename(possible) %in% nm]
  if (length(possible)!=1L) stop("Expected ONE MHU-1 file for ",stem,
    ". Found ",length(possible),". Put the expected file in PROJECT_DIR.")
  normalizePath(possible,winslash="/")
}

mhu1_count_file <- locate_ibSLS("20260331_COUNTtable_for_MHU1")
mhu1_meta_file  <- locate_ibSLS("20260331_Metadata_for_MHU1")
message("S1 count matrix: ",mhu1_count_file)
message("S1 metadata: ",mhu1_meta_file)

mhu1_meta <- data.table::fread(mhu1_meta_file,data.table=FALSE,
                               check.names=FALSE)
need_mhu1_meta <- c("ibSLSDataReleaseID","GroupID","TissueName")
require_columns(mhu1_meta,need_mhu1_meta,"MHU-1 metadata")
mhu1_meta <- mhu1_meta[
  grepl("SOLEUS",toupper(trimws(as.character(mhu1_meta$TissueName)))),
  ,drop=FALSE]
if (nrow(mhu1_meta)!=9L) stop("Expected 9 soleus samples in MHU-1 metadata; found ",
                             nrow(mhu1_meta))
mhu1_meta$sample <- as.character(mhu1_meta$ibSLSDataReleaseID)
mhu1_meta$group <- toupper(trimws(as.character(mhu1_meta$GroupID)))
if (!all(mhu1_meta$group %in% c("MHU01_MG","MHU01_AG","MHU01_GC")))
  stop("Unexpected/other-mission MHU-1 GroupID: ",
       paste(unique(mhu1_meta$group),collapse=", "))
mhu1_meta$group <- sub("^MHU01_","",mhu1_meta$group)
if (!all(table(factor(mhu1_meta$group,levels=c("MG","AG","GC")))==3L))
  stop("MHU-1 must have 3 biological samples/group in soleus. Group counts: ",
       paste(names(table(mhu1_meta$group)),table(mhu1_meta$group),collapse="; "))
if (!"D000017" %in% mhu1_meta$sample ||
    mhu1_meta$group[match("D000017",mhu1_meta$sample)] != "GC")
  stop("MHU-1 expected atypical GC sample D000017 is missing or misclassified")
if (anyDuplicated(mhu1_meta$sample)) stop("Duplicated ibSLS biological sample IDs")

# Read exactly the gene-level count table supplied by ibSLS; no inferred groups.
raw_mhu1 <- data.table::fread(mhu1_count_file,data.table=FALSE,check.names=FALSE)
if (!"GeneID" %in% names(raw_mhu1)) stop("MHU-1 count table lacks GeneID")
if (!all(mhu1_meta$sample %in% names(raw_mhu1))) {
  stop("MHU-1 soleus sample IDs missing from count matrix: ",
       paste(setdiff(mhu1_meta$sample,names(raw_mhu1)),collapse=", "))
}
mhu1_meta <- mhu1_meta[order(match(mhu1_meta$group,c("AG","MG","GC")),
                             mhu1_meta$sample),,drop=FALSE]
raw_ids <- sub("\\.[0-9]+$","",as.character(raw_mhu1$GeneID))
mhu1_counts <- as.matrix(raw_mhu1[,mhu1_meta$sample,drop=FALSE])
storage.mode(mhu1_counts) <- "double"
if (anyNA(raw_ids) || any(!nzchar(raw_ids)) ||
    any(!is.finite(mhu1_counts)) || any(mhu1_counts<0))
  stop("MHU-1 count input has invalid IDs, missing, negative or non-finite values")
if (anyDuplicated(raw_ids)) {
  mhu1_counts <- rowsum(mhu1_counts,group=raw_ids,reorder=FALSE)
} else rownames(mhu1_counts) <- raw_ids
stopifnot(identical(colnames(mhu1_counts),mhu1_meta$sample))

# S1A/B: normalise using the full 3-group, 9-sample specimen collection.
all_design <- stats::model.matrix(~0+factor(mhu1_meta$group,
                                 levels=c("GC","AG","MG")))
y_qc <- edgeR::DGEList(counts=mhu1_counts)
keep_qc <- edgeR::filterByExpr(y_qc,design=all_design)
if (sum(keep_qc)<1000L) stop("Not enough expressed MHU-1 genes for QC")
y_qc <- edgeR::calcNormFactors(y_qc[keep_qc,,keep.lib.sizes=FALSE],method="TMM")
logcpm <- edgeR::cpm(y_qc,log=TRUE,prior.count=2)
variances <- apply(logcpm,1,stats::var)
var_ids <- head(order(variances,decreasing=TRUE),500L)
pca <- stats::prcomp(t(logcpm[var_ids,,drop=FALSE]),center=TRUE,scale.=FALSE)
pca_points <- data.frame(sample=rownames(pca$x),
                          PC1=pca$x[,1],PC2=pca$x[,2],stringsAsFactors=FALSE)
pca_points$group <- mhu1_meta$group[match(pca_points$sample,mhu1_meta$sample)]
pca_points$atypical <- pca_points$sample=="D000017"
explained <- 100*(pca$sdev^2)/sum(pca$sdev^2)

qc_cols <- c(MG="#3679B0",AG="#459B7B",GC="#C68B43")
fig_s1a <- ggplot2::ggplot(pca_points,ggplot2::aes(PC1,PC2,colour=group))+
  ggplot2::geom_point(size=3.5)+
  ggplot2::geom_point(data=pca_points[pca_points$atypical,,drop=FALSE],
                      mapping=ggplot2::aes(PC1,PC2),
                      shape=21,fill=NA,colour="#A31D25",size=6.0,stroke=1.2,
                      inherit.aes=FALSE)+
  ggplot2::geom_text(data=pca_points[pca_points$atypical,,drop=FALSE],
                      ggplot2::aes(label=sample),hjust=-0.10,vjust=-0.75,
                      size=4.0,colour="#A31D25",nudge_x=-25,show.legend=FALSE)+
  ggplot2::scale_colour_manual(values=qc_cols,breaks=c("MG","AG","GC"))+
  ggplot2::labs(title="A  MHU-1 soleus: sample PCA",
       x=sprintf("PC1 (%.1f%%)",explained[1]),
       y=sprintf("PC2 (%.1f%%)",explained[2]),colour="Group")+
  ggplot2::theme_bw(base_size=13)+
  ggplot2::theme(legend.position="bottom",plot.title=ggplot2::element_text(face="bold",size=14),
       legend.text=ggplot2::element_text(size=11),panel.grid.minor=ggplot2::element_blank())+
  ggplot2::coord_cartesian(clip="off")

# Correlations across the same expressed genes and nine normalized specimens.
r_matrix <- stats::cor(logcpm,method="pearson")
corr_long <- as.data.frame(as.table(r_matrix),stringsAsFactors=FALSE)
names(corr_long) <- c("sample_x","sample_y","pearson_r")
axis_order <- mhu1_meta$sample
axis_labels <- setNames(paste0(mhu1_meta$group,"-",
                                sub("^D0000","",mhu1_meta$sample),
                                ifelse(mhu1_meta$sample=="D000017","*","")),axis_order)
corr_long$sample_x <- factor(as.character(corr_long$sample_x),levels=axis_order)
corr_long$sample_y <- factor(as.character(corr_long$sample_y),levels=rev(axis_order))
corr_long$atypical <- as.character(corr_long$sample_x)=="D000017" |
                      as.character(corr_long$sample_y)=="D000017"
fig_s1b <- ggplot2::ggplot(corr_long,ggplot2::aes(sample_x,sample_y,fill=pearson_r))+
  ggplot2::geom_tile(ggplot2::aes(colour=atypical),linewidth=0.9)+
  ggplot2::scale_colour_manual(values=c(`TRUE`="#A31D25",`FALSE`="white"),
                               guide="none")+
  ggplot2::scale_fill_gradient(low="#ECF4FA",high="#17476D",
               limits=c(min(corr_long$pearson_r)-0.005,1),name="Pearson r")+
  ggplot2::scale_x_discrete(labels=axis_labels,drop=FALSE)+
  ggplot2::scale_y_discrete(labels=axis_labels,drop=FALSE)+
  ggplot2::coord_fixed()+
  ggplot2::labs(title="B  MHU-1 soleus: sample correlations",x=NULL,y=NULL)+
  ggplot2::theme_minimal(base_size=13)+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=14),
       panel.grid=ggplot2::element_blank(),
       axis.text.x=ggplot2::element_text(angle=55,hjust=1,size=9.5),
       axis.text.y=ggplot2::element_text(size=9.5),
       legend.position="bottom",legend.key.width=grid::unit(1.15,"cm"),
       legend.text=ggplot2::element_text(size=10))

# Full MSigDB mouse gene sets: identical registry to UNIFIED pipeline.
make_mhu1_sets <- function() {
  k <- msigdbr::msigdbr(species="Mus musculus",collection="C2",
                        subcollection="CP:KEGG_LEGACY")
  h <- msigdbr::msigdbr(species="Mus musculus",collection="H")
  z <- split(c(as.character(k$gene_symbol),as.character(h$gene_symbol)),
             c(as.character(k$gs_name),as.character(h$gs_name)))
  lapply(z,unique)
}
all_mhu1_sets <- make_mhu1_sets()
if (!all(pathway_ids %in% names(all_mhu1_sets)))
  stop("S1: missing original KEGG or Hallmark gene sets: ",
       paste(setdiff(pathway_ids,names(all_mhu1_sets)),collapse=", "))

# Match Ensembl mouse genes to official gene symbols using the same mapping.
valid_ids <- grepl("^ENSMUSG[0-9]+$",rownames(mhu1_counts))
if (sum(valid_ids)<1000L) stop("MHU-1 gene identifiers unexpectedly non-Ensembl")
syms <- rep(NA_character_,nrow(mhu1_counts))
mi <- rownames(mhu1_counts)[valid_ids]
known <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
       keys=unique(mi),keytype="ENSEMBL",column="SYMBOL",multiVals="first")
syms[valid_ids] <- unname(known[mi])
syms[!valid_ids] <- rownames(mhu1_counts)[!valid_ids]
ann_mhu1 <- data.frame(gene_id=rownames(mhu1_counts),symbol=syms,
                        stringsAsFactors=FALSE)

# Refit every two-group contrast separately, matching the unified fit_contrast.
fit_mhu1_sensitivity <- function(sample_ids,num,den,run_name) {
  meta_fit <- mhu1_meta[match(sample_ids,mhu1_meta$sample),,drop=FALSE]
  if (anyNA(meta_fit$sample) || any(!meta_fit$group %in% c(num,den)))
    stop("S1 invalid sample subset for ",run_name)
  tab <- table(factor(meta_fit$group,levels=c(den,num)))
  if (any(tab<2)) stop("S1 contrast requires at least two biological mice/group: ",run_name)
  xx <- mhu1_counts[,sample_ids,drop=FALSE]
  d <- stats::model.matrix(~0+factor(meta_fit$group,levels=c(den,num)))
  colnames(d) <- c("CONTROL","EXPOSURE")
  yy <- edgeR::DGEList(counts=xx)
  keep <- edgeR::filterByExpr(yy,design=d)
  if (sum(keep)<1000L) stop("Too few expressed genes in ",run_name)
  yy <- edgeR::calcNormFactors(yy[keep,,keep.lib.sizes=FALSE],method="TMM")
  vv <- limma::voom(yy,d,plot=FALSE)
  cm <- limma::makeContrasts(EXPOSURE-CONTROL,levels=d)
  lf <- limma::lmFit(vv,d)
  lf <- limma::eBayes(limma::contrasts.fit(lf,cm),robust=TRUE)
  tt <- limma::topTable(lf,coef=1,sort.by="none",number=Inf)
  tt$gene_id <- rownames(tt)
  tt$symbol <- ann_mhu1$symbol[match(tt$gene_id,ann_mhu1$gene_id)]
  # Deduplicate by average expression, not by the observed effect or significance.
  tt <- tt[!is.na(tt$symbol) & nzchar(tt$symbol) & is.finite(tt$t),,drop=FALSE]
  tt <- tt[order(-tt$AveExpr,tt$gene_id),,drop=FALSE]
  tt <- tt[!duplicated(tt$symbol),,drop=FALSE]
  if (nrow(tt)<1000L) stop("S1: too few annotated mouse genes for ",run_name)
  ranks <- stats::setNames(tt$t,tt$symbol)
  ranks <- sort(ranks,decreasing=TRUE)
  paths <- all_mhu1_sets[lengths(lapply(all_mhu1_sets,intersect,names(ranks)))>=10L]
  # IMPORTANT (Windows): never inherit a multicore/SOCK BiocParallel backend.
  # Explicit SerialParam avoids serialization errors such as
  # "Error in serialize(data, node$con): error writing to connection".
  # This changes execution mode ONLY; gene rankings, contrasts and FDR stay the same.
  ee <- as.data.frame(fgsea::fgseaMultilevel(
    pathways = paths, stats = ranks,
    minSize = 10L, maxSize = 500L, eps = 0,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  ee <- ee[ee$pathway %in% pathway_ids,c("pathway","NES","pval","padj"),drop=FALSE]
  if (nrow(ee)!=length(pathway_ids)) stop("S1 GSEA: missing selected pathways for ",run_name)
  ee$comparison <- paste(num,"vs",den)
  ee$contrast <- paste(num,"vs",den,sep="_")
  ee$run <- run_name
  ee$n_numerator <- unname(tab[num]);ee$n_denominator <- unname(tab[den])
  ee
}

# The complete 9-animal model uses only two relevant groups per contrast,
# exactly as in the unified analysis. Omission fits exclude the named GC mouse.
comparison_specs <- list(c("MG","GC"),c("AG","GC"),c("MG","AG"))
fit_rows <- list()
for (contrast in comparison_specs) {
  num <- contrast[1]; den <- contrast[2]
  full_ids <- mhu1_meta$sample[mhu1_meta$group %in% c(num,den)]
  full_name <- paste(num,"vs",den)
  message(sprintf("S1: fitting %s vs %s — full dataset", num, den))
  fit_rows[[length(fit_rows)+1L]] <- fit_mhu1_sensitivity(full_ids,num,den,"Full")
  if (den=="GC") {
    for (omit in mhu1_meta$sample[mhu1_meta$group=="GC"]) {
      ids <- setdiff(full_ids,omit)
      message(sprintf("S1: fitting %s vs %s — omit %s", num, den, omit))
      fit_rows[[length(fit_rows)+1L]] <- fit_mhu1_sensitivity(
        ids,num,den,paste0("Exclude ",omit))
    }
  }
}
s1_enrichment <- do.call(rbind,fit_rows)
s1_enrichment <- s1_enrichment[
  order(match(s1_enrichment$comparison,c("MG vs GC","AG vs GC","MG vs AG")),
        s1_enrichment$pathway,s1_enrichment$run),,drop=FALSE]
utils::write.csv(s1_enrichment,
  file.path(FIGURE_DIR,"Supplementary_Figure_S1_MHU1_sensitivity_GSEA.csv"),
  row.names=FALSE)
utils::write.csv(mhu1_meta[,c("sample","group","TissueName","ibSLSSampleID")],
  file.path(FIGURE_DIR,"Supplementary_Figure_S1_MHU1_sample_registry.csv"),
  row.names=FALSE)
utils::write.csv(pca_points,
  file.path(FIGURE_DIR,"Supplementary_Figure_S1_MHU1_PCA_coordinates.csv"),
  row.names=FALSE)
utils::write.csv(corr_long,
  file.path(FIGURE_DIR,"Supplementary_Figure_S1_MHU1_sample_correlations.csv"),
  row.names=FALSE)

# Optional reproducibility audit of the full, non-omission contrast versus the
# existing UNIFIED CSV. Warn rather than silently claiming identical pipelines.
p_all <- utils::read.csv(PATHWAYS_FILE,check.names=FALSE,stringsAsFactors=FALSE)
truth <- p_all[p_all$dataset=="MHU1" & p_all$tissue=="SOLEUS" &
               p_all$cohort=="ALL" & p_all$pathway %in% pathway_ids,
               c("contrast","pathway","NES","padj"),drop=FALSE]
full_s1 <- s1_enrichment[s1_enrichment$run=="Full",
                     c("contrast","pathway","NES","padj"),drop=FALSE]
audit <- merge(full_s1,truth,by=c("contrast","pathway"),
               suffixes=c("_S1_refit","_unified"))
audit$NES_difference <- audit$NES_S1_refit-audit$NES_unified
utils::write.csv(audit,
  file.path(FIGURE_DIR,"Supplementary_Figure_S1_unified_NES_reproducibility_audit.csv"),
  row.names=FALSE)
if (nrow(audit)!=3L*length(pathway_ids) ||
    any(abs(audit$NES_difference)>0.05,na.rm=TRUE)) {
  warning("S1 refitted full-model enrichment differs from the saved unified ",
          "results (or some MHU1 comparisons are absent). Check ",
          "Supplementary_Figure_S1_unified_NES_reproducibility_audit.csv ",
          "before including S1 in the manuscript.")
}

s1_enrichment$comparison <- factor(s1_enrichment$comparison,
                                    levels=c("MG vs GC","AG vs GC","MG vs AG"))
s1_enrichment$pathway <- factor(s1_enrichment$pathway,
                                levels=rev(pathway_ids))
# Separate full-data values from the omit-one-control estimates for clear legend.
base_values <- s1_enrichment[s1_enrichment$run=="Full",,drop=FALSE]
omit_values <- s1_enrichment[s1_enrichment$run!="Full",,drop=FALSE]
range_values <- stats::aggregate(NES~comparison+pathway,omit_values,
          function(x) c(lo=min(x),hi=max(x)))
range_values <- data.frame(comparison=range_values$comparison,
           pathway=range_values$pathway,
           xmin=range_values$NES[,"lo"],xmax=range_values$NES[,"hi"])
fig_s1c <- ggplot2::ggplot()+
  ggplot2::geom_vline(xintercept=0,linetype="dashed",colour="#909090",
                       linewidth=0.45)+
  ggplot2::geom_segment(data=range_values,
       ggplot2::aes(x=xmin,xend=xmax,y=pathway,yend=pathway),
       colour="#A6ADB5",linewidth=1.2)+
  ggplot2::geom_point(data=omit_values,
       ggplot2::aes(x=NES,y=pathway,colour=run),size=2.5,alpha=0.95)+
  ggplot2::geom_point(data=base_values,
       ggplot2::aes(x=NES,y=pathway,shape=run),size=3.1,
       fill="#20252B",colour="#20252B")+
  ggplot2::facet_wrap(~comparison,nrow=1,scales="fixed")+
  ggplot2::scale_y_discrete(labels=pathway_labels)+
  ggplot2::scale_colour_manual(values=c(`Exclude D000016`="#9299A3",
             `Exclude D000017`="#BC393D",`Exclude D000018`="#C1C7CD"),
             name="Sensitivity fit")+
  ggplot2::scale_shape_manual(values=c(18),name="Reference")+
  ggplot2::labs(title="C  Sensitivity of metabolic enrichment to GC sample removal",
    x="Normalised enrichment score (NES)",y=NULL)+
  ggplot2::theme_bw(base_size=13)+
  ggplot2::theme(
    plot.title=ggplot2::element_text(face="bold",size=14),
    axis.text.y=ggplot2::element_text(size=10.9),
    axis.text.x=ggplot2::element_text(size=10.5),
    axis.title.x=ggplot2::element_text(size=12.2),
    strip.text=ggplot2::element_text(face="bold",size=12),
    panel.grid.minor=ggplot2::element_blank(),
    legend.position="bottom",legend.text=ggplot2::element_text(size=10),
    legend.title=ggplot2::element_text(size=10.5))

# Combine 3 panels. SVG remains editable/vector; PNG is high resolution.
supp_fig1 <- (fig_s1a | fig_s1b) / fig_s1c +
  patchwork::plot_layout(heights=c(1.0,1.50))
write_fig(supp_fig1,"Supplementary_Figure_S1_MHU1_QC_and_sensitivity",
          S1_FIG_WIDTH,S1_FIG_HEIGHT)
write_fig(fig_s1a,"Supplementary_Figure_S1A_MHU1_PCA",7.5,5.5)
write_fig(fig_s1b,"Supplementary_Figure_S1B_MHU1_sample_correlations",7.5,6.2)
write_fig(fig_s1c,"Supplementary_Figure_S1C_MHU1_leave_one_GC_out",10.7,5.8)

cat(paste0("\nSUPPLEMENTARY FIGURE S1\n", 
  "Input counts: ",mhu1_count_file,"\n",
  "Input metadata: ",mhu1_meta_file,"\n",
  "QC: all 9 MHU-1 soleus biological samples; 3 per group\n",
  "D000017 is retained in the complete analysis and highlighted for inspection.\n",
  "Sensitivity: MG_vs_GC and AG_vs_GC independently refitted after omitting ",
  "D000016/D000017/D000018; MG_vs_AG onboard comparison does not involve GC.\n",
  "RNA-seq model: edgeR filterByExpr, TMM, voom, limma robust eBayes.\n",
  "Enrichment: moderated t rank in mouse symbols; MSigDB KEGG Legacy + ",
  "Hallmarks; fgseaMultilevel, 10..500 genes, eps=0.\n",
  "Supplementary sensitivity estimates are descriptive, not extra replicates.\n"),
  file=file.path(FIGURE_DIR,"Supplementary_Figure_S1_provenance.txt"))
message("\nCOMPLETE: Figures 2 and 3 and Supplementary Figure S1 are saved to: ",
        FIGURE_DIR)
