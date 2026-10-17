# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

###############################################################################
# FIGURE 11 — Denervation recapitulates the ALS-directional muscle signature
#
# Reads ONLY results already produced by 4_denervation_hypothesis_analysis.R.
# Does NOT rerun differential expression, CAMERA, fgsea, or t-tests.
#
# Source workbook: denervation_hypothesis_summary.xlsx
# Required sheets: Direction concordance, CAMERA, fgsea,
#                  Module score tests, Sample scores
#
# Panels
# A: Gene-level concordance at |log2FC| >= 0.10 (all three comparisons)
# B: fgsea NES plus CAMERA FDR (two fixed-direction module components)
# C: Every individual sample's signed module score (n=3 per group)
# D: Sample-level score differences and their actual 95% CIs
#
# Outputs: Figure_11_denervation.png and .svg (400 dpi PNG).
###############################################################################

# OPTIONAL: set this explicitly to the original analysis's results workbook.
# If left as NA, the script searches the two typical working-directory paths.
WORKBOOK <- NA_character_
OUTDIR <- NA_character_  # default: publication_figures next to WORKBOOK

CRAN_PACKAGES <- c("readxl", "dplyr", "tidyr", "ggplot2",
                   "patchwork", "ragg", "svglite")
missing_packages <- CRAN_PACKAGES[
  !vapply(CRAN_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork)
})

if (is.na(WORKBOOK)) {
  candidates <- c(
    file.path(getwd(), "denervation_hypothesis_analysis",
              "denervation_hypothesis_summary.xlsx"),
    file.path(getwd(), "denervation_hypothesis_summary.xlsx")
  )
  available <- candidates[file.exists(candidates)]
  if (!length(available)) {
    stop("Cannot locate denervation_hypothesis_summary.xlsx. ",
         "Set WORKBOOK to its full path at the beginning of this script.\n",
         "Paths checked:\n", paste(candidates, collapse = "\n"), call. = FALSE)
  }
  WORKBOOK <- available[[1L]]
}
stopifnot(file.exists(WORKBOOK))
if (is.na(OUTDIR)) OUTDIR <- file.path(dirname(WORKBOOK), "publication_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# CRITICAL: the summary workbook exported by openxlsx has an A1-only
# worksheet dimension recorded in XML, even though real data extend beyond A1.
# We ALWAYS supply explicit cell limits so readxl does not truncate the tables.
read_result <- function(sheet, ncols, rows = 1000L) {
  x <- as.data.frame(readxl::read_excel(
    WORKBOOK, sheet = sheet,
    range = readxl::cell_limits(c(1L, 1L), c(rows, ncols)),
    guess_max = rows,
    .name_repair = "minimal"
  ), check.names = FALSE)
  x <- x[rowSums(!is.na(x)) > 0L, , drop = FALSE]
  if (!nrow(x)) stop("Empty results sheet: ", sheet, call. = FALSE)
  x
}

need <- function(dat, cols, label) {
  absent <- setdiff(cols, names(dat))
  if (length(absent)) {
    stop(label, ": missing columns: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

concordance <- read_result("Direction concordance", 12L)
camera <- read_result("CAMERA", 11L)
fgsea <- read_result("fgsea", 14L)
score_tests <- read_result("Module score tests", 16L)
scores <- read_result("Sample scores", 4L)

need(concordance, c("dataset", "contrast", "module_component", "n_evaluable",
  "n_concordant_at_abs_log2FC_0_10", "n_discordant_at_abs_log2FC_0_10",
  "n_small_effect"), "Direction concordance")
need(camera, c("dataset", "contrast", "gene_set", "Direction",
               "FDR_all_tests"), "CAMERA")
need(fgsea, c("dataset", "contrast", "pathway", "NES"), "fgsea")
need(score_tests, c("dataset", "contrast", "numerator", "denominator",
                   "score_difference", "ci_low", "ci_high", "FDR_all_tests"),
     "Module score tests")
need(scores, c("dataset", "sample", "group", "module_score"), "Sample scores")

# Explicit contrast order. All data values are read from the workbook.
keys <- data.frame(
  dataset = c("GSE217576", "GSE117440", "GSE117440"),
  contrast = c("DEN_vs_CTL", "DEN_vs_SHAM",
               "Electrical_stimulation_vs_DEN"),
  comparison = c("Mouse DEN vs CTL", "Rat DEN vs SHAM",
                 "Rat STIM vs DEN"),
  stringsAsFactors = FALSE
)
pretty_order <- keys$comparison
as_comparison <- function(df) {
  df <- dplyr::inner_join(df, keys, by = c("dataset", "contrast"))
  df$comparison <- factor(df$comparison, levels = rev(pretty_order))
  df
}
fmt_fdr <- function(p) {
  vapply(p, function(z) {
    if (!is.finite(z)) return("NA")
    if (z < 0.001) formatC(z, format = "e", digits = 2)
    else formatC(z, format = "f", digits = 3)
  }, character(1))
}

# Validate one result per expected comparison before plotting.
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  match_n <- function(df) sum(df$dataset == k$dataset &
                             df$contrast == k$contrast, na.rm = TRUE)
  stopifnot(match_n(score_tests) == 1L)
  stopifnot(match_n(concordance[concordance$module_component == "ALL", ]) == 1L)
  stopifnot(match_n(camera) == 2L, match_n(fgsea) == 2L)
}
stopifnot(nrow(scores) == 15L, !anyDuplicated(scores$sample))

# ============================================================================
# A: fraction concordant, discordant, or below effect threshold
# For stimulation, "concordant" means OPPOSITE to the ALS directions.
# ============================================================================
ca <- concordance |>
  dplyr::filter(module_component == "ALL") |>
  as_comparison() |>
  dplyr::mutate(
    n_evaluable = as.numeric(n_evaluable),
    concordant = as.numeric(n_concordant_at_abs_log2FC_0_10),
    discordant = as.numeric(n_discordant_at_abs_log2FC_0_10),
    small_effect = as.numeric(n_small_effect)
  )
stopifnot(all(ca$concordant + ca$discordant + ca$small_effect == ca$n_evaluable))
ca_long <- ca |>
  dplyr::select(comparison, n_evaluable, concordant, discordant, small_effect) |>
  tidyr::pivot_longer(cols = c(concordant, discordant, small_effect),
                      names_to = "status", values_to = "count") |>
  dplyr::mutate(percent = 100 * count / n_evaluable,
                status = factor(status,
                     levels = c("concordant", "discordant", "small_effect")))
ca_labels <- ca |>
  dplyr::mutate(label = sprintf("%d/%d (%.1f%%)",
    concordant, n_evaluable, 100 * concordant/n_evaluable))

col_status <- c(concordant = "#274C75", discordant = "#C76360",
                small_effect = "#BAC5D0")

base_theme <- function() {
  ggplot2::theme_minimal(base_size = 12.5, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = margin(b = 7)),
      plot.subtitle = ggplot2::element_text(size = 10.5, colour = "#475569"),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(colour = "#263444", size = 10.5),
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = margin(8, 10, 8, 7),
      strip.text = ggplot2::element_text(face = "bold", size = 10.5),
      strip.background = ggplot2::element_rect(fill = "#F2F5F8", colour = NA)
    )
}

pA <- ggplot(ca_long, aes(x = percent, y = comparison, fill = status)) +
  geom_col(width = 0.64, colour = "white", linewidth = 0.3) +
  geom_text(data = ca_labels,
            aes(x = 102, y = comparison, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3.4,
            colour = "#172433", fontface = "bold") +
  scale_fill_manual(values = col_status,
    labels = c(concordant = "Expected direction",
               discordant = "Opposite direction",
               small_effect = "|log2 fold-change| < 0.10")) +
  scale_x_continuous(limits = c(0, 139), breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "A   Gene-level directional concordance",
       x = "Proportion of evaluable module genes (%)", y = NULL) +
  base_theme() +
  theme(legend.position = "bottom", legend.text = element_text(size = 9),
        legend.key.width = grid::unit(0.38, "cm"))

# ============================================================================
# B: fgsea NES (direction and strength) + CAMERA global FDR (annotations)
# These are different statistics; explicitly identify each in the panel.
# ============================================================================
cb <- camera |>
  dplyr::select(dataset, contrast, gene_set, Direction, FDR_all_tests) |>
  dplyr::inner_join(
    fgsea |>
      dplyr::select(dataset, contrast, pathway, NES),
    by = c("dataset", "contrast", "gene_set" = "pathway")
  ) |>
  as_comparison() |>
  dplyr::mutate(
    component = factor(gene_set, levels = c("ALS_DOWN", "ALS_UP"),
                       labels = c("ALS-down genes", "ALS-up genes")),
    NES = as.numeric(NES),
    FDR_all_tests = as.numeric(FDR_all_tests),
    label = paste0("NES ", sprintf("%+.2f", NES), "\nFDR ", fmt_fdr(FDR_all_tests)),
    label_colour = ifelse(abs(NES) >= 2.1, "white", "#1D2A38")
  )
stopifnot(nrow(cb) == 6L, all(is.finite(cb$NES)),
          all(is.finite(cb$FDR_all_tests)))

nes_bound <- max(abs(cb$NES)) * 1.08
pB <- ggplot(cb, aes(x = component, y = comparison, fill = NES)) +
  geom_tile(width = 0.96, height = 0.82, colour = "white", linewidth = 1.1) +
  geom_text(aes(label = label, colour = label_colour),
            size = 3.6, lineheight = 1.15, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#356B9B", mid = "#F9FAFB", high = "#C5644C",
    midpoint = 0, limits = c(-nes_bound, nes_bound),
    name = "fgsea NES") +
  labs(title = "B   Directional module enrichment",
       x = NULL, y = NULL) +
  base_theme() +
  theme(panel.grid = element_blank(), legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(face = "bold"))

# ============================================================================
# C: exact individual scores, with overlay showing group means
# ============================================================================
scores <- scores |>
  dplyr::mutate(
    module_score = as.numeric(module_score),
    dataset_label = factor(dataset,
      levels = c("GSE217576", "GSE117440"),
      labels = c("Mouse · GSE217576", "Rat · GSE117440")),
    group = factor(group, levels = c("CTL", "DENV", "SHAM", "DEN", "DN_SM"))
  )
stopifnot(all(is.finite(scores$module_score)))

score_colour <- c(CTL = "#7B8794", DENV = "#C45F52",
                  SHAM = "#7B8794", DEN = "#C45F52", DN_SM = "#2E4B83")
set.seed(20261014)
pC <- ggplot(scores, aes(x = group, y = module_score, colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.42,
             colour = "#A3ABB6") +
  geom_point(position = position_jitter(width = 0.075, height = 0),
             size = 3.05, alpha = 0.95) +
  stat_summary(fun = mean, geom = "point", colour = "#142638",
               shape = 95, size = 10, show.legend = FALSE) +
  facet_grid(. ~ dataset_label, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = score_colour, guide = "none") +
  scale_x_discrete(labels = c(CTL = "Control", DENV = "Denervated",
                              SHAM = "Sham", DEN = "Denervated",
                              DN_SM = "Denerv. +\nstimulation")) +
  labs(title = "C   Sample-level signed module scores",
       x = NULL, y = "Mean signed gene z-score") +
  base_theme() +
  theme(panel.grid.major.x = element_blank())

# ============================================================================
# D: reported mean-score differences, 95% CI and global FDR
# ============================================================================
sd <- score_tests |>
  as_comparison() |>
  dplyr::mutate(
    score_difference = as.numeric(score_difference),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    FDR_all_tests = as.numeric(FDR_all_tests),
    fdr_label = paste0("FDR ", fmt_fdr(FDR_all_tests)),
    type = ifelse(contrast == "Electrical_stimulation_vs_DEN",
                  "Stimulation", "Denervation")
  )
stopifnot(nrow(sd) == 3L, all(is.finite(sd$score_difference)),
          all(is.finite(sd$ci_low)), all(is.finite(sd$ci_high)))

label_x <- max(sd$ci_high) + 0.19
pD <- ggplot(sd, aes(y = comparison, x = score_difference)) +
  geom_vline(xintercept = 0, linewidth = 0.65, colour = "#A5ADB6",
             linetype = "dashed") +
  geom_segment(aes(x = ci_low, xend = ci_high,
                   y = comparison, yend = comparison, colour = type),
               linewidth = 1.1) +
  geom_point(aes(colour = type), size = 3.4) +
  geom_text(aes(x = label_x, label = fdr_label),
            hjust = 0, colour = "#263444", size = 3.7) +
  scale_colour_manual(values = c(Denervation = "#C45F52",
                                 Stimulation = "#2E4B83"), guide = "none") +
  scale_x_continuous(limits = c(min(-1.45, min(sd$ci_low) - 0.16),
                                  label_x + 1.04),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "D   Signed module-score differences",
       x = "Mean difference (numerator − denominator)", y = NULL) +
  base_theme()

# ============================================================================
# Compose and export; separate datasets are NEVER merged for scoring/analysis.
# ============================================================================
main <- (pA | pB) /
  (patchwork::free(pC, side = "l") | pD) +
  patchwork::plot_layout(heights = c(1, 1.12), widths = c(1.04, 1))

stem <- file.path(OUTDIR, "Figure_11_denervation")
ggsave(paste0(stem, ".png"), main, width = 13.5, height = 9.8,
       units = "in", dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem, ".svg"), main, width = 13.5, height = 9.8,
       units = "in", bg = "white", device = svglite::svglite)

# Audit log: derived from the workbook, not from manually entered results.
audit <- c(
  paste0("Source workbook: ", normalizePath(WORKBOOK, winslash = "/")),
  "A: Gene-level concordance (n_concordant / n_evaluable):",
  capture.output(print(ca[, c("dataset", "contrast", "concordant", "n_evaluable")])),
  "B: fgsea NES + CAMERA FDR:",
  capture.output(print(cb[, c("dataset", "contrast", "gene_set", "NES", "FDR_all_tests")])),
  "D: Module-score effects and intervals:",
  capture.output(print(sd[, c("dataset", "contrast", "score_difference",
                               "ci_low", "ci_high", "FDR_all_tests")]))
)
writeLines(audit, con = paste0(stem, "_source_audit.txt"))
message("Figure 11 created directly from: ", normalizePath(WORKBOOK, winslash = "/"))
message("Figure saved as PNG and SVG at: ", stem)
