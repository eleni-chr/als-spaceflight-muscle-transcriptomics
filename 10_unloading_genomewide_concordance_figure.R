#!/usr/bin/env Rscript
# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# =============================================================================
# ALS vs spaceflight | FIGURE 6 - Genome-wide concordance during unloading
# Replots final, already-computed human and mouse-to-human results from 09_unloading_genomewide_concordance_analysis.R.
# Includes ONE shared Genes per bin colourbar for A/B (identical scales).
# Layout fix: independent width allocations for the top and bottom rows prevent
# panel C's y-axis text from pushing panels A/B toward the right margin.
# No differential-expression analyses, orthologue mapping or correlations are
# refitted. Validate every summary statistic against its source gene-level data.
#
# A. Human: OSD-540 ULLS vs GSE271607 dry immersion (vastus lateralis)
# B. Human: GSE271607 dry immersion (vastus lateralis vs soleus)
# C. OSD-21 mouse gastrocnemius HLU vs three human unloading responses
#
# Expected input: CROSS_SPECIES_UNLOADING_CONCORDANCE_results.zip (or extracted
#                directory), produced by 09_unloading_genomewide_concordance_analysis.R.
# Saves PNG + SVG at publication size, plus individual figure panels.
# =============================================================================

# ------------------------------ USER OPTIONS --------------------------------
PROJECT_DIR <- getwd()
RESULTS_ZIP <- file.path(PROJECT_DIR, "CROSS_SPECIES_UNLOADING_CONCORDANCE_results.zip")
RESULTS_DIR <- file.path(PROJECT_DIR, "CROSS_SPECIES_UNLOADING_CONCORDANCE_results")
OUT_DIR <- file.path(PROJECT_DIR, "MANUSCRIPT_FIGURE_6_UNLOADING")
FIG_WIDTH_IN <- 8.8
FIG_HEIGHT_IN <- 8.2  # taller canvas accommodates one shared colourbar under A/B
DPI <- 450
PLOT_QUANTILE <- 0.995  # plotting only: rho is calculated on ALL shared genes
N_BINS <- 40L
RHO_CHECK_TOL <- 0.0006
# -----------------------------------------------------------------------------

for (pkg in c("ggplot2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install ", pkg, " with install.packages('", pkg, "').")
  }
}
if (!requireNamespace("svglite", quietly = TRUE)) {
  message("Package svglite not found. For the best SVG output, install.packages('svglite').")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Input discovery supports unpacked files as well as the original ZIP.
find_input <- function(relative_path) {
  subfolder <- "CROSS_SPECIES_UNLOADING_CONCORDANCE_results"
  candidates <- c(
    file.path(RESULTS_DIR, relative_path),
    file.path(RESULTS_DIR, subfolder, relative_path),
    file.path(PROJECT_DIR, subfolder, relative_path),
    file.path(PROJECT_DIR, "cross_species_review", subfolder, relative_path),
    file.path(PROJECT_DIR, relative_path)
  )
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(normalizePath(found[1L], winslash = "/", mustWork = TRUE))

  if (file.exists(RESULTS_ZIP)) {
    members <- utils::unzip(RESULTS_ZIP, list = TRUE)$Name
    member <- paste0(subfolder, "/", relative_path)
    if (!member %in% members) {
      stop("Required file not found in results ZIP: ", member)
    }
    extract_dir <- file.path(OUT_DIR, "_inputs_from_zip")
    dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(RESULTS_ZIP, files = member, exdir = extract_dir, overwrite = TRUE)
    path <- file.path(extract_dir, member)
    if (file.exists(path)) return(path)
  }
  stop("Cannot locate ", relative_path,
       ". Put CROSS_SPECIES_UNLOADING_CONCORDANCE_results.zip in PROJECT_DIR ",
       "or set RESULTS_DIR to the extracted analysis output folder.")
}

CORRELATIONS_FILE <- find_input("tables/all_five_genomewide_correlations.csv")
MATCHED_FILE <- find_input("tables/all_five_matched_gene_effects.csv")
SENSITIVITY_FILE <- find_input("tables/sensitivity_excluding_frozen_174_module.csv")

summary <- utils::read.csv(CORRELATIONS_FILE, check.names = FALSE)
matched <- utils::read.csv(MATCHED_FILE, check.names = FALSE)
sensitivity <- utils::read.csv(SENSITIVITY_FILE, check.names = FALSE)

require_cols <- function(d, cols, label) {
  missing <- setdiff(cols, colnames(d))
  if (length(missing)) stop("Missing columns in ", label, ": ", paste(missing, collapse=", "))
}
require_cols(summary, c("comparison", "n_shared", "spearman_rho"), "correlations")
require_cols(matched, c("comparison", "human_symbol", "log2FC_first", "log2FC_second"), "matched genes")
require_cols(sensitivity, c("comparison", "rho_without_174", "n_remaining"), "sensitivity")

expected_keys <- c(
  "HUMAN_ULLS_vs_DI_VL",
  "HUMAN_DI_VL_vs_DI_SOLEUS",
  "MOUSE_HLU_vs_HUMAN_ULLS_VL",
  "MOUSE_HLU_vs_HUMAN_DI_VL",
  "MOUSE_HLU_vs_HUMAN_DI_SOLEUS"
)
if (nrow(summary) != 5L || !setequal(summary$comparison, expected_keys)) {
  stop("Expected exactly the five prespecified human-human/mouse-human comparisons.")
}
if (!setequal(unique(matched$comparison), expected_keys)) {
  stop("The matched-gene file has missing or unexpected comparisons.")
}
if (!setequal(sensitivity$comparison, expected_keys)) {
  stop("The sensitivity file has missing or unexpected comparisons.")
}
if (anyDuplicated(paste(matched$comparison, matched$human_symbol, sep = "||"))) {
  stop("Duplicated comparison/gene entries: cannot safely compute or plot concordance.")
}

# QC: verify the published rho and n from the exported raw matched gene pairs.
qc <- do.call(rbind, lapply(expected_keys, function(key) {
  d <- matched[matched$comparison == key &
                 is.finite(matched$log2FC_first) &
                 is.finite(matched$log2FC_second), , drop = FALSE]
  r <- summary[summary$comparison == key, , drop = FALSE]
  recalculated <- stats::cor(d$log2FC_first, d$log2FC_second, method = "spearman")
  data.frame(comparison = key,
             n_reported = r$n_shared, n_from_matched_table = nrow(d),
             rho_reported = r$spearman_rho, rho_recalculated = recalculated,
             abs_rho_difference = abs(r$spearman_rho - recalculated),
             stringsAsFactors = FALSE)
}))
utils::write.csv(qc, file.path(OUT_DIR, "Figure_6_correlation_input_QC.csv"), row.names = FALSE)
if (any(qc$n_reported != qc$n_from_matched_table) ||
    any(!is.finite(qc$rho_recalculated)) ||
    any(qc$abs_rho_difference > RHO_CHECK_TOL)) {
  stop("Figure 6 source-table validation failed. Check Figure_6_correlation_input_QC.csv")
}
if (any(!is.finite(matched$log2FC_first)) || any(!is.finite(matched$log2FC_second))) {
  stop("Matched gene-effect table contains non-finite values.")
}
message("Source QC passed for all five comparisons.")

# ---- VISUAL STYLE ----------------------------------------------------------
font_family <- "sans"
ink <- "#22303B"
blue <- "#286991"

# All titles intentionally short so A/B don't collide at ~17-cm print width.
make_density <- function(key, title, xlabel, ylabel) {
  d <- matched[matched$comparison == key, , drop = FALSE]
  r <- summary[summary$comparison == key, , drop = FALSE]

  # 99.5th percentile of |log2FC| among BOTH paired contrast vectors.
  # This is only a visualization window. Spearman rho + n use every matched gene.
  lim <- as.numeric(stats::quantile(abs(c(d$log2FC_first, d$log2FC_second)),
                                    probs = PLOT_QUANTILE, na.rm = TRUE))
  lim <- max(0.5, lim)
  view <- d[abs(d$log2FC_first) <= lim & abs(d$log2FC_second) <= lim, , drop=FALSE]

  ggplot2::ggplot(view, ggplot2::aes(x = log2FC_first, y = log2FC_second)) +
    # A/B will receive the SAME count-to-colour scale after both plots
    # are constructed. This permits one correctly shared colourbar.
    ggplot2::geom_bin_2d(bins = N_BINS, show.legend = TRUE) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed",
                         linewidth=0.4, colour="#7A8085") +
    ggplot2::geom_hline(yintercept=0, linewidth=0.25, colour="#A4A9AD") +
    ggplot2::geom_vline(xintercept=0, linewidth=0.25, colour="#A4A9AD") +
    ggplot2::coord_equal(xlim=c(-lim, lim), ylim=c(-lim, lim), expand=FALSE) +
    ggplot2::labs(
      title=title,
      subtitle=sprintf("\u03c1 = %.3f   |   n = %s genes",
                       r$spearman_rho, format(r$n_shared, big.mark=",", trim=TRUE)),
      x=xlabel, y=ylabel
    ) +
    ggplot2::theme_bw(base_size=12.6, base_family=font_family) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(size=14.2, face="bold", colour=ink,
                                       margin=ggplot2::margin(b=5)),
      plot.subtitle=ggplot2::element_text(size=11.6, colour=ink,
                                          margin=ggplot2::margin(b=9)),
      axis.title.x=ggplot2::element_text(size=11.9, face="plain", colour=ink,
                                         margin=ggplot2::margin(t=7)),
      axis.title.y=ggplot2::element_text(size=11.9, face="plain", colour=ink,
                                         margin=ggplot2::margin(r=6)),
      axis.text=ggplot2::element_text(size=11, colour=ink),
      panel.grid.major=ggplot2::element_line(colour="#EDEFF1", linewidth=0.25),
      panel.grid.minor=ggplot2::element_blank(),
      plot.title.position="plot",
      plot.subtitle.position="plot",
      plot.margin=ggplot2::margin(9,16,4,9)
    )
}

fig6A <- make_density(
  "HUMAN_ULLS_vs_DI_VL", "A. Human ULLS vs dry immersion",
  "ULLS, vastus lateralis (log2 fold change)", "DI, vastus lateralis (log2 fold change)")
fig6B <- make_density(
  "HUMAN_DI_VL_vs_DI_SOLEUS", "B. Human dry immersion: VL vs soleus",
  "Vastus lateralis (log2 fold change)", "Soleus (log2 fold change)")

# SHARED A/B LEGEND ----------------------------------------------------------
# geom_bin_2d counts the number of genes in each bin; the density colours are
# NOT probabilities. Compute the maximum observed bin count in either panel,
# then apply exactly the same count-to-colour limits and breaks to both.
# Without fixed limits, ggplot would scale A/B separately, making a single
# collected colourbar scientifically misleading.
count_A <- ggplot2::ggplot_build(fig6A)$data[[1L]]$count
count_B <- ggplot2::ggplot_build(fig6B)$data[[1L]]$count
common_max_count <- max(c(count_A, count_B), na.rm = TRUE)
if (!is.finite(common_max_count) || common_max_count < 1) {
  stop("No non-empty 2D bins for Figure 6A/B; cannot create colour legend.")
}

# Always include the shared maximum in the legend; use integer counts.
count_breaks <- unique(as.integer(round(seq(0, common_max_count, length.out = 5))))
count_breaks <- count_breaks[count_breaks >= 0 & count_breaks <= common_max_count]
shared_bin_fill <- function() {
  ggplot2::scale_fill_gradient(
    low = "#E1EBF2", high = "#245D86",
    limits = c(0, common_max_count),
    breaks = count_breaks,
    name = "Genes per bin",
    guide = ggplot2::guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(65, "mm"),
      barheight = grid::unit(3.5, "mm")
    )
  )
}
fig6A <- fig6A + shared_bin_fill() +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 11.5, colour=ink),
                 legend.text = ggplot2::element_text(size = 10.8, colour=ink))
fig6B <- fig6B + shared_bin_fill() +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 11.5, colour=ink),
                 legend.text = ggplot2::element_text(size = 10.8, colour=ink))

c_keys <- expected_keys[3:5]
c_labels <- c(
  "Human ULLS (VL)",
  "Human DI (VL)",
  "Human DI (soleus)"
)
m <- summary[match(c_keys, summary$comparison), , drop = FALSE]
m$label <- factor(c_labels, levels=rev(c_labels))
m$annotation <- sprintf("\u03c1 = %.3f  |  n = %s",
                        m$spearman_rho,
                        format(m$n_shared, big.mark=",", trim=TRUE))

# Keep C on its own full-width row; text column is to the right of the plotted rho.
fig6C <- ggplot2::ggplot(m, ggplot2::aes(y=label, x=spearman_rho)) +
  ggplot2::geom_vline(xintercept=0, linewidth=0.65, colour="#7E8B96") +
  ggplot2::geom_segment(ggplot2::aes(x=0, xend=spearman_rho, yend=label),
                        linewidth=1.1, colour="#92B3C9") +
  ggplot2::geom_point(shape=21, size=4.5, stroke=0.65,
                      fill=blue, colour="#184B70") +
  ggplot2::geom_text(ggplot2::aes(label=annotation), x=0.215, hjust=0,
                     colour=ink, size=4.15, family=font_family) +
  ggplot2::scale_x_continuous(limits=c(-0.012, 0.48),
                              breaks=c(0,0.1,0.2,0.3,0.4),
                              labels=c("0", "0.1", "0.2", "0.3", "0.4"),
                              expand=ggplot2::expansion(mult=0)) +
  ggplot2::labs(title="C  Mouse HLU vs human unloading",
                x="Spearman's \u03c1 across shared genes", y=NULL) +
  ggplot2::theme_minimal(base_size=13, base_family=font_family) +
  ggplot2::theme(
    plot.title=ggplot2::element_text(size=14.5, face="bold", colour=ink),
    plot.subtitle=ggplot2::element_text(size=11.9, colour=ink),
    axis.text.y=ggplot2::element_text(size=12.0, colour=ink,
                                      margin=ggplot2::margin(r=8)),
    axis.text.x=ggplot2::element_text(size=11.5, colour=ink),
    axis.title.x=ggplot2::element_text(size=12.7, colour=ink,
                                       margin=ggplot2::margin(t=8)),
    panel.grid.major.y=ggplot2::element_blank(),
    panel.grid.minor=ggplot2::element_blank(),
    panel.grid.major.x=ggplot2::element_line(colour="#E5E9EC"),
    plot.title.position="plot",
    plot.subtitle.position="plot",
    plot.margin=ggplot2::margin(7,16,10,9)
  )

# IMPORTANT LAYOUT FIX ------------------------------------------------------
# Panel C's y-axis labels must NOT set the common left-axis width for the
# top row. Otherwise ggplot/patchwork reserves a huge left blank margin above
# A/B and their titles can collide. Let panel C keep its own outer text area.
# wrap_elements() isolates the bottom plot's axes from the top-row gtable
# widths. It is compatible with older patchwork releases too.
fig6C_independent <- patchwork::wrap_elements(full=fig6C)

# The top row uses two equal-width panels with shortened titles. The bottom
# uses its own alignment so the C labels do not shift A/B to the right.
# Collect the IDENTICAL A/B fill guides within the top row, and show one
# horizontal colourbar below A and B. The wrapped lower row is isolated and
# contributes no legend. Do not collect guides at the outer patchwork level.
fig6_top <- (fig6A | fig6B) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
fig6 <- (fig6_top / fig6C_independent) +
  patchwork::plot_layout(heights=c(1.35, 0.86)) &
  ggplot2::theme(plot.background=ggplot2::element_rect(fill="white", colour=NA))

save_plot <- function(plot, stem, width, height) {
  png_name <- file.path(OUT_DIR, paste0(stem, ".png"))
  svg_name <- file.path(OUT_DIR, paste0(stem, ".svg"))
  ggplot2::ggsave(png_name, plot, width=width, height=height, units="in",
                  dpi=DPI, bg="white", limitsize=FALSE)
  if (requireNamespace("svglite", quietly=TRUE)) {
    ggplot2::ggsave(svg_name, plot, width=width, height=height, units="in",
                    device=svglite::svglite, bg="white", limitsize=FALSE)
  } else {
    grDevices::svg(svg_name, width=width, height=height, bg="white")
    print(plot)
    grDevices::dev.off()
  }
  message("Saved: ", png_name, " and ", svg_name)
}

save_plot(fig6, "Figure_6_UNLOADING_GENOMEWIDE_CONCORDANCE_REVISED",
          FIG_WIDTH_IN, FIG_HEIGHT_IN)
save_plot(fig6A, "Figure_6A_human_ULLS_vs_DI_VL", 6.0, 5.5)
save_plot(fig6B, "Figure_6B_human_DI_VL_vs_soleus", 6.0, 5.5)
save_plot(fig6C, "Figure_6C_mouse_vs_human_unloading", 8.8, 3.7)

utils::write.csv(summary[match(expected_keys,summary$comparison),],
                 file.path(OUT_DIR,"Figure_6_all_five_correlations.csv"), row.names=FALSE)
utils::write.csv(sensitivity[match(expected_keys,sensitivity$comparison),],
                 file.path(OUT_DIR,"Figure_6_module_exclusion_sensitivity.csv"), row.names=FALSE)
writeLines(c(
  "FIGURE 6 - genome-wide gene-level Spearman correlations during terrestrial unloading",
  "A/B: 2D bin counts of paired gene-level log2FC. One shared 'Genes per bin' colourbar.",
  paste0("A/B common fill range: 0 to ", common_max_count, " genes per bin."),
  "Only central plotting window of effects is shown (99.5th percentile of absolute effects).",
  "All reported rho and n are from ALL shared evaluable genes BEFORE plotting restriction.",
  "C: points and horizontal connectors indicate observed rho, NOT confidence intervals.",
  "Figure layout: abbreviated A/B titles; C labels kept independent of top-row widths.",
  "VL = vastus lateralis; ULLS = unilateral lower-limb suspension; DI = dry immersion.",
  "Mouse effects: OSD-21 hindlimb unloading in mouse gastrocnemius, mapped to human orthologues.",
  "The five tests are descriptive gene-level comparisons; nominal P values do not assess",
  "variation among biological replicates or independence of genes.",
  "All 174-module exclusion results are in the CSV sensitivity table, not plotted in Figure 6.",
  paste0("Correlation input: ", CORRELATIONS_FILE),
  paste0("Matched gene-level input: ", MATCHED_FILE),
  paste0("Module sensitivity input: ", SENSITIVITY_FILE)
), file.path(OUT_DIR,"Figure_6_PROVENANCE.txt"))
message("DONE. Figure 6 PNG and SVG saved in: ", OUT_DIR)
