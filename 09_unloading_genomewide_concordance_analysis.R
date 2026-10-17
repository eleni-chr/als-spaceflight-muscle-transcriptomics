#!/usr/bin/env Rscript
# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# =============================================================================
# ALS vs spaceflight | Figure 6: genome-wide unloading concordance
# Analysis using the gene-level differential-expression tables from script 02.
#
# A: OSD-540 human VL (ULLS vs baseline) vs GSE271607 human VL (DI after vs before)
# B: GSE271607 human VL vs human soleus (DI after vs before)
# C: OSD-21 mouse gastrocnemius (HLU vs control) vs each of the three human
#    unloading contrasts above (OSD-540 VL; GSE271607 VL; GSE271607 soleus).
#
# No differential-expression or module-selection model is rerun. For cross-
# species comparison mouse gene symbols are harmonised to HUMAN symbols using
# org.Mm.eg.db + babelgene, matching the established module workflow.
# Correlations are over all common evaluable genes, NOT DEGs or just 174 genes.
# For transparency, a second analysis excludes the frozen 174-gene module.
#
# Results are exploratory. Gene-level Spearman P-values do NOT test animal-level
# replication and cannot resolve platform, species, muscle or intervention effects.
# =============================================================================

# -------------------------- USER SETTINGS ----------------------------------
PROJECT_DIR <- getwd()            # Set to the analysis project root if needed
OUT_DIR <- file.path(PROJECT_DIR, "CROSS_SPECIES_UNLOADING_CONCORDANCE_results")

# Point directly to files if your directory structure differs; otherwise leave
# as NA and the script will find the exact basenames in known output folders.
DE_PATHS <- c(
  OSD21_HLU      = NA_character_,
  OSD540_ULLS    = NA_character_,
  GSE271607_VL   = NA_character_,
  GSE271607_SOL  = NA_character_
)
DE_BASENAMES <- c(
  OSD21_HLU     = "OSD21_HLU_vs_CTRL_Gene_symbol_level.csv",
  OSD540_ULLS   = "OSD540_ULLS_vs_Baseline_Gene_symbol_level.csv",
  GSE271607_VL  = "GSE271607_VL_After_vs_Before_Gene_symbol_level.csv",
  GSE271607_SOL = "GSE271607_Soleus_After_vs_Before_Gene_symbol_level.csv"
)

# Optional explicit FROZEN module path. If NA, search the unified results
# directory and results.zip; no module is defined from these new correlations.
MODULE_PATH <- NA_character_
EXCLUDE_174_SENSITIVITY <- TRUE

# Matches the published module workflow's orthologue source and evidence rule.
ORTHOLOGY_MIN_SUPPORT <- 3L
# Conservatively exclude human symbols reached from >1 measured mouse symbol.
# This is 'unique among observed mappings', not a proof of reciprocal orthology.
REQUIRE_UNIQUE_MOUSE_TO_HUMAN <- TRUE

# Reproducibility check against the two original manuscript human comparisons.
EXPECTED_HUMAN <- data.frame(
  comparison = c("HUMAN_ULLS_vs_DI_VL", "HUMAN_DI_VL_vs_DI_SOLEUS"),
  expected_n = c(14637L, 14725L),
  expected_rho = c(0.396, 0.600),
  stringsAsFactors = FALSE
)
HUMAN_N_TOLERANCE <- 0.05       # >5% change in matched genes requires review
HUMAN_RHO_TOLERANCE <- 0.03    # >0.03 change in rho requires review

FIG_WIDTH <- 10.5               # Compact canvas; keeps labels legible at ~17 cm
FIG_HEIGHT <- 9.6
DPI <- 450
PLOT_QUANTILE <- 0.995         # Extreme effects clipped in plot ONLY
# ---------------------------------------------------------------------------

required <- c("dplyr", "tibble", "ggplot2", "patchwork", "babelgene",
              "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db")
missing_packages <- required[!vapply(required, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing packages before running:\n  ",
       paste(missing_packages, collapse = ", "),
       "\nFor Bioconductor packages use BiocManager::install(...).")
}
if (!requireNamespace("svglite", quietly = TRUE)) {
  warning("svglite not installed: SVG output will use grDevices::svg().")
}

DIR_TABLE <- file.path(OUT_DIR, "tables")
DIR_FIG <- file.path(OUT_DIR, "figures")
DIR_QC <- file.path(OUT_DIR, "QC")
for (x in c(DIR_TABLE, DIR_FIG, DIR_QC)) {
  dir.create(x, showWarnings = FALSE, recursive = TRUE)
}

message("\n=== Genome-wide unloading concordance: 5 prespecified comparisons ===")
message("PROJECT_DIR: ", normalizePath(PROJECT_DIR, winslash = "/"))

# ------------------------- FILE DISCOVERY ----------------------------------
search_roots <- unique(c(
  file.path(PROJECT_DIR, "NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL"),
  file.path(PROJECT_DIR, "NASA_ALS_metabolic_pathway_results_ALS_DERIVED_HUMAN_BIOPSIES"),
  file.path(PROJECT_DIR, "NASA_ALS_metabolic_pathway_results_ALS_DERIVED"),
  file.path(PROJECT_DIR, "nasa_als_results"),
  file.path(PROJECT_DIR, "reversible_ALS_unloading_results"),
  PROJECT_DIR
))
search_roots <- search_roots[dir.exists(search_roots)]

find_exact_basename <- function(filename) {
  # Direct locations are much faster than walking every raw-data subdirectory.
  subdirs <- c("", "DE", "tables", "Results", "gene_level", "DE_results")
  for (root in search_roots) {
    probes <- file.path(root, subdirs, filename)
    found <- probes[file.exists(probes)]
    if (length(found)) {
      return(normalizePath(found[1], mustWork = TRUE, winslash = "/"))
    }
  }
  for (root in search_roots) {
    matches <- list.files(root, pattern = "\\.csv$", full.names = TRUE,
                          recursive = TRUE, ignore.case = TRUE)
    exact <- matches[tolower(basename(matches)) == tolower(filename)]
    if (length(exact)) {
      exact <- exact[order(nchar(exact))]
      if (length(exact) > 1) warning("Multiple copies of ", filename,
                                   "; using: ", exact[1])
      return(normalizePath(exact[1], mustWork = TRUE, winslash = "/"))
    }
  }
  NA_character_
}

for (nm in names(DE_BASENAMES)) {
  explicit <- DE_PATHS[[nm]]
  if (!is.na(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) stop("Input ", nm, " does not exist: ", explicit)
    DE_PATHS[[nm]] <- normalizePath(explicit, mustWork = TRUE, winslash = "/")
  } else {
    DE_PATHS[[nm]] <- find_exact_basename(DE_BASENAMES[[nm]])
  }
}
if (anyNA(DE_PATHS)) {
  missing <- names(DE_PATHS)[is.na(DE_PATHS)]
  stop("\nMissing existing gene-level differential-expression result(s):\n",
       paste(paste0("  ", missing, " : ", DE_BASENAMES[missing]), collapse = "\n"),
       "\n\nCopy the four final DE CSV files to PROJECT_DIR, or set their ",
       "full paths in DE_PATHS at the top of this script.\n",
       "Do not substitute actual-flight OSD-21 output for HLU_vs_CTRL.")
}
input_manifest <- data.frame(dataset = names(DE_PATHS), file = unname(DE_PATHS),
                             stringsAsFactors = FALSE)
write.csv(input_manifest, file.path(DIR_QC, "DE_input_manifest.csv"), row.names = FALSE)
for (i in seq_along(DE_PATHS)) message("  ", names(DE_PATHS)[i], ": ", DE_PATHS[i])

# ------------------------ STANDARDISE DE FILES -----------------------------
normalize_name <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))
clean_symbol <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[out %in% c("", "NA", "N/A") | is.na(out)] <- NA_character_
  out
}
pick_col <- function(df, candidates, required = TRUE, source = "") {
  norm <- normalize_name(names(df))
  idx <- match(normalize_name(candidates), norm)
  idx <- idx[!is.na(idx)]
  if (length(idx)) return(names(df)[idx[1]])
  if (required) stop("Cannot find gene ID or log2FC column in ", source,
                     ". Available columns: ", paste(names(df), collapse = ", "))
  NULL
}

read_de <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  sym_col <- pick_col(x, c("symbol", "gene_symbol", "genesymbol", "gene"),
                      source = basename(path))
  lfc_col <- pick_col(x, c("log2FC", "logFC", "log2FoldChange", "LFC"),
                      source = basename(path))
  p_col <- pick_col(x, c("p_value", "pvalue", "p_val", "pval", "P.Value"),
                    required = FALSE)
  # The comparison effect direction is inherited unchanged from the exact DE
  # filename. No arbitrary sign inversion or DE-significance filtering.
  out <- data.frame(
    source_symbol = trimws(as.character(x[[sym_col]])),
    log2FC = suppressWarnings(as.numeric(x[[lfc_col]])),
    p_value = if (!is.null(p_col)) suppressWarnings(as.numeric(x[[p_col]]))
              else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$log2FC) & !is.na(out$source_symbol) &
               nzchar(out$source_symbol), , drop = FALSE]
  if (!nrow(out)) stop("No evaluable genes in: ", path)
  out
}

# Collapse repeated symbols conservatively. Never double-count multiple
# Ensembl IDs that map to one human symbol. Discard conflicting signed effects.
collapse_effects <- function(df, group_col, tag) {
  d <- df[!is.na(df[[group_col]]) & nzchar(df[[group_col]]) &
            is.finite(df$log2FC), , drop = FALSE]
  names(d)[names(d) == group_col] <- "human_symbol"
  count_before <- nrow(d)
  d <- tibble::as_tibble(d) |>
    dplyr::group_by(.data$human_symbol) |>
    dplyr::filter(dplyr::n_distinct(sign(.data$log2FC[.data$log2FC != 0])) <= 1) |>
    dplyr::summarise(log2FC = stats::median(.data$log2FC),
                     source_rows = dplyr::n(), .groups = "drop")
  message("  ", tag, ": ", count_before, " mapped rows -> ", nrow(d),
          " unique human gene symbols (ambiguous duplicate effects excluded)")
  as.data.frame(d)
}

human_harmonise <- function(raw, tag) {
  ids <- trimws(raw$source_symbol)
  ens <- sub("\\.[0-9]+$", "", ids)
  is_ensg <- grepl("^ENSG[0-9]+$", toupper(ens))
  human <- clean_symbol(ids)
  if (any(is_ensg)) {
    keys <- unique(ens[is_ensg])
    ann <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                 keys = keys, keytype = "ENSEMBL",
                                 column = "SYMBOL", multiVals = "first")
    human[is_ensg] <- clean_symbol(unname(ann[ens[is_ensg]]))
  }
  raw$human_symbol <- human
  collapse_effects(raw, "human_symbol", tag)
}

# ------------------------ MOUSE ORTHOLOGUE MAP ------------------------------
mouse_harmonise <- function(raw) {
  input <- unique(clean_symbol(raw$source_symbol))
  input <- input[!is.na(input)]
  mm_official <- AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db,
                                     keytype = "SYMBOL")
  official_lut <- stats::setNames(mm_official, toupper(mm_official))
  source <- data.frame(
    mouse_symbol = input,
    mouse_symbol_official = unname(official_lut[input]),
    stringsAsFactors = FALSE
  )
  source <- source[!is.na(source$mouse_symbol_official), , drop = FALSE]
  if (!nrow(source)) stop("No mouse official SYMBOL mappings. Check that OSD-21 ",
                         "file is a gene-symbol-level file, not Ensembl-only.")
  message("  Mouse input: ", length(input), " unique symbols; ", nrow(source),
          " recognised as official mouse SYMBOLs")
  ortho <- babelgene::orthologs(
    genes = unique(source$mouse_symbol_official),
    species = "mouse", human = FALSE,
    min_support = ORTHOLOGY_MIN_SUPPORT, top = TRUE
  )
  needed <- c("symbol", "human_symbol", "support_n")
  if (!all(needed %in% names(ortho))) {
    stop("Unexpected babelgene columns: ", paste(names(ortho), collapse = ", "))
  }
  ortho <- ortho[, needed, drop = FALSE]
  names(ortho) <- c("mouse_symbol_official", "human_symbol", "support_n")
  ortho$human_symbol <- clean_symbol(ortho$human_symbol)
  ortho <- ortho[!is.na(ortho$human_symbol), , drop = FALSE]
  ortho <- ortho[order(ortho$mouse_symbol_official,
                       -as.numeric(ortho$support_n)), , drop = FALSE]
  ortho <- ortho[!duplicated(ortho$mouse_symbol_official), , drop = FALSE]
  map <- merge(source, ortho, by = "mouse_symbol_official", all.x = TRUE, sort = FALSE)
  map <- map[!is.na(map$human_symbol), , drop = FALSE]
  n_mapped_before <- nrow(map)
  if (REQUIRE_UNIQUE_MOUSE_TO_HUMAN) {
    n_per_human <- table(map$human_symbol)
    map <- map[unname(n_per_human[map$human_symbol]) == 1L, , drop = FALSE]
  }
  if (anyDuplicated(map$human_symbol)) stop("Orthologue map not unique by human gene")
  message("  Orthologues: ", n_mapped_before, " mapped mouse genes; ", nrow(map),
          " retained after uniqueness QC")
  write.csv(map, file.path(DIR_QC, "OSD21_mouse_to_human_orthologue_map.csv"),
            row.names = FALSE)

  raw$mouse_symbol <- clean_symbol(raw$source_symbol)
  joined <- merge(raw, map[, c("mouse_symbol", "human_symbol")],
                  by = "mouse_symbol", all.x = FALSE, all.y = FALSE,
                  sort = FALSE)
  collapse_effects(joined, "human_symbol", "OSD21 mouse HLU")
}

message("\nLoading FINAL gene-level DE effects...")
effects <- list(
  OSD21_HLU     = mouse_harmonise(read_de(DE_PATHS[["OSD21_HLU"]])),
  OSD540_ULLS   = human_harmonise(read_de(DE_PATHS[["OSD540_ULLS"]]),
                                "OSD540 human VL"),
  GSE271607_VL  = human_harmonise(read_de(DE_PATHS[["GSE271607_VL"]]),
                                "GSE271607 human VL"),
  GSE271607_SOL = human_harmonise(read_de(DE_PATHS[["GSE271607_SOL"]]),
                                "GSE271607 human soleus")
)

qc_counts <- data.frame(dataset = names(effects),
                        n_unique_human_symbols = vapply(effects, nrow, integer(1)))
write.csv(qc_counts, file.path(DIR_QC, "gene_universe_counts.csv"), row.names = FALSE)
print(qc_counts, row.names = FALSE)

# ------------------------ FIVE PLANNED COMPARISONS --------------------------
comparison_spec <- data.frame(
  comparison = c("HUMAN_ULLS_vs_DI_VL", "HUMAN_DI_VL_vs_DI_SOLEUS",
                 "MOUSE_HLU_vs_HUMAN_ULLS_VL", "MOUSE_HLU_vs_HUMAN_DI_VL",
                 "MOUSE_HLU_vs_HUMAN_DI_SOLEUS"),
  first = c("OSD540_ULLS", "GSE271607_VL", "OSD21_HLU", "OSD21_HLU", "OSD21_HLU"),
  second = c("GSE271607_VL", "GSE271607_SOL", "OSD540_ULLS",
             "GSE271607_VL", "GSE271607_SOL"),
  group = c("Human vs human", "Human vs human", rep("Mouse vs human", 3)),
  plot_label = c("ULLS (VL) vs dry immersion (VL)",
                 "Dry immersion: VL vs soleus",
                 "Mouse HLU vs human ULLS (VL)",
                 "Mouse HLU vs human dry immersion (VL)",
                 "Mouse HLU vs human dry immersion (soleus)"),
  stringsAsFactors = FALSE
)

pairwise_genes <- list()
cor_results <- vector("list", nrow(comparison_spec))

for (i in seq_len(nrow(comparison_spec))) {
  row <- comparison_spec[i, ]
  x <- effects[[row$first]][, c("human_symbol", "log2FC")]
  y <- effects[[row$second]][, c("human_symbol", "log2FC")]
  names(x)[2] <- "log2FC_first"
  names(y)[2] <- "log2FC_second"
  paired <- merge(x, y, by = "human_symbol", sort = FALSE)
  paired <- paired[is.finite(paired$log2FC_first) &
                     is.finite(paired$log2FC_second), , drop = FALSE]
  if (anyDuplicated(paired$human_symbol)) {
    stop("Duplicate human gene keys after join: ", row$comparison)
  }
  if (nrow(paired) < 100L) {
    stop("Too few shared genes (", nrow(paired), ") for ", row$comparison,
         ". Inspect orthologue mapping and source columns.")
  }
  if (stats::sd(paired$log2FC_first) == 0 ||
      stats::sd(paired$log2FC_second) == 0) {
    stop("Constant fold-change vector in ", row$comparison)
  }
  test <- suppressWarnings(stats::cor.test(paired$log2FC_first,
                    paired$log2FC_second, method = "spearman", exact = FALSE))
  rho <- unname(test$estimate)
  paired$comparison <- row$comparison
  paired$first <- row$first
  paired$second <- row$second
  pairwise_genes[[i]] <- paired
  cor_results[[i]] <- data.frame(
    comparison = row$comparison, first = row$first, second = row$second,
    group = row$group, n_shared = nrow(paired), spearman_rho = rho,
    spearman_p_approx = test$p.value,
    stringsAsFactors = FALSE
  )
  message(sprintf("  %-37s n=%6d  rho=% .4f  p=%s", row$comparison,
                  nrow(paired), rho, format.pval(test$p.value, digits = 3)))
}

cor_results <- do.call(rbind, cor_results)
cor_results$BH_FDR_across_5_tests <- p.adjust(cor_results$spearman_p_approx,
                                            method = "BH")
all_pairs <- do.call(rbind, pairwise_genes)
write.csv(cor_results, file.path(DIR_TABLE, "all_five_genomewide_correlations.csv"),
          row.names = FALSE)
write.csv(all_pairs, file.path(DIR_TABLE, "all_five_matched_gene_effects.csv"),
          row.names = FALSE)

# --------------------------- REPRODUCIBILITY -------------------------------
expected <- merge(EXPECTED_HUMAN, cor_results, by = "comparison", all.x = TRUE)
expected$n_difference <- expected$n_shared - expected$expected_n
expected$n_percent_difference <- 100 * expected$n_difference / expected$expected_n
expected$rho_difference <- expected$spearman_rho - expected$expected_rho
expected$within_tolerance <- abs(expected$n_percent_difference) <=
                               (100 * HUMAN_N_TOLERANCE) &
                             abs(expected$rho_difference) <= HUMAN_RHO_TOLERANCE
write.csv(expected, file.path(DIR_QC, "human_human_original_values_check.csv"),
          row.names = FALSE)
if (any(!expected$within_tolerance)) {
  warning("Human–human correlations differ materially from manuscript values ",
          "(n=14637/rho=.396; n=14725/rho=.600). See QC/human_human_original_values_check.csv.",
          " Reconcile preprocessing before updating the manuscript.")
} else {
  message("Human–human comparisons agree with the original manuscript within tolerance.")
}

# ----------------- FROZEN 174-GENE EXCLUSION SENSITIVITY -------------------
read_module <- function() {
  if (!is.na(MODULE_PATH) && nzchar(MODULE_PATH)) {
    if (!file.exists(MODULE_PATH)) stop("MODULE_PATH missing: ", MODULE_PATH)
    z <- read.csv(MODULE_PATH, stringsAsFactors = FALSE, check.names = FALSE)
    src <- MODULE_PATH
  } else {
    candidates <- c(
      file.path(PROJECT_DIR, "UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results", "QC",
                "FROZEN_module_174.csv"),
      file.path(PROJECT_DIR, "results", "QC", "FROZEN_module_174.csv"),
      file.path(PROJECT_DIR, "QC", "FROZEN_module_174.csv"),
      file.path(PROJECT_DIR, "FROZEN_module_174.csv"))
    found <- candidates[file.exists(candidates)]
    if (length(found)) {
      src <- found[1]
      z <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      zip_candidates <- c(file.path(PROJECT_DIR, "results.zip"),
                          file.path(PROJECT_DIR,
                                    "UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results", "results.zip"))
      found_zip <- zip_candidates[file.exists(zip_candidates)]
      if (!length(found_zip)) return(NULL)
      zz <- found_zip[1]
      members <- unzip(zz, list = TRUE)$Name
      member <- members[members == "QC/FROZEN_module_174.csv"]
      if (!length(member)) return(NULL)
      dest <- file.path(DIR_QC, "_frozen_module_input")
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      unzip(zz, files = member[1], exdir = dest, overwrite = TRUE)
      src <- file.path(dest, member[1])
      z <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  scol <- pick_col(z, c("human_symbol", "symbol", "gene_symbol"),
                   source = src)
  genes <- unique(clean_symbol(z[[scol]]))
  genes <- genes[!is.na(genes)]
  if (length(genes) != 174L) {
    stop("Frozen module contains ", length(genes), " genes rather than 174: ",
         src, ". Do not use an alternative module without documenting it.")
  }
  message("Loaded frozen 174-gene module from: ", src)
  genes
}

module_genes <- if (EXCLUDE_174_SENSITIVITY) read_module() else NULL
if (EXCLUDE_174_SENSITIVITY && is.null(module_genes)) {
  warning("Frozen 174-gene list not found. Primary five comparisons and Fig 6 ",
          "will still be generated; module-excluded sensitivity is omitted. ",
          "Place the unified results.zip in PROJECT_DIR or set MODULE_PATH.")
}

if (!is.null(module_genes)) {
  sensitivity <- vector("list", length(pairwise_genes))
  for (i in seq_along(pairwise_genes)) {
    x <- pairwise_genes[[i]]
    remaining <- x[!x$human_symbol %in% module_genes, , drop = FALSE]
    test <- suppressWarnings(stats::cor.test(
      remaining$log2FC_first, remaining$log2FC_second,
      method = "spearman", exact = FALSE))
    sensitivity[[i]] <- data.frame(
      comparison = comparison_spec$comparison[i],
      n_original = nrow(x), n_module_removed = nrow(x) - nrow(remaining),
      n_remaining = nrow(remaining),
      original_rho = cor_results$spearman_rho[i],
      rho_without_174 = unname(test$estimate),
      delta_rho = unname(test$estimate) - cor_results$spearman_rho[i],
      p_without_174_approx = test$p.value,
      stringsAsFactors = FALSE
    )
  }
  sensitivity <- do.call(rbind, sensitivity)
  sensitivity$BH_FDR_across_5_tests <- p.adjust(
    sensitivity$p_without_174_approx, method = "BH")
  write.csv(sensitivity,
            file.path(DIR_TABLE, "sensitivity_excluding_frozen_174_module.csv"),
            row.names = FALSE)
  message("\n174-gene-excluded sensitivity:")
  print(sensitivity[, c("comparison", "n_remaining", "rho_without_174")],
        row.names = FALSE)
}

# ----------------------------- FIGURE 6 -------------------------------------
# A and B: full shared gene universe shown as 2-D binned scatter (not DEG-only).
# C: three cross-species rho values; no confidence intervals are invented.
# Per-gene Spearman rho and P calculated on all points BEFORE clipping for display.

format_p <- function(p) if (is.na(p)) "NA" else if (p < 2.2e-16) {
  "<2.2e-16"
} else {
  format(p, digits = 2, scientific = TRUE)
}

make_density_panel <- function(comparison, title, x_lab, y_lab) {
  d <- all_pairs[all_pairs$comparison == comparison, , drop = FALSE]
  r <- cor_results[cor_results$comparison == comparison, , drop = FALSE]
  axis_limit <- stats::quantile(abs(c(d$log2FC_first, d$log2FC_second)),
                                PLOT_QUANTILE, na.rm = TRUE, names = FALSE)
  axis_limit <- max(axis_limit, 0.5)
  # Plotting-only truncation of outliers is explicit in the figure legend.
  plotdat <- d[abs(d$log2FC_first) <= axis_limit &
               abs(d$log2FC_second) <= axis_limit, , drop = FALSE]
  ggplot2::ggplot(plotdat,
                  ggplot2::aes(x = .data$log2FC_first,
                               y = .data$log2FC_second)) +
    ggplot2::geom_bin_2d(bins = 48) +
    ggplot2::scale_fill_gradient(low = "#C7DCE9", high = "#245D86",
                                name = "Genes per bin") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "#828282") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "#828282") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         linewidth = 0.3, colour = "#747474") +
    ggplot2::coord_fixed(xlim = c(-axis_limit, axis_limit),
                         ylim = c(-axis_limit, axis_limit), expand = FALSE) +
    ggplot2::labs(title = title,
      subtitle = sprintf("rho = %.3f  |  n = %s genes  |  P %s",
                         r$spearman_rho, format(r$n_shared, big.mark = ","),
                         ifelse(r$spearman_p_approx < 2.2e-16, "< 2.2e-16",
                                paste("=", format_p(r$spearman_p_approx)))),
      x = x_lab, y = y_lab) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12.5),
      axis.title = ggplot2::element_text(size = 13.2),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11.5),
      panel.grid.minor = ggplot2::element_blank()
    )
}

figA <- make_density_panel(
  "HUMAN_ULLS_vs_DI_VL", "A  Independent human VL unloading paradigms",
  "OSD-540 ULLS vs baseline (log2FC)",
  "GSE271607 dry immersion vs before (log2FC)")
figB <- make_density_panel(
  "HUMAN_DI_VL_vs_DI_SOLEUS", "B  Human dry immersion: two muscles",
  "Vastus lateralis, after vs before (log2FC)",
  "Soleus, after vs before (log2FC)")

m <- cor_results[cor_results$group == "Mouse vs human", , drop = FALSE]
m$label <- factor(
  c("vs human ULLS (VL)", "vs human dry immersion (VL)",
    "vs human dry immersion (soleus)"),
  levels = rev(c("vs human ULLS (VL)", "vs human dry immersion (VL)",
                 "vs human dry immersion (soleus)")))
m$annotation <- sprintf("rho = %+.3f  |  n = %s", m$spearman_rho,
                        format(m$n_shared, big.mark = ","))
lo <- min(-0.35, min(m$spearman_rho) - 0.13)
hi <- max(1.2, max(m$spearman_rho) + 0.7)  # room for text in a separate right column
annotation_x <- hi - 0.03
figC <- ggplot2::ggplot(m, ggplot2::aes(y = .data$label,
                                        x = .data$spearman_rho)) +
  ggplot2::geom_vline(xintercept = 0, colour = "#8F969D", linewidth = 0.65) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$spearman_rho,
                                    yend = .data$label),
                        linewidth = 1.2, colour = "#B7C4D0") +
  ggplot2::geom_point(size = 4.7, shape = 21, fill = "#357DA5",
                      colour = "#1F4965", stroke = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = .data$annotation),
                     x = annotation_x, hjust = 1, size = 4.4,
                     colour = "#293847") +
  ggplot2::scale_x_continuous(limits = c(lo, hi)) +
  ggplot2::labs(
    title = "C  Mouse HLU versus human unloading responses",
    subtitle = "OSD-21 mouse gastrocnemius HLU vs control; orthologue-matched comparisons",
    x = "Spearman correlation of shared-gene log2 fold changes", y = NULL) +
  ggplot2::theme_minimal(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 17, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 13),
    axis.text.y = ggplot2::element_text(size = 14, colour = "#293847"),
    axis.text.x = ggplot2::element_text(size = 12.5),
    axis.title.x = ggplot2::element_text(size = 13.5),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 18, 12, 14)
  )

fig6 <- patchwork::wrap_plots(figA, figB, figC,
                              design = "AABB\nAABB\nCCCC\nCCCC",
                              guides = "keep") &
  ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))

save_figure <- function(plot, stem, width, height) {
  pngfile <- file.path(DIR_FIG, paste0(stem, ".png"))
  svgfile <- file.path(DIR_FIG, paste0(stem, ".svg"))
  ggplot2::ggsave(pngfile, plot, width = width, height = height,
                  units = "in", dpi = DPI, limitsize = FALSE, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(svgfile, plot, width = width, height = height,
                    units = "in", device = svglite::svglite,
                    limitsize = FALSE, bg = "white")
  } else {
    grDevices::svg(filename = svgfile, width = width, height = height,
                   bg = "white")
    print(plot)
    grDevices::dev.off()
  }
  message("  Figure: ", pngfile)
  message("          ", svgfile)
}

save_figure(fig6, "Figure_6_UNLOADING_human_and_mouse_genomewide_concordance",
            FIG_WIDTH, FIG_HEIGHT)
save_figure(figA, "Figure_6A_human_ULLS_vs_DI_VL", 7.0, 6.0)
save_figure(figB, "Figure_6B_human_DI_VL_vs_soleus", 7.0, 6.0)
save_figure(figC, "Figure_6C_mouse_HLU_vs_human_unloading", 11.0, 4.7)

# -------------------------- AUDIT TRAIL ------------------------------------
audit <- c(
  "Genome-wide mechanical unloading concordance | script 13",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste("ggplot2:", as.character(utils::packageVersion("ggplot2"))),
  paste("babelgene:", as.character(utils::packageVersion("babelgene"))),
  paste("One-to-one among mapped observed mouse genes:",
        REQUIRE_UNIQUE_MOUSE_TO_HUMAN),
  paste("Orthologue min support:", ORTHOLOGY_MIN_SUPPORT),
  "Mouse OSD-21 comparisons use HLU_vs_CTRL (NOT actual flight).",
  "All input comparisons use log2FC in experimental condition minus reference.",
  "All correlations used shared evaluable genes without DEG filtering.",
  "P-values: asymptotic Spearman tests treating genes as observations; descriptive only.",
  "Figure A/B bins display the central 99.5% of |log2FC| values for readability,",
  "but statistical results use the full, untrimmed shared-gene universe.",
  "Cross-species effects mix microarray/RNA-seq platforms; interpret cautiously.",
  if (is.null(module_genes)) "174 module exclusion NOT RUN: list absent."
  else "174 module exclusion RUN using the frozen human-symbol list.",
  "Original human-human correlations are recorded and checked in QC.",
  "The script does NOT redefine or select the ALS-associated module.",
  "Input files:", paste(names(DE_PATHS), DE_PATHS, sep = " = "),
  "Session info:", capture.output(utils::sessionInfo())
)
writeLines(audit, file.path(OUT_DIR, "METHODS_AND_QC_AUDIT_TRAIL.txt"))

message("\nDone. Tables: ", DIR_TABLE, "\nFigure 6: ", DIR_FIG,
        "\nQC: ", DIR_QC)
