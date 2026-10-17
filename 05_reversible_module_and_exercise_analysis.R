# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

###############################################################################
# Reversible ALS–unloading module and human ALS exercise analysis
#
# This script defines the directional module from the SOD1G93A disease
# trajectory and independent unloading/recovery studies. It evaluates
# disease-stage robustness and tests the module in the paired human resistance
# exercise dataset GSE250455.
#
# GSE250455 is evaluated after module definition and does not affect module
# membership. Required data and prior analysis outputs are configured below.
###############################################################################
# =============================================================================
# 0. USER SETTINGS
# =============================================================================

PROJECT_DIR <- "."

# ---- Previous ALS/spaceflight analysis outputs ----
# Leave as NA to auto-detect one of the known output folders below.
PREVIOUS_RESULTS_DIR <- NA_character_

PREVIOUS_RESULTS_CANDIDATES <- c(
  file.path(PROJECT_DIR, "NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL"),
  file.path(PROJECT_DIR, "NASA_ALS_metabolic_pathway_results_ALS_DERIVED_HUMAN_BIOPSIES"),
  file.path(PROJECT_DIR, "NASA_ALS_metabolic_pathway_results_ALS_DERIVED"),
  file.path(PROJECT_DIR, "nasa_als_results")
)

# ---- SOD1 significant mRNA workbook ----
SOD1_WORKBOOK <- file.path(PROJECT_DIR, "SOD_data.xlsx")

# These are intentionally fixed because they are the mutation-specific
# hSOD1G93A comparisons in the supplied supplementary workbook.
SOD1_SHEET_VS_WT  <- "Table S3"
SOD1_SHEET_VS_NTG <- "Table S4"

# ---- Human ALS resistance-exercise dataset ----
GSE250455_FILE <- file.path(PROJECT_DIR, "GSE250455_RawCounts.txt")
if (!file.exists(GSE250455_FILE)) {
  gz_alt <- paste0(GSE250455_FILE, ".gz")
  if (file.exists(gz_alt)) GSE250455_FILE <- gz_alt
}

# ---- Output ----
OUTDIR <- file.path(PROJECT_DIR, "reversible_ALS_unloading_results")
TABLE_DIR <- file.path(OUTDIR, "tables")
FIG_DIR <- file.path(OUTDIR, "figures")
QC_DIR <- file.path(OUTDIR, "QC")

# ---- Module definition ----
# A very small log2FC can have an unstable sign. Effects smaller than this are
# treated as "no directional support" for module selection.
MIN_ABS_LFC <- 0.10

# Full module:
# - SOD1: same direction in >=2 of 3 stages
# - unloading: same direction in >=2 independent studies
# - recovery: reversed in >=1 study and NOT reinforced in another recovery study
MIN_SOD1_SUPPORT_FULL <- 2
MIN_UNLOADING_SUPPORT <- 2
MIN_RECOVERY_REVERSE <- 1

# High-confidence subset requires reversal in BOTH recovery datasets when both
# are available for that gene.
REQUIRE_BOTH_RECOVERY_FOR_HIGH_CONF <- TRUE

# Permutations for the formal convergence test.
# 5,000 is a sensible working value. Use 10,000+ for final manuscript analyses.
N_PERM <- 20000
SEED <- 20261004

# GSE250455 expression filtering.
MIN_CPM <- 1
MIN_SAMPLES_CPM <- 4

# Significance thresholds are NOT required for primary module membership,
# because replication/directionality is the central criterion. They are used
# to create an additional "statistically supported" subset.
SUPPORT_NOMINAL_P <- 0.05
SUPPORT_FDR <- 0.10

set.seed(SEED)

# =============================================================================
# 1. PACKAGES
# =============================================================================

cran_pkgs <- c(
  "tidyverse",
  "openxlsx",
  "readxl",
  "ggrepel",
  "scales",
  "babelgene"
)

bioc_pkgs <- c(
  "edgeR",
  "limma",
  "fgsea",
  "AnnotationDbi",
  "org.Mm.eg.db",
  "org.Hs.eg.db"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) == 0) return(invisible(NULL))

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing)
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(readxl)
  library(ggrepel)
  library(scales)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
})

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. GENERAL HELPERS
# =============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  toupper(x)
}

norm_id <- function(x) {
  # Uppercase FIRST, then remove punctuation/separators.
  # The previous order removed lowercase letters entirely, which caused
  # mixed-case columns such as "Wald", "p-value", "p-adjusted", "symbol",
  # and "statistic" to be misidentified.
  gsub("[^A-Z0-9]+", "", toupper(as.character(x)))
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

first_matching_col <- function(nms, patterns, required = FALSE, label = "column") {
  nn <- norm_id(nms)
  for (p in patterns) {
    hit <- which(str_detect(nn, p))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  if (required) {
    stop(
      "Could not identify ", label, ". Available columns:\n",
      paste(nms, collapse = ", ")
    )
  }
  NA_character_
}

infer_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- case_when(
    str_detect(x, "(^|[^0-9])90([^0-9]|$)|PRESYM") ~ "90",
    str_detect(x, "(^|[^0-9])120([^0-9]|$)|EARLY") ~ "120",
    str_detect(x, "(^|[^0-9])150([^0-9]|$)|LATE") ~ "150",
    TRUE ~ NA_character_
  )
  out
}

signed_direction <- function(x, min_abs = MIN_ABS_LFC) {
  # Exact zero must remain directionless even when min_abs = 0.
  # For non-zero thresholds, values exactly on the threshold retain the
  # previous inclusive >= / <= behaviour.
  case_when(
    is.na(x) | !is.finite(x) ~ 0L,
    x == 0 ~ 0L,
    x > 0 & x >= min_abs ~ 1L,
    x < 0 & x <= -min_abs ~ -1L,
    TRUE ~ 0L
  )
}

safe_binom_p <- function(k, n, alternative = "greater") {
  if (is.na(n) || n < 1) return(NA_real_)
  binom.test(k, n, p = 0.5, alternative = alternative)$p.value
}

safe_wilcox_p <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  if (length(x) < 3 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

safe_t_p <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(NA_real_)
  t.test(x, mu = 0, alternative = alternative)$p.value
}

read_delim_auto <- function(path) {
  ext <- tolower(tools::file_ext(gsub("\\.gz$", "", path)))
  if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else if (ext %in% c("txt", "tsv")) {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    stop("Unsupported text format: ", path)
  }
}

standardise_effect_table <- function(df, stage_hint = NA_character_, source = NA_character_) {
  nms <- colnames(df)

  symbol_col <- first_matching_col(
    nms,
    c("^SYMBOL$", "^GENESYMBOL$", "^GENENAME$", "^GENE$"),
    required = TRUE,
    label = "gene-symbol column"
  )

  lfc_col <- first_matching_col(
    nms,
    c("^LOG2FC$", "^LOG2FOLDCHANGE$", "^LOGFC$", "^LFC$"),
    required = TRUE,
    label = "log2-fold-change column"
  )

  stat_col <- first_matching_col(
    nms,
    c("^STAT$", "^T$", "^WALDSTAT$", "^WALD$", "^MODERATEDT$"),
    required = FALSE
  )

  p_col <- first_matching_col(
    nms,
    c("^PVALUE$", "^PVAL$", "^PVALUE", "^PVALUEUNADJUSTED$", "^P\\.VALUE$"),
    required = FALSE
  )

  fdr_col <- first_matching_col(
    nms,
    c("^FDR$", "^PADJ$", "^ADJPVALUE$", "^ADJPVAL$", "^ADJ\\.P\\.VAL$"),
    required = FALSE
  )

  stage_col <- first_matching_col(
    nms,
    c("^STAGE$", "^TIMEPOINT$", "^AGE$", "^DAY$", "^DAYS$"),
    required = FALSE
  )

  out <- tibble(
    symbol = clean_symbol(df[[symbol_col]]),
    log2FC = safe_num(df[[lfc_col]]),
    stat = if (!is.na(stat_col)) safe_num(df[[stat_col]]) else NA_real_,
    p_value = if (!is.na(p_col)) safe_num(df[[p_col]]) else NA_real_,
    FDR = if (!is.na(fdr_col)) safe_num(df[[fdr_col]]) else NA_real_,
    stage_raw = if (!is.na(stage_col)) as.character(df[[stage_col]]) else stage_hint,
    source_file = source
  ) %>%
    mutate(
      stage = infer_stage(stage_raw),
      stage = if_else(is.na(stage) & !is.na(stage_hint), stage_hint, stage)
    ) %>%
    filter(!is.na(symbol), is.finite(log2FC))

  # Resolve duplicate symbols by retaining the row with the strongest absolute
  # statistic, then smallest p-value, then largest absolute log2FC.
  out %>%
    mutate(
      abs_stat_sort = if_else(is.finite(stat), abs(stat), -Inf),
      p_sort = if_else(is.finite(p_value), p_value, Inf),
      abs_lfc_sort = abs(log2FC)
    ) %>%
    arrange(symbol, desc(abs_stat_sort), p_sort, desc(abs_lfc_sort)) %>%
    distinct(symbol, stage, .keep_all = TRUE) %>%
    dplyr::select(-abs_stat_sort, -p_sort, -abs_lfc_sort, -stage_raw)
}

read_effect_file <- function(path, stage_hint = NA_character_) {
  if (!file.exists(path)) stop("File not found: ", path)

  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv", "txt", "tsv") || str_detect(tolower(path), "\\.(csv|txt|tsv)\\.gz$")) {
    x <- read_delim_auto(path)
    return(standardise_effect_table(x, stage_hint = stage_hint, source = basename(path)))
  }

  if (ext %in% c("xlsx", "xlsm", "xls")) {
    sheets <- openxlsx::getSheetNames(path)
    parts <- list()

    for (s in sheets) {
      x <- tryCatch(
        openxlsx::read.xlsx(path, sheet = s, check.names = FALSE),
        error = function(e) NULL
      )
      if (is.null(x) || nrow(x) == 0) next

      sh <- stage_hint
      if (is.na(sh)) sh <- infer_stage(s)

      z <- tryCatch(
        standardise_effect_table(x, stage_hint = sh, source = paste0(basename(path), "::", s)),
        error = function(e) NULL
      )

      if (!is.null(z) && nrow(z) > 0) parts[[s]] <- z
    }

    if (length(parts) == 0) {
      stop("Could not identify usable gene-level DE sheets in: ", path)
    }

    return(bind_rows(parts))
  }

  stop("Unsupported effect-file type: ", path)
}

# =============================================================================
# 3. LOAD SOD1 SIGNIFICANT mRNAs DIRECTLY FROM TABLES S3 + S4
# =============================================================================

signed_fold_change_to_log2 <- function(fc) {
  # The supplementary tables use a signed fold-change convention:
  #   +2 means 2-fold up; -2 means 2-fold down.
  # All reported significant values in Tables S3/S4 have |FC| >= 2.
  #
  # Convert magnitude to log2 scale while preserving the explicitly reported
  # sign. abs(log2(abs(fc))) is used so the sign can never be accidentally
  # inverted if a future table contains a magnitude below 1.
  fc <- safe_num(fc)

  out <- rep(NA_real_, length(fc))
  ok <- is.finite(fc) & fc != 0

  out[ok] <- sign(fc[ok]) * abs(log2(abs(fc[ok])))
  out
}

read_sod1_sheet_stage <- function(path, sheet, stage) {
  # Table S3 / S4 layout:
  #   columns 1:4   = 90 d
  #   columns 6:9   = 120 d
  #   columns 11:14 = 150 d
  # Row 4 contains headers; row 5 onward contains significant genes.
  #
  # This parser deliberately uses BASE-R column access and naming. It does not
  # use dplyr::rename()/transmute() for p-value/FDR fields, avoiding the
  # non-standard-evaluation error "object 'p_value' not found".

  col_idx <- switch(
    as.character(stage),
    `90`  = 1:4,
    `120` = 6:9,
    `150` = 11:14,
    stop("Unsupported SOD1 stage: ", stage)
  )

  raw <- readxl::read_excel(
    path,
    sheet = sheet,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  raw <- as.data.frame(raw, check.names = FALSE)

  if (nrow(raw) < 5 || ncol(raw) < max(col_idx)) {
    stop("Unexpected layout in ", sheet, " for ", stage, " d.")
  }

  dat <- raw[5:nrow(raw), col_idx, drop = FALSE]

  out <- data.frame(
    ensembl_id  = sub("\\..*$", "", trimws(as.character(dat[[1]]))),
    p           = safe_num(dat[[2]]),
    fdr         = safe_num(dat[[3]]),
    fold_change = safe_num(dat[[4]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out$log2FC <- signed_fold_change_to_log2(out$fold_change)

  keep <- !is.na(out$ensembl_id) &
          out$ensembl_id != "" &
          is.finite(out$log2FC)

  out <- out[keep, , drop = FALSE]
  out <- out[!duplicated(out$ensembl_id), , drop = FALSE]

  tibble::as_tibble(out)
}


map_mouse_ensembl_symbols <- function(ids) {
  ids <- unique(sub("\\..*$", "", as.character(ids)))

  syms <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )

  data.frame(
    ensembl_id = ids,
    symbol = clean_symbol(unname(syms[ids])),
    stringsAsFactors = FALSE
  )
}


load_sod1_stage_effects <- function() {
  if (!file.exists(SOD1_WORKBOOK)) {
    stop(
      "SOD1 supplementary workbook not found: ", SOD1_WORKBOOK, "\n",
      "Place SOD_data.xlsx in PROJECT_DIR or update SOD1_WORKBOOK."
    )
  }

  available_sheets <- readxl::excel_sheets(SOD1_WORKBOOK)
  needed <- c(SOD1_SHEET_VS_WT, SOD1_SHEET_VS_NTG)
  missing_sheets <- setdiff(needed, available_sheets)

  if (length(missing_sheets) > 0) {
    stop(
      "Required SOD1 sheet(s) missing from workbook: ",
      paste(missing_sheets, collapse = ", ")
    )
  }

  stages <- c("90", "120", "150")
  expected_counts <- c(`90` = 558L, `120` = 544L, `150` = 685L)
  out <- vector("list", length(stages))
  names(out) <- stages

  for (st in stages) {
    message("Reading SOD1 ", st, " d significant mRNAs from Tables S3 + S4...")

    vs_wt <- as.data.frame(
      read_sod1_sheet_stage(SOD1_WORKBOOK, SOD1_SHEET_VS_WT, st),
      stringsAsFactors = FALSE
    )
    names(vs_wt) <- c(
      "ensembl_id", "p_WT", "FDR_WT", "fold_change_WT", "log2FC_WT"
    )

    vs_ntg <- as.data.frame(
      read_sod1_sheet_stage(SOD1_WORKBOOK, SOD1_SHEET_VS_NTG, st),
      stringsAsFactors = FALSE
    )
    names(vs_ntg) <- c(
      "ensembl_id", "p_NTg", "FDR_NTg", "fold_change_NTg", "log2FC_NTg"
    )

    message(
      "  Table S3 rows: ", nrow(vs_wt),
      "; Table S4 rows: ", nrow(vs_ntg)
    )

    # Inner join = significant in BOTH hSOD1G93A-vs-hSOD1WT and
    # hSOD1G93A-vs-NTg.
    both <- merge(
      vs_wt,
      vs_ntg,
      by = "ensembl_id",
      all = FALSE,
      sort = FALSE
    )

    both$direction_WT  <- sign(both$log2FC_WT)
    both$direction_NTg <- sign(both$log2FC_NTg)

    n_opposite <- sum(
      both$direction_WT != both$direction_NTg,
      na.rm = TRUE
    )

    if (n_opposite > 0) {
      warning(
        "At ", st, " d, ", n_opposite,
        " genes were significant against both controls but had opposite directions. ",
        "They will be excluded."
      )
    }

    same_direction <- !is.na(both$direction_WT) &
                      !is.na(both$direction_NTg) &
                      both$direction_WT == both$direction_NTg

    both <- both[same_direction, , drop = FALSE]

    message(
      "  Significant in BOTH comparisons, same direction: ", nrow(both)
    )

    if (nrow(both) != expected_counts[[st]]) {
      warning(
        "Expected ", expected_counts[[st]],
        " genes significant in BOTH SOD1 comparisons at ", st,
        " d, but extracted ", nrow(both), ". Check workbook layout/version."
      )
    }

    # Both control comparisons have the same direction. Use their mean log2
    # effect as the mutation-specific SOD1 effect.
    both$log2FC <- rowMeans(
      cbind(both$log2FC_WT, both$log2FC_NTg),
      na.rm = TRUE
    )

    # Conservative annotations: use the less significant result of the two
    # comparisons. Module membership does not depend on these because the source
    # sheets already contain significant genes.
    both$p_value <- pmax(both$p_WT, both$p_NTg, na.rm = TRUE)
    both$FDR <- pmax(both$FDR_WT, both$FDR_NTg, na.rm = TRUE)

    mapping <- map_mouse_ensembl_symbols(both$ensembl_id)

    stage_pre <- merge(
      both,
      mapping,
      by = "ensembl_id",
      all.x = TRUE,
      sort = FALSE
    )

    stage_pre <- stage_pre[
      !is.na(stage_pre$symbol) & stage_pre$symbol != "",
      ,
      drop = FALSE
    ]

    n_unmapped <- nrow(both) - nrow(stage_pre)
    if (n_unmapped > 0) {
      warning(
        "At ", st, " d, ", n_unmapped,
        " mutation-specific Ensembl IDs did not map to an official gene symbol."
      )
    }

    # Collapse rare duplicate Ensembl-to-symbol mappings using explicit base-R
    # indexing; no bare p_value/FDR symbols are evaluated anywhere in this block.
    by_symbol <- split(stage_pre, stage_pre$symbol)

    collapse_one_symbol <- function(z) {
      data.frame(
        symbol = z$symbol[[1]],
        ensembl_id = paste(sort(unique(z$ensembl_id)), collapse = ";"),
        stage = st,
        log2FC = mean(z$log2FC, na.rm = TRUE),
        p_value = max(z$p_value, na.rm = TRUE),
        FDR = max(z$FDR, na.rm = TRUE),
        fold_change_WT = mean(z$fold_change_WT, na.rm = TRUE),
        fold_change_NTg = mean(z$fold_change_NTg, na.rm = TRUE),
        log2FC_WT = mean(z$log2FC_WT, na.rm = TRUE),
        log2FC_NTg = mean(z$log2FC_NTg, na.rm = TRUE),
        p_WT = max(z$p_WT, na.rm = TRUE),
        p_NTg = max(z$p_NTg, na.rm = TRUE),
        FDR_WT = max(z$FDR_WT, na.rm = TRUE),
        FDR_NTg = max(z$FDR_NTg, na.rm = TRUE),
        dataset = "SOD1G93A",
        source_file = basename(SOD1_WORKBOOK),
        stringsAsFactors = FALSE
      )
    }

    stage_dat <- do.call(
      rbind,
      lapply(by_symbol, collapse_one_symbol)
    )
    rownames(stage_dat) <- NULL

    n_symbol_collapses <- nrow(stage_pre) - nrow(stage_dat)
    if (n_symbol_collapses > 0) {
      message(
        "  Collapsed ", n_symbol_collapses,
        " duplicate Ensembl-to-symbol mapping row(s)."
      )
    }

    stage_dat <- stage_dat[, c(
      "symbol",
      "ensembl_id",
      "stage",
      "log2FC",
      "p_value",
      "FDR",
      "fold_change_WT",
      "fold_change_NTg",
      "log2FC_WT",
      "log2FC_NTg",
      "p_WT",
      "p_NTg",
      "FDR_WT",
      "FDR_NTg",
      "dataset",
      "source_file"
    ), drop = FALSE]

    stage_dat <- stage_dat[order(stage_dat$symbol), , drop = FALSE]

    message("  Official gene symbols retained: ", nrow(stage_dat))
    out[[st]] <- tibble::as_tibble(stage_dat)
  }

  ans <- dplyr::bind_rows(out)

  message("\nSOD1 mutation-specific significant sets extracted:")
  print(
    ans %>%
      dplyr::count(stage, name = "n_mapped_unique_symbols") %>%
      dplyr::arrange(factor(stage, levels = c("90", "120", "150")))
  )

  ans
}


# =============================================================================
# 3B. CROSS-SPECIES GENE HARMONISATION
# =============================================================================
#
# All cross-dataset analyses are carried out in HUMAN ORTHOLOG GENE-SYMBOL
# SPACE:
#   * SOD1G93A mouse genes -> human ortholog symbols
#   * OSD-21 mouse genes   -> human ortholog symbols
#   * OSD-540 / GSE271607 human Ensembl IDs -> human gene symbols
#   * GSE250455 already provides human gene symbols
#
# Mouse-to-human orthologs are obtained with babelgene using its highest-
# supported match per mouse gene (top = TRUE) and at least three supporting
# orthology resources (min_support = 3). This avoids relying on capitalization
# or assuming that identical mouse/human symbols are always orthologous.

mouse_to_human_ortholog_map <- function(mouse_symbols) {
  # The rest of this workflow stores gene symbols in upper case for consistent
  # joins across datasets. babelgene, however, matches non-human gene SYMBOLS
  # case-sensitively (e.g. mouse "Pgm1", not "PGM1"). Therefore first restore
  # official mouse symbol case using org.Mm.eg.db, then query babelgene.
  mouse_symbols_input <- unique(clean_symbol(mouse_symbols))
  mouse_symbols_input <- mouse_symbols_input[!is.na(mouse_symbols_input)]

  if (length(mouse_symbols_input) == 0) {
    return(tibble(
      mouse_symbol = character(),
      mouse_symbol_official = character(),
      human_symbol = character(),
      ortholog_support_n = integer(),
      ortholog_support = character()
    ))
  }

  # Build a case-insensitive lookup from the official mouse SYMBOL vocabulary.
  mm_official <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
  mm_lookup <- setNames(mm_official, toupper(mm_official))

  mouse_symbol_official <- unname(mm_lookup[mouse_symbols_input])

  input_map <- tibble(
    mouse_symbol = mouse_symbols_input,
    mouse_symbol_official = mouse_symbol_official
  ) %>%
    filter(!is.na(mouse_symbol_official))

  if (nrow(input_map) == 0) {
    stop(
      "None of the supplied mouse symbols could be restored to official ",
      "mouse SYMBOL case using org.Mm.eg.db."
    )
  }

  ortho <- babelgene::orthologs(
    genes = unique(input_map$mouse_symbol_official),
    species = "mouse",
    human = FALSE,
    min_support = 3,
    top = TRUE
  )

  if (nrow(ortho) == 0) {
    stop("babelgene returned no mouse-to-human ortholog mappings.")
  }

  # babelgene returns columns including:
  # human_symbol, symbol (mouse), support, support_n.
  ortho_small <- tibble(
    mouse_symbol_official = as.character(ortho$symbol),
    human_symbol = clean_symbol(ortho$human_symbol),
    ortholog_support_n = safe_num(ortho$support_n),
    ortholog_support = as.character(ortho$support)
  ) %>%
    filter(
      !is.na(mouse_symbol_official),
      !is.na(human_symbol)
    ) %>%
    arrange(mouse_symbol_official, desc(ortholog_support_n)) %>%
    distinct(mouse_symbol_official, .keep_all = TRUE)

  out <- input_map %>%
    left_join(
      ortho_small,
      by = "mouse_symbol_official"
    ) %>%
    filter(!is.na(human_symbol)) %>%
    arrange(mouse_symbol, desc(ortholog_support_n)) %>%
    distinct(mouse_symbol, .keep_all = TRUE)

  message(
    "    Mouse ortholog mapping: ",
    length(mouse_symbols_input), " input symbols; ",
    nrow(input_map), " recognized as official mouse symbols; ",
    nrow(out), " mapped to human orthologs."
  )

  out
}


human_ensembl_to_symbol_map <- function(ids) {
  ids_raw <- as.character(ids)
  ids_clean <- sub("\\..*$", "", ids_raw)
  ids_clean <- unique(ids_clean[grepl("^ENSG[0-9]+$", toupper(ids_clean))])

  if (length(ids_clean) == 0) {
    return(tibble(
      human_ensembl = character(),
      human_symbol = character()
    ))
  }

  syms <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ids_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )

  tibble(
    human_ensembl = ids_clean,
    human_symbol = clean_symbol(unname(syms[ids_clean]))
  ) %>%
    filter(!is.na(human_symbol))
}


collapse_human_gene_effects <- function(x) {
  # Multiple source identifiers can occasionally resolve to the same human gene
  # symbol. Collapse them to one row per human gene. If mapped source genes have
  # conflicting effect directions, remove that ambiguous human symbol rather
  # than arbitrarily choosing one effect.
  x %>%
    filter(!is.na(symbol), symbol != "", is.finite(log2FC)) %>%
    mutate(.effect_sign = sign(log2FC)) %>%
    group_by(symbol) %>%
    filter(n_distinct(.effect_sign[.effect_sign != 0]) <= 1L) %>%
    summarise(
      log2FC = median(log2FC, na.rm = TRUE),
      stat = if (all(is.na(stat))) NA_real_ else median(stat, na.rm = TRUE),
      p_value = if (all(is.na(p_value))) NA_real_ else max(p_value, na.rm = TRUE),
      FDR = if (all(is.na(FDR))) NA_real_ else max(FDR, na.rm = TRUE),
      AveExpr = if (all(is.na(AveExpr))) NA_real_ else median(AveExpr, na.rm = TRUE),
      source_id = paste(sort(unique(source_id)), collapse = ";"),
      source_file = dplyr::first(source_file),
      .groups = "drop"
    )
}


harmonize_human_de_to_symbol <- function(x) {
  raw_id <- as.character(x$symbol)
  stripped <- sub("\\..*$", "", raw_id)
  is_ensg <- grepl("^ENSG[0-9]+$", toupper(stripped))

  # Start by treating non-Ensembl identifiers as already being gene symbols.
  mapped_symbol <- clean_symbol(raw_id)

  if (any(is_ensg)) {
    mp <- human_ensembl_to_symbol_map(stripped[is_ensg])
    mp_vec <- mp$human_symbol
    names(mp_vec) <- mp$human_ensembl
    mapped_symbol[is_ensg] <- unname(mp_vec[stripped[is_ensg]])
  }

  n_input <- nrow(x)
  n_ensg <- sum(is_ensg)
  n_mapped <- sum(!is.na(mapped_symbol) & mapped_symbol != "")

  message(
    "    Human ID harmonisation: ", n_input, " rows; ",
    n_ensg, " Ensembl IDs; ",
    n_mapped, " rows mapped to human symbols."
  )

  x %>%
    mutate(
      source_id = raw_id,
      symbol = mapped_symbol
    ) %>%
    collapse_human_gene_effects()
}


harmonize_mouse_de_to_human <- function(x) {
  if (!"symbol" %in% names(x)) {
    stop(
      "Mouse DE table has no 'symbol' column. Available columns: ",
      paste(names(x), collapse = ", ")
    )
  }

  source_mouse <- clean_symbol(x$symbol)
  mp <- mouse_to_human_ortholog_map(source_mouse)

  out <- x %>%
    dplyr::mutate(source_id = source_mouse) %>%
    dplyr::left_join(
      mp %>% dplyr::select(
        mouse_symbol,
        human_symbol,
        ortholog_support_n,
        ortholog_support
      ),
      by = c("source_id" = "mouse_symbol")
    ) %>%
    dplyr::mutate(symbol = human_symbol)

  n_mapped <- sum(!is.na(out$symbol))

  message(
    "    Mouse-to-human ortholog harmonisation: ",
    nrow(x), " mouse rows; ",
    n_mapped, " mapped rows."
  )

  if (n_mapped == 0) {
    stop(
      "Mouse-to-human ortholog harmonisation produced zero mapped rows. ",
      "Check mouse gene-symbol formatting."
    )
  }

  out %>%
    dplyr::select(
      symbol,
      log2FC,
      stat,
      p_value,
      FDR,
      AveExpr,
      source_id,
      source_file
    ) %>%
    collapse_human_gene_effects()
}


harmonize_sod1_to_human <- function(sod1_mouse) {
  if (!"symbol" %in% names(sod1_mouse)) {
    stop(
      "SOD1 table has no 'symbol' column. Available columns: ",
      paste(names(sod1_mouse), collapse = ", ")
    )
  }

  mp <- mouse_to_human_ortholog_map(sod1_mouse$symbol)

  # Avoid rename() entirely here because Bioconductor packages can mask tidyverse
  # generics in some R sessions. Construct the mouse-symbol field explicitly.
  x0 <- as.data.frame(sod1_mouse, stringsAsFactors = FALSE)
  x0$mouse_symbol <- clean_symbol(x0$symbol)
  x0$symbol <- NULL
  x0 <- tibble::as_tibble(x0)

  x <- x0 %>%
    dplyr::left_join(
      mp %>%
        dplyr::select(
          mouse_symbol,
          mouse_symbol_official,
          human_symbol,
          ortholog_support_n,
          ortholog_support
        ),
      by = "mouse_symbol"
    ) %>%
    dplyr::filter(!is.na(human_symbol)) %>%
    dplyr::mutate(
      symbol = human_symbol,
      .effect_sign = sign(log2FC)
    )

  # Resolve rare many-mouse-to-one-human cases conservatively within each SOD1
  # stage. If mapped mouse genes disagree in direction, drop that human symbol
  # for that stage. Otherwise average the mapped effect and use the larger
  # p/FDR as conservative annotations.
  out <- x %>%
    dplyr::group_by(stage, symbol) %>%
    dplyr::filter(
      dplyr::n_distinct(.effect_sign[.effect_sign != 0]) <= 1L
    ) %>%
    dplyr::summarise(
      mouse_symbol = paste(sort(unique(mouse_symbol)), collapse = ";"),
      mouse_symbol_official =
        paste(sort(unique(mouse_symbol_official)), collapse = ";"),
      ensembl_id = paste(sort(unique(ensembl_id)), collapse = ";"),
      log2FC = mean(log2FC, na.rm = TRUE),
      p_value = if (all(is.na(p_value))) NA_real_ else max(p_value, na.rm = TRUE),
      FDR = if (all(is.na(FDR))) NA_real_ else max(FDR, na.rm = TRUE),
      n_mouse_source_genes = dplyr::n(),
      ortholog_support_n =
        if (all(is.na(ortholog_support_n))) NA_real_
        else max(ortholog_support_n, na.rm = TRUE),
      dataset = "SOD1G93A",
      source_file = dplyr::first(source_file),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      factor(stage, levels = c("90", "120", "150")),
      symbol
    )

  message("\nSOD1 mouse-to-human ortholog mapping:")
  qc <- sod1_mouse %>%
    dplyr::count(stage, name = "n_mouse_genes") %>%
    dplyr::left_join(
      out %>% dplyr::count(stage, name = "n_human_ortholog_genes"),
      by = "stage"
    ) %>%
    dplyr::mutate(
      n_human_ortholog_genes =
        tidyr::replace_na(n_human_ortholog_genes, 0L),
      mapping_percent = 100 * n_human_ortholog_genes / n_mouse_genes
    )

  print(qc)

  write_csv(
    qc,
    file.path(QC_DIR, "SOD1_mouse_to_human_ortholog_mapping_summary.csv")
  )

  write_csv(
    mp,
    file.path(QC_DIR, "mouse_to_human_ortholog_map_SOD1.csv")
  )

  out
}


# =============================================================================
# 4. LOAD THE EXACT UNLOADING / RECOVERY GENE-LEVEL DE FILES
# =============================================================================
#
# Do NOT try to infer these comparisons from dataset/contrast text stored inside
# arbitrary CSV files. The previous analyses already produced the exact
# gene-level files we need, and OSD-21 and the human datasets live in DIFFERENT
# result directories. We therefore locate the files by their known basenames.
#
# Required comparisons:
#   OSD-21 HLU:          OSD21_HLU_vs_CTRL_Gene_symbol_level.csv
#   OSD-21 reloading:    OSD21_HLU_RL_vs_HLU_Gene_symbol_level.csv
#   OSD-540 ULLS:        OSD540_ULLS_vs_Baseline_Gene_symbol_level.csv
#   OSD-540 recovery:    OSD540_AR_vs_ULLS_Gene_symbol_level.csv
#   GSE271607 VL:        GSE271607_VL_After_vs_Before_Gene_symbol_level.csv
#   GSE271607 soleus:    GSE271607_Soleus_After_vs_Before_Gene_symbol_level.csv
#
# GSE271607 VL + soleus are from the same dry-immersion study/cohort, so they
# are combined into ONE study-level effect below and are not counted as two
# independent unloading replications.

DE_FILE_BASENAMES <- c(
  OSD21_unloading =
    "OSD21_HLU_vs_CTRL_Gene_symbol_level.csv",
  OSD21_recovery =
    "OSD21_HLU_RL_vs_HLU_Gene_symbol_level.csv",
  OSD540_unloading =
    "OSD540_ULLS_vs_Baseline_Gene_symbol_level.csv",
  OSD540_recovery =
    "OSD540_AR_vs_ULLS_Gene_symbol_level.csv",
  GSE271607_VL =
    "GSE271607_VL_After_vs_Before_Gene_symbol_level.csv",
  GSE271607_Soleus =
    "GSE271607_Soleus_After_vs_Before_Gene_symbol_level.csv"
)


candidate_de_roots <- function() {
  roots <- character()

  # Explicitly supplied previous-results folder gets first priority.
  if (!is.na(PREVIOUS_RESULTS_DIR) && nzchar(PREVIOUS_RESULTS_DIR)) {
    roots <- c(roots, PREVIOUS_RESULTS_DIR)
  }

  # Then search each known result directory separately. This is important:
  # OSD-21 normally lives in nasa_als_results, whereas the human unloading
  # outputs live in the human-biopsies pathway-analysis directory.
  roots <- c(
    roots,
    PREVIOUS_RESULTS_CANDIDATES,
    PROJECT_DIR
  )

  roots <- roots[dir.exists(roots)]

  if (length(roots) == 0) {
    stop(
      "None of the expected previous-analysis directories exists.\n",
      "Expected one or more of:\n",
      paste(PREVIOUS_RESULTS_CANDIDATES, collapse = "\n")
    )
  }

  # Keep the stated priority while removing duplicate paths.
  roots_norm <- normalizePath(roots, winslash = "/", mustWork = FALSE)
  roots[!duplicated(roots_norm)]
}


find_exact_de_file <- function(target_basename) {
  roots <- candidate_de_roots()

  for (root in roots) {
    # Fast direct checks first.
    direct_candidates <- c(
      file.path(root, target_basename),
      file.path(root, "DE", target_basename)
    )

    direct_hit <- direct_candidates[file.exists(direct_candidates)]
    if (length(direct_hit) > 0) {
      return(normalizePath(direct_hit[[1]], winslash = "/", mustWork = TRUE))
    }

    # Then recursively search this result directory. Match the BASENAME exactly,
    # case-insensitively, rather than trying to infer the biological contrast.
    csvs <- list.files(
      root,
      pattern = "\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(csvs) > 0) {
      hit <- csvs[
        tolower(basename(csvs)) == tolower(target_basename)
      ]

      if (length(hit) > 0) {
        # If duplicate copies exist in a result tree, prefer the shortest path.
        hit <- hit[order(nchar(hit))]
        return(normalizePath(hit[[1]], winslash = "/", mustWork = TRUE))
      }
    }
  }

  stop(
    "\nCould not find required DE file:\n  ", target_basename, "\n\n",
    "Searched these directories:\n  ",
    paste(candidate_de_roots(), collapse = "\n  "),
    "\n\nThis analysis uses the exact previously generated gene-level DE files."
  )
}


read_gene_level_de_file <- function(path) {
  x <- readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )

  nms <- names(x)
  nms_norm <- norm_id(nms)

  # Direct normalized-name lookup. This intentionally does NOT use
  # first_matching_col(), because that helper is designed around regex matching
  # and was normalizing the regex strings themselves (e.g. "^SYMBOL$"), causing
  # exact columns such as "symbol" not to be found.
  pick_col <- function(candidates, required = FALSE, label = "column") {
    candidates_norm <- norm_id(candidates)
    idx <- match(candidates_norm, nms_norm)
    idx <- idx[!is.na(idx)]

    if (length(idx) == 0) {
      if (required) {
        stop(
          "Could not identify ", label, " in ", basename(path), ".\n",
          "Available columns: ", paste(nms, collapse = ", ")
        )
      }
      return(NULL)
    }

    nms[[idx[[1]]]]
  }

  symbol_col <- pick_col(
    c("symbol", "gene_symbol", "genesymbol", "gene"),
    required = TRUE,
    label = "gene-symbol column"
  )

  lfc_col <- pick_col(
    c("log2FC", "logFC", "log2FoldChange"),
    required = TRUE,
    label = "log2FC column"
  )

  stat_col <- pick_col(
    c("statistic", "stat", "t", "WaldStatistic", "QLF", "F"),
    required = FALSE
  )

  p_col <- pick_col(
    c("p_value", "pvalue", "p_val", "pval", "P.Value"),
    required = FALSE
  )

  fdr_col <- pick_col(
    c("FDR", "padj", "adj.P.Val", "adj_p_value", "adjusted_p_value"),
    required = FALSE
  )

  ave_col <- pick_col(
    c("AveExpr", "baseMean", "logCPM"),
    required = FALSE
  )

  symbol_vec <- clean_symbol(x[[symbol_col]])
  lfc_vec <- safe_num(x[[lfc_col]])

  stat_vec <- if (!is.null(stat_col)) {
    safe_num(x[[stat_col]])
  } else {
    rep(NA_real_, nrow(x))
  }

  p_vec <- if (!is.null(p_col)) {
    safe_num(x[[p_col]])
  } else {
    rep(NA_real_, nrow(x))
  }

  fdr_vec <- if (!is.null(fdr_col)) {
    safe_num(x[[fdr_col]])
  } else {
    rep(NA_real_, nrow(x))
  }

  ave_vec <- if (!is.null(ave_col)) {
    v <- safe_num(x[[ave_col]])

    if (norm_id(ave_col) == "BASEMEAN") {
      v <- log2(v + 1)
    }
    v
  } else {
    rep(NA_real_, nrow(x))
  }

  out <- tibble(
    symbol = symbol_vec,
    log2FC = lfc_vec,
    stat = stat_vec,
    p_value = p_vec,
    FDR = fdr_vec,
    AveExpr = ave_vec,
    source_file = normalizePath(path, winslash = "/", mustWork = TRUE)
  ) %>%
    filter(
      !is.na(symbol),
      symbol != "",
      is.finite(log2FC)
    )

  # Gene-symbol files should already contain one row per gene. If duplicates are
  # present, collapse them conservatively rather than double-counting.
  out %>%
    group_by(symbol) %>%
    summarise(
      log2FC = median(log2FC, na.rm = TRUE),
      stat = if (all(is.na(stat))) NA_real_ else median(stat, na.rm = TRUE),
      p_value = if (all(is.na(p_value))) NA_real_ else max(p_value, na.rm = TRUE),
      FDR = if (all(is.na(FDR))) NA_real_ else max(FDR, na.rm = TRUE),
      AveExpr = if (all(is.na(AveExpr))) NA_real_ else median(AveExpr, na.rm = TRUE),
      source_file = dplyr::first(source_file),
      .groups = "drop"
    )
}


load_one_loading_effect <- function(
  path,
  study_id,
  role,
  source_species = c("mouse", "human")
) {
  source_species <- match.arg(source_species)

  message(
    "  Harmonising ", study_id, " / ", role,
    " (", source_species, ")..."
  )

  x <- read_gene_level_de_file(path)

  x <- if (source_species == "mouse") {
    harmonize_mouse_de_to_human(x)
  } else {
    harmonize_human_de_to_symbol(x)
  }

  x %>%
    mutate(
      study_id = study_id,
      role = role,
      source_species = source_species
    )
}


combine_gse271607_tissues <- function(vl_path, soleus_path) {
  message("  Harmonising GSE271607 vastus lateralis (human)...")
  vl <- read_gene_level_de_file(vl_path) %>%
    harmonize_human_de_to_symbol() %>%
    mutate(tissue = "Vastus lateralis")

  message("  Harmonising GSE271607 soleus (human)...")
  sol <- read_gene_level_de_file(soleus_path) %>%
    harmonize_human_de_to_symbol() %>%
    mutate(tissue = "Soleus")

  bind_rows(vl, sol) %>%
    group_by(symbol) %>%
    summarise(
      # One study-level effect: median of the two muscle-specific effects.
      log2FC = median(log2FC, na.rm = TRUE),
      stat = if (all(is.na(stat))) NA_real_ else median(stat, na.rm = TRUE),

      # Conservative support annotation because the two tissues are NOT
      # independent biological studies.
      p_value = if (all(is.na(p_value))) NA_real_ else max(p_value, na.rm = TRUE),
      FDR = if (all(is.na(FDR))) NA_real_ else max(FDR, na.rm = TRUE),

      AveExpr = if (all(is.na(AveExpr))) NA_real_ else median(AveExpr, na.rm = TRUE),
      n_tissues_available = n_distinct(tissue),
      source_id = paste(sort(unique(source_id)), collapse = ";"),
      source_file = paste(sort(unique(source_file)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      study_id = "GSE271607",
      role = "unloading",
      source_species = "human"
    )
}


load_unloading_recovery_effects <- function() {
  paths <- vapply(
    DE_FILE_BASENAMES,
    find_exact_de_file,
    FUN.VALUE = character(1)
  )

  message("\nExact gene-level DE files used:")
  for (nm in names(paths)) {
    message("  ", nm, ": ", paths[[nm]])
  }

  manifest <- tibble(
    input = names(paths),
    basename = unname(DE_FILE_BASENAMES[names(paths)]),
    resolved_path = unname(paths)
  )

  write_csv(
    manifest,
    file.path(QC_DIR, "unloading_recovery_DE_file_manifest.csv")
  )

  osd21_unloading <- load_one_loading_effect(
    paths[["OSD21_unloading"]],
    study_id = "OSD21",
    role = "unloading",
    source_species = "mouse"
  )

  osd21_recovery <- load_one_loading_effect(
    paths[["OSD21_recovery"]],
    study_id = "OSD21",
    role = "recovery",
    source_species = "mouse"
  )

  osd540_unloading <- load_one_loading_effect(
    paths[["OSD540_unloading"]],
    study_id = "OSD540",
    role = "unloading",
    source_species = "human"
  )

  osd540_recovery <- load_one_loading_effect(
    paths[["OSD540_recovery"]],
    study_id = "OSD540",
    role = "recovery",
    source_species = "human"
  )

  gse271607_unloading <- combine_gse271607_tissues(
    paths[["GSE271607_VL"]],
    paths[["GSE271607_Soleus"]]
  )

  out <- bind_rows(
    osd21_unloading,
    osd21_recovery,
    osd540_unloading,
    osd540_recovery,
    gse271607_unloading
  )

  if (nrow(out) == 0) {
    stop("The required DE files were located, but no usable gene-level rows were read.")
  }

  got <- out %>%
    distinct(study_id, role) %>%
    arrange(role, study_id)

  message("\nMechanical-loading effects loaded:")
  print(got)

  expected_pairs <- tribble(
    ~study_id, ~role,
    "OSD21", "unloading",
    "OSD21", "recovery",
    "OSD540", "unloading",
    "OSD540", "recovery",
    "GSE271607", "unloading"
  )

  missing_pairs <- anti_join(
    expected_pairs,
    got,
    by = c("study_id", "role")
  )

  if (nrow(missing_pairs) > 0) {
    stop(
      "Some required unloading/recovery comparisons were not loaded:\n",
      paste(
        paste(missing_pairs$study_id, missing_pairs$role, sep = " / "),
        collapse = "\n"
      )
    )
  }

  message("\nGenes available per loading comparison:")
  print(
    out %>%
      count(study_id, role, name = "n_genes") %>%
      arrange(role, study_id)
  )

  out
}

# =============================================================================
# 5. DERIVE THE REVERSIBLE ALS–UNLOADING MODULE
# =============================================================================

build_sod1_consensus <- function(sod1, stages_use) {
  x <- sod1 %>%
    filter(stage %in% stages_use) %>%
    mutate(direction = signed_direction(log2FC)) %>%
    group_by(symbol) %>%
    summarise(
      n_stage_available = n_distinct(stage),
      n_stage_up = sum(direction == 1L),
      n_stage_down = sum(direction == -1L),
      n_stage_nonzero = sum(direction != 0L),
      median_sod1_log2FC = median(log2FC, na.rm = TRUE),
      min_sod1_p = if (all(is.na(p_value))) NA_real_ else min(p_value, na.rm = TRUE),
      min_sod1_FDR = if (all(is.na(FDR))) NA_real_ else min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sod1_direction = case_when(
        n_stage_up > n_stage_down ~ 1L,
        n_stage_down > n_stage_up ~ -1L,
        TRUE ~ 0L
      ),
      n_sod1_concordant = pmax(n_stage_up, n_stage_down)
    )

  x
}

build_loading_consensus <- function(load_effects) {
  unload <- load_effects %>%
    filter(role == "unloading") %>%
    mutate(direction = signed_direction(log2FC)) %>%
    dplyr::select(symbol, study_id, log2FC, direction, p_value, FDR) %>%
    pivot_wider(
      names_from = study_id,
      values_from = c(log2FC, direction, p_value, FDR),
      names_sep = "__"
    )

  recovery <- load_effects %>%
    filter(role == "recovery") %>%
    mutate(direction = signed_direction(log2FC)) %>%
    dplyr::select(symbol, study_id, log2FC, direction, p_value, FDR) %>%
    pivot_wider(
      names_from = study_id,
      values_from = c(log2FC, direction, p_value, FDR),
      names_sep = "__"
    )

  full_join(unload, recovery, by = "symbol", suffix = c("_unload", "_recovery"))
}

derive_module <- function(
  sod1,
  load_effects,
  stages_use = c("90", "120", "150"),
  min_sod1_support = MIN_SOD1_SUPPORT_FULL,
  min_unloading_support = MIN_UNLOADING_SUPPORT,
  min_recovery_reverse = MIN_RECOVERY_REVERSE
) {
  s <- build_sod1_consensus(sod1, stages_use)

  unload_long <- load_effects %>%
    filter(role == "unloading") %>%
    mutate(direction = signed_direction(log2FC)) %>%
    dplyr::select(symbol, study_id, unload_lfc = log2FC, unload_direction = direction,
           unload_p = p_value, unload_FDR = FDR)

  rec_long <- load_effects %>%
    filter(role == "recovery") %>%
    mutate(direction = signed_direction(log2FC)) %>%
    dplyr::select(symbol, study_id, recovery_lfc = log2FC, recovery_direction = direction,
           recovery_p = p_value, recovery_FDR = FDR)

  unload_sum <- unload_long %>%
    inner_join(s %>% dplyr::select(symbol, sod1_direction), by = "symbol") %>%
    group_by(symbol, sod1_direction) %>%
    summarise(
      n_unloading_available = n_distinct(study_id),
      n_unloading_concordant = sum(unload_direction == sod1_direction & unload_direction != 0L),
      n_unloading_discordant = sum(unload_direction == -sod1_direction & unload_direction != 0L),
      min_unloading_p = if (all(is.na(unload_p))) NA_real_ else min(unload_p, na.rm = TRUE),
      min_unloading_FDR = if (all(is.na(unload_FDR))) NA_real_ else min(unload_FDR, na.rm = TRUE),
      median_unloading_log2FC = median(unload_lfc, na.rm = TRUE),
      .groups = "drop"
    )

  rec_sum <- rec_long %>%
    inner_join(s %>% dplyr::select(symbol, sod1_direction), by = "symbol") %>%
    group_by(symbol, sod1_direction) %>%
    summarise(
      n_recovery_available = n_distinct(study_id),
      n_recovery_reverse = sum(recovery_direction == -sod1_direction & recovery_direction != 0L),
      n_recovery_reinforce = sum(recovery_direction == sod1_direction & recovery_direction != 0L),
      min_recovery_p = if (all(is.na(recovery_p))) NA_real_ else min(recovery_p, na.rm = TRUE),
      min_recovery_FDR = if (all(is.na(recovery_FDR))) NA_real_ else min(recovery_FDR, na.rm = TRUE),
      median_recovery_log2FC = median(recovery_lfc, na.rm = TRUE),
      .groups = "drop"
    )

  ann <- s %>%
    left_join(unload_sum, by = c("symbol", "sod1_direction")) %>%
    left_join(rec_sum, by = c("symbol", "sod1_direction")) %>%
    mutate(
      across(
        c(
          n_unloading_available, n_unloading_concordant, n_unloading_discordant,
          n_recovery_available, n_recovery_reverse, n_recovery_reinforce
        ),
        ~replace_na(.x, 0L)
      ),
      als_direction = sod1_direction,
      expected_exercise_direction = -als_direction,

      directional_module =
        sod1_direction != 0L &
        n_sod1_concordant >= min_sod1_support &
        n_unloading_concordant >= min_unloading_support &
        n_recovery_reverse >= min_recovery_reverse &
        n_recovery_reinforce == 0L,

      # The SOD1 input itself is already restricted to significant genes.
      # Therefore this stricter subset asks for additional statistical support
      # from the unloading and recovery datasets; it does not require SOD1
      # p/FDR columns to be present in the supplied significant-only tables.
      statistically_supported =
        directional_module &
        (
          (!is.na(min_unloading_FDR) & min_unloading_FDR <= SUPPORT_FDR) |
          (!is.na(min_unloading_p) & min_unloading_p <= SUPPORT_NOMINAL_P)
        ) &
        (
          (!is.na(min_recovery_FDR) & min_recovery_FDR <= SUPPORT_FDR) |
          (!is.na(min_recovery_p) & min_recovery_p <= SUPPORT_NOMINAL_P)
        ),

      # Keep the recovery-stringency criterion separate from statistical
      # support. This avoids calling a gene "high confidence" solely because it
      # reverses in both recovery datasets.
      two_recovery_models =
        directional_module &
        n_recovery_available >= 2L &
        n_recovery_reverse >= 2L,

      # High-confidence = directional module + additional statistical support
      # + the configured recovery-stringency criterion.
      high_confidence =
        statistically_supported &
        if (isTRUE(REQUIRE_BOTH_RECOVERY_FOR_HIGH_CONF)) {
          n_recovery_available >= 2L & n_recovery_reverse >= 2L
        } else {
          n_recovery_reverse >= 1L
        }
    )

  # Add per-study effects for transparent inspection.
  per_study <- load_effects %>%
    dplyr::select(symbol, study_id, role, log2FC, p_value, FDR) %>%
    unite("key", role, study_id, sep = "__") %>%
    pivot_wider(
      names_from = key,
      values_from = c(log2FC, p_value, FDR),
      names_sep = "__"
    )

  ann %>%
    left_join(per_study, by = "symbol") %>%
    arrange(
      desc(high_confidence),
      desc(two_recovery_models),
      desc(statistically_supported),
      desc(directional_module),
      desc(n_unloading_concordant),
      desc(n_recovery_reverse),
      desc(n_sod1_concordant),
      symbol
    )
}

# =============================================================================
# 6. FORMAL PERMUTATION TEST OF CROSS-DATASET CONVERGENCE
# =============================================================================

module_permutation_test <- function(
  sod1,
  load_effects,
  stages_use = c("90", "120", "150"),
  min_sod1_support = MIN_SOD1_SUPPORT_FULL,
  min_unloading_support = MIN_UNLOADING_SUPPORT,
  min_recovery_reverse = MIN_RECOVERY_REVERSE,
  n_perm = N_PERM
) {
  s <- build_sod1_consensus(sod1, stages_use) %>%
    filter(
      sod1_direction != 0L,
      n_sod1_concordant >= min_sod1_support
    ) %>%
    dplyr::select(symbol, sod1_direction)

  studies <- load_effects %>%
    mutate(direction = signed_direction(log2FC)) %>%
    dplyr::select(symbol, study_id, role, direction) %>%
    distinct() %>%
    mutate(col = paste(role, study_id, sep = "__")) %>%
    dplyr::select(symbol, col, direction) %>%
    pivot_wider(names_from = col, values_from = direction, values_fill = 0L)

  dat <- inner_join(s, studies, by = "symbol")

  unload_cols <- grep("^unloading__", colnames(dat), value = TRUE)
  recovery_cols <- grep("^recovery__", colnames(dat), value = TRUE)

  if (length(unload_cols) < min_unloading_support) {
    stop("Not enough unloading studies were found for the permutation test.")
  }
  if (length(recovery_cols) < min_recovery_reverse) {
    stop("Not enough recovery studies were found for the permutation test.")
  }

  count_module <- function(d) {
    als <- d$sod1_direction

    u_mat <- as.matrix(d[, unload_cols, drop = FALSE])
    r_mat <- as.matrix(d[, recovery_cols, drop = FALSE])

    u_conc <- rowSums(u_mat == als & u_mat != 0L)
    r_rev <- rowSums(r_mat == -als & r_mat != 0L)
    r_reinforce <- rowSums(r_mat == als & r_mat != 0L)

    sum(
      u_conc >= min_unloading_support &
      r_rev >= min_recovery_reverse &
      r_reinforce == 0L
    )
  }

  observed <- count_module(dat)

  null_counts <- numeric(n_perm)
  n <- nrow(dat)

  for (b in seq_len(n_perm)) {
    dperm <- dat

    # Keep the recurrently significant SOD1 genes and their ALS directions fixed.
    # Independently permute gene labels for each unloading and recovery study.
    # This preserves each loading dataset's marginal direction distribution while
    # destroying gene-specific ALS-loading-recovery alignment.
    for (cc in c(unload_cols, recovery_cols)) {
      dperm[[cc]] <- sample(dperm[[cc]], size = n, replace = FALSE)
    }

    null_counts[b] <- count_module(dperm)
  }

  empirical_p <- (1 + sum(null_counts >= observed)) / (n_perm + 1)

  summary <- tibble(
    stages_used = paste(stages_use, collapse = ","),
    n_recurrently_significant_sod1_genes_tested = nrow(dat),
    observed_module_size = observed,
    null_mean = mean(null_counts),
    null_sd = sd(null_counts),
    null_95th_percentile = unname(quantile(null_counts, 0.95)),
    empirical_p_one_sided = empirical_p,
    n_permutations = n_perm
  )

  list(summary = summary, null = null_counts)
}

# =============================================================================
# 7. SOD1 LEAVE-ONE-STAGE-OUT TRAJECTORY
# =============================================================================

test_module_in_sod1_stage <- function(stage_effects, module, heldout_stage) {
  # IMPORTANT FOR SIGNIFICANT-ONLY SOD1 INPUT:
  # The held-out stage table contains only significant genes.
  # Therefore:
  #   - a matching row means the module gene is significantly dysregulated there;
  #   - no matching row means "not significantly detected", NOT log2FC = 0.
  # We test recurrence/direction, not a whole-module effect-size distribution.

  m <- module %>%
    filter(directional_module) %>%
    dplyr::select(symbol, als_direction)

  heldout_sig <- stage_effects %>%
    filter(stage == heldout_stage) %>%
    # The SOD1 workbook provides significant-gene effect sizes plus p/FDR, but
    # no test-statistic column. LOSO recurrence does not need a statistic, so
    # select only fields that are actually present.
    dplyr::select(symbol, log2FC, p_value, FDR)

  x <- inner_join(
    heldout_sig,
    m,
    by = "symbol"
  ) %>%
    mutate(
      signed_log2FC = log2FC * als_direction,
      aligned = signed_log2FC > 0,
      recurrence_class = if_else(
        aligned,
        "Significant in expected ALS direction",
        "Significant in opposite direction"
      )
    )

  n_module <- nrow(m)
  n_recurrent <- nrow(x)
  n_aligned <- sum(x$aligned)
  n_opposite <- n_recurrent - n_aligned

  summary <- tibble(
    stage = heldout_stage,
    n_independently_defined_module_genes = n_module,
    n_significant_at_heldout_stage = n_recurrent,
    n_significant_expected_direction = n_aligned,
    n_significant_opposite_direction = n_opposite,

    # Fraction of the independently defined module that is significant at the
    # held-out stage, regardless of direction.
    significant_recurrence_fraction =
      if_else(n_module > 0, n_recurrent / n_module, NA_real_),

    # Primary trajectory metric: fraction of the independently defined module
    # that is BOTH significant and in the predicted ALS direction at that stage.
    expected_direction_recurrence_fraction =
      if_else(n_module > 0, n_aligned / n_module, NA_real_),

    # Conditional directional test among recurrent significant genes only.
    direction_concordance_among_recurrent =
      if_else(n_recurrent > 0, n_aligned / n_recurrent, NA_real_),

    direction_sign_test_p =
      if (n_recurrent > 0)
        safe_binom_p(n_aligned, n_recurrent, alternative = "greater")
      else NA_real_,

    # Descriptive only: this is calculated among recurrent significant genes,
    # NOT across the whole module, because missing genes have unknown log2FC.
    median_signed_log2FC_among_recurrent =
      if (n_recurrent > 0) median(x$signed_log2FC, na.rm = TRUE) else NA_real_,
    mean_signed_log2FC_among_recurrent =
      if (n_recurrent > 0) mean(x$signed_log2FC, na.rm = TRUE) else NA_real_
  )

  list(summary = summary, genes = x)
}

run_sod1_loso <- function(sod1, load_effects) {
  stages <- c("90", "120", "150")
  summaries <- list()
  gene_tables <- list()
  modules <- list()

  for (heldout in stages) {
    train <- setdiff(stages, heldout)

    # With only two training stages, require agreement in BOTH.
    mod <- derive_module(
      sod1,
      load_effects,
      stages_use = train,
      min_sod1_support = length(train),
      min_unloading_support = MIN_UNLOADING_SUPPORT,
      min_recovery_reverse = MIN_RECOVERY_REVERSE
    )

    modules[[heldout]] <- mod %>%
      filter(directional_module) %>%
      mutate(heldout_stage = heldout, training_stages = paste(train, collapse = ","))

    tst <- test_module_in_sod1_stage(sod1, mod, heldout)

    summaries[[heldout]] <- tst$summary %>%
      mutate(
        training_stages = paste(train, collapse = ","),
        n_training_module_genes = sum(mod$directional_module)
      )

    gene_tables[[heldout]] <- tst$genes %>%
      mutate(
        heldout_stage = heldout,
        training_stages = paste(train, collapse = ",")
      )
  }

  list(
    summary = bind_rows(summaries),
    genes = bind_rows(gene_tables),
    modules = bind_rows(modules)
  )
}

# =============================================================================
# 8. GSE250455: PAIRED HUMAN ALS PRE/POST RESISTANCE EXERCISE
# =============================================================================

read_gse250455 <- function(path) {
  if (!file.exists(path)) stop("GSE250455 file not found: ", path)

  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )

  if (!all(c("gene_id", "gene_name") %in% colnames(x))) {
    stop("GSE250455 file must contain gene_id and gene_name columns.")
  }

  sample_cols <- setdiff(colnames(x), c("gene_id", "gene_name"))

  meta <- tibble(sample = sample_cols) %>%
    mutate(
      participant = str_remove(sample, "-(Pre|Post)$"),
      condition = case_when(
        str_detect(sample, "-Pre$") ~ "Pre",
        str_detect(sample, "-Post$") ~ "Post",
        TRUE ~ NA_character_
      ),
      condition = factor(condition, levels = c("Pre", "Post"))
    )

  if (any(is.na(meta$condition))) {
    stop("Could not parse Pre/Post status from all GSE250455 sample names.")
  }

  pairing <- meta %>%
    count(participant, condition) %>%
    pivot_wider(names_from = condition, values_from = n, values_fill = 0)

  if (!all(pairing$Pre == 1 & pairing$Post == 1)) {
    stop("GSE250455 samples do not form complete one-to-one Pre/Post pairs.")
  }

  # Collapse duplicate gene symbols by SUMMING RSEM expected counts.
  # This avoids arbitrary transcript/gene-symbol duplicate selection.
  counts_df <- x %>%
    mutate(gene_name = clean_symbol(gene_name)) %>%
    filter(!is.na(gene_name)) %>%
    group_by(gene_name) %>%
    summarise(
      across(all_of(sample_cols), ~sum(safe_num(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  counts <- as.matrix(counts_df[, sample_cols, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- counts_df$gene_name

  list(counts = counts, meta = meta)
}

run_gse250455_edger <- function(counts, meta) {
  # Reorder metadata to count matrix columns.
  # Reorder metadata to match the count-matrix columns.
  # Use an explicit dplyr namespace because IRanges/S4Vectors can mask slice().
  meta <- meta %>%
    dplyr::slice(match(colnames(counts), sample))

  participant <- factor(meta$participant)
  condition <- factor(meta$condition, levels = c("Pre", "Post"))

  design <- model.matrix(~ participant + condition)
  colnames(design) <- make.names(colnames(design))

  y <- edgeR::DGEList(counts = counts)

  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")

  y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)

  coef_name <- grep("^conditionPost$", colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    stop(
      "Could not identify conditionPost coefficient. Design columns: ",
      paste(colnames(design), collapse = ", ")
    )
  }

  qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table %>%
    rownames_to_column("symbol") %>%
    as_tibble() %>%
    transmute(
      symbol = clean_symbol(symbol),
      log2FC = logFC,
      logCPM = logCPM,
      F_stat = F,
      p_value = PValue,
      FDR = FDR
    ) %>%
    arrange(FDR, p_value)

  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)

  list(
    y = y,
    design = design,
    fit = fit,
    qlf = qlf,
    results = tt,
    logCPM = logcpm,
    meta = meta,
    keep = keep
  )
}

save_gse250455_qc <- function(gse_fit) {
  logcpm <- gse_fit$logCPM
  meta <- gse_fit$meta

  pca <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)

  var_exp <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

  pdat <- as_tibble(pca$x[, 1:2, drop = FALSE], rownames = "sample") %>%
    left_join(meta, by = "sample")

  p <- ggplot(pdat, aes(PC1, PC2, label = participant, shape = condition)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = Inf) +
    labs(
      title = "GSE250455 paired ALS muscle RNA-seq",
      subtitle = "TMM-normalized logCPM",
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2])
    ) +
    theme_bw(base_size = 11)

  ggsave(
    file.path(FIG_DIR, "GSE250455_PCA.svg"),
    p, width = 7.5, height = 6, device = "svg"
  )

  lib <- tibble(
    sample = colnames(gse_fit$y$counts),
    library_size = gse_fit$y$samples$lib.size,
    norm_factor = gse_fit$y$samples$norm.factors
  ) %>%
    left_join(meta, by = "sample")

  write_csv(lib, file.path(QC_DIR, "GSE250455_library_sizes_and_norm_factors.csv"))
}

module_expression_score <- function(logcpm, meta, module) {
  m <- module %>%
    dplyr::filter(directional_module) %>%
    dplyr::transmute(
      symbol,
      als_direction,
      expected_exercise_direction = -als_direction
    ) %>%
    dplyr::filter(symbol %in% rownames(logcpm))

  if (nrow(m) < 2) {
    stop("Too few module genes are present in GSE250455 to calculate a module score.")
  }

  expr <- logcpm[m$symbol, , drop = FALSE]

  # Gene-wise z score across samples, then orient every gene so a POSITIVE
  # change means movement in the predicted recovery/exercise direction.
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0

  orient <- m$expected_exercise_direction
  names(orient) <- m$symbol
  z_oriented <- z * orient[rownames(z)]

  score <- colMeans(z_oriented, na.rm = TRUE)

  score_df <- tibble(
    sample = names(score),
    recovery_module_score = unname(score)
  ) %>%
    dplyr::left_join(meta, by = "sample")

  wide <- score_df %>%
    dplyr::select(participant, condition, recovery_module_score) %>%
    tidyr::pivot_wider(names_from = condition, values_from = recovery_module_score) %>%
    dplyr::mutate(delta_Post_minus_Pre = Post - Pre)

  paired_t <- function(alternative) {
    if (nrow(wide) < 3) return(NA_real_)
    stats::t.test(
      wide$Post, wide$Pre,
      paired = TRUE,
      alternative = alternative
    )$p.value
  }

  paired_wilcox <- function(alternative) {
    if (nrow(wide) < 3) return(NA_real_)
    suppressWarnings(
      stats::wilcox.test(
        wide$Post, wide$Pre,
        paired = TRUE,
        alternative = alternative,
        exact = FALSE
      )$p.value
    )
  }

  summary <- tibble(
    n_participants = nrow(wide),
    n_module_genes_measured = nrow(m),
    mean_delta = mean(wide$delta_Post_minus_Pre),
    median_delta = median(wide$delta_Post_minus_Pre),

    # Positive Post-Pre = predicted recovery direction.
    paired_t_p_predicted_recovery = paired_t("greater"),
    paired_t_p_opposite_direction = paired_t("less"),
    paired_t_p_two_sided = paired_t("two.sided"),

    paired_wilcox_p_predicted_recovery = paired_wilcox("greater"),
    paired_wilcox_p_opposite_direction = paired_wilcox("less"),
    paired_wilcox_p_two_sided = paired_wilcox("two.sided")
  )

  list(score_df = score_df, paired = wide, summary = summary, genes = m)
}

test_exercise_module_gene_level <- function(gse_results, module) {
  m <- module %>%
    dplyr::filter(directional_module) %>%
    dplyr::select(symbol, als_direction) %>%
    dplyr::mutate(expected_exercise_direction = -als_direction)

  x <- dplyr::inner_join(gse_results, m, by = "symbol") %>%
    dplyr::mutate(
      signed_exercise_log2FC = log2FC * expected_exercise_direction,
      direction_class = dplyr::case_when(
        signed_exercise_log2FC > 0 ~ "Predicted recovery direction",
        signed_exercise_log2FC < 0 ~ "Opposite / ALS-associated direction",
        TRUE ~ "No directional change"
      ),
      aligned_with_predicted_recovery = signed_exercise_log2FC > 0,
      aligned_with_opposite_direction = signed_exercise_log2FC < 0,
      signed_exercise_stat = sign(log2FC) * sqrt(pmax(F_stat, 0)) *
        expected_exercise_direction
    )

  n_nonzero <- sum(x$signed_exercise_log2FC != 0, na.rm = TRUE)
  n_predicted <- sum(x$signed_exercise_log2FC > 0, na.rm = TRUE)
  n_opposite <- sum(x$signed_exercise_log2FC < 0, na.rm = TRUE)

  summary <- tibble(
    n_module_genes = nrow(m),
    n_module_genes_measured = nrow(x),
    n_nonzero_direction = n_nonzero,
    n_predicted_recovery_direction = n_predicted,
    n_opposite_ALS_direction = n_opposite,
    proportion_predicted_recovery =
      if (n_nonzero > 0) n_predicted / n_nonzero else NA_real_,
    proportion_opposite_ALS_direction =
      if (n_nonzero > 0) n_opposite / n_nonzero else NA_real_,

    sign_test_p_predicted_recovery =
      safe_binom_p(n_predicted, n_nonzero, alternative = "greater"),
    sign_test_p_opposite_direction =
      safe_binom_p(n_opposite, n_nonzero, alternative = "greater"),
    sign_test_p_two_sided =
      safe_binom_p(n_predicted, n_nonzero, alternative = "two.sided"),

    median_signed_log2FC = median(x$signed_exercise_log2FC, na.rm = TRUE),
    mean_signed_log2FC = mean(x$signed_exercise_log2FC, na.rm = TRUE),

    wilcox_p_predicted_recovery =
      safe_wilcox_p(x$signed_exercise_log2FC, alternative = "greater"),
    wilcox_p_opposite_direction =
      safe_wilcox_p(x$signed_exercise_log2FC, alternative = "less"),
    wilcox_p_two_sided =
      safe_wilcox_p(x$signed_exercise_log2FC, alternative = "two.sided"),

    t_test_p_predicted_recovery =
      safe_t_p(x$signed_exercise_log2FC, alternative = "greater"),
    t_test_p_opposite_direction =
      safe_t_p(x$signed_exercise_log2FC, alternative = "less"),
    t_test_p_two_sided =
      safe_t_p(x$signed_exercise_log2FC, alternative = "two.sided")
  )

  list(summary = summary, genes = x)
}

exercise_expression_matched_permutation <- function(
  gse_results,
  module,
  n_perm = N_PERM
) {
  bg <- gse_results %>%
    dplyr::filter(is.finite(logCPM), is.finite(F_stat), !is.na(symbol)) %>%
    dplyr::mutate(
      pseudo_t = sign(log2FC) * sqrt(pmax(F_stat, 0)),
      expr_bin = dplyr::ntile(logCPM, 10)
    )

  m <- module %>%
    dplyr::filter(directional_module) %>%
    dplyr::select(symbol, als_direction) %>%
    dplyr::mutate(expected_exercise_direction = -als_direction) %>%
    dplyr::inner_join(
      bg %>% dplyr::select(symbol, pseudo_t, expr_bin),
      by = "symbol"
    )

  if (nrow(m) < 3) {
    return(tibble(
      n_module_genes = nrow(m),
      observed_mean_signed_stat = NA_real_,
      null_mean = NA_real_,
      null_sd = NA_real_,
      empirical_p_predicted_recovery = NA_real_,
      empirical_p_opposite_direction = NA_real_,
      empirical_p_two_sided = NA_real_,
      n_permutations = n_perm
    ))
  }

  observed <- mean(m$pseudo_t * m$expected_exercise_direction)

  null <- numeric(n_perm)

  # Match random genes to the module's expression-decile structure.
  module_by_bin <- split(m, m$expr_bin)
  bg_by_bin <- split(bg, bg$expr_bin)

  for (b in seq_len(n_perm)) {
    vals <- numeric(0)

    for (bn in names(module_by_bin)) {
      mm <- module_by_bin[[bn]]
      bb <- bg_by_bin[[bn]]

      if (is.null(bb) || nrow(bb) == 0) next

      n_take <- nrow(mm)
      replace_flag <- nrow(bb) < n_take

      sampled <- bb %>%
        dplyr::slice_sample(n = n_take, replace = replace_flag)

      # Preserve the module's expected-direction composition within each
      # expression bin.
      dirs <- sample(mm$expected_exercise_direction, size = n_take, replace = FALSE)
      vals <- c(vals, sampled$pseudo_t * dirs)
    }

    null[b] <- mean(vals)
  }

  null_mean <- mean(null)
  obs_dev <- abs(observed - null_mean)
  null_dev <- abs(null - null_mean)

  tibble(
    n_module_genes = nrow(m),
    observed_mean_signed_stat = observed,
    null_mean = null_mean,
    null_sd = sd(null),
    null_2.5th_percentile = unname(quantile(null, 0.025)),
    null_5th_percentile = unname(quantile(null, 0.05)),
    null_95th_percentile = unname(quantile(null, 0.95)),
    null_97.5th_percentile = unname(quantile(null, 0.975)),

    # Positive oriented statistic = predicted recovery direction.
    empirical_p_predicted_recovery =
      (1 + sum(null >= observed)) / (n_perm + 1),

    # Negative oriented statistic = opposite / ALS-associated direction.
    empirical_p_opposite_direction =
      (1 + sum(null <= observed)) / (n_perm + 1),

    # Two-sided test around the empirical null centre.
    empirical_p_two_sided =
      (1 + sum(null_dev >= obs_dev)) / (n_perm + 1),

    n_permutations = n_perm
  )
}

exercise_fgsea <- function(gse_results, module) {
  ranks <- gse_results %>%
    filter(!is.na(symbol), is.finite(F_stat)) %>%
    mutate(rank_stat = sign(log2FC) * sqrt(pmax(F_stat, 0))) %>%
    arrange(desc(abs(rank_stat))) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(symbol, rank_stat) %>%
    deframe()

  m <- module %>%
    filter(directional_module)

  pathways <- list(
    ALS_module_down_expected_to_increase_with_exercise =
      m %>% filter(als_direction == -1L) %>% pull(symbol) %>% intersect(names(ranks)),
    ALS_module_up_expected_to_decrease_with_exercise =
      m %>% filter(als_direction == 1L) %>% pull(symbol) %>% intersect(names(ranks))
  )

  pathways <- pathways[lengths(pathways) >= 3]

  if (length(pathways) == 0) {
    return(tibble())
  }

  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 3,
    maxSize = 500,
    eps = 0
  ) %>%
    as_tibble() %>%
    arrange(padj, pval)
}

save_exercise_module_plot <- function(score_obj) {
  d <- score_obj$score_df

  p <- ggplot(
    d,
    aes(
      x = condition,
      y = recovery_module_score,
      group = participant
    )
  ) +
    geom_line(alpha = 0.65) +
    geom_point(size = 2.6) +
    labs(
      title = "Reversible ALS–unloading module in human ALS muscle",
      subtitle = "Positive Post–Pre change = movement in the predicted recovery direction",
      x = NULL,
      y = "Direction-oriented module expression score"
    ) +
    theme_bw(base_size = 11)

  ggsave(
    file.path(FIG_DIR, "GSE250455_recovery_module_paired_score.svg"),
    p, width = 7, height = 5.5, device = "svg"
  )
}

save_exercise_gene_plot <- function(gene_results) {
  if (nrow(gene_results) == 0) return(invisible(NULL))

  top <- gene_results %>%
    dplyr::arrange(dplyr::desc(abs(signed_exercise_log2FC))) %>%
    # slice_head() already returns all available rows when fewer than n exist;
    # n must be a constant, so do not use n() inside its n argument.
    dplyr::slice_head(n = 40) %>%
    dplyr::mutate(
      symbol = factor(symbol, levels = rev(symbol))
    )

  p <- ggplot(top, aes(signed_exercise_log2FC, symbol)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(size = 2.3) +
    labs(
      title = "Human ALS exercise response of reversible module genes",
      subtitle = "Positive values indicate movement in the predicted recovery direction",
      x = "Direction-oriented exercise log2FC",
      y = NULL
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(FIG_DIR, "GSE250455_module_gene_effects.svg"),
    p, width = 8, height = max(6, 0.22 * nrow(top) + 2), device = "svg"
  )
}

# =============================================================================
# 9. SENSITIVITY ANALYSIS ACROSS MODULE LFC THRESHOLDS
# =============================================================================

module_threshold_sensitivity <- function(sod1, load_effects, thresholds = c(0, 0.05, 0.10, 0.20)) {
  old <- MIN_ABS_LFC
  on.exit(assign("MIN_ABS_LFC", old, envir = .GlobalEnv), add = TRUE)

  out <- list()

  for (thr in thresholds) {
    assign("MIN_ABS_LFC", thr, envir = .GlobalEnv)

    mod <- derive_module(
      sod1,
      load_effects,
      stages_use = c("90", "120", "150"),
      min_sod1_support = MIN_SOD1_SUPPORT_FULL,
      min_unloading_support = MIN_UNLOADING_SUPPORT,
      min_recovery_reverse = MIN_RECOVERY_REVERSE
    )

    out[[as.character(thr)]] <- tibble(
      min_abs_lfc = thr,
      n_directional_module = sum(mod$directional_module),
      n_statistically_supported = sum(mod$statistically_supported),
      n_two_recovery_models = sum(mod$two_recovery_models),
      n_high_confidence = sum(mod$high_confidence)
    )
  }

  bind_rows(out)
}


# =============================================================================
# 10. TARGETED GENES OF INTEREST
# =============================================================================

TARGET_GENES <- c(
  "PGM1", "TPI1", "PGK1", "ENO3", "LDHA", "FBP2", "PRKAA2",
  "PDK4", "PVALB"
)

# =============================================================================
# 11. MAIN ANALYSIS
# =============================================================================

message("============================================================")
message("REVERSIBLE ALS–UNLOADING MODULE ANALYSIS")
message("============================================================")
message("SOD1 mode: mutation-specific significant genes = intersection of Tables S3 + S4")
message("Table S3: hSOD1G93A vs hSOD1WT; Table S4: hSOD1G93A vs NTg.")
message("Absent SOD1 genes are treated as not significantly detected, never as zero effect.")

# ---- SOD1 ----
# First read the mutation-specific significant mouse genes exactly as reported
# in Tables S3/S4, then convert them to human ortholog space for all
# cross-species analyses.
sod1_mouse <- load_sod1_stage_effects()

write_csv(
  sod1_mouse,
  file.path(TABLE_DIR, "SOD1_stage_effects_MOUSE_original_symbols.csv")
)

sod1 <- harmonize_sod1_to_human(sod1_mouse)

write_csv(
  sod1,
  file.path(TABLE_DIR, "SOD1_stage_effects_HUMAN_ORTHOLOG_space.csv")
)

writeLines(
  c(
    "SOD1 INPUT MODE: significant-only mutation-specific gene sets read directly from SOD_data.xlsx.",
    "Table S3 contains hSOD1G93A vs hSOD1WT significant mRNAs.",
    "Table S4 contains hSOD1G93A vs NTg significant mRNAs.",
    "At each age, only genes significant in BOTH comparisons and with the same direction are retained.",
    "Expected mouse intersections before ortholog mapping: 558 genes at 90 d, 544 at 120 d, 685 at 150 d.",
    "For cross-species analyses, mouse SOD1 genes are mapped to human ortholog symbols with babelgene.",
    "babelgene settings: species='mouse', human=FALSE, min_support=3, top=TRUE.",
    "OSD-21 mouse genes are mapped with the same method.",
    "OSD-540 and GSE271607 human Ensembl IDs are mapped to official human SYMBOLs using org.Hs.eg.db.",
    "GSE250455 already supplies human gene symbols.",
    "Thus all cross-dataset module analyses are performed in human ortholog gene-symbol space.",
    "A gene absent from a stage-specific SOD1 intersection is interpreted as not significantly dysregulated at that stage."
  ),
  file.path(QC_DIR, "cross_species_harmonisation_note.txt")
)

message("\nSOD1 rows by stage BEFORE and AFTER mouse-to-human mapping:")
sod1_stage_counts <- sod1_mouse %>%
  count(stage, name = "n_mouse_significant_genes") %>%
  left_join(
    sod1 %>% count(stage, name = "n_human_ortholog_genes"),
    by = "stage"
  ) %>%
  mutate(
    percent_retained = 100 * n_human_ortholog_genes / n_mouse_significant_genes
  )

print(sod1_stage_counts)

write_csv(
  sod1_stage_counts,
  file.path(QC_DIR, "SOD1_input_coverage_and_ortholog_mapping.csv")
)

# ---- Previous unloading/recovery DE outputs ----
# Load the exact six previously generated gene-level DE files by basename.
# OSD-21 and the human datasets may reside in different result directories.
load_effects <- load_unloading_recovery_effects()

write_csv(
  load_effects,
  file.path(TABLE_DIR, "unloading_recovery_gene_effects_HUMAN_ORTHOLOG_space.csv")
)

# QC: all three unloading studies should now share genes with the SOD1 human-
# ortholog universe. This specifically guards against the previous problem in
# which OSD-540/GSE271607 Ensembl IDs were mislabeled as "symbol".
sod1_human_universe <- unique(sod1$symbol)

cross_species_overlap_qc <- load_effects %>%
  filter(role == "unloading") %>%
  group_by(study_id) %>%
  summarise(
    n_genes = n(),
    n_overlap_with_SOD1_human_ortholog_space =
      sum(symbol %in% sod1_human_universe),
    .groups = "drop"
  )

message("\nCross-species harmonisation QC:")
print(cross_species_overlap_qc)

if (any(cross_species_overlap_qc$n_overlap_with_SOD1_human_ortholog_space == 0)) {
  stop(
    "Cross-species harmonisation failed: at least one unloading study has zero ",
    "overlap with the SOD1 human-ortholog gene space. See QC table."
  )
}

write_csv(
  cross_species_overlap_qc,
  file.path(QC_DIR, "cross_species_overlap_QC.csv")
)

# ---- Full reversible module ----
full_module <- derive_module(
  sod1,
  load_effects,
  stages_use = c("90", "120", "150"),
  min_sod1_support = MIN_SOD1_SUPPORT_FULL,
  min_unloading_support = MIN_UNLOADING_SUPPORT,
  min_recovery_reverse = MIN_RECOVERY_REVERSE
)

write_csv(
  full_module,
  file.path(TABLE_DIR, "reversible_ALS_unloading_module_all_genes.csv")
)

directional_module <- full_module %>% dplyr::filter(directional_module)
supported_module <- full_module %>% dplyr::filter(statistically_supported)
two_recovery_module <- full_module %>% dplyr::filter(two_recovery_models)
high_conf_module <- full_module %>% dplyr::filter(high_confidence)

write_csv(
  directional_module,
  file.path(TABLE_DIR, "reversible_ALS_unloading_module_DIRECTIONAL.csv")
)
write_csv(
  supported_module,
  file.path(TABLE_DIR, "reversible_ALS_unloading_module_SUPPORTED.csv")
)
write_csv(
  two_recovery_module,
  file.path(TABLE_DIR, "reversible_ALS_unloading_module_TWO_RECOVERY_MODELS.csv")
)
write_csv(
  high_conf_module,
  file.path(TABLE_DIR, "reversible_ALS_unloading_module_HIGH_CONFIDENCE.csv")
)

module_counts <- tibble(
  module = c(
    "Directional",
    "Statistically_supported",
    "Two_recovery_models",
    "High_confidence_supported_and_two_recovery"
  ),
  n_genes = c(
    nrow(directional_module),
    nrow(supported_module),
    nrow(two_recovery_module),
    nrow(high_conf_module)
  )
)

write_csv(module_counts, file.path(TABLE_DIR, "module_counts.csv"))

message("\nModule counts:")
print(module_counts)

# ---- Formal convergence permutation test ----
perm <- module_permutation_test(
  sod1,
  load_effects,
  stages_use = c("90", "120", "150"),
  min_sod1_support = MIN_SOD1_SUPPORT_FULL,
  min_unloading_support = MIN_UNLOADING_SUPPORT,
  min_recovery_reverse = MIN_RECOVERY_REVERSE,
  n_perm = N_PERM
)

write_csv(
  perm$summary,
  file.path(TABLE_DIR, "module_convergence_permutation_test.csv")
)

null_df <- tibble(module_size = perm$null)

p_null <- ggplot(null_df, aes(module_size)) +
  geom_histogram(bins = 40) +
  geom_vline(
    xintercept = perm$summary$observed_module_size,
    linetype = 2,
    linewidth = 1
  ) +
  labs(
    title = "Permutation test of reversible ALS–unloading convergence",
    subtitle = paste0(
      "Observed module size = ", perm$summary$observed_module_size,
      "; empirical P = ", signif(perm$summary$empirical_p_one_sided, 3)
    ),
    x = "Module size under gene-label permutation",
    y = "Permutations"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(FIG_DIR, "module_convergence_permutation_null.svg"),
  p_null, width = 7, height = 5.5, device = "svg"
)

# ---- Threshold robustness ----
sens <- module_threshold_sensitivity(sod1, load_effects)
write_csv(
  sens,
  file.path(TABLE_DIR, "module_LFC_threshold_sensitivity.csv")
)

# ---- SOD1 leave-one-stage-out trajectory ----
loso <- run_sod1_loso(sod1, load_effects)

write_csv(
  loso$summary,
  file.path(TABLE_DIR, "SOD1_LOSO_module_trajectory_summary.csv")
)
write_csv(
  loso$genes,
  file.path(TABLE_DIR, "SOD1_LOSO_module_trajectory_gene_results.csv")
)
write_csv(
  loso$modules,
  file.path(TABLE_DIR, "SOD1_LOSO_modules.csv")
)

trajectory_plot_data <- loso$summary %>%
  mutate(
    stage = factor(stage, levels = c("90", "120", "150")),
    stage_label = factor(
      stage,
      levels = c("90", "120", "150"),
      labels = c(
        "90 d\nPresymptomatic",
        "120 d\nEarly symptomatic",
        "150 d\nLate symptomatic"
      )
    )
  )

# Primary significant-only trajectory figure:
# genes absent from the held-out significant list count as "not recapitulated",
# never as zero-effect genes.
p_traj <- ggplot(
  trajectory_plot_data,
  aes(stage_label, expected_direction_recurrence_fraction, group = 1)
) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Independent recurrence of the reversible ALS–unloading module",
    subtitle = paste0(
      "Each stage is tested using a module derived without that stage; ",
      "SOD1 inputs contain significant genes only"
    ),
    x = NULL,
    y = "Module genes significant in the expected ALS direction"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(FIG_DIR, "SOD1_LOSO_expected_direction_recurrence.svg"),
  p_traj, width = 7.5, height = 5.5, device = "svg"
)

p_recurrence <- ggplot(
  trajectory_plot_data,
  aes(stage_label, significant_recurrence_fraction, group = 1)
) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Significant recurrence across SOD1 disease stages",
    subtitle = "Genes absent from a significant-only stage table are classified as not significantly detected",
    x = NULL,
    y = "Independently defined module genes significant at held-out stage"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(FIG_DIR, "SOD1_LOSO_significant_recurrence.svg"),
  p_recurrence, width = 7.5, height = 5.5, device = "svg"
)

# Conditional direction concordance among only those module genes that recur
# significantly at each held-out stage.
p_direction <- ggplot(
  trajectory_plot_data,
  aes(stage_label, direction_concordance_among_recurrent, group = 1)
) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Direction of recurrent significant module genes",
    subtitle = "Dashed line = 50% directional concordance",
    x = NULL,
    y = "Recurrent genes in expected ALS direction"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(FIG_DIR, "SOD1_LOSO_direction_concordance_among_recurrent.svg"),
  p_direction, width = 7.5, height = 5.5, device = "svg"
)

# =============================================================================
# 12. HUMAN ALS EXERCISE VALIDATION — GSE250455
# =============================================================================

message("\n============================================================")
message("GSE250455 HUMAN ALS RESISTANCE-EXERCISE VALIDATION")
message("============================================================")

gse <- read_gse250455(GSE250455_FILE)

message(
  "GSE250455: ", nrow(gse$counts), " unique gene symbols; ",
  ncol(gse$counts), " samples; ",
  n_distinct(gse$meta$participant), " paired participants."
)

gse_fit <- run_gse250455_edger(gse$counts, gse$meta)

write_csv(
  gse_fit$results,
  file.path(TABLE_DIR, "GSE250455_edgeR_paired_Post_vs_Pre.csv")
)

save_gse250455_qc(gse_fit)

# Primary validation uses the directional module.
exercise_gene <- test_exercise_module_gene_level(
  gse_fit$results,
  full_module
)

write_csv(
  exercise_gene$summary,
  file.path(TABLE_DIR, "GSE250455_module_gene_level_summary.csv")
)
write_csv(
  exercise_gene$genes,
  file.path(TABLE_DIR, "GSE250455_module_gene_level_results.csv")
)

# Secondary sensitivity analysis using the corrected high-confidence subset:
# statistically supported AND stringent recovery support.
high_conf_for_validation <- full_module %>%
  dplyr::mutate(directional_module = high_confidence)

exercise_gene_high_conf <- test_exercise_module_gene_level(
  gse_fit$results,
  high_conf_for_validation
)

write_csv(
  exercise_gene_high_conf$summary,
  file.path(TABLE_DIR, "GSE250455_HIGH_CONF_module_gene_level_summary.csv")
)
write_csv(
  exercise_gene_high_conf$genes,
  file.path(TABLE_DIR, "GSE250455_HIGH_CONF_module_gene_level_results.csv")
)

exercise_score <- module_expression_score(
  gse_fit$logCPM,
  gse_fit$meta,
  full_module
)

write_csv(
  exercise_score$summary,
  file.path(TABLE_DIR, "GSE250455_module_sample_score_summary.csv")
)
write_csv(
  exercise_score$paired,
  file.path(TABLE_DIR, "GSE250455_module_sample_score_pairs.csv")
)

save_exercise_module_plot(exercise_score)
save_exercise_gene_plot(exercise_gene$genes)

exercise_perm <- exercise_expression_matched_permutation(
  gse_fit$results,
  full_module,
  n_perm = N_PERM
)

write_csv(
  exercise_perm,
  file.path(TABLE_DIR, "GSE250455_module_expression_matched_permutation.csv")
)

exercise_perm_high_conf <- exercise_expression_matched_permutation(
  gse_fit$results,
  high_conf_for_validation,
  n_perm = N_PERM
)

write_csv(
  exercise_perm_high_conf,
  file.path(TABLE_DIR, "GSE250455_HIGH_CONF_expression_matched_permutation.csv")
)

exercise_gsea <- exercise_fgsea(
  gse_fit$results,
  full_module
)

write_csv(
  exercise_gsea,
  file.path(TABLE_DIR, "GSE250455_module_fgsea.csv")
)

# ---- Targeted genes of interest ----
targeted <- gse_fit$results %>%
  filter(symbol %in% TARGET_GENES) %>%
  mutate(
    target_group = case_when(
      symbol %in% c("PGM1", "TPI1", "PGK1", "ENO3", "LDHA", "FBP2", "PRKAA2") ~
        "Original 7-gene ALS core",
      symbol == "PDK4" ~ "PDK4",
      symbol == "PVALB" ~ "PVALB",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(target_group, FDR, p_value)

write_csv(
  targeted,
  file.path(TABLE_DIR, "GSE250455_targeted_ALS_metabolic_genes.csv")
)

# =============================================================================
# 13. COMPACT CROSS-ANALYSIS SUMMARY
# =============================================================================

cross_summary <- tibble(
  item = c(
    "Directional reversible SOD1-unloading module genes",
    "Statistically supported module genes",
    "Genes reversing in both recovery models",
    "High-confidence genes (statistically supported + both recovery models)",
    "Module convergence permutation P",
    "SOD1 90d LOSO significant expected-direction recurrence",
    "SOD1 120d LOSO significant expected-direction recurrence",
    "SOD1 150d LOSO significant expected-direction recurrence",

    "GSE250455 module genes measured",
    "GSE250455 proportion moving in predicted recovery direction",
    "GSE250455 proportion moving in opposite/ALS direction",
    "GSE250455 sign-test P predicted recovery",
    "GSE250455 sign-test P opposite direction",
    "GSE250455 sign-test P two-sided",
    "GSE250455 paired module-score t-test P predicted recovery",
    "GSE250455 paired module-score t-test P opposite direction",
    "GSE250455 paired module-score t-test P two-sided",
    "GSE250455 paired module-score Wilcoxon P predicted recovery",
    "GSE250455 paired module-score Wilcoxon P opposite direction",
    "GSE250455 paired module-score Wilcoxon P two-sided",
    "GSE250455 expression-matched permutation P predicted recovery",
    "GSE250455 expression-matched permutation P opposite direction",
    "GSE250455 expression-matched permutation P two-sided",

    "GSE250455 HIGH-CONF genes measured",
    "GSE250455 HIGH-CONF proportion predicted recovery",
    "GSE250455 HIGH-CONF sign-test P opposite direction",
    "GSE250455 HIGH-CONF permutation P opposite direction"
  ),
  value = c(
    nrow(directional_module),
    nrow(supported_module),
    nrow(two_recovery_module),
    nrow(high_conf_module),
    perm$summary$empirical_p_one_sided,
    loso$summary$expected_direction_recurrence_fraction[match("90", loso$summary$stage)],
    loso$summary$expected_direction_recurrence_fraction[match("120", loso$summary$stage)],
    loso$summary$expected_direction_recurrence_fraction[match("150", loso$summary$stage)],

    exercise_gene$summary$n_module_genes_measured,
    exercise_gene$summary$proportion_predicted_recovery,
    exercise_gene$summary$proportion_opposite_ALS_direction,
    exercise_gene$summary$sign_test_p_predicted_recovery,
    exercise_gene$summary$sign_test_p_opposite_direction,
    exercise_gene$summary$sign_test_p_two_sided,
    exercise_score$summary$paired_t_p_predicted_recovery,
    exercise_score$summary$paired_t_p_opposite_direction,
    exercise_score$summary$paired_t_p_two_sided,
    exercise_score$summary$paired_wilcox_p_predicted_recovery,
    exercise_score$summary$paired_wilcox_p_opposite_direction,
    exercise_score$summary$paired_wilcox_p_two_sided,
    exercise_perm$empirical_p_predicted_recovery,
    exercise_perm$empirical_p_opposite_direction,
    exercise_perm$empirical_p_two_sided,

    exercise_gene_high_conf$summary$n_module_genes_measured,
    exercise_gene_high_conf$summary$proportion_predicted_recovery,
    exercise_gene_high_conf$summary$sign_test_p_opposite_direction,
    exercise_perm_high_conf$empirical_p_opposite_direction
  )
)

write_csv(
  cross_summary,
  file.path(TABLE_DIR, "PRIMARY_cross_analysis_summary.csv")
)

# =============================================================================
# 14. EXCEL WORKBOOK
# =============================================================================

wb_file <- file.path(OUTDIR, "reversible_ALS_unloading_results.xlsx")

wb_sheets <- list(
  PRIMARY_summary = cross_summary,
  Module_counts = module_counts,
  Module_permutation = perm$summary,
  Module_directional = directional_module,
  Module_supported = supported_module,
  Module_two_recovery = two_recovery_module,
  Module_high_confidence = high_conf_module,
  Module_threshold_sensitivity = sens,
  SOD1_LOSO_summary = loso$summary,
  SOD1_LOSO_gene_results = loso$genes,

  GSE250455_DE = gse_fit$results,
  GSE250455_gene_summary = exercise_gene$summary,
  GSE250455_module_genes = exercise_gene$genes,
  GSE250455_score_summary = exercise_score$summary,
  GSE250455_score_pairs = exercise_score$paired,
  GSE250455_permutation = exercise_perm,
  GSE250455_HC_gene_summary = exercise_gene_high_conf$summary,
  GSE250455_HC_genes = exercise_gene_high_conf$genes,
  GSE250455_HC_permutation = exercise_perm_high_conf,
  GSE250455_fgsea = exercise_gsea,
  GSE250455_targeted = targeted
)

openxlsx::write.xlsx(
  wb_sheets,
  wb_file,
  overwrite = TRUE
)

# =============================================================================
# 15. SESSION INFO
# =============================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(QC_DIR, "sessionInfo.txt")
)

message("\n============================================================")
message("ANALYSIS COMPLETE")
message("============================================================")
message("Primary outputs:")
message("  ", file.path(TABLE_DIR, "reversible_ALS_unloading_module_DIRECTIONAL.csv"))
message("  ", file.path(TABLE_DIR, "module_convergence_permutation_test.csv"))
message("  ", file.path(TABLE_DIR, "SOD1_LOSO_module_trajectory_summary.csv"))
message("  ", file.path(TABLE_DIR, "GSE250455_edgeR_paired_Post_vs_Pre.csv"))
message("  ", file.path(TABLE_DIR, "GSE250455_module_gene_level_summary.csv"))
message("  ", file.path(TABLE_DIR, "GSE250455_module_sample_score_summary.csv"))
message("  ", file.path(TABLE_DIR, "PRIMARY_cross_analysis_summary.csv"))
message("  ", wb_file)
message("============================================================")
