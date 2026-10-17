# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

###############################################################################
# Denervation hypothesis analysis for the predefined ALS-unloading module
###############################################################################

# Purpose:
#   Test whether denervation shifts the frozen, directional 174-gene
#   ALS-unloading module in the ALS-associated direction in two independent
#   GEO datasets, and whether electrical stimulation moves it in reverse.
#   The module is never redefined using these datasets.
#
# Place the downloaded GEO input files in DATA_DIR. Results are written
# beneath DATA_DIR/denervation_hypothesis_analysis.
#
# GSE217576: six-sample expression matrix with no raw-count column;
#             log2(x + 1), limma.
# GSE117440: rat expected counts; edgeR filtering/TMM and limma-voom. Fractional
#             expected counts are retained and are not rounded.
#
# CAMERA and fgsea are the primary directional gene-set tests. Gene-level
# concordance and sample-level module scores are supporting summaries; gene
# signs are not treated as independent biological replicates.

DATA_DIR <- getwd()
OUT_DIR <- file.path(DATA_DIR, "denervation_hypothesis_analysis")
MODULE_FILE <- file.path(
  DATA_DIR, "reversible_ALS_unloading_results", "tables",
  "reversible_ALS_unloading_module_DIRECTIONAL.csv"
)
GSE217576_FILE <- file.path(DATA_DIR, "GSE217576_bulk_RNA-seq_DENVvsCTL.xlsx")
GSE117440_DIR <- file.path(DATA_DIR, "GSE117440")

INSTALL_MISSING <- TRUE
MIN_ORTHOLOG_SUPPORT <- 3L
MIN_ABS_LFC <- 0.10
MIN_GENESET_SIZE <- 10L
N_PERM_SIMPLE <- 10000L
SEED <- 20261013L
options(timeout = max(600, getOption("timeout")))
options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(SEED)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(NULL))
  if (!INSTALL_MISSING) {
    stop("Install missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing)
  }
}

install_if_missing(
  c("readxl", "readr", "dplyr", "tibble", "ggplot2", "openxlsx", "babelgene")
)
install_if_missing(
  c("limma", "edgeR", "fgsea", "AnnotationDbi", "org.Mm.eg.db", "org.Rn.eg.db"),
  bioc = TRUE
)
suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tibble)
  library(ggplot2); library(openxlsx); library(limma); library(edgeR)
  library(fgsea); library(AnnotationDbi); library(org.Mm.eg.db)
  library(org.Rn.eg.db); library(babelgene)
})

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")
QC_DIR <- file.path(OUT_DIR, "QC")
for (p in c(TABLE_DIR, FIGURE_DIR, QC_DIR)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

assert(dir.exists(DATA_DIR), paste("Data folder not found:", DATA_DIR))
for (f in c(GSE217576_FILE, MODULE_FILE)) {
  assert(file.exists(f), paste("Required input file not found:", f))
}
assert(dir.exists(GSE117440_DIR), paste("Dataset folder not found:", GSE117440_DIR))

write_csv <- function(x, path) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    if (is.list(x[[nm]])) {
      x[[nm]] <- vapply(x[[nm]], function(z) paste(as.character(z), collapse = ";"), character(1))
    }
  }
  readr::write_csv(x, path, na = "")
  invisible(path)
}

# Load the frozen human-orthologue module from the prior analysis.
module_raw <- readr::read_csv(MODULE_FILE, show_col_types = FALSE)
assert(all(c("symbol", "als_direction") %in% names(module_raw)),
       "Module file must contain symbol and als_direction columns.")
if ("directional_module" %in% names(module_raw)) {
  keep <- tolower(as.character(module_raw$directional_module)) %in% c("true", "1", "yes")
  module_raw <- module_raw[keep, , drop = FALSE]
}
module <- data.frame(
  symbol = toupper(trimws(as.character(module_raw$symbol))),
  als_direction = suppressWarnings(as.integer(module_raw$als_direction)),
  stringsAsFactors = FALSE
)
module <- module[!is.na(module$symbol) & module$symbol != "", , drop = FALSE]
module$component <- ifelse(module$als_direction < 0, "ALS_DOWN", "ALS_UP")
assert(
  nrow(module) == 174L && !anyDuplicated(module$symbol) &&
    sum(module$als_direction == -1L) == 128L &&
    sum(module$als_direction == 1L) == 46L,
  "Expected the frozen module with 174 genes: 128 ALS_DOWN and 46 ALS_UP."
)

clean_token <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\|.*$", "", x)
  x <- sub("\\.\\d+$", "", x)
  x[x == ""] <- NA_character_
  x
}

build_symbol_lookup <- function(orgdb) {
  official <- as.character(AnnotationDbi::keys(orgdb, keytype = "SYMBOL"))
  official_lookup <- setNames(official, toupper(official))
  alias_lookup <- character()
  if ("ALIAS" %in% AnnotationDbi::keytypes(orgdb)) {
    alias_keys <- AnnotationDbi::keys(orgdb, keytype = "ALIAS")
    alias_tbl <- suppressMessages(AnnotationDbi::select(
      orgdb, keys = alias_keys, keytype = "ALIAS", columns = "SYMBOL"
    ))
    alias_tbl <- alias_tbl[
      !is.na(alias_tbl$ALIAS) & !is.na(alias_tbl$SYMBOL),
      c("ALIAS", "SYMBOL"), drop = FALSE
    ]
    alias_tbl$key <- toupper(as.character(alias_tbl$ALIAS))
    alias_tbl <- alias_tbl[order(alias_tbl$key, alias_tbl$SYMBOL), , drop = FALSE]
    alias_tbl <- alias_tbl[!duplicated(alias_tbl$key), , drop = FALSE]
    alias_lookup <- setNames(as.character(alias_tbl$SYMBOL), alias_tbl$key)
  }
  list(official = official_lookup, alias = alias_lookup)
}

restore_official_symbols <- function(x, lookup) {
  token <- toupper(clean_token(x))
  result <- unname(lookup$official[token])
  needs_alias <- is.na(result) & !is.na(token)
  if (any(needs_alias) && length(lookup$alias)) {
    result[needs_alias] <- unname(lookup$alias[token[needs_alias]])
  }
  result
}

map_ids_to_symbol <- function(ids, orgdb, keytype) {
  ids <- clean_token(ids)
  result <- rep(NA_character_, length(ids))
  available <- as.character(AnnotationDbi::keys(orgdb, keytype = keytype))
  valid <- unique(ids[!is.na(ids) & ids %in% available])
  if (length(valid)) {
    mapped <- AnnotationDbi::mapIds(
      orgdb, keys = valid, keytype = keytype, column = "SYMBOL", multiVals = "first"
    )
    use <- !is.na(ids) & ids %in% valid
    result[use] <- as.character(mapped[ids[use]])
  }
  result
}

make_ortholog_lookup <- function(source_symbols, species) {
  source_symbols <- unique(as.character(source_symbols))
  source_symbols <- source_symbols[!is.na(source_symbols) & source_symbols != ""]
  assert(length(source_symbols) > 0L, paste("No valid", species, "symbols to map."))
  orth <- babelgene::orthologs(
    genes = source_symbols, species = species, human = FALSE,
    min_support = MIN_ORTHOLOG_SUPPORT, top = TRUE
  )
  assert(nrow(orth) > 0L, paste("No supported", species, "to human orthologues found."))
  assert(all(c("symbol", "human_symbol", "support_n") %in% names(orth)),
         "Unexpected columns returned by babelgene::orthologs().")
  out <- data.frame(
    source_key = toupper(as.character(orth$symbol)),
    human_symbol = toupper(as.character(orth$human_symbol)),
    ortholog_support_n = suppressWarnings(as.integer(orth$support_n)),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$source_key) & !is.na(out$human_symbol) & out$human_symbol != "",
             , drop = FALSE]
  out <- out[order(out$source_key, -out$ortholog_support_n), , drop = FALSE]
  out[!duplicated(out$source_key), , drop = FALSE]
}

harmonize_to_human <- function(expr, feature_ids, source_symbols, species, weights = NULL) {
  assert(nrow(expr) == length(feature_ids), "Feature ID count does not match matrix rows.")
  assert(nrow(expr) == length(source_symbols), "Symbol count does not match matrix rows.")
  assert(all(is.finite(expr)), "Expression matrix contains non-finite values.")
  feature_tbl <- data.frame(
    row_index = seq_len(nrow(expr)),
    feature_id = as.character(feature_ids),
    source_symbol = as.character(source_symbols),
    source_key = toupper(as.character(source_symbols)),
    mean_expression = rowMeans(expr),
    stringsAsFactors = FALSE
  )
  orth <- make_ortholog_lookup(source_symbols, species)
  mapping <- merge(feature_tbl, orth, by = "source_key", all.x = TRUE, sort = FALSE)
  mapping$status <- ifelse(
    is.na(mapping$source_symbol) | mapping$source_symbol == "",
    "no_source_symbol",
    ifelse(is.na(mapping$human_symbol), "no_supported_human_orthologue", "candidate")
  )
  candidate <- which(mapping$status == "candidate")
  assert(length(candidate) > 0L, paste("No", species, "features mapped to human orthologues."))
  ordered <- candidate[order(
    mapping$source_key[candidate], -mapping$mean_expression[candidate],
    mapping$feature_id[candidate]
  )]
  duplicate_source <- duplicated(mapping$source_key[ordered])
  mapping$status[ordered[duplicate_source]] <- "duplicate_source_symbol_excluded"
  selected <- ordered[!duplicate_source]
  human_symbols <- mapping$human_symbol[selected]
  duplicate_human <- duplicated(human_symbols) | duplicated(human_symbols, fromLast = TRUE)
  mapping$status[selected[duplicate_human]] <- "duplicate_human_orthologue_excluded"
  selected <- selected[!duplicate_human]
  mapping$status[selected] <- "retained_one_to_one"
  assert(length(selected) > 0L, "No one-to-one human orthologues remain after mapping.")
  mapping <- mapping[order(mapping$row_index), c(
    "row_index", "feature_id", "source_symbol", "human_symbol",
    "ortholog_support_n", "mean_expression", "status"
  ), drop = FALSE]
  kept <- mapping[mapping$status == "retained_one_to_one", , drop = FALSE]
  kept <- kept[order(kept$row_index), , drop = FALSE]
  expr_human <- expr[kept$row_index, , drop = FALSE]
  rownames(expr_human) <- kept$human_symbol
  weights_human <- NULL
  if (!is.null(weights)) {
    weights_human <- weights[kept$row_index, , drop = FALSE]
    rownames(weights_human) <- kept$human_symbol
  }
  list(expr = expr_human, weights = weights_human, mapping_qc = mapping)
}

numeric_matrix <- function(x, sample_cols, row_ids) {
  m <- do.call(cbind, lapply(sample_cols, function(nm) as.numeric(x[[nm]])))
  colnames(m) <- sample_cols
  rownames(m) <- as.character(row_ids)
  storage.mode(m) <- "numeric"
  assert(all(is.finite(m)), "Input expression matrix contains missing/non-finite values.")
  assert(all(m >= 0), "Expression values must be non-negative.")
  m
}

make_metadata <- function(samples, groups, levels) {
  assert(length(samples) == length(groups), "Sample and group vector lengths differ.")
  assert(!anyNA(groups), "Could not assign a group to every sample.")
  data.frame(
    sample = as.character(samples),
    group = factor(as.character(groups), levels = levels),
    stringsAsFactors = FALSE
  )
}

read_gse217576 <- function(mm_lookup) {
  message("Reading GSE217576 workbook...")
  x <- as.data.frame(readxl::read_excel(
    GSE217576_FILE, sheet = 1, .name_repair = "minimal"
  ), check.names = FALSE)
  required <- c("gene id", "gene name", "CTL1", "CTL2", "CTL3", "DENV1", "DENV2", "DENV3")
  assert(all(required %in% names(x)), "GSE217576 workbook columns differ from expected layout.")
  samples <- c("CTL1", "CTL2", "CTL3", "DENV1", "DENV2", "DENV3")
  expr <- numeric_matrix(x, samples, x[["gene id"]])
  symbols <- restore_official_symbols(x[["gene name"]], mm_lookup)
  needs_id <- is.na(symbols)
  if (any(needs_id)) {
    ens <- map_ids_to_symbol(x[["gene id"]][needs_id], org.Mm.eg.db, "ENSEMBL")
    symbols[needs_id] <- restore_official_symbols(ens, mm_lookup)
  }
  groups <- ifelse(grepl("^CTL", samples), "CTL", "DENV")
  meta <- make_metadata(samples, groups, c("CTL", "DENV"))
  assert(sum(meta$group == "CTL") == 3L && sum(meta$group == "DENV") == 3L,
         "GSE217576 must contain three controls and three denervated samples.")
  list(
    name = "GSE217576", species = "mouse", expr = log2(expr + 1),
    feature_ids = rownames(expr), source_symbols = symbols, meta = meta,
    contrast_defs = data.frame(
      contrast = "DEN_vs_CTL", expression = "DENV-CTL",
      numerator = "DENV", denominator = "CTL", phase = "denervation",
      stringsAsFactors = FALSE
    )
  )
}

read_gse117440 <- function() {
  message("Reading GSE117440 per-sample expected-count files...")
  files <- sort(list.files(
    GSE117440_DIR,
    pattern = "^GSM[0-9]+_(SHAM|DN|DN-SM)[0-9]+\\.gene\\.FPKM\\.txt$",
    full.names = TRUE,
    ignore.case = TRUE
  ))
  assert(length(files) == 9L, "GSE117440 should contain nine SHAM, DN, and DN-SM files.")
  samples <- sub("\\.gene\\.FPKM\\.txt$", "", basename(files), ignore.case = TRUE)
  groups <- ifelse(
    grepl("_SHAM", samples, ignore.case = TRUE), "SHAM",
    ifelse(grepl("_DN-SM", samples, ignore.case = TRUE), "DN_SM", "DEN")
  )
  meta <- make_metadata(samples, groups, c("SHAM", "DEN", "DN_SM"))
  observed_n <- table(meta$group)
  assert(all(as.integer(observed_n[c("SHAM", "DEN", "DN_SM")]) == c(3L, 3L, 3L)),
         "GSE117440 must contain three samples in each SHAM, DN, and DN-SM group.")
  file_data <- lapply(files, function(f) {
    z <- as.data.frame(readr::read_tsv(
      f, show_col_types = FALSE, progress = FALSE
    ), check.names = FALSE)
    assert(all(c("gene_id", "expected_count", "FPKM") %in% names(z)),
           paste("Unexpected GSE117440 file layout:", basename(f)))
    z$gene_id <- clean_token(z$gene_id)
    z$expected_count <- as.numeric(z$expected_count)
    z <- z[!is.na(z$gene_id) & is.finite(z$expected_count), , drop = FALSE]
    if (anyDuplicated(z$gene_id)) {
      z <- stats::aggregate(expected_count ~ gene_id, data = z, FUN = sum)
    } else {
      z <- z[, c("gene_id", "expected_count"), drop = FALSE]
    }
    z
  })
  all_ids <- sort(unique(unlist(lapply(file_data, function(z) z$gene_id), use.names = FALSE)))
  counts <- matrix(0, nrow = length(all_ids), ncol = length(files),
                   dimnames = list(all_ids, samples))
  for (j in seq_along(file_data)) {
    z <- file_data[[j]]
    counts[z$gene_id, j] <- z$expected_count
  }
  assert(all(is.finite(counts)) && all(counts >= 0), "GSE117440 counts are invalid.")
  list(
    name = "GSE117440", species = "rat", counts = counts,
    feature_ids = rownames(counts),
    meta = meta,
    contrast_defs = data.frame(
      contrast = c("DEN_vs_SHAM", "Electrical_stimulation_vs_DEN"),
      expression = c("DEN-SHAM", "DN_SM-DEN"),
      numerator = c("DEN", "DN_SM"), denominator = c("SHAM", "DEN"),
      phase = c("denervation", "intervention"),
      stringsAsFactors = FALSE
    )
  )
}

###############################################################################
# 4. Differential expression and frozen-module tests
###############################################################################

build_design <- function(meta) {
  group_levels <- levels(meta$group)
  assert(length(group_levels) >= 2L, "At least two groups are required.")
  design <- stats::model.matrix(~ 0 + group, data = meta)
  colnames(design) <- group_levels
  design
}

expected_component_direction <- function(component, phase) {
  if (phase == "denervation") {
    if (component == "ALS_DOWN") "Down" else "Up"
  } else {
    if (component == "ALS_DOWN") "Up" else "Down"
  }
}

expected_nes_sign <- function(component, phase) {
  if (phase == "denervation") {
    if (component == "ALS_DOWN") -1 else 1
  } else {
    if (component == "ALS_DOWN") 1 else -1
  }
}

test_module_score <- function(score, meta, def, dataset) {
  numerator_values <- score[meta$group == def$numerator]
  denominator_values <- score[meta$group == def$denominator]
  tt <- tryCatch(
    stats::t.test(numerator_values, denominator_values, var.equal = FALSE),
    error = function(e) NULL
  )
  effect <- mean(numerator_values, na.rm = TRUE) - mean(denominator_values, na.rm = TRUE)
  expected_sign <- if (def$phase == "denervation") 1 else -1
  data.frame(
    dataset = dataset,
    contrast = def$contrast,
    phase = def$phase,
    numerator = def$numerator,
    denominator = def$denominator,
    n_numerator = length(numerator_values),
    n_denominator = length(denominator_values),
    mean_score_numerator = mean(numerator_values, na.rm = TRUE),
    mean_score_denominator = mean(denominator_values, na.rm = TRUE),
    score_difference = effect,
    expected_score_sign = expected_sign,
    direction_matches_expectation = is.finite(effect) && sign(effect) == expected_sign,
    p_value = if (is.null(tt)) NA_real_ else tt$p.value,
    ci_low = if (is.null(tt)) NA_real_ else unname(tt$conf.int[1]),
    ci_high = if (is.null(tt)) NA_real_ else unname(tt$conf.int[2]),
    stringsAsFactors = FALSE
  )
}

analyze_dataset <- function(ds, harmonized) {
  dataset <- ds$name
  expr <- harmonized$expr
  weights <- harmonized$weights
  meta <- ds$meta
  design <- build_design(meta)
  contrast_matrix <- limma::makeContrasts(
    contrasts = ds$contrast_defs$expression,
    levels = design
  )
  colnames(contrast_matrix) <- ds$contrast_defs$contrast

  fit0 <- if (is.null(weights)) {
    limma::lmFit(expr, design)
  } else {
    limma::lmFit(expr, design, weights = weights)
  }
  fit <- limma::contrasts.fit(fit0, contrast_matrix)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  module_genes <- split(module$symbol, module$component)
  module_indices <- lapply(module_genes, function(g) which(rownames(expr) %in% g))
  module_indices <- module_indices[lengths(module_indices) >= MIN_GENESET_SIZE]

  score_symbols <- intersect(module$symbol, rownames(expr))
  assert(length(score_symbols) >= 20L, paste(dataset, "has too few module genes for a score."))
  score_module <- module[match(score_symbols, module$symbol), , drop = FALSE]
  score_expr <- expr[score_symbols, , drop = FALSE]
  z <- t(scale(t(score_expr)))
  z[!is.finite(z)] <- 0
  signed_z <- sweep(z, 1, score_module$als_direction, function(a, b) a * b)
  sample_score <- colMeans(signed_z)
  sample_scores <- data.frame(
    dataset = dataset,
    sample = meta$sample,
    group = as.character(meta$group),
    module_score = as.numeric(sample_score[meta$sample]),
    stringsAsFactors = FALSE
  )

  de_tables <- list()
  module_effects <- list()
  camera_tables <- list()
  fgsea_tables <- list()
  score_tests <- list()

  for (j in seq_len(nrow(ds$contrast_defs))) {
    def <- ds$contrast_defs[j, , drop = FALSE]
    contrast_name <- def$contrast
    tt <- limma::topTable(
      fit, coef = contrast_name, number = Inf,
      sort.by = "none", adjust.method = "BH"
    )
    de <- data.frame(
      symbol = rownames(tt),
      log2FC = as.numeric(tt$logFC),
      average_expression = as.numeric(tt$AveExpr),
      moderated_t = as.numeric(tt$t),
      p_value = as.numeric(tt$P.Value),
      FDR = as.numeric(tt$adj.P.Val),
      B = as.numeric(tt$B),
      dataset = dataset,
      contrast = contrast_name,
      phase = def$phase,
      stringsAsFactors = FALSE
    )
    de_tables[[contrast_name]] <- de
    write_csv(de, file.path(TABLE_DIR, paste0("DE_", dataset, "_", contrast_name, ".csv")))

    me <- merge(module, de, by = "symbol", all = FALSE, sort = FALSE)
    # `def$phase` is a single value. Using it as the test argument to
    # ifelse() would return a one-element result and recycle the first gene's
    # sign across the whole module. Branch on phase, then retain the per-gene
    # direction vector.
    if (def$phase == "denervation") {
      me$expected_effect_sign <- me$als_direction
    } else if (def$phase == "intervention") {
      me$expected_effect_sign <- -me$als_direction
    } else {
      stop("Unexpected contrast phase: ", def$phase, call. = FALSE)
    }
    me$effect_status <- ifelse(
      abs(me$log2FC) < MIN_ABS_LFC,
      "small_effect",
      ifelse(sign(me$log2FC) == me$expected_effect_sign, "concordant", "discordant")
    )
    module_effects[[contrast_name]] <- me

    if (length(module_indices)) {
      camera_args <- list(
        y = expr,
        index = module_indices,
        design = design,
        contrast = contrast_matrix[, contrast_name],
        use.ranks = FALSE
      )
      if (!is.null(weights)) camera_args$weights <- weights
      cam <- as.data.frame(do.call(limma::camera, camera_args))
      cam$gene_set <- rownames(cam)
      cam$dataset <- dataset
      cam$contrast <- contrast_name
      cam$phase <- def$phase
      cam$expected_direction <- vapply(
        cam$gene_set, expected_component_direction, character(1), phase = def$phase
      )
      cam$direction_matches_expectation <-
        tolower(as.character(cam$Direction)) == tolower(cam$expected_direction)
      camera_tables[[contrast_name]] <- cam
    }

    pathway_sets <- lapply(module_genes, function(g) intersect(g, rownames(expr)))
    pathway_sets <- pathway_sets[lengths(pathway_sets) >= MIN_GENESET_SIZE]
    if (length(pathway_sets)) {
      rank_stat <- fit$t[, contrast_name]
      names(rank_stat) <- rownames(expr)
      rank_stat <- sort(rank_stat[is.finite(rank_stat)], decreasing = TRUE)
      fg <- as.data.frame(fgsea::fgseaMultilevel(
        pathways = pathway_sets,
        stats = rank_stat,
        minSize = MIN_GENESET_SIZE,
        maxSize = Inf,
        nPermSimple = N_PERM_SIMPLE,
        eps = 1e-10
      ))
      if (nrow(fg)) {
        fg$leadingEdge <- vapply(
          fg$leadingEdge, function(z) paste(as.character(z), collapse = ";"), character(1)
        )
        fg$dataset <- dataset
        fg$contrast <- contrast_name
        fg$phase <- def$phase
        fg$expected_NES_sign <- vapply(
          as.character(fg$pathway), expected_nes_sign, numeric(1), phase = def$phase
        )
        fg$NES_matches_expectation <- sign(fg$NES) == fg$expected_NES_sign
        fgsea_tables[[contrast_name]] <- fg
      }
    }

    score_tests[[contrast_name]] <- test_module_score(
      sample_score[meta$sample], meta, def, dataset
    )
  }

  list(
    dataset = dataset,
    expr = expr,
    meta = meta,
    de = de_tables,
    module_effects = module_effects,
    camera = camera_tables,
    fgsea = fgsea_tables,
    sample_scores = sample_scores,
    score_tests = score_tests
  )
}

###############################################################################
# 5. Prepare all datasets and fit the primary and intervention contrasts
###############################################################################

mm_lookup <- build_symbol_lookup(org.Mm.eg.db)
rn_lookup <- build_symbol_lookup(org.Rn.eg.db)

ds217576 <- read_gse217576(mm_lookup)
ds117440 <- read_gse117440()

write_csv(ds217576$meta, file.path(QC_DIR, "GSE217576_sample_metadata.csv"))
write_csv(ds117440$meta, file.path(QC_DIR, "GSE117440_sample_metadata.csv"))

harm217576 <- harmonize_to_human(
  expr = ds217576$expr,
  feature_ids = ds217576$feature_ids,
  source_symbols = ds217576$source_symbols,
  species = "mouse"
)
write_csv(harm217576$mapping_qc, file.path(QC_DIR, "GSE217576_orthologue_mapping.csv"))

result217576 <- analyze_dataset(ds217576, harm217576)

# GSE117440 provides fractional expected counts. Keep them as supplied, filter
# low-expression genes, apply TMM normalization, and use voom sample weights.
design117440 <- build_design(ds117440$meta)
y117440 <- edgeR::DGEList(counts = ds117440$counts, group = ds117440$meta$group)
keep117440 <- edgeR::filterByExpr(y117440, group = ds117440$meta$group)
assert(sum(keep117440) > 1000L, "GSE117440 filtering retained unexpectedly few genes.")
y117440 <- y117440[keep117440, , keep.lib.sizes = FALSE]
y117440 <- edgeR::calcNormFactors(y117440, method = "TMM")
voom117440 <- limma::voom(y117440, design = design117440, plot = FALSE)
rat_symbols <- restore_official_symbols(
  map_ids_to_symbol(rownames(voom117440$E), org.Rn.eg.db, "ENTREZID"),
  rn_lookup
)
harm117440 <- harmonize_to_human(
  expr = voom117440$E,
  feature_ids = rownames(voom117440$E),
  source_symbols = rat_symbols,
  species = "rat",
  weights = voom117440$weights
)
write_csv(harm117440$mapping_qc, file.path(QC_DIR, "GSE117440_orthologue_mapping.csv"))
result117440 <- analyze_dataset(ds117440, harm117440)
results <- list(result217576, result117440)

###############################################################################
# 6. Combined summaries, figures, and output files
###############################################################################

all_module_effects <- dplyr::bind_rows(
  unlist(lapply(results, function(x) x$module_effects), recursive = FALSE)
)
all_scores <- dplyr::bind_rows(lapply(results, function(x) x$sample_scores))
all_score_tests <- dplyr::bind_rows(
  unlist(lapply(results, function(x) x$score_tests), recursive = FALSE)
)
all_camera <- dplyr::bind_rows(
  unlist(lapply(results, function(x) x$camera), recursive = FALSE)
)
all_fgsea <- dplyr::bind_rows(
  unlist(lapply(results, function(x) x$fgsea), recursive = FALSE)
)

if (nrow(all_camera)) {
  all_camera$FDR_all_tests <- p.adjust(all_camera$PValue, method = "BH")
}
if (nrow(all_fgsea)) {
  all_fgsea$FDR_all_tests <- p.adjust(all_fgsea$pval, method = "BH")
}
if (nrow(all_score_tests)) {
  all_score_tests$FDR_all_tests <- p.adjust(all_score_tests$p_value, method = "BH")
}

# Count directional agreement at the gene level as a descriptive summary.
effect_rows <- dplyr::bind_rows(
  all_module_effects %>% dplyr::mutate(module_component = component),
  all_module_effects %>% dplyr::mutate(module_component = "ALL")
)
concordance <- effect_rows %>%
  dplyr::group_by(dataset, contrast, phase, module_component) %>%
  dplyr::summarise(
    n_evaluable = sum(is.finite(log2FC)),
    n_concordant_at_abs_log2FC_0_10 =
      sum(effect_status == "concordant", na.rm = TRUE),
    n_discordant_at_abs_log2FC_0_10 =
      sum(effect_status == "discordant", na.rm = TRUE),
    n_small_effect = sum(effect_status == "small_effect", na.rm = TRUE),
    n_sign_concordant_without_threshold = sum(
      is.finite(log2FC) & sign(log2FC) == expected_effect_sign
    ),
    n_concordant_FDR05 = sum(
      effect_status == "concordant" & FDR < 0.05,
      na.rm = TRUE
    ),
    pct_concordant_at_abs_log2FC_0_10 =
      100 * n_concordant_at_abs_log2FC_0_10 / pmax(n_evaluable, 1),
    pct_sign_concordant_without_threshold =
      100 * n_sign_concordant_without_threshold / pmax(n_evaluable, 1),
    .groups = "drop"
  )

coverage <- dplyr::bind_rows(lapply(results, function(res) {
  present <- rownames(res$expr)
  dplyr::bind_rows(lapply(c("ALS_DOWN", "ALS_UP", "ALL"), function(comp) {
    expected <- if (comp == "ALL") {
      module$symbol
    } else {
      module$symbol[module$component == comp]
    }
    measured <- intersect(expected, present)
    data.frame(
      dataset = res$dataset,
      module_component = comp,
      n_expected = length(expected),
      n_measured = length(measured),
      n_missing = length(setdiff(expected, measured)),
      percent_measured = 100 * length(measured) / length(expected),
      missing_symbols = paste(setdiff(expected, measured), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}))

write_csv(all_module_effects, file.path(TABLE_DIR, "module_gene_effects_by_contrast.csv"))
write_csv(concordance, file.path(TABLE_DIR, "module_directional_concordance.csv"))
write_csv(coverage, file.path(TABLE_DIR, "module_gene_coverage.csv"))
write_csv(all_camera, file.path(TABLE_DIR, "CAMERA_directional_module_tests.csv"))
write_csv(all_fgsea, file.path(TABLE_DIR, "fgsea_directional_module_tests.csv"))
write_csv(all_scores, file.path(TABLE_DIR, "sample_level_module_scores.csv"))
write_csv(all_score_tests, file.path(TABLE_DIR, "sample_level_module_score_contrasts.csv"))

score_plot <- ggplot2::ggplot(
  all_scores,
  ggplot2::aes(x = group, y = module_score, color = group)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = group),
    outlier.shape = NA,
    alpha = 0.15,
    width = 0.58
  ) +
  ggplot2::geom_jitter(width = 0.08, size = 2.2, alpha = 0.9) +
  ggplot2::facet_wrap(~dataset, scales = "free_x") +
  ggplot2::labs(
    x = NULL,
    y = "ALS-direction module score (mean signed gene z-score)",
    title = "Denervation and intervention effects on the frozen ALS-unloading module"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )
ggplot2::ggsave(
  file.path(FIGURE_DIR, "module_scores_by_dataset.png"),
  score_plot, width = 10, height = 6, dpi = 300
)
ggplot2::ggsave(
  file.path(FIGURE_DIR, "module_scores_by_dataset.pdf"),
  score_plot, width = 10, height = 6
)

# PCA plots use the 500 most variable retained human-orthologue genes in each
# study. They are study-specific QC views; samples are never pooled.
for (res in results) {
  expr <- res$expr
  vars <- apply(expr, 1, stats::var)
  top <- order(vars, decreasing = TRUE)
  top <- top[seq_len(min(500L, length(top)))]
  pc <- stats::prcomp(t(expr[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  pca_data <- data.frame(
    sample = rownames(pc$x),
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    group = as.character(res$meta$group[match(rownames(pc$x), res$meta$sample)]),
    stringsAsFactors = FALSE
  )
  pca_plot <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = PC1, y = PC2, color = group, label = sample)
  ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(
      check_overlap = TRUE, vjust = -0.8, size = 2.7, show.legend = FALSE
    ) +
    ggplot2::labs(
      title = paste(res$dataset, "PCA"),
      x = sprintf("PC1 (%.1f%%)", pct[1]),
      y = sprintf("PC2 (%.1f%%)", pct[2]),
      color = "Group"
    ) +
    ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(
    file.path(FIGURE_DIR, paste0(res$dataset, "_PCA.png")),
    pca_plot, width = 7, height = 5.5, dpi = 300
  )
  ggplot2::ggsave(
    file.path(FIGURE_DIR, paste0(res$dataset, "_PCA.pdf")),
    pca_plot, width = 7, height = 5.5
  )
}

summary_sheets <- list(
  "Module coverage" = coverage,
  "Direction concordance" = concordance,
  "CAMERA" = all_camera,
  "fgsea" = all_fgsea,
  "Module score tests" = all_score_tests,
  "Sample scores" = all_scores
)
summary_file <- file.path(OUT_DIR, "denervation_hypothesis_summary.xlsx")
summary_wb <- openxlsx::createWorkbook()
for (sheet_name in names(summary_sheets)) {
  sheet_data <- as.data.frame(summary_sheets[[sheet_name]], stringsAsFactors = FALSE)
  openxlsx::addWorksheet(summary_wb, sheet_name)
  openxlsx::writeData(
    summary_wb,
    sheet = sheet_name,
    x = sheet_data,
    colNames = TRUE,
    rowNames = FALSE,
    withFilter = nrow(sheet_data) > 0L
  )
  openxlsx::freezePane(summary_wb, sheet = sheet_name, firstRow = TRUE)
}
openxlsx::saveWorkbook(summary_wb, file = summary_file, overwrite = TRUE)

# Confirm every summary tab contains the expected header and number of rows.
for (sheet_name in names(summary_sheets)) {
  exported <- openxlsx::read.xlsx(summary_file, sheet = sheet_name)
  expected <- as.data.frame(summary_sheets[[sheet_name]], stringsAsFactors = FALSE)
  assert(
    identical(names(exported), names(expected)) && nrow(exported) == nrow(expected),
    paste("Summary workbook export check failed for sheet:", sheet_name)
  )
}

writeLines(
  c(
    "Denervation hypothesis analysis",
    paste("Data folder:", normalizePath(DATA_DIR, winslash = "/", mustWork = TRUE)),
    paste("Frozen module:", normalizePath(MODULE_FILE, winslash = "/", mustWork = TRUE)),
    "",
    "Primary denervation contrasts:",
    "  GSE217576: DENV vs CTL",
    "  GSE117440: DN vs SHAM",
    "",
    "Intervention contrasts:",
    "  GSE117440: electrical stimulation (DN-SM) vs DN",
    "",
    "The module was not reselected from these datasets. Gene-level direction",
    "counts are descriptive; CAMERA, fgsea, and sample-level scores are the",
    "prespecified inferential summaries. Samples from distinct datasets were",
    "not pooled.",
    "",
    "GSE117440 has no sham-plus-stimulation group, so stimulation is compared",
    "with DN only. GSE217576 and GSE117440 were analyzed separately.",
    "",
    capture.output(sessionInfo())
  ),
  con = file.path(QC_DIR, "analysis_notes_and_session_info.txt")
)

message("Analysis complete. Results are in: ", normalizePath(OUT_DIR, winslash = "/"))
