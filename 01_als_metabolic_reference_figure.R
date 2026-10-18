# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# Figure 1 | SOD1G93A ALS metabolic reference
# ------------------------------------------------------------
# Uses the original DAVID KEGG enrichment tables in SOD_data.xlsx:
#   Table S8  = downregulated mRNAs at 90 days
#   Table S10 = downregulated mRNAs at 120 days
#   Table S12 = downregulated mRNAs at 150 days
#
# Figure: dot-matrix of nine KEGG pathways across the three stages.
#   Dot size   = DAVID Count (downregulated genes mapping to pathway)
#   Dot colour = -log10(Benjamini-adjusted enrichment P)
#   Grey dot   = Benjamini-adjusted P >= 0.05 (NOT 'absent')
#   GSEA       = pathway has exact MSigDB gene set used downstream
#
# CRITICAL: Use Benjamini-adjusted enrichment P < 0.05 for recurrence.
# Raw P < 0.05 would incorrectly classify all nine pathways as
# enriched at all three stages in these source tables.
#
# From the directory containing this script and the workbook, run:
#   source("01_als_metabolic_reference_figure.R")
#
# Requires: install.packages(c("readxl", "ggplot2", "ragg", "svglite"))
# Writes 600-dpi PNG, vector SVG, and audit CSV in the working directory.

# --------------- Configuration ----------------------------
INPUT_FILE <- "SOD_data.xlsx"
# The uploaded workbook may have a numbered name; this is a fallback.
if (!file.exists(INPUT_FILE) && file.exists("SOD_data(2).xlsx")) {
  INPUT_FILE <- "SOD_data(2).xlsx"
}
OUTPUT_PREFIX <- "Figure_1_ALS_metabolic_reference"
Q_CUTOFF <- 0.05

for (pkg in c("readxl", "ggplot2", "ragg", "svglite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing R package: ", pkg,
         ". Install with install.packages('", pkg, "').")
  }
}
if (!file.exists(INPUT_FILE)) {
  stop("Workbook not found in ", getwd(),
       ". Please put SOD_data.xlsx in the R working directory.")
}

sheets <- c("90 days" = "Table S8",
            "120 days" = "Table S10",
            "150 days" = "Table S12")
stopifnot(all(unname(sheets) %in% readxl::excel_sheets(INPUT_FILE)))

reference <- data.frame(
  pathway = c(
    "Glycolysis/gluconeogenesis",
    "Fructose and mannose metabolism",
    "Carbon metabolism",
    "AMPK signaling pathway",
    "Insulin signaling pathway",
    "Glucagon signaling pathway",
    "Purine metabolism",
    "Pentose phosphate pathway",
    "Starch and sucrose metabolism"
  ),
  tested_in_GSEA = c(TRUE, TRUE, FALSE, FALSE, TRUE,
                     FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Normalize typography and UK/US spellings, but require an EXACT
# pathway-term match; a partial match would incorrectly conflate
# "Carbon metabolism" with "Central carbon metabolism in cancer".
normalize_term <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("signaling", "signaling", x, fixed = TRUE)
  gsub("[^a-z0-9]", "", x)
}

read_stage <- function(age, sheet) {
  # Row 1 = sheet title; row 2 blank; row 3 = real header.
  tab <- as.data.frame(readxl::read_excel(
    INPUT_FILE, sheet = sheet, skip = 2, .name_repair = "unique"
  ), stringsAsFactors = FALSE)
  names(tab) <- trimws(names(tab))
  need <- c("Category", "Term", "Count", "P-Value",
            "Benjamini", "Fold Enrichment", "FDR")
  missing_columns <- setdiff(need, names(tab))
  if (length(missing_columns)) {
    stop(sheet, " is missing columns: ",
         paste(missing_columns, collapse = ", "))
  }

  # Only actual KEGG rows, not other functional annotation sections.
  tab <- tab[!is.na(tab$Category) &
               trimws(as.character(tab$Category)) == "KEGG_PATHWAY",,
             drop = FALSE]
  tab$normalized_term <- normalize_term(tab$Term)

  results <- lapply(seq_len(nrow(reference)), function(i) {
    name <- reference$pathway[i]
    matches <- which(tab$normalized_term == normalize_term(name))
    if (length(matches) != 1L) {
      stop("Expected exactly one KEGG row for '", name, "' in ",
           sheet, "; found ", length(matches), ".")
    }
    hit <- tab[matches, , drop = FALSE]
    data.frame(
      pathway = name,
      age = age,
      source_sheet = sheet,
      term_in_source = as.character(hit$Term),
      count = as.numeric(hit$Count),
      p_raw = as.numeric(hit[["P-Value"]]),
      p_benjamini = as.numeric(hit$Benjamini),
      p_fdr_david = as.numeric(hit$FDR),
      fold_enrichment = as.numeric(hit[["Fold Enrichment"]]),
      tested_in_GSEA = reference$tested_in_GSEA[i],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

# ---------------- Extract + validate ---------------------
all_data <- do.call(rbind, Map(read_stage, names(sheets),
                               unname(sheets)))
if (nrow(all_data) != 27L ||
    any(!is.finite(all_data$p_benjamini)) ||
    any(!is.finite(all_data$count))) {
  stop("Unexpected missing/non-numeric values in the 9 x 3 KEGG data.")
}
all_data$significant <- all_data$p_benjamini < Q_CUTOFF
all_data$minus_log10_q <- -log10(pmax(all_data$p_benjamini,
                                      .Machine$double.xmin))

# Audit values. These are the exact numbers used to make the plot.
write.csv(all_data,
          paste0(OUTPUT_PREFIX, "_extracted_values.csv"),
          row.names = FALSE)

recurrence <- aggregate(as.integer(all_data$significant),
                        by = list(pathway = all_data$pathway),
                        FUN = sum)
names(recurrence)[2] <- "stages_with_BH_q_below_0.05"
recurrence <- recurrence[match(reference$pathway,
                               recurrence$pathway), , drop = FALSE]
print(recurrence, row.names = FALSE)
expected <- c(3L, 3L, 3L, 3L, 3L, 3L, 3L, 1L, 1L)
if (!identical(as.integer(recurrence$stages_with_BH_q_below_0.05),
               expected)) {
  warning("Source data do not reproduce 7 pathways at all stages +",
          " 2 pathways at one stage. Check the audit CSV.")
}

# -------------------- Figure -----------------------------
all_data$pathway <- factor(
  all_data$pathway, levels = rev(reference$pathway)
)
all_data$age <- factor(all_data$age,
                       levels = c(names(sheets), "GSEA"))

sig <- all_data[all_data$significant, , drop = FALSE]
not_sig <- all_data[!all_data$significant, , drop = FALSE]
gsea <- data.frame(
  pathway = factor(reference$pathway,
                   levels = levels(all_data$pathway)),
  age = factor(rep("GSEA", nrow(reference)),
               levels = levels(all_data$age)),
  text = ifelse(reference$tested_in_GSEA, "Yes", "–"),
  tested = reference$tested_in_GSEA,
  stringsAsFactors = FALSE
)

p <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 3.5,
                      linewidth = 0.45, colour = "#C6CDD3") +
  # Grey circles: pathway observed, but enrichment NOT significant
  # after Benjamini correction.
  ggplot2::geom_point(
    data = not_sig,
    ggplot2::aes(x = age, y = pathway, size = count),
    shape = 21, fill = "#ECEFF1", colour = "#AFBAC2",
    stroke = 0.45
  ) +
  # Coloured circles: Benjamini-adjusted enrichment P < 0.05.
  ggplot2::geom_point(
    data = sig,
    ggplot2::aes(x = age, y = pathway,
                 size = count, fill = minus_log10_q),
    shape = 21, colour = "#24495B", stroke = 0.35
  ) +
  ggplot2::geom_text(
    data = gsea,
    ggplot2::aes(x = age, y = pathway, label = text),
    colour = ifelse(gsea$tested, "#205F83", "#9AA5AD"),
    fontface = ifelse(gsea$tested, "bold", "plain"),
    size = 4.1
  ) +
  ggplot2::scale_x_discrete(drop = FALSE,
                            expand = ggplot2::expansion(add = 0.55)) +
  ggplot2::scale_y_discrete(drop = FALSE,
                            expand = ggplot2::expansion(add = 0.7)) +
  ggplot2::scale_size_continuous(
    name = "Downregulated genes\nin pathway",
    range = c(3.2, 8.5),
    breaks = c(4, 8, 12, 16, 20)
  ) +
  ggplot2::scale_fill_gradient(
    name = expression(-log[10]("Benjamini-adjusted P")),
    low = "#B9D9E8", high = "#155475"
  ) +
  ggplot2::labs(
    x = NULL, y = NULL,
    caption = paste0(
      "Coloured dots: Benjamini-adjusted P < ", Q_CUTOFF,
      "; grey dots: adjusted P >= ", Q_CUTOFF,
      ". Dot size indicates KEGG gene count. GSEA: pathway tested downstream."
    )
  ) +
  ggplot2::guides(
    size = ggplot2::guide_legend(
      order = 1, override.aes = list(fill = "#84BAD2")
    ),
    fill = ggplot2::guide_colourbar(order = 2)
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      linewidth = 0.3, colour = "#E8EDF0"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 12, face = "bold", colour = "#263B4B"
    ),
    axis.text.y = ggplot2::element_text(
      size = 11.5, colour = "#263B4B"
    ),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 10.5),
    legend.text = ggplot2::element_text(size = 10),
    plot.caption = ggplot2::element_text(
      size = 9, colour = "#586A78", hjust = 0,
      margin = ggplot2::margin(t = 13)
    ),
    plot.margin = ggplot2::margin(10, 18, 12, 10)
  )

print(p)

# PNG: high-resolution raster image with reliable text rendering.
ggplot2::ggsave(
  filename = paste0(OUTPUT_PREFIX, ".png"),
  plot = p,
  device = ragg::agg_png,
  width = 10.5, height = 5.7, units = "in",
  dpi = 600, bg = "white"
)

# SVG: scalable vector image, suitable for editing in Inkscape/Illustrator.
ggplot2::ggsave(
  filename = paste0(OUTPUT_PREFIX, ".svg"),
  plot = p,
  device = svglite::svglite,
  width = 10.5, height = 5.7, units = "in",
  bg = "white"
)

message("Saved PNG, SVG and values CSV to: ", normalizePath(getwd()))
