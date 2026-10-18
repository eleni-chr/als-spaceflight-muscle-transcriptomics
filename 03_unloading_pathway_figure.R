# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

###############################################################################
# FIGURE 4
# Metabolic pathway enrichment across mechanical-unloading models
#
# Single integrated dot plot:
#   - Six ALS-reference metabolic pathways
#   - Three broader Hallmark bioenergetic pathways
#
# Columns = unloading experiments
# Dot fill = NES
# Dot size = -log10(FDR), capped at 10
# * = FDR < 0.05
###############################################################################

library(tidyverse)
library(scales)
library(svglite)
library(grid)


# =============================================================================
# 1. INPUT FILES
# =============================================================================

# Inputs are written by 02_spaceflight_unloading_reference_analysis.R.

# OSD-21 hindlimb unloading
f_osd21 <- file.path("NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL", "GSEA",
  "OSD21_HLU_vs_CTRL_primary_GSEA_all.csv")

# OSD-540 unilateral lower-limb suspension
f_osd540 <- file.path("NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL", "GSEA", "OSD540_ULLS_vs_Baseline_primary_GSEA_all.csv")

# GSE271607 dry immersion
f_gse271607_vl <- file.path("NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL", "GSEA",
  "GSE271607_VL_After_vs_before_primary_GSEA_all.csv")

f_gse271607_sol <- file.path("NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL", "GSEA",
  "GSE271607_Soleus_After_vs_before_primary_GSEA_all.csv")


# =============================================================================
# 2. READ GSEA RESULTS
# =============================================================================

read_gsea <- function(file, comparison_name) {
  
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  x <- read_csv(file, show_col_types = FALSE)
  
  required <- c("pathway", "NES", "padj")
  
  if (!all(required %in% colnames(x))) {
    stop(
      "\nThe following file does not contain pathway, NES and padj columns:\n",
      file,
      "\n\nColumns present:\n",
      paste(colnames(x), collapse = ", ")
    )
  }
  
  x %>%
    dplyr::select(pathway, NES, padj) %>%
    mutate(
      comparison = comparison_name,
      pathway_upper = toupper(pathway)
    )
}


gsea_osd21 <- read_gsea(
  f_osd21,
  "OSD-21\nHLU"
)

gsea_osd540 <- read_gsea(
  f_osd540,
  "OSD-540\nULLS"
)

gsea_vl <- read_gsea(
  f_gse271607_vl,
  "GSE271607\nVastus lateralis"
)

gsea_sol <- read_gsea(
  f_gse271607_sol,
  "GSE271607\nSoleus"
)


gsea_unloading <- bind_rows(
  gsea_osd21,
  gsea_osd540,
  gsea_vl,
  gsea_sol
)


# =============================================================================
# 3. DEFINE PATHWAYS
# =============================================================================

# Six formally tested ALS-reference pathways
als_reference_sets <- c(
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "KEGG_PURINE_METABOLISM",
  "KEGG_STARCH_AND_SUCROSE_METABOLISM"
)

als_reference_labels <- c(
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS" =
    "Glycolysis / gluconeogenesis",
  
  "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM" =
    "Fructose and mannose metabolism",
  
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY" =
    "Pentose-phosphate pathway",
  
  "KEGG_INSULIN_SIGNALING_PATHWAY" =
    "Insulin signalling",
  
  "KEGG_PURINE_METABOLISM" =
    "Purine metabolism",
  
  "KEGG_STARCH_AND_SUCROSE_METABOLISM" =
    "Starch and sucrose metabolism"
)


# Three broader Hallmark bioenergetic pathways
hallmark_sets <- c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM"
)

hallmark_labels <- c(
  "HALLMARK_GLYCOLYSIS" =
    "Glycolysis",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" =
    "Oxidative phosphorylation",
  
  "HALLMARK_FATTY_ACID_METABOLISM" =
    "Fatty-acid metabolism"
)


comparison_levels <- c(
  "OSD-21\nHLU",
  "OSD-540\nULLS",
  "GSE271607\nVastus lateralis",
  "GSE271607\nSoleus"
)


# =============================================================================
# 4. BUILD INTEGRATED FIGURE DATA
# =============================================================================

als_data <- gsea_unloading %>%
  filter(pathway_upper %in% als_reference_sets) %>%
  mutate(
    pathway_label = recode(
      pathway_upper,
      !!!als_reference_labels
    ),
    pathway_group = "ALS-reference\nmetabolic pathways"
  )


hallmark_data <- gsea_unloading %>%
  filter(pathway_upper %in% hallmark_sets) %>%
  mutate(
    pathway_label = recode(
      pathway_upper,
      !!!hallmark_labels
    ),
    pathway_group = "Broader\nbioenergetic\npathways"
  )


figure4_data <- bind_rows(
  als_data,
  hallmark_data
) %>%
  mutate(
    
    comparison = factor(
      comparison,
      levels = comparison_levels
    ),
    
    pathway_group = factor(
      pathway_group,
      levels = c(
        "ALS-reference\nmetabolic pathways",
        "Broader\nbioenergetic\npathways"
      )
    ),
    
    pathway_label = factor(
      pathway_label,
      levels = rev(
        c(
          unname(als_reference_labels),
          unname(hallmark_labels)
        )
      )
    ),
    
    # Protect against padj = 0
    padj_for_plot = pmax(padj, 1e-300),
    
    logFDR = -log10(padj_for_plot),
    
    # Cap very small FDR values so extremely significant results
    # do not create disproportionately large dots
    logFDR_plot = pmin(logFDR, 10),
    
    significant = padj < 0.05
  )


# =============================================================================
# 5. CHECK EXPECTED RESULTS
# =============================================================================

cat(
  "\nFigure 4 rows:",
  nrow(figure4_data),
  "\n"
)

cat(
  "Expected rows: 36 (9 pathways x 4 comparisons)\n\n"
)


expected_pathways <- c(
  als_reference_sets,
  hallmark_sets
)

expected_combinations <- expand_grid(
  pathway_upper = expected_pathways,
  comparison = comparison_levels
)

observed_combinations <- figure4_data %>%
  transmute(
    pathway_upper,
    comparison = as.character(comparison)
  )

missing_combinations <- anti_join(
  expected_combinations,
  observed_combinations,
  by = c("pathway_upper", "comparison")
)

if (nrow(missing_combinations) > 0) {
  
  message(
    "\nMissing pathway/comparison combinations:"
  )
  
  print(missing_combinations)
  
} else {
  
  message(
    "\nAll 36 expected pathway/comparison combinations were found."
  )
}


# =============================================================================
# 6. COMMON NES COLOUR SCALE
# =============================================================================

nes_limit <- max(
  abs(figure4_data$NES),
  na.rm = TRUE
)

# Round upwards to nearest 0.5
nes_limit <- ceiling(nes_limit * 2) / 2

cat(
  "\nNES colour scale:",
  -nes_limit,
  "to",
  nes_limit,
  "\n"
)


# =============================================================================
# 7. CREATE FIGURE 4
# =============================================================================

figure4 <- ggplot(
  figure4_data,
  aes(
    x = comparison,
    y = pathway_label
  )
) +
  
  # Light horizontal guides
  geom_hline(
    yintercept = seq_along(
      levels(figure4_data$pathway_label)
    ),
    colour = "grey92",
    linewidth = 0.35
  ) +
  
  # Light vertical guides
  geom_vline(
    xintercept = seq_along(
      levels(figure4_data$comparison)
    ),
    colour = "grey92",
    linewidth = 0.35
  ) +
  
  # Main pathway dots
  geom_point(
    aes(
      fill = NES,
      size = logFDR_plot
    ),
    shape = 21,
    colour = "grey25",
    stroke = 0.35
  ) +
  
  # Significance marker
  geom_text(
    data = figure4_data %>%
      filter(significant),
    
    aes(label = "*"),
    
    size = 5.5,
    fontface = "bold",
    colour = "black",
    nudge_y = 0.20,
    show.legend = FALSE
  ) +
  
  # Separate the pathway classes while retaining a single integrated figure
  facet_grid(
    rows = vars(pathway_group),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  # NES colour scale
  scale_fill_gradient2(
    low = "#E98B7D",
    mid = "#FFFDF7",
    high = "#70B7D6",
    midpoint = 0,
    limits = c(
      -nes_limit,
      nes_limit
    ),
    oob = scales::squish,
    name = "NES"
  ) +
  
  # FDR size scale
  scale_size_continuous(
    name = expression(-log[10]("FDR")),
    range = c(3.5, 8),
    limits = c(0, 10),
    breaks = c(
      2,
      4,
      6,
      8,
      10
    ),
    labels = c(
      "2",
      "4",
      "6",
      "8",
      "\u226510"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    
    # X-axis labels
    axis.text.x = element_text(
      size = 13,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 10)
    ),
    
    # Pathway labels
    axis.text.y = element_text(
      size = 13,
      colour = "black"
    ),
    
    axis.ticks = element_blank(),
    
    # Put pathway-class headings on the left
    strip.placement = "outside",
    
    strip.background = element_blank(),
    
    # Vertical pathway-class headings
    strip.text.y.left = element_text(
      size = 13,
      face = "bold",
      colour = "black",
      angle = 90,
      lineheight = 0.9,
      margin = margin(
        r = 14,
        l = 14
      )
    ),
    
    # Larger visual separation between the two pathway groups
    panel.spacing.y = unit(
      1.5,
      "lines"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    # Extra left margin prevents vertical group headings being clipped
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 10,
      l = 55
    )
  )


figure4


# =============================================================================
# 8. SAVE FIGURE
# =============================================================================

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)


ggsave(
  filename = file.path(
    "figures",
    "Figure4_mechanical_unloading_pathway_enrichment.svg"
  ),
  plot = figure4,
  width = 10.5,
  height = 7,
  device = svglite::svglite
)


ggsave(
  filename = file.path(
    "figures",
    "Figure4_mechanical_unloading_pathway_enrichment.pdf"
  ),
  plot = figure4,
  width = 10.5,
  height = 7
)


ggsave(
  filename = file.path(
    "figures",
    "Figure4_mechanical_unloading_pathway_enrichment.png"
  ),
  plot = figure4,
  width = 10.5,
  height = 7,
  dpi = 600
)


# =============================================================================
# 9. SAVE EXACT SOURCE DATA USED IN FIGURE 4
# =============================================================================

figure4_source_data <- figure4_data %>%
  dplyr::select(
    pathway_group,
    comparison,
    pathway,
    pathway_label,
    NES,
    padj,
    logFDR
  )


write_csv(
  figure4_source_data,
  file.path(
    "figures",
    "Figure4_source_data.csv"
  )
)


message("\nFigure 4 complete.")
