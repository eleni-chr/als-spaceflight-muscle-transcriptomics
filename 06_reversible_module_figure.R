# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

###############################################################################
# FIGURE 7
# A recurrent and reversible gene-level programme links ALS muscle
# with mechanical unloading
#
# Run AFTER:
# 05_reversible_module_and_exercise_analysis.R
#
# Required objects/functions:
#   sod1
#   load_effects
#   full_module
#   perm
#   loso
#   build_sod1_consensus()
#   FIG_DIR
#
# Panel A:
#   i    SOD1G93A disease stages
#   ii   Mechanical unloading
#   iii  Restoration of loading
#   iv   Module features
#
# Panel B:
#   Empirical permutation test
#
# Panel C:
#   SOD1G93A leave-one-stage-out robustness
###############################################################################


# =============================================================================
# 1. PACKAGES
# =============================================================================

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "patchwork",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(scales)
})


# =============================================================================
# 2. CHECK REQUIRED OBJECTS / FUNCTIONS
# =============================================================================

required_objects <- c(
  "sod1",
  "load_effects",
  "full_module",
  "perm",
  "loso"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    inherits = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing: ",
    paste(missing_objects, collapse = ", "),
    "\nRun the final reversible ALS-unloading analysis script first."
  )
}

if (!exists("build_sod1_consensus")) {
  stop(
    "Function build_sod1_consensus() is missing.\n",
    "Run the final reversible ALS-unloading analysis script first."
  )
}

if (!exists("FIG_DIR")) {
  FIG_DIR <- "."
}

if (!dir.exists(FIG_DIR)) {
  dir.create(
    FIG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


MIN_LFC_FIGURE <- if (exists("MIN_ABS_LFC")) {
  MIN_ABS_LFC
} else {
  0.10
}

MIN_SOD1_SUPPORT_FIG <- if (exists("MIN_SOD1_SUPPORT_FULL")) {
  MIN_SOD1_SUPPORT_FULL
} else {
  2L
}

# Visual gap between downregulated and upregulated gene blocks
HEATMAP_GAP <- 4


# =============================================================================
# 3. CHECK FINAL PERMUTATION RUN
# =============================================================================

if (length(perm$null) != 20000) {
  stop(
    "\nThe manuscript reports the final 20,000-permutation analysis.\n",
    "Current perm$null contains ",
    length(perm$null),
    " permutations.\n\n",
    "Set N_PERM <- 20000 in the main analysis script, rerun it, ",
    "then rerun this figure code."
  )
}


# =============================================================================
# 4. PUBLICATION THEME — ENLARGED TEXT
# =============================================================================

theme_figure6 <- ggplot2::theme_classic(
  base_size = 14
) +
  ggplot2::theme(
    
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 17,
      hjust = 0
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 12.5,
      hjust = 0,
      margin = ggplot2::margin(b = 10)
    ),
    
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text = ggplot2::element_text(
      size = 11.5
    ),
    
    strip.background = ggplot2::element_blank(),
    
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 12.5
    ),
    
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = ggplot2::element_text(
      size = 11
    ),
    
    plot.margin = ggplot2::margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )


# =============================================================================
# 5. PREPARE MODULE
# =============================================================================

module_fig <- full_module %>%
  dplyr::filter(
    directional_module
  ) %>%
  dplyr::mutate(
    
    all_three_unloading =
      n_unloading_available >= 3L &
      n_unloading_concordant >= 3L,
    
    both_recoveries =
      two_recovery_models,
    
    statistical_support =
      statistically_supported,
    
    high_confidence_gene =
      high_confidence
  )


# =============================================================================
# 6. MODULE COUNTS
# =============================================================================

n_module <- nrow(module_fig)

n_down <- sum(
  module_fig$als_direction == -1L,
  na.rm = TRUE
)

n_up <- sum(
  module_fig$als_direction == 1L,
  na.rm = TRUE
)

n_all3_unloading <- sum(
  module_fig$all_three_unloading,
  na.rm = TRUE
)

n_both_recovery <- sum(
  module_fig$both_recoveries,
  na.rm = TRUE
)

n_stat_support <- sum(
  module_fig$statistical_support,
  na.rm = TRUE
)

n_high_conf <- sum(
  module_fig$high_confidence_gene,
  na.rm = TRUE
)


# =============================================================================
# 7. TOTAL RECURRENT SOD1G93A GENES VS EVALUABLE GENES
# =============================================================================
#
# 529 = all recurrent SOD1G93A genes
# 526 = recurrent SOD1G93A genes evaluable in the loading-data universe
# 174 = final reversible ALS-unloading module
# =============================================================================

sod1_recurrent_all <- build_sod1_consensus(
  sod1,
  c(
    "90",
    "120",
    "150"
  )
) %>%
  dplyr::filter(
    sod1_direction != 0L,
    n_sod1_concordant >= MIN_SOD1_SUPPORT_FIG
  )

n_recurrent_sod1_total <-
  nrow(sod1_recurrent_all)

n_recurrent_sod1_evaluable <-
  as.integer(
    perm$summary$n_recurrently_significant_sod1_genes_tested[1]
  )

n_recurrent_not_evaluable <-
  n_recurrent_sod1_total -
  n_recurrent_sod1_evaluable


message("\nFigure 7 counts:")
message("Total recurrent SOD1G93A genes: ", n_recurrent_sod1_total)
message("Evaluable recurrent genes: ", n_recurrent_sod1_evaluable)
message("Not evaluable in loading universe: ", n_recurrent_not_evaluable)
message("Directional reversible module: ", n_module)
message("SOD1G93A downregulated: ", n_down)
message("SOD1G93A upregulated: ", n_up)
message("Concordant in all 3 unloading studies: ", n_all3_unloading)
message("Reversed in both recovery studies: ", n_both_recovery)
message("Statistically supported: ", n_stat_support)
message("High confidence: ", n_high_conf)


# =============================================================================
# 8. CHECK COUNTS
# =============================================================================

expected_counts <- c(
  recurrent_sod1_total = 529,
  recurrent_sod1_evaluable = 526,
  module = 174,
  down = 128,
  up = 46,
  all3_unloading = 59,
  both_recovery = 85,
  statistical_support = 102,
  high_confidence = 55
)

observed_counts <- c(
  recurrent_sod1_total = n_recurrent_sod1_total,
  recurrent_sod1_evaluable = n_recurrent_sod1_evaluable,
  module = n_module,
  down = n_down,
  up = n_up,
  all3_unloading = n_all3_unloading,
  both_recovery = n_both_recovery,
  statistical_support = n_stat_support,
  high_confidence = n_high_conf
)

if (!all(observed_counts == expected_counts)) {
  
  warning(
    "\nSome Figure 7 counts differ from the expected final values:\n",
    paste(
      names(observed_counts),
      ": observed = ",
      observed_counts,
      ", expected = ",
      expected_counts,
      collapse = "\n"
    )
  )
}


# =============================================================================
# 9. ORDER GENES FOR PANEL A
# =============================================================================

module_fig <- module_fig %>%
  
  dplyr::mutate(
    
    direction_order = dplyr::case_when(
      als_direction == -1L ~ 1L,
      als_direction == 1L ~ 2L,
      TRUE ~ 3L
    )
  ) %>%
  
  dplyr::arrange(
    
    direction_order,
    
    dplyr::desc(
      high_confidence_gene
    ),
    
    dplyr::desc(
      both_recoveries
    ),
    
    dplyr::desc(
      statistical_support
    ),
    
    dplyr::desc(
      all_three_unloading
    ),
    
    dplyr::desc(
      n_unloading_concordant
    ),
    
    dplyr::desc(
      n_recovery_reverse
    ),
    
    dplyr::desc(
      n_sod1_concordant
    ),
    
    symbol
  ) %>%
  
  dplyr::mutate(
    
    row_index =
      dplyr::row_number(),
    
    plot_y_base =
      n_module + 1 - row_index,
    
    # Shift the downregulated block upward to create a visual gap
    plot_y =
      dplyr::if_else(
        als_direction == -1L,
        plot_y_base + HEATMAP_GAP,
        plot_y_base
      )
  )


# =============================================================================
# 10. DIRECTION-BLOCK LABEL POSITIONS
# =============================================================================

up_block_midpoint <-
  (1 + n_up) / 2

down_block_midpoint <-
  (
    (n_up + 1) +
      n_module
  ) / 2 +
  HEATMAP_GAP


row_block_breaks <- c(
  down_block_midpoint,
  up_block_midpoint
)

# Single-line labels because they are rotated vertically
row_block_labels <- c(
  
  paste0(
    "Downregulated\nin SOD1G93A\n(n=",
    n_down,
    ")"
  ),
  
  paste0(
    "Upregulated\nin SOD1G93A\n(n=",
    n_up,
    ")"
  )
)


panel_A_y_limits <- c(
  0.5,
  n_module + HEATMAP_GAP + 0.5
)


# =============================================================================
# 11. CONDITION DEFINITIONS
# =============================================================================

condition_key <- tibble::tribble(
  
  ~condition_code,
  ~condition_group,
  ~condition_label,
  ~condition_order,
  
  "SOD90",
  "SOD1G93A",
  "90 d",
  1,
  
  "SOD120",
  "SOD1G93A",
  "120 d",
  2,
  
  "SOD150",
  "SOD1G93A",
  "150 d",
  3,
  
  "OSD21_UNLOAD",
  "Mechanical unloading",
  "OSD-21\nHLU",
  4,
  
  "OSD540_UNLOAD",
  "Mechanical unloading",
  "OSD-540\nULLS",
  5,
  
  "GSE271607_UNLOAD",
  "Mechanical unloading",
  "GSE271607\nDry immersion",
  6,
  
  "OSD21_RECOVERY",
  "Restoration of loading",
  "OSD-21\nReloading",
  7,
  
  "OSD540_RECOVERY",
  "Restoration of loading",
  "OSD-540\nActive recovery",
  8
) %>%
  
  dplyr::mutate(
    
    condition_group = factor(
      condition_group,
      levels = c(
        "SOD1G93A",
        "Mechanical unloading",
        "Restoration of loading"
      )
    )
  )


# =============================================================================
# 12. FACET LABELS: i-iii
# =============================================================================

facet_labels <- c(
  
  SOD1G93A =
    "'i'~SOD1^G93A",
  
  `Mechanical unloading` =
    "'ii'~'Mechanical unloading'",
  
  `Restoration of loading` =
    "'iii'~'Restoration of loading'"
)


# =============================================================================
# 13. SOD1 EFFECTS
# =============================================================================
#
# SOD1 input is significant-gene-only.
# Therefore an absent gene at a stage means:
# "Not significantly detected"
# =============================================================================

sod1_heat <- sod1 %>%
  
  dplyr::mutate(
    stage_chr = as.character(stage)
  ) %>%
  
  dplyr::filter(
    
    symbol %in%
      module_fig$symbol,
    
    stage_chr %in%
      c(
        "90",
        "120",
        "150"
      )
  ) %>%
  
  dplyr::transmute(
    
    symbol,
    
    condition_code =
      dplyr::case_when(
        
        stage_chr == "90" ~
          "SOD90",
        
        stage_chr == "120" ~
          "SOD120",
        
        stage_chr == "150" ~
          "SOD150",
        
        TRUE ~
          NA_character_
      ),
    
    log2FC
  ) %>%
  
  dplyr::filter(
    !is.na(
      condition_code
    )
  )


# =============================================================================
# 14. UNLOADING / RECOVERY EFFECTS
# =============================================================================
#
# Missing gene-level effect:
# "Not available"
# =============================================================================

loading_heat <- load_effects %>%
  
  dplyr::filter(
    symbol %in%
      module_fig$symbol
  ) %>%
  
  dplyr::mutate(
    
    condition_code =
      dplyr::case_when(
        
        study_id == "OSD21" &
          role == "unloading" ~
          "OSD21_UNLOAD",
        
        study_id == "OSD540" &
          role == "unloading" ~
          "OSD540_UNLOAD",
        
        study_id == "GSE271607" &
          role == "unloading" ~
          "GSE271607_UNLOAD",
        
        study_id == "OSD21" &
          role == "recovery" ~
          "OSD21_RECOVERY",
        
        study_id == "OSD540" &
          role == "recovery" ~
          "OSD540_RECOVERY",
        
        TRUE ~
          NA_character_
      )
  ) %>%
  
  dplyr::filter(
    !is.na(
      condition_code
    )
  ) %>%
  
  dplyr::transmute(
    symbol,
    condition_code,
    log2FC
  ) %>%
  
  dplyr::distinct(
    symbol,
    condition_code,
    .keep_all = TRUE
  )


# =============================================================================
# 15. COMBINE EFFECT TABLES
# =============================================================================

heat_effects <- dplyr::bind_rows(
  sod1_heat,
  loading_heat
)


# =============================================================================
# 16. COMPLETE HEATMAP MATRIX
# =============================================================================
#
# Red:
#   concordant with recurrent SOD1G93A consensus
#
# Beige:
#   measured effect, but |log2FC| < threshold
#
# Blue:
#   opposite to recurrent SOD1G93A consensus
#
# Light grey:
#   not significantly detected in significant-only SOD1 stage dataset
#
# Dark grey:
#   gene/effect not available in unloading/recovery comparison
# =============================================================================

heatmap_data <- tidyr::crossing(
  
  symbol =
    module_fig$symbol,
  
  condition_code =
    condition_key$condition_code
  
) %>%
  
  dplyr::left_join(
    
    heat_effects,
    
    by = c(
      "symbol",
      "condition_code"
    )
  ) %>%
  
  dplyr::left_join(
    
    module_fig %>%
      
      dplyr::select(
        symbol,
        als_direction,
        plot_y
      ),
    
    by =
      "symbol"
  ) %>%
  
  dplyr::left_join(
    
    condition_key,
    
    by =
      "condition_code"
  ) %>%
  
  dplyr::mutate(
    
    state =
      dplyr::case_when(
        
        is.na(log2FC) &
          condition_group == "SOD1G93A" ~
          "Not significantly detected",
        
        is.na(log2FC) &
          condition_group != "SOD1G93A" ~
          "Not available",
        
        abs(log2FC) <
          MIN_LFC_FIGURE ~
          "Below directional threshold",
        
        sign(log2FC) ==
          als_direction ~
          "Concordant with consensus",
        
        sign(log2FC) ==
          -als_direction ~
          "Opposite to consensus",
        
        TRUE ~
          "Below directional threshold"
      ),
    
    state =
      factor(
        state,
        levels = c(
          "Concordant with consensus",
          "Below directional threshold",
          "Opposite to consensus",
          "Not significantly detected",
          "Not available"
        )
      ),
    
    condition_code =
      factor(
        condition_code,
        levels =
          condition_key$condition_code
      )
  )


# =============================================================================
# 17. CONDITION LABELS
# =============================================================================

condition_labels <- stats::setNames(
  condition_key$condition_label,
  condition_key$condition_code
)


# =============================================================================
# 18. PANEL A — i-iii DIRECTIONAL HEATMAPS
# =============================================================================

panel_A_subtitle <- paste0(
  
  n_recurrent_sod1_total,
  " recurrent SOD1G93A genes identified; ",
  
  n_recurrent_sod1_evaluable,
  " were evaluable across the loading datasets; ",
  
  n_module,
  " met the module criteria"
)


pA_heat <- ggplot2::ggplot(
  
  heatmap_data,
  
  ggplot2::aes(
    x = condition_code,
    y = plot_y,
    fill = state
  )
  
) +
  
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.08
  ) +
  
  ggplot2::facet_grid(
    
    . ~ condition_group,
    
    scales = "free_x",
    
    space = "free_x",
    
    labeller =
      ggplot2::labeller(
        condition_group =
          ggplot2::as_labeller(
            facet_labels,
            label_parsed
          )
      )
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = c(
      
      "Concordant with consensus" =
        "#B84A48",
      
      "Below directional threshold" =
        "#E8DED1",
      
      "Opposite to consensus" =
        "#416FA6",
      
      "Not significantly detected" =
        "#D9D9D9",
      
      "Not available" =
        "#9EA4AA"
    ),
    
    drop = FALSE
  ) +
  
  ggplot2::scale_x_discrete(
    labels = condition_labels,
    drop = TRUE
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits =
      panel_A_y_limits,
    
    breaks =
      row_block_breaks,
    
    labels =
      row_block_labels,
    
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::labs(
    
    title =
      "A   Recurrent reversible ALS–unloading programme",
    
    subtitle =
      panel_A_subtitle,
    
    x =
      NULL,
    
    y =
      NULL,
    
    fill =
      "Direction relative to recurrent\nSOD1G93A consensus"
  ) +
  
  theme_figure6 +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_text(
        size = 11.5,
        lineheight = 0.98,
        margin = ggplot2::margin(
          t = 5
        )
      ),
    
    # Vertical labels alongside the downregulated/upregulated blocks
    axis.text.y =
      ggplot2::element_text(
        size = 10.5,
        colour = "black",
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        margin = ggplot2::margin(
          r = 8
        )
      ),
    
    axis.ticks.y =
      ggplot2::element_blank(),
    
    axis.line.y =
      ggplot2::element_blank(),
    
    panel.spacing.x =
      grid::unit(
        8,
        "pt"
      ),
    
    legend.position =
      "bottom",
    
    legend.key.width =
      grid::unit(
        0.56,
        "cm"
      ),
    
    legend.key.height =
      grid::unit(
        0.44,
        "cm"
      )
  ) +
  
  ggplot2::guides(
    
    fill =
      ggplot2::guide_legend(
        title.position = "left",
        title.hjust = 0,
        nrow = 2,
        byrow = TRUE
      )
  )


# =============================================================================
# 19. PANEL A — iv MODULE FEATURE MATRIX
# =============================================================================

feature_labels <- c(
  
  all_three_unloading =
    paste0(
      "All 3\nunloading\n(n=",
      n_all3_unloading,
      ")"
    ),
  
  both_recoveries =
    paste0(
      "Both\nrecoveries\n(n=",
      n_both_recovery,
      ")"
    ),
  
  statistical_support =
    paste0(
      "Statistical\nsupport\n(n=",
      n_stat_support,
      ")"
    ),
  
  high_confidence_gene =
    paste0(
      "High\nconfidence\n(n=",
      n_high_conf,
      ")"
    )
)


feature_data <- module_fig %>%
  
  dplyr::select(
    symbol,
    plot_y,
    all_three_unloading,
    both_recoveries,
    statistical_support,
    high_confidence_gene
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      all_three_unloading,
      both_recoveries,
      statistical_support,
      high_confidence_gene
    ),
    
    names_to =
      "feature",
    
    values_to =
      "present"
  ) %>%
  
  dplyr::mutate(
    
    feature =
      factor(
        feature,
        levels = c(
          "all_three_unloading",
          "both_recoveries",
          "statistical_support",
          "high_confidence_gene"
        )
      ),
    
    present =
      factor(
        present,
        levels = c(
          FALSE,
          TRUE
        ),
        labels = c(
          "Not present",
          "Present"
        )
      ),
    
    feature_group =
      factor(
        "iv   Module features",
        levels = "iv   Module features"
      )
  )


pA_features <- ggplot2::ggplot(
  
  feature_data,
  
  ggplot2::aes(
    x = feature,
    y = plot_y,
    fill = present
  )
  
) +
  
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.08
  ) +
  
  ggplot2::facet_grid(
    . ~ feature_group
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = c(
      "Not present" = "#E3E5E7",
      "Present" = "#444444"
    ),
    
    drop = FALSE
  ) +
  
  ggplot2::scale_x_discrete(
    labels = feature_labels
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits =
      panel_A_y_limits,
    
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::labs(
    
    title = NULL,
    
    subtitle = NULL,
    
    x = NULL,
    
    y = NULL,
    
    fill = "Module feature"
  ) +
  
  theme_figure6 +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_text(
        size = 10.5,
        lineheight = 0.92,
        margin = ggplot2::margin(
          t = 5
        )
      ),
    
    axis.text.y =
      ggplot2::element_blank(),
    
    axis.ticks =
      ggplot2::element_blank(),
    
    axis.line.x =
      ggplot2::element_blank(),
    
    axis.line.y =
      ggplot2::element_blank(),
    
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    legend.position =
      "bottom",
    
    legend.title =
      ggplot2::element_text(
        face = "bold",
        size = 10.5
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 10
      ),
    
    legend.key.width =
      grid::unit(
        0.48,
        "cm"
      ),
    
    legend.key.height =
      grid::unit(
        0.40,
        "cm"
      ),
    
    legend.key =
      ggplot2::element_rect(
        colour = "#B8B8B8",
        linewidth = 0.3
      )
  ) +
  
  ggplot2::guides(
    
    fill =
      ggplot2::guide_legend(
        title.position = "left",
        title.hjust = 0,
        nrow = 1,
        byrow = TRUE
      )
  )


# =============================================================================
# 20. COMBINE PANEL A
# =============================================================================

panel_A <- pA_heat + pA_features +
  
  patchwork::plot_layout(
    widths = c(
      5.2,
      2.0
    ),
    guides = "keep"
  )


# =============================================================================
# 21. PANEL B — PERMUTATION TEST
# =============================================================================

perm_summary <-
  perm$summary[1, ]


observed_module_size <-
  as.numeric(
    perm_summary$observed_module_size
  )


null_mean <-
  as.numeric(
    perm_summary$null_mean
  )


null_95 <-
  as.numeric(
    perm_summary$null_95th_percentile
  )


perm_p <-
  as.numeric(
    perm_summary$empirical_p_one_sided
  )


null_plot_data <- tibble::tibble(
  module_size =
    perm$null
)


if (perm_p < 0.001) {
  
  p_text <- formatC(
    perm_p,
    format = "e",
    digits = 1
  )
  
} else {
  
  p_text <- formatC(
    perm_p,
    format = "f",
    digits = 3
  )
}


hist_breaks <- seq(
  
  floor(
    min(
      null_plot_data$module_size
    )
  ) - 0.5,
  
  ceiling(
    max(
      null_plot_data$module_size
    )
  ) + 0.5,
  
  by = 1
)


hist_counts <- graphics::hist(
  
  null_plot_data$module_size,
  
  breaks =
    hist_breaks,
  
  plot =
    FALSE
)


hist_ymax <- max(
  hist_counts$counts
)


panel_B_xmin <- floor(
  min(
    null_plot_data$module_size
  ) / 10
) * 10


panel_B_xmax <- ceiling(
  (
    observed_module_size + 5
  ) / 10
) * 10


panel_B_breaks <- seq(
  panel_B_xmin,
  panel_B_xmax,
  by = 20
)


pB <- ggplot2::ggplot(
  
  null_plot_data,
  
  ggplot2::aes(
    x = module_size
  )
  
) +
  
  ggplot2::geom_histogram(
    
    binwidth = 1,
    
    boundary = 0.5,
    
    fill = "#BFBFBF",
    
    colour = "white",
    
    linewidth = 0.15
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = null_mean,
    
    linewidth = 0.8,
    
    linetype = "dotted",
    
    colour = "#444444"
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = null_95,
    
    linewidth = 0.8,
    
    linetype = "dashed",
    
    colour = "#555555"
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = observed_module_size,
    
    linewidth = 1.15,
    
    colour = "#B84A48"
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = null_mean,
    
    y = hist_ymax * 0.97,
    
    label =
      paste0(
        "Null mean = ",
        round(
          null_mean,
          1
        )
      ),
    
    angle = 90,
    
    hjust = 1,
    
    vjust = -0.5,
    
    size = 4.0,
    
    colour = "#444444"
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = null_95,
    
    y = hist_ymax * 0.97,
    
    label =
      paste0(
        "95th percentile = ",
        round(
          null_95,
          0
        )
      ),
    
    angle = 90,
    
    hjust = 1,
    
    vjust = -0.5,
    
    size = 4.0,
    
    colour = "#444444"
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = observed_module_size,
    
    y = hist_ymax * 0.97,
    
    label =
      paste0(
        "Observed = ",
        observed_module_size
      ),
    
    angle = 90,
    
    hjust = 1,
    
    vjust = -0.55,
    
    size = 4.2,
    
    fontface = "bold",
    
    colour = "#B84A48"
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      panel_B_xmin,
      panel_B_xmax
    ),
    
    breaks =
      panel_B_breaks,
    
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::scale_y_continuous(
    
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.08
        )
      )
  ) +
  
  ggplot2::labs(
    
    title =
      "B   Cross-dataset convergence exceeds chance",
    
    subtitle =
      paste0(
        "20,000 gene-label permutations; empirical P = ",
        p_text
      ),
    
    x =
      "Reversible module size under permutation",
    
    y =
      "Number of permutations"
  ) +
  
  theme_figure6 +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 15
      )
  )


# =============================================================================
# 22. PANEL C — SOD1G93A LEAVE-ONE-STAGE-OUT ROBUSTNESS
# =============================================================================

loso_plot_data <- loso$summary %>%
  
  dplyr::mutate(
    
    stage =
      factor(
        as.character(
          stage
        ),
        levels = c(
          "90",
          "120",
          "150"
        )
      ),
    
    stage_label =
      factor(
        as.character(
          stage
        ),
        levels = c(
          "90",
          "120",
          "150"
        ),
        labels = c(
          "90 d\nPresymptomatic",
          "120 d\nEarly symptomatic",
          "150 d\nLate symptomatic"
        )
      ),
    
    recurrence_percent =
      expected_direction_recurrence_fraction *
      100,
    
    point_label =
      paste0(
        
        sprintf(
          "%.1f%%",
          recurrence_percent
        ),
        
        "\n",
        
        "(",
        
        n_significant_expected_direction,
        
        "/",
        
        n_independently_defined_module_genes,
        
        ")"
      )
  )


pC_title <- expression(
  bold(
    paste(
      "C   Module recurrence across the ",
      SOD1^G93A,
      " disease course"
    )
  )
)


pC <- ggplot2::ggplot(
  
  loso_plot_data,
  
  ggplot2::aes(
    x = stage_label,
    y = recurrence_percent,
    group = 1
  )
  
) +
  
  ggplot2::geom_line(
    linewidth = 1.0,
    colour = "#444444"
  ) +
  
  ggplot2::geom_point(
    
    size = 4.6,
    
    shape = 21,
    
    stroke = 1.0,
    
    fill = "#416FA6",
    
    colour = "white"
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        point_label
    ),
    
    vjust = 1.8,
    
    size = 4.1,
    
    fontface = "bold",
    
    lineheight = 0.92
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits = c(
      0,
      100
    ),
    
    breaks = seq(
      0,
      100,
      by = 25
    ),
    
    labels =
      function(x) {
        paste0(
          x,
          "%"
        )
      },
    
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.02
        )
      )
  ) +
  
  ggplot2::labs(
    
    title =
      pC_title,
    
    subtitle =
      "Each stage was tested using a module defined without that stage",
    
    x =
      NULL,
    
    y =
      "Module genes recurring\nin expected direction (%)"
  ) +
  
  theme_figure6 +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 15
      ),
    
    axis.text.x =
      ggplot2::element_text(
        size = 11,
        lineheight = 0.98
      )
  )


# =============================================================================
# 23. COMBINE PANELS B + C
# =============================================================================

bottom_row <- pB + pC +
  
  patchwork::plot_layout(
    widths = c(
      1.08,
      1
    )
  )


# =============================================================================
# 24. FINAL FIGURE
# =============================================================================

Figure_7 <- panel_A / bottom_row +
  
  patchwork::plot_layout(
    heights = c(
      1.55,
      1
    )
  )


print(
  Figure_7
)


# =============================================================================
# 25. SAVE
# =============================================================================

ggplot2::ggsave(
  
  filename = file.path(
    FIG_DIR,
    "Figure_7_recurrent_reversible_ALS_unloading_module.svg"
  ),
  
  plot = Figure_7,
  
  width = 15,
  
  height = 10.8,
  
  units = "in",
  
  bg = "white"
)


ggplot2::ggsave(
  
  filename = file.path(
    FIG_DIR,
    "Figure_7_recurrent_reversible_ALS_unloading_module.pdf"
  ),
  
  plot = Figure_7,
  
  width = 15,
  
  height = 10.8,
  
  units = "in",
  
  bg = "white"
)


ggplot2::ggsave(
  
  filename = file.path(
    FIG_DIR,
    "Figure_7_recurrent_reversible_ALS_unloading_module.png"
  ),
  
  plot = Figure_7,
  
  width = 15,
  
  height = 10.8,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white"
)


message(
  "\nFigure 7 saved to:\n",
  normalizePath(
    FIG_DIR,
    winslash = "/",
    mustWork = FALSE
  )
)
