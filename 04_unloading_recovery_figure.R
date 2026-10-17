# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# ============================================================
# FIGURE 5
# Reversal of unloading-associated metabolic pathway changes
#
# Panels:
#   A. OSD-21
#   B. OSD-540
#
# x-axis = pathway enrichment during unloading
# y-axis = pathway enrichment during restoration/recovery
#
# Opposite-sign quadrants indicate reversal.
# Dashed y = -x line indicates equal-magnitude reversal.
#
# Point fill = residual NES after restoration vs baseline/control
#
# * = FDR < 0.05 during BOTH unloading and restoration,
#     with opposite NES directions
#
# Short pathway labels are used directly on the figure to
# prevent overcrowding.
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)
library(grid)


# ============================================================
# 2. Paths
# ============================================================

base_dir <- file.path(getwd(), "NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL")

summary_dir <- file.path(base_dir, "summary")

figure_dir <- file.path(
  base_dir,
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ref_file <- file.path(
  summary_dir,
  "SOD1_reference_pathways_all_datasets.csv"
)

bio_file <- file.path(
  summary_dir,
  "Bioenergetic_exact_pathways_all_datasets.csv"
)


if (!file.exists(ref_file)) {
  stop("Cannot find: ", ref_file)
}

if (!file.exists(bio_file)) {
  stop("Cannot find: ", bio_file)
}


# ============================================================
# 3. Read data
# ============================================================

ref_raw <- read.csv(
  ref_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bio_raw <- read.csv(
  bio_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. Keep primary analysis if appropriate
# ============================================================

keep_primary <- function(x) {
  
  if ("variant" %in% names(x)) {
    
    if (any(
      x$variant == "primary",
      na.rm = TRUE
    )) {
      
      x <- x %>%
        filter(
          variant == "primary"
        )
    }
  }
  
  x
}


ref_raw <- keep_primary(ref_raw)
bio_raw <- keep_primary(bio_raw)


# ============================================================
# 5. Cleaning functions
# ============================================================

clean_name <- function(x) {
  
  x %>%
    toupper() %>%
    str_replace_all(
      "[^A-Z0-9]+",
      "_"
    ) %>%
    str_replace_all(
      "^_|_$",
      ""
    )
}


clean_dataset <- function(x) {
  
  x %>%
    toupper() %>%
    str_replace_all(
      "[^A-Z0-9]",
      ""
    )
}


# ============================================================
# 6. ALS-reference pathways
# ============================================================

ref_dat <- ref_raw %>%
  
  mutate(
    
    dataset_id =
      clean_dataset(dataset),
    
    pathway_id =
      clean_name(reference_pathway)
    
  ) %>%
  
  mutate(
    
    pathway = case_when(
      
      str_detect(
        pathway_id,
        "GLYCOLYSIS.*GLUCONEOGENESIS"
      ) ~
        "Glycolysis / gluconeogenesis",
      
      str_detect(
        pathway_id,
        "FRUCTOSE.*MANNOSE"
      ) ~
        "Fructose / mannose metabolism",
      
      str_detect(
        pathway_id,
        "PENTOSE.*PHOSPHATE"
      ) ~
        "Pentose-phosphate pathway",
      
      str_detect(
        pathway_id,
        "INSULIN.*SIGNAL"
      ) ~
        "Insulin signalling",
      
      str_detect(
        pathway_id,
        "PURINE.*METABOLISM"
      ) ~
        "Purine metabolism",
      
      str_detect(
        pathway_id,
        "STARCH.*SUCROSE"
      ) ~
        "Starch / sucrose metabolism",
      
      TRUE ~
        NA_character_
    )
  ) %>%
  
  filter(
    !is.na(pathway)
  ) %>%
  
  transmute(
    
    dataset_id,
    contrast,
    
    pathway_set =
      "ALS-reference pathways",
    
    pathway,
    
    NES =
      as.numeric(NES),
    
    padj =
      as.numeric(padj)
  )


# ============================================================
# 7. Hallmark pathways
#
# More flexible matching than before so it works whether the
# CSV contains "HALLMARK_GLYCOLYSIS", "Glycolysis", etc.
# ============================================================

bio_dat <- bio_raw %>%
  
  mutate(
    
    dataset_id =
      clean_dataset(dataset),
    
    pathway_id =
      clean_name(bioenergetic_pathway)
    
  ) %>%
  
  mutate(
    
    pathway = case_when(
      
      str_detect(
        pathway_id,
        "GLYCOLYSIS"
      ) ~
        "Hallmark glycolysis",
      
      str_detect(
        pathway_id,
        "OXIDATIVE.*PHOSPHORYLATION"
      ) ~
        "Oxidative phosphorylation",
      
      str_detect(
        pathway_id,
        "FATTY.*ACID.*METABOLISM"
      ) ~
        "Fatty-acid metabolism",
      
      TRUE ~
        NA_character_
    )
  ) %>%
  
  filter(
    !is.na(pathway)
  ) %>%
  
  transmute(
    
    dataset_id,
    contrast,
    
    pathway_set =
      "Hallmark bioenergetic pathways",
    
    pathway,
    
    NES =
      as.numeric(NES),
    
    padj =
      as.numeric(padj)
  )


# ============================================================
# 8. Combine
# ============================================================

all_pathways <- bind_rows(
  ref_dat,
  bio_dat
)


# ============================================================
# 9. Compact pathway labels
# ============================================================

pathway_short <- c(
  
  "Glycolysis / gluconeogenesis" =
    "Gly/GNG",
  
  "Fructose / mannose metabolism" =
    "Fru/Man",
  
  "Pentose-phosphate pathway" =
    "PPP",
  
  "Insulin signalling" =
    "Insulin",
  
  "Purine metabolism" =
    "Purine",
  
  "Starch / sucrose metabolism" =
    "Starch/Suc",
  
  "Hallmark glycolysis" =
    "H-Gly",
  
  "Oxidative phosphorylation" =
    "OXPHOS",
  
  "Fatty-acid metabolism" =
    "FA metab"
)


# ============================================================
# 10. Build reversal dataset
# ============================================================

make_reversal_data <- function(
    dat,
    dataset_name,
    unloading_contrast,
    restoration_contrast,
    final_contrast) {
  
  
  target_dataset <-
    clean_dataset(dataset_name)
  
  
  d0 <- dat %>%
    
    filter(
      dataset_id == target_dataset
    )
  
  
  if (nrow(d0) == 0) {
    
    stop(
      "Dataset not found: ",
      dataset_name,
      "\nAvailable datasets: ",
      paste(
        sort(unique(dat$dataset_id)),
        collapse = ", "
      )
    )
  }
  
  
  requested_contrasts <- c(
    unloading_contrast,
    restoration_contrast,
    final_contrast
  )
  
  
  missing_contrasts <- setdiff(
    requested_contrasts,
    unique(d0$contrast)
  )
  
  
  if (length(missing_contrasts) > 0) {
    
    stop(
      "\nMissing contrasts for ",
      dataset_name,
      ":\n",
      paste(
        missing_contrasts,
        collapse = "\n"
      ),
      "\n\nAvailable contrasts:\n",
      paste(
        sort(unique(d0$contrast)),
        collapse = "\n"
      )
    )
  }
  
  
  d <- d0 %>%
    
    mutate(
      
      phase = case_when(
        
        contrast ==
          unloading_contrast ~
          "unloading",
        
        contrast ==
          restoration_contrast ~
          "restoration",
        
        contrast ==
          final_contrast ~
          "final",
        
        TRUE ~
          NA_character_
      )
    ) %>%
    
    filter(
      !is.na(phase)
    )
  
  
  # ----------------------------------------------------------
  # Check duplicates
  # ----------------------------------------------------------
  
  duplicates <- d %>%
    
    count(
      pathway_set,
      pathway,
      phase
    ) %>%
    
    filter(
      n > 1
    )
  
  
  if (nrow(duplicates) > 0) {
    
    print(duplicates)
    
    stop(
      "Duplicate pathway/contrast rows found."
    )
  }
  
  
  # ----------------------------------------------------------
  # Wide format
  # ----------------------------------------------------------
  
  d <- d %>%
    
    select(
      pathway_set,
      pathway,
      phase,
      NES,
      padj
    ) %>%
    
    pivot_wider(
      
      names_from =
        phase,
      
      values_from =
        c(
          NES,
          padj
        ),
      
      names_sep =
        "_"
    )
  
  
  # ----------------------------------------------------------
  # Interpretation
  # ----------------------------------------------------------
  
  d <- d %>%
    
    mutate(
      
      reversed_direction =
        NES_unloading *
        NES_restoration < 0,
      
      significant_unloading =
        !is.na(padj_unloading) &
        padj_unloading < 0.05,
      
      significant_restoration =
        !is.na(padj_restoration) &
        padj_restoration < 0.05,
      
      significant_reversal =
        reversed_direction &
        significant_unloading &
        significant_restoration,
      
      significant_final =
        !is.na(padj_final) &
        padj_final < 0.05,
      
      short_label =
        unname(
          pathway_short[pathway]
        ),
      
      plot_label =
        paste0(
          short_label,
          ifelse(
            significant_restoration,
            "*",
            ""
          )
        )
      
    )
  
  
  message(
    dataset_name,
    " rows: ",
    nrow(d)
  )
  
  
  message(
    "  ALS-reference: ",
    sum(
      d$pathway_set ==
        "ALS-reference pathways"
    )
  )
  
  
  message(
    "  Hallmark: ",
    sum(
      d$pathway_set ==
        "Hallmark bioenergetic pathways"
    )
  )
  
  
  d
}


# ============================================================
# 11. OSD-21
# ============================================================

osd21 <- make_reversal_data(
  
  dat =
    all_pathways,
  
  dataset_name =
    "OSD21",
  
  unloading_contrast =
    "HLU_vs_CTRL",
  
  restoration_contrast =
    "HLU_RL_vs_HLU",
  
  final_contrast =
    "HLU_RL_vs_CTRL"
  
)


# ============================================================
# 12. OSD-540
# ============================================================

osd540 <- make_reversal_data(
  
  dat =
    all_pathways,
  
  dataset_name =
    "OSD540",
  
  unloading_contrast =
    "ULLS_vs_Baseline",
  
  restoration_contrast =
    "AR_vs_ULLS",
  
  final_contrast =
    "AR_vs_Baseline"
  
)


# ============================================================
# 13. Check pathway counts
# ============================================================

cat("\n")
cat("OSD-21 pathways:\n")

print(
  osd21 %>%
    select(
      pathway_set,
      pathway
    )
)


cat("\n")
cat("OSD-540 pathways:\n")

print(
  osd540 %>%
    select(
      pathway_set,
      pathway
    )
)


# ============================================================
# 14. Save source data
# ============================================================

source_dat <- bind_rows(
  
  osd21 %>%
    mutate(
      dataset = "OSD-21"
    ),
  
  osd540 %>%
    mutate(
      dataset = "OSD-540"
    )
  
) %>%
  
  select(
    
    dataset,
    pathway_set,
    pathway,
    
    NES_unloading,
    padj_unloading,
    
    NES_restoration,
    padj_restoration,
    
    NES_final,
    padj_final,
    
    reversed_direction,
    significant_reversal,
    significant_final
    
  )


write.csv(
  
  source_dat,
  
  file.path(
    summary_dir,
    "Figure4_reversal_plot_source_data.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 15. Common x/y limits
# ============================================================

xy_values <- c(
  
  osd21$NES_unloading,
  osd21$NES_restoration,
  
  osd540$NES_unloading,
  osd540$NES_restoration
  
)


xy_limit <- max(
  abs(xy_values),
  na.rm = TRUE
)


xy_limit <- ceiling(
  (xy_limit + 0.35) * 2
) / 2


# ============================================================
# 16. Residual NES colour limits
# ============================================================

residual_values <- c(
  osd21$NES_final,
  osd540$NES_final
)


residual_limit <- max(
  abs(residual_values),
  na.rm = TRUE
)


residual_limit <- ceiling(
  residual_limit * 2
) / 2


if (residual_limit == 0) {
  residual_limit <- 1
}


# ============================================================
# 17. Plotting function
# ============================================================

make_reversal_panel <- function(
    dat,
    panel_title,
    panel_subtitle,
    x_title,
    y_title,
    xy_limit,
    residual_limit) {
  
  
  ggplot(
    
    dat,
    
    aes(
      x = NES_unloading,
      y = NES_restoration
    )
    
  ) +
    
    
    # --------------------------------------------------------
  # Zero lines
  # --------------------------------------------------------
  
  geom_hline(
    
    yintercept = 0,
    
    linewidth = 0.45,
    
    colour = "grey70"
    
  ) +
    
    
    geom_vline(
      
      xintercept = 0,
      
      linewidth = 0.45,
      
      colour = "grey70"
      
    ) +
    
    
    # --------------------------------------------------------
  # Equal-magnitude reversal line
  # --------------------------------------------------------
  
  geom_abline(
    
    intercept = 0,
    
    slope = -1,
    
    linetype = "dashed",
    
    linewidth = 0.8,
    
    colour = "#4C78A8"
    
  ) +
    
    
    # --------------------------------------------------------
  # Points
  # --------------------------------------------------------
  
  geom_point(
    
    aes(
      shape = pathway_set,
      fill = NES_final
    ),
    
    size = 6.0,
    
    colour = "black",
    
    stroke = 0.7
    
  ) +
    
    
    # --------------------------------------------------------
  # Compact labels
  #
  # Much shorter labels + stronger repulsion.
  # --------------------------------------------------------
  
  geom_text_repel(
    
    aes(
      label = plot_label
    ),
    
    seed = 42,
    
    size = 4.6,
    
    fontface = "plain",
    
    box.padding = 0.8,
    
    point.padding = 0.9,
    
    min.segment.length = 0,
    
    segment.colour = "grey55",
    
    segment.size = 0.4,
    
    max.overlaps = Inf,
    
    force = 3.2,
    
    force_pull = 0.2,
    
    max.time = 5,
    
    direction = "both"
    
  ) +
    
    
    # --------------------------------------------------------
  # Shape
  # --------------------------------------------------------
  
  scale_shape_manual(
    
    values = c(
      
      "ALS-reference pathways" =
        21,
      
      "Hallmark bioenergetic pathways" =
        24
      
    ),
    
    name = NULL
    
  ) +
    
    
    # --------------------------------------------------------
  # Fill
  # --------------------------------------------------------
  
  scale_fill_gradient2(
    
    low =
      "#D6604D",
    
    mid =
      "#F7F7F7",
    
    high =
      "#4393C3",
    
    midpoint =
      0,
    
    limits = c(
      -residual_limit,
      residual_limit
    ),
    
    oob =
      scales::squish,
    
    name =
      paste0(
        "NES after restoration\n",
        "vs baseline/control"
      )
    
  ) +
    
    
    # --------------------------------------------------------
  # Axes
  # --------------------------------------------------------
  
  coord_equal(
    
    xlim = c(
      -xy_limit,
      xy_limit
    ),
    
    ylim = c(
      -xy_limit,
      xy_limit
    ),
    
    clip = "off"
    
  ) +
    
    
    scale_x_continuous(
      
      breaks =
        pretty(
          c(
            -xy_limit,
            xy_limit
          ),
          n = 5
        )
      
    ) +
    
    
    scale_y_continuous(
      
      breaks =
        pretty(
          c(
            -xy_limit,
            xy_limit
          ),
          n = 5
        )
      
    ) +
    
    
    # --------------------------------------------------------
  # Labels
  # --------------------------------------------------------
  
  labs(
    
    title =
      panel_title,
    
    subtitle =
      panel_subtitle,
    
    x =
      x_title,
    
    y =
      y_title
    
  ) +
    
    
    # --------------------------------------------------------
  # Theme
  # --------------------------------------------------------
  
  theme_classic(
    base_size = 16
  ) +
    
    theme(
      
      plot.title =
        element_text(
          size = 19,
          face = "bold",
          hjust = 0
        ),
      
      plot.subtitle =
        element_text(
          size = 15,
          hjust = 0,
          margin = margin(
            b = 14
          )
        ),
      
      axis.title.x =
        element_text(
          size = 16,
          face = "bold",
          margin = margin(
            t = 12
          )
        ),
      
      axis.title.y =
        element_text(
          size = 16,
          face = "bold",
          margin = margin(
            r = 12
          )
        ),
      
      axis.text =
        element_text(
          size = 13.5,
          colour = "black"
        ),
      
      axis.ticks =
        element_line(
          colour = "black"
        ),
      
      legend.title =
        element_text(
          size = 14,
          face = "bold"
        ),
      
      legend.text =
        element_text(
          size = 13
        ),
      
      legend.key.height =
        unit(
          0.75,
          "cm"
        ),
      
      legend.key.width =
        unit(
          1.0,
          "cm"
        ),
      
      plot.margin =
        margin(
          t = 14,
          r = 35,
          b = 18,
          l = 28
        )
      
    )
}


# ============================================================
# 18. Panel A
# ============================================================

p_A <- make_reversal_panel(
  
  dat =
    osd21,
  
  panel_title =
    "A   OSD-21 - Mouse",
  
  panel_subtitle =
    "12 days hindlimb unloading \u2192 3.5 h reloading",
  
  x_title =
    "NES during unloading\nHLU vs control",
  
  y_title =
    "NES during restoration\nReloading vs HLU",
  
  xy_limit =
    xy_limit,
  
  residual_limit =
    residual_limit
  
)


# ============================================================
# 19. Panel B
# ============================================================

p_B <- make_reversal_panel(
  
  dat =
    osd540,
  
  panel_title =
    "B   OSD-540 - Human",
  
  panel_subtitle =
    "10 days unilateral lower-limb suspension \u2192 active recovery",
  
  x_title =
    "NES during unloading\nULLS vs baseline",
  
  y_title =
    "NES during restoration\nRecovery vs ULLS",
  
  xy_limit =
    xy_limit,
  
  residual_limit =
    residual_limit
  
)


# ============================================================
# 20. Assemble
# ============================================================

fig4 <- p_A + p_B +
  
  plot_layout(
    
    ncol = 2,
    
    guides = "collect"
    
  ) &
  
  theme(
    
    legend.position =
      "bottom",
    
    legend.box =
      "horizontal"
    
  )


# ============================================================
# 21. Display
# ============================================================

print(fig4)


# ============================================================
# 22. Save
# ============================================================

ggsave(
  
  file.path(
    figure_dir,
    "Figure4_reversal_of_unloading_pathways.pdf"
  ),
  
  fig4,
  
  width = 14.5,
  height = 7.7,
  units = "in"
  
)


ggsave(
  
  file.path(
    figure_dir,
    "Figure4_reversal_of_unloading_pathways.png"
  ),
  
  fig4,
  
  width = 14.5,
  height = 7.7,
  units = "in",
  
  dpi = 600
  
)


ggsave(
  
  file.path(
    figure_dir,
    "Figure4_reversal_of_unloading_pathways.svg"
  ),
  
  fig4,
  
  width = 14.5,
  height = 7.7,
  units = "in"
  
)


# ============================================================
# 23. Summary
# ============================================================

reversal_summary <- bind_rows(
  
  osd21 %>%
    mutate(
      dataset = "OSD-21"
    ),
  
  osd540 %>%
    mutate(
      dataset = "OSD-540"
    )
  
) %>%
  
  group_by(
    dataset
  ) %>%
  
  summarise(
    
    n_pathways =
      n(),
    
    n_directional_reversal =
      sum(
        reversed_direction,
        na.rm = TRUE
      ),
    
    n_significant_reversal =
      sum(
        significant_reversal,
        na.rm = TRUE
      ),
    
    n_significant_vs_final_baseline =
      sum(
        significant_final,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  )


print(
  reversal_summary
)


write.csv(
  
  reversal_summary,
  
  file.path(
    summary_dir,
    "Figure4_reversal_summary.csv"
  ),
  
  row.names = FALSE
  
)
