# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# ============================================================================
# FIGURE 10 — Human ALS resistance exercise and the recurrent ALS–unloading module
# ============================================================================
# This script uses the GSE250455 outputs from
# 05_reversible_module_and_exercise_analysis.R
#
# Required for Panel C:
#   reversible_ALS_unloading_results/tables/
#     GSE250455_module_sample_score_pairs.csv
# Columns: participant, Pre, Post, delta_Post_minus_Pre
#
# Optional verification file (uploaded previously):
#   GSE250455_module_sample_score_summary.csv
#
# IMPORTANT: The analysis script already defines scores such that higher values
# indicate movement in the PREDICTED RECOVERY direction. DO NOT flip signs here.
# A negative Post - Pre score change is therefore the NON-RECOVERY direction.
# A: 119 of 167 module genes show further ALS-associated changes
# B: SOD1-down module GSEA: NES = -2.50, FDR = 2.7e-12
# C: Seven actual pairs of participant module scores (no simulated values)
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)

# ---------------- INPUTS -----------------------------------------------------
# Works from the project root OR from its results/tables folder.
PAIR_CANDIDATES <- c(
  file.path("reversible_ALS_unloading_results", "tables",
            "GSE250455_module_sample_score_pairs.csv"),
  file.path("tables", "GSE250455_module_sample_score_pairs.csv"),
  "GSE250455_module_sample_score_pairs.csv"
)
SUM_CANDIDATES <- c(
  "GSE250455_module_sample_score_summary.csv",
  file.path("reversible_ALS_unloading_results", "tables",
            "GSE250455_module_sample_score_summary.csv"),
  file.path("tables", "GSE250455_module_sample_score_summary.csv")
)

locate_file <- function(candidates, required = TRUE) {
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(found[[1]])
  if (required) stop(
    "Cannot find GSE250455_module_sample_score_pairs.csv.\n",
    "Run the v16 analysis first. It exports the file under ",
    "reversible_ALS_unloading_results/tables/.\n",
    "Current working directory: ", getwd(), "\nSearched:\n  ",
    paste(candidates, collapse = "\n  "), call. = FALSE
  )
  NULL
}

PAIRS_FILE <- locate_file(PAIR_CANDIDATES)
SUMMARY_FILE <- locate_file(SUM_CANDIDATES, required = FALSE)

pairs_raw <- read_csv(PAIRS_FILE, show_col_types = FALSE)
required <- c("participant", "Pre", "Post", "delta_Post_minus_Pre")
missing <- setdiff(required, names(pairs_raw))
if (length(missing)) stop(
  "The paired-score CSV is missing column(s): ", paste(missing, collapse = ", "),
  ".\nObserved: ", paste(names(pairs_raw), collapse = ", "), call. = FALSE
)

# Read precisely the data columns exported by module_expression_score() in v16.
pairs <- pairs_raw %>%
  transmute(
    participant = as.character(participant),
    Pre  = as.numeric(Pre),
    Post = as.numeric(Post),
    delta_reported = as.numeric(delta_Post_minus_Pre)
  ) %>%
  arrange(participant) %>%
  mutate(
    delta = Post - Pre,
    response = case_when(
      delta < 0 ~ "Non-recovery direction",
      delta > 0 ~ "Predicted recovery",
      TRUE ~ "No change"
    ),
    participant_label = paste0("P", row_number())
  )

if (nrow(pairs) != 7L || anyDuplicated(pairs$participant) ||
    anyNA(pairs) || !all(is.finite(pairs$Pre) & is.finite(pairs$Post))) {
  stop("The paired-score table must contain exactly seven unique, complete participants.")
}
if (any(abs(pairs$delta - pairs$delta_reported) > 1e-7)) {
  stop("Post - Pre does not agree with delta_Post_minus_Pre in paired-score CSV.")
}

n_nonrecovery <- sum(pairs$delta < 0)
n_recovery <- sum(pairs$delta > 0)
n_unchanged <- sum(pairs$delta == 0)
if (n_nonrecovery != 5L || n_recovery != 2L || n_unchanged != 0L) {
  stop("The measured paired scores do NOT reproduce the reported 5/7 non-recovery ",
       "and 2/7 recovery directions. Check that you used the v16 module scores.")
}

# Calculate actual summary values and p-values from the paired measurements.
mean_delta <- mean(pairs$delta)
median_delta <- median(pairs$delta)
t_p <- t.test(pairs$Post, pairs$Pre, paired = TRUE,
              alternative = "two.sided")$p.value
wilcox_p <- suppressWarnings(
  wilcox.test(pairs$Post, pairs$Pre, paired = TRUE,
              alternative = "two.sided", exact = FALSE)$p.value
)

# Confirm this really is the same run represented by the previously supplied
# one-row summary CSV, where available.
if (!is.null(SUMMARY_FILE)) {
  summary <- read_csv(SUMMARY_FILE, show_col_types = FALSE)
  fields <- c("n_participants", "n_module_genes_measured", "mean_delta",
              "median_delta", "paired_t_p_two_sided",
              "paired_wilcox_p_two_sided")
  if (nrow(summary) != 1L || !all(fields %in% names(summary))) {
    stop("Unexpected one-row summary CSV structure.")
  }
  checks <- c(
    nrow(pairs) == summary$n_participants[[1]],
    summary$n_module_genes_measured[[1]] == 167,
    isTRUE(all.equal(mean_delta, summary$mean_delta[[1]], tolerance = 1e-7)),
    isTRUE(all.equal(median_delta, summary$median_delta[[1]], tolerance = 1e-7)),
    isTRUE(all.equal(t_p, summary$paired_t_p_two_sided[[1]], tolerance = 1e-6)),
    isTRUE(all.equal(wilcox_p, summary$paired_wilcox_p_two_sided[[1]], tolerance = 1e-6))
  )
  if (!all(checks)) stop(
    "The participant-pair file and summary CSV do not match. ",
    "Use outputs from the same v16 analysis run."
  )
}

cat("\nPanel C — observed participant changes:\n")
print(pairs %>% dplyr::select(participant_label, Pre, Post, delta, response), n = Inf)
cat(sprintf(
  "\n%d non-recovery; %d predicted recovery; mean delta %.4f; ",
  n_nonrecovery, n_recovery, mean_delta
))
cat(sprintf("median delta %.4f; paired t P %.4f; Wilcoxon P %.4f\n\n",
            median_delta, t_p, wilcox_p))

# ---------------- MANUSCRIPT STATISTICS --------------------------------------
N_MODULE_TOTAL <- 174L
N_MEASURED <- 167L
N_NONRECOVERY_GENES <- 119L
N_RECOVERY_GENES <- 48L
SIGN_P <- 1.9e-8  # one-sided (opposite to predicted recovery)
MATCHED_P <- 0.00020
GSEA_NES <- -2.50
GSEA_FDR <- 2.7e-12
stopifnot(N_NONRECOVERY_GENES + N_RECOVERY_GENES == N_MEASURED)

# ---------------- STYLE ------------------------------------------------------
BASE_SIZE <- 16
COL_NONRECOVERY <- "#C95858"
COL_RECOVERY <- "#4C78A8"
COL_TEXT <- "#222222"
COL_MUTED <- "#686868"

th <- theme_classic(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    text = element_text(colour = COL_TEXT),
    plot.title = element_text(size = 17, face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 13.5, colour = COL_MUTED,
                                 margin = margin(b = 12)),
    plot.caption = element_text(size = 12.5, colour = COL_TEXT,
                                hjust = 0, margin = margin(t = 12)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 13, colour = COL_TEXT),
    plot.margin = margin(12, 18, 14, 14),
    legend.position = "none"
  )

# ---------------- PANEL A ----------------------------------------------------
pct_non <- 100 * N_NONRECOVERY_GENES / N_MEASURED
pct_rec <- 100 * N_RECOVERY_GENES / N_MEASURED

genes <- tibble(
  direction = c("Non-recovery", "Predicted recovery"),
  n = c(N_NONRECOVERY_GENES, N_RECOVERY_GENES),
  pct = c(pct_non, pct_rec),
  xmin = c(0, pct_non),
  xmax = c(pct_non, 100)
) %>% mutate(xmid = (xmin + xmax) / 2)

pA <- ggplot(genes) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = .82, ymax = 1.18,
                fill = direction), colour = "white", linewidth = 1.1) +
  geom_text(aes(x = xmid, y = 1,
                label = sprintf("%d / %d\n%.1f%%", n, N_MEASURED, pct)),
            colour = "white", fontface = "bold", size = 5.1,
            lineheight = 1.05) +
  annotate("segment", x = 50, xend = 50, y = .74, yend = 1.26,
           colour = COL_TEXT, linewidth = .7, linetype = "dashed") +
  annotate("text", x = 50, y = 1.34, label = "50%", size = 3.9,
           colour = COL_MUTED) +
  annotate("text", x = pct_non / 2, y = .65,
           label = "Further SOD1G93A/unloading\nchange", size = 4.1,
           fontface = "bold", colour = COL_NONRECOVERY, lineheight = 1.05) +
  annotate("text", x = pct_non + pct_rec / 2, y = .65,
           label = "Predicted\nreversal", size = 4.1,
           fontface = "bold", colour = COL_RECOVERY, lineheight = 1.05) +
  scale_fill_manual(values = c("Non-recovery" = COL_NONRECOVERY,
                               "Predicted recovery" = COL_RECOVERY)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"),
                     limits = c(0, 100), expand = c(0, 0)) +
  coord_cartesian(ylim = c(.4, 1.44), clip = "off") +
  labs(title = "Gene-level directional response",
       x = "Proportion of measurable module genes", y = NULL,
       caption = paste0(
         "Directional sign test (119 vs 48): one-sided P = 1.9 × 10⁻⁸\n",
         "Expression-matched random-gene-set test (mean signed statistic):\n",
         "empirical P = 0.00020"
       )) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank())

# ---------------- PANEL B ----------------------------------------------------
# NES summary/lollipop (not a running GSEA enrichment curve).
pB <- ggplot() +
  annotate("rect", xmin = -3.3, xmax = 0, ymin = .6, ymax = 1.62,
           fill = COL_NONRECOVERY, alpha = .075) +
  annotate("rect", xmin = 0, xmax = 3.3, ymin = .6, ymax = 1.62,
           fill = COL_RECOVERY, alpha = .075) +
  annotate("segment", x = 0, xend = 0, y = .60, yend = 1.62,
           linetype = "dashed", colour = COL_MUTED, linewidth = .75) +
  annotate("segment", x = 0, xend = GSEA_NES, y = 1, yend = 1,
           linewidth = 2.0, colour = COL_NONRECOVERY) +
  annotate("point", x = GSEA_NES, y = 1, size = 7.5,
           colour = COL_NONRECOVERY) +
  annotate("text", x = -2.13, y = .79, label = "NES = −2.50",
           fontface = "bold", size = 4.7, colour = COL_NONRECOVERY) +
  annotate("text", x = -2.13, y = .68, label = "FDR = 2.7 × 10⁻¹²",
           fontface = "bold", size = 4.0, colour = COL_TEXT) +
  annotate("text", x = -1.6, y = 1.47,
           label = "Further ALS/unloading\nchange", size = 3.7,
           fontface = "bold", colour = COL_NONRECOVERY, lineheight = 1.05) +
  annotate("text", x = 1.6, y = 1.47,
           label = "Predicted\nrecovery", size = 3.7,
           fontface = "bold", colour = COL_RECOVERY, lineheight = 1.05) +
  scale_x_continuous(limits = c(-3.3, 3.3), breaks = -3:3,
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(.6, 1.62), clip = "off") +
  labs(title = "Coordinated module response",
       x = "Normalised enrichment score (NES)", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank())

# ---------------- PANEL C ----------------------------------------------------
# Actual matched pre/post values for each participant. These scores were
# gene-wise standardised and oriented in v16, so HIGHER = predicted recovery.
# Change in score (Post - Pre), not absolute score relative to zero, defines
# the direction of each participant's response.
score_long <- pairs %>%
  dplyr::select(participant_label, Pre, Post, response) %>%
  pivot_longer(cols = c(Pre, Post),
               names_to = "timepoint", values_to = "score") %>%
  mutate(x_position = if_else(timepoint == "Pre", 1, 2))

pC <- ggplot(
  score_long,
  aes(x = x_position, y = score,
      group = participant_label, colour = response)
) +
  geom_line(linewidth = 1.25, alpha = .85) +
  geom_point(size = 3.8) +
  scale_colour_manual(
    values = c("Non-recovery direction" = COL_NONRECOVERY,
               "Predicted recovery" = COL_RECOVERY),
    limits = c("Non-recovery direction", "Predicted recovery"),
    labels = c("Non-recovery", "Predicted recovery")
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Before exercise", "After 3 months"),
    limits = c(0.85, 2.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(.12, .12))) +
  labs(
    title = "Individual participant module scores",
    x = NULL,
    y = "Recovery-oriented module score\n(mean gene-wise z score)",
    colour = NULL
  ) +
  th +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12.5),
    legend.key.width = grid::unit(1.1, "cm"),
    plot.margin = margin(12, 18, 14, 14)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.4, size = 3)))

# ---------------- COMBINE & EXPORT -------------------------------------------
figure8 <- ((pA | pB) / (pC | plot_spacer())) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 21, face = "bold", family = "sans"),
      plot.tag.position = c(0, 1)
    )
  )

print(figure8)

ggsave("Figure_10_GSE250455.pdf", figure8,
       width = 14, height = 11.2, units = "in")
ggsave("Figure_10_GSE250455.png", figure8,
       width = 14, height = 11.2, units = "in", dpi = 600, bg = "white")
ggsave("Figure_10_GSE250455.tiff", figure8,
       width = 14, height = 11.2, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

write_csv(pairs %>% dplyr::select(participant_label, participant, Pre, Post, delta, response),
          "Figure_10_GSE250455_PanelC_data.csv")
message("Figure 10 saved. Actual paired-score source: ", PAIRS_FILE)

# ------------------------------------------------------------
# Print and save publication figure-legend draft. Participant-level
# inferential tests belong here rather than cluttering Panel C.
# ------------------------------------------------------------
legend_text <- sprintf(paste0(
  "Figure 10. Resistance exercise in human ALS skeletal muscle does not ",
  "reverse the recurrent ALS–unloading programme. ",
  "(A) Among 167 measurable module genes, 119 (71.3%%) changed in the ",
  "ALS/unloading-associated direction and 48 (28.7%%) in the predicted ",
  "recovery direction after exercise. The one-sided binomial sign test ",
  "assessed this directional imbalance (P = 1.9 × 10⁻⁸); the separate ",
  "expression-matched permutation test evaluated the mean oriented ",
  "exercise statistic against random gene sets (empirical P = 0.00020). ",
  "(B) Pre-ranked GSEA for the SOD1-downregulated module component ",
  "showed negative enrichment following exercise (NES = −2.50; ",
  "FDR = 2.7 × 10⁻¹²), opposite to the predicted recovery response. ",
  "(C) Recovery-oriented module scores for seven paired pre-exercise ",
  "and three-month post-exercise ALS muscle biopsies. Each line ",
  "represents one participant; positive Post − Pre changes indicate ",
  "the predicted recovery direction. Five participants showed a ",
  "negative change and two showed a positive change. Mean Δ = %.3f; ",
  "median Δ = %.3f. Paired two-sided t-test P = %.3f; ",
  "two-sided Wilcoxon signed-rank test P = %.3f. " ,
  "These participant-level tests were not statistically significant."
), mean_delta, median_delta, t_p, wilcox_p)

writeLines(legend_text, "Figure_10_GSE250455_legend.txt", useBytes = TRUE)
cat("\n", legend_text, "\n", sep = "")
