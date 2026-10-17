# Analysis of Muscle Transcriptomic Signatures in Spaceflight and Mechanical Unloading

Analysis code and documentation covering skeletal-muscle transcriptomic and proteomic analyses across ALS, spaceflight, mechanical unloading, recovery, exercise, and denervation.

## Manuscript

This repository contains all the data analysis code and figure generation code accompanying a manuscript submitted for publication.

**Publication citation:** _Coming soon._

**Code author:** Dr Eleni Christoforidou

## Contents

The scripts are numbered in the order of the analysis workflow.

| Script | Analysis or figure |
|---|---|
| `01_als_metabolic_reference_figure.R` | Figure 1: SOD1G93A metabolic pathway reference |
| `02_spaceflight_unloading_reference_analysis.R` | Core differential-expression and pathway analyses supporting Figures 4–5 and downstream module analyses |
| `03_unloading_pathway_figure.R` | Figure 4: metabolic pathway responses to mechanical unloading |
| `04_unloading_recovery_figure.R` | Figure 5: pathway responses during restoration of mechanical loading |
| `05_reversible_module_and_exercise_analysis.R` | Definition and robustness analysis of the ALS–unloading module; analysis of human ALS resistance-exercise data |
| `06_reversible_module_figure.R` | Figure 7: recurrent and reversible ALS–unloading module |
| `07_actual_spaceflight_analysis.R` | Integrated NASA and JAXA spaceflight analyses |
| `08_spaceflight_figures_2_3_and_s1.R` | Figures 2–3 and Supplementary Figure S1 |
| `09_unloading_genomewide_concordance_analysis.R` | Genome-wide concordance across human and mouse unloading datasets |
| `10_unloading_genomewide_concordance_figure.R` | Figure 6: unloading-associated genome-wide concordance |
| `11_spaceflight_module_figure_8_and_s2.R` | Figure 8 and Supplementary Figure S2: fixed ALS module in spaceflight |
| `12_human_als_module_validation_analysis.R` | Human ALS skeletal-muscle validation in GSE215364 and GSE41414 |
| `13_human_als_validation_figures_9_and_s5.R` | Figure 9 and Supplementary Figure S5 |
| `14_als_resistance_exercise_figure_10.R` | Figure 10: resistance exercise in human ALS muscle |
| `15_denervation_module_analysis.R` | Denervation analysis in GSE217576 and GSE117440 |
| `16_denervation_figure_11.R` | Figure 11: denervation and electrical stimulation |
| `17_astronaut_muscle_proteomics_analysis.R` | Exploratory analysis of PXD028435 and PXD081616 astronaut muscle proteomics |
| `18_astronaut_proteomics_figures_s3_s4.R` | Supplementary Figures S3–S4 |

Some scripts also write quality-control tables or sensitivity outputs that are not shown as figures in the manuscript.

## Data and execution

The input datasets and intermediate analysis outputs are not included here. Source data should be obtained from the repositories and supplementary files cited in the manuscript. Dataset identifiers used by these scripts include NASA OSDR studies (including OSD-21, OSD-101, OSD-104, OSD-401, OSD-419, OSD-540, OSD-665, OSD-666, OSD-684, OSD-770, and OSD-111), JAXA MHU-1/MHU-4/MHU-5, GEO series GSE271607, GSE215364, GSE41414, GSE250455, GSE217576, and GSE117440, and ProteomeXchange datasets PXD028435 and PXD081616.

Before running a script, read its configuration section and place the required source files in the specified project folders. Paths are relative to the working directory by default; update the configuration values if your local layout differs. Run upstream analysis scripts before their figure scripts:

- Run 02 before 03 and 04.
- Run 02 and 05 before the downstream module analyses that need their outputs.
- Run 07 before 08 and 11.
- Run 09 before 10.
- Run 12 before 13.
- Run 15 before 16.
- Run 17 before 18.

`06_reversible_module_figure.R` uses objects created by script 05; source it in the same R session after script 05.

The scripts use the R packages listed in their configuration or package-loading sections. Package versions and operating-system differences may affect reproducibility. No analysis was run as part of preparing this release folder because the required input datasets are not included here.

## Peer review and use

The manuscript will undergo peer review; this repository has not been peer-reviewed. Review of the manuscript does not validate this code. Users should treat the code and derived results cautiously, inspect the data-processing choices, and independently verify any result they rely on.

## Disclaimer

The scripts and associated outputs are provided “as is” for research and educational use, without warranty of any kind. The author accepts no responsibility or liability for errors, omissions, interpretations, decisions, or consequences arising from use of this repository. Users are responsible for assessing the code’s suitability for their own work and for verifying results independently.

## License

This code is released under the Apache License 2.0. Please see the LICENSE file for details. By using this code, you agree to the terms of the Apache License 2.0, which grants explicit patent rights and requires that any redistributed or derivative works include the original copyright and attribution notices.
