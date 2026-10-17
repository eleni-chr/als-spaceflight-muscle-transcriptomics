# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# Figure 9 and Supplementary Figure S5: human ALS module validation.
# Exports PNG and SVG only; retains all sensitivity and CAMERA results.
# Run after human_ALS_module_validation_GSE215364_GSE41414.R.
# Uses exported results only; does not rerun or change the analysis.

RESULTS_DIR <- file.path(getwd(), 'human_ALS_module_validation')
WORKBOOK <- file.path(RESULTS_DIR, 'human_ALS_module_validation_results.xlsx')
OUTDIR <- file.path(RESULTS_DIR, 'publication_figures')
# Optional: change WORKBOOK above if you moved the results workbook.
FIGURE_NUMBER <- '9'
BASE_SIZE <- 14

packages <- c('readxl', 'dplyr', 'tidyr', 'ggplot2', 'patchwork', 'ragg', 'svglite')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos='https://cloud.r-project.org')
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
if (!file.exists(WORKBOOK)) stop('Workbook not found: ', WORKBOOK)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# Read explicit cell ranges: the openxlsx export can contain incorrect A1
# worksheet dimension metadata, despite having the full data in its XML.
index <- readxl::read_excel(WORKBOOK, sheet='INDEX', range='A1:B1000') |>
  filter(!is.na(sheet), !is.na(table))
read_table <- function(name, cols, rows=20000L) {
  hit <- index$sheet[index$table == name]
  if (length(hit) != 1L) stop('Expected one indexed table: ', name)
  x <- readxl::read_excel(WORKBOOK, sheet=hit,
    range=readxl::cell_limits(c(1L,1L), c(rows,cols)),
    guess_max=20000L)
  x[rowSums(!is.na(x)) > 0L, , drop=FALSE]
}
direction <- read_table('directional_summary',10,1000)
gsea <- read_table('fgsea_all_models',12,1000)
camera <- read_table('CAMERA_all_models',7,1000)
perm <- read_table('expression_matched_competitive_tests',12,1000)
cross <- read_table('cross_cohort_gene_effects',15)
cross_summary <- read_table('cross_cohort_summary',5,1000)
primary <- c('GSE215364_primary','GSE41414_primary')
models <- c(primary, 'GSE215364_covariate_adjusted',
  'GSE41414_covariate_adjusted','GSE41414_muscle_adjusted','GSE41414_deltoid_only')
model_labels <- c('GSE215364: primary','GSE41414: primary',
  'GSE215364: covariate-adjusted','GSE41414: covariate-adjusted',
  'GSE41414: muscle-adjusted','GSE41414: deltoid only')
component_labels <- c(ALL='Full module', SOD1_DOWN='SOD1G93A-down component',
  SOD1_UP='SOD1G93A-up component')
cohort_labels <- c(GSE215364_primary='GSE215364\nRNA-seq',
  GSE41414_primary='GSE41414\nMicroarray')
colours <- c(SOD1_DOWN='#2878B5',SOD1_UP='#C54B40')
fmt <- function(x) {
  vapply(x,function(z) {
    if (!is.finite(z)) return('NA')
    if (z < 0.001) formatC(z,format='e',digits=2)
    else formatC(z,format='f',digits=3)
  },character(1))
}
theme_fig <- function() theme_classic(base_size=BASE_SIZE) + theme(
  plot.title=element_text(face='bold',size=BASE_SIZE+2),
  plot.subtitle=element_text(size=BASE_SIZE-2,margin=margin(b=10)),
  axis.title=element_text(size=BASE_SIZE),
  strip.background=element_blank(),strip.text=element_text(face='bold'),
  legend.title=element_blank(),legend.position='bottom',
  plot.margin=margin(12,16,12,12))

# A: concordance, discordance and small effects among evaluable genes.
da <- direction |> filter(analysis %in% primary, threshold==0.1) |>
  mutate(cohort=factor(analysis,levels=primary,labels=unname(cohort_labels)),
    component=factor(component,levels=c('ALL','SOD1_DOWN','SOD1_UP'),
      labels=unname(component_labels)),
    label=sprintf('%d/%d\n%.1f%%',n_concordant,n_evaluable,100*fraction_concordant))
bars <- da |> pivot_longer(c(n_concordant,n_discordant,n_neutral),
  names_to='status', values_to='n') |>
  mutate(status=factor(status,levels=c('n_concordant','n_discordant','n_neutral'),
    labels=c('Concordant','Discordant','Below threshold')),percent=100*n/n_evaluable)
pA <- ggplot(bars,aes(cohort,percent,fill=status)) +
  geom_col(width=0.65,position=position_stack(reverse=TRUE)) +
  geom_text(data=da,aes(cohort,100*fraction_concordant/2,label=label),
    inherit.aes=FALSE,colour='white',fontface='bold',size=4.1) +
  facet_wrap(~component,nrow=1) +
  scale_fill_manual(values=c(Concordant='#347B70',Discordant='#B75D63',
    'Below threshold'='#D8D8D8'),drop=FALSE) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=expansion(mult=c(0,0.02))) +
  labs(title='A   Directional concordance',subtitle='|log2 fold-change| ≥ 0.10; denominators are evaluable genes',
    x=NULL,y='Evaluable genes (%)') + theme_fig()

# B: signed NES; labels use BH correction across four primary tests.
gb <- gsea |> filter(analysis %in% primary) |>
  mutate(cohort=factor(analysis,levels=primary,labels=unname(cohort_labels)),
    label=paste0('FDR ',fmt(FDR_four_primary_tests),
      ifelse(FDR_four_primary_tests<0.05,' *','')))
pB <- ggplot(gb,aes(cohort,NES,fill=pathway)) +
  geom_hline(yintercept=0,colour='grey45',linewidth=0.4) +
  geom_col(position=position_dodge(width=0.75),width=0.6) +
  geom_text(aes(y=NES+ifelse(NES>0,0.25,-0.25),label=label),
    position=position_dodge(width=0.75),size=3.7) +
  scale_fill_manual(values=colours,breaks=names(colours),
    labels=unname(component_labels[names(colours)])) +
  scale_y_continuous(limits=c(-4.1,3.6),breaks=-4:3) +
  labs(title='B   Directional gene-set enrichment',
    subtitle='* FDR < 0.05 across four primary tests',x=NULL,y='Normalised enrichment score (NES)') + theme_fig()

# C: observed versus mean concordance of expression-matched random sets.
pc <- perm |> filter(analysis %in% primary) |>
  mutate(cohort=factor(analysis,levels=primary,labels=unname(cohort_labels)))
pl <- pc |> dplyr::select(cohort,Observed=observed_concordance_0.1,
  'Matched null mean'=null_concordance_mean) |>
  pivot_longer(-cohort,names_to='series',values_to='fraction')
pC <- ggplot(pl,aes(cohort,100*fraction,fill=series)) +
  geom_col(position=position_dodge(width=0.75),width=0.65) +
  geom_text(aes(label=sprintf('%.1f%%',100*fraction)),
    position=position_dodge(width=0.75),vjust=-0.5,size=4.2) +
  geom_text(data=pc,aes(cohort,98,label=paste0('Empirical P=',fmt(p_concordance_expected))),
    inherit.aes=FALSE,size=3.9) +
  scale_fill_manual(values=c(Observed='#347B70','Matched null mean'='#BABABA')) +
  scale_y_continuous(limits=c(0,105),breaks=seq(0,100,25),expand=expansion(mult=c(0,0))) +
  labs(title='C   Agreement exceeds matched background',
    subtitle='20,000 expression-matched permutations; |log2FC| ≥ 0.10',
    x=NULL,y='Concordant genes (%)') + theme_fig()

# D: primary ALS-versus-control log2FC in the 99 shared module genes.
module_cross <- cross |> filter(module_gene==TRUE,
  is.finite(log2FC_RNA),is.finite(log2FC_array))
cs <- cross_summary |> filter(scope=='shared_module_genes')
if (nrow(cs)!=1L || nrow(module_cross)!=cs$n) stop('Cross-cohort coverage mismatch.')
rho <- cor(module_cross$log2FC_RNA,module_cross$log2FC_array,method='spearman')
if (abs(rho-cs$rho)>1e-8) stop('Cross-cohort correlation mismatch.')
pD <- ggplot(module_cross,aes(log2FC_RNA,log2FC_array,colour=component)) +
  geom_hline(yintercept=0,linetype='dashed',colour='grey65',linewidth=0.4) +
  geom_vline(xintercept=0,linetype='dashed',colour='grey65',linewidth=0.4) +
  geom_point(size=2.8,alpha=0.85) +
  scale_colour_manual(values=colours,breaks=names(colours),
    labels=unname(component_labels[names(colours)])) +
  labs(title='D   Cross-cohort module agreement',
    subtitle=sprintf('n=%d; Spearman rho=%.3f; same direction %.1f%%',
      cs$n,cs$rho,100*cs$fraction_same_direction),
    x='GSE215364 ALS vs control\n(log2 fold-change)',y='GSE41414 ALS vs control\n(log2 fold-change)') + theme_fig()
main <- (pA | pB) / (pC | pD) + plot_layout(heights=c(1,1))

# Supplement: keep all six models and nonsignificant CAMERA findings visible.
model_factor <- function(x) factor(x,levels=rev(models),labels=rev(model_labels))
sg <- gsea |> mutate(model=model_factor(analysis),
  component=factor(pathway,levels=names(colours),labels=unname(component_labels[names(colours)])),
  label=paste0('NES ',sprintf('%.2f',NES),'\nFDR ',fmt(padj)))
pSA <- ggplot(sg,aes(component,model,fill=NES)) +
  geom_tile(colour='white',linewidth=1) + geom_text(aes(label=label),size=3.6) +
  scale_fill_gradient2(low='#6CA4CD',mid='white',high='#D98E87',midpoint=0,limits=c(-3.5,3.5)) +
  labs(title='A   GSEA across primary and sensitivity models',
    subtitle='FDR corrected across the two components within each model',x=NULL,y=NULL,fill='NES') +
  theme_fig() + theme(axis.text.x=element_text(size=12),legend.position='none')
sc <- camera |> mutate(model=model_factor(analysis),
  component=factor(component,levels=names(colours),labels=unname(component_labels[names(colours)])),
  significance=ifelse(FDR<0.05,'FDR < 0.05','FDR ≥ 0.05'),
  label=paste0(Direction,'\nFDR ',fmt(FDR)))
pSB <- ggplot(sc,aes(component,model,fill=significance)) +
  geom_tile(colour='white',linewidth=1) + geom_text(aes(label=label),size=3.7) +
  scale_fill_manual(values=c('FDR < 0.05'='#B8D7CD','FDR ≥ 0.05'='#ECECEC')) +
  labs(title='B   Correlation-aware CAMERA tests',
    subtitle='FDR corrected within each model; Up/Down is the observed direction',x=NULL,y=NULL) + theme_fig()
sd <- direction |> filter(component=='ALL',threshold==0.1) |>
  mutate(model=model_factor(analysis),label=sprintf('%d/%d (%.1f%%)',
    n_concordant,n_evaluable,100*fraction_concordant))
pSC <- ggplot(sd,aes(100*fraction_concordant,model)) +
  geom_col(fill='#347B70',width=0.65) +
  geom_text(aes(label=label),hjust=-0.1,size=4) +
  scale_x_continuous(limits=c(0,110),breaks=seq(0,100,25),expand=expansion(mult=c(0,0))) +
  labs(title='C   Full-module directional concordance',
    subtitle='|log2 fold-change| ≥ 0.10; denominators are evaluable genes',x='Concordant genes (%)',y=NULL) + theme_fig()
supplement <- pSA / pSB / pSC + plot_layout(heights=c(1,1,1))

save_figure <- function(plot,stem,width,height) {
  ggsave(file.path(OUTDIR,paste0(stem,'.png')),plot,width=width,height=height,
    units='in',dpi=400,bg='white',limitsize=FALSE,device=ragg::agg_png)
  ggsave(file.path(OUTDIR,paste0(stem,'.svg')),plot,width=width,height=height,
    units='in',device=svglite::svglite,bg='white')
}
save_figure(main,paste0('Figure_',FIGURE_NUMBER,'_human_ALS_validation'),17,11)
save_figure(supplement,'Supplementary_Figure_S5_human_ALS_sensitivity',12,15)
write.csv(da,file.path(OUTDIR,'Figure8_directional_concordance_source.csv'),row.names=FALSE)
write.csv(gb,file.path(OUTDIR,'Figure8_GSEA_source.csv'),row.names=FALSE)
write.csv(pc,file.path(OUTDIR,'Figure8_permutation_source.csv'),row.names=FALSE)
write.csv(module_cross,file.path(OUTDIR,'Figure8_cross_cohort_source.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(OUTDIR,'figure_sessionInfo.txt'))
print(main)
print(supplement)
message('PNG and SVG figures saved in: ',normalizePath(OUTDIR,winslash='/'))
