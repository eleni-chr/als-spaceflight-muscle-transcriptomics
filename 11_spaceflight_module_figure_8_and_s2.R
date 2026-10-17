#!/usr/bin/env Rscript
# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# =============================================================================
# FIGURE 8 + SUPPLEMENTARY FIGURE S2 | ALS vs actual spaceflight
# Frozen 174-gene ALS-associated load-sensitive module
#
# Input: unified actual-spaceflight pipeline v9 (JAXA fixed) outputs.
#   UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results/module/
#      frozen_module_GSEA.csv
#      frozen_module_gene_directionality.csv
#      sample_scores.csv
#      sample_score_group_tests.csv
# or a results.zip that contains these four files in its module/ directory.
#
# Figure 8: 16 primary comparisons
#   A | NES for ALS_DOWN and ALS_UP; * = global BH FDR < 0.05 over 32 tests
#   B | Whole-module directional concordance (174-gene frozen reference)
#   C | Mean difference in signed per-sample module score and Welch 95% CI;
#       * = BH FDR < 0.05 over 16 primary score tests
# Supplementary Figure S2: MHU-1 MG/AG/GC comparisons
#   A | GSEA (primary and two separately adjusted secondary comparisons)
#   B | Whole-module directional concordance
#   C | Sample-score differences and Welch 95% CI
#
# Symbols in both figures: A circles ALS_DOWN / triangles ALS_UP;
# B diamonds show whole-module concordance; C squares show score differences.
# The script does NOT rerun differential expression, fgsea, or gene selection.
# It checks the full comparison inventory and refuses truncated / old results.
# =============================================================================

# ------------------------- CHANGE THESE IF NEEDED ----------------------------
PROJECT_DIR <- getwd()   # folder containing results.zip or unified results folder
RESULTS_DIR <- file.path(PROJECT_DIR, "UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results")
ZIP_FILE <- file.path(PROJECT_DIR, "results.zip")
OUTPUT_DIR <- file.path(PROJECT_DIR, "MANUSCRIPT_FIGURE_8_AND_S2")
MAIN_W <- 19.2; MAIN_H <- 11.0
SUPP_W <- 15.8; SUPP_H <- 6.2
DPI <- 400
BASE_FONT <- 13.6
ALPHA <- 0.05
# -----------------------------------------------------------------------------

for (pkg in c("ggplot2", "patchwork", "svglite")) {
  if (!requireNamespace(pkg, quietly=TRUE))
    stop("Please install ",pkg," first: install.packages('",pkg,"')",call.=FALSE)
}
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)

# Results may reside in the unified output folder OR in its original ZIP.
# Reject ambiguity: never quietly choose one of two similarly named analyses.
find_csv <- function(filename) {
  rel <- file.path("module", filename)
  candidates <- c(file.path(RESULTS_DIR,rel),file.path(PROJECT_DIR,rel),
                  file.path(PROJECT_DIR,filename))
  present <- candidates[file.exists(candidates)]
  if (length(present)) return(normalizePath(present[1],winslash="/"))
  # Search for a nested unified folder, but exclude any generated figure files.
  nested <- list.files(PROJECT_DIR,pattern=paste0("^",filename,"$"),
                       recursive=TRUE,full.names=TRUE)
  nested <- nested[!grepl("MANUSCRIPT_FIGURE_8_AND_S2",nested,fixed=TRUE)]
  if (length(nested)==1L) return(normalizePath(nested,winslash="/"))
  if (length(nested)>1L) stop("Ambiguous source: ",filename," appears in:\n",
                             paste(nested,collapse="\n"),call.=FALSE)
  zips <- c(ZIP_FILE,file.path(RESULTS_DIR,"results.zip"))
  zips <- unique(zips[file.exists(zips)])
  if (!length(zips)) stop("Missing input: ",rel,
      "\nPut unified results folder or its results.zip in PROJECT_DIR.",call.=FALSE)
  hits <- list()
  for (z in zips) {
    contents <- utils::unzip(z,list=TRUE)$Name
    contents <- gsub("\\\\","/",contents)
    matches <- contents[endsWith(contents,paste0("/",rel)) | contents==rel]
    for (m in matches) hits[[length(hits)+1L]] <- list(zip=z,member=m)
  }
  if (length(hits)!=1L) stop("Expected ONE ZIP entry ending module/",filename,
      "; found ",length(hits),". Remove duplicate ZIPs or specify RESULTS_DIR.",call.=FALSE)
  cache <- file.path(OUTPUT_DIR,"_input_cache")
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  utils::unzip(hits[[1]]$zip,files=hits[[1]]$member,
               exdir=cache,overwrite=TRUE)
  extracted <- file.path(cache,hits[[1]]$member)
  if (!file.exists(extracted)) stop("ZIP extraction failed: ",filename)
  normalizePath(extracted,winslash="/")
}

paths <- setNames(lapply(c("frozen_module_GSEA.csv",
                            "frozen_module_gene_directionality.csv",
                            "sample_scores.csv",
                            "sample_score_group_tests.csv"),find_csv),
                  c("gsea","genes","samples","tests"))
for (nm in names(paths)) message(nm,": ",paths[[nm]])
read_one <- function(p) utils::read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
M <- read_one(paths$gsea); G <- read_one(paths$genes)
SC <- read_one(paths$samples); ST <- read_one(paths$tests)

require_cols <- function(d,cols,name) {
  absent <- setdiff(cols,names(d))
  if (length(absent)) stop(name," is missing column(s): ",
                           paste(absent,collapse=", "),call.=FALSE)
}
require_cols(M,c("dataset","tissue","cohort","contrast","primary",
                 "pathway","NES","pval","FDR_primary_across_ALL_spaceflight",
                 "FDR_secondary_across_ALL_spaceflight"),"GSEA table")
require_cols(G,c("dataset","tissue","cohort","contrast","human_symbol",
                 "als_sign","log2FoldChange"),"Gene table")
require_cols(ST,c("dataset","tissue","cohort","contrast","primary",
                  "mean_FLT_minus_control","pvalue","FDR_primary_scores"),"Score tests")
require_cols(SC,c("dataset","tissue","cohort","contrast","sample","group",
                  "ALS_oriented_module_score"),"Sample scores")

as_bool <- function(x) {
  y <- toupper(trimws(as.character(x)))
  if (any(!is.na(x) & !y %in% c("TRUE","FALSE","T","F","1","0")))
    stop("Unexpected TRUE/FALSE encoding in a primary flag",call.=FALSE)
  ifelse(is.na(x),NA,y %in% c("TRUE","T","1"))
}
M$primary <- as_bool(M$primary); ST$primary <- as_bool(ST$primary)
if ("primary" %in% names(G)) G$primary <- as_bool(G$primary)

# Exact inventory, consistent with the unified actual-flight pipeline and
# existing manuscript Figure 2 (dataset + tissue + cohort + contrast).
registry <- data.frame(
  dataset=c("OSD101","OSD419","OSD104","OSD401","OSD401",
            "OSD770","OSD665","OSD666","OSD111","OSD111",
            "OSD21","OSD684","OSD684","MHU1","MHU4","MHU5"),
  tissue=c("GASTROCNEMIUS","GASTROCNEMIUS","SOLEUS",
           "GASTROCNEMIUS","QUADRICEPS","SOLEUS","EDL","QUADRICEPS",
           "SOLEUS","EDL","GASTROCNEMIUS","MYOBUNDLE","MYOBUNDLE",
           "SOLEUS","SOLEUS","SOLEUS"),
  cohort=c(rep("ALL",11),"YOUNG","OLD","ALL","ALL","ALL"),
  contrast=c(rep("FLT_vs_HGC",3),rep("FLT_vs_GC",2),
             rep("FLT_vs_HGC",5),"FLT_vs_GC","FLT_vs_GC","FLT_vs_GC",
             "MG_vs_AG","PG_vs_GC","PG_vs_GC"),
  label=c("RR-1  |  Gastrocnemius L (OSD-101)",
          "RR-1  |  Gastrocnemius R (OSD-419)",
          "RR-1  |  Soleus (OSD-104)",
          "RR-5  |  Gastrocnemius","RR-5  |  Quadriceps",
          "RR-23  |  Soleus","RR-23  |  EDL","RR-23  |  Quadriceps",
          "BION-M1  |  Soleus","BION-M1  |  EDL",
          "STS-108  |  Gastrocnemius",
          "CRS-21  |  Young human myobundles",
          "CRS-21  |  Older human myobundles",
          "MHU-1  |  Soleus, MG vs onboard 1 g",
          "MHU-4  |  Soleus, PG vs ground",
          "MHU-5  |  Soleus, PG vs ground"),
  archive=c(rep("NASA OSDR",13),rep("JAXA ibSLS",3)),
  stringsAsFactors=FALSE)
keyfun <- function(d) paste(d$dataset,d$tissue,d$cohort,d$contrast,sep="||")
registry$key <- keyfun(registry)
for (tab in list(M=M,G=G,SC=SC,ST=ST)) {
  if (anyNA(tab[,c("dataset","tissue","cohort","contrast")]))
    stop("Missing contrast-identification metadata",call.=FALSE)
}
M$key <- keyfun(M); G$key <- keyfun(G)
SC$key <- keyfun(SC); ST$key <- keyfun(ST)
M$pathway <- as.character(M$pathway)
components <- c("ALS_DOWN","ALS_UP")
MP <- M[M$primary,,drop=FALSE]
observed <- unique(MP$key)
if (length(observed)!=16L || !setequal(observed,registry$key))
  stop("The primary GSEA contrasts are NOT the required 16. Missing:\n",
       paste(setdiff(registry$key,observed),collapse="\n"),
       "\nExtra:\n",paste(setdiff(observed,registry$key),collapse="\n"),
       "\nUse complete, recent UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results.",call.=FALSE)
if (nrow(MP)!=32L || anyDuplicated(paste(MP$key,MP$pathway)) ||
    !setequal(unique(MP$pathway),components) ||
    any(table(MP$key)!=2L)) stop("Expected exactly 32 primary GSEA component rows")
if (any(!is.finite(MP$NES)) || any(!is.finite(MP$pval)) ||
    any(!is.finite(MP$FDR_primary_across_ALL_spaceflight)))
  stop("Non-finite primary GSEA NES/pvalue/FDR",call.=FALSE)
# Confirm the FDR column really belongs to the current 32-test family.
check_q <- p.adjust(MP$pval,method="BH")
if (any(abs(check_q-MP$FDR_primary_across_ALL_spaceflight)>1e-7))
  stop("Primary module FDR is not BH-adjusted over these 32 tests",call.=FALSE)

GP <- G[G$key %in% registry$key,,drop=FALSE]
SP <- ST[ST$primary,,drop=FALSE]
if (!setequal(unique(GP$key),registry$key) ||
    !setequal(unique(SP$key),registry$key) ||
    nrow(SP)!=16L || anyDuplicated(SP$key))
  stop("Gene directionality and score tables must cover all 16 primary contrasts",call.=FALSE)
if (any(!is.finite(SP$pvalue)) || any(!is.finite(SP$FDR_primary_scores)) ||
    any(abs(p.adjust(SP$pvalue,"BH")-SP$FDR_primary_scores)>1e-7))
  stop("Score-test FDR must be BH-adjusted over all 16 primary tests",call.=FALSE)

# Whole 174-gene module concordance: each measured symbol counts once, using
# direction of its finite log2FC, with NO per-gene significance threshold.
whole_concordance <- function(df) {
  k <- unique(df$key)
  do.call(rbind,lapply(k,function(ki) {
    g <- df[df$key==ki,,drop=FALSE]
    if (anyDuplicated(g$human_symbol))
      stop("Duplicate module symbol in ",ki,"; fix gene harmonisation",call.=FALSE)
    if (nrow(g)!=174L) stop("Expected 174 frozen module symbols for ",ki,
                           "; got ",nrow(g),call.=FALSE)
    if (any(!g$als_sign %in% c(-1,1))) stop("Invalid frozen ALS gene direction")
    good <- is.finite(g$log2FoldChange)
    if (!any(good)) stop("No module genes evaluable for ",ki)
    matched <- sign(g$log2FoldChange[good])==g$als_sign[good]
    data.frame(key=ki,n_evaluable=sum(good),n_concordant=sum(matched),
               pct_concordant=100*mean(matched),stringsAsFactors=FALSE)
  }))
}
CON <- whole_concordance(G)
if (any(CON$n_evaluable>174L | CON$n_evaluable<130L))
  stop("Implausible module gene coverage in directional concordance table")

# Rebuild sample-score confidence intervals from the *biological-sample* scores,
# not from the 174 gene observations. Welch t-test agrees with pipeline.
score_intervals <- function(tests,samples) {
  do.call(rbind,lapply(seq_len(nrow(tests)),function(i) {
    d <- tests[i,,drop=FALSE]; ss <- samples[samples$key==d$key,,drop=FALSE]
    if (anyDuplicated(ss$sample))
      stop("Duplicate biological sample in module scores: ",d$key)
    groups <- strsplit(d$contrast,"_vs_",fixed=TRUE)[[1]]
    if (length(groups)!=2L) stop("Unrecognised comparison: ",d$contrast)
    pos <- as.numeric(ss$ALS_oriented_module_score[ss$group==groups[1]])
    neg <- as.numeric(ss$ALS_oriented_module_score[ss$group==groups[2]])
    pos <- pos[is.finite(pos)];neg <- neg[is.finite(neg)]
    if (length(pos)<2L || length(neg)<2L)
      stop("Insufficient biological samples for score CI: ",d$key)
    t <- stats::t.test(pos,neg,conf.level=0.95,var.equal=FALSE)
    if (abs(unname(t$estimate[1]-t$estimate[2])-d$mean_FLT_minus_control)>1e-5)
      stop("Score difference in sample scores disagrees with group test: ",d$key)
    if (is.finite(d$pvalue) && abs(t$p.value-d$pvalue)>1e-5)
      stop("Score p-value disagrees with sample-level Welch test: ",d$key)
    data.frame(key=d$key,ci_low=t$conf.int[1],ci_high=t$conf.int[2],
               n_exposed=length(pos),n_control=length(neg),stringsAsFactors=FALSE)
  }))
}
# CI reconstruction is needed only for 16 primary comparisons plus the two
# supplementary MHU-1 controls; unrelated secondary contrasts must not prevent
# the main figure from rendering if they lack sample-level Welch estimates.
mhu1_secondary_keys <- paste("MHU1","SOLEUS","ALL",
                             c("MG_vs_GC","AG_vs_GC"),sep="||")
needed_score_keys <- c(registry$key,mhu1_secondary_keys)
ST <- ST[ST$key %in% needed_score_keys,,drop=FALSE]
SC <- SC[SC$key %in% needed_score_keys,,drop=FALSE]
if (anyDuplicated(ST$key)) stop("Duplicate sample-score tests for a contrast")
CI <- score_intervals(ST,SC)
ST <- merge(ST,CI,by="key",sort=FALSE)

# Propagate display metadata, preserving exact registry order.
attach_info <- function(x,reg=registry) {
  ix <- match(x$key,reg$key)
  if (anyNA(ix)) stop("Unregistered comparison in plotting data",call.=FALSE)
  x$label <- reg$label[ix];x$archive <- reg$archive[ix]
  x$y <- nrow(reg) - ix + 1
  x
}
MP <- attach_info(MP)
CP <- attach_info(CON[CON$key %in% registry$key,,drop=FALSE])
SP <- attach_info(ST[ST$primary,,drop=FALSE])
MP$component <- factor(MP$pathway,levels=components)
MP$y_plot <- MP$y+ifelse(MP$pathway=="ALS_DOWN",0.16,-0.16)
MP$asterisk <- ifelse(MP$FDR_primary_across_ALL_spaceflight<ALPHA,"*","")
SP$asterisk <- ifelse(SP$FDR_primary_scores<ALPHA,"*","")

COL_OSDR <- "#E69F00"; COL_JAXA <- "#CC79A7"
COL_INK <- "#30353B"; COL_GRID <- "#E8E8E8"
cols <- c("NASA OSDR"=COL_OSDR,"JAXA ibSLS"=COL_JAXA)
shapes <- c("ALS_DOWN"=21,"ALS_UP"=24)
sep_main <- nrow(registry)-c(3,5,8,10,11,13,14,15)+0.5

base_theme <- ggplot2::theme_classic(base_size=BASE_FONT) + ggplot2::theme(
  panel.grid.major.x=ggplot2::element_line(colour=COL_GRID,linewidth=.26),
  panel.grid.minor=ggplot2::element_blank(),
  axis.text.x=ggplot2::element_text(colour=COL_INK,size=BASE_FONT*.88),
  axis.text.y=ggplot2::element_text(colour=COL_INK,size=BASE_FONT*.84),
  axis.title.x=ggplot2::element_text(size=BASE_FONT*.92,margin=ggplot2::margin(t=11)),
  plot.title=ggplot2::element_text(face="bold",size=BASE_FONT*1.03,
                                 margin=ggplot2::margin(b=12)),
  legend.title=ggplot2::element_text(size=BASE_FONT*.9,face="bold"),
  legend.text=ggplot2::element_text(size=BASE_FONT*.84),
  legend.spacing.y=grid::unit(4,"pt"),
  plot.margin=ggplot2::margin(12,8,10,7))

horizontal_lines <- function(seps) {
  ggplot2::geom_hline(yintercept=seps,col=COL_GRID,linewidth=.53)
}
yscale <- function(n,labels=NULL) {
  ggplot2::scale_y_continuous(breaks=seq_len(n),
     labels=if (is.null(labels)) NULL else rev(labels),
     limits=c(.42,n+.58),expand=ggplot2::expansion(mult=0))
}
clean_y <- ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                          axis.ticks.y=ggplot2::element_blank(),
                          axis.line.y=ggplot2::element_blank())

# ------------------------------- MAIN PANEL A -------------------------------
neslimit <- max(3.25,ceiling(2*max(abs(MP$NES)))/2+0.3)
# Attach star to the right of each GSEA symbol; paired components occupy
# separate y offsets so their stars cannot be mistaken for another component.
MP$star_x <- MP$NES + 0.20
figA <- ggplot2::ggplot(MP,ggplot2::aes(x=NES,y=y_plot))+
  horizontal_lines(sep_main)+
  ggplot2::geom_vline(xintercept=0,linetype="dotted",linewidth=.6,col="#777777")+
  ggplot2::geom_point(ggplot2::aes(shape=component,colour=archive),
                      size=3.75,stroke=1.0,fill="white")+
  ggplot2::geom_text(data=MP[MP$asterisk=="*",,drop=FALSE],
       ggplot2::aes(x=star_x,y=y_plot,label=asterisk),inherit.aes=FALSE,
       colour=COL_INK,size=5.5,fontface="bold",hjust=0,vjust=.4)+
  ggplot2::scale_shape_manual(values=shapes,name="Module component",
      labels=c("ALS-downregulated","ALS-upregulated"))+
  ggplot2::scale_colour_manual(values=cols,name="Data repository")+
  yscale(16,registry$label)+
  ggplot2::scale_x_continuous(limits=c(-neslimit,neslimit+.42),
     breaks=pretty(c(-neslimit,neslimit),n=7))+
  ggplot2::labs(title="A   Module-component enrichment",x="Normalised enrichment score (NES)",y=NULL)+
  base_theme+
  ggplot2::guides(shape=ggplot2::guide_legend(order=1,ncol=1,override.aes=list(fill="white")),
                  colour=ggplot2::guide_legend(order=2,ncol=1,
                      override.aes=list(shape=16,size=3.9)))

# ------------------------------- MAIN PANEL B -------------------------------
figB <- ggplot2::ggplot(CP,ggplot2::aes(y=y,x=pct_concordant))+
  horizontal_lines(sep_main)+
  ggplot2::geom_vline(xintercept=50,linetype="dotted",linewidth=.6,col="#777777")+
  ggplot2::geom_segment(ggplot2::aes(x=50,xend=pct_concordant,yend=y,colour=archive),
                        linewidth=.95,alpha=.70)+
  ggplot2::geom_point(ggplot2::aes(fill=archive),shape=23,colour=COL_INK,
                      size=3.9,stroke=.7)+
  ggplot2::geom_text(ggplot2::aes(x=129,label=paste0(n_concordant,"/",n_evaluable)),
                     hjust=1,size=4.15,colour=COL_INK)+
  ggplot2::scale_colour_manual(values=cols,guide="none")+
  ggplot2::scale_fill_manual(values=cols,guide="none")+
  yscale(16)+
  ggplot2::scale_x_continuous(breaks=c(0,25,50,75,100),limits=c(0,135),
                               expand=ggplot2::expansion(mult=0))+
  ggplot2::labs(title="B   Whole-module gene direction",
               x="Genes matching ALS direction (%)",y=NULL)+base_theme+clean_y

# ------------------------------- MAIN PANEL C -------------------------------
score_rng <- range(c(SP$ci_low,SP$ci_high,SP$mean_FLT_minus_control))
score_ext <- max(.30,diff(score_rng)*.11)
SP$star_x <- SP$mean_FLT_minus_control+score_ext*.19
figC <- ggplot2::ggplot(SP,ggplot2::aes(y=y,x=mean_FLT_minus_control))+
  horizontal_lines(sep_main)+
  ggplot2::geom_vline(xintercept=0,linetype="dotted",linewidth=.6,col="#777777")+
  ggplot2::geom_segment(ggplot2::aes(x=ci_low,xend=ci_high,yend=y),
                        linewidth=.76,colour="#676D73")+
  ggplot2::geom_point(ggplot2::aes(fill=archive),shape=22,colour=COL_INK,
                      size=3.8,stroke=.7)+
  ggplot2::geom_text(data=SP[SP$asterisk=="*",,drop=FALSE],
        ggplot2::aes(x=star_x,y=y,label=asterisk),inherit.aes=FALSE,
        size=5.5,fontface="bold",hjust=0,vjust=.4,col=COL_INK)+
  ggplot2::scale_fill_manual(values=cols,guide="none")+
  yscale(16)+
  ggplot2::scale_x_continuous(limits=c(min(score_rng[1]-.15,-.05),
       max(score_rng[2]+score_ext,.05)),breaks=pretty(score_rng,n=6))+
  ggplot2::labs(title="C   Biological-sample module scores",
     x="Exposure minus control (signed z-score)",y=NULL)+base_theme+clean_y

fig8 <- (figA | figB | figC)+
   patchwork::plot_layout(widths=c(1.77,1.09,1.18),guides="collect")+
    ggplot2::theme(legend.position="right",legend.direction="vertical",
                   legend.box="vertical",legend.box.just="left")

# ------------------------- SUPPLEMENTARY FIGURE S2 --------------------------
# MHU-1 onboard MG/AG primary; MG/GC and AG/GC are *secondary*. All
# secondary GSEA FDRs come from the pipeline's separately defined family
# of secondary comparisons (including eligible other missions), not BH3/BH6.
reg_sx <- data.frame(dataset=rep("MHU1",3),tissue=rep("SOLEUS",3),
  cohort=rep("ALL",3),contrast=c("MG_vs_AG","MG_vs_GC","AG_vs_GC"),
  label=c("MG vs onboard 1 g  [primary]",
          "MG vs ground  [secondary]","Onboard 1 g vs ground  [secondary]"),
  archive=rep("JAXA ibSLS",3),stringsAsFactors=FALSE)
reg_sx$key <- keyfun(reg_sx)
MS <- M[M$key %in% reg_sx$key,,drop=FALSE]
GS <- G[G$key %in% reg_sx$key,,drop=FALSE]
SS <- ST[ST$key %in% reg_sx$key,,drop=FALSE]
if (nrow(MS)!=6L || anyDuplicated(paste(MS$key,MS$pathway)) ||
    !setequal(unique(MS$key),reg_sx$key) ||
    !setequal(unique(SS$key),reg_sx$key) || nrow(SS)!=3L)
  stop("Supplementary figure needs all 3 MHU-1 contrasts and 6 GSEA rows",call.=FALSE)
if (any(!MS$primary[MS$contrast=="MG_vs_AG"]) ||
    any(MS$primary[MS$contrast!="MG_vs_AG"]))
  stop("MHU-1 primary/secondary flags are incorrect",call.=FALSE)
if (any(!is.finite(MS$FDR_secondary_across_ALL_spaceflight[
                  MS$contrast!="MG_vs_AG"])))
  stop("Missing globally adjusted secondary comparison FDR for MHU-1",call.=FALSE)
cs <- whole_concordance(GS)
MS <- attach_info(MS,reg_sx);cs <- attach_info(cs,reg_sx)
SS <- attach_info(SS,reg_sx)
MS$component <- factor(MS$pathway,levels=components)
MS$y_plot <- MS$y+ifelse(MS$pathway=="ALS_DOWN",.16,-.16)
MS$plot_FDR <- ifelse(MS$primary,MS$FDR_primary_across_ALL_spaceflight,
                        MS$FDR_secondary_across_ALL_spaceflight)
MS$asterisk <- ifelse(MS$plot_FDR<ALPHA,"*","")
MS$star_x <- MS$NES+.19
# Mark significance for the primary score test only. The secondary score
# comparisons were not BH-adjusted in the source pipeline.
SS$asterisk <- ifelse(SS$primary & is.finite(SS$FDR_primary_scores) &
                        SS$FDR_primary_scores<ALPHA,"*","")

shade_sx <- function() ggplot2::annotate("rect",xmin=-Inf,xmax=Inf,ymin=.5,
                           ymax=2.5,fill="#F3F3F3",alpha=.75)
sx_theme <- base_theme + ggplot2::theme(plot.title=ggplot2::element_text(
                                         face="bold",size=BASE_FONT*1.00))
sx_nes <- max(2.9,ceiling(2*max(abs(MS$NES)))/2+.3)
sxA <- ggplot2::ggplot(MS,ggplot2::aes(x=NES,y=y_plot))+
  shade_sx()+ggplot2::geom_vline(xintercept=0,linetype="dotted",colour="#777777")+
  ggplot2::geom_point(ggplot2::aes(shape=component),size=4,stroke=1,
                      colour=COL_JAXA,fill="white")+
  ggplot2::geom_text(data=MS[MS$asterisk=="*",,drop=FALSE],
     ggplot2::aes(x=star_x,y=y_plot,label=asterisk),inherit.aes=FALSE,
     size=5.6,col=COL_INK,fontface="bold",hjust=0)+
  ggplot2::scale_shape_manual(values=shapes,name="Component",
                        labels=c("ALS down","ALS up"))+
  yscale(3,reg_sx$label)+
  ggplot2::scale_x_continuous(limits=c(-sx_nes,sx_nes+.35))+
  ggplot2::labs(title="A   MHU-1 enrichment",x="Normalised enrichment score",y=NULL)+
  sx_theme
sxB <- ggplot2::ggplot(cs,ggplot2::aes(x=pct_concordant,y=y))+
  shade_sx()+ggplot2::geom_vline(xintercept=50,linetype="dotted",colour="#777777")+
  ggplot2::geom_segment(ggplot2::aes(x=50,xend=pct_concordant,yend=y),
                        linewidth=.9,colour=COL_JAXA)+
  ggplot2::geom_point(size=3.9,shape=23,fill=COL_JAXA,colour=COL_INK)+
  ggplot2::geom_text(ggplot2::aes(x=125,label=paste0(n_concordant,"/",n_evaluable)),
                     size=4.1,hjust=1)+
  yscale(3)+ggplot2::scale_x_continuous(limits=c(0,130),breaks=c(0,25,50,75,100),
       expand=ggplot2::expansion(mult=0))+
  ggplot2::labs(title="B   Whole-module gene direction",x="ALS-direction genes (%)",y=NULL)+
  sx_theme+clean_y
sx_range <- range(c(SS$ci_low,SS$ci_high))
sx_pad <- max(.25,diff(sx_range)*.10)
sxC <- ggplot2::ggplot(SS,ggplot2::aes(x=mean_FLT_minus_control,y=y))+
  shade_sx()+ggplot2::geom_vline(xintercept=0,linetype="dotted",colour="#777777")+
  ggplot2::geom_segment(ggplot2::aes(x=ci_low,xend=ci_high,yend=y),
                        linewidth=.75,colour="#676D73")+
  ggplot2::geom_point(shape=22,size=3.8,fill=COL_JAXA,colour=COL_INK)+
  ggplot2::geom_text(data=SS[SS$asterisk=="*",,drop=FALSE],
        ggplot2::aes(x=mean_FLT_minus_control+sx_pad*.18,y=y,label=asterisk),
        inherit.aes=FALSE,hjust=0,size=5.4,fontface="bold")+
  yscale(3)+ggplot2::scale_x_continuous(limits=c(sx_range[1]-sx_pad,
                                             sx_range[2]+sx_pad))+
  ggplot2::labs(title="C   Module-score contrast",x="Difference (signed z-score)",y=NULL)+
  sx_theme+clean_y

figS <- (sxA | sxB | sxC)+
  patchwork::plot_layout(widths=c(1.55,1.00,1.07),guides="collect")+
  ggplot2::theme(legend.position="right",legend.direction="vertical",
                 legend.box="vertical")

# ----------------------------- SAVE OUTPUTS ---------------------------------
save_pair <- function(p,stem,w,h) {
  png <- file.path(OUTPUT_DIR,paste0(stem,".png"))
  svg <- file.path(OUTPUT_DIR,paste0(stem,".svg"))
  ggplot2::ggsave(png,p,width=w,height=h,units="in",dpi=DPI,bg="white",
                  limitsize=FALSE)
  ggplot2::ggsave(svg,p,width=w,height=h,units="in",bg="white",
                  device=svglite::svglite,limitsize=FALSE)
  message("Saved: ",png,"\nSaved: ",svg)
}
save_pair(fig8,"Figure8_spaceflight_frozen_module_PRIMARY_16",MAIN_W,MAIN_H)
save_pair(figS,"Supplementary_Figure_SX_MHU1_module_controls",SUPP_W,SUPP_H)

# Numeric audit tables: never hide the exact contrasts or FDR families.
mainA <- MP[,c("key","label","archive","pathway","NES","pval",
                "FDR_primary_across_ALL_spaceflight")]
mainB <- CP[,c("key","n_evaluable","n_concordant","pct_concordant")]
mainC <- SP[,c("key","n_exposed","n_control","mean_FLT_minus_control",
                "ci_low","ci_high","pvalue","FDR_primary_scores")]
write.csv(mainA,file.path(OUTPUT_DIR,"Figure8_panelA_primary_GSEA.csv"),row.names=FALSE)
write.csv(mainB,file.path(OUTPUT_DIR,"Figure8_panelB_primary_concordance.csv"),row.names=FALSE)
write.csv(mainC,file.path(OUTPUT_DIR,"Figure8_panelC_primary_sample_scores.csv"),row.names=FALSE)
write.csv(MS[,c("key","label","pathway","primary","NES","pval","plot_FDR")],
   file.path(OUTPUT_DIR,"Supplementary_S2_panelA_MHU1_GSEA.csv"),row.names=FALSE)
write.csv(cs[,c("key","n_evaluable","n_concordant","pct_concordant")],
   file.path(OUTPUT_DIR,"Supplementary_S2_panelB_MHU1_concordance.csv"),row.names=FALSE)
write.csv(SS[,c("key","primary","mean_FLT_minus_control","ci_low",
                 "ci_high","pvalue","FDR_primary_scores")],
   file.path(OUTPUT_DIR,"Supplementary_S2_panelC_MHU1_sample_scores.csv"),row.names=FALSE)
writeLines(c(
  "Figure 8 and Supplementary Figure S2: provenance and interpretation",
  "Figure 8: exactly 16 primary comparisons; 32 GSEA tests, 16 sample-score tests.",
  "Panel A: circles=ALS_DOWN, triangles=ALS_UP; NES and global BH-adjusted FDR from v9.",
  "Panel B: diamonds=whole-module concordance; fraction of genes with finite log2FC",
  "and sign agreeing with fixed ALS sign; no fold-change or per-gene FDR threshold.",
  "Panel C: squares=mean score differences; Welch 95% CIs; unadjusted CI",
  "coverage, significance asterisk by BH over 16 primary score tests.",
  "Supplementary SX: MHU-1 MG vs AG primary; MG vs GC and AG vs GC secondary.",
  "Secondary GSEA FDR is separate from primary FDR and taken from complete",
  "secondary comparison test family in the unified pipeline.",
  "Secondary sample score p-values are NOT BH-adjusted; no secondary stars.",
  "The graph does not itself establish the cause of an atypical GC sample.",
  "For MHU-1 ground-control QC/leave-one-out sensitivity, use separate S1 QC.",
  paste("Inputs:",paste(unlist(paths),collapse=" | ")),
  paste("Generated:",format(Sys.time(),tz="UTC",usetz=TRUE))),
  file.path(OUTPUT_DIR,"FIGURE8_S2_PROVENANCE.txt"))
message("DONE: outputs in ",normalizePath(OUTPUT_DIR,winslash="/"))
