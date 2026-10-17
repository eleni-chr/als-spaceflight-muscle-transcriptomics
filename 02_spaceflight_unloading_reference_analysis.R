# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

###############################################################################
# NASA / HUMAN UNLOADING / ALS SKELETAL-MUSCLE TRANSCRIPTOMICS
# FULL RE-RUN SCRIPT
#
# Put this script and ALL input files in the SAME folder and run from the top.
#
# NASA RNA-seq primary rule:
#   ONLY per-sample RSEM files containing rRNArm and ending .genes.results
#   (or .genes.results.gz) are used for OSD-101, OSD-419, OSD-684 and OSD-540.
#
# OSD-21: normalized microarray probeset matrix -> limma
# GSE271607: supplied raw Read_counts -> paired DESeq2, VL/soleus separately
#
# ALS reference:
#   SOD_data.xlsx Tables S7-S12: directional SOD1G93A mRNA KEGG reference
###############################################################################

input_dir <- "."
output_dir <- file.path(input_dir,"NASA_ALS_metabolic_pathway_results_FIXED_CANONICAL")
FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 1
MIN_COUNT <- 10
MIN_SAMPLES <- 3
set.seed(20260930)

cran_pkgs <- c("tidyverse","openxlsx","readxl","ggrepel","svglite","msigdbr")
bioc_pkgs <- c("DESeq2","tximport","limma","AnnotationDbi","org.Mm.eg.db","org.Hs.eg.db","fgsea")

install_if_missing <- function(pkgs, bioc=FALSE){
  m <- pkgs[!vapply(pkgs,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))]
  if(!length(m)) return(invisible(NULL))
  if(bioc){
    if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
    BiocManager::install(m,ask=FALSE,update=FALSE)
  } else install.packages(m)
}
install_if_missing(cran_pkgs)
install_if_missing(bioc_pkgs,TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(openxlsx); library(readxl); library(ggrepel)
  library(svglite); library(msigdbr); library(DESeq2); library(tximport)
  library(limma); library(AnnotationDbi); library(org.Mm.eg.db)
  library(org.Hs.eg.db); library(fgsea)
})

dir.create(output_dir,showWarnings=FALSE,recursive=TRUE)
for(d in c("DE","GSEA","ALS_reference","robustness","figures","QC","summary"))
  dir.create(file.path(output_dir,d),showWarnings=FALSE,recursive=TRUE)

ALS_CORE_MOUSE <- c("Pgm1","Tpi1","Pgk1","Eno3","Ldha","Fbp2","Prkaa2")
ALS_CORE_HUMAN <- toupper(ALS_CORE_MOUSE)

# Primary ALS-derived pathway registry.
#
# IMPORTANT: GSEA extraction below uses EXACT MSigDB set names only. There is
# deliberately no fuzzy pathway matching. Three SOD1 pathways do not have an
# exact full-pathway counterpart in the MSigDB collections used here and are
# therefore retained in the reference table but marked as not formally tested.
SOD1_REFERENCE_PATHWAYS <- tribble(
  ~reference_pathway, ~sod1_terms, ~msigdb_set, ~formal_test, ~note,
  "Glycolysis / gluconeogenesis",
    "GLYCOLYSIS_GLUCONEOGENESIS",
    "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
    TRUE, "Exact KEGG pathway",
  "Fructose and mannose metabolism",
    "FRUCTOSE_AND_MANNOSE_METABOLISM",
    "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",
    TRUE, "Exact KEGG pathway",
  "Pentose phosphate pathway",
    "PENTOSE_PHOSPHATE_PATHWAY",
    "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
    TRUE, "Exact KEGG pathway",
  "Carbon metabolism",
    "CARBON_METABOLISM",
    NA_character_,
    FALSE, "No exact full KEGG carbon-metabolism set in the loaded MSigDB collections; not substituted",
  "AMPK signalling",
    "AMPK_SIGNALING_PATHWAY;AMPK_SIGNALLING_PATHWAY",
    NA_character_,
    FALSE, "No exact full KEGG AMPK-signalling set in the loaded MSigDB collections; not substituted",
  "Insulin signalling",
    "INSULIN_SIGNALING_PATHWAY;INSULIN_SIGNALLING_PATHWAY",
    "KEGG_INSULIN_SIGNALING_PATHWAY",
    TRUE, "Exact KEGG pathway",
  "Glucagon signalling",
    "GLUCAGON_SIGNALING_PATHWAY;GLUCAGON_SIGNALLING_PATHWAY",
    NA_character_,
    FALSE, "No exact full KEGG glucagon-signalling set in the loaded MSigDB collections; not substituted",
  "Purine metabolism",
    "PURINE_METABOLISM",
    "KEGG_PURINE_METABOLISM",
    TRUE, "Exact KEGG pathway",
  "Starch and sucrose metabolism",
    "STARCH_AND_SUCROSE_METABOLISM",
    "KEGG_STARCH_AND_SUCROSE_METABOLISM",
    TRUE, "Exact KEGG pathway"
)

# Separate, exact secondary bioenergetic panel used to describe the broader
# unloading phenotype. These are NOT labelled as SOD1-specific pathway tests.
BIOENERGETIC_EXACT_PATHWAYS <- tribble(
  ~bioenergetic_pathway, ~msigdb_set,
  "Glycolysis", "HALLMARK_GLYCOLYSIS",
  "Oxidative phosphorylation", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "Fatty acid metabolism", "HALLMARK_FATTY_ACID_METABOLISM"
)

clean_ens <- function(x){
  x <- as.character(x); x <- sub("\\|.*$","",x); sub("\\..*$","",x)
}
norm_path <- function(x){
  x %>% toupper() %>%
    str_replace_all("^HALLMARK_|^KEGG_LEGACY_|^KEGG_MEDICUS_|^KEGG_|^REACTOME_","") %>%
    str_replace_all("[^A-Z0-9]+","_") %>% str_replace_all("_+","_") %>%
    str_replace_all("^_|_$","")
}
find_one <- function(pattern,required=TRUE){
  z <- list.files(input_dir,full.names=TRUE)
  h <- z[str_detect(basename(z),regex(pattern,ignore_case=TRUE))]
  if(!length(h)){
    if(required) stop("Could not find file matching: ",pattern)
    return(NA_character_)
  }
  if(length(h)>1) message("Multiple matches for ",pattern,"; using ",basename(h[1]))
  h[1]
}
map_ens <- function(ids,species=c("mouse","human")){
  species <- match.arg(species); ids <- unique(clean_ens(ids))
  db <- if(species=="mouse") org.Mm.eg.db else org.Hs.eg.db
  sy <- AnnotationDbi::mapIds(db,keys=ids,keytype="ENSEMBL",column="SYMBOL",multiVals="first")
  tibble(ensembl_id=ids,symbol=unname(sy[ids]))
}
collapse_de <- function(x){
  x %>% filter(!is.na(symbol),symbol!="",is.finite(statistic)) %>%
    group_by(symbol) %>% slice_max(abs(statistic),n=1,with_ties=FALSE) %>% ungroup()
}
collapse_expr <- function(m,map){
  as.data.frame(m) %>% rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id=clean_ens(ensembl_id)) %>% left_join(map,by="ensembl_id") %>%
    filter(!is.na(symbol),symbol!="") %>% dplyr::select(-ensembl_id) %>%
    group_by(symbol) %>% summarise(across(where(is.numeric),~mean(.x,na.rm=TRUE)),.groups="drop") %>%
    column_to_rownames("symbol") %>% as.matrix()
}
sample_id_rsem <- function(x){
  b <- basename(x); b <- str_remove(b,"\\.gz$"); b <- str_remove(b,"\\.genes\\.results$")
  str_remove(b,"_rRNArm.*$")
}
save_pca <- function(expr,meta,title,file){
  x <- expr[,meta$sample,drop=FALSE]; v <- apply(x,1,var,na.rm=TRUE); v <- v[is.finite(v)]
  k <- names(sort(v,decreasing=TRUE))[seq_len(min(1000,length(v)))]
  pc <- prcomp(t(x[k,,drop=FALSE]),center=TRUE,scale.=FALSE)
  pct <- 100*pc$sdev^2/sum(pc$sdev^2)
  d <- as_tibble(pc$x[,1:2,drop=FALSE],rownames="sample") %>% left_join(meta,by="sample")
  g <- ggplot(d,aes(PC1,PC2,colour=group,label=sample))+geom_point(size=3)+
    ggrepel::geom_text_repel(size=2.2,max.overlaps=Inf)+theme_bw()+
    labs(title=title,x=paste0("PC1 (",round(pct[1],1),"%)"),
         y=paste0("PC2 (",round(pct[2],1),"%)"),colour="Group")
  ggsave(file,g,width=8,height=6,device="svg")
}
save_volcano <- function(de,genes,title,file){
  d <- de %>% mutate(nlf=-log10(pmax(FDR,.Machine$double.xmin)),
    status=case_when(FDR<FDR_CUTOFF&log2FC>=LFC_CUTOFF~"Up",
                     FDR<FDR_CUTOFF&log2FC<=-LFC_CUTOFF~"Down",TRUE~"NS"),
    core=toupper(symbol)%in%toupper(genes))
  g <- ggplot(d,aes(log2FC,nlf))+geom_point(aes(shape=status),alpha=.4,size=1.2)+
    geom_vline(xintercept=c(-1,1),linetype=2)+geom_hline(yintercept=-log10(.05),linetype=2)+
    geom_point(data=filter(d,core),size=2.5)+
    ggrepel::geom_text_repel(data=filter(d,core),aes(label=symbol),size=2.8,max.overlaps=Inf)+
    theme_bw()+labs(title=title,x="log2 fold change",y="-log10(FDR)")
  ggsave(file,g,width=8,height=6,device="svg")
}
save_scatter <- function(a,b,xlab,ylab,title,file){
  d <- inner_join(dplyr::select(a,symbol,x=log2FC),dplyr::select(b,symbol,y=log2FC),by="symbol") %>%
    filter(is.finite(x),is.finite(y))
  rho <- cor(d$x,d$y,method="spearman",use="complete.obs")
  ct <- suppressWarnings(cor.test(d$x,d$y,method="spearman",exact=FALSE))
  g <- ggplot(d,aes(x,y))+geom_hline(yintercept=0,linetype=2)+geom_vline(xintercept=0,linetype=2)+
    geom_point(alpha=.25,size=1.1)+theme_bw()+labs(title=title,
    subtitle=paste0("Spearman rho=",round(rho,3),"; P=",signif(ct$p.value,3)),x=xlab,y=ylab)
  ggsave(file,g,width=7,height=6.5,device="svg")
  tibble(comparison=title,n_genes=nrow(d),spearman_rho=rho,p_value=ct$p.value)
}

safe_msig <- function(sp,c,sc=NULL){
  tryCatch(if(is.null(sc)) msigdbr(species=sp,collection=c)
           else msigdbr(species=sp,collection=c,subcollection=sc),
           error=function(e){message("Skipping MSigDB ",c,"/",sc,": ",e$message); tibble()})
}
get_pathways <- function(sp){
  bind_rows(safe_msig(sp,"H"),safe_msig(sp,"C2","CP:KEGG_LEGACY"),
            safe_msig(sp,"C2","CP:KEGG_MEDICUS"),safe_msig(sp,"C2","CP:REACTOME")) %>%
    dplyr::select(gs_name,gene_symbol) %>% filter(!is.na(gene_symbol),gene_symbol!="") %>%
    distinct() %>% {split(.$gene_symbol,.$gs_name)}
}
run_gsea <- function(de,pathways,dataset,contrast,variant="primary"){
  r <- de %>% filter(!is.na(symbol),symbol!="",is.finite(statistic)) %>%
    group_by(symbol) %>% slice_max(abs(statistic),n=1,with_ties=FALSE) %>% ungroup()
  s <- r$statistic; names(s) <- r$symbol; s <- sort(s[is.finite(s)],decreasing=TRUE)
  fgseaMultilevel(pathways=pathways,stats=s,minSize=5,maxSize=500,eps=0) %>%
    as_tibble() %>% arrange(padj,desc(abs(NES))) %>%
    mutate(dataset=dataset,contrast=contrast,variant=variant,pathway_norm=norm_path(pathway),
           leadingEdge=vapply(leadingEdge,paste,collapse=";",FUN.VALUE=character(1))) %>%
    relocate(dataset,contrast,variant)
}
filter_metabolic <- function(x){
  kw <- "GLYCOL|GLUCOSE|CARBOHYDRATE|PENTOSE|FRUCTOSE|MANNOSE|PYRUVATE|AMPK|OXIDATIVE_PHOSPHORYLATION|MITOCHONDR|FATTY_ACID|LIPID_METABOL|INSULIN|GLUCAGON|TCA|CITRIC_ACID|CARBON_METABOL|PURINE|STARCH|SUCROSE"
  filter(x,str_detect(pathway,regex(kw,ignore_case=TRUE)))
}
extract_exact_gsea_sets <- function(x, spec, label_col){
  map_dfr(seq_len(nrow(spec)), function(i){
    ref <- spec[i,]
    label <- as.character(ref[[label_col]][1])
    set_name <- as.character(ref$msigdb_set[1])

    if(is.na(set_name) || set_name == ""){
      return(tibble(
        dataset=unique(x$dataset)[1], contrast=unique(x$contrast)[1],
        variant=unique(x$variant)[1], label=label, msigdb_set=NA_character_,
        test_status="Not tested: no exact canonical MSigDB set specified",
        NES=NA_real_, pval=NA_real_, padj=NA_real_, size=NA_integer_
      ))
    }

    h <- x %>% dplyr::filter(pathway == set_name)

    if(nrow(h) == 0){
      return(tibble(
        dataset=unique(x$dataset)[1], contrast=unique(x$contrast)[1],
        variant=unique(x$variant)[1], label=label, msigdb_set=set_name,
        test_status="Not tested: exact set absent from this GSEA result",
        NES=NA_real_, pval=NA_real_, padj=NA_real_, size=NA_integer_
      ))
    }

    h <- h %>% dplyr::slice(1)
    tibble(
      dataset=h$dataset, contrast=h$contrast, variant=h$variant,
      label=label, msigdb_set=set_name, test_status="Tested",
      NES=h$NES, pval=h$pval, padj=h$padj, size=h$size
    )
  })
}

core_from_gsea <- function(x){
  extract_exact_gsea_sets(x, SOD1_REFERENCE_PATHWAYS, "reference_pathway") %>%
    dplyr::rename(reference_pathway=label) %>%
    dplyr::left_join(
      SOD1_REFERENCE_PATHWAYS %>%
        dplyr::select(reference_pathway, formal_test, note),
      by="reference_pathway"
    )
}

bioenergetic_from_gsea <- function(x){
  extract_exact_gsea_sets(x, BIOENERGETIC_EXACT_PATHWAYS, "bioenergetic_pathway") %>%
    dplyr::rename(bioenergetic_pathway=label)
}



# ----- DESeq2 / RSEM helpers --------------------------------------------------
find_rnarm_rsem <- function(){
  z <- list.files(input_dir,full.names=TRUE)
  allr <- z[str_detect(basename(z),regex("\\.genes\\.results(\\.gz)?$",ignore_case=TRUE))]
  use <- allr[str_detect(basename(allr),regex("rRNArm.*\\.genes\\.results(\\.gz)?$",ignore_case=TRUE))]
  write_csv(tibble(file=basename(allr),used_primary=allr%in%use,
                   reason=if_else(allr%in%use,"rRNArm RSEM","Ignored: non-rRNArm")),
            file.path(output_dir,"QC","RSEM_file_selection_audit.csv"))
  use
}
classify_rsem <- function(files){
  b <- basename(files)
  ds <- case_when(
    str_detect(b,regex("GLDS[-_]?101|OSD[-_]?101",ignore_case=TRUE))~"OSD101",
    str_detect(b,regex("GLDS[-_]?419|OSD[-_]?419",ignore_case=TRUE))~"OSD419",
    str_detect(b,regex("GLDS[-_]?615|OSD[-_]?684|GLDS[-_]?684",ignore_case=TRUE))~"OSD684",
    str_detect(b,regex("OSD[-_]?540|GLDS[-_]?540|GSE211204",ignore_case=TRUE))~"OSD540",
    TRUE~NA_character_)
  gsm <- suppressWarnings(as.integer(str_remove(str_extract(b,"GSM\\d+"),"GSM")))
  ds[is.na(ds)&gsm>=6456083&gsm<=6456113] <- "OSD540"
  tibble(file=files,dataset=ds)
}
import_rsem <- function(files,meta,design){
  names(files) <- rownames(meta)
  txi <- tximport(files,type="rsem",txIn=FALSE,txOut=FALSE)
  stopifnot(identical(colnames(txi$counts),rownames(meta)))

  # DESeq2 requires every effective length used in the tximport offset to be
  # finite and strictly >0. Some RSEM gene rows have effective_length=0.
  valid_length <- apply(txi$length,1,function(z) all(is.finite(z) & z>0))
  message("RSEM: removing ",sum(!valid_length),
          " genes with zero/non-finite effective length; ",
          sum(valid_length)," retained before count filtering.")
  if(!any(valid_length)) stop("No RSEM genes remain after effective-length validation.")

  txi$counts <- txi$counts[valid_length,,drop=FALSE]
  txi$abundance <- txi$abundance[valid_length,,drop=FALSE]
  txi$length <- txi$length[valid_length,,drop=FALSE]

  dds <- DESeqDataSetFromTximport(txi,colData=meta,design=design)
  keep <- rowSums(counts(dds)>=MIN_COUNT,na.rm=TRUE)>=MIN_SAMPLES
  dds <- dds[keep,]
  message("RSEM: ",nrow(dds)," genes retained after expression filtering.")
  dds
}
res_table <- function(rr,map,dataset,contrast){
  t <- as.data.frame(rr) %>% rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id=clean_ens(ensembl_id),statistic=stat) %>% left_join(map,by="ensembl_id") %>%
    transmute(ensembl_id,symbol,baseMean,log2FC=log2FoldChange,lfcSE,statistic,
              p_value=pvalue,FDR=padj,dataset=dataset,contrast=contrast,
              significant_FDR=!is.na(FDR)&FDR<FDR_CUTOFF,
              significant_FDR_LFC=!is.na(FDR)&FDR<FDR_CUTOFF&abs(log2FC)>=LFC_CUTOFF)
  list(ensembl=t,symbol=collapse_de(t))
}
fit_deseq <- function(dds,map,dataset,contrasts){
  dds <- DESeq(dds); vsd <- vst(dds,blind=FALSE); ev <- assay(vsd); rownames(ev) <- clean_ens(rownames(ev))
  out <- list()
  for(nm in names(contrasts)){
    sp <- contrasts[[nm]]
    rr <- if(is.character(sp)) results(dds,contrast=sp,alpha=FDR_CUTOFF)
          else results(dds,contrast=sp,alpha=FDR_CUTOFF)
    out[[nm]] <- res_table(rr,map,dataset,nm)
  }
  list(results=out,dds=dds,vsd=vsd,expr_vst=ev)
}
save_de <- function(fit,prefix){
  for(nm in names(fit$results)){
    r <- fit$results[[nm]]
    if(!is.null(r$ensembl)) write_csv(r$ensembl,file.path(output_dir,"DE",paste0(prefix,"_",nm,"_Ensembl_level.csv")))
    write_csv(r$symbol,file.path(output_dir,"DE",paste0(prefix,"_",nm,"_Gene_symbol_level.csv")))
  }
}

# ----- OSD-21 microarray helpers ---------------------------------------------
prep21 <- function(df){
  sc <- grep("^Mmus_",colnames(df),value=TRUE)
  df %>% filter(!is.na(SYMBOL),SYMBOL!="") %>% dplyr::select(SYMBOL,all_of(sc)) %>%
    group_by(SYMBOL) %>% summarise(across(all_of(sc),~median(.x,na.rm=TRUE)),.groups="drop") %>%
    column_to_rownames("SYMBOL") %>% as.matrix()
}
fit_limma <- function(expr,meta,contrasts,dataset){
  e <- expr[,meta$sample,drop=FALSE]
  design <- model.matrix(~0+group,data=meta); colnames(design) <- levels(meta$group); rownames(design)<-meta$sample
  cm <- makeContrasts(contrasts=unname(contrasts),levels=design); colnames(cm)<-names(contrasts)
  f <- eBayes(contrasts.fit(lmFit(e,design),cm),robust=TRUE)
  out <- list()
  for(nm in colnames(cm)){
    t <- topTable(f,coef=nm,number=Inf,sort.by="P") %>% rownames_to_column("symbol") %>%
      transmute(symbol,log2FC=logFC,statistic=t,p_value=P.Value,FDR=adj.P.Val,B,
                dataset=dataset,contrast=nm,
                significant_FDR=FDR<FDR_CUTOFF,
                significant_FDR_LFC=FDR<FDR_CUTOFF&abs(log2FC)>=LFC_CUTOFF)
    out[[nm]] <- list(ensembl=NULL,symbol=t)
  }
  list(results=out,expr_vst=e,fit=f)
}

# ----- Flexible Excel enrichment-table reader --------------------------------
find_sheet <- function(path,n){
  sh <- excel_sheets(path)
  hit <- sh[str_detect(sh,regex(paste0("(^|[^0-9])S?",n,"([^0-9]|$)"),ignore_case=TRUE))]
  if(!length(hit)) stop("Cannot find Table S",n," in ",basename(path),". Sheets: ",paste(sh,collapse=" | "))
  hit[1]
}
read_sheet_auto <- function(path,sheet){
  pre <- suppressMessages(read_excel(path,sheet=sheet,col_names=FALSE,n_max=25))
  if(!nrow(pre)) return(tibble())
  score <- apply(pre,1,function(z){
    s <- tolower(paste(unlist(z),collapse=" "))
    sum(str_detect(s,c("pathway","term","description","p.value","p-value","fdr","benjamini","genes","count")))
  })
  h <- if(max(score,na.rm=TRUE)>0) which.max(score) else 1
  suppressMessages(read_excel(path,sheet=sheet,skip=h-1,.name_repair="unique"))
}
find_col <- function(df,pats){
  nm <- colnames(df)
  sc <- map_int(nm,~sum(map_lgl(pats,function(p) str_detect(.x,regex(p,ignore_case=TRUE)))))
  if(!length(sc)||max(sc)==0) return(NA_character_)
  nm[which.max(sc)]
}
extract_enrich <- function(path,n,source,age,direction){
  sh <- find_sheet(path,n); x <- read_sheet_auto(path,sh); if(!nrow(x)) return(tibble())
  tc <- find_col(x,c("pathway","term","description","name"))
  fc <- find_col(x,c("adjust","adj","fdr","benjamini","q.value","q-value"))
  pc <- find_col(x,c("^p$","p.value","p-value","p value"))
  if(is.na(tc)){
    cc <- names(x)[vapply(x,function(z)is.character(z)||is.factor(z),logical(1))]
    tc <- if(length(cc)) cc[1] else colnames(x)[1]
  }
  tibble(source=source,table=paste0("S",n),sheet=sh,age=age,direction=direction,
         term=as.character(x[[tc]]),
         p_value=if(!is.na(pc)) suppressWarnings(as.numeric(x[[pc]])) else NA_real_,
         FDR=if(!is.na(fc)) suppressWarnings(as.numeric(x[[fc]])) else NA_real_) %>%
    filter(!is.na(term),term!="") %>% mutate(term_norm=norm_path(term))
}

# ----- Locate inputs ----------------------------------------------------------
message("Locating files...")
f21 <- find_one("GLDS[-_]?21.*array.*normalized.*probeset.*\\.csv$")
f_sod <- find_one("^SOD_data\\.xlsx$")
f271 <- find_one("GSE271607.*Read_counts.*\\.txt(\\.gz)?$")
f419_star <- find_one("GLDS[-_]?419.*STAR.*(count|unnormal).*\\.csv$",FALSE)

rall <- find_rnarm_rsem()
rc <- classify_rsem(rall)
files101 <- filter(rc,dataset=="OSD101")$file
files419 <- filter(rc,dataset=="OSD419")$file
files684 <- filter(rc,dataset=="OSD684")$file
files540 <- filter(rc,dataset=="OSD540")$file
if(!length(files540)){
  u <- filter(rc,is.na(dataset))$file
  if(length(u)==31){ message("Using 31 unassigned rRNArm files as OSD-540."); files540 <- u }
}
if(!length(files101)) stop("No OSD-101 rRNArm .genes.results detected.")
if(!length(files419)) stop("No OSD-419 rRNArm .genes.results detected.")
if(!length(files684)) stop("No OSD-684/GLDS-615 rRNArm .genes.results detected.")
if(length(files540)!=31) stop("Expected 31 OSD-540 rRNArm .genes.results files; found ",length(files540))
write_csv(mutate(rc,file=basename(file)),file.path(output_dir,"QC","RSEM_dataset_classification.csv"))

# ----- Parse ALS reference ----------------------------------------------------
message("Parsing SOD1 reference tables...")
sod <- bind_rows(
  extract_enrich(f_sod,7,"SOD1G93A","90d","Up"),
  extract_enrich(f_sod,8,"SOD1G93A","90d","Down"),
  extract_enrich(f_sod,9,"SOD1G93A","120d","Up"),
  extract_enrich(f_sod,10,"SOD1G93A","120d","Down"),
  extract_enrich(f_sod,11,"SOD1G93A","150d","Up"),
  extract_enrich(f_sod,12,"SOD1G93A","150d","Down")
)
write_csv(sod,file.path(output_dir,"ALS_reference","SOD1_mRNA_KEGG_tables_S7_S12_parsed.csv"))
sod_down <- filter(sod,direction=="Down")
sod_ref <- map_dfr(seq_len(nrow(SOD1_REFERENCE_PATHWAYS)),function(i){
  r <- SOD1_REFERENCE_PATHWAYS[i,]
  accepted <- stringr::str_split(as.character(r$sod1_terms),";",simplify=FALSE)[[1]] %>% stringr::str_trim()
  h <- sod_down %>% dplyr::filter(term_norm %in% accepted)
  tibble(
    reference_pathway=as.character(r$reference_pathway),
    detected_in_SOD1_down_tables=nrow(h)>0,
    ages=if(nrow(h)) paste(sort(unique(h$age)),collapse=";") else NA_character_,
    n_ages=dplyr::n_distinct(h$age),
    matched_SOD1_terms=if(nrow(h)) paste(unique(h$term),collapse=" | ") else NA_character_,
    exact_msigdb_set=as.character(r$msigdb_set),
    formal_GSEA_test=as.logical(r$formal_test),
    test_note=as.character(r$note)
  )
})
write_csv(sod_ref,file.path(output_dir,"ALS_reference","SOD1_mRNA_metabolic_reference_summary.csv"))

# ----- OSD-101 ---------------------------------------------------------------
message("OSD-101...")
meta101 <- tibble(file=files101,sample=sample_id_rsem(files101),dataset="OSD-101",
  group=case_when(str_detect(basename(file),regex("_FLT_|-FLT_",ignore_case=TRUE))~"Flight",
                  str_detect(basename(file),regex("_GC_|-GC_",ignore_case=TRUE))~"Ground",
                  TRUE~NA_character_))
if(any(is.na(meta101$group))){ print(meta101); stop("Cannot parse OSD-101 groups.") }
meta101 <- mutate(meta101,group=factor(group,levels=c("Ground","Flight")))
cd <- meta101 %>% dplyr::select(sample,group) %>% as.data.frame(); rownames(cd)<-cd$sample; cd$sample<-NULL
dds101 <- import_rsem(meta101$file,cd,~group); map101 <- map_ens(rownames(dds101),"mouse")
fit101 <- fit_deseq(dds101,map101,"OSD-101",list(Flight_vs_Ground=c("group","Flight","Ground")))
save_de(fit101,"OSD101")
save_pca(fit101$expr_vst,meta101,"OSD-101: mouse gastrocnemius",
         file.path(output_dir,"figures","OSD101_PCA.svg"))

# ----- OSD-419 ---------------------------------------------------------------
message("OSD-419...")
meta419 <- tibble(file=files419,sample=sample_id_rsem(files419),dataset="OSD-419",
  group=case_when(str_detect(basename(file),regex("_FLT_|-FLT_",ignore_case=TRUE))~"Flight",
                  str_detect(basename(file),regex("_GC_|-GC_",ignore_case=TRUE))~"Ground",
                  TRUE~NA_character_))
if(any(is.na(meta419$group))){ print(meta419); stop("Cannot parse OSD-419 groups.") }
meta419 <- mutate(meta419,group=factor(group,levels=c("Ground","Flight")))
cd <- meta419 %>% dplyr::select(sample,group) %>% as.data.frame(); rownames(cd)<-cd$sample; cd$sample<-NULL
dds419 <- import_rsem(meta419$file,cd,~group); map419 <- map_ens(rownames(dds419),"mouse")
fit419 <- fit_deseq(dds419,map419,"OSD-419",list(Flight_vs_Ground=c("group","Flight","Ground")))
save_de(fit419,"OSD419")
save_pca(fit419$expr_vst,meta419,"OSD-419: mouse gastrocnemius",
         file.path(output_dir,"figures","OSD419_PCA.svg"))

# ----- OSD-21 ---------------------------------------------------------------
message("OSD-21...")
x21 <- read_csv(f21,show_col_types=FALSE); expr21 <- prep21(x21); s21 <- colnames(expr21)
meta21 <- tibble(sample=s21,dataset="OSD-21",
 group=case_when(str_detect(s21,"_HLU-RL_")~"HLU_RL",str_detect(s21,"_HLU_")~"HLU",
                 str_detect(s21,"_CTRL_")~"CTRL",str_detect(s21,"_GC_")~"GC",
                 str_detect(s21,"_FLT_")~"FLT",TRUE~NA_character_)) %>%
 mutate(group=factor(group,levels=c("CTRL","HLU","HLU_RL","GC","FLT")))
if(any(is.na(meta21$group))) stop("Cannot parse all OSD-21 groups.")
fit21 <- fit_limma(expr21,meta21,c(Flight_vs_Ground="FLT-GC",HLU_vs_CTRL="HLU-CTRL",
                                   HLU_RL_vs_HLU="HLU_RL-HLU",HLU_RL_vs_CTRL="HLU_RL-CTRL"),"OSD-21")
save_de(fit21,"OSD21")
save_pca(expr21,meta21,"OSD-21: mouse gastrocnemius",file.path(output_dir,"figures","OSD21_PCA.svg"))

# ----- OSD-684 / GLDS-615 ----------------------------------------------------
message("OSD-684...")
meta684 <- tibble(file=files684,sample=sample_id_rsem(files684),
                  gsm=str_extract(basename(files684),"GSM\\d+"),dataset="OSD-684") %>%
 mutate(gsm_num=suppressWarnings(as.integer(str_remove(gsm,"GSM"))),
 group=case_when(
   gsm_num%in%7469191:7469193~"Young_Ground",gsm_num%in%7469194:7469196~"Old_Ground",
   gsm_num%in%7469197:7469199~"Young_Flight",gsm_num%in%7469200:7469202~"Old_Flight",
   str_detect(basename(file),regex("young.*ground|ground.*young",ignore_case=TRUE))~"Young_Ground",
   str_detect(basename(file),regex("old.*ground|older.*ground|ground.*old",ignore_case=TRUE))~"Old_Ground",
   str_detect(basename(file),regex("young.*flight|flight.*young",ignore_case=TRUE))~"Young_Flight",
   str_detect(basename(file),regex("old.*flight|older.*flight|flight.*old",ignore_case=TRUE))~"Old_Flight",
   TRUE~NA_character_),
 group=factor(group,levels=c("Young_Ground","Young_Flight","Old_Ground","Old_Flight")))
if(any(is.na(meta684$group))){ print(meta684); stop("Cannot parse all OSD-684 groups.") }
cd <- meta684 %>% dplyr::select(sample,group) %>% as.data.frame(); rownames(cd)<-cd$sample; cd$sample<-NULL
dds684 <- import_rsem(meta684$file,cd,~group); map684 <- map_ens(rownames(dds684),"human")
dds684 <- DESeq(dds684); rn <- resultsNames(dds684)
yf <- rn[str_detect(rn,regex("Young_Flight.*Young_Ground",ignore_case=TRUE))][1]
og <- rn[str_detect(rn,regex("Old_Ground.*Young_Ground",ignore_case=TRUE))][1]
of <- rn[str_detect(rn,regex("Old_Flight.*Young_Ground",ignore_case=TRUE))][1]
if(any(is.na(c(yf,og,of)))) stop("Cannot construct OSD-684 composite contrasts. Coefficients: ",paste(rn,collapse=", "))
mk684 <- function(spec,nm){
  rr <- if(is.character(spec)) results(dds684,contrast=spec,alpha=FDR_CUTOFF) else {
    v <- rep(0,length(rn)); names(v)<-rn; v[names(spec)]<-spec
    results(dds684,contrast=unname(v),alpha=FDR_CUTOFF)
  }
  res_table(rr,map684,"OSD-684",nm)
}
fit684 <- list(results=list(
 Young_Flight_vs_Ground=mk684(c("group","Young_Flight","Young_Ground"),"Young_Flight_vs_Ground"),
 Old_Flight_vs_Ground=mk684(c("group","Old_Flight","Old_Ground"),"Old_Flight_vs_Ground"),
 Average_Flight_vs_Ground=mk684(setNames(c(.5,-.5,.5),c(yf,og,of)),"Average_Flight_vs_Ground"),
 Age_by_Flight_Interaction=mk684(setNames(c(-1,-1,1),c(yf,og,of)),"Age_by_Flight_Interaction")
),dds=dds684)
fit684$vsd <- vst(dds684,blind=FALSE); fit684$expr_vst<-assay(fit684$vsd); rownames(fit684$expr_vst)<-clean_ens(rownames(fit684$expr_vst))
save_de(fit684,"OSD684")
save_pca(fit684$expr_vst,meta684,"OSD-684: human engineered skeletal muscle",
         file.path(output_dir,"figures","OSD684_PCA.svg"))

# ----- OSD-540 / GSE211204 ---------------------------------------------------
message("OSD-540 human vastus lateralis...")
geo540 <- tibble(
 gsm=c(paste0("GSM",6456083:6456093),paste0("GSM",6456094:6456104),paste0("GSM",6456105:6456113)),
 subject=c(as.character(1:11),as.character(1:11),as.character(c(1,2,3,6,7,8,9,10,11))),
 condition=c(rep("Baseline",11),rep("ULLS",11),rep("AR",9)),idx=as.character(1:31))
meta540 <- tibble(file=files540,raw_sample=sample_id_rsem(files540),
                  gsm=str_extract(basename(files540),"GSM\\d+"),
                  idx=str_match(basename(files540),"(?:^|[_\\-])S(\\d+)(?:[_\\.\\-]|$)")[,2]) %>%
 left_join(dplyr::select(geo540,gsm,subject,condition),by="gsm")
if(any(is.na(meta540$subject))){
  fb <- tibble(idx=meta540$idx) %>% left_join(dplyr::select(geo540,idx,sub2=subject,cond2=condition),by="idx")
  meta540 <- mutate(meta540,subject=coalesce(subject,fb$sub2),condition=coalesce(condition,fb$cond2))
}
if(any(is.na(meta540$subject))||any(is.na(meta540$condition))){ print(meta540); stop("Cannot map all OSD-540 files.") }
meta540 <- meta540 %>% mutate(sample=if_else(!is.na(gsm),gsm,paste0("OSD540_",row_number())),
 subject=factor(subject),condition=factor(condition,levels=c("Baseline","ULLS","AR")),group=condition,dataset="OSD-540")
cd <- meta540 %>% dplyr::select(sample,subject,condition) %>% as.data.frame(); rownames(cd)<-cd$sample; cd$sample<-NULL
dds540 <- import_rsem(meta540$file,cd,~subject+condition); map540 <- map_ens(rownames(dds540),"human")
fit540 <- fit_deseq(dds540,map540,"OSD-540",list(
 ULLS_vs_Baseline=c("condition","ULLS","Baseline"),
 AR_vs_ULLS=c("condition","AR","ULLS"),
 AR_vs_Baseline=c("condition","AR","Baseline")))
save_de(fit540,"OSD540")
save_pca(fit540$expr_vst,meta540,"OSD-540: human vastus lateralis",
         file.path(output_dir,"figures","OSD540_PCA.svg"))

# ----- GSE271607 --------------------------------------------------------------
message("GSE271607 human dry immersion...")

x271 <- readr::read_tsv(f271,show_col_types=FALSE)
ids271 <- clean_ens(x271[[1]])

# The count file also contains annotation columns (gene_name, gene_biotype).
# Only the 40 Sol/Vas before/after columns are biological samples.
sample_cols271 <- colnames(x271)[
  stringr::str_detect(
    colnames(x271),
    stringr::regex("^(Sol|Vas)_(before|after)_\\d+$",ignore_case=TRUE)
  )
]
if(length(sample_cols271)!=40){
  stop("Expected 40 GSE271607 sample columns but found ",length(sample_cols271),
       ". Detected: ",paste(sample_cols271,collapse=", "))
}

cnt271 <- as.matrix(x271[,sample_cols271,drop=FALSE])
storage.mode(cnt271) <- "numeric"
rownames(cnt271) <- ids271
if(anyNA(cnt271)) stop("NA values detected in GSE271607 count matrix.")
if(any(cnt271<0)) stop("Negative values detected in GSE271607 count matrix.")
if(any(abs(cnt271-round(cnt271))>1e-6)){
  warning("GSE271607 contains non-integer count values; rounding for DESeq2.")
}
cnt271 <- round(cnt271)
if(anyDuplicated(rownames(cnt271))){
  cnt271 <- rowsum(cnt271,group=rownames(cnt271),reorder=FALSE)
}

meta271 <- tibble(sample=colnames(cnt271)) %>%
  mutate(
    muscle=case_when(
      str_detect(sample,regex("^Sol_",ignore_case=TRUE))~"Soleus",
      str_detect(sample,regex("^Vas_",ignore_case=TRUE))~"VastusLateralis",
      TRUE~NA_character_),
    condition=case_when(
      str_detect(sample,regex("_before_",ignore_case=TRUE))~"Before",
      str_detect(sample,regex("_after_",ignore_case=TRUE))~"After",
      TRUE~NA_character_),
    subject=str_extract(sample,"\\d+$")
  )
if(any(is.na(meta271$muscle))||any(is.na(meta271$condition))||any(is.na(meta271$subject))){
  print(meta271,n=Inf)
  stop("Cannot parse GSE271607 biological sample columns.")
}
print(meta271 %>% dplyr::count(muscle,condition))
map271 <- map_ens(rownames(cnt271),"human")

fit271muscle <- function(which_muscle,dataset){
  m <- meta271 %>%
    dplyr::filter(muscle==which_muscle) %>%
    mutate(subject=factor(subject),
           condition=factor(condition,levels=c("Before","After")),
           group=condition,dataset=dataset)

  pair_check <- m %>% dplyr::count(subject,condition)
  if(nrow(m)!=20 || any(pair_check$n!=1)){
    stop("Unexpected pairing in GSE271607 ",which_muscle," samples.")
  }

  cd271 <- m %>% dplyr::select(sample,subject,condition) %>% as.data.frame()
  rownames(cd271) <- cd271$sample
  cd271$sample <- NULL

  d <- DESeqDataSetFromMatrix(
    countData=cnt271[,m$sample,drop=FALSE],
    colData=cd271,
    design=~subject+condition
  )
  d <- d[rowSums(counts(d)>=MIN_COUNT)>=MIN_SAMPLES,]
  message(dataset,": ",nrow(d)," genes retained after expression filtering.")

  f <- fit_deseq(d,map271,dataset,
                 list(After_vs_Before=c("condition","After","Before")))
  list(fit=f,meta=m)
}

tmp <- fit271muscle("VastusLateralis","GSE271607-VL")
fit271_vl <- tmp$fit
meta271_vl <- tmp$meta

tmp <- fit271muscle("Soleus","GSE271607-Soleus")
fit271_sol <- tmp$fit
meta271_sol <- tmp$meta

save_de(fit271_vl,"GSE271607_VL")
save_de(fit271_sol,"GSE271607_Soleus")
save_pca(fit271_vl$expr_vst,meta271_vl,"GSE271607: human vastus lateralis",
         file.path(output_dir,"figures","GSE271607_VL_PCA.svg"))
save_pca(fit271_sol$expr_vst,meta271_sol,"GSE271607: human soleus",
         file.path(output_dir,"figures","GSE271607_Soleus_PCA.svg"))

# ----- Volcanoes --------------------------------------------------------------
fits <- list(OSD101=list(f=fit101,g=ALS_CORE_MOUSE),OSD419=list(f=fit419,g=ALS_CORE_MOUSE),
 OSD21=list(f=fit21,g=ALS_CORE_MOUSE),OSD684=list(f=fit684,g=ALS_CORE_HUMAN),
 OSD540=list(f=fit540,g=ALS_CORE_HUMAN),GSE271607_VL=list(f=fit271_vl,g=ALS_CORE_HUMAN),
 GSE271607_Soleus=list(f=fit271_sol,g=ALS_CORE_HUMAN))
for(px in names(fits)) for(cn in names(fits[[px]]$f$results))
 save_volcano(fits[[px]]$f$results[[cn]]$symbol,fits[[px]]$g,paste0(px,": ",cn),
              file.path(output_dir,"figures",paste0(px,"_",cn,"_volcano.svg")))


# ----- Primary contrast registry ---------------------------------------------
primary_contrasts <- tribble(
 ~dataset,~contrast,~display_name,~evidence_class,
 "OSD-101","Flight_vs_Ground","Mouse flight: OSD-101","Actual spaceflight",
 "OSD-419","Flight_vs_Ground","Mouse flight: OSD-419","Actual spaceflight",
 "OSD-21","Flight_vs_Ground","Mouse flight: OSD-21","Actual spaceflight",
 "OSD-21","HLU_vs_CTRL","Mouse HLU: OSD-21","Ground unloading",
 "OSD-21","HLU_RL_vs_HLU","Mouse reloading: OSD-21","Reloading",
 "OSD-684","Young_Flight_vs_Ground","Human engineered muscle: young","Engineered human muscle in flight",
 "OSD-684","Old_Flight_vs_Ground","Human engineered muscle: older","Engineered human muscle in flight",
 "OSD-540","ULLS_vs_Baseline","Human VL: ULLS","Human biopsy unloading",
 "GSE271607-VL","After_vs_Before","Human VL: dry immersion","Human biopsy unloading",
 "GSE271607-Soleus","After_vs_Before","Human soleus: dry immersion","Human biopsy unloading"
)
write_csv(primary_contrasts,file.path(output_dir,"summary","primary_contrasts.csv"))

# ----- GSEA -------------------------------------------------------------------
message("Running GSEA...")
mouse_pathways <- get_pathways("Mus musculus")
human_pathways <- get_pathways("Homo sapiens")
allg <- list(); allm <- list(); allcore <- list(); allbio <- list()

run_fit_gsea <- function(fit,dataset,prefix,pathways,variant="primary"){
  for(cn in names(fit$results)){
    k <- paste(dataset,cn,variant,sep="__")
    fg <- run_gsea(fit$results[[cn]]$symbol,pathways,dataset,cn,variant)
    allg[[k]] <<- fg; allm[[k]] <<- filter_metabolic(fg); allcore[[k]] <<- core_from_gsea(fg); allbio[[k]] <<- bioenergetic_from_gsea(fg)
    write_csv(fg,file.path(output_dir,"GSEA",paste0(prefix,"_",cn,"_",variant,"_GSEA_all.csv")))
    write_csv(filter_metabolic(fg),file.path(output_dir,"GSEA",paste0(prefix,"_",cn,"_",variant,"_GSEA_metabolic.csv")))
  }
}
run_fit_gsea(fit101,"OSD-101","OSD101",mouse_pathways)
run_fit_gsea(fit419,"OSD-419","OSD419",mouse_pathways)
run_fit_gsea(fit21,"OSD-21","OSD21",mouse_pathways)
run_fit_gsea(fit684,"OSD-684","OSD684",human_pathways)
run_fit_gsea(fit540,"OSD-540","OSD540",human_pathways)
run_fit_gsea(fit271_vl,"GSE271607-VL","GSE271607_VL",human_pathways)
run_fit_gsea(fit271_sol,"GSE271607-Soleus","GSE271607_Soleus",human_pathways)

gsea_all <- bind_rows(allg) %>% mutate(global_FDR=p.adjust(pval,"BH"))
gsea_met <- bind_rows(allm) %>% left_join(dplyr::select(gsea_all,dataset,contrast,variant,pathway,global_FDR),
                                           by=c("dataset","contrast","variant","pathway"))
sod_core_all <- bind_rows(allcore) %>%
  mutate(
    global_FDR=if_else(test_status=="Tested",p.adjust(pval,"BH"),NA_real_),
    significant_within_contrast=test_status=="Tested" & !is.na(padj) & padj<FDR_CUTOFF,
    significant_global=test_status=="Tested" & !is.na(global_FDR) & global_FDR<FDR_CUTOFF,
    concordant_with_SOD1=case_when(
      test_status!="Tested"~NA,
      is.na(NES)~NA,
      TRUE~NES<0),
    interpretation=case_when(
      test_status!="Tested"~test_status,
      is.na(NES)~"Exact set present in registry but no GSEA result",
      NES<0~"Concordant with SOD1 down-enrichment",
      NES>0~"Opposite to SOD1 direction",
      TRUE~"No directional shift")
  )

bioenergetic_all <- bind_rows(allbio) %>%
  mutate(
    significant_within_contrast=test_status=="Tested" & !is.na(padj) & padj<FDR_CUTOFF,
    direction=case_when(is.na(NES)~NA_character_,NES<0~"Down",NES>0~"Up",TRUE~"No shift")
  )

write_csv(gsea_all,file.path(output_dir,"summary","GSEA_all_datasets.csv"))
write_csv(gsea_met,file.path(output_dir,"summary","GSEA_metabolic_all_datasets.csv"))
write_csv(sod_core_all,file.path(output_dir,"summary","SOD1_reference_pathways_all_datasets.csv"))
write_csv(bioenergetic_all,file.path(output_dir,"summary","Bioenergetic_exact_pathways_all_datasets.csv"))
write_csv(SOD1_REFERENCE_PATHWAYS,file.path(output_dir,"ALS_reference","SOD1_fixed_pathway_registry.csv"))

# ----- SOD1 pathway heatmap ---------------------------------------------------
ch <- sod_core_all %>% inner_join(primary_contrasts,by=c("dataset","contrast")) %>%
 mutate(display_name=factor(display_name,levels=primary_contrasts$display_name),
        reference_pathway=factor(reference_pathway,levels=rev(SOD1_REFERENCE_PATHWAYS$reference_pathway)))
g <- ggplot(ch,aes(display_name,reference_pathway,fill=NES))+geom_tile()+
 geom_text(aes(label=case_when(test_status!="Tested"~"NT",is.na(NES)~"NA",
                               padj<.05~paste0(sprintf("%.2f",NES),"*"),
                               TRUE~sprintf("%.2f",NES))),size=2.4)+
 theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+
 labs(title="ALS-derived metabolic pathway response across spaceflight and unloading",
      subtitle="Exact canonical sets only. NES < 0 = SOD1-concordant; * FDR < 0.05; NT = no exact set available",
      x=NULL,y=NULL,fill="NES")
ggsave(file.path(output_dir,"figures","SOD1_derived_metabolic_pathway_NES_heatmap.svg"),g,
       width=14,height=7,device="svg")

# ----- Exact bioenergetic pathway heatmap ------------------------------------
bh <- bioenergetic_all %>%
  inner_join(primary_contrasts,by=c("dataset","contrast")) %>%
  mutate(
    display_name=factor(display_name,levels=primary_contrasts$display_name),
    bioenergetic_pathway=factor(
      bioenergetic_pathway,
      levels=rev(BIOENERGETIC_EXACT_PATHWAYS$bioenergetic_pathway)))

g <- ggplot(bh,aes(display_name,bioenergetic_pathway,fill=NES))+geom_tile()+
  geom_text(aes(label=case_when(is.na(NES)~"NA",padj<.05~paste0(sprintf("%.2f",NES),"*"),
                                TRUE~sprintf("%.2f",NES))),size=2.5)+
  theme_bw(base_size=10)+
  theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+
  labs(title="Exact bioenergetic pathway response across spaceflight and unloading",
       subtitle="Hallmark gene sets; * within-contrast FDR < 0.05",
       x=NULL,y=NULL,fill="NES")
ggsave(file.path(output_dir,"figures","Bioenergetic_exact_pathways_NES_heatmap.svg"),g,
       width=14,height=4.8,device="svg")

# ----- Exploratory metabolic heatmap -----------------------------------------
mp <- gsea_met %>% inner_join(primary_contrasts,by=c("dataset","contrast"))
tops <- mp %>% group_by(pathway) %>%
 summarise(n=n_distinct(paste(dataset,contrast)),best=min(padj,na.rm=TRUE),.groups="drop") %>%
 mutate(best=if_else(is.infinite(best),NA_real_,best)) %>% filter(n>=2) %>%
 arrange(best) %>% slice_head(n=25) %>% pull(pathway)
mh <- filter(mp,pathway%in%tops) %>%
 mutate(display_name=factor(display_name,levels=primary_contrasts$display_name))
if(nrow(mh)){
 g <- ggplot(mh,aes(display_name,fct_reorder(pathway,NES),fill=NES))+geom_tile()+
  geom_text(aes(label=ifelse(padj<.05,"*","")),size=4)+theme_bw(base_size=9)+
  theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+
  labs(title="Exploratory metabolic pathway responses",subtitle="Pre-ranked GSEA; * FDR < 0.05",
       x=NULL,y=NULL,fill="NES")
 ggsave(file.path(output_dir,"figures","Exploratory_metabolic_GSEA_NES_heatmap.svg"),g,
        width=14,height=10,device="svg")
}

# ----- OSD-21 HLU robustness -------------------------------------------------
message("OSD-21 HLU robustness...")
M11 <- "Mmus_C57-6J_GST-M11_HLU_Rep1"
core_noM11 <- tibble()
if(M11%in%meta21$sample){
 m <- filter(meta21,sample!=M11) %>% droplevels()
 f <- fit_limma(expr21[,m$sample,drop=FALSE],m,c(HLU_vs_CTRL="HLU-CTRL"),"OSD-21")
 fg <- run_gsea(f$results$HLU_vs_CTRL$symbol,mouse_pathways,"OSD-21","HLU_vs_CTRL","M11_excluded")
 core_noM11 <- core_from_gsea(fg)
 write_csv(f$results$HLU_vs_CTRL$symbol,file.path(output_dir,"robustness","OSD21_HLU_vs_CTRL_M11_excluded_DE.csv"))
 write_csv(fg,file.path(output_dir,"robustness","OSD21_HLU_vs_CTRL_M11_excluded_GSEA.csv"))
 write_csv(core_noM11,file.path(output_dir,"robustness","OSD21_HLU_vs_CTRL_M11_excluded_SOD1_reference.csv"))
}
hlu <- filter(meta21,group=="HLU")$sample
loo <- list()
for(drop in hlu){
 m <- filter(meta21,sample!=drop) %>% droplevels()
 f <- fit_limma(expr21[,m$sample,drop=FALSE],m,c(HLU_vs_CTRL="HLU-CTRL"),"OSD-21")
 fg <- run_gsea(f$results$HLU_vs_CTRL$symbol,mouse_pathways,"OSD-21","HLU_vs_CTRL",paste0("drop_",drop))
 loo[[drop]] <- core_from_gsea(fg) %>% mutate(omitted_sample=drop)
}
loo_all <- bind_rows(loo)
write_csv(loo_all,file.path(output_dir,"robustness","HLU_leave_one_out_SOD1_reference_results.csv"))
hlu_rob <- loo_all %>%
 dplyr::filter(test_status=="Tested") %>%
 group_by(reference_pathway,msigdb_set) %>%
 summarise(n_leave_one_out_runs=n(),median_NES=median(NES,na.rm=TRUE),
           min_NES=min(NES,na.rm=TRUE),max_NES=max(NES,na.rm=TRUE),
           proportion_negative=mean(NES<0,na.rm=TRUE),
           proportion_padj_lt_0_05=mean(padj<.05,na.rm=TRUE),
           proportion_concordant_significant=mean(NES<0&padj<.05,na.rm=TRUE),.groups="drop")
write_csv(hlu_rob,file.path(output_dir,"robustness","HLU_SOD1_reference_robustness_summary.csv"))

# ----- Optional OSD-419 STAR sensitivity -------------------------------------
star_summary <- tibble()
if(!is.na(f419_star)&&file.exists(f419_star)){
 message("OSD-419 STAR sensitivity using ",basename(f419_star))
 sx <- read_csv(f419_star,show_col_types=FALSE); colnames(sx)[1]<-"ensembl_id"; sx$ensembl_id<-clean_ens(sx$ensembl_id)
 ss <- setdiff(colnames(sx),"ensembl_id")
 sm <- tibble(sample=ss,group=case_when(
   str_detect(sample,regex("_FLT_|-FLT_",ignore_case=TRUE))~"Flight",
   str_detect(sample,regex("_GC_|-GC_",ignore_case=TRUE))~"Ground",TRUE~NA_character_))
 if(!any(is.na(sm$group))){
  cm <- as.matrix(sx[,ss,drop=FALSE]); storage.mode(cm)<-"numeric"; rownames(cm)<-sx$ensembl_id
  sm <- mutate(sm,group=factor(group,levels=c("Ground","Flight")))
  cd <- as.data.frame(sm); rownames(cd)<-cd$sample; cd$sample<-NULL
  dd <- DESeqDataSetFromMatrix(round(cm),cd,~group); dd<-dd[rowSums(counts(dd)>=MIN_COUNT)>=MIN_SAMPLES,]
  mp2 <- map_ens(rownames(dd),"mouse")
  fs <- fit_deseq(dd,mp2,"OSD-419-STAR",list(Flight_vs_Ground=c("group","Flight","Ground")))
  comp <- inner_join(dplyr::select(fit419$results$Flight_vs_Ground$symbol,symbol,RSEM=log2FC),
                     dplyr::select(fs$results$Flight_vs_Ground$symbol,symbol,STAR=log2FC),by="symbol")
  rho <- cor(comp$RSEM,comp$STAR,method="spearman",use="complete.obs")
  star_summary <- tibble(comparison="OSD-419 RSEM vs STAR",n_genes=nrow(comp),spearman_rho=rho)
  write_csv(comp,file.path(output_dir,"robustness","OSD419_RSEM_vs_STAR_gene_effects.csv"))
  write_csv(star_summary,file.path(output_dir,"robustness","OSD419_RSEM_vs_STAR_summary.csv"))
  fg <- run_gsea(fs$results$Flight_vs_Ground$symbol,mouse_pathways,"OSD-419-STAR","Flight_vs_Ground","STAR_sensitivity")
  write_csv(core_from_gsea(fg),file.path(output_dir,"robustness","OSD419_STAR_SOD1_reference.csv"))
 }
} else message("No STAR-labelled OSD-419 merged count matrix found; STAR sensitivity skipped.")

# ----- Cross-dataset gene-effect comparisons ---------------------------------
cors <- bind_rows(
 save_scatter(fit21$results$Flight_vs_Ground$symbol,fit21$results$HLU_vs_CTRL$symbol,
              "OSD-21 flight log2FC","OSD-21 HLU log2FC","OSD-21: flight vs HLU",
              file.path(output_dir,"figures","OSD21_Flight_vs_HLU_gene_effects.svg")),
 save_scatter(fit101$results$Flight_vs_Ground$symbol,fit21$results$Flight_vs_Ground$symbol,
              "OSD-101 flight log2FC","OSD-21 flight log2FC","Independent mouse flight: OSD-101 vs OSD-21",
              file.path(output_dir,"figures","OSD101_vs_OSD21_Flight_gene_effects.svg")),
 save_scatter(fit101$results$Flight_vs_Ground$symbol,fit419$results$Flight_vs_Ground$symbol,
              "OSD-101 left gastroc log2FC","OSD-419 right gastroc log2FC","RR-1 left/right concordance",
              file.path(output_dir,"figures","OSD101_vs_OSD419_gene_effects.svg")),
 save_scatter(fit540$results$ULLS_vs_Baseline$symbol,fit271_vl$results$After_vs_Before$symbol,
              "Human VL ULLS log2FC","Human VL dry-immersion log2FC","Human vastus lateralis unloading",
              file.path(output_dir,"figures","Human_VL_OSD540_vs_GSE271607_gene_effects.svg")),
 save_scatter(fit271_vl$results$After_vs_Before$symbol,fit271_sol$results$After_vs_Before$symbol,
              "Human VL dry-immersion log2FC","Human soleus dry-immersion log2FC","Dry immersion: VL vs soleus",
              file.path(output_dir,"figures","GSE271607_VL_vs_Soleus_gene_effects.svg"))
)
write_csv(cors,file.path(output_dir,"summary","cross_dataset_gene_effect_correlations.csv"))

# ----- DE counts --------------------------------------------------------------
sumde <- function(f,dataset) map_dfr(names(f$results),function(cn){
 x<-f$results[[cn]]$symbol
 tibble(dataset=dataset,contrast=cn,genes_tested=nrow(x),
 FDR_lt_0_05=sum(x$FDR<.05,na.rm=TRUE),
 FDR_lt_0_05_log2FC_ge_1=sum(x$FDR<.05&abs(x$log2FC)>=1,na.rm=TRUE),
 up_FDR_LFC=sum(x$FDR<.05&x$log2FC>=1,na.rm=TRUE),
 down_FDR_LFC=sum(x$FDR<.05&x$log2FC<=-1,na.rm=TRUE))
})
deg <- bind_rows(sumde(fit101,"OSD-101"),sumde(fit419,"OSD-419"),sumde(fit21,"OSD-21"),
                 sumde(fit684,"OSD-684"),sumde(fit540,"OSD-540"),
                 sumde(fit271_vl,"GSE271607-VL"),sumde(fit271_sol,"GSE271607-Soleus"))
write_csv(deg,file.path(output_dir,"summary","DEG_summary_counts.csv"))

# ----- SOD1 concordance summary ----------------------------------------------
sod_conc <- sod_core_all %>%
 inner_join(primary_contrasts,by=c("dataset","contrast")) %>%
 group_by(dataset,contrast,display_name,evidence_class) %>%
 summarise(
   SOD1_reference_pathways_total=n(),
   exact_pathways_tested=sum(test_status=="Tested"),
   exact_pathways_not_tested=sum(test_status!="Tested"),
   pathways_negative=sum(test_status=="Tested" & NES<0,na.rm=TRUE),
   pathways_positive=sum(test_status=="Tested" & NES>0,na.rm=TRUE),
   pathways_negative_FDR05=sum(test_status=="Tested" & NES<0 & padj<.05,na.rm=TRUE),
   pathways_positive_FDR05=sum(test_status=="Tested" & NES>0 & padj<.05,na.rm=TRUE),
   median_NES=ifelse(sum(test_status=="Tested")>0,median(NES[test_status=="Tested"],na.rm=TRUE),NA_real_),
   proportion_concordant=ifelse(sum(test_status=="Tested")>0,mean(NES[test_status=="Tested"]<0,na.rm=TRUE),NA_real_),
   .groups="drop")
write_csv(sod_conc,file.path(output_dir,"summary","SOD1_reference_concordance_summary.csv"))

# ----- Metadata / QC ----------------------------------------------------------
mq <- bind_rows(mutate(map101,dataset="OSD-101",species="mouse"),
                mutate(map419,dataset="OSD-419",species="mouse"),
                mutate(map684,dataset="OSD-684",species="human"),
                mutate(map540,dataset="OSD-540",species="human"),
                mutate(map271,dataset="GSE271607",species="human")) %>%
 mutate(mapped=!is.na(symbol)&symbol!="")
mqs <- mq %>% group_by(dataset,species) %>%
 summarise(n_Ensembl=n(),n_mapped=sum(mapped),percent_mapped=100*mean(mapped),.groups="drop")
write_csv(mqs,file.path(output_dir,"QC","Ensembl_mapping_summary.csv"))
write_csv(mutate(meta101,file=basename(file),group=as.character(group)),file.path(output_dir,"QC","OSD101_sample_metadata_used.csv"))
write_csv(mutate(meta419,file=basename(file),group=as.character(group)),file.path(output_dir,"QC","OSD419_sample_metadata_used.csv"))
write_csv(mutate(meta21,group=as.character(group)),file.path(output_dir,"QC","OSD21_sample_metadata_used.csv"))
write_csv(mutate(meta684,file=basename(file),group=as.character(group)),file.path(output_dir,"QC","OSD684_sample_metadata_used.csv"))
write_csv(mutate(meta540,file=basename(file),group=as.character(group),condition=as.character(condition)),file.path(output_dir,"QC","OSD540_sample_metadata_used.csv"))
write_csv(meta271,file.path(output_dir,"QC","GSE271607_sample_metadata_used.csv"))

# ----- Summary workbook -------------------------------------------------------
message("Writing summary workbook...")
wb <- createWorkbook()
tabs <- list(
 SOD1_reference=sod_ref,SOD1_fixed_sets=SOD1_REFERENCE_PATHWAYS,SOD1_KEGG_raw=sod,
 SOD1_pathways_all=sod_core_all,SOD1_concordance=sod_conc,Bioenergetic_exact=bioenergetic_all,DEG_counts=deg,
 Metabolic_GSEA=gsea_met,Effect_correlations=cors,HLU_LOO_robustness=hlu_rob,
 Mapping_QC=mqs,OSD101_metadata=mutate(meta101,file=basename(file)),
 OSD419_metadata=mutate(meta419,file=basename(file)),OSD21_metadata=meta21,
 OSD684_metadata=mutate(meta684,file=basename(file)),OSD540_metadata=mutate(meta540,file=basename(file)),
 GSE271607_metadata=meta271
)
if(nrow(star_summary)) tabs$OSD419_STAR_sens <- star_summary
for(n in names(tabs)){ addWorksheet(wb,n); writeData(wb,n,tabs[[n]]) }
saveWorkbook(wb,file.path(output_dir,"summary","NASA_ALS_pathway_summary.xlsx"),overwrite=TRUE)

params <- tibble(
 parameter=c("NASA_RNAseq_primary_input","NASA_RNAseq_DE_method","OSD540_design",
             "GSE271607_input","GSE271607_design","OSD21_method","ALS_primary_reference",
             "MIN_COUNT","MIN_SAMPLES","FDR_CUTOFF","LFC_CUTOFF","M11_sample"),
 value=c("rRNArm per-sample RSEM .genes.results only","tximport(RSEM)+DESeq2",
         "~ subject + condition","GSE271607 raw gene Read_counts",
         "~ subject + condition; VL and soleus separately",
         "limma on normalized probeset expression collapsed to symbol",
         "SOD1G93A mRNA KEGG S7-S12; exact canonical MSigDB sets only for formal cross-dataset pathway tests",
         
         MIN_COUNT,MIN_SAMPLES,FDR_CUTOFF,LFC_CUTOFF,M11)
)
write_csv(params,file.path(output_dir,"QC","analysis_parameters.csv"))
capture.output(sessionInfo(),file=file.path(output_dir,"QC","sessionInfo.txt"))

message("")
message("===============================================================")
message("FULL ANALYSIS COMPLETE")
message("Results: ",normalizePath(output_dir,winslash="/",mustWork=FALSE))
message("Start with:")
message("  summary/NASA_ALS_pathway_summary.xlsx")
message("  summary/SOD1_reference_concordance_summary.csv")
message("  summary/SOD1_reference_pathways_all_datasets.csv")
message("  summary/Bioenergetic_exact_pathways_all_datasets.csv")
message("  robustness/HLU_SOD1_reference_robustness_summary.csv")
message("  figures/SOD1_derived_metabolic_pathway_NES_heatmap.svg")
message("===============================================================")
