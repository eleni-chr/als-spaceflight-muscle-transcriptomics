# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# Independent human ALS validation of the FROZEN 174-gene ALS-unloading module
# Version 1.0, 2026-10-10
# Run in RStudio: source("12_human_als_module_validation_analysis.R")
# Set the working directory to the project folder containing the downloaded
# GEO data. Internet access may be needed for package and GEO metadata retrieval.
# Reads the previous v16 module CSV, or sheet Module_directional of its workbook.
# That module is ALREADY in human orthologue symbol space: no new mapping is done.
# Human cohorts are NEVER used to select/redefine the module.
#
# Outputs: human_ALS_module_validation/{tables,figures,QC,cache} plus an XLSX.
# Positive log2FC/stat/NES = higher in ALS than control. Expected: SOD1_DOWN
# negative enrichment; SOD1_UP positive enrichment. Mixed module not used as an
# unsigned GSEA set. fgsea is competitive; CAMERA adds correlation-aware testing.
# Primary models: ~condition. Available age/sex and muscle adjustment are separate
# sensitivities, never silently substituted for the primary model. Muscle-adjusted
# effects require an additive model; deltoid-only microarray is also analysed.
# Cross-gene direction counts/permutations are competitive, exploratory evidence,
# not independent-patient tests. No pooling of expression across platforms.

PROJECT_DIR <- getwd()
RNA_DIR <- file.path(PROJECT_DIR, 'GSE215364')
ARRAY_DIR <- file.path(PROJECT_DIR, 'GSE41414')
OUTDIR <- file.path(PROJECT_DIR, 'human_ALS_module_validation')
# Set explicitly only if your earlier output is in another folder. CSV or XLSX.
MODULE_FILE <- NA_character_
MODULE_SHEET <- 'Module_directional'
# Optional CSVs with sample, condition (Control/ALS), muscle, age, sex columns.
# sample MUST be a GSM accession. These override GEO metadata by sample ID.
RNA_METADATA_OVERRIDE <- file.path(RNA_DIR, 'sample_metadata.csv')
ARRAY_METADATA_OVERRIDE <- file.path(ARRAY_DIR, 'sample_metadata.csv')
INSTALL_MISSING <- TRUE
ARRAY_INPUT <- 'auto' # auto: prefer CEL, otherwise GEO deposited log2 RMA matrix
N_PERM <- 20000L
SEED <- 20261010L
MIN_GSEA_SIZE <- 10L
LFC_THRESHOLDS <- c(0, 0.10)
options(timeout = max(600, getOption('timeout')), repos = c(CRAN='https://cloud.r-project.org'))
set.seed(SEED)

install_needed <- function(x, bioc=FALSE) {
  missing <- x[!vapply(x, requireNamespace, logical(1), quietly=TRUE)]
  if (!length(missing)) return(invisible(NULL))
  if (!INSTALL_MISSING) stop('Install packages: ', paste(missing, collapse=', '))
  if (bioc) {
    if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
    BiocManager::install(missing, ask=FALSE, update=FALSE)
  } else install.packages(missing)
}
install_needed(c('readr','readxl','openxlsx','ggplot2'))
install_needed(c('DESeq2','limma','edgeR','fgsea','GEOquery','Biobase',
                 'AnnotationDbi','org.Hs.eg.db','hgfocus.db'), TRUE)
for (p in c('tables','figures','QC','cache')) dir.create(file.path(OUTDIR,p), recursive=TRUE, showWarnings=FALSE)
TD <- file.path(OUTDIR,'tables'); FD <- file.path(OUTDIR,'figures')
QD <- file.path(OUTDIR,'QC'); CD <- file.path(OUTDIR,'cache')
writeLines(capture.output(sessionInfo()), file.path(QD,'sessionInfo_start.txt'))
tables <- list(); model_log <- list()
emit <- function(x, name, qc=FALSE) {
  x <- as.data.frame(x)
  for (j in seq_along(x)) if (is.list(x[[j]])) x[[j]] <- vapply(x[[j]], paste, character(1), collapse=';')
  readr::write_csv(x, file.path(if (qc) QD else TD, paste0(name,'.csv')), na='')
  if (!qc) tables[[name]] <<- x
  invisible(x)
}
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)
clean <- function(x) toupper(trimws(as.character(x)))
gsm <- function(x) {
  z <- regmatches(x, regexpr('GSM[0-9]+', x, ignore.case=TRUE))
  if (length(z) != length(x)) stop('Could not extract a GSM ID from every filename.')
  toupper(z)
}
saveplot <- function(p, name, w=8, h=5) {
  ggplot2::ggsave(file.path(FD,paste0(name,'.pdf')),p,width=w,height=h)
  ggplot2::ggsave(file.path(FD,paste0(name,'.png')),p,width=w,height=h,dpi=300)
}

# 1. Frozen module from the EXACT previous human-orthologue output.
load_module <- function() {
  f <- MODULE_FILE
  if (is.na(f)) {
    candidates <- list.files(PROJECT_DIR, recursive=TRUE, full.names=TRUE,
      pattern='^reversible_ALS_unloading_module_DIRECTIONAL\\.csv$')
    if (!length(candidates)) candidates <- list.files(PROJECT_DIR, recursive=TRUE,
      full.names=TRUE, pattern='^reversible_ALS_unloading_results.*\\.xlsx$')
    assert(length(candidates)==1L, paste0('Set MODULE_FILE at the top to your previous ',
      'reversible_ALS_unloading_module_DIRECTIONAL.csv or results workbook. Found ',
      length(candidates),' candidates:\n',paste(candidates,collapse='\n')))
    f <- candidates[1]
  }
  assert(file.exists(f),paste('Module file not found:',f))
  x <- if (grepl('\\.xlsx$',f,ignore.case=TRUE)) as.data.frame(readxl::read_excel(f,sheet=MODULE_SHEET)) else
    as.data.frame(readr::read_csv(f,show_col_types=FALSE))
  if ('directional_module' %in% names(x)) x <- x[tolower(as.character(x$directional_module)) %in% c('true','1'),,drop=FALSE]
  assert(all(c('symbol','als_direction') %in% names(x)),
    'Expected v16 module columns symbol and als_direction (already human symbols).')
  m <- data.frame(symbol=clean(x$symbol), als_direction=as.integer(x$als_direction))
  assert(!anyNA(m) && !anyDuplicated(m$symbol),'Module has missing/duplicate symbols or directions.')
  assert(nrow(m)==174L && sum(m$als_direction==-1L)==128L && sum(m$als_direction==1L)==46L,
    'Frozen module must contain exactly 174 unique human genes: 128 down and 46 up. Check the selected input.')
  m$component <- ifelse(m$als_direction==-1,'SOD1_DOWN','SOD1_UP')
  emit(data.frame(file=normalizePath(f,winslash='/'),md5=unname(tools::md5sum(f))), 'module_source',TRUE)
  emit(m,'frozen_module'); m
}
module <- load_module()
sets <- split(module$symbol,module$component)

# 2. GEO metadata and explicit sample matching (never infer group from order).
EXPECTED_RNA <- paste0('GSM',6634377:6634386)
EXPECTED_ARRAY <- paste0('GSM',c(1016766:1016775,1016777,1016778,1016780,1016782))
extract_field <- function(pd, pattern) {
  vapply(seq_len(nrow(pd)),function(i) {
    row <- as.character(unlist(pd[i,,drop=FALSE],use.names=FALSE))
    hit <- row[grepl(paste0('^(',pattern,')\\s*:'),row,ignore.case=TRUE)]
    if (length(hit)) trimws(sub('^[^:]*:','',hit[1])) else NA_character_
  },character(1))
}
get_geo <- function(acc, folder) {
  cache <- file.path(CD,paste0(acc,'_GEO.rds'))
  if (file.exists(cache)) return(readRDS(cache))
  local <- list.files(folder,pattern='series_matrix.*\\.txt(\\.gz)?$',full.names=TRUE,recursive=TRUE,ignore.case=TRUE)
  z <- if (length(local)==1) GEOquery::getGEO(filename=local[1], getGPL=FALSE) else
    GEOquery::getGEO(acc,GSEMatrix=TRUE,getGPL=FALSE,destdir=CD)
  if (is.list(z)) {
    assert(length(z)==1L,paste(acc,'has multiple matrices; select the correct platform explicitly.'))
    z <- z[[1]]
  }
  assert(inherits(z,'ExpressionSet'),paste(acc,'did not return a GEO expression matrix.'))
  saveRDS(z,cache); z
}
metadata <- function(es, acc, override) {
  pd <- Biobase::pData(es)
  emit(cbind(sample=rownames(pd),pd),paste0(acc,'_GEO_metadata_raw'),TRUE)
  ids <- rownames(pd)
  condition <- extract_field(pd,'disease state|disease status|disease|diagnosis')
  title <- as.character(pd$title)
  condition[grepl('ALS',title,ignore.case=TRUE)] <- 'ALS'
  condition[grepl('control|CTR',title,ignore.case=TRUE)] <- 'Control'
  condition <- ifelse(grepl('^ALS$|amyotrophic',condition,ignore.case=TRUE),'ALS',
    ifelse(grepl('healthy|normal|control',condition,ignore.case=TRUE),'Control',NA_character_))
  muscle <- extract_field(pd,'muscle|biopsied muscle|tissue')
  source <- if ('source_name_ch1' %in% names(pd)) as.character(pd$source_name_ch1) else rep('',nrow(pd))
  for (site in c('deltoid','quadriceps','biceps','tibialis anterior')) {
    take <- grepl(site,paste(muscle,source),ignore.case=TRUE)
    muscle[take] <- site
  }
  muscle[is.na(muscle) | clean(muscle) %in% c('MUSCLE','SKELETAL MUSCLE','')] <- NA_character_
  sex <- clean(extract_field(pd,'sex|gender'))
  sex <- ifelse(sex %in% c('M','MALE'),'M',ifelse(sex %in% c('F','FEMALE'),'F',NA))
  age <- suppressWarnings(as.numeric(extract_field(pd,'age|age \\(years\\)')))
  out <- data.frame(sample=ids,condition=condition,muscle=muscle,age=age,sex=sex)
  if (file.exists(override)) {
    ov <- as.data.frame(readr::read_csv(override,show_col_types=FALSE))
    assert('sample' %in% names(ov) && !anyDuplicated(ov$sample),'Metadata override needs unique sample GSM IDs.')
    assert(all(ov$sample %in% ids),'Unknown GSM ID in metadata override.')
    for (nm in intersect(names(out),setdiff(names(ov),'sample'))) {
      ix <- match(ov$sample,out$sample); good <- !is.na(ov[[nm]]) & ov[[nm]]!=''
      out[ix[good],nm] <- ov[[nm]][good]
    }
  }
  expected <- if (acc=='GSE215364') EXPECTED_RNA else EXPECTED_ARRAY
  assert(setequal(out$sample,expected),paste(acc,'sample inventory differs from the expected cohort.'))
  assert(!anyNA(out$condition) && all(out$condition %in% c('Control','ALS')),'Unresolved condition labels. Edit sample_metadata.csv.')
  assert(sum(out$condition=='ALS')==length(expected)/2 && sum(out$condition=='Control')==length(expected)/2,
    paste(acc,'does not have the expected equal ALS/control groups.'))
  control_ids <- if (acc=='GSE215364') EXPECTED_RNA[1:5] else EXPECTED_ARRAY[1:7]
  assert(setequal(out$sample[out$condition=='Control'],control_ids),
    paste(acc,'condition assignments conflict with the deposited GSM cohort labels.'))
  out$condition <- factor(out$condition,levels=c('Control','ALS'))
  rownames(out) <- out$sample
  emit(out,paste0(acc,'_sample_metadata'),TRUE)
  # A template is saved for missing muscle/covariates; no invented covariates.
  readr::write_csv(out,file.path(QD,paste0(acc,'_metadata_override_template.csv')))
  out
}
unpack <- function(folder,acc) {
  assert(dir.exists(folder),paste('Directory not found:',folder))
  tarfiles <- list.files(folder,pattern=paste0('^',acc,'_RAW\\.tar$'),full.names=TRUE,ignore.case=TRUE)
  dest <- file.path(folder,paste0(acc,'_RAW_extracted'))
  if (length(tarfiles)) {
    marker <- file.path(dest,'extraction_complete.txt')
    if (!file.exists(marker)) {
      dir.create(dest,showWarnings=FALSE)
      entries <- utils::untar(tarfiles[1],list=TRUE)
      assert(!any(grepl('(^[/\\\\])|(^|/)\\.\\.(/|$)',entries)), 'Unsafe paths in archive.')
      utils::untar(tarfiles[1],exdir=dest)
      writeLines(unname(tools::md5sum(tarfiles[1])),marker)
    } else assert(readLines(marker,n=1)==unname(tools::md5sum(tarfiles[1])),
      'Archive changed: delete the RAW_extracted folder and rerun.')
  }
}
unpack(RNA_DIR,'GSE215364'); unpack(ARRAY_DIR,'GSE41414')
rna_es <- get_geo('GSE215364',RNA_DIR); array_es <- get_geo('GSE41414',ARRAY_DIR)
rna_meta <- metadata(rna_es,'GSE215364',RNA_METADATA_OVERRIDE)
array_meta <- metadata(array_es,'GSE41414',ARRAY_METADATA_OVERRIDE)

# 3. Gene annotation: drop ambiguous ID-to-symbol links, log every mapping.
map_ids <- function(ids, db, type) {
  valid <- intersect(unique(ids),AnnotationDbi::keys(db,keytype=type))
  assert(length(valid)>0,paste('No IDs found in annotation:',type))
  z <- AnnotationDbi::select(db,keys=valid,keytype=type,columns='SYMBOL')
  z <- unique(data.frame(id=as.character(z[[type]]),symbol=clean(z$SYMBOL)))
  z <- z[!is.na(z$symbol) & nzchar(z$symbol),,drop=FALSE]
  n <- table(z$id); z <- z[z$id %in% names(n)[n==1],,drop=FALSE]
  out <- data.frame(id=ids,symbol=z$symbol[match(ids,z$id)])
  out
}
read_counts <- function(f) {
  con <- if (grepl('\\.gz$',f)) gzfile(f,'rt') else file(f,'rt')
  preview <- readLines(con,n=30,warn=FALSE); close(con)
  preview <- preview[!grepl('^#|^\\s*$',preview)]
  assert(length(preview)>0,paste('Empty count file:',f))
  first <- strsplit(preview[1],'\t')[[1]]
  header <- !grepl('^ENSG[0-9]+(\\.[0-9]+)?$',first[1]) &&
    !(length(first)==2 && is.finite(suppressWarnings(as.numeric(first[2]))))
  con <- if (grepl('\\.gz$',f)) gzfile(f,'rt') else file(f,'rt')
  on.exit(close(con))
  x <- read.delim(con,header=header,comment.char='#',check.names=FALSE,stringsAsFactors=FALSE)
  assert(ncol(x)>=2,paste('Expected tab-separated gene counts:',f))
  if (header && 'Geneid' %in% names(x)) {
    idcol <- match('Geneid',names(x)); values <- setdiff(seq_len(ncol(x)),match(c('Geneid','Chr','Start','End','Strand','Length'),names(x)))
  } else { idcol <- 1L; values <- setdiff(seq_len(ncol(x)),idcol) }
  assert(length(values)==1L,paste('Expected one sample count column; inspect:',f))
  ids <- sub('\\.[0-9]+$','',trimws(as.character(x[[idcol]])))
  count <- suppressWarnings(as.numeric(x[[values]]))
  assert(all(is.finite(count)) && all(count>=0) && all(abs(count-round(count))<1e-7),
    paste('Counts must be finite nonnegative integers; no TPM/FPKM/fractional rounding:',f))
  assert(!anyDuplicated(ids),paste('Duplicate gene IDs:',f))
  data.frame(id=ids,count=count)
}
rf <- list.files(RNA_DIR,pattern='geneCOUNT.*\\.txt(\\.gz)?$',recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
rid <- gsm(basename(rf))
assert(length(rf)==10L && !anyDuplicated(rid) && setequal(rid,EXPECTED_RNA),
  'Need exactly ten GSM*_geneCOUNT.txt(.gz) files. Avoid duplicate compressed/uncompressed copies.')
rf <- rf[match(rna_meta$sample,rid)]
rc <- lapply(rf,read_counts)
assert(all(vapply(rc,function(x) identical(x$id,rc[[1]]$id),logical(1))),
  'Gene rows differ across count files; align/inspect inputs explicitly before analysis.')
counts <- do.call(cbind,lapply(rc,`[[`,'count')); rownames(counts) <- rc[[1]]$id
colnames(counts) <- rna_meta$sample
if (all(grepl('^ENSG[0-9]+$',rownames(counts)))) {
  rmap <- map_ids(rownames(counts),org.Hs.eg.db::org.Hs.eg.db,'ENSEMBL')
} else if (all(grepl('^[0-9]+$',rownames(counts)))) {
  rmap <- map_ids(rownames(counts),org.Hs.eg.db::org.Hs.eg.db,'ENTREZID')
} else {
  # Only confirmed official symbols are accepted, not a case-converted mouse list.
  sy <- clean(rownames(counts)); valid <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype='SYMBOL')
  rmap <- data.frame(id=rownames(counts),symbol=ifelse(sy %in% valid,sy,NA))
}
emit(rmap,'GSE215364_ID_mapping',TRUE)
ok <- !is.na(rmap$symbol)
counts <- rowsum(counts[ok,,drop=FALSE],rmap$symbol[ok],reorder=FALSE)
assert(all(is.finite(counts)) && max(counts)<=.Machine$integer.max,
  'Aggregated counts exceed integer range or are non-finite.')
storage.mode(counts) <- 'integer'
emit(data.frame(sample=colnames(counts),library_size=colSums(counts)),'GSE215364_library_sizes',TRUE)
assert(all(colSums(counts)>0),'Zero RNA-seq library size.')

# 4. Microarray: raw CEL -> RMA, otherwise use deposited log2 RMA unchanged.
cf <- list.files(ARRAY_DIR,pattern='\\.CEL(\\.gz)?$',full.names=TRUE,recursive=TRUE,ignore.case=TRUE)
assert(ARRAY_INPUT %in% c('auto','CEL','matrix'),'ARRAY_INPUT must be auto, CEL or matrix.')
if (ARRAY_INPUT=='CEL' || (ARRAY_INPUT=='auto' && length(cf))) {
  install_needed(c('affy','hgfocuscdf'),TRUE)
  cid <- gsm(basename(cf))
  assert(length(cf)==14L && !anyDuplicated(cid) && setequal(cid,EXPECTED_ARRAY),
    'Need exactly 14 uniquely GSM-labelled CEL files. Avoid duplicate compressed/uncompressed copies.')
  cf <- cf[match(array_meta$sample,cid)]
  raw <- affy::ReadAffy(filenames=cf)
  chip_name <- affy::cdfName(raw)
  message("Detected CEL chip type: ", chip_name)
  
  chip_key <- gsub("[^a-z0-9]", "", tolower(chip_name))
  
  assert(
    length(chip_key) == 1L && chip_key %in% c("hgfocus", "focus"),
    paste0("Unexpected CEL chip type: ", chip_name,
           ". Expected HG-Focus.")
  )
  aes <- affy::rma(raw); array <- Biobase::exprs(aes)
  colnames(array) <- array_meta$sample
  array_source <- 'Raw CEL: affy RMA jointly across 14 arrays'
} else {
  array <- Biobase::exprs(array_es)[,array_meta$sample,drop=FALSE]
  assert(nrow(array)==8793L && all(is.finite(array)) && max(array)<30 && min(array)>-10,
    'Deposited matrix does not resemble the expected 8793-probeset log2 RMA matrix. Inspect before proceeding.')
  array_source <- 'Deposited log2 RMA matrix: no second log transform or normalization'
}
assert(all(is.finite(array)),'Non-finite microarray expression.')
amap <- map_ids(rownames(array),hgfocus.db::hgfocus.db,'PROBEID')
amap$mean_expression <- rowMeans(array)
# One probe per symbol: highest mean signal, chosen independently of case-control DE.
ord <- order(amap$symbol,-amap$mean_expression,amap$id,na.last=TRUE)
selected <- ord[!is.na(amap$symbol[ord]) & !duplicated(amap$symbol[ord])]
amap$selected <- seq_len(nrow(amap)) %in% selected
emit(amap,'GSE41414_probe_mapping',TRUE)
array <- array[selected,,drop=FALSE]; rownames(array) <- amap$symbol[selected]
emit(data.frame(source=array_source,n_gene_symbols=nrow(array)),'GSE41414_preprocessing',TRUE)
inputs <- c(rf,cf)
emit(data.frame(file=inputs,md5=unname(tools::md5sum(inputs))),'local_input_manifest',TRUE)

# 5. Primary DE and explicitly labelled covariate/site sensitivities.
check_design <- function(meta, formula) {
  mm <- model.matrix(formula,meta)
  assert(nrow(mm)==nrow(meta),'Missing covariate values would drop samples.')
  assert(qr(mm)$rank==ncol(mm),'Design is not full rank; condition is confounded with covariates.')
  assert(nrow(mm)-ncol(mm)>=3L,'Fewer than three residual degrees of freedom.')
  mm
}
fits <- list()
fit_one <- function(name,mat,meta,formula,rna=FALSE) {
  message('Fitting ',name,' ',deparse(formula))
  mat <- mat[,meta$sample,drop=FALSE]
  design <- check_design(meta,formula)
  if (rna) {
    y <- edgeR::DGEList(mat)
    keep <- edgeR::filterByExpr(y,design=design)
    assert(sum(keep)>1000,'Too few expressed RNA-seq genes.')
    dds <- DESeq2::DESeqDataSetFromMatrix(mat[keep,,drop=FALSE],meta,formula)
    dds <- DESeq2::DESeq(dds,quiet=TRUE)
    rr <- as.data.frame(DESeq2::results(dds,contrast=c('condition','ALS','Control'),alpha=0.05))
    de <- data.frame(symbol=rownames(rr),log2FC=rr$log2FoldChange,stat=rr$stat,
      p_value=rr$pvalue,FDR=rr$padj,mean_expression=rr$baseMean)
    expr <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds,blind=FALSE))
    # CAMERA sensitivity uses voom on the same filtered genes/design.
    y <- edgeR::calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
    v <- limma::voom(y,design,plot=FALSE)
    camera_input <- v
    emit(data.frame(symbol=rownames(mat),expression_pass=keep),paste0(name,'_expression_filter'),TRUE)
    saveRDS(dds,file.path(CD,paste0(name,'_dds.rds')))
  } else {
    fit <- limma::eBayes(limma::lmFit(mat,design),trend=FALSE)
    coef <- match('conditionALS',colnames(design))
    rr <- limma::topTable(fit,coef=coef,number=Inf,sort.by='none')
    de <- data.frame(symbol=rownames(rr),log2FC=rr$logFC,stat=rr$t,
      p_value=rr$P.Value,FDR=rr$adj.P.Val,mean_expression=rr$AveExpr)
    expr <- mat; camera_input <- mat
    saveRDS(fit,file.path(CD,paste0(name,'_limma.rds')))
  }
  coef <- match('conditionALS',colnames(design))
  assert(!is.na(coef),'ALS-vs-Control coefficient missing.')
  emit(de,paste0(name,'_DE'))
  emit(data.frame(sample=meta$sample,design,check.names=FALSE),paste0(name,'_design'),TRUE)
  model_log[[name]] <<- data.frame(analysis=name,formula=deparse(formula),n_ALS=sum(meta$condition=='ALS'),
    n_control=sum(meta$condition=='Control'),n_genes=nrow(de),status='completed')
  list(de=de,expr=expr,meta=meta,design=design,camera_input=camera_input,coef=coef)
}
add_fit <- function(name,mat,meta,formula,rna=FALSE,optional=FALSE) {
  if (!optional) fits[[name]] <<- fit_one(name,mat,meta,formula,rna) else {
    z <- tryCatch(fit_one(name,mat,meta,formula,rna),error=function(e) {
      model_log[[name]] <<- data.frame(analysis=name,formula=deparse(formula),n_ALS=sum(meta$condition=='ALS'),
        n_control=sum(meta$condition=='Control'),n_genes=NA,status=paste('not run:',conditionMessage(e)))
      warning(name,': ',conditionMessage(e)); NULL
    })
    if (!is.null(z)) fits[[name]] <<- z
  }
}
add_fit('GSE215364_primary',counts,rna_meta,~condition,TRUE)
add_fit('GSE41414_primary',array,array_meta,~condition)
for (acc in c('GSE215364','GSE41414')) {
  meta <- if (acc=='GSE215364') rna_meta else array_meta
  mat <- if (acc=='GSE215364') counts else array
  isrna <- acc=='GSE215364'
  covars <- c('age','sex')
  use <- covars[vapply(covars,function(nm) !anyNA(meta[[nm]]) && length(unique(meta[[nm]]))>1,logical(1))]
  if (length(use)) add_fit(paste0(acc,'_covariate_adjusted'),mat,meta,
    reformulate(c(use,'condition')),isrna,TRUE)
  if (!anyNA(meta$muscle) && length(unique(meta$muscle))>1) add_fit(paste0(acc,'_muscle_adjusted'),mat,meta,
    ~muscle+condition,isrna,TRUE)
}
assert(!anyNA(array_meta$muscle),paste('GSE41414 muscle metadata missing. Copy the QC metadata override template',
  'to GSE41414/sample_metadata.csv, fill biopsy muscle from the GEO records and rerun. Primary DE is saved.'))
dm <- array_meta[tolower(array_meta$muscle)=='deltoid',,drop=FALSE]
assert(sum(dm$condition=='ALS')==4L && sum(dm$condition=='Control')==7L,
  'Expected deltoid sensitivity cohort: 4 ALS and 7 controls. Check muscle assignments.')
add_fit('GSE41414_deltoid_only',array,dm,~condition)
emit(do.call(rbind,model_log),'model_inventory',TRUE)

# 6. Frozen-module evaluation: complete coverage, fgsea, CAMERA and concordance.
summaries <- list(); gsea_all <- list(); cameras <- list(); permutations <- list()
evaluate <- function(name,obj) {
  de <- obj$de
  finite <- is.finite(de$stat) & is.finite(de$log2FC) & is.finite(de$p_value)
  bg <- de[finite,,drop=FALSE]
  coverage <- merge(module,de,by='symbol',all.x=TRUE,sort=FALSE)
  coverage$present_in_expression <- coverage$symbol %in% rownames(obj$expr)
  coverage$rank_evaluable <- is.finite(coverage$stat) & is.finite(coverage$log2FC) & is.finite(coverage$p_value)
  coverage$reason <- ifelse(!coverage$present_in_expression,'not annotated / not on platform / expression filtered',
    ifelse(!coverage$rank_evaluable,'DE statistic unavailable (e.g. Cook outlier)','evaluable'))
  coverage$aligned_log2FC <- coverage$log2FC*coverage$als_direction
  emit(coverage,paste0(name,'_module_genes'))
  assert(sum(coverage$rank_evaluable)>=MIN_GSEA_SIZE,
    paste(name,'has too few evaluable module genes; inspect annotation and filtering.'))
  ranks <- sort(setNames(bg$stat,bg$symbol),decreasing=TRUE)
  emit(data.frame(symbol=names(ranks),stat=unname(ranks)),paste0(name,'_ranks'))
  sizes <- vapply(sets,function(s) sum(s %in% names(ranks)),integer(1))
  valid <- sets[sizes>=MIN_GSEA_SIZE]
  if (length(valid)) {
    g <- as.data.frame(fgsea::fgseaMultilevel(pathways=valid,stats=ranks,minSize=MIN_GSEA_SIZE,
      maxSize=500,eps=0,nPermSimple=10000,nproc=1))
    g$analysis <- name; g$expected_direction <- ifelse(g$pathway=='SOD1_DOWN',-1,1)
    g$expected_NES <- sign(g$NES)==g$expected_direction
    gsea_all[[name]] <<- g
    emit(g,paste0(name,'_fgsea'))
    for (s in names(valid)) saveplot(fgsea::plotEnrichment(valid[[s]],ranks)+
      ggplot2::labs(title=paste(name,s),subtitle='Ranking: ALS versus control; positive = higher in ALS'),
      paste0(name,'_GSEA_',s))
  }
  emit(data.frame(component=names(sizes),n_frozen=vapply(sets,length,integer(1)),n_ranked=sizes,
    tested=sizes>=MIN_GSEA_SIZE),paste0(name,'_GSEA_coverage'),TRUE)
  ci <- lapply(sets,function(s) which(rownames(obj$expr) %in% s))
  ci <- ci[lengths(ci)>=MIN_GSEA_SIZE]
  if (length(ci)) {
    ca <- limma::camera(obj$camera_input,index=ci,design=obj$design,contrast=obj$coef,inter.gene.cor=NA)
    ca$component <- rownames(ca); ca$analysis <- name
    cameras[[name]] <<- ca; emit(ca,paste0(name,'_CAMERA'))
  }
  # Descriptive threshold counts: denominator includes neutral effects.
  for (component in c('ALL',names(sets))) for (threshold in LFC_THRESHOLDS) {
    q <- coverage[coverage$rank_evaluable & (component=='ALL' | coverage$component==component),,drop=FALSE]
    aligned <- q$aligned_log2FC
    k <- sum(aligned>0 & aligned>=threshold); opposite <- sum(aligned<0 & aligned<=-threshold)
    summaries[[paste(name,component,threshold)]] <<- data.frame(analysis=name,component=component,
      threshold=threshold,n_frozen=if(component=='ALL')174L else length(sets[[component]]),n_evaluable=nrow(q),
      n_concordant=k,n_discordant=opposite,n_neutral=nrow(q)-k-opposite,
      fraction_concordant=if(nrow(q))k/nrow(q) else NA_real_,
      n_concordant_FDR05=sum(aligned>0 & aligned>=threshold & q$FDR<0.05,na.rm=TRUE))
  }
  # Matched gene-set null: nonmodule background, 10 expression strata;
  # one unique random gene per module slot, preserving the up/down split.
  # Reports finite competitive statistics; this does not model patient sampling.
  bg <- bg[order(bg$symbol),,drop=FALSE]
  bg$bin <- pmin(10L,ceiling(rank(bg$mean_expression,ties.method='first')/nrow(bg)*10))
  q <- coverage[coverage$rank_evaluable,,drop=FALSE]
  q$bin <- bg$bin[match(q$symbol,bg$symbol)]
  pool <- split(which(!bg$symbol %in% module$symbol),bg$bin[!bg$symbol %in% module$symbol])
  need <- table(q$bin)
  assert(all(vapply(names(need),function(b) length(pool[[b]])>=need[[b]],logical(1))),
    paste(name,'insufficient nonmodule genes for expression-matched null.'))
  obs <- mean(q$stat*q$als_direction)
  obs_count <- mean(q$aligned_log2FC>=0.1 & q$aligned_log2FC>0)
  null <- numeric(N_PERM); null_count <- numeric(N_PERM)
  for (b in seq_len(N_PERM)) {
    idx <- integer(nrow(q))
    for (bin in names(need)) {
      slots <- which(q$bin==as.integer(bin)); pp <- pool[[bin]]
      idx[slots] <- pp[sample.int(length(pp),length(slots),replace=FALSE)]
    }
    null[b] <- mean(bg$stat[idx]*q$als_direction)
    a <- bg$log2FC[idx]*q$als_direction
    null_count[b] <- mean(a>=0.1 & a>0)
  }
  pm <- data.frame(analysis=name,n=nrow(q),n_perm=N_PERM,observed_aligned_stat=obs,
    null_mean=mean(null),null_q95=unname(quantile(null,0.95)),
    p_expected=(1+sum(null>=obs))/(N_PERM+1),p_opposite=(1+sum(null<=obs))/(N_PERM+1),
    observed_concordance_0.1=obs_count,null_concordance_mean=mean(null_count),
    p_concordance_expected=(1+sum(null_count>=obs_count))/(N_PERM+1))
  permutations[[name]] <<- pm
  emit(data.frame(iteration=seq_len(N_PERM),aligned_stat=null,concordance_0.1=null_count),paste0(name,'_matched_null'),TRUE)
  # PCA and descriptive signed module score; no gene-wise tests of patient means.
  expr <- obj$expr
  vary <- apply(expr,1,sd)>0
  pc <- prcomp(t(expr[vary,,drop=FALSE]),center=TRUE,scale.=FALSE)
  pd <- data.frame(sample=rownames(pc$x),PC1=pc$x[,1],PC2=pc$x[,2],condition=obj$meta$condition,muscle=obj$meta$muscle)
  emit(pd,paste0(name,'_PCA'),TRUE)
  saveplot(ggplot2::ggplot(pd,ggplot2::aes(PC1,PC2,colour=condition,label=sample))+
    ggplot2::geom_point(size=3)+ggplot2::geom_text(vjust=-0.6,size=3)+ggplot2::theme_bw()+
    ggplot2::labs(title=name),paste0(name,'_PCA'))
  ix <- match(module$symbol,rownames(expr)); good <- !is.na(ix)
  mm <- expr[ix[good],,drop=FALSE]; dirs <- module$als_direction[good]
  sdv <- apply(mm,1,sd); nz <- is.finite(sdv) & sdv>0
  zz <- sweep(sweep(mm[nz,,drop=FALSE],1,rowMeans(mm[nz,,drop=FALSE]),'-'),1,sdv[nz],'/')
  scores <- data.frame(sample=colnames(expr),condition=obj$meta$condition,
    score=colMeans(sweep(zz,1,dirs[nz],'*')),n_genes=sum(nz))
  emit(scores,paste0(name,'_sample_scores'))
  saveplot(ggplot2::ggplot(scores,ggplot2::aes(condition,score,colour=condition))+
    ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(position=ggplot2::position_jitter(width=0.08),size=3)+
    ggplot2::theme_bw()+ggplot2::labs(title=name,y='Signed within-cohort gene z-score (descriptive)'),paste0(name,'_module_score'))
}
for (nm in names(fits)) evaluate(nm,fits[[nm]])
summary <- do.call(rbind,summaries); emit(summary,'directional_summary')
g <- do.call(rbind,gsea_all)
if (!is.null(g)) {
  # Predefined primary family = 2 components x 2 full cohorts. Sensitivity q-values
  # are reported separately; they are not additional independent replications.
  primary <- g$analysis %in% c('GSE215364_primary','GSE41414_primary')
  g$FDR_four_primary_tests <- NA_real_
  g$FDR_four_primary_tests[primary] <- p.adjust(g$pval[primary],method='BH',n=4)
  emit(g,'fgsea_all_models')
}
if(length(cameras)) emit(do.call(rbind,cameras),'CAMERA_all_models')
pm <- do.call(rbind,permutations)
pm$FDR_expected_all_models <- p.adjust(pm$p_expected,'BH')
emit(pm,'expression_matched_competitive_tests')

# 7. Independent cohort effect agreement. P-values are exploratory across genes.
a <- fits$GSE215364_primary$de; b <- fits$GSE41414_primary$de
cross <- merge(a,b,by='symbol',suffixes=c('_RNA','_array'))
cross <- cross[is.finite(cross$log2FC_RNA) & is.finite(cross$log2FC_array),,drop=FALSE]
cross <- merge(cross,module,by='symbol',all.x=TRUE)
cross$module_gene <- !is.na(cross$als_direction)
cross$same_direction <- cross$log2FC_RNA*cross$log2FC_array>0
emit(cross,'cross_cohort_gene_effects')
cs <- lapply(c('all_shared_genes','shared_module_genes'),function(s) {
  z <- if(s=='shared_module_genes') cross[cross$module_gene,,drop=FALSE] else cross
  ct <- if(nrow(z)>=3 && sd(z$log2FC_RNA)>0 && sd(z$log2FC_array)>0)
    cor.test(z$log2FC_RNA,z$log2FC_array,method='spearman',exact=FALSE) else NULL
  data.frame(scope=s,n=nrow(z),rho=if(is.null(ct))NA else unname(ct$estimate),
    exploratory_p=if(is.null(ct))NA else ct$p.value,fraction_same_direction=mean(z$same_direction))
})
emit(do.call(rbind,cs),'cross_cohort_summary')
cm <- cross[cross$module_gene,,drop=FALSE]
if(nrow(cm)) saveplot(ggplot2::ggplot(cm,ggplot2::aes(log2FC_RNA,log2FC_array,colour=component))+
  ggplot2::geom_hline(yintercept=0,colour='grey70')+ggplot2::geom_vline(xintercept=0,colour='grey70')+
  ggplot2::geom_point(size=2)+ggplot2::theme_bw()+ggplot2::labs(x='GSE215364 ALS vs control log2FC',
  y='GSE41414 ALS vs control log2FC',title='Shared frozen-module genes'), 'cross_cohort_module_effects')
if(!is.null(g)) saveplot(ggplot2::ggplot(g,ggplot2::aes(NES,analysis,colour=pathway))+
  ggplot2::geom_vline(xintercept=0,colour='grey70')+ggplot2::geom_point(size=3)+
  ggplot2::theme_bw()+ggplot2::labs(title='Directional module enrichment',subtitle='Expected: SOD1_DOWN negative; SOD1_UP positive'),
  'module_NES_all_models',9,6)

# 8. Compact workbook plus complete CSVs and cached objects for figure scripts.
# Short unique worksheet names avoid Excel's 31-character limit.
wb <- openxlsx::createWorkbook()
sheet_index <- data.frame(sheet=sprintf('T%03d',seq_along(tables)),table=names(tables))
openxlsx::addWorksheet(wb,'INDEX'); openxlsx::writeData(wb,'INDEX',sheet_index)
for (i in seq_along(tables)) {
  nm <- sheet_index$sheet[i]; openxlsx::addWorksheet(wb,nm)
  openxlsx::writeData(wb,nm,tables[[i]],withFilter=TRUE)
  openxlsx::freezePane(wb,nm,firstRow=TRUE)
}
openxlsx::saveWorkbook(wb,file.path(OUTDIR,'human_ALS_module_validation_results.xlsx'),overwrite=TRUE)
saveRDS(list(module=module,fits=fits,summary=summary,gsea=g,cross=cross),
  file.path(CD,'human_ALS_validation_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(QD,'sessionInfo.txt'))
writeLines(c(
  'Independent frozen-module evaluation. Human data never redefine module.',
  'Primary: GSE215364 RNA-seq DESeq2 ~condition and GSE41414 RMA limma ~condition.',
  'GSE41414 deltoid only: 4 ALS vs 7 controls; normalization and probe choice retained from full cohort.',
  'RNA: unambiguous Ensembl/Entrez -> official human symbol; counts summed by symbol before fitting.',
  'Array: unambiguous probe -> symbol; largest mean-expression probe retained before DE.',
  'RNA filtering: edgeR filterByExpr with the corresponding model design.',
  'Ranking: finite DESeq2 Wald or limma moderated t statistics with finite raw P-values; no DE-significance preselection.',
  'fgsea two-sided P-values; inspect NES sign. Four primary tests receive a joint BH FDR.',
  'CAMERA estimates intergene correlation (inter.gene.cor=NA); RNA uses voom as a sensitivity.',
  'Matched random-set P-values are competitive/exploratory and do not account for all intergene correlation.',
  'Missing genes are not zeros and never count as discordant; coverage records all frozen genes.',
  'Module scores are descriptive within-cohort z-scores; magnitudes are not comparable across cohorts.',
  'Covariate/site models are sensitivities. Inspect QC metadata and model inventories before interpretation.',
  'GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215364',
  'GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41414',
  'Methods: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html',
  'Methods: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html'
),file.path(QD,'analysis_notes.txt'))
message('\nCOMPLETE: ',file.path(OUTDIR,'human_ALS_module_validation_results.xlsx'))
message('Share the workbook, QC/model_inventory.csv and QC/*sample_metadata.csv for interpretation.')
