#!/usr/bin/env Rscript
# Code author: Dr Eleni Christoforidou
# Copyright (c) 2026 Dr Eleni Christoforidou
# SPDX-License-Identifier: Apache-2.0
# See the repository LICENSE file for the full license terms.

# =============================================================================
# UNIFIED ACTUAL-SPACEFLIGHT ANALYSIS | ALS vs spaceflight manuscript
# Processes each NASA and JAXA comparison separately; counts are not pooled across missions.
# Computational integration, NOT an unloading/analogue analysis
#
# Independent predefined comparisons, NEVER pool counts between missions.
# NASA RR-1: OSD-101 left GST, OSD-419 right GST, OSD-104 SOL
# NASA RR-5: OSD-401 GST and QUAD
# NASA RR-23: OSD-770 SOL, OSD-665 EDL, OSD-666 QUAD
# NASA BION-M1: OSD-111 SOL and EDL (NORMALIZED microarray probesets required)
# NASA STS-108: OSD-21 flight GST vs contemporaneous flight ground controls
# NASA tissue chip: OSD-684 young/old human myobundles (GLDS-615)
# JAXA: MHU-1 MG vs onboard artificial 1g (plus ground sensitivity),
#       MHU-4/5 ~1/6g vs ground; soleus.
#
# FIXED ALS REFERENCE: six original SOD1G93A KEGG terms + three Hallmarks;
# and the original frozen 174-gene module (128 down, 46 up).
# The original nine KEGG pathways include three without exact MSigDB legacy
# counterparts: DO NOT replace those with approximately related terms.
#
# NASA RSEM aggregated expected counts -> edgeR TMM / limma-voom (fractional
# counts allowed; NOT tximport with effective-length offsets).
# Microarray -> existing background-corrected, RMA-normalized probeset values
# -> limma; NEVER use *_raw_intensities_probe* as normalized expression.
# JAXA count tables -> same voom framework, fitted separately within mission.
#
# All individual genes ranked genome-wide; gene sets are evaluated independently.
# Primary group comparisons BH-adjusted together for fixed-module enrichment.
# Same-mission muscles are NOT independent missions/replications.
# Missing datasets are written to INCLUSION.csv, NOT silently included.
# =============================================================================

# ============================= CONFIGURATION =================================
DATA_DIR <- getwd()
OUT_DIR <- file.path(DATA_DIR,"UNIFIED_ACTUAL_SPACEFLIGHT_ALS_results")
MODULE_FILE <- file.path(DATA_DIR,"reversible_ALS_unloading_results","tables",
                         "reversible_ALS_unloading_module_DIRECTIONAL.csv")
# Strict = stop at the end if any registered experiment was missing/failed,
# but preserve the results of successfully completed experiments.
STRICT_ALL_STUDIES <- FALSE
AUTO_INSTALL_PACKAGES <- TRUE
SEED <- 20261015L
ALPHA <- 0.05
MIN_MODULE_COVERAGE <- 0.60
MIN_GENES_FOR_GSEA <- 10L
# For studies with normalized arrays, limma uses log2 expression directly.
# For NASA RNA-seq, detection is based on the exact GLDS identifier below.
# Override any unusual local path by putting it in OVERRIDE below.
OVERRIDE <- list(
  # OSD104 = list(data="C:/.../GLDS-104_...csv", isa="C:/.../s_OSD-104.txt"),
  # MHU1 = list(data="C:/.../COUNTtable_for_MHU1.tsv", metadata="...tsv")
)
# Users may restrict to studies to reduce test runtime; BUT such a run must
# not be reported as the COMPLETE spaceflight analysis. Default = all studies.
RUN_ONLY <- NULL
# =============================================================================

set.seed(SEED)
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
for (d in c("QC","DE","GSEA","module","figures","concordance"))
  dir.create(file.path(OUT_DIR,d),showWarnings=FALSE,recursive=TRUE)
LOG <- file.path(OUT_DIR,"execution_log.txt")
cat("UNIFIED ACTUAL SPACEFLIGHT ANALYSIS\n",file=LOG)
msg <- function(...) {x <- paste0(...,collapse="");message(x)
  cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),x,"\n",file=LOG,append=TRUE)}
stopf <- function(...) stop(paste0(...,collapse=""),call.=FALSE)
wcsv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na="")
safe <- function(x) gsub("[^A-Za-z0-9_-]+","_",x)

bioc <- c("edgeR","limma","fgsea","org.Mm.eg.db","org.Hs.eg.db","AnnotationDbi")
cran <- c("data.table","msigdbr","babelgene","ggplot2","svglite","statmod")
missing <- function(x) x[!vapply(x,requireNamespace,logical(1),quietly=TRUE)]
if (AUTO_INSTALL_PACKAGES && length(c(missing(bioc),missing(cran)))) {
  if (length(missing(cran))) install.packages(missing(cran),repos="https://cloud.r-project.org")
  if (length(missing(bioc))) {
    if (!requireNamespace("BiocManager",quietly=TRUE))
      install.packages("BiocManager",repos="https://cloud.r-project.org")
    BiocManager::install(missing(bioc),ask=FALSE,update=FALSE)
  }
}
if (length(c(missing(bioc),missing(cran)))) stopf("Missing packages: ",
  paste(c(missing(bioc),missing(cran)),collapse=", "))

# -------------------- Registry of actual-flight specimens -------------------
# GLDS and OSD identifiers differ in RR-23 (604/605/673) and OSD-684 (615).
registry <- list(
  OSD101=list(platform="rnaseq",glds="101",isa="OSD-101_metadata_GLDS-101-ISA",mission="RR-1",species="mouse",target="GASTROCNEMIUS"),
  OSD419=list(platform="rnaseq",glds="419",isa="OSD-419_metadata_GLDS-419-ISA",mission="RR-1",species="mouse",target="GASTROCNEMIUS"),
  OSD104=list(platform="rnaseq",glds="104",isa="OSD-104_metadata_OSD-104-ISA",mission="RR-1",species="mouse",target="SOLEUS"),
  OSD401=list(platform="rnaseq",glds="401",isa="OSD-401_metadata_GLDS-401-ISA",mission="RR-5",species="mouse",target=c("GASTROCNEMIUS","QUADRICEPS")),
  OSD770=list(platform="rnaseq",glds="673",isa="OSD-770_metadata_OSD-770-ISA",mission="RR-23",species="mouse",target="SOLEUS"),
  OSD665=list(platform="rnaseq",glds="604",isa="OSD-665_metadata_OSD-665-ISA",mission="RR-23",species="mouse",target="EDL"),
  OSD666=list(platform="rnaseq",glds="605",isa="OSD-666_metadata_OSD-666-ISA",mission="RR-23",species="mouse",target="QUADRICEPS"),
  OSD111=list(platform="array",glds="111",isa="OSD-111_metadata_OSD-111-ISA",mission="BION-M1",species="mouse",target=c("SOLEUS","EDL")),
  OSD21=list(platform="array",glds="21",isa="OSD-21_metadata_OSD-21-ISA",mission="STS-108",species="mouse",target="GASTROCNEMIUS"),
  OSD684=list(platform="rnaseq",glds="615",isa="OSD-684_metadata_GLDS-615-ISA",mission="CRS-21",species="human",target="MYOBUNDLE"),
  MHU1=list(platform="jaxa",mission="MHU-1",species="mouse",target="SOLEUS"),
  MHU4=list(platform="jaxa",mission="MHU-4",species="mouse",target="SOLEUS"),
  MHU5=list(platform="jaxa",mission="MHU-5",species="mouse",target="SOLEUS")
)
if (!is.null(RUN_ONLY)) {
  if (!all(RUN_ONLY %in% names(registry))) stopf("Unknown studies in RUN_ONLY")
  registry <- registry[RUN_ONLY]
}
files <- list.files(DATA_DIR,recursive=TRUE,full.names=TRUE,all.files=FALSE)
files <- files[file.exists(files) & !dir.exists(files)]
files <- files[!startsWith(normalizePath(files,winslash="/"),
  paste0(normalizePath(OUT_DIR,winslash="/"),"/"))]
# Prefer exact basename, then (1) Windows duplicate, never another GLDS number.
find_named <- function(study,expected_name,target="data") {
  v <- OVERRIDE[[study]][[target]]
  if (!is.null(v)) {
    if (!file.exists(v)) stopf(study," explicitly configured ",target," missing: ",v)
    return(normalizePath(v))
  }
  hits <- files[base::basename(files)==expected_name]
  if (!length(hits)) {
    # Accept Windows duplicate-download suffix before either CSV or TSV extension.
    alt <- sub("(\\.[^.]+)$", "(1)\\1", expected_name)
    hits <- files[base::basename(files)==alt]
  }
  if (!length(hits)) return(NA_character_)
  if (length(hits)>1L) stopf(study," ambiguous ",target,": ",paste(hits,collapse=" | "),
    "\nSet OVERRIDE[['",study,"']]$",target," to the exact file.")
  normalizePath(hits[[1]])
}
find_isa <- function(study,cfg) {
  v <- OVERRIDE[[study]][["isa"]]
  if (!is.null(v)) {if (!file.exists(v)) stopf("Missing ISA override ",v); return(v)}
  f <- paste0("s_",sub("^OSD","OSD-",study),".txt")
  # Prefer ISA folder, including changed legacy prefix conventions.
  candidate <- file.path(DATA_DIR,cfg$isa,f)
  if (file.exists(candidate)) return(candidate)
  h <- files[basename(files)==f]
  if (length(h)>1L) {
    matchfolder <- grepl(cfg$isa,h,fixed=TRUE)
    if (sum(matchfolder)==1) h <- h[matchfolder]
  }
  if (length(h)==1L) return(h)
  # Can read zipped ISA archive, using its s_*.txt entry.
  z <- files[grepl(paste0("OSD-",sub("^OSD","",study),"_metadata.*ISA\\.zip$"),
                   basename(files),ignore.case=TRUE)]
  if (length(z)==1L) return(z)
  NA_character_
}
read_isa <- function(path) {
  if (is.na(path)) stopf("Missing ISA sample metadata")
  if (grepl("\\.zip$",path,ignore.case=TRUE)) {
    members <- utils::unzip(path,list=TRUE)$Name
    s <- members[grepl("(^|/)s_[^/]+\\.txt$",members)]
    if (length(s)!=1) stopf("Expected one s_*.txt inside ISA ZIP: ",path)
    conn <- unz(path,s[[1]],open="r");on.exit(close(conn))
    x <- utils::read.delim(conn,check.names=FALSE,quote="",fill=TRUE,
                           stringsAsFactors=FALSE,comment.char="")
  } else x <- utils::read.delim(path,check.names=FALSE,quote="",fill=TRUE,
                             stringsAsFactors=FALSE,comment.char="")
  names(x) <- make.unique(names(x))
  x
}
getcol <- function(x,pattern,required=TRUE) {
  n <- grep(pattern,names(x),value=TRUE,ignore.case=TRUE)
  if (!length(n)) {if (required) stopf("Column not found for pattern '",pattern,
                            "'. Available: ",paste(head(names(x),35),collapse=", "))
                  return(NULL)}
  n[[1]]
}
key <- function(x) {x <- sub("^GLDS[-_][0-9]+_rna[-_]seq_","",as.character(x),ignore.case=TRUE)
  tolower(gsub("[^a-zA-Z0-9]","",trimws(x)))}
metadata_assign <- function(isa,sample_names,study,cfg) {
  sn <- getcol(isa,"^Sample Name$")
  sc <- getcol(isa,"^Source Name$",required=FALSE)
  m <- match(key(sample_names),key(isa[[sn]]))
  if (anyNA(m)) stopf(study,": ",sum(is.na(m))," samples unmatched to ISA; examples: ",
                      paste(head(sample_names[is.na(m)],6),collapse=", "),
                      ". Confirm the exact study metadata and assay sample names.")
  x <- isa[m,,drop=FALSE]
  sf <- grep("^Factor Value\\[Spaceflight",names(x),value=TRUE,ignore.case=TRUE)
  if (!length(sf)) sf <- grep("^Factor Value\\[.*(spaceflight|flight)",names(x),value=TRUE,ignore.case=TRUE)
  if (!length(sf)) stopf(study," lacks ISA Factor Value[Spaceflight] annotation. No guessing group from filename.")
  raw <- toupper(trimws(as.character(x[[sf[[1]]]])))
  group <- rep(NA_character_,length(raw))
  group[grepl("SPACE[ -]?FLIGHT|^FLIGHT$",raw)] <- "FLT"
  group[grepl("GROUND CONTROL|^GROUND$",raw)] <- "GC"
  group[grepl("VIVARIUM",raw)] <- "VC"
  housing <- getcol(x,"^Parameter Value\\[(habitat|housing)",required=FALSE)
  if (!is.null(housing)) {
    h <- toupper(as.character(x[[housing]]))
    group[group=="GC" & grepl("HABITAT|RODENT FLIGHT HARDWARE|TRANSPORTER",h)] <- "HGC"
    group[group=="GC" & grepl("VIVARIUM|STANDARD VIVARIUM",h)] <- "VC"
  }
  tissue_col <- getcol(x,"^Factor Value\\[Tissue",required=FALSE)
  if (is.null(tissue_col)) tissue_col <- getcol(x,"^Characteristics\\[(Material Type|Tissue|Organ)",required=FALSE)
  t <- if (!is.null(tissue_col)) toupper(as.character(x[[tissue_col]])) else rep("",nrow(x))
  t <- paste(t,toupper(as.character(x[[sn]])))
  tissue <- rep(NA_character_,length(t))
  tissue[grepl("SOLEUS|(^|[^A-Z])SLS([^A-Z]|$)",t)] <- "SOLEUS"
  tissue[grepl("GASTROCNEMIUS|(^|[^A-Z])GST([^A-Z]|$)",t)] <- "GASTROCNEMIUS"
  tissue[grepl("QUADRICEPS|QUAD",t)] <- "QUADRICEPS"
  tissue[grepl("EXTENSOR DIGITORUM|(^|[^A-Z])EDL([^A-Z]|$)",t)] <- "EDL"
  if (study=="OSD684") tissue[] <- "MYOBUNDLE"
  if (length(cfg$target)==1L && is.na(tissue[1]) &&
      !any(!is.na(tissue))) tissue[] <- cfg$target
  age_col <- getcol(x,"^Factor Value\\[Donor Age",required=FALSE)
  if (is.null(age_col)) age_col <- getcol(x,"^Characteristics\\[Donor Age",required=FALSE)
  age_group <- rep("ALL",nrow(x))
  if (study=="OSD684") {
    if (is.null(age_col)) stopf("OSD684 missing donor-age metadata. Cannot stratify young/old.")
    a <- toupper(as.character(x[[age_col]]))
    age_group <- ifelse(grepl("21.?40|YOUNG|18.?40",a),"YOUNG",
                        ifelse(grepl("OLDER|OLD|60|65|70",a),"OLD",NA_character_))
    if (anyNA(age_group)) stopf("OSD684 donor age mapping unresolved: ",paste(unique(a),collapse=" | "))
  }
  sex <- getcol(x,"^Characteristics\\[Sex",required=FALSE)
  age <- getcol(x,"^Characteristics\\[Age at Launch",required=FALSE)
  data.frame(sample=sample_names,group=group,tissue=tissue,age_group=age_group,
             source=if (!is.null(sc)) as.character(x[[sc]]) else as.character(x[[sn]]),
             sex=if (is.null(sex)) NA_character_ else as.character(x[[sex]]),
             age_at_launch=if (is.null(age)) NA_character_ else as.character(x[[age]]),
             isa_sample=as.character(x[[sn]]),stringsAsFactors=FALSE)
}

# JAXA metadata reader: map sample IDs and gravity using actual column content;
# if ambiguous, STOP that study with diagnostics instead of fabricating groups.
read_jaxa_metadata <- function(path,cols,study) {
  z <- data.table::fread(path,data.table=FALSE,check.names=FALSE)
  possible <- names(z)[grepl("sample|name|id",names(z),ignore.case=TRUE)]
  if (!length(possible)) possible <- names(z)
  matches <- sapply(possible,function(n) sum(key(cols) %in% key(z[[n]])))
  if (max(matches)==0) stopf(study,": cannot match sample IDs to ibSLS metadata.\nMetadata columns: ",
      paste(names(z),collapse=", "),"\nData columns: ",paste(head(cols,12),collapse=", "))
  sc <- possible[which.max(matches)]
  idx <- match(key(cols),key(z[[sc]]))
  if (anyNA(idx)) stopf(study,": not all count samples matched metadata; ",
                        paste(cols[is.na(idx)],collapse=", "))
  z <- z[idx,,drop=FALSE]
  # ibSLS GroupID is mission-prefixed, e.g. MHU01_MG, MHU04_PG.
  # Match the EXACT registered mission and remove its prefix before modelling.
  # Never allow the metadata from one MHU mission to be applied to another.
  if (!("GroupID" %in% names(z))) stopf(study,": ibSLS GroupID column missing. Available: ",
                                        paste(names(z),collapse=", "))
  expected_prefix <- switch(study, MHU1="MHU01", MHU4="MHU04", MHU5="MHU05")
  if (is.null(expected_prefix)) stopf(study,": no explicit JAXA group schema")
  full_group <- toupper(trimws(as.character(z$GroupID)))
  if (any(!startsWith(full_group,paste0(expected_prefix,"_"))))
    stopf(study,": metadata includes a different mission's GroupID; observed: ",
          paste(unique(full_group),collapse=", "))
  grp <- substring(full_group,nchar(expected_prefix)+2L)
  allowed <- if (study=="MHU1") c("MG","AG","GC") else c("PG","GC")
  if (anyNA(grp) || any(!grp %in% allowed))
    stopf(study,": invalid JAXA gravity group: ",paste(unique(full_group),collapse=" | "),
          ". Expected mission ",expected_prefix," and groups ",paste(allowed,collapse=", "))
  msg(study,": ibSLS GroupID remapped: ",
      paste(names(table(grp)),as.integer(table(grp)),sep="=",collapse="; "))
  tissue <- rep("SOLEUS",nrow(z))
  tcols <- grep("(tissue|muscle|organ)",names(z),value=TRUE,ignore.case=TRUE)
  if (length(tcols)) {
    tt <- toupper(as.character(z[[tcols[[1]]]]))
    if (any(grepl("SOLEUS",tt))) tissue <- ifelse(grepl("SOLEUS",tt),"SOLEUS","OTHER")
  }
  data.frame(sample=cols,group=grp,tissue=tissue,age_group="ALL",
             source=cols,sex=NA_character_,age_at_launch=NA_character_,
             isa_sample=cols,stringsAsFactors=FALSE)
}

read_counts <- function(path) {
  x <- data.table::fread(path,data.table=FALSE,check.names=FALSE)
  ids <- sub("\\.[0-9]+$","",as.character(x[[1]]))
  cols <- names(x)[-1L]
  ok <- vapply(x[-1L],function(v) is.numeric(v) && all(is.finite(v)) &&
                   all(v>=0),logical(1))
  cols <- cols[ok]
  if (length(cols)<4L) stopf("Counts expected >=4 sample columns: ",path,
                             "; columns ",paste(names(x),collapse=", "))
  if (anyDuplicated(cols)) stopf("Duplicate sample names in ",path)
  mat <- as.matrix(x[,cols,drop=FALSE]);storage.mode(mat) <- "double"
  if (anyNA(ids) || any(!nzchar(ids))) stopf("Missing gene IDs in ",path)
  if (anyDuplicated(ids)) mat <- rowsum(mat,group=ids,reorder=FALSE)
  else rownames(mat) <- ids
  list(mat=mat,kind="counts",annotations=NULL)
}
read_normalized_array <- function(path) {
  if (grepl("raw_intensities|normalized_intensities_probe",basename(path),ignore.case=TRUE))
    stopf("DO NOT fit limma to raw/unsummarized probe intensities: ",path,
      "\nDownload GLDS-111_array_normalized_expression_probeset_GLmicroarray.csv")
  x <- data.table::fread(path,data.table=FALSE,check.names=FALSE)
  id_col <- getcol(x,"^(ProbesetID|probeset_id|probeset|ID_REF)$",required=FALSE)
  symbol_col <- getcol(x,"^SYMBOL$",required=FALSE)
  ens_col <- getcol(x,"^ENSEMBL$",required=FALSE)
  if (is.null(id_col) && is.null(symbol_col) && is.null(ens_col))
    stopf("Microarray matrix lacks ProbesetID/SYMBOL/ENSEMBL annotation: ",path)
  n <- names(x)
  sam <- n[vapply(x,is.numeric,logical(1))]
  sam <- sam[grepl("(FLT|_GC_|_VIV_|CONTROL|GSM|RR|Mmus)",sam,ignore.case=TRUE)]
  if (length(sam)<4) stopf("Cannot identify normalized array sample columns: ",path)
  m <- as.matrix(x[,sam,drop=FALSE]);storage.mode(m) <- "double"
  if (any(!is.finite(m))) stopf("Non-finite normalized expression values in ",path)
  # NASA processed probeset values are already log2 expression.
  if (max(m,na.rm=TRUE)>100) stopf("Expression exceeds typical log2 range; ",
                                  "check that this is NORMALIZED probeset matrix, not raw intensities.")
  ids <- if (!is.null(id_col)) as.character(x[[id_col]]) else paste0("probeset_",seq_len(nrow(x)))
  if (anyDuplicated(ids)) ids <- make.unique(ids)
  rownames(m) <- ids
  anno <- data.frame(gene_id=ids,
    mouse_symbol=if (is.null(symbol_col)) NA_character_ else as.character(x[[symbol_col]]),
    ensembl=if (is.null(ens_col)) NA_character_ else as.character(x[[ens_col]]),
    stringsAsFactors=FALSE)
  # Keep uniquely mapped Ensembl rows only if NASA provides mapping count.
  map_col <- getcol(x,"^count_ENSEMBL_mappings$",required=FALSE)
  if (!is.null(map_col)) {
    bad <- !is.na(x[[map_col]]) & as.numeric(x[[map_col]]) != 1
    anno$mouse_symbol[bad] <- NA_character_;anno$ensembl[bad] <- NA_character_
  }
  list(mat=m,kind="array",annotations=anno)
}

pool_technical <- function(m,meta,study) {
  # ISA Source Name repeated only for technical sequencing libraries.
  # Same animal in different muscles remains DIFFERENT measurements.
  k <- paste(meta$source,meta$tissue,meta$group,meta$age_group,sep="||")
  grp <- split(seq_len(nrow(meta)),factor(k,levels=unique(k)))
  if (any(lengths(grp)>1L) && any(vapply(grp,function(i) {
    length(i)>1 && !all(grepl("techrep",meta$sample[i],ignore.case=TRUE))
  },logical(1)))) stopf(study,": duplicated biological source lacking techrep labels. Inspect ISA.")
  m2 <- vapply(grp,function(i) rowSums(m[,i,drop=FALSE]),numeric(nrow(m)))
  if (is.vector(m2)) m2 <- matrix(m2,ncol=1)
  rownames(m2) <- rownames(m)
  new <- meta[vapply(grp,function(i) i[[1L]],integer(1)),,drop=FALSE]
  new$sample <- sub("_techrep[0-9]+$","",new$sample,ignore.case=TRUE)
  if (anyDuplicated(new$sample)) stopf(study,": pooled IDs not unique")
  colnames(m2) <- new$sample
  wcsv(data.frame(assay_column=meta$sample,pooled_sample=new$sample[match(k,names(grp))],
                   biological_source=meta$source,tissue=meta$tissue,group=meta$group),
       file.path(OUT_DIR,"QC",paste0(study,"_technical_replicate_map.csv")))
  list(mat=m2,meta=new)
}

annotate_genes <- function(row_ids,annotations,species) {
  core <- sub("\\.[0-9]+$","",row_ids)
  ans <- data.frame(gene_id=row_ids,mouse_symbol=NA_character_,human_symbol=NA_character_,
                    stringsAsFactors=FALSE)
  if (species=="mouse") {
    if (!is.null(annotations)) {
      ans$mouse_symbol <- annotations$mouse_symbol[match(row_ids,annotations$gene_id)]
      extra <- annotations$ensembl[match(row_ids,annotations$gene_id)]
      ind <- which(is.na(ans$mouse_symbol) & !is.na(extra) & grepl("^ENSMUSG",extra))
      if (length(ind)) {
        e <- sub("\\.[0-9]+$","",extra[ind])
        ans$mouse_symbol[ind] <- as.character(AnnotationDbi::mapIds(
            org.Mm.eg.db::org.Mm.eg.db,keys=e,keytype="ENSEMBL",
            column="SYMBOL",multiVals="first"))
      }
    } else {
      yes <- grepl("^ENSMUSG[0-9]+$",core)
      if (any(yes)) ans$mouse_symbol[yes] <- as.character(AnnotationDbi::mapIds(
         org.Mm.eg.db::org.Mm.eg.db,keys=unique(core[yes]),keytype="ENSEMBL",
         column="SYMBOL",multiVals="first")[core[yes]])
      ans$mouse_symbol[!yes] <- core[!yes]
    }
    names_q <- unique(ans$mouse_symbol[!is.na(ans$mouse_symbol) & nzchar(ans$mouse_symbol)])
    if (length(names_q)<1000) stopf("Poor mouse gene annotation: ",length(names_q)," genes")
    mp <- babelgene::orthologs(genes=names_q,species="mouse",human=FALSE,
                              min_support=3,top=TRUE)
    if (!all(c("symbol","human_symbol","support_n") %in% names(mp)))
      stopf("Unexpected babelgene ortholog columns: ",paste(names(mp),collapse=", "))
    mp <- mp[order(-mp$support_n,mp$symbol),,drop=FALSE]
    mp <- mp[!duplicated(mp$symbol),,drop=FALSE]
    mp <- mp[!duplicated(toupper(mp$human_symbol)),,drop=FALSE]
    ans$human_symbol <- toupper(mp$human_symbol[match(ans$mouse_symbol,mp$symbol)])
  } else {
    yes <- grepl("^ENSG[0-9]+$",core)
    if (any(yes)) ans$human_symbol[yes] <- as.character(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,keys=unique(core[yes]),keytype="ENSEMBL",
      column="SYMBOL",multiVals="first")[core[yes]])
    ans$human_symbol[!yes] <- core[!yes]
  }
  ans$human_symbol <- toupper(ans$human_symbol)
  ans$human_symbol[!nzchar(ans$human_symbol)] <- NA_character_
  ans
}

load_module <- function(f) {
  if (!file.exists(f)) stopf("Frozen reference missing: ",f,
    "\nCopy original reversible_ALS_unloading_module_DIRECTIONAL.csv or set MODULE_FILE.")
  d <- data.table::fread(f,data.table=FALSE)
  clean <- tolower(gsub("[^A-Za-z0-9]","",names(d)))
  i <- match("symbol",clean)
  if (is.na(i)) i <- match("genesymbol",clean)
  j <- match("alsdirection",clean)
  if (is.na(j)) j <- match("sod1direction",clean)
  if (is.na(i)||is.na(j)) stopf("Module needs symbol and als_direction; found ",paste(names(d),collapse=", "))
  x <- toupper(as.character(d[[i]])); z <- toupper(trimws(as.character(d[[j]])))
  sign <- ifelse(z %in% c("-1","DOWN","DOWNREGULATED"),-1L,
                 ifelse(z %in% c("+1","1","UP","UPREGULATED"),1L,NA_integer_))
  if (length(x)!=174L || anyDuplicated(x) || anyNA(sign) ||
      sum(sign==-1L)!=128L || sum(sign==1L)!=46L) stopf(
      "Frozen module invalid. Must be 174 unique genes: 128 DOWN, 46 UP.")
  data.frame(human_symbol=x,als_sign=sign,stringsAsFactors=FALSE)
}
module <- load_module(MODULE_FILE)
wcsv(module,file.path(OUT_DIR,"QC","FROZEN_module_174.csv"))

make_gene_sets <- function(species) {
  name <- if (species=="mouse") "Mus musculus" else "Homo sapiens"
  k <- msigdbr::msigdbr(species=name,collection="C2",subcollection="CP:KEGG_LEGACY")
  h <- msigdbr::msigdbr(species=name,collection="H")
  paths <- split(c(as.character(k$gene_symbol),as.character(h$gene_symbol)),
                 c(as.character(k$gs_name),as.character(h$gs_name)))
  paths <- lapply(paths,unique)
  ref <- c("KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
           "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",
           "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
           "KEGG_INSULIN_SIGNALING_PATHWAY",
           "KEGG_PURINE_METABOLISM",
           "KEGG_STARCH_AND_SUCROSE_METABOLISM")
  hall <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
            "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_GLYCOLYSIS")
  if (!all(c(ref,hall) %in% names(paths))) stopf("Missing exact original KEGG/Hallmark gene sets for ",species,
         ": ",paste(setdiff(c(ref,hall),names(paths)),collapse=", "),
         ". Do not silently substitute another collection.")
  list(all=paths,ref=ref,hall=hall,selected=c(ref,hall))
}
sets_mouse <- make_gene_sets("mouse")
sets_human <- make_gene_sets("human")

# Fit all contrasts for each study/tissue using its own appropriately processed
# expression data. Design includes no cross-mission or cross-muscle pooling.
fit_contrast <- function(mat,meta,numerator,denominator,kind) {
  keep_s <- meta$group %in% c(numerator,denominator)
  md <- meta[keep_s,,drop=FALSE]
  x <- mat[,md$sample,drop=FALSE]
  tab <- table(md$group)
  if (!all(c(numerator,denominator) %in% names(tab)) || any(tab[c(numerator,denominator)]<2L))
    stopf("Need >=2 independent biological samples/group, got: ",
          paste(names(tab),tab,collapse="; "))
  f <- factor(md$group,levels=c(denominator,numerator))
  design <- stats::model.matrix(~0+f);colnames(design) <- c("CONTROL","EXPOSURE")
  if (kind=="counts") {
    y <- edgeR::DGEList(counts=x)
    expressed <- edgeR::filterByExpr(y,design=design)
    if (sum(expressed)<1000) stopf("Too few expressed genes after filterByExpr")
    y <- edgeR::calcNormFactors(y[expressed,,keep.lib.sizes=FALSE],method="TMM")
    v <- limma::voom(y,design,plot=FALSE)
    es <- v$E
    fit <- limma::lmFit(v,design)
  } else {
    if (nrow(x)<1000) stopf("Too few microarray probesets")
    es <- x
    fit <- limma::lmFit(x,design)
  }
  cm <- limma::makeContrasts(EXPOSURE-CONTROL,levels=design)
  fit <- limma::eBayes(limma::contrasts.fit(fit,cm),robust=TRUE)
  out <- limma::topTable(fit,coef=1,sort.by="none",number=Inf)
  out$gene_id <- rownames(out)
  names(out)[names(out)=="logFC"] <- "log2FoldChange"
  names(out)[names(out)=="t"] <- "stat"
  names(out)[names(out)=="P.Value"] <- "pvalue"
  names(out)[names(out)=="adj.P.Val"] <- "padj"
  names(out)[names(out)=="AveExpr"] <- "avgExpr"
  list(de=out,expression=es,sample_meta=md)
}

# De-duplicate by mean expression measured WITHOUT using differential-expression
# direction/statistic. This avoids picking the probe/gene with largest |t|.
unique_by_gene <- function(df,symcol) {
  x <- df[!is.na(df[[symcol]]) & nzchar(df[[symcol]]) & is.finite(df$stat),,drop=FALSE]
  x <- x[order(-x$avgExpr,x$gene_id),,drop=FALSE]
  x[!duplicated(x[[symcol]]),,drop=FALSE]
}
fgsea_once <- function(paths,stats) {
  paths <- paths[lengths(lapply(paths,intersect,names(stats)))>=MIN_GENES_FOR_GSEA]
  if (!length(paths)) stopf("No gene sets have sufficient coverage")
  r <- as.data.frame(fgsea::fgseaMultilevel(pathways=paths,stats=stats,
                       minSize=MIN_GENES_FOR_GSEA,maxSize=500,eps=0))
  if ("leadingEdge" %in% names(r)) r$leadingEdge <- vapply(r$leadingEdge,
                                                   paste,collapse=";",FUN.VALUE=character(1))
  r
}

all_de <- list();all_pathways <- list();all_modules <- list();all_gene <- list()
all_summary <- list();all_samples <- list();all_scoretests <- list();all_qc <- list()
inclusion <- list(); all_model_symbols <- list(); fitted <- list()

append_result <- function(study,mission,tissue,cohort,cmp,primary,fit,ann,species) {
  id <- paste(study,tissue,cohort,cmp,sep="__")
  d <- fit$de
  i <- match(d$gene_id,ann$gene_id)
  d$mouse_symbol <- ann$mouse_symbol[i]
  d$human_symbol <- ann$human_symbol[i]
  dm <- if (species=="mouse") unique_by_gene(d,"mouse_symbol") else NULL
  dh <- unique_by_gene(d,"human_symbol")
  ranks_h <- setNames(dh$stat,dh$human_symbol)
  ranks_h <- sort(ranks_h,decreasing=TRUE)
  ranks_native <- if (species=="mouse") setNames(dm$stat,dm$mouse_symbol) else ranks_h
  ranks_native <- sort(ranks_native,decreasing=TRUE)
  if (length(ranks_h)<1000 || length(ranks_native)<1000)
    stopf(id,": fewer than 1000 unique mapped genes")
  chosen <- if (species=="mouse") sets_mouse else sets_human
  enr <- fgsea_once(chosen$all,ranks_native)
  enr$dataset<-study;enr$mission<-mission;enr$tissue<-tissue;enr$cohort<-cohort
  enr$contrast<-cmp;enr$primary<-primary;enr$species<-species
  selected <- enr[enr$pathway %in% chosen$selected,,drop=FALSE]
  if (nrow(selected)!=9) stopf(id,": only ",nrow(selected),
                               "/9 planned metabolic pathways tested")
  component <- list(ALS_DOWN=module$human_symbol[module$als_sign==-1L],
                    ALS_UP=module$human_symbol[module$als_sign==1L])
  m <- fgsea_once(component,ranks_h)
  if (nrow(m)!=2) stopf(id,": both frozen module components required")
  m$expected_NES_sign <- ifelse(m$pathway=="ALS_DOWN",-1L,1L)
  m$concordant_with_ALS <- sign(m$NES)==m$expected_NES_sign
  m$dataset<-study;m$mission<-mission;m$tissue<-tissue;m$cohort<-cohort
  m$contrast<-cmp;m$primary<-primary;m$species<-species
  idx <- match(module$human_symbol,dh$human_symbol)
  g <- cbind(module,data.frame(log2FoldChange=dh$log2FoldChange[idx],
         stat=dh$stat[idx],pvalue=dh$pvalue[idx],padj=dh$padj[idx]))
  g$direction_match <- ifelse(is.na(g$log2FoldChange),NA,
          sign(g$log2FoldChange)==g$als_sign)
  g$direction_match_LFC_ge_0_10 <- !is.na(g$direction_match) &
           g$direction_match & abs(g$log2FoldChange)>=0.10
  g$dataset<-study;g$mission<-mission;g$tissue<-tissue;g$cohort<-cohort
  g$contrast<-cmp;g$primary<-primary
  n_eval <- sum(!is.na(g$direction_match))
  if (n_eval<174L*MIN_MODULE_COVERAGE) stopf(id,": frozen module coverage only ",n_eval,"/174")
  summ <- data.frame(dataset=study,mission=mission,tissue=tissue,cohort=cohort,
    contrast=cmp,primary=primary,species=species,
    n_numerator=sum(fit$sample_meta$group==strsplit(cmp,"_vs_",fixed=TRUE)[[1]][1]),
    n_denominator=sum(fit$sample_meta$group==strsplit(cmp,"_vs_",fixed=TRUE)[[1]][2]),
    n_gene_tested=nrow(d),n_DEGs_FDR05=sum(d$padj<ALPHA,na.rm=TRUE),
    n_module_evaluable=n_eval,
    n_module_concordant=sum(g$direction_match,na.rm=TRUE),
    pct_module_concordant=100*mean(g$direction_match,na.rm=TRUE),
    n_concordant_absLFC01=sum(g$direction_match_LFC_ge_0_10),
    stringsAsFactors=FALSE)
  # Sample-level independent validation, using a fixed gene set and per-gene
  # within-study standardized expression; biological samples are the replicates.
  # Only fit EXPOSE vs CONTROL samples; z transform within all samples of that
  # tissue/cohort. A score is NOT used to select module genes.
  es <- fit$expression
  gene_to_expr <- data.frame(gene_id=rownames(es),avgExpr=rowMeans(es))
  ann_e <- merge(ann,gene_to_expr,by="gene_id",all=FALSE)
  ann_e <- ann_e[!is.na(ann_e$human_symbol) &
      ann_e$human_symbol %in% module$human_symbol,,drop=FALSE]
  ann_e <- ann_e[order(-ann_e$avgExpr,ann_e$gene_id),,drop=FALSE]
  ann_e <- ann_e[!duplicated(ann_e$human_symbol),,drop=FALSE]
  ix <- match(ann_e$gene_id,rownames(es))
  zz <- es[ix,,drop=FALSE]
  sdgene <- apply(zz,1,stats::sd)
  take <- is.finite(sdgene) & sdgene>0
  zz <- zz[take,,drop=FALSE]
  ss <- module$als_sign[match(ann_e$human_symbol[take],module$human_symbol)]
  z <- t(scale(t(zz)))
  score <- colMeans(z*ss)
  sc <- data.frame(dataset=study,mission=mission,tissue=tissue,cohort=cohort,
                   contrast=cmp,sample=names(score),
                   group=fit$sample_meta$group[match(names(score),fit$sample_meta$sample)],
                   ALS_oriented_module_score=as.numeric(score),
                   n_genes_used=nrow(z),stringsAsFactors=FALSE)
  pos <- sc$ALS_oriented_module_score[sc$group==strsplit(cmp,"_vs_",fixed=TRUE)[[1]][1]]
  neg <- sc$ALS_oriented_module_score[sc$group==strsplit(cmp,"_vs_",fixed=TRUE)[[1]][2]]
  pv <- tryCatch(stats::t.test(pos,neg)$p.value,error=function(e) NA_real_)
  test <- data.frame(dataset=study,mission=mission,tissue=tissue,cohort=cohort,
                     contrast=cmp,primary=primary,n_genes=nrow(z),
                     mean_FLT_minus_control=mean(pos)-mean(neg),pvalue=pv)
  # PCA quality check: computed only within this tissue, cohort and contrast;
  # no cross-mission sample-level batch mixing.
  if (ncol(es)>=4 && nrow(es)>=500) {
    vars <- apply(es,1,stats::var)
    take_pca <- head(order(vars,decreasing=TRUE),min(500L,nrow(es)))
    pca <- stats::prcomp(t(es[take_pca,,drop=FALSE]),center=TRUE,scale.=FALSE)
    pc <- as.data.frame(pca$x[,seq_len(min(2L,ncol(pca$x))),drop=FALSE])
    if (ncol(pc)==2L) {
      pc$sample <- rownames(pc)
      pc$group <- fit$sample_meta$group[match(pc$sample,fit$sample_meta$sample)]
      wcsv(pc,file.path(OUT_DIR,"QC",paste0(safe(id),"_PCA.csv")))
      plot_pc <- ggplot2::ggplot(pc,ggplot2::aes(PC1,PC2,color=group))+
        ggplot2::geom_point(size=2.8)+
        ggplot2::labs(title=paste(id,"PCA QC"),
                      subtitle="Top 500 variable features within this contrast")+
        ggplot2::theme_bw(base_size=11)
      ggplot2::ggsave(file.path(OUT_DIR,"figures",paste0(safe(id),"_PCA_QC.png")),
                      plot_pc,width=7,height=5.5,dpi=250,bg="white")
    }
  }
  out_de <- d;out_de$dataset<-study;out_de$tissue<-tissue;out_de$cohort<-cohort
  out_de$contrast<-cmp
  wcsv(out_de,file.path(OUT_DIR,"DE",paste0(safe(id),"_gene_effects.csv")))
  wcsv(enr,file.path(OUT_DIR,"GSEA",paste0(safe(id),"_allMSigDB.csv")))
  all_de[[id]] <<- out_de
  all_pathways[[id]] <<- selected
  all_modules[[id]] <<- m
  all_gene[[id]] <<- g
  all_summary[[id]] <<- summ
  all_samples[[id]] <<- sc
  all_scoretests[[id]] <<- test
  all_model_symbols[[id]] <<- dh[,c("human_symbol","log2FoldChange"),drop=FALSE]
  msg("FINISHED ",id,"  genes=",nrow(d)," DEGs=",summ$n_DEGs_FDR05,
      "  module coverage=",n_eval,"  match=",round(summ$pct_module_concordant,1),"%")
  invisible(NULL)
}

analyse_study <- function(study,cfg) {
  is_jaxa <- cfg$platform=="jaxa"
  if (is_jaxa) {
    dataset <- sub("^MHU","MHU",study)
    datafile <- find_named(study,paste0("20260331_COUNTtable_for_",dataset,".tsv"))
    if (is.na(datafile)) return(list(status="MISSING_DATA",detail="JAXA ibSLS COUNTtable TSV not found"))
    metafile <- find_named(study,paste0("20260331_Metadata_for_",dataset,".tsv"),"metadata")
    if (is.na(metafile)) return(list(status="MISSING_METADATA",detail="JAXA ibSLS Metadata TSV not found"))
    obj <- read_counts(datafile)
    meta <- read_jaxa_metadata(metafile,colnames(obj$mat),study)
    # JAXA tables may contain multiple tissues, selected by tissue metadata.
    if (!any(meta$tissue=="SOLEUS")) stopf(study," lacks identified soleus samples")
    age_groups <- "ALL"
  } else {
    if (cfg$platform=="rnaseq") {
      f <- paste0("GLDS-",cfg$glds,
        "_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv")
      datafile <- find_named(study,f)
    } else {
      f <- paste0("GLDS-",cfg$glds,
        "_array_normalized_expression_probeset_GLmicroarray.csv")
      datafile <- find_named(study,f)
      if (is.na(datafile) && study=="OSD111") return(list(status="NEEDS_NORMALIZED_ARRAY",
         detail="Uploaded GLDS-111_array_raw_intensities_probe_GLmicroarray.csv is not suitable. Download GLDS-111_array_normalized_expression_probeset_GLmicroarray.csv"))
    }
    if (is.na(datafile)) return(list(status="MISSING_DATA",detail=paste("Data matrix not found:",f)))
    metafile <- find_isa(study,cfg)
    if (is.na(metafile)) return(list(status="MISSING_METADATA",detail=paste("ISA s_*.txt or zip not found for",study)))
    obj <- if (cfg$platform=="rnaseq") read_counts(datafile) else read_normalized_array(datafile)
    isa <- read_isa(metafile)
    meta <- metadata_assign(isa,colnames(obj$mat),study,cfg)
    age_groups <- unique(meta$age_group[!is.na(meta$age_group)])
  }
  wcsv(meta,file.path(OUT_DIR,"QC",paste0(study,"_original_sample_assignments.csv")))
  if (obj$kind=="counts") {
    pp <- pool_technical(obj$mat,meta,study)
    obj$mat<-pp$mat;meta<-pp$meta
  }
  wcsv(meta,file.path(OUT_DIR,"QC",paste0(study,"_final_biological_samples.csv")))
  q <- unique(meta[,c("source","sex","age_at_launch","tissue","age_group","group"),drop=FALSE])
  all_qc[[study]] <<- data.frame(dataset=study,mission=cfg$mission,
     tissue=q$tissue,cohort=q$age_group,sex=q$sex,age_at_launch=q$age_at_launch,
     group=q$group,biological_source=q$source,stringsAsFactors=FALSE)
  if (anyDuplicated(meta$sample) || !identical(colnames(obj$mat),meta$sample))
    stopf(study,": sample ordering mismatch after pooling")
  ann <- annotate_genes(rownames(obj$mat),obj$annotations,cfg$species)
  wcsv(ann,file.path(OUT_DIR,"QC",paste0(study,"_gene_annotation.csv")))
  done <- 0L
  for (tissue in cfg$target) {
    for (age in age_groups) {
      idx <- which(meta$tissue==tissue & meta$age_group==age & !is.na(meta$group))
      if (!length(idx)) {msg(study," ",tissue," ",age,": no mapped samples; skipped");next}
      mt <- meta[idx,,drop=FALSE];mm <- obj$mat[,mt$sample,drop=FALSE]
      groups <- unique(mt$group)
      contrasts <- if (is_jaxa) {
        if (study=="MHU1") list(c("MG","AG","primary"),c("MG","GC","secondary"),
                                  c("AG","GC","secondary"))
        else list(c("PG","GC","primary"))
      } else {
        baseline <- if ("HGC" %in% groups) "HGC" else "GC"
        list(c("FLT",baseline,"primary"),
             if ("VC" %in% groups) c("FLT","VC","secondary") else NULL)
      }
      contrasts <- Filter(Negate(is.null),contrasts)
      for (v in contrasts) {
        if (!all(v[1:2] %in% groups)) {
          if (v[3]=="primary") stopf(study," ",tissue,": required comparison missing: ",
                                      paste(v[1:2],collapse=" vs "),
                                      ". Groups=",paste(groups,collapse=", "))
          next
        }
        id <- paste0(v[1],"_vs_",v[2])
        fit <- fit_contrast(mm,mt,v[1],v[2],obj$kind)
        append_result(study,cfg$mission,tissue,age,id,v[3]=="primary",fit,ann,cfg$species)
        done <- done+1L
      }
    }
  }
  if (!done) stopf(study,": NO comparisons fitted: check metadata tissue, group, cohort")
  list(status="ANALYSED",detail=paste0(done," contrast(s); input=",basename(datafile),
                                         "; metadata=",basename(metafile)))
}

for (name in names(registry)) {
  msg("=== ",name," [",registry[[name]]$mission,"] ===")
  r <- tryCatch(analyse_study(name,registry[[name]]),error=function(e) {
    list(status="FAILED_QC",detail=conditionMessage(e))})
  inclusion[[name]] <- data.frame(dataset=name,mission=registry[[name]]$mission,
    platform=registry[[name]]$platform,species=registry[[name]]$species,
    expected_tissues=paste(registry[[name]]$target,collapse=";"),
    status=r$status,details=r$detail,stringsAsFactors=FALSE)
  msg(name," status: ",r$status,"  ",r$detail)
}
INCL <- do.call(rbind,inclusion)
wcsv(INCL,file.path(OUT_DIR,"QC","STUDY_INCLUSION_MANIFEST.csv"))
if (!length(all_summary)) stopf("No studies analysed. Check QC/STUDY_INCLUSION_MANIFEST.csv ",
                                "and the required file list in the README.")

merge_table <- function(v) {
  if (!length(v)) return(data.frame())
  data.table::rbindlist(v,fill=TRUE,use.names=TRUE) |> as.data.frame()
}
S <- merge_table(all_summary)
P <- merge_table(all_pathways)
M <- merge_table(all_modules)
G <- merge_table(all_gene)
SC <- merge_table(all_samples)
ST <- merge_table(all_scoretests)
# One fixed-module primary-testing family across *all* actual-flight contrasts
# assessed in this run; separately correct secondary/sensitivity contrasts.
M$FDR_primary_across_ALL_spaceflight <- NA_real_
M$FDR_secondary_across_ALL_spaceflight <- NA_real_
ix <- which(M$primary & is.finite(M$pval))
M$FDR_primary_across_ALL_spaceflight[ix] <- p.adjust(M$pval[ix],method="BH")
ix <- which(!M$primary & is.finite(M$pval))
M$FDR_secondary_across_ALL_spaceflight[ix] <- p.adjust(M$pval[ix],method="BH")
# pathways: report FDR from the complete MSigDB universe PER contrast,
# plus panel-wide BH among 9 * all primary spaceflight contrasts.
P$FDR_primary_9panel_all_contrasts <- NA_real_
ix <- which(P$primary & is.finite(P$pval))
P$FDR_primary_9panel_all_contrasts[ix] <- p.adjust(P$pval[ix],method="BH")
ST$FDR_primary_scores <- NA_real_
ix <- which(ST$primary & is.finite(ST$pvalue))
ST$FDR_primary_scores[ix] <- p.adjust(ST$pvalue[ix],method="BH")
if (length(all_qc)) wcsv(merge_table(all_qc),
    file.path(OUT_DIR,"QC","sex_age_tissue_sample_registry.csv"))
wcsv(S,file.path(OUT_DIR,"QC","flight_contrast_summary.csv"))
wcsv(P,file.path(OUT_DIR,"GSEA","six_SOD1_KEGG_and_three_Hallmarks.csv"))
wcsv(M,file.path(OUT_DIR,"module","frozen_module_GSEA.csv"))
wcsv(G,file.path(OUT_DIR,"module","frozen_module_gene_directionality.csv"))
wcsv(SC,file.path(OUT_DIR,"module","sample_scores.csv"))
wcsv(ST,file.path(OUT_DIR,"module","sample_score_group_tests.csv"))
MS <- merge(S,M[,c("dataset","mission","tissue","cohort","contrast","pathway",
                     "NES","pval","padj","concordant_with_ALS",
                     "FDR_primary_across_ALL_spaceflight","FDR_secondary_across_ALL_spaceflight")],
          by=c("dataset","mission","tissue","cohort","contrast"),all.x=TRUE)
wcsv(MS,file.path(OUT_DIR,"module","manuscript_module_summary.csv"))

# Correlations computed only across shared mapped human-symbol gene effects;
# never report individual genes as biological replicates, or infer causation.
all_ids <- names(all_model_symbols)
pair_rows <- list()
if (length(all_ids)>=2) for (j in seq_len(length(all_ids)-1L))
  for (k in seq.int(j+1L,length(all_ids))) {
    a<-all_model_symbols[[j]]; b<-all_model_symbols[[k]]
    common<-merge(a,b,by="human_symbol",suffixes=c("_a","_b"))
    valid<-is.finite(common$log2FoldChange_a) & is.finite(common$log2FoldChange_b)
    common<-common[valid,,drop=FALSE]
    pair_rows[[length(pair_rows)+1L]]<-data.frame(
      contrast1=all_ids[j],contrast2=all_ids[k],n_shared=nrow(common),
      spearman_rho=if (nrow(common)>50) cor(common$log2FoldChange_a,
                    common$log2FoldChange_b,method="spearman") else NA_real_)
  }
if (length(pair_rows)) wcsv(merge_table(pair_rows),
  file.path(OUT_DIR,"concordance","genomewide_cross_study_correlations.csv"))

save_fig <- function(p,name,w=14,h=8) {
  base <- file.path(OUT_DIR,"figures",name)
  ggplot2::ggsave(paste0(base,".png"),p,width=w,height=h,dpi=350,bg="white")
  ggplot2::ggsave(paste0(base,".svg"),p,width=w,height=h,bg="white",device=svglite::svglite)
}
# Sort by mission for final consolidated figure. Asterisk marks GLOBAL primary
# module FDR significance, not per-study padj.
mp <- M[M$primary,,drop=FALSE]
mp$label <- paste(mp$mission,mp$dataset,mp$tissue,mp$cohort,sep=" | ")
mp$label <- factor(mp$label,levels=rev(unique(mp$label)))
mp$component <- factor(mp$pathway,levels=c("ALS_DOWN","ALS_UP"))
mp$mark <- ifelse(!is.na(mp$FDR_primary_across_ALL_spaceflight) &
                  mp$FDR_primary_across_ALL_spaceflight<ALPHA,"*","")
fig_m <- ggplot2::ggplot(mp,ggplot2::aes(NES,label,color=component))+
  ggplot2::geom_vline(xintercept=0,linetype="dashed",linewidth=.3)+
  ggplot2::geom_point(position=ggplot2::position_dodge(width=.45),size=2.8)+
  ggplot2::geom_text(ggplot2::aes(label=mark),position=ggplot2::position_dodge(width=.45),
                    hjust=-0.65,color="black",size=5,show.legend=FALSE)+
  ggplot2::labs(title="ALS-associated load-sensitive module in actual spaceflight",
    subtitle="Fixed 174-gene reference. * BH FDR < 0.05 across all primary module tests.",
    x="Normalized enrichment score (flight or partial gravity vs indicated control)",
    y=NULL,color="Module component")+
  ggplot2::theme_bw(base_size=12)+ggplot2::theme(legend.position="bottom")
save_fig(fig_m,"Figure_UNIFIED_frozen_module_flight",15,max(6,3+0.42*length(unique(mp$label))))
pp <- P[P$primary,,drop=FALSE]
pp$label <- paste(pp$mission,pp$dataset,pp$tissue,pp$cohort,sep=" | ")
pp$label <- factor(pp$label,levels=rev(unique(pp$label)))
pp$mark <- ifelse(pp$padj<ALPHA,"*","")
pp$pathway <- factor(pp$pathway,levels=c(sets_mouse$ref,sets_mouse$hall))
fig_p <- ggplot2::ggplot(pp,ggplot2::aes(pathway,label,color=NES))+
  ggplot2::geom_point(size=2.5)+
  ggplot2::geom_text(ggplot2::aes(label=mark),nudge_y=.27,size=4,color="black")+
  ggplot2::scale_color_gradient2(low="#3B6BA5",mid="white",high="#BC5858",midpoint=0)+
  ggplot2::labs(title="Prespecified ALS-reference metabolic pathways in actual spaceflight",
    subtitle="* FDR < 0.05 in full MSigDB enrichment universe within each contrast",
    x=NULL,y=NULL,color="NES")+
  ggplot2::theme_bw(base_size=11)+ggplot2::theme(axis.text.x=ggplot2::element_text(
        angle=50,hjust=1,size=9))
save_fig(fig_p,"Figure_UNIFIED_metabolic_pathway_panel",15,max(6,3+0.4*length(unique(pp$label))))
# Report session versions, study inclusion, testing family and limits.
capture.output(sessionInfo(),file=file.path(OUT_DIR,"QC","sessionInfo.txt"))
writeLines(c("UNIFIED ACTUAL SPACEFLIGHT: PROVENANCE AND ANALYSIS BOUNDARIES",
  paste0("Date: ",Sys.Date()),
  "Only included original actual-flight or in-orbit partial-gravity samples and matched controls.",
  "Excluded hindlimb unloading, reloading, dry immersion, exercise and denervation.",
  "One voom/limma or normalized-microarray/limma DE model per experimental contrast.",
  "RSEM combined expected counts lack effective-length offsets; analysis is not tximport.",
  "Module was frozen from SOD1/unloading before spaceflight validation.",
  "Pathways: exact six original KEGG reference terms plus three Hallmarks.",
  "Pathway padj = GSEA correction across full KEGG legacy + Hallmark universe per contrast.",
  "FDR_primary_9panel_all_contrasts = supplementary panel-wide BH sensitivity.",
  "Module FDR_primary_across_ALL_spaceflight = BH over all completed primary flight comparisons x both components.",
  "Secondary vivarium or alternative-gravity comparisons adjusted separately.",
  "Within-mission muscle contrasts share animals/flight context and are not independent missions.",
  "Age/sex/mission effects confounded; no across-study causation inference.",
  "Prior manuscript results may differ because platform, normalization, mapping or GSEA background changed.",
  "JAXA ibSLS mission-prefixed GroupID labels explicitly mapped: MHU01_*; MHU04_*; MHU05_*.",
  paste0("INCLUSION status: ",paste(INCL$dataset,INCL$status,sep="=",collapse="; "))),
  file.path(OUT_DIR,"QC","METHODS_AUDIT_TRAIL.txt"))
msg("COMPLETE analyses=",nrow(S)," primary=",sum(S$primary),
  " registered studies=",nrow(INCL)," analysed=",sum(INCL$status=="ANALYSED"))
if (any(INCL$status!="ANALYSED")) {
  warning("PARTIAL RUN. Missing/failed studies listed in QC/STUDY_INCLUSION_MANIFEST.csv. ",
          "Do not label consolidated figures as the complete spaceflight collection.")
  if (STRICT_ALL_STUDIES) stopf("STRICT_ALL_STUDIES: incomplete study inclusion; see QC manifest")
}
