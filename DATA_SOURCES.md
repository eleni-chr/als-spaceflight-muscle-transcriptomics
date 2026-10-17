# Data sources and input manifest

This manifest documents the external and supplementary data inputs used by the analysis scripts. **Input data are not part of the public code release.** Obtain them from the original source and review that source's current terms before use. This manifest is for provenance and reproducibility; it does not grant rights to redistribute third-party material.

Local filenames and paths below describe the analysis bundle used to prepare the scripts. The source repository may revise or reorganize files. Where an exact download filename is not stable, the accession and local pattern identify the expected input.

## ALS reference workbook

| Source | Local input | Used by | Source and citation |
|---|---|---|---|
| Supplementary tables from the hSOD1G93A skeletal-muscle study | SOD_data.xlsx (Tables S7–S12 are used for the KEGG reference) | Scripts 01, 02, 05 | Moody et al. *Dysregulated microRNA-mRNA networks in SOD1G93A skeletal muscle reveal metabolic impairments in amyotrophic lateral sclerosis*. Disease Models & Mechanisms (2026). [Article DOI](https://doi.org/10.1242/dmm.052963). Retrieve the supplementary workbook from the publisher record and cite the article. |

## NASA Open Science Data Repository

For each OSDR study, use the dataset's current citation and license information from its OSDR record. RSEM and microarray files are not interchangeable; use the file type specified below.

| OSDR study | GLDS record | Local inputs used | Scripts | Source |
|---|---|---|---|---|
| OSD-101 | GLDS-101 | Per-sample RSEM *.genes.results files; unnormalized count matrix and ISA sample metadata for the integrated analysis | 02, 07 | [OSD-101](https://osdr.nasa.gov/bio/repo/data/studies/OSD-101) |
| OSD-104 | GLDS-104 | Unnormalized RSEM count matrix and ISA sample metadata | 07 | [OSD-104](https://osdr.nasa.gov/bio/repo/data/studies/OSD-104) |
| OSD-111 | GLDS-111 | Normalized microarray probeset matrix and ISA sample metadata | 07 | [OSD-111](https://osdr.nasa.gov/bio/repo/data/studies/OSD-111) |
| OSD-21 | GLDS-21 | Normalized microarray probeset matrix and ISA sample metadata | 02, 07 | [OSD-21](https://osdr.nasa.gov/bio/repo/data/studies/OSD-21) |
| OSD-401 | GLDS-401 | Unnormalized RSEM count matrix and ISA sample metadata | 07 | [OSD-401](https://osdr.nasa.gov/bio/repo/data/studies/OSD-401) |
| OSD-419 | GLDS-419 | Per-sample RSEM *.genes.results files, unnormalized RSEM count matrix, optional STAR count matrix, and ISA sample metadata | 02, 07 | [OSD-419](https://osdr.nasa.gov/bio/repo/data/studies/OSD-419) |
| OSD-540 | GLDS-540 | Per-sample RSEM *.genes.results files | 02 | [OSD-540](https://osdr.nasa.gov/bio/repo/data/studies/OSD-540) |
| OSD-665 | GLDS-604 | Unnormalized RSEM count matrix and ISA sample metadata | 07 | [OSD-665](https://osdr.nasa.gov/bio/repo/data/studies/OSD-665) |
| OSD-666 | GLDS-605 | Unnormalized RSEM count matrix and ISA sample metadata | 07 | [OSD-666](https://osdr.nasa.gov/bio/repo/data/studies/OSD-666) |
| OSD-684 | GLDS-615 | Per-sample RSEM *.genes.results files, unnormalized RSEM count matrix, and ISA sample metadata | 02, 07 | [OSD-684](https://osdr.nasa.gov/bio/repo/data/studies/OSD-684) |
| OSD-770 | GLDS-673 | Unnormalized RSEM count matrix and ISA sample metadata | 07 | [OSD-770](https://osdr.nasa.gov/bio/repo/data/studies/OSD-770) |

OSDR licenses and use conditions can vary by dataset and file. Check the study record before reuse; public access alone should not be treated as a redistribution license.

## JAXA / ibSLS

| Source | Local input | Used by | Source and citation |
|---|---|---|---|
| Mouse Habitat Unit missions MHU-1, MHU-4, and MHU-5 | 20260331_COUNTtable_for_MHU1.tsv, 20260331_COUNTtable_for_MHU4.tsv, 20260331_COUNTtable_for_MHU5.tsv, and matching 20260331_Metadata_for_MHU*.tsv files | 07 | Download from the [ibSLS downloads page](https://ibsls.megabank.tohoku.ac.jp/downloads). Follow the current [ibSLS conditions of use](https://ibsls.megabank.tohoku.ac.jp/condition_of_use), including the mission-specific and database citations and acknowledgement. |

## NCBI Gene Expression Omnibus

Cite each GEO Series accession and its linked source publication. NCBI's [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) notes that submitters may hold rights in submitted materials; check the record and linked publication for any additional conditions.

| Accession | Local input | Used by | Source |
|---|---|---|---|
| GSE271607 | GSE271607_Read_counts.txt | 02 | [GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271607) |
| GSE250455 | GSE250455_RawCounts.txt | 05 | [GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250455) |
| GSE215364 | Ten sample-level GSM*_geneCOUNT.txt files in GSE215364/; GEO sample metadata are retrieved by script 12 when needed | 12 | [GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215364) |
| GSE41414 | Series matrix and the study's CEL files in GSE41414/; GEO metadata are retrieved by script 12 when needed | 12 | [GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41414) |
| GSE217576 | GSE217576_bulk_RNA-seq_DENVvsCTL.xlsx | 15 | [GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217576) |
| GSE117440 | Nine sample-level processed tables in GSE117440/ | 15 | [GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117440) |

## Proteomics supplementary materials

| Source | Local input | Used by | Source and citation |
|---|---|---|---|
| PXD028435 / Murgia et al. | Supplementary workbooks MOESM3–MOESM10 and MOESM12 in PXD028435/ | 17 | [ProteomeXchange record](https://www.ebi.ac.uk/pride/archive/projects/PXD028435); [Murgia et al. article](https://www.nature.com/articles/s41526-024-00406-3). |
| PXD081616 / Blottner et al. | Supplementary workbook and SDRF metadata in PXD081616/ | 17 | [ProteomeXchange record](https://www.ebi.ac.uk/pride/archive/projects/PXD081616); [Blottner et al. article](https://www.mdpi.com/1422-0067/24/4/4095). |

The article pages identify the cited papers as open access under CC BY 4.0. Check each downloaded supplement for a separate rights notice, and follow the current source-record terms for repository data. For PRIDE data, EMBL-EBI's [terms of use](https://www.ebi.ac.uk/about/terms-of-use/) state that original data-owner terms continue to apply.

## Reproduction notes

- Download inputs from the source links above and place them at the paths or with the filename patterns listed here. The .gitignore file excludes these local inputs from Git while leaving them available for local analysis.
- Run analysis scripts in numeric order. Figure scripts may depend on tables or workbooks generated by earlier analysis scripts; those derived files are not source inputs and are not included here.
- Script 12 uses GEOquery to obtain GEO metadata when no local cache is available, so that step requires internet access. Some annotation and pathway operations in the scripts may also access online resources.
- The repository's Apache-2.0 license applies to the R code only. It does not apply to data obtained from NASA, JAXA/ibSLS, GEO, PRIDE, or publishers.