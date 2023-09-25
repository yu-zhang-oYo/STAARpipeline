[![R build status](https://github.com/xihaoli/STAARpipeline/workflows/R-CMD-check/badge.svg)](https://github.com/xihaoli/STAARpipeline/actions)
[![Build Status](https://travis-ci.com/xihaoli/STAARpipeline.svg?branch=main)](https://app.travis-ci.com/github/xihaoli/STAARpipeline)
[![Build status](https://ci.appveyor.com/api/projects/status/ltr225p13idh2934/branch/main?svg=true)](https://ci.appveyor.com/project/xihaoli/staarpipeline/branch/main)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# For MR

**Step 0**: install the STAARpipeline package with modified code in R so that the output also includes IDs

**Note**: 

1. I only change the code for gene_centric_coding and Individual_analysis
2. the format of the output is different from the original one.

```
devtools::install_github("yu-zhang-oYo/STAARpipeline", ref = "STAARpipeline_MR")
```


**Step 1**: Fit STAAR null model
```
cd GA
Rscript STAARpipeline_Null_Model.r
```

**Step 2**: Individual (single variant) analysis
```
mkdir single_output
Rscript STAARpipeline_Individual_Analysis.r 22
```

**Step 3.1**: Gene-centric coding analysis
```
mkdir gene_centric_output
Rscript STAARpipeline_Gene_Centric_Coding.r 22
```

**Step 4**: simple LASSO for all the variables in R

```
library(SeqArray)
library(dplyr)
library(glmnet)
library(AER)

setwd("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/GA")

## aGDS directory
agds_dir <- get(load("../agds_dir.Rdata"))
agds_dir <- paste0(".", agds_dir)
## Null model
obj_nullmodel <- get(load("./outputobj_nullmodel.Rdata"))

agds.path <- agds_dir[1]
genofile <- seqOpen(agds.path)

phenotype.id <- as.character(obj_nullmodel$id_include)

# load individual_analysis results
single <- read.table(file="single_output/individual_analysis_22.out", header=TRUE, stringsAsFactors=F)
single_sig <- single[single$pvalue < 1e-8,]
variantIDs_select_single <- single_sig$variantIDs


# load gene_centric_coding results
load("gene_centric_output/gene_centric_coding_22.Rdata")


# choose significant results for gene_centric_coding results
# Note: choose suitable alpha, so that we can have suitable numbers of variants for LASSO
#.      if too much or too less variants, the prediction is only one values
staar_o_values_sig <- lapply(results_coding, function(x, alpha=1e-4) {
  if ((length(x)>0) && x$'STAAR-O' < alpha) {
    return(x)
  }
})
# Remove NULL values
staar_o_values_sig <- Filter(Negate(is.null), staar_o_values_sig)
variantIDs_select_gene <- unique(unlist(sapply(staar_o_values_sig, function(x) x$variantIDs)))

variantIDs_select_total <- c(variantIDs_select_single, variantIDs_select_gene)
# there may be duplicated IDs in single and gene
variantIDs_select_total <- unique(variantIDs_select_total)

# Function to extract SNP values for a given set of IDs
extract_SNP_values <- function(IDs, genofile, phenotype.id, variant=TRUE) {
  # Extract all available IDs from the genofile
  # different varinat.id may have the same annotation/id
  if(variant){
    seqSetFilter(genofile, variant.id = IDs, sample.id = phenotype.id)
  } else{
    seqSetFilter(genofile, sample.id = phenotype.id)
    seqSetFilterAnnotID(genofile, id = IDs)
  }
  SNP_values <- seqGetData(genofile, "$dosage_alt")  
  SNP_values <- as.matrix(SNP_values)
  # Retrieve the sample IDs
  sample_ids <- seqGetData(genofile, "sample.id")
  variant_ids <- seqGetData(genofile, "variant.id")
  df <- cbind(sample_ids, SNP_values)
  colnames(df) <- c("sample.id", variant_ids)
  df <- as.data.frame(df)
  
  # Reset the filter
  seqResetFilter(genofile)

  return(df)
}

df_total <- extract_SNP_values(IDs=variantIDs_select_total, genofile=genofile, phenotype.id=phenotype.id, variant=TRUE)
seqClose(genofile)

phenotype <- read.csv("./pheno/pheno.csv")
df <- merge(df_total, phenotype, by="sample.id")
var_expo <- "GA"
df <- df[,c(variantIDs_select_total, var_expo)]
df <- na.omit(df)
X_select <- df %>%
  mutate(across(all_of(as.character(variantIDs_select_total)), as.numeric)) %>%
  select(all_of(as.character(variantIDs_select_total))) %>%
  as.matrix()
y_select <- df[ , var_expo]

# simple LASSO for all the variables
cv_model <- cv.glmnet(X_select, y_select, alpha = 1)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(X_select, y_select, alpha = 1, lambda = best_lambda)
coefficients <- coef(best_model)
df$predicted_expo <- predict(best_model, s = best_lambda, newx = X_select)
hist(df$predicted_expo)

# use the prediction for the second stage directly
var_outc <- "GA_fake"
df[[var_outc]] <- df[[var_expo]]
formula_str <- paste(var_outc, " ~ ", "predicted_expo")
second_stage_model <- lm(formula_str, data = df)
# Display the summary of the second stage model
summary(second_stage_model)


# use the prediction of the exposure as IV
var_outc <- "GA_fake"
df[[var_outc]] <- df[[var_expo]]
formula_str <- paste(var_outc, " ~ ", var_expo, "|", "predicted_expo")
model_ivreg <- ivreg(as.formula(formula_str), data = df)
summary(model_ivreg, diagnostics = TRUE)

```




# STAARpipeline
This is an R package for performing association analysis of whole-genome/whole-exome sequencing (WGS/WES) studies using STAAR pipeline.
## Description
**STAARpipeline** is an R package for phenotype-genotype association analyses of biobank-scale WGS/WES data, including single variant analysis and variant set analysis. The single variant analysis in STAARpipeline provides individual *P* values of variants given a MAF or MAC cut-off. The variant set analysis in STAARpipeline includes gene-centric analysis and non-gene-centric analysis of rare variants. The gene-centric coding analysis provides five genetic categories: putative loss of function (pLoF), missense, disruptive missense, pLoF and disruptive missense, and synonymous. The gene-centric noncoding analysis provides eight genetic categories: promoter or enhancer overlaid with CAGE or DHS sites, UTR, upstream, downstream, and noncoding RNA genes. The non-gene-centric analysis includes sliding window analysis with fixed sizes and dynamic window analysis with data-adaptive sizes. STAARpipeline also provides analytical follow-up of dissecting association signals independent of known variants via conditional analysis using <a href="https://github.com/xihaoli/STAARpipelineSummary">**STAARpipelineSummary**</a>.

**STAARpipeline** and **STAARpipelineSummary** are implemented as a collection of apps. Please see the apps <a href="https://github.com/xihaoli/staarpipeline-rap">**staarpipeline**</a>, <a href="https://github.com/xihaoli/staarpipelinesummary_varset-rap">**staarpipelinesummary_varset**</a> and <a href="https://github.com/xihaoli/staarpipelinesummary_indvar-rap">**staarpipelinesummary_indvar**</a> that run on the UK Biobank Research Analysis Platform for more details.
## Workflow Overview
![STAARpipeline_workflow](docs/STAARpipeline_workflow.jpg)
## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 3.5.1)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.
## Dependencies
STAARpipeline links to R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a> and <a href="https://cran.r-project.org/web/packages/RcppArmadillo/index.html">RcppArmadillo</a>, and also imports R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://github.com/xihaoli/STAAR">STAAR</a>, <a href="https://github.com/zilinli1988/SCANG">SCANG</a>, <a href="https://cran.r-project.org/web/packages/dplyr/index.html">dplyr</a>, <a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">SeqArray</a>, <a href="https://bioconductor.org/packages/release/bioc/html/SeqVarTools.html">SeqVarTools</a>, <a href="https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html">GenomicFeatures</a>, <a href="https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html">TxDb.Hsapiens.UCSC.hg38.knownGene</a>, <a href="https://cran.r-project.org/web/packages/GMMAT/index.html">GMMAT</a>, <a href="https://bioconductor.org/packages/release/bioc/html/GENESIS.html">GENESIS</a>, <a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>. These dependencies should be installed before installing STAARpipeline.
## Installation
```
library(devtools)
devtools::install_github("xihaoli/STAARpipeline",ref="main")
```
## Docker Image
A [docker image for STAARpipeline](https://hub.docker.com/repository/docker/zilinli/staarpipeline), including R (version 3.6.1) built with Intel MKL and all STAAR-related packages (STAAR, SCANG, STAARpipeline, STAARpipelineSummary) pre-installed, is located in the Docker Hub. The docker image can be pulled using
```
docker pull zilinli/staarpipeline:0.9.6
```
## Usage
Please see the <a href="docs/STAARpipeline_manual.pdf">**STAARpipeline** user manual</a> for detailed usage of STAARpipeline package. Please see the <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">**STAARpipeline** tutorial</a> for a detailed example of analyzing sequencing data using STAARpipeline.
## Data Availability
The whole-genome functional annotation data assembled from a variety of sources and the precomputed annotation principal components are available at the [Functional Annotation of Variant - Online Resource (FAVOR)](https://favor.genohub.org) site and [FAVOR Essential Database](https://doi.org/10.7910/DVN/1VGTJI).
## Version
The current version is 0.9.6 (January 1, 2023).
## Citation
If you use **STAARpipeline** and **STAARpipelineSummary** for your work, please cite:

Zilin Li*, Xihao Li*, Hufeng Zhou, Sheila M. Gaynor, Margaret Sunitha Selvaraj, Theodore Arapoglou, Corbin Quick, Yaowu Liu, Han Chen, Ryan Sun, Rounak Dey, Donna K. Arnett, Paul L. Auer, Lawrence F. Bielak, Joshua C. Bis, Thomas W. Blackwell, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jennifer A. Brody, Brian E. Cade, Matthew P. Conomos, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Paul S. de Vries, Ravindranath Duggirala, Nora Franceschini, Barry I. Freedman, Harald H. H. Göring, Xiuqing Guo, Rita R. Kalyani, Charles Kooperberg, Brian G. Kral, Leslie A. Lange, Bridget M. Lin, Ani Manichaikul, Alisa K. Manning, Lisa W. Martin, Rasika A. Mathias, James B. Meigs, Braxton D. Mitchell, May E. Montasser, Alanna C. Morrison, Take Naseri, Jeffrey R. O’Connell, Nicholette D. Palmer, Patricia A. Peyser, Bruce M. Psaty, Laura M. Raffield, Susan Redline, Alexander P. Reiner, Muagututi’a Sefuiva Reupena, Kenneth M. Rice, Stephen S. Rich, Jennifer A. Smith, Kent D. Taylor, Margaret A. Taub, Ramachandran S. Vasan, Daniel E. Weeks, James G. Wilson, Lisa R. Yanek, Wei Zhao, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Jerome I. Rotter, Cristen J. Willer, Pradeep Natarajan, Gina M. Peloso, & Xihong Lin. (2022). **A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies**. _Nature Methods_, _19_(12), 1599-1611. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/36303018">36303018</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10008172/">PMC10008172</a>. DOI: <a href="https://doi.org/10.1038/s41592-022-01640-x">10.1038/s41592-022-01640-x</a>.

Xihao Li*, Zilin Li*, Hufeng Zhou, Sheila M. Gaynor, Yaowu Liu, Han Chen, Ryan Sun, Rounak Dey, Donna K. Arnett, Stella Aslibekyan, Christie M. Ballantyne, Lawrence F. Bielak, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jai G. Broome, Matthew P. Conomos, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Barry I. Freedman, Xiuqing Guo, George Hindy, Marguerite R. Irvin, Sharon L. R. Kardia, Sekar Kathiresan, Alyna T. Khan, Charles L. Kooperberg, Cathy C. Laurie, X. Shirley Liu, Michael C. Mahaney, Ani W. Manichaikul, Lisa W. Martin, Rasika A. Mathias, Stephen T. McGarvey, Braxton D. Mitchell, May E. Montasser, Jill E. Moore, Alanna C. Morrison, Jeffrey R. O'Connell, Nicholette D. Palmer, Akhil Pampana, Juan M. Peralta, Patricia A. Peyser, Bruce M. Psaty, Susan Redline, Kenneth M. Rice, Stephen S. Rich, Jennifer A. Smith, Hemant K. Tiwari, Michael Y. Tsai, Ramachandran S. Vasan, Fei Fei Wang, Daniel E. Weeks, Zhiping Weng, James G. Wilson, Lisa R. Yanek, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Benjamin M. Neale, Shamil R. Sunyaev, Gonçalo R. Abecasis, Jerome I. Rotter, Cristen J. Willer, Gina M. Peloso, Pradeep Natarajan, & Xihong Lin. (2020). **Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale**. _Nature Genetics_, _52_(9), 969-983. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/32839606">32839606</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7483769/">PMC7483769</a>. DOI: <a href="https://doi.org/10.1038/s41588-020-0676-4">10.1038/s41588-020-0676-4</a>.
## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
