# QUAIL: an framework to estimate genetic effects on the variance of quantitative traits

## Introduction

QUAIL (**qua**ntile **i**ntegral **l**inear model) is a quantile regression-based framework to estimate genetic effects on the variance of quantitative traits. QUAIL can be used in 

* **Genome-wide vQTL analysis** - QUAIL constructs a quantile integral phenotype which aggregates information from all quantile levels, and only requires fitting two linear regressions per SNP in genome-wide analysis.
* **Evaluating the vPGS performance** - QUAIL can be extended to continuous predictors such as vPGS and quantify the performance of vPGS in predicting the phenotypic variability.

![QUAIL workflow](https://github.com/qlu-lab/QUAIL/blob/main/Fig/QUAIL_Workflow.png)

## Updates
- Apr 13, 2021: Initial release. Release the codes for Genome-wide vQTL analysis and evaluating the vPGS performance.

## Prerequisites

The software is developed using R and tested in Linux environments. The statistical computing software R (>=3.5.1) and the following R packages are required:

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) (>=1.11.8)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html) (>=1.6.6)
* [quantreg](https://cran.r-project.org/web/packages/quantreg/index.html) (>=5.85)
* [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) (>=3.5.1)
* [BEDMatrix](https://cran.r-project.org/web/packages/BEDMatrix/index.html) (>=2.0.3)

## Genome-wide vQTL analysis
There are two steps to conduct the Genome-wide vQTL analysis:

### Step1: Obtain the quantile integrated rank score
The following analysis transform the phenotype in pheno.txt into the quantile integrated rank score after adjusting the covaraites in covar.txt using 2000 quantile levels and 5 cores are used for parallel computing.
```bash
$ Rscript Obtain_Rank_Score.R \
  --pheno pheno.txt \
  --covar covar.txt \
  --output pheno_rank_score.txt \
  --num_levels 2000 \
  --num_cores 5
```
where the inputs are

| Flag | Description |
|-----|------------------------------------------------------------------------|
| pheno      | The path to the phenotype file |
| covar         | The path to the covariate file |
| output        | The path to the output phenotypic rank score file |                                                    
| num_levels     | Number of quantile levels to fit |
| num_cores        | Number of cores for parellel computing |

### Step2: Perform Genome-wide vQTL analysis
```bash
$ Rscript QUAIL_vQTL.R \
  --pheno_rs pheno_rank_score.txt \
  --geno test \
  --covar covar.txt \
  --output output_1-100.txt \
  --num_cores 5 \
  --start 1 \
  --end 1000
```
where the inputs are

| Flag | Description |
|-----|------------------------------------------------------------------------|
| pheno_rs      | The path to output phenotypic rank score file from step1|
| geno         | The path of the genotype file following the [PLINK format](https://www.cog-genomics.org/plink/1.9). |
| covar        | The path to the covariate file |                                                    
| output     | The path to the output summary statistics file |
| num_cores        | Number of cores for parellel computing |
| start         | (Optinoal) Index of SNP that starts computing |
| end       | (Optinoal) Index of SNP that ends computing |


#### Explanation of Output

The final result has the following fields:

| Column | Description |
|-----|-------------|
| CHR | The chromosomal location of the gene |
| Nsnps | Number of SNPs in the gene expression imputation model |
| Nsnps.used | Number of SNPs used in building the association test |                                                 
| Gene | The name of the gene |
| Matching | The matching method for analysing the trios |
| Beta | The estimated effect size |
| SE | The estimated standard error of Beta |
| Z | The Z test statistic for testing transmission disequilibrium |
| P | The P-value for testing transmission disequilibrium |

## Evaluate the predictive performance of vPGS

```bash
$ Rscript QUAIL_vPGS.R \
  --pheno pheno.txt \
  --vpgs  test.all.score \
  --covar covar.txt \
  --output output_vpgs.txt \
  --num_levels 500 \
  --num_cores 5
```

The final result has the following fields:

| Column | Description |
|-----|-------------|
| CHR | The chromosomal location of the gene |
| Nsnps | Number of SNPs in the gene expression imputation model |
| Nsnps.used | Number of SNPs used in building the association test |                                                 
| Gene | The name of the gene |
| Matching | The matching method for analysing the trios |
| Beta | The estimated effect size |
| SE | The estimated standard error of Beta |
| Z | The Z test statistic for testing transmission disequilibrium |
| P | The P-value for testing transmission disequilibrium |

## Credits
