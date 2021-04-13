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

### Explanation of Output

## Evaluate the predictive performance of vPGS

## Credits
