# QUAIL: a unified framework to estimate genetic effects on the variance of quantitative traits

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

## Download QUAIL

You can download QUAIL by:

```
$ git clone https://github.com/qlu-lab/QUAIL
$ cd ./QUAIL
```


## Input Data Format
### Phenotype file

The input phenotype file need to be an n x 3 table, where n is the sample size. The columns in order are the FID, IID, and the phenotype. Here is an example of the phenotype file `pheno.txt`:
```
FID	IID	bmi
1	1	31.0920070855422
2	2	34.9722947258896
3	3	26.651358296193
4	4	11.3265202767703
5	5	28.1784063613473
...
...
```


### Covariate file

The input Covariate file needs to be an n x (m+2) table, where m is the number of covariates to include. The columns in order are the FID, IID, and the covariates. Here is an example of the covariate file `covar.txt`:

```
FID	IID	sex	age	pc1
1	1	0	40	0.0496566810517587
2	2	1	41	-0.30689382604118
3	3	1	42	0.0662231532694345
4	4	0	43	0.653681021333332
5	5	1	44	0.50856163868585
...
...
```

### Genotype file

The input genotype file needs to be in the plink bed/bim/fam format. The path only includes the prefix, not the suffix. For example, the path to the input genotype file is `geno` where the genotype files are `geno.bed, geno.bim, geno.fam`.

### vPGS file

The input vPGS file needs to be an n x 3 table, where n is the sample size. The columns in order are the FID, IID, and the vPGS. Here is an example of the vPGS file `vpgs.txt`:
```
FID	IID	vpgs
1000	1000	-0.166044638683817
1001	1001	0.325459663014065
1002	1002	-0.0589986740920006
1003	1003	0.176636938025438
1004	1004	0.266151552296392
1005	1005	0.376450240002708
...
...
```

## Genome-wide vQTL analysis
There are two steps to conduct the Genome-wide vQTL analysis:

### Step1: Obtain the quantile integrated rank score
#### Example:
The following script transform the phenotype in `pheno.txt` into the quantile integrated rank score using `2000` quantile levels. It adjusted the covaraites in `covar.txt` and `5` cores are used for parallel computing. The output phenotypic rank score file will be written to `pheno_rank_score.txt`.
```bash
Rscript Obtain_Rank_Score.R \
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
| num_cores        | Number of cores for parallel computing |

### Step2: Perform genome-wide vQTL analysis

#### Example:
The following script perform genome-wide vQTL analysis from the `1`-`1000` SNPs in plink genotype file `geno` using output of quantile integrated rank score in step1 `pheno_rank_score.txt`. This analysis adjusted the covaraites in `covar.txt` and `5` cores are used for parallel computing. 
```bash
Rscript QUAIL_vQTL.R \
  --pheno_rs pheno_rank_score.txt \
  --geno geno \
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
| num_cores        | Number of cores for parallel computing |
| start          | (Optional) Index of SNP that starts computing |
| end       | (Optional) Index of SNP that ends computing |

Note that 
* For both genotype vQTL analysis and vPGS analysis(described following), the `--num_levels` depends on the sample size in the analytic sample. We recommend using 500 levels when the sample size is around 10000 and 2000 levels when the sample size > 10000. However, you can check the robustness of your results by using more levels.

* The script is designed to run on chromosome segments to facilitate parallel computation on the cluster. If `--start` or `--end` is not specified, the script will perform the analysis on all SNPs in the plink `test` file.



#### Explanation of Output

The final result has the following fields:

| Column | Description |
|-----|-------------|
| CHR | Chromosome |
| SNP | rs ID |
| BP | Base pair position |                                                 
| A1 | Allele 1 (effect allele) |
| A2 | Allele 2 (non-effect allele) |
| FREQ | allele frequency of allele 1 |
| BETA | The estimated effect size |
| SE | The estimated standard error of BETA |
| P | The P-value for testing variance effects |
| N | Sample size |

Note that the BETA in the output of QUAIL is the standardized effect size. You can calculate the allele effect size by BETA_allele = BETA/sqrt(2 * FREQ * (1-FREQ)), where BETA and FREQ come from the output of QUAIL.

## Evaluate the predictive performance of vPGS
#### Example:
The following script evaluate the performance of vPGS from the `test.all.score` in predicting the variability of phenotype from `pheno.txt` using `500` quantile levels. The analysis adjusted the covaraites in `covar.txt` and `5` cores are used for parallel computing. The performance is written out to `output_vpgs.txt`.
```bash
Rscript QUAIL_vPGS.R \
  --pheno pheno.txt \
  --vpgs  vpgs.txt \
  --covar covar.txt \
  --output output_vpgs.txt \
  --num_levels 500 \
  --num_cores 5
```
Note that 
* The choice of `--num_levels` can be found in the Genome-wide vQTL analysis part.

The final result has the following fields:

| Column | Description |
|-----|-------------|
| BETA | The estimated effect size of vPGS to phenotypic variability|
| SE | The estimated standard error of BETA |
| P | The P-value for testing variance effects for vPGS |
| N | Sample size |

## Example
```bash
cd ./QUAIL
Rscript Obtain_Rank_Score.R \
  --pheno ./test_data/pheno_test.txt \
  --covar ./test_data/covar_test.txt \
  --output pheno_rank_score.txt \
  --num_levels 2000 \
  --num_cores 5
```


## Citation

If you use QUAIL, please cite

Miao, J., Lin, Y., Wu, Y., Zheng, B., Schmitz, L. L., Fletcher, J. M., & Lu, Q. (2021). [A quantile integral linear model to quantify genetic effects on phenotypic variability.](https://www.biorxiv.org/content/10.1101/2021.04.14.439847v1) bioRxiv, 2021.2004.2014.439847

## License

All rights reserved for [Lu-Laboratory](http://qlu-lab.org/)
