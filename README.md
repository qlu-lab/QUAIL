# QUAIL: a unified framework to estimate genetic effects on the variance of quantitative traits

## QUAIL

QUAIL (**qua**ntile **i**ntegral **l**inear model) is a quantile regression-based framework to estimate genetic effects on the variance of quantitative traits. QUAIL can be used in 

* **Genome-wide vQTL analysis** - QUAIL constructs a quantile integral phenotype which aggregates information from all quantile levels, and only requires fitting two linear regressions per SNP in genome-wide analysis.
* **Evaluating the vPGS performance** - QUAIL can be extended to continuous predictors such as vPGS and quantify the performance of vPGS in predicting the phenotypic variability.

![QUAIL workflow](https://github.com/qlu-lab/QUAIL/blob/main/Fig/QUAIL_Workflow.png)

## Manual

`QUAIL` can be downloaded via `git clone https://github.com/qlu-lab/QUAIL`

Please see the [wiki](https://github.com/qlu-lab/QUAIL/wiki) for the short tutorials describing the two basic functions (Genome-wide vQTL analysis and Evaluating the vPGS performance), as well as the detailed manual of `QUAIL`.

## Version History
- Mar 14, 2023: Speed up the step2 and move the tutorials into wiki.
- Aug 22, 2022: Add the simulation codes.
- Feb 25, 2022: Add the dispersion effects.
- Jan 13, 2022: Add the test data part.
- Apr 13, 2021: Initial release. Release the codes for Genome-wide vQTL analysis and evaluating the vPGS performance.


## Citation

If you use QUAIL, please cite

Miao, J., Lin, Y., Wu, Y., Zheng, B., Schmitz, L. L., Fletcher, J. M., & Lu, Q. (2022). A quantile integral linear model to quantify genetic effects on phenotypic variability. Proceedings of the National Academy of Sciences, 119(39), e2212959119. https://doi.org/doi:10.1073/pnas.2212959119

## Contact

For questions and comments, please open a GitHub issue or contact Jiacheng Miao at jmiao24@wisc.edu.

## "Birds" familial links
* [PIGEON](https://github.com/qlu-lab/PIGEON)  (**P**olygen**I**c **G**ene-**E**nvironment interacti**ON**) is unified statistical framework to estimate polygenic gene-environment (GxE) interactions using GWIS (and GWAS) summary statistics.
