#!/bin/bash

# Step1: Obtain the quantile integrated rank score
cd QUAIL
Rscript Step1_QUAIL_rank_score.R \
--pheno test_data/pheno_test.txt \
--covar test_data/covar_test.txt \
--output test_data/pheno_rank_score.txt \
--num_levels 2000 \
--num_cores 5

cd QUAIL
Rscript Step1_QUAIL_rank_score.R \
--pheno test_data/pheno_test.txt \
--covar test_data/covar_test.txt \
--output test_data/pheno_rank_score.txt \
--num_levels 2000

### Step2: Perform genome-wide vQTL analysis using plink2
cd QUAIL
chmod a+x plink/plink2

Rscript Step2_QUAIL_vQTL.R \
--pheno_rs test_data/pheno_rank_score.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_scan \
--plink_path plink/plink2 \
--dispersion \
--pheno test_data/pheno_test.txt

Rscript Step2_QUAIL_vQTL.R \
--pheno_rs test_data/pheno_rank_score.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_scan \
--freq test_data/test.afreq \
--plink_path plink/plink2

Rscript Step2_QUAIL_vQTL.R \
--pheno_rs test_data/pheno_rank_score.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_scan \
--plink_path plink/plink2

/z/Comp/lu_group/Members/jmiao24/Software/QUAIL/plink/plink2 \
--bfile /z/Comp/lu_group/Members/jmiao24/Software/QUAIL/test_data/test \
--freq 

