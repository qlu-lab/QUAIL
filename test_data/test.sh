#!/bin/bash

# Step1: Obtain the quantile integrated rank score
cd QUAIL
/s/bin/Rscript Step1_QUAIL_rank_score.R \
--pheno test_data/pheno_test.txt \
--covar test_data/covar_test.txt \
--output test_data/pheno_rank_score.txt \
--num_levels 2000 \
--num_cores 5

cd QUAIL
/s/bin/Rscript Step1_QUAIL_rank_score.R \
--pheno test_data/pheno_test.txt \
--covar test_data/covar_test.txt \
--output test_data/pheno_rank_score.txt \
--num_levels 2000

### Step2: Perform genome-wide vQTL analysis using plink2
cd QUAIL
chmod a+x plink/plink2

/s/bin/Rscript Step2_QUAIL_vQTL.R \
--pheno_rs test_data/pheno_rank_score.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_scan \
--plink_path plink/plink2 \
--dispersion \
--pheno test_data/pheno_test.txt

/s/bin/Rscript Step2_QUAIL_vQTL.R \
--pheno_rs test_data/pheno_rank_score.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_scan \
--plink_path plink/plink2

# Old QUAIL
cd QUAIL
/s/pkg/linux64/R/3.5.1/bin/Rscript QUAIL_vQTL.R \
--pheno_rs test_data/pheno_rank_score.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_sca_old \
--num_cores 5

cd QUAIL
/s/pkg/linux64/R/3.5.1/bin/Rscript QUAIL_Disp.R \
--pheno_rs test_data/pheno_rank_score.txt \
--pheno test_data/pheno_test.txt \
--covar test_data/covar_test.txt \
--geno test_data/test \
--output test_data/vQTL_sca_old \
--analysis both \
--num_cores 5