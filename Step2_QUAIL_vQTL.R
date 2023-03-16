# Packages
suppressMessages(library(data.table))
suppressMessages(library(quantreg))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(quadprog))
suppressMessages(library(progress))
suppressMessages(library(MASS))

options(stringsAsFactors = F)
option_list = list(
  make_option("--pheno_rs", action = "store", default = NA, type = "character", help = "The path of the phenotype rank score file"),
  make_option("--geno", action = "store", default = NA, type = "character", help = "The path of the genotype file"),
  make_option("--covar", action = "store", default = NA, type = "character", help = "The path of the covariate file"),
  make_option("--output", action = "store", default = NA, type = "character", help = "The prefix of of the output file"),
  make_option("--plink_path", action = "store", default = NA, type = "character", help = "The path to the plink files"),
  make_option("--dispersion", action = "store_true", default = FALSE, help = "whether to perform the dispersion test"),
  make_option("--pheno", action = "store", default = NA, type = "character", help = "The path to the phenotype file"),
  make_option("--freq", action = "store", default = NA, type = "character", help = "The path of an external allele frequency file")
)

opt = parse_args(OptionParser(option_list=option_list))

cat("********************************************************************* \n")
cat("* Quantile integral linear model (QUAIL) \n")
cat("* Version 1.0.0 \n")
cat("* Genome-wide vQTL analysis \n")
cat("* Step2: Performing genome-wide vQTL analysis \n")
cat("* (C) Jiacheng Miao \n")
cat("* University of Wisconsin–Madison \n")
cat("* https://github.com/qlu-lab/QUAIL \n")
cat("* GNU General Public License v3 \n")
cat("********************************************************************* \n \n")

# --- 0. I/O check
if (is.na(opt$pheno_rs) | is.na(opt$geno) |  is.na(opt$covar) | is.na(opt$output)) {
    cat("ERROR: Missing essential inputs.\n")
    q("no")
}

cat("Options in effect: \n")
cat("Rscript Step2_QUAIL_vQTL.R \\ \n")
cat(paste0("--pheno ", opt$pheno_rs, " \\ \n"))
cat(paste0("--pheno_rs ", opt$pheno_rs, " \\ \n"))
cat(paste0("--geno ", opt$geno, " \\ \n"))
cat(paste0("--covar ", opt$covar, " \\ \n"))
cat(paste0("--output ", opt$output, " \\ \n"))
cat(paste0("--plink_path ", opt$plink_path, " \\ \n"))
if (!is.na(opt$dispersion)){
  cat(paste0("--dispersion ", opt$dispersion, " \\ \n"))
}
if (!is.na(opt$pheno)){
  cat(paste0("--pheno ", opt$pheno, " \\ \n"))
}
if (!is.na(opt$freq)){
  cat(paste0("--freq ", opt$freq, " \\ \n"))
}


# Input parameters
pheno_rank_score <- opt$pheno_rs
genotype <- opt$geno
covariate <- opt$covar
output <- opt$output
plink_path <- opt$plink_path
dispersion <- opt$dispersion
pheno <- opt$pheno
freq <- opt$freq

rank_score <- fread(pheno_rank_score)
rank_score_name <- colnames(rank_score)[3]

# Add the raw phenotype into the covaraites table
if (dispersion){
  pheno_table <- fread(pheno)
  pheno_name <- colnames(pheno_table)[3]
  covar_table <- fread(covariate)
  covar_table$pheno <- pheno_table[, 3]
  fwrite(covar_table, paste0(output, "_covar_tmp.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
}
# Plinlk2 to run genome-wide vQTL analysis
cat("\n### Use plink2 to run the Genome-wide vQTL analysis ###\n")
if (!is.na(freq)){
  job_vqtl <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ", " --read-freq ", freq,  " --pheno ", pheno_rank_score , " --covar ", covariate, " --out ", output, "_QUAIL_vQTL --linear --no-psam-pheno")
  system(job_vqtl)
}else{
  job_vqtl <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ", " --pheno ", pheno_rank_score , " --covar ", covariate, " --out ", output, "_QUAIL_vQTL --linear --no-psam-pheno")
  system(job_vqtl)
}

# Format the  Genome-wide vQTL summary statistics
cat("\n### Format the Genome-wide vQTL summary statistics ###\n")
df <- fread(paste0(output, "_QUAIL_vQTL.", rank_score_name, ".glm.linear"))
df <- df[df$TEST == "ADD", ]
df$A2 <- ifelse(df$A1 == df$REF, df$ALT, df$REF)
colnames(df) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "BETA", "SE", "Z", "P", "A2")
col_out <- c("CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "Z", "P", "N")
df <- as.data.frame(df)[!(is.na(df$P)), col_out]
cat("--- Saving the QUAIL vQTL  results with ", nrow(df), " SNPs\n")
fwrite(df, paste0(output, "_QUAIL_vQTL.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
cat("--- The QUAIL vQTL results are saved in ", paste0(output, "_QUAIL_vQTL.txt") ,"!\n")

if (dispersion){
  cat("\n### Use plink2 to run the Genome-wide dispersion analysis ###\n")
  # Additive GWAS

  if (!is.na(freq)){
    job_add <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ", " --read-freq ", freq, " --pheno ", pheno , " --covar ", covariate, " --out ", output, "_GWAS_add --linear --no-psam-pheno")
    system(job_add)
  }else{
    job_add <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ", " --pheno ", pheno , " --covar ", covariate, " --out ", output, "_GWAS_add --linear --no-psam-pheno")
    system(job_add)
  }
  # Read the additive GWAS
  df_add <- fread(paste0(output, "_GWAS_add.", pheno_name, ".glm.linear"))
  df_add <- df_add[df_add$TEST == "ADD", ]
  df_add$A2 <- ifelse(df_add$A1 == df_add$REF, df_add$ALT, df_add$REF)
  colnames(df_add) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "BETA", "SE", "Z", "P", "A2")
  col_out <- c("CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "Z", "P")
  df_add <- as.data.frame(df_add)[, col_out]

  # vQTL GWAS with raw phenotype with covaraites

  if (!is.na(freq)){
    job_disp <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ", " --read-freq ", freq,  " --pheno ", pheno_rank_score , " --covar ", paste0(output, "_covar_tmp.txt"), " --out ", output, "_QUAIL_disp --linear --no-psam-pheno")
    system(job_disp)
  }else{
    job_disp <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ",  " --pheno ", pheno_rank_score , " --covar ", paste0(output, "_covar_tmp.txt"), " --out ", output, "_QUAIL_disp --linear --no-psam-pheno")
    system(job_disp)
  }

  # Read the vQTL GWAS with raw phenotype with covaraites
  df_vqtl <- fread(paste0(output, "_QUAIL_disp.", rank_score_name, ".glm.linear"))
  df_vqtl <- df_vqtl[df_vqtl$TEST == "ADD", ]
  df_vqtl$A2 <- ifelse(df_vqtl$A1 == df_vqtl$REF, df_vqtl$ALT, df_vqtl$REF)
  colnames(df_vqtl) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "BETA", "SE", "Z", "P", "A2")
  col_out <- c("CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "Z", "P", "N")
  df_vqtl <- as.data.frame(df_vqtl)[, col_out]


  # Dispersion test
  df_add_MHC <- df_add[df_add$CHR !=6, ]
  df_vqtl_MHC <- df_vqtl[df_vqtl$CHR !=6, ]
  snp_ovp_MHC <- intersect(df_add_MHC$SNP, df_vqtl_MHC$SNP)
  df_add_MHC_ovp <- df_add_MHC[match(snp_ovp_MHC, df_add_MHC$SNP), ]
  df_vqtl_MHC_ovp <- df_vqtl_MHC[match(snp_ovp_MHC, df_vqtl_MHC$SNP), ]

  snp_ovp <- intersect(df_add$SNP, df_vqtl$SNP)
  df_add_ovp <- df_add[match(snp_ovp, df_add$SNP), ]
  df_vqtl_ovp <- df_vqtl[match(snp_ovp, df_vqtl$SNP), ]

  # Only use no MHC
  r_disp <- rlm(df_vqtl_MHC_ovp$BETA ~ 0 + df_add_MHC_ovp$BETA)
  add_noise <- mean(df_add_MHC_ovp$SE^2 , na.rm=T)
  add_noise_adj <- 1+add_noise/(var(df_add_MHC_ovp$BETA,na.rm=T)-add_noise)
  r_disp_adj <- r_disp$coefficients[1]*add_noise_adj
  # Use MHC
  dispersion <- df_vqtl_ovp$BETA-r_disp_adj*df_add_ovp$BETA
  dispersion_se <- sqrt(df_vqtl_ovp$SE^2+(r_disp_adj^2)*df_add_ovp$SE^2)
  dispersion_t <- dispersion/dispersion_se
  dispersion_pval <- pchisq(dispersion_t^2,1,lower.tail=F)
  df_disp_tmp <- data.frame(BETA_disp = dispersion, SE_disp = dispersion_se, Z_disp = dispersion_t, P_disp = dispersion_pval, N = df_vqtl_MHC$N)
  cat("\n")
  df_disp <- df_add_MHC
  df_raw_vqtl <- df[match(df_disp$SNP, df$SNP), ]
  df_disp <- cbind(df_disp, as.data.frame(df_raw_vqtl)[, c("BETA", "SE", "Z", "P")])
  df_disp <- cbind(df_disp, df_disp_tmp)
  df_disp <- df_disp[!(is.na(df_disp$P_disp)), ]
  colnames(df_disp) <- c("CHR", "BP", "SNP", "A1", "A2", "BETA_add", "SE_add", "Z_add", "P_add", "BETA_var", "SE_var", "Z_var", "P_var", "BETA_disp", "SE_disp", "Z_disp", "P_disp", "N")
  cat("--- Saving the QUAIL dispersion results with ", nrow(df_disp), " SNPs\n")
  fwrite(df_disp, paste0(output, "_QUAIL_Disp.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  cat("--- The QUAIL dispersion results are saved in ", paste0(output, "_QUAIL_Disp.txt") ,"!\n")

  # Remove all the redundant files
  system(paste0("rm -rf ", paste0(output, "_covar_tmp.txt")))
  system(paste0("rm -rf ", paste0(output, "_GWAS_add.", pheno_name, ".glm.linear")))
  system(paste0("rm -rf ", paste0(output, "_QUAIL_disp.", rank_score_name, ".glm.linear")))

}
system(paste0("rm -rf ", paste0(output, "_QUAIL_vQTL.", rank_score_name, ".glm.linear")))
cat("\n### Step2 of QUAIL is finished! ###\n")