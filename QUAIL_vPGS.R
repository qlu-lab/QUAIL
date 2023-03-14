# Packages
rm(list = ls())
suppressMessages(require(data.table))
suppressMessages(require(quantreg))
suppressMessages(require(parallel))
suppressMessages(library(optparse))
suppressMessages(library(BEDMatrix))

options(stringsAsFactors = F)
option_list = list(
  make_option("--pheno", action = "store", default = NA, type = "character", help = "The path of the phenotype file"),
  make_option("--vpgs", action = "store", default = NA, type = "character", help = "The path of the vpgs file"),
  make_option("--covar", action = "store", default = NA, type = "character", help = "The path of the covariate file"),
  make_option("--output", action = "store", default = NA, type = "character", help = "The path of the output file"),
  make_option("--num_levels", action = "store", default = NA, type = "integer", help = "Number of quantile levels to fit")
)

opt = parse_args(OptionParser(option_list=option_list))

# --- 0. I/O check
if (is.na(opt$pheno) | is.na(opt$vpgs) |  is.na(opt$covar) | is.na(opt$output) | is.na(opt$num_levels) | is.na(opt$num_cores) ) {
    cat("ERROR: Missing essential inputs.\n")
    q("no")
}

cat("********************************************************************* \n")
cat("* Quantile integral linear model (QUAIL) \n")
cat("* Version 1.0.0 \n")
cat("* Assessing vPGS predictive performance \n")
cat("* (C) Jiacheng Miao \n")
cat("* University of Wisconsin–Madison \n")
cat("* https://github.com/qlu-lab/QUAIL \n")
cat("* GNU General Public License v3 \n")
cat("********************************************************************* \n \n")

# Input parameters
phenotype <- opt$pheno
vpgs_path <- opt$vpgs
covariate <- opt$covar
output <- opt$output
num_levels <- opt$num_levels

# --- 1. Read the input files

# read the phenotype files
cat("Preparing the data:\n")
pheno <- as.data.frame(fread(phenotype,stringsAsFactors = F))
pheno <- subset(pheno,!duplicated(pheno$IID))
pheno[pheno[,3]==-9, 3] <- NA  # Plink format

# read the covariates files
covar <- as.data.frame(fread(covariate,data.table = F,stringsAsFactors = F))
covar <- subset(covar,!duplicated(covar$IID))

# read the vPGS files
vpgs <- as.data.frame(fread(vpgs_path,data.table = F,stringsAsFactors = F))
vpgs <- subset(vpgs, !duplicated(vpgs$IID))

# align phenotype and covariates files
id_intersect <- intersect(intersect(pheno$IID, covar$IID), vpgs$IID)
pheno <- pheno[match(id_intersect, pheno$IID), ]
covar <- covar[match(id_intersect, covar$IID), ]
vpgs <- vpgs[match(id_intersect, vpgs$IID), 3]

# standardized the vpgs
vpgs <- scale(vpgs)

# prepare files for the analysis
covar_vpgs <- covar[, 3:ncol(covar)]
covar_vpgs$vpgs <- vpgs
covar_no_vpgs <- covar[, 3:ncol(covar)]
pheno_qr <- lm(pheno[, 3] ~ vpgs)$residuals


# --- 2. Obtain the quantile integrated rank score

cat("Fitting the quantile regression to obtain quantile rank score:\n")

# Step1: fit the null model for k quantiles
Get_a_i_tau_vpgs <- function(i){
    tau_curr <- i/(num_levels + 1)
    # Fit quantile regression between phenotype and covarites at quantile i/k
    Qreg <- rq(pheno_qr ~ ., data = covar_no_vpgs, tau = tau_curr, method = "fn")
    Qreg_vpgs <- rq(pheno_qr ~ ., data = covar_vpgs, tau = tau_curr, method = "fn")
    # Constrcuct the quantile rank score for each quantile 
    coeff <-  summary(Qreg_vpgs , se = "ker")$coefficients
    SE_vpgs <- coeff[nrow(coeff),2]
    a_i_tau <- (tau_curr - ifelse(residuals(Qreg) < 0 , 1, 0))*SE_vpgs/sqrt(-tau_curr^2 + tau_curr)
    return(a_i_tau)
}

Fit_a_i_tau_vpgs <- function(start = 1, end = num_levels){
    df_one_q <- lapply(start:end, Get_a_i_tau_vpgs)
    return(df_one_q)
}

# To suppress some Warnings that are harmless: https://github.com/HenrikBengtsson/future/issues/218
df_all_q <- suppressWarnings(Fit_a_i_tau_vpgs(start=1, end =num_levels))

# Step2: Add the quantile rank score across quantiles
cat("Obtaining the integrated rank score\n")
int_rank_score <- 0
for (i in 1:num_levels){
    ## upper quantile weight = 1, lower quantile weight = -1
    weight <- ifelse( (i > ((num_levels)/2)), 1, -1)
    int_rank_score <- int_rank_score + weight*df_all_q[[i]]
}

int_rank_score <- int_rank_score/((num_levels)/2)

# --- 3. Use QUAIL to evaluate the vpgs performance

# Obtain the vpgs_star
cat("Begin evaluating vPGS performance:\n")
m_vpgs_covar <- lm(vpgs ~ ., data = covar_no_vpgs)
vpgs_star <- m_vpgs_covar$residual

# Run regression between Y_QI and G_star
Y_QI <- sqrt(length(vpgs_star))*int_rank_score
coeff_vpgs <- summary(lm(Y_QI ~ vpgs_star))$coefficients
df_out <- data.frame(BETA = coeff_vpgs[2, 1], SE = coeff_vpgs[2, 2], P =coeff_vpgs[2, 4], N = length(vpgs_star))
cat("The evaluation result is .\n")
print(df_out)

# Write out the results
cat("Write out the evaluation results.\n")
fwrite(df_out, output, col.names = T, row.names = F, quote = F, sep = "\t", na = "NA")
cat("Finish!\n")
