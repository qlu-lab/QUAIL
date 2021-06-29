# Packages
suppressMessages(require(data.table))
suppressMessages(require(quantreg))
suppressMessages(require(parallel))
suppressMessages(library(optparse))

options(stringsAsFactors = F)
option_list = list(
  make_option("--pheno", action = "store", default = NA, type = "character", help = "The path of the phenotype file"),
  make_option("--covar", action = "store", default = NA, type = "character", help = "The path of the covariate file"),
  make_option("--output", action = "store", default = NA, type = "character", help = "The path of the output file"),
  make_option("--num_levels", action = "store", default = NA, type = "integer", help = "Number of quantile levels to fit"),
  make_option("--num_cores", action = "store", default = NA, type = "integer", help = "Number of cores for parellel computing")
)

opt = parse_args(OptionParser(option_list=option_list))

# --- 0. I/O check
if(length(opt) < 5){
  cat("ERROR: Missing essential inputs\n")
  q("no")
}

if(opt$num_levels %% 2 == 1) {
    cat("ERROR: The num_levels need to be a even number\n")
    q("no")
}

# Input parameters
phenotype <- opt$pheno
covariate <- opt$covar
output <- opt$output
num_levels <- opt$num_levels
num_cores <- opt$num_cores

# --- 1. Read the input files

# read the phenotype files
cat("Preparing the data:\n")
pheno <- as.data.frame(fread(phenotype,stringsAsFactors = F))
pheno <- subset(pheno,!duplicated(pheno$IID))
pheno[pheno[,3]==-9, 3] <- NA  # Plink format

# read the covariates files
covar <- as.data.frame(fread(covariate,data.table = F,stringsAsFactors = F))
covar <- subset(covar,!duplicated(covar$IID))

# align phenotype and covariates files
id_intersect <- intersect(pheno$IID, covar$IID)
pheno <- pheno[match(id_intersect, pheno$IID), ]
covar <- covar[match(id_intersect, covar$IID), ]

pheno_qr <- pheno[, 3]
covar_qr <- as.data.frame(covar[, 3:ncol(covar)])
covar_qr$d_rv <- rnorm(nrow(covar_qr))

# Progress bar of the mclapply; Source: https://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply/26892969#26892969
mclapply2 <- function(X, FUN, ..., 
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
    mc.cleanup = TRUE, mc.allow.recursive = TRUE,
    mc.progress=TRUE, mc.style=3) {
    if (!is.vector(X) || is.object(X)) X <- as.list(X)

    if (mc.progress) {
        f <- fifo(tempfile(), open="w+b", blocking=T)
        p <- parallel:::mcfork()
        pb <- txtProgressBar(0, length(X), style=mc.style)
        setTxtProgressBar(pb, 0) 
        progress <- 0
        if (inherits(p, "masterProcess")) {
            while (progress < length(X)) {
                readBin(f, "double")
                progress <- progress + 1
                setTxtProgressBar(pb, progress) 
            }
            cat("\n")
            parallel:::mcexit()
        }
    }
    tryCatch({
        result <- mclapply(X, function(...) {
                res <- FUN(...)
                if (mc.progress) writeBin(1, f)
                res
            }, 
            mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
            mc.silent = mc.silent, mc.cores = mc.cores,
            mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
        )

    }, finally = {
        if (mc.progress) close(f)
    })
    result
}

# --- 2. Obtain the quantile integrated rank score

# Step1: fit the null model for k quantiles
Get_a_i_tau <- function(i){
    tau_curr <- i/(num_levels + 1)
    # Fit quantile regression between phenotype and covarites at quantile i/k
    Qreg <- rq(pheno_qr ~ ., data = covar_qr, tau = tau_curr, method = "fn")
    # Constrcuct the quantile rank score for each quantile 
    coeff <-  summary(Qreg , se = "ker")$coefficients
    SE_d_rv <- coeff[nrow(coeff),2]
    a_i_tau <- (tau_curr - ifelse(residuals(Qreg) < 0 , 1, 0))*SE_d_rv/(sqrt(-tau_curr^2 + tau_curr))
    return(a_i_tau)
}

Fit_a_i_tau <- function(start = 1, end = num_levels){
    df_one_q <- mclapply2(start:end, Get_a_i_tau, mc.cores = num_cores)
    return(df_one_q)
}

cat("Fitting the quantile regression to obtain quantile rank score:\n")

# To suppress some Warnings that are harmless: https://github.com/HenrikBengtsson/future/issues/218
df_all_q <- suppressWarnings(Fit_a_i_tau(start=1, end =num_levels))

# Step2: Add the quantile rank score across quantiles
cat("Obtaining the integrated rank score\n")
int_rank_score <- 0
for (i in 1:num_levels){
    ## upper quantile weight = 1, lower quantile weight = -1
    weight <- ifelse( (i > (num_levels)/2), 1, -1)
    int_rank_score <- int_rank_score + weight*df_all_q[[i]]
}

int_rank_score <- int_rank_score/((num_levels)/2)
out_int_rank_score <- data.frame(FID = pheno$FID, IID = pheno$IID, int_rank_score = int_rank_score)

cat("Write out the integrated rank score\n")
fwrite(out_int_rank_score, output, col.names = T, row.names = F, quote = F, sep = "\t", na = "NA")
cat("Finish!")
