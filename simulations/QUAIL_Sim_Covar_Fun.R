## Functions used in QUAIL 
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(dglm))
suppressPackageStartupMessages(library(quantreg))

numscale = function(x) {
  if(sum(x==0)==length(x)){x}else{(x - mean(x))/sd(x)}
}

Sim_Covar <- function(seed) {
  set.seed(42 * as.numeric(seed))
  # Simulate SNP
  N <- 10^3 # Sample size
  maf <- runif(1, 0.05, 0.5)
  geno <- rbinom(N, 2, maf)
  
  # Simulate error terms
  e <- rnorm(N) # Standard normal distribution
  
  # Simulate Environment
  E <- rnorm(N)
  
  # Simulate covariates
  covar <- rep(0, N)
  covar[which(geno == 0)] <- rbinom(length(which(geno == 0)), 1, 0.2)
  covar[which(geno == 1)] <- rbinom(length(which(geno == 1)), 1, 0.5)
  covar[which(geno == 2)] <- rbinom(length(which(geno == 2)), 1, 0.8)
  
  ##  Variance explianed by mean effects and variance efects
  rsq_covar_main <- 0.2
  rsq_covarxe <- 0.2
  rsq_e <- 1 - rsq_covar_main - rsq_covarxe
  
  # Simulate outcome
  covar_main <- numscale(covar) * rnorm(1, 0, 1)
  covarxe <- numscale(covar) * E  * rnorm(1, 0, 1)
  y <- numscale(covar_main) * sqrt(rsq_covar_main) +
    numscale(covarxe) * sqrt(rsq_covarxe) +
    numscale(e) * sqrt(rsq_e)
  return(list(geno = geno, y = y, covar = covar))
}


### INT
int = function(x){return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}

QUAIL_Covar <- function(x, y, covar, quantiles = 100){
  x <- (x - mean(x))/sd(x)
  covar <- (covar - mean(covar))/sd(covar)
  y <- lm(y~ x + covar)$residuals
  a_i_tau_list <- rep()
  k <- quantiles - 1 ### number of quantiles to fit the model
  d <- rnorm(length(x))
  for (i in 1:k){
    # print(i)
    tau_curr <- i/(k+1)
    Qreg <- rq(y ~ d + covar, tau = tau_curr, method = "fn")
    SE_d_cur <- summary(Qreg , se = "ker")$coefficients[2,2]
    balance <-   sqrt(-tau_curr^2 + tau_curr)/SE_d_cur    
    a_i_tau <- ifelse(residuals(Qreg) < 0 ,1,0)
    a_i_tau_list[[i]] <- (-a_i_tau + tau_curr)/balance
  }
  
  int_rank_score <- data.frame(score = rep(0, length(y)))
  for (i in 1:(k-1)){
    if ( i > quantiles/2){
      term2 <- 1
    }else {
      term2 <- -1
    }
    int_rank_score <- int_rank_score + term2*a_i_tau_list[[i]]
  }
  
  int_rank_score <- int_rank_score/(quantiles/2)
  Xstar <-  lm(x ~ covar)$residuals
  null_phi_1 <- int_rank_score$score*sqrt(sum(Xstar^2))
  p <- summary(lm(null_phi_1 ~ Xstar))$coefficients[2,4]
  return(p)
}

DRM_Covar <- function(x, y, covar){
  y <- lm(y~covar)$residuals
  x <- as.factor(x)
  y.i <- tapply(y, x, median)
  z.ij <- abs(y - y.i[x])
  res <- summary(lm(z.ij~x))$coef[2,]
  p <- res[4]
  return(p)
}

Lev_Covar <- function(x, y, covar){
  y <- lm(y~covar)$residuals
  P <- leveneTest(y~as.factor(x), center = median)$"P"[1]
  return(P)
}

HLMM_Covar <- function(x, y, covar){
  P <- coefficients(summary(dglm(y~ x + covar,dformula=~ x + covar))$"dispersion.summary")[2,4]
  return(P)
}

HLMM_INT_Covar <- function(x, y, covar){
  P <- coefficients(summary(dglm(int(y)~x + covar ,dformula=~ x + covar))$"dispersion.summary")[2,4]
  return(P)
}