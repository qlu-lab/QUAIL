## Functions used in QUAIL 
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(dglm))
suppressPackageStartupMessages(library(quantreg))

numscale = function(x) {
  if(sum(x==0)==length(x)){x}else{(x - mean(x))/sd(x)}
}

Sim <- function(seed, g_type = "mean"){
  set.seed(42 * as.numeric(seed))
  # Simulate SNP
  N <- 10^3 # Sample size
  maf <- runif(1, 0.05, 0.5)
  geno <- rbinom(N, 2, maf)
  
  # Simulate error terms
  e <- rt(N, 6) # Chisq distribution with df = 6
  e <- (e - 6)/sqrt(2*6) # standardized to have a mean of 0 and variance of 1
  
  # Simulate Environment
  E <- rnorm(N)
  
  # Simulate the effect
  ##  b_g is the mean effect, f_g is the variance effect
  if(g_type == "nei"){b_g = 0; f_g = 0
  }else if(g_type =="both"){b_g = rnorm(1,0,1); f_g = rnorm(1,0,1)
  }else if(g_type == "mean"){b_g = rnorm(1,0,1); f_g = 0
  }else if(g_type == "var"){b_g = 0; f_g =  rnorm(1,0,1)}
  
  ##  Variance explianed by mean effects and variance efects
  rsq_g <- 0.005
  rsq_gxe <- 0.005
  rsq_e <- 1 - rsq_g * ifelse(b_g^2 >0, 1, 0) - rsq_gxe * ifelse(f_g^2 >0, 1, 0)
  
  # Simulate outcome
  g = geno * b_g
  gxe = geno * E * f_g
  y = numscale(g)*sqrt(rsq_g) + numscale(gxe)*sqrt(rsq_gxe) + numscale(e)*sqrt(rsq_e)
  return(list(geno = geno, y = y))
}

### INT
int = function(x){return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}

QUAIL <- function(x, y, quantiles = 100){
  x <- (x - mean(x))/sd(x)
  y <- lm(y~x)$residuals
  a_i_tau_list <- rep()
  k <- quantiles - 1 ### number of quantiles to fit the model
  d <- rnorm(length(x))
  for (i in 1:k){
    # print(i)
    tau_curr <- i/(k+1)
    Qreg <- rq(y ~ d, tau = tau_curr, method = "fn")
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
  Xstar <-  lm(x ~ 1)$residuals
  null_phi_1 <- int_rank_score$score*sqrt(sum(Xstar^2))
  p <- summary(lm(null_phi_1 ~ Xstar))$coefficients[2,4]
  return(p)
}

DRM <- function(x, y){
  x <- as.factor(x)
  y.i <- tapply(y, x, median)
  z.ij <- abs(y - y.i[x])
  res <- summary(lm(z.ij~x))$coef[2,]
  p <- res[4]
  return(p)
}

Lev <- function(x, y){
  P <- leveneTest(y~as.factor(x), center = median)$"P"[1]
  return(P)
}

HLMM <- function(x, y){
  P <- coefficients(summary(dglm(y~x,dformula=~ x))$"dispersion.summary")[2,4]
  return(P)
}

HLMM_INT <- function(x, y){
  P <- coefficients(summary(dglm(int(y)~x ,dformula=~ x))$"dispersion.summary")[2,4]
  return(P)
}