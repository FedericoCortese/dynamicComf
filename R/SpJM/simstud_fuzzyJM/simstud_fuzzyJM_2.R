library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_fuzzyJM_2.R")

# K=2 soft ---------------------------------------------------------------------

max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

K_grid=2:3
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,2,length.out=5)

TT=c(1000,2000)
P=c(5,10)

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  TT=TT,
                  P=P)

mu=1
Sigma_rho=0
ar_rho = 0.9
tau = .2
ncores=detectCores()-1

start=Sys.time()
res_list_soft_K2 <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  mi    <- hp$m[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  
  soft_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = set.seed(i)
  )
  
  ground_truth=data.frame(p1=soft_scen$pi_1,p2=1-soft_scen$pi_1)
  
  Yinput=soft_scen[,(1:Pi)+1]
  
  attempt <- 1
  last_loss <- NA_real_
  
  repeat {
    result <- tryCatch({
      fit <- fuzzy_jump_cpp(
        Yinput,
        K        = Ki,
        lambda   = li,
        m        = mi,
        max_iter = max_iter,
        n_init   = n_init,
        tol      = tol,
        verbose  = FALSE
      )
      list(success = TRUE, loss = fit$loss,
           PE=fit$PE,PB=fit$PB,PB_lambda=fit$PB_lambda,XB=fit$XB,
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
      last_loss <- result$loss
      PE <- result$PE
      PB <- result$PB
      PB_lambda <- result$PB_lambda
      XB <- result$XB
      break
    }
    if (attempt >= max_retries) {
      warning(sprintf("Row %d: reached max attempts (%d), setting loss=NA", i, max_retries))
      break
    }
    attempt <- attempt + 1
    Sys.sleep(0.1)
  }
  
  list(
    S=S,
    ground_truth=ground_truth,
    K         = Ki,
    lambda    = li,
    m         = mi,
    loss      = last_loss,
    attempts  = attempt,
    PE=PE,
    PB=PB,
    PB_lambda=PB_lambda,
    XB=XB
  )
}, mc.cores = ncores)
end=Sys.time()
end-start


# K=2 hard ----------------------------------------------------------------



# K=3 ---------------------------------------------------------------------


