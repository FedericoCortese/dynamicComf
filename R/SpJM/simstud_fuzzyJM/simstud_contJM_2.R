library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

###

source("Utils_fuzzyJM_2.R")
sourceCpp("cont_jump.cpp")
# K=2 soft contJM ---------------------------------------------------------------------

max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

K_grid=2:3
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,2,length.out=5)
seed=1:50

TT=c(1000,2000)
P=c(5,10)

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

mu=1
Sigma_rho=0
ar_rho = 0.9
tau = .2
ncores=detectCores()-1

start=Sys.time()
res_list_soft_K2_cont <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  mi    <- hp$m[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  soft_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=cbind(soft_scen$pi_1,1-soft_scen$pi_1)
  
  Yinput=soft_scen[,(1:Pi)+1]
  
  attempt <- 1
  last_loss <- NA_real_
  
  repeat {
    result <- tryCatch({
      fit <- cont_jump(as.matrix(Yinput), 
                             Ki, 
                             jump_penalty= li, 
                             alpha=2,
                             initial_states=NULL,
                             max_iter=max_iter, 
                             n_init=n_init, 
                             tol=tol, 
                             mode_loss=T,
                             #method="euclidean",
                             grid_size=.05
      )
      list(success = TRUE, 
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
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
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start


# K=2 hard contJM ---------------------------------------------------------------------

max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

K_grid=2:3
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,2,length.out=5)
seed=1:50

TT=c(1000,2000)
P=c(5,10)

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

mu=1
Sigma_rho=0
ar_rho = 0.9
tau = 5
ncores=detectCores()-1


start=Sys.time()
res_list_hard_K2_cont <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  mi    <- hp$m[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  hard_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=cbind(hard_scen$pi_1,1-hard_scen$pi_1)
  
  Yinput=hard_scen[,(1:Pi)+1]
  
  attempt <- 1
  last_loss <- NA_real_
  
  repeat {
    result <- tryCatch({
      fit <- cont_jump(as.matrix(Yinput), 
                       Ki, 
                       jump_penalty= li, 
                       alpha=2,
                       initial_states=NULL,
                       max_iter=max_iter, 
                       n_init=n_init, 
                       tol=tol, 
                       mode_loss=T,
                       #method="euclidean",
                       grid_size=.05
      )
      list(success = TRUE, 
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
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
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start