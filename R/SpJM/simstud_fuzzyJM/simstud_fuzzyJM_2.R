library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_fuzzyJM_2.R")

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

# K=2 soft ---------------------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = .2
ncores=detectCores()-1
Ktrue=2

start=Sys.time()
res_list_soft_K2 <- mclapply(seq_len(nrow(hp)), function(i) {
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
    K  = Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=soft_scen[,c("pi_1","pi_2")]
  
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
    TT = TTi,
    P = Pi,
    seed= seedi,
    mu=mu, 
    Sigma_rho=Sigma_rho,
    ar_rho=ar_rho,
    Ktrue=Ktrue,
    tau=tau,
    loss      = last_loss,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

save(res_list_soft_K2,file='res_list_soft_K2.Rdata')

lista_risultati <- lapply(res_list_soft_K2, function(el) {
  # estraggo S e ground_truth come vettori
  cS  <- el$S
  cGT <- el$ground_truth
  
  # calcolo la distanza di Hellinger
  hd <- hellinger_distance_matrix(cS, cGT)
  
  # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
  data.frame(
    jhellinger_dist = hd,
    lambda          = el$lambda,
    K               = el$K,
    m               = el$m,
    TT              = el$TT,
    P               = el$P,
    seed            = el$seed,
    loss            = el$loss,
    attempts        = el$attempts,
    stringsAsFactors = FALSE
  )
})

# 3. Combino tutti i data.frame in uno solo
results_df_soft_K2_fuzzy <- do.call(rbind, lista_risultati)
save(results_df_soft_K2_fuzzy,file='hellinger_df_soft_K2_fuzzy.Rdata')

# K=2 hard ----------------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = 5
Ktrue=2
ncores=detectCores()-1

start=Sys.time()
res_list_hard_K2 <- mclapply(seq_len(nrow(hp)), function(i) {
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
    K  = Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=hard_scen[,c("pi_1","pi_2")]
  
  Yinput=hard_scen[,(1:Pi)+1]
  
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
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
      last_loss <- result$loss
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
    TT = TTi,
    P = Pi,
    seed= seedi,
    mu=mu, 
    Ktrue=Ktrue,
    Sigma_rho=Sigma_rho,
    ar_rho=ar_rho,
    tau=tau,
    loss      = last_loss,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

lista_risultati <- lapply(res_list_hard_K2, function(el) {
  # estraggo S e ground_truth come vettori
  cS  <- el$S
  cGT <- el$ground_truth
  
  # calcolo la distanza di Hellinger
  hd <- hellinger_distance_matrix(cS, cGT)
  
  # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
  data.frame(
    jhellinger_dist = hd,
    lambda          = el$lambda,
    K               = el$K,
    m               = el$m,
    TT              = el$TT,
    P               = el$P,
    seed            = el$seed,
    loss            = el$loss,
    attempts        = el$attempts,
    stringsAsFactors = FALSE
  )
})

# 3. Combino tutti i data.frame in uno solo
results_df_hard_K2_fuzzy <- do.call(rbind, lista_risultati)
save(results_df_hard_K2_fuzzy,file='hellinger_df_hard_K2_fuzzy.Rdata')


# K=3 soft ----------------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = .2
ncores=detectCores()-1
Ktrue=3

start=Sys.time()
res_list_soft_K3 <- mclapply(seq_len(nrow(hp)), function(i) {
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
    K  = Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=soft_scen[,c("pi_1","pi_2","pi_3")]
  
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
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
      last_loss <- result$loss
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
    TT = TTi,
    P = Pi,
    seed= seedi,
    Ktrue=Ktrue,
    mu=mu, 
    Sigma_rho=Sigma_rho,
    ar_rho=ar_rho,
    tau=tau,
    loss      = last_loss,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

save(res_list_soft_K2,file='res_list_soft_K3.Rdata')

lista_risultati <- lapply(res_list_soft_K3, function(el) {
  # estraggo S e ground_truth come vettori
  cS  <- el$S
  cGT <- el$ground_truth
  
  # calcolo la distanza di Hellinger
  hd <- hellinger_distance_matrix(cS, cGT)
  
  # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
  data.frame(
    jhellinger_dist = hd,
    lambda          = el$lambda,
    K               = el$K,
    m               = el$m,
    TT              = el$TT,
    P               = el$P,
    seed            = el$seed,
    loss            = el$loss,
    attempts        = el$attempts,
    PE              = el$PE,
    PB              = el$PB,
    PB_lambda       = el$PB_lambda,
    XB              = el$XB,
    stringsAsFactors = FALSE
  )
})

# 3. Combino tutti i data.frame in uno solo
results_df_soft_K3_fuzzy <- do.call(rbind, lista_risultati)
save(results_df_soft_K3_fuzzy,file='hellinger_df_soft_K3_fuzzy.Rdata')


# K=3 hard ----------------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = 5
Ktrue=3
ncores=detectCores()-1

start=Sys.time()
res_list_hard_K3 <- mclapply(seq_len(nrow(hp)), function(i) {
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
    K  = Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=hard_scen[,c("pi_1","pi_2","pi_3")]
  
  Yinput=hard_scen[,(1:Pi)+1]
  
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
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
      last_loss <- result$loss
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
    TT = TTi,
    P = Pi,
    seed= seedi,
    mu=mu, 
    Ktrue=Ktrue,
    Sigma_rho=Sigma_rho,
    ar_rho=ar_rho,
    tau=tau,
    loss      = last_loss,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

lista_risultati <- lapply(res_list_hard_K3, function(el) {
  # estraggo S e ground_truth come vettori
  cS  <- el$S
  cGT <- el$ground_truth
  
  # calcolo la distanza di Hellinger
  hd <- hellinger_distance_matrix(cS, cGT)
  
  # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
  data.frame(
    jhellinger_dist = hd,
    lambda          = el$lambda,
    K               = el$K,
    m               = el$m,
    TT              = el$TT,
    P               = el$P,
    seed            = el$seed,
    loss            = el$loss,
    attempts        = el$attempts,
    PE              = el$PE,
    PB              = el$PB,
    PB_lambda       = el$PB_lambda,
    XB              = el$XB,
    stringsAsFactors = FALSE
  )
})

# 3. Combino tutti i data.frame in uno solo
results_df_hard_K3_fuzzy <- do.call(rbind, lista_risultati)
save(results_df_hard_K3_fuzzy,file='hellinger_df_hard_K3_fuzzy.Rdata')

