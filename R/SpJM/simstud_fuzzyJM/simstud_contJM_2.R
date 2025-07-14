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

max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

K_grid=2:3
#lambda_grid=c(0,.5,1,5,10,25,50,100)
lambda_grid=c(seq(0,.5,length.out=11),5,10,50,100)
seed=1:50

TT=c(1000,2000)
P=c(5,10)


hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

# K=2 soft cont JM --------------------------------------------------------



mu=1
Sigma_rho=0
ar_rho = 0.99
tau = .2
ncores=detectCores()-1
Ktrue=2

start=Sys.time()
res_list_soft_K2_cont <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  soft_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    K=Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=soft_scen[,c("pi_1", "pi_2")]
  
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
                             mode_loss=F,
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
      
      # Reorder states
      # Compute MAP and re‐order states
      old_MAP <- apply(result$S, 1, which.max)
      MAP <- order_states_condMed(Yinput[, 1], old_MAP)
      tab <- table(MAP, old_MAP)
      new_order <- apply(tab, 1, which.max)
      result$S <- result$S[, new_order]
      
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
    TT        = TTi,
    P         = Pi,
    seed      = seedi,
    mu         = mu,
    Sigma_rho = Sigma_rho,
    ar_rho    = ar_rho,
    tau       = tau,
    Ktrue     = Ktrue,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

save(res_list_soft_K2_cont,file='res_cont_soft_K2.Rdata')

# lista_risultati <- mclapply(res_list_soft_K2_cont, function(el) {
#   # estraggo S e ground_truth come vettori
#   cS  <- as.matrix(el$S)
#   cGT <- as.matrix(el$ground_truth)
#   
#   MAP_true=apply(cS,1,which.max)
#   MAP_est=apply(cGT,1,which.max)
#   confusion <-table(MAP_true, MAP_est)
#   mapping <- apply(confusion, 1, which.max)
#   cS=cS[,mapping]
#   
#   # calcolo la distanza di Hellinger
#   hd <- mean(vapply(
#     seq_len(nrow(cS)),
#     function(i) hellinger_distance_vec(cS[i, ], cGT[i, ]),
#     numeric(1)
#   ))
#   
#   # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
#   data.frame(
#     jhellinger_dist = hd,
#     lambda          = el$lambda,
#     K               = el$K,
#     stringsAsFactors = FALSE
#   )
# },mc.cores = ncores)
# 
# # 3. Combino tutti i data.frame in uno solo
# results_df_soft_K2_cont <- do.call(rbind, lista_risultati)
# 
# save(results_df_soft_K2_cont,file='hellinger_df_soft_K2_cont.Rdata')

# K=2 hard contJM ---------------------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = 5
ncores=detectCores()-1
Ktrue=2

start=Sys.time()
res_list_hard_K2_cont <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  hard_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    K=Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=hard_scen[,c("pi_1", "pi_2")]
  
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
                       mode_loss=F,
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
      
      # Reorder states
      # Compute MAP and re‐order states
      old_MAP <- apply(result$S, 1, which.max)
      MAP <- order_states_condMed(Yinput[, 1], old_MAP)
      tab <- table(MAP, old_MAP)
      new_order <- apply(tab, 1, which.max)
      result$S <- result$S[, new_order]
      
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
    TT        = TTi,
    P         = Pi,
    seed      = seedi,
    mu         = mu,
    Sigma_rho = Sigma_rho,
    ar_rho    = ar_rho,
    tau       = tau,
    Ktrue     = Ktrue,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

save(res_list_hard_K2_cont,file='res_cont_hard_K2.Rdata')

# lista_risultati <- mclapply(res_list_hard_K2_cont, function(el) {
#   # estraggo S e ground_truth come vettori
#   cS  <- as.matrix(el$S)
#   cGT <- as.matrix(el$ground_truth)
#   
#   MAP_true=apply(cS,1,which.max)
#   MAP_est=apply(cGT,1,which.max)
#   confusion <-table(MAP_true, MAP_est)
#   mapping <- apply(confusion, 1, which.max)
#   cS=cS[,mapping]
#   
#   # calcolo la distanza di Hellinger
#   hd <- mean(vapply(
#     seq_len(nrow(cS)),
#     function(i) hellinger_distance_vec(cS[i, ], cGT[i, ]),
#     numeric(1)
#   ))
#   
#   # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
#   data.frame(
#     jhellinger_dist = hd,
#     lambda          = el$lambda,
#     K               = el$K,
#     stringsAsFactors = FALSE
#   )
# },mc.cores = ncores)
# 
# # 3. Combino tutti i data.frame in uno solo
# results_df_hard_K2_cont <- do.call(rbind, lista_risultati)
# 
# save(results_df_hard_K2_cont,file='hellinger_df_hard_K2_cont.Rdata')


# K=3 soft cont JM --------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = .2
ncores=detectCores()-1
Ktrue=3

start=Sys.time()
res_list_soft_K3_cont <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  soft_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    K=Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=soft_scen[,c("pi_1", "pi_2", "pi_3")]
  
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
                       mode_loss=F,
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
      
      # Reorder states
      # Compute MAP and re‐order states
      old_MAP <- apply(result$S, 1, which.max)
      MAP <- order_states_condMed(Yinput[, 1], old_MAP)
      tab <- table(MAP, old_MAP)
      new_order <- apply(tab, 1, which.max)
      result$S <- result$S[, new_order]
      
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
    TT        = TTi,
    P         = Pi,
    seed      = seedi,
    mu         = mu,
    Sigma_rho = Sigma_rho,
    ar_rho    = ar_rho,
    tau       = tau,
    Ktrue     = Ktrue,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

save(res_list_soft_K3_cont,file='res_cont_soft_K3.Rdata')

# lista_risultati <- mclapply(res_list_soft_K3_cont, function(el) {
#   # estraggo S e ground_truth come vettori
#   cS  <- as.matrix(el$S)
#   cGT <- as.matrix(el$ground_truth)
#   
#   MAP_true=apply(cS,1,which.max)
#   MAP_est=apply(cGT,1,which.max)
#   confusion <-table(MAP_true, MAP_est)
#   mapping <- apply(confusion, 1, which.max)
#   cS=cS[,mapping]
#   
#   # calcolo la distanza di Hellinger
#   hd <- mean(vapply(
#     seq_len(nrow(cS)),
#     function(i) hellinger_distance_vec(cS[i, ], cGT[i, ]),
#     numeric(1)
#   ))
#   
#   # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
#   data.frame(
#     jhellinger_dist = hd,
#     lambda          = el$lambda,
#     K               = el$K,
#     stringsAsFactors = FALSE
#   )
# },mc.cores = ncores)
# 
# # 3. Combino tutti i data.frame in uno solo
# results_df_soft_K3_cont <- do.call(rbind, lista_risultati)
# 
# save(results_df_soft_K3_cont,file='hellinger_df_soft_K3_cont.Rdata')


# K= 3 hard cont ----------------------------------------------------------

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = 5
ncores=detectCores()-1
Ktrue=3

start=Sys.time()
res_list_hard_K3_cont <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  hard_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    K=Ktrue,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=hard_scen[,c("pi_1", "pi_2", "pi_3")]
  
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
                       mode_loss=F,
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
      
      # Reorder states
      # Compute MAP and re‐order states
      old_MAP <- apply(result$S, 1, which.max)
      MAP <- order_states_condMed(Yinput[, 1], old_MAP)
      tab <- table(MAP, old_MAP)
      new_order <- apply(tab, 1, which.max)
      result$S <- result$S[, new_order]
      
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
    TT        = TTi,
    P         = Pi,
    seed      = seedi,
    mu         = mu,
    Sigma_rho = Sigma_rho,
    ar_rho    = ar_rho,
    tau       = tau,
    Ktrue     = Ktrue,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

save(res_list_hard_K3_cont,file='res_cont_hard_K3.Rdata')

# lista_risultati <- mclapply(res_list_hard_K3_cont, function(el) {
#   # estraggo S e ground_truth come vettori
#   cS  <- as.matrix(el$S)
#   cGT <- as.matrix(el$ground_truth)
#   
#   MAP_true=apply(cS,1,which.max)
#   MAP_est=apply(cGT,1,which.max)
#   confusion <-table(MAP_true, MAP_est)
#   mapping <- apply(confusion, 1, which.max)
#   cS=cS[,mapping]
#   
#   # calcolo la distanza di Hellinger
#   hd <- mean(vapply(
#     seq_len(nrow(cS)),
#     function(i) hellinger_distance_vec(cS[i, ], cGT[i, ]),
#     numeric(1)
#   ))
#   
#   # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
#   data.frame(
#     jhellinger_dist = hd,
#     lambda          = el$lambda,
#     K               = el$K,
#     stringsAsFactors = FALSE
#   )
# },mc.cores = ncores)
# 
# # 3. Combino tutti i data.frame in uno solo
# results_df_hard_K3_cont <- do.call(rbind, lista_risultati)
# 
# save(results_df_hard_K3_cont,file='hellinger_df_hard_K3_cont.Rdata')

