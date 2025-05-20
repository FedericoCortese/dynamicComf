library(cluster)
library(StatMatch)
library(poliscidata) # for weighted mode 
library(pdfCluster) # for ARI
library(parallel)
library(foreach)
library(doParallel)

hellinger_distance_vec <- function(p, q) {
  sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)
}


compute_entropy <- function(prob_matrix, base = exp(1)) {
  # Ensure the input is a matrix
  prob_matrix <- as.matrix(prob_matrix)
  
  # Replace zeros with NA to avoid log(0); these will be treated as zero in entropy
  prob_matrix[prob_matrix == 0] <- NA
  
  # Compute entropy for each row
  entropy_values <- apply(prob_matrix, 1, function(row) {
    # Remove NA values (originally zeros)
    valid_probs <- row[!is.na(row)]
    -sum(valid_probs * log(valid_probs, base = base))
  })
  
  return(sum(entropy_values))
}

weighted_median <- function(x, weights) {
  # Ensure x and weights are of the same length
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length.")
  }
  
  # Sort x and weights by x
  sorted_indices <- order(x)
  x <- x[sorted_indices]
  weights <- weights[sorted_indices]
  
  # Compute cumulative weights
  cumulative_weights <- cumsum(weights)
  total_weight <- sum(weights)
  
  # Find the smallest x such that the cumulative weight is >= 50% of the total weight
  weighted_median <- x[which(cumulative_weights >= total_weight / 2)[1]]
  
  return(weighted_median)
}


weighted_mode <- function(x, weights) {
  # Ensure x and weights are of the same length
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length.")
  }
  
  # Aggregate weights for each unique value of x
  unique_x <- unique(x)
  aggregated_weights <- sapply(unique_x, function(val) sum(weights[x == val]))
  
  # Find the value of x with the maximum aggregated weight
  weighted_mode <- unique_x[which.max(aggregated_weights)]
  
  return(weighted_mode)
}


Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

get_cat=function(y,mc,mu_val,phi){
  # Function to simulate categorical data
  
  # Arguments:
  # y: continuous variable 
  # mc: Markov chain states
  # mu: numeric mean value
  # phi: conditional probability for the categorical outcome k in state k
  
  library(dplyr)
  
  mu_val=c(-mu_val,0,mu_val)
  #K=length(unique(mc))
  #mu=seq(-mu,mu,length.out=K)
  #phi1=(1-phi)/(K-1)
  phi1=(1-phi)/2
  
  TT=length(y)
  for(i in 1:TT){
    k=mc[i]
    switch(k,
           "1"={
             threshold=c(qnorm(phi1,mu_val[1]),qnorm(phi+phi1,mu_val[1]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=1
             }
             else if(y[i]<threshold[1]){
               y[i]=2
             }
             else{
               y[i]=3
             }
           },
           "2"={
             threshold=c(qnorm(phi1,mu_val[2]),qnorm(phi+phi1,mu_val[2]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=2
             }
             else if(y[i]<threshold[1]){
               y[i]=3
             }
             else{
               y[i]=1
             }
           },
           "3"={
             threshold=c(qnorm(phi1,mu_val[3]),qnorm(phi+phi1,mu_val[3]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=3
             }
             else if(y[i]<threshold[1]){
               y[i]=1
             }
             else{
               y[i]=2
             }
           }
    )
  }
  return(y)
  
}

order_states_condMean=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional mean for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMeans=sort(tapply(y,s,mean,na.rm=T))
  
  states_temp=match(s,names(condMeans))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

order_states_condMed=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional median for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMed=sort(tapply(y,s,median,na.rm=T))
  
  states_temp=match(s,names(condMed))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

punct=function(x,pNAs,typeNA){
  
  # x is a vector (column of the dataset)
  # pNAs is the percentage of missing values
  # typeNA is the type of missing values (0 for random, 1 for continuous, all other values will turn into no missing imputation)
  
  TT=length(x)
  pTT=round(TT*pNAs)
  if(typeNA==0){
    NAindx=sample(1:TT,pTT,replace = F)
    x[NAindx]=NA
  }
  else if(typeNA==1){
    NAindx=sample(1:(TT-pTT),1,replace = F)
    NAindx=seq(NAindx,NAindx+pTT)
    x[NAindx]=NA
  }
  
  return(x)
  
}

initialize_states <- function(Y, K) {
  n <- nrow(Y)
  
  ### Repeat the following few times?
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(daisy(Y, metric = "gower"))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  ###
  
  # init_stats=rep(0,n)
  # For cycle
  # for(i in 1:n){
  #   init_stats[i]=which.min(gower.dist(Y[i,],centroids))
  # }
  
  # Using sapply and vapply
  # init_stats2 <- sapply(1:n, function(i) which.min(gower.dist(Y[i,], centroids)))
  # init_stats3 <- vapply(1:n, function(i) which.min(gower.dist(Y[i,], centroids)), integer(1))
  
  # Faster solution 
  dist_matrix <- gower.dist(Y, centroids)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}
######

get_cat_t <- function(y, mc, mu_val, phi, df = 5) {
  # y: vettore continuo (una variabile)
  # mc: stati latenti (valori da 1 a K)
  # mu_val: vettore di K valori medi (es. c(-mu, 0, mu))
  # phi: probabilità condizionata
  # df: gradi di libertà della t-Student
  
  TT <- length(y)
  K <- length(mu_val)
  phi1 <- (1 - phi) / (K - 1)  # non serve direttamente, solo se vuoi usare probabilità
  
  for (i in 1:TT) {
    k <- mc[i]
    mu_k <- mu_val[k]
    
    # Calcola l'intervallo centrato in mu_k che contiene prob = phi sotto la t-Student
    half_width <- qt((1 + phi) / 2, df = df)
    lower <- mu_k - half_width
    upper <- mu_k + half_width
    
    # Verifica se y[i] cade nell'intervallo
    if (y[i] >= lower && y[i] <= upper) {
      y[i] <- k  # corretto: valore coerente con lo stato
    } else {
      # Scegli casualmente un altro k tra gli altri K-1
      y[i] <- sample(setdiff(1:K, k), 1)
    }
  }
  return(y)
}

sim_data_stud_t=function(seed=123,
                         TT,
                         P,
                         Pcat,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=4,
                         phi=.8,
                         pers=.95){
  
  
  MU=seq(mu, -mu, length.out=Ktrue)
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    # u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    u = mvtnorm::rmvt(TT, sigma = (nu-2)*Sigma/nu, df = nu, delta = rep(MU[k],P))
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  SimData=data.frame(SimData)
  
  if(!is.null(Pcat)){
    for (j in 1:Pcat) {
      SimData[, j] <- get_cat_t(SimData[, j], x, MU, phi=phi, df = nu)
      SimData[, j]=factor(SimData[, j],levels=1:Ktrue)
    }  
  }
  
  
  
  return(list(
    SimData=SimData,
    mchain=x,
    TT=TT,
    P=P,
    K=Ktrue,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}


simstud_fuzzyJM_coord=function(seed,
                               lambda,
                               TT,
                               P,
                               K,
                               mu=1.5,
                               phi=.8,
                               rho=0,
                               Pcat=NULL,
                               pers=.9,
                               m=1.5,
                               nu=4){
  # Simulate
  simDat=sim_data_stud_t(seed=123,
                         TT=TT,
                         P=P,
                         Pcat=Pcat,
                         Ktrue=K,
                         mu=mu,
                         rho=rho,
                         nu=nu,
                         phi=phi,
                         pers=pers)
  
  Y=simDat$SimData
  
  success <- FALSE
  trials=1
  while (!success&trials<10) {
    est <- try(fuzzy_jump_coord(Y=Y, 
                                K=K, 
                                lambda=lambda, 
                                m=m,
                                max_iter=10, 
                                n_init=10, tol=1e-8, 
                                verbose=FALSE
                                
    ), silent = TRUE)
    trials=trials+1
    
    if (!inherits(est, "try-error")) {
      success <- TRUE  # Exit the loop if no error
    } else {
      message("Retrying fuzzy_jump() due to an error...")
    }
  }
  
  
  
  old_MAP=factor(est$MAP,levels=1:K)
  MAP=factor(relabel_clusters(est$MAP,simDat$mchain),levels=1:K)
  #est$MAP=factor(relabel_clusters(est$MAP,simDat$mchain),levels=1:K)
  simDat$mchain=factor(simDat$mchain,levels=1:K)
  
  BAC=caret::confusionMatrix(MAP,simDat$mchain)$overall[1]
  ARI=adj.rand.index(MAP,simDat$mchain)
  
  tab <- table(MAP, old_MAP)
  new_order <- apply(tab, 1, which.max)
  best_S=est$best_S[, new_order]
  
  fRI=fclust::ARI.F(simDat$mchain,best_S)
  
  # Return
  return(list(
    S=best_S,
    MAP=MAP ,
    ARI=ARI,
    BAC=BAC,
    fRI=fRI,
    seed=seed,
    lambda=lambda,
    TT=TT,
    P=P,
    K=K,
    Pcat=Pcat,
    true_seq=simDat$mchain
  ))
  
}

######

fuzzy_jump_cpp <- function(Y, 
                             K, 
                             lambda=1e-5, 
                             m=1,
                             max_iter=5, 
                             n_init=10, tol=1e-16, 
                             verbose=FALSE
                             
) {
  # Fit jump model for mixed type data 
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # K: number of states
  # lambda: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  # time_vec is a vector of time points, needed if times are not equally sampled
  
  # Value:
  # best_s: estimated state sequence
  # Y: imputed data
  # Y.orig: original data
  # condMM: state-conditional medians and modes
  
  library(Rcpp)
  Rcpp::sourceCpp("simplex_pgd.cpp")
  
  K=as.integer(K)
  
  TT <- nrow(Y)
  P <- ncol(Y)
  best_loss <- NULL
  best_S <- NULL
  
  res_list <- lapply(1:n_init, function(init) {
    # k-prot++
    s <- initialize_states(Y, K)
    S <- matrix(0, nrow = TT, ncol = K)
    row_indices <- rep(1:TT)  
    S[cbind(row_indices, s)] <- 1 
    
    # inizializzazione mu
    mu <- sapply(1:K, function(k) {
      apply(Y, 2, function(x) poliscidata::wtd.median(x, weights = S[, k]^m))
    })
    mu <- t(mu)
    colnames(mu) <- colnames(Y)
    
    S_old <- S
    V <- gower.dist(Y, mu)
    loss_old <- sum(V * S^m) + lambda * sum(abs(S[1:(TT-1), ] - S[2:TT, ]))^2
    
    
    for (it in 1:max_iter) {
      Rcpp::sourceCpp("simplex_pgd.cpp")
      
      # S(1)
      S[1, ] <- optimize_pgd_1T(
        init   = rep(1/K, K),
        g      = V[1, ],
        s_t1   = S[2, ],
        lambda = lambda,
        m      = m
      )$par
      
      # S(2) to S(T-1)
      for (t in 2:(TT - 1)) {
        S[t, ] <- optimize_pgd_2T(
          init    = rep(1/K, K),
          g       = V[t, ],
          s_prec  = S[t - 1, ],
          s_succ  = S[t + 1, ],
          lambda  = lambda,
          m       = m
        )$par
      }
      
      # S(T)
      S[TT, ] <- optimize_pgd_1T(
        init   = rep(1/K, K),
        g      = V[TT, ],
        s_t1   = S[TT - 1, ],
        lambda = lambda,
        m      = m
      )$par
      
      # aggiorna mu
      mu <- sapply(1:K, function(k) {
        apply(Y, 2, function(x) poliscidata::wtd.median(x, weights = S[, k]^m))
      })
      mu <- t(mu)
      colnames(mu) <- colnames(Y)
      
      V <- gower.dist(Y, mu)
      loss <- sum(V * S^m) + lambda * sum(abs(S[1:(TT-1), ] - S[2:TT, ]))^2
      
      if (verbose) cat(sprintf('Iteration %d: %.6e\n', it, loss))
      
      if (!is.null(tol)) {
        epsilon <- loss_old - loss
        if (epsilon < tol) break
      } else if (all(S == S_old)) {
        break
      }
      loss_old <- loss
      S_old <- S
    }
    
    list(S = S, loss = loss_old)
  })
  
  best_idx <- which.min(sapply(res_list, function(x) x$loss))
  best_S <- res_list[[best_idx]]$S
  best_loss <- res_list[[best_idx]]$loss
  
  MAP=apply(best_S, 1, which.max)

  old_MAP=apply(best_S,1,which.max)
  MAP=order_states_condMed(Y[,1],old_MAP)
  
  tab <- table(MAP, old_MAP)
  new_order <- apply(tab, 1, which.max)
  
  # Reorder the columns of S accordingly
  best_S <- best_S[, new_order]
  
  
  return(list(best_S=best_S,
              loss=best_loss,
              MAP=MAP,
              Y=Y
  ))
}

library(Rcpp)        # for your optimize_pgd_1T / _2T
library(foreach)
library(doParallel)
# (also poliscidata, gower, etc.)

# make sure your pgd functions are already loaded:
# Rcpp::sourceCpp("simplex_pgd.cpp")

fuzzy_jump_cpp_parallel <- function(Y, 
                                    K, 
                                    lambda   = 1e-5, 
                                    m        = 1,
                                    max_iter = 5, 
                                    n_init   = 10, 
                                    tol      = 1e-16, 
                                    verbose  = FALSE) {
  # 1) Compile/load C++ solvers once in the master
  Rcpp::sourceCpp("simplex_pgd.cpp")
  
  # 2) Required R packages & functions in master
  require(poliscidata); require(gower)
  # (assume initialize_states() and order_states_condMed() are in your GlobalEnv)
  
  # 3) Launch cluster
  ncores <- max(1, parallel::detectCores() - 1)
  cl     <- parallel::makeCluster(ncores)
  
  # 4) On *each* worker: load Rcpp, source the same C++ file, load packages
  parallel::clusterEvalQ(cl, {
    library(Rcpp)
    sourceCpp("simplex_pgd.cpp")
    library(poliscidata)
    library(gower)
  })
  
  # 5) Export any pure-R helpers into the workers
  parallel::clusterExport(cl, c("initialize_states", "order_states_condMed"), envir = environment())
  
  # 6) Register and run foreach
  doParallel::registerDoParallel(cl)
  best <- foreach(init = seq_len(n_init),
                  .packages = c("poliscidata","gower","StatMatch","cluster"),
                  .combine  = function(a,b) if (a$loss <= b$loss) a else b
  ) %dopar% {
    TT <- nrow(Y); P <- ncol(Y); K <- as.integer(K)
    
    # init
    s <- initialize_states(Y, K)
    S <- matrix(0, TT, K)
    S[cbind(seq_len(TT), s)] <- 1
    
    mu <- matrix(NA, K, P)
    for (k in 1:K) {
      mu[k, ] <- apply(Y, 2, function(x) poliscidata::wtd.median(x, weights = S[,k]^m))
    }
    V        <- gower.dist(Y, mu)
    loss_old <- sum(V * S^m) + lambda * sum(abs(S[-TT,] - S[-1,]))^2
    
    # coordinate updates
    for (it in seq_len(max_iter)) {
      # t=1
      r1 <- optimize_pgd_1T(rep(1/K, K), V[1, ], S[2, ], lambda, m)
      S[1,] <- r1$par
      # t=2..TT-1
      for (t in 2:(TT-1)) {
        r2 <- optimize_pgd_2T(rep(1/K,K), V[t,], S[t-1,], S[t+1,], lambda, m)
        S[t,] <- r2$par
      }
      # t=TT
      rT <- optimize_pgd_1T(rep(1/K,K), V[TT,], S[TT-1,], lambda, m)
      S[TT,] <- rT$par
      
      # recompute
      for (k in 1:K) {
        mu[k, ] <- apply(Y, 2, function(x) poliscidata::wtd.median(x, weights = S[,k]^m))
      }
      V    <- gower.dist(Y, mu)
      loss <- sum(V * S^m) + lambda * sum(abs(S[-TT,] - S[-1,]))^2
      
      if (verbose) message(sprintf("init %d, iter %d: loss=%.6e", init, it, loss))
      if (loss_old - loss < tol) break
      loss_old <- loss
    }
    
    list(loss = loss_old, S = S)
  }
  
  parallel::stopCluster(cl)
  
  # pick & reorder
  loss=best$loss
  best_S  <- best$S
  oldMAP  <- apply(best_S,1,which.max)
  MAP     <- order_states_condMed(Y[,1], oldMAP)
  tab     <- table(MAP, oldMAP)
  new_ord <- apply(tab,1,which.max)
  best_S  <- best_S[, new_ord]
  
  list(best_S = best_S, MAP = MAP, Y = Y,loss=loss)
}

fuzzyJM_gap <- function(Y,
                     K_grid    = 2:3,
                     lambda_grid = seq(0.01, .5, .05),
                     m_grid     = seq(1.01,2,length.out=3),
                     tol       = NULL,
                     max_iter   = 10,
                     verbose   = FALSE,
                     n_cores   = NULL,
                     B         = 10,
                     n_init=10) {
  
  require(foreach)
  require(doParallel)
  require(dplyr)
  
  # 1) build your parameter grid
  grid <- expand.grid(K = K_grid, lambda = lambda_grid,m=m_grid, b = 0:B)
  
  # 2) set up cores
  if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  results_list <- foreach(i = seq_len(nrow(grid)),
                          .packages = c("poliscidata","Rcpp","StatMatch"),
                          .export   = c("Y","fuzzy_jump_cpp",
                                        "weighted_median",
                                        "order_states_condMed",
                                        "initialize_states",
                                        "tol","max_iter",
                                        "n_init"),
                          .errorhandling = "pass") %dopar% {
                            params <- grid[i, ]
                            K_val  <- params$K
                            lambda  <- params$lambda
                            b      <- params$b
                            m = params$m
                            
                            set.seed(1000*i + b)
                            
                            Y_input <- if (b == 0) Y else Y[sample(1:nrow(Y),nrow(Y)),]
                            permuted <- (b != 0)
                            
                            out <- tryCatch({
                              res <- fuzzy_jump_cpp(Y_input,
                                           lambda   = lambda,
                                           K       = K_val,
                                           m =m,
                                           tol     = tol,
                                           n_init = n_init,
                                           verbose = FALSE,
                                           max_iter  = max_iter)
                              
                              list(
                                meta = data.frame(K        = K_val,
                                                  lambda    = lambda,
                                                  m=m,
                                                  loss     = res$loss,
                                                  permuted = permuted,
                                                  i=i),
                                cosa = if (!permuted)
                                  list(K       = K_val,
                                       lambda   = lambda,
                                       m=m,
                                       s       = res$best_S,
                                       i=i
                                  )
                                else NULL,
                                error = NA_character_
                              )
                              
                            }, error = function(e) {
                              # on error, record the params and the message
                              list(
                                meta  = data.frame(K        = K_val,
                                                   lambda    = lambda,
                                                   loss     = NA_real_,
                                                   permuted = permuted),
                                cosa  = NULL,
                                error = e$message
                              )
                            })
                            out
                          }
  stopCluster(cl)
  
  # 4) pull out all the meta–data and cosa results
  meta_df <- do.call(rbind, lapply(results_list, `[[`, "meta"))
  # cosa_results <- Filter(Negate(is.null),
  #                        lapply(results_list, `[[`, "cosa"))
  
  # 5) inspect errors, if any
  errors <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(x$meta, error = x$error, stringsAsFactors = FALSE)
  })) %>% filter(!is.na(error))
  
  if (nrow(errors) && verbose) {
    message("Some iterations failed. Inspect 'errors' data frame.")
    print(errors)
  }
  
  # 6) compute GAP statistics
  gap_stats <- meta_df %>%
    group_by(K, lambda, m) %>%
    summarise(
      log_O           = log(loss[!permuted]),
      log_O_star_mean = mean(log(loss[permuted]), na.rm = TRUE),
      se_log_O_star   = sd(log(loss[permuted]),   na.rm = TRUE),
      GAP             = log_O_star_mean - log_O,
      .groups         = "drop"
    )
  
  list(
    gap_stats    = gap_stats,
    errors       = if (nrow(errors)) errors else NULL
  )
  
  
}

# simstud gaussian AR(1) ---------------------------------------------------------------------


simulate_fuzzy_mixture_mv <- function(
    TT = 1000,
    P = 2,
    mu = 1,
    Sigma_rho = 0.5,
    ar_rho = 0.9,
    tau = 0.5,
    seed = NULL
) {
  
  # The higher tau, the 'harder' the clustering
  
  if (!is.null(seed)) set.seed(seed)
  # carica MASS per mvrnorm
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required but not installed.")
  }
  
  # vettori di media
  mu1 <- rep(-mu, P)
  mu2 <- rep(mu, P)
  
  # matrice di covarianza comune
  Sigma <- matrix(Sigma_rho, nrow = P, ncol = P)
  diag(Sigma) <- 1
  
  # pre-allocazioni
  alpha <- numeric(TT)
  pi_t  <- numeric(TT)
  y_mat <- matrix(0, nrow = TT, ncol = P)
  
  # inizializza latente
  alpha[1] <- rnorm(1, 0, tau)
  pi_t[1]  <- pnorm(alpha[1])
  
  # simula AR(1) e pesi
  for (t in 2:TT) {
    alpha[t] <- ar_rho * alpha[t - 1] + rnorm(1, 0, tau)
    pi_t[t]  <- pnorm(alpha[t])
  }
  
  # estrai y_t dal mix
  for (t in 1:TT) {
    if (runif(1) < pi_t[t]) {
      y_mat[t, ] <- MASS::mvrnorm(1, mu1, Sigma)
    } else {
      y_mat[t, ] <- MASS::mvrnorm(1, mu2, Sigma)
    }
  }
  
  MAP=I(pi_t<.5)+1
  
  # restituisci data.frame
  df <- data.frame(time = 1:TT, 
                   as.data.frame(y_mat),
                   alpha = alpha, pi_1 = pi_t, MAP)
  
  
  names(df)[(1:P)+1] <- paste0("Y", 1:P)
  return(df)
}



