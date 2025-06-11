library(cluster)
library(StatMatch)
library(poliscidata) # for weighted mode 
library(pdfCluster) # for ARI
library(parallel)
library(foreach)
library(doParallel)
library(gower)

hellinger_distance_vec <- function(p, q) {
  sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)
}

hellinger_distance_matrix=function(S,true_distr){
  if(!is.matrix(S)){
    # But it works only when K=2
    S=cbind(S,1-S)
  }
  hellinger_ts <- apply(cbind(S, true_distr), 1, function(row) {
    p <- row[1:2]
    q <- row[3:4]
    hellinger_distance_vec(p, q)
  })
  
  return(mean(hellinger_ts))
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

library(Rcpp)        # for your optimize_pgd_1T / _2T
library(foreach)
library(doParallel)

fuzzy_jump_cpp <- function(Y, 
                           K, 
                           lambda = 1e-5, 
                           m = 1,
                           max_iter = 5, 
                           n_init = 10, 
                           tol = 1e-16, 
                           verbose = FALSE,
                           parallel = FALSE,
                           n_cores = NULL
) {
  # Fit jump model for mixed‐type data with optional parallel initializations
  #
  # Arguments:
  #   Y            - data.frame with mixed data types (categorical vars must be factors)
  #   K            - number of states
  #   lambda       - penalty for the number of jumps
  #   m            - fuzziness exponent (for soft membership)
  #   max_iter     - maximum number of iterations per initialization
  #   n_init       - number of random initializations
  #   tol          - convergence tolerance
  #   verbose      - print progress per iteration (TRUE/FALSE)
  #   parallel     - if TRUE, use mclapply for parallel initializations
  #   n_cores      - number of cores for mclapply; if NULL, uses detectCores() - 1
  #
  # Value:
  #   List with:
  #     best_S      - TT×K soft‐membership matrix from best initialization
  #     loss        - best objective value
  #     MAP         - re‐ordered state sequence (1..K)
  #     Y           - imputed data (NA replaced)
  #     PE          - partitioning entropy
  #     PB          - PB index
  #     PB_lambda   - PB index using loss in denominator
  #     XB          - Xie–Beni index
  
  K <- as.integer(K)
  TT <- nrow(Y)
  P  <- ncol(Y)
  
  # Replace missing values with column medians initially
  for (j in seq_len(P)) {
    if (any(is.na(Y[[j]]))) {
      medj <- median(Y[[j]], na.rm = TRUE)
      Y[[j]][is.na(Y[[j]])] <- medj
    }
  }
  
  run_one <- function(init) {
    # Single initialization
    # 1) Initialize hard states via k‐prototypes++ (or custom routine)
    s <- initialize_states(Y, K)
    S <- matrix(0, nrow = TT, ncol = K)
    row_idx <- seq_len(TT)
    S[cbind(row_idx, s)] <- 1
    
    # 2) Initialize mu: state‐conditional medians/modes
    mu <- sapply(seq_len(K), function(k) {
      apply(Y, 2, function(x) poliscidata::wtd.median(x, weights = S[, k]^m))
    })
    mu <- t(mu)
    colnames(mu) <- colnames(Y)
    
    S_old <- S
    V <- gower.dist(Y, mu)
    loss_old <- sum(V * (S^m)) + lambda * sum(abs(S[1:(TT-1), ] - S[2:TT, ])^2)
    
    for (it in seq_len(max_iter)) {
      # Source Rcpp optimization routines (they must be in PATH)
      Rcpp::sourceCpp("simplex_pgd.cpp")
      
      # Update S[1, ]
      S[1, ] <- optimize_pgd_1T(
        init   = rep(1/K, K),
        g      = V[1, ],
        s_t1   = S[2, ],
        lambda = lambda,
        m      = m
      )$par
      
      # Update S[2:(TT-1), ]
      if (TT > 2) {
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
      }
      
      # Update S[TT, ]
      S[TT, ] <- optimize_pgd_1T(
        init   = rep(1/K, K),
        g      = V[TT, ],
        s_t1   = S[TT - 1, ],
        lambda = lambda,
        m      = m
      )$par
      
      # Recompute mu
      mu <- sapply(seq_len(K), function(k) {
        apply(Y, 2, function(x) poliscidata::wtd.median(x, weights = S[, k]^m))
      })
      mu <- t(mu)
      colnames(mu) <- colnames(Y)
      
      # Recompute distances and loss
      V <- gower.dist(Y, mu)
      loss <- sum(V * (S^m)) + lambda * sum(abs(S[1:(TT-1), ] - S[2:TT, ])^2)
      
      if (verbose) cat(sprintf("Initialization %d, Iteration %d: loss = %.6e\n", init, it, loss))
      
      # Check convergence
      if (!is.null(tol)) {
        if ((loss_old - loss) < tol) break
      } else if (all(S == S_old)) {
        break
      }
      loss_old <- loss
      S_old <- S
    }
    
    list(S = S, loss = loss_old, mu = mu)
  }
  
  # Choose apply function based on parallel flag
  if (parallel) {
    library(parallel)
    if (is.null(n_cores)) {
      n_cores <- max(detectCores() - 1, 1)
    }
    res_list <- mclapply(seq_len(n_init), run_one, mc.cores = n_cores)
  } else {
    res_list <- lapply(seq_len(n_init), run_one)
  }
  
  # Find best initialization
  losses <- vapply(res_list, function(x) x$loss, numeric(1))
  best_idx <- which.min(losses)
  best_run <- res_list[[best_idx]]
  best_S   <- best_run$S
  best_loss<- best_run$loss
  best_mu  <- best_run$mu
  
  # Compute MAP and re‐order states
  old_MAP <- apply(best_S, 1, which.max)
  MAP <- order_states_condMed(Y[, 1], old_MAP)
  tab <- table(MAP, old_MAP)
  new_order <- apply(tab, 1, which.max)
  best_S <- best_S[, new_order]
  
  # Compute cluster validity indices
  PE <- compute_entropy(best_S, base = exp(1))
  
  unc_med <- apply(Y, 2, median)
  ref <- as.data.frame(t(unc_med))
  colnames(ref) <- colnames(Y)
  E1 <- sum(gower_dist(Y, ref))
  
  Dmat <- gower.dist(best_mu)
  DK <- sum(Dmat[lower.tri(Dmat)])
  
  Jm <- sum((gower.dist(Y, best_mu)) * (best_S^m))
  
  PB <- (DK * (1/K) * E1 / Jm)^2
  PB_lambda <- (DK * (1/K) * E1 / best_loss)^2
  
  XB <- Jm / (TT * min(Dmat[lower.tri(Dmat)]))
  
  return(list(
    best_S    = best_S,
    loss      = best_loss,
    MAP       = MAP,
    Y         = Y,
    PE        = PE,
    PB        = PB,
    PB_lambda = PB_lambda,
    XB        = XB
  ))
}






# simstud gaussian AR(1) ---------------------------------------------------------------------


# simulate_fuzzy_mixture_mv <- function(
#     TT = 1000,
#     P = 2,
#     mu = 1,
#     Sigma_rho = 0.5,
#     ar_rho = 0.9,
#     tau = 0.5,
#     seed = NULL
# ) {
#   
#   # The higher tau, the 'harder' the clustering
#   
#   if (!is.null(seed)) set.seed(seed)
#   # carica MASS per mvrnorm
#   if (!requireNamespace("MASS", quietly = TRUE)) {
#     stop("Package 'MASS' is required but not installed.")
#   }
#   
#   # vettori di media
#   mu1 <- rep(-mu, P)
#   mu2 <- rep(mu, P)
#   
#   # matrice di covarianza comune
#   Sigma <- matrix(Sigma_rho, nrow = P, ncol = P)
#   diag(Sigma) <- 1
#   
#   # pre-allocazioni
#   alpha <- numeric(TT)
#   pi_t  <- numeric(TT)
#   y_mat <- matrix(0, nrow = TT, ncol = P)
#   
#   # inizializza latente
#   alpha[1] <- rnorm(1, 0, tau)
#   pi_t[1]  <- pnorm(alpha[1])
#   
#   # simula AR(1) e pesi
#   for (t in 2:TT) {
#     alpha[t] <- ar_rho * alpha[t - 1] + rnorm(1, 0, tau)
#     pi_t[t]  <- pnorm(alpha[t])
#   }
#   
#   # estrai y_t dal mix
#   for (t in 1:TT) {
#     if (runif(1) < pi_t[t]) {
#       y_mat[t, ] <- MASS::mvrnorm(1, mu1, Sigma)
#     } else {
#       y_mat[t, ] <- MASS::mvrnorm(1, mu2, Sigma)
#     }
#   }
#   
#   MAP=I(pi_t<.5)+1
#   
#   # restituisci data.frame
#   df <- data.frame(time = 1:TT, 
#                    as.data.frame(y_mat),
#                    alpha = alpha, pi_1 = pi_t, MAP)
#   
#   
#   names(df)[(1:P)+1] <- paste0("Y", 1:P)
#   return(df)
# }

simulate_fuzzy_mixture_mv<- function(
    TT = 1000,
    P  = 2,
    K  = 2,
    mu = 1,
    Sigma_rho = 0,
    ar_rho    = 0.9,
    tau       = 0.5,
    seed      = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required but not installed.")
  }
  
  # --- Generate K centroids spaced from -mu to +mu ---
  mu_vals <- seq(-mu, mu, length.out = K)
  mus     <- matrix(mu_vals, nrow = K, ncol = P, byrow = FALSE)
  
  # --- Common covariance ---
  Sigma <- matrix(Sigma_rho, nrow = P, ncol = P)
  diag(Sigma) <- 1
  
  # --- Storage ---
  alpha_mat <- matrix(0, nrow = TT, ncol = K)
  pi_mat    <- matrix(0, nrow = TT, ncol = K)
  y_mat     <- matrix(0, nrow = TT, ncol = P)
  
  # --- Initialize latent scores ---
  alpha_mat[1, ] <- rnorm(K, 0, tau)
  pi_mat[1, ]    <- exp(alpha_mat[1, ]) / sum(exp(alpha_mat[1, ]))
  
  # --- Simulate AR(1) + softmax weights ---
  for (t in 2:TT) {
    alpha_mat[t, ] <- ar_rho * alpha_mat[t - 1, ] + rnorm(K, 0, tau)
    ealpha <- exp(alpha_mat[t, ])
    pi_mat[t, ] <- ealpha / sum(ealpha)
  }
  
  # --- Draw from mixture ---
  for (t in 1:TT) {
    k_t       <- which.max(rmultinom(1, size = 1, prob = pi_mat[t, ]))
    y_mat[t,] <- MASS::mvrnorm(1, mu = mus[k_t, ], Sigma = Sigma)
  }
  
  # --- Build output ---
  MAP <- max.col(pi_mat)
  df <- data.frame(
    time      = seq_len(TT),
    y_mat,
    alpha_mat,
    pi_mat,
    MAP       = MAP
  )
  
  # Name columns
  names(df)[2:(1+P)]          <- paste0("Y",   1:P)
  names(df)[(2+P):(1+P+K)]    <- paste0("alpha_",1:K)
  names(df)[(2+P+K):(1+P+2*K)]<- paste0("pi_",   1:K)
  
  return(df)
}


# Cross-val ---------------------------------------------------------------

cv_sparse_jump <- function(
    Y,
    true_states,
    K_grid=NULL,
    m_grid=NULL,
    lambda_grid=NULL,
    n_folds = 5,
    parallel=F,
    n_cores=NULL
) {
  
  # cv_sparse_jump: Cross-validate Sparse Jump Model parameters (K, kappa, lambda)
  
  # Arguments:
  #   Y           - data matrix (N × P)
  #   true_states - vector of true states for ARI computation
  #   K_grid      - vector of candidate numbers of states
  #   m_grid  - vector of candidate fuzzyness parameters
  #   lambda_grid - vector of candidate lambdas
  #   n_folds     - number of folds for cross-validation (default: 5)
  #   parallel    - logical; TRUE for parallel execution (default: FALSE)
  
  # Value:
  #   A data.frame with one row per (K, kappa, lambda) combination, containing:
  #     K      - number of states tested
  #     m  - fuzzyness hyperparameter value
  #     lambda - jump penalty value
  #     fARI    - mean fuzzy Adjusted Rand Index across folds
  
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  if(is.null(m_grid)) {
    m_grid <- c(1.01,1.25,1.5,1.75,2)  # Default range for kappa
  }
  if(is.null(lambda_grid)) {
    lambda_grid <- seq(0,1,.1)  # Default range for lambda
  }
  
  # Libreria per fARI
  library(fclust)
  
  N <- nrow(Y)
  P <- ncol(Y)
  
  # Suddivido gli N campioni in n_folds blocchi contigui
  fold_indices <- split(
    1:N,
    rep(1:n_folds, each = ceiling(N / n_folds), length.out = N)
  )
  
  # Funzione che, per una tripla (K, m, lambda) e un fold (train_idx, val_idx),
  # calcola l’fARI sui punti di validazione
  fold_fari <- function(K, m, lambda, train_idx, val_idx) {
    # 1) Fit del modello sparse_jump su soli dati di TRAIN
    res <- fuzzy_jump_cpp(Y=as.matrix(Y[train_idx, , drop = FALSE]), 
                     K=as.integer(K), 
                     lambda=lambda, 
                     m=m,
                     max_iter=10, 
                     n_init=10, tol=1e-8, 
                     verbose=F
                     
      )
    
    states_train <- res$best_S
    
    
    # 2) Calcolo dei centroidi (medi anche solo sulle feature selezionate)
    #    Ogni riga di "centroids" è il centro per uno stato k = 1..K
    centroids <- res$mu
    
    # 3) A ciascun punto assegna la probabilità di appartenere a ciascuno stato
    #    (usando il reciproco della distanza di Gower tra il punto e i centroidi)
    Y_val_feats     <- as.matrix(Y[val_idx, feat_idx, drop = FALSE])
    pred_val_states <- matrix(0,nrow=length(val_idx),ncol=K)
    
    for (i in seq_along(val_idx)) {
      dists <- rep(Inf, K)
      for (k in 1:K) {
        if (!any(is.na(centroids[k, ]))) {
          # Essentially fuzzy cmeans distance
          dists[k]= gower.dist(Y_val_feats[i, ], centroids[k, ])^(-2/(m-1))
          #dists[k] <- sum((Y_val_feats[i, ] - centroids[k, ])^2)
        }
      }
      pred_val_states[i,] <- dists/sum(dists)
    }
    
    # ARI.F(VC=true_states[val_idx],U=clust$U)
    
    # 4) Calcolo ARI tra etichette vere e quelle predette sul blocco VAL
    return(adjustedRandIndex(true_states[val_idx], pred_val_states))
  }
  
  # Costruisco la griglia di tutte le combinazioni di (K, kappa, lambda)
  grid <- expand.grid(
    K      = K_grid,
    m  = m_grid,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Data.frame in cui raccogliere (K, kappa, lambda, media_ARI)
  results <- data.frame(
    K      = integer(0),
    m  = integer(0),
    lambda = numeric(0),
    fARI    = numeric(0)
  )
  
  # Loop su ciascuna riga della griglia
  if (!parallel) {
    # applichiamo una funzione su ogni riga di 'grid'
    results_list <- lapply(seq_len(nrow(grid)), function(row) {
      K_val     <- as.integer(grid$K[row])
      m_val <- as.integer(grid$m[row])
      lambda_val<-        grid$lambda[row]
      
      # calcolo ARI su ciascun fold
      fari_vals <- numeric(n_folds)
      for (f in seq_len(n_folds)) {
        val_idx   <- fold_indices[[f]]
        train_idx <- setdiff(seq_len(N), val_idx)
        fari_vals[f] <- fold_fari(K_val, m_val, lambda_val, train_idx, val_idx)
      }
      mean_fari <- mean(fari_vals)
      
      # ritorno un data.frame di una sola riga
      data.frame(
        K      = K_val,
        m  = m_val,
        lambda = lambda_val,
        fARI    = mean_fari,
        stringsAsFactors = FALSE
      )
    })
    
    # combino tutti i data.frame in un unico data.frame
    results <- do.call(rbind, results_list)
  }
  
  
  # 2) VERSIONE PARALLELA: usare mclapply()
  if (parallel) {
    if(is.null(n_cores)){
      n_cores <- detectCores() - 1
    }
    
    results_list <- mclapply(
      seq_len(nrow(grid)),
      function(row) {
        K_val     <- as.integer(grid$K[row])
        m_val <- as.integer(grid$m[row])
        lambda_val<-        grid$lambda[row]
        
        # calcolo ARI su ciascun fold
        fari_vals <- numeric(n_folds)
        for (f in seq_len(n_folds)) {
          val_idx   <- fold_indices[[f]]
          train_idx <- setdiff(seq_len(N), val_idx)
          fari_vals[f] <- fold_ari(K_val, m_val, lambda_val, train_idx, val_idx)
        }
        mean_fari <- mean(fari_vals)
        
        # ritorno un data.frame di una sola riga
        data.frame(
          K      = K_val,
          m  = m_val,
          lambda = lambda_val,
          fARI    = mean_fari,
          stringsAsFactors = FALSE
        )
      },
      mc.cores = n_cores
    )
    
    # combino tutti i data.frame in un unico data.frame
    results <- do.call(rbind, results_list)
  }
  
  return(results)
}

# continuous JM -----------------------------------------------------------

# See "cont_jump.cpp"
# sourceCpp("cont_jump.cpp")

# discretize_prob_simplex <- function(n_c, grid_size) {
#   # Sample grid points on a probability simplex.
#   N <- as.integer(1 / grid_size)
# 
#   # Generate all combinations and filter those that sum to N
#   tuples <- expand.grid(rep(list(0:N), n_c))
#   valid_tuples <- tuples[rowSums(tuples) == N, ]
# 
#   # Reverse the order and scale to get the simplex points
#   lst <- as.matrix(valid_tuples[nrow(valid_tuples):1, ]) / N
#   rownames(lst) <- NULL
#   return(lst)
# }
# 
# onerun_contJM=function(Y,K,
#                        jump_penalty,alpha,grid_size,mode_loss=T,
#                        max_iter,tol,initial_states,M){
#   tryCatch({
#     TT=nrow(Y)
#     loss_old <- 1e10
# 
#     S_old=matrix(0,nrow=TT,ncol=K)
# 
#     # State initialization through kmeans++
#     if (!is.null(initial_states)) {
#       s <- initial_states
#     } else {
#       #s <- initialize_states_jumpR(Y, K)
#       s=maotai::kmeanspp(Y,k=K)
#     }
#     mu=matrix(NA,nrow=K,ncol=P)
#     for (i in unique(s)) {
# 
#       if(sum(s==i)>1){
#         mu[i,] <- colMeans(Y[s==i,])
#       }
#       else{
#         mu[i,]=mean(Y[s==i,])
#       }
#     }
# 
#     for (it in 1:max_iter) {
# 
#       # E Step
# 
#       # Compute loss matrix
#       loss_mx <- matrix(NA, nrow(Y), nrow(mu))
#       for (t in 1:nrow(Y)) {
#         for (k in 1:nrow(mu)) {
#           loss_mx[t, k] <- .5*sqrt(sum((Y[t, ] - mu[k, ])^2))
#         }
#       }
# 
#       # Continuous model
#       prob_vecs <- discretize_prob_simplex(K, grid_size)
#       pairwise_l1_dist <- as.matrix(dist(prob_vecs, method = "manhattan")) / 2
#       jump_penalty_mx <- jump_penalty * (pairwise_l1_dist ^ alpha)
# 
#       if (mode_loss) {
#         # Adding mode loss
#         m_loss <- log(rowSums(exp(-jump_penalty_mx)))
#         m_loss <- m_loss - m_loss[1]  # Offset a constant
#         jump_penalty_mx <- jump_penalty_mx + m_loss
#       }
# 
#       LARGE_FLOAT=1e1000
#       # Handle continuous model with probability vectors
#       if (!is.null(prob_vecs)) {
#         loss_mx[is.nan(loss_mx)] <- LARGE_FLOAT
#         loss_mx <- loss_mx %*% t(prob_vecs)
#       }
# 
# 
#       N <- ncol(loss_mx)
# 
#       loss_mx[is.nan(loss_mx)] <- Inf
# 
#       # DP algorithm initialization
#       values <- matrix(NA, TT, N)
#       assign <- integer(TT)
# 
#       # Initial step
#       values[1, ] <- loss_mx[1, ]
# 
#       # DP iteration (bottleneck)
#       for (t in 2:TT) {
#         values[t, ] <- loss_mx[t, ] + apply(values[t - 1, ] + jump_penalty_mx, 2, min)
#       }
# 
#       S=matrix(NA,nrow=TT,ncol=K)
# 
#       # Find optimal path backwards
#       assign[TT] <- which.min(values[TT, ])
#       value_opt <- values[TT, assign[TT]]
# 
#       S[TT,]=prob_vecs[assign[TT],]
# 
#       # Traceback
#       for (t in (TT - 1):1) {
#         assign[t] <- which.min(values[t, ] + jump_penalty_mx[, assign[t + 1]])
#         S[t,]=prob_vecs[assign[t],]
#       }
# 
#       # M Step
# 
#       for(k in 1:K){
#         # What if the denominator is 0??
#         mu[k,]=apply(Y, 2, function(x) weighted.mean(x, w = S[,k]))
#       }
# 
#       # Re-fill-in missings
#       for(i in 1:P){
#         Y[,i]=ifelse(M[,i],mu[apply(S,1,which.max),i],Y[,i])
#       }
# 
#       if (!is.null(tol)) {
#         epsilon <- loss_old - value_opt
#         if (epsilon < tol) {
#           break
#         }
#       }
# 
#       else if (all(S == S_old)) {
#         break
#       }
# 
#       S_old=S
# 
#       loss_old <- value_opt
#     }
# 
# 
# 
#     return(list(S=S,value_opt=value_opt,mu=mu))}, error = function(e) {
#       # Return a consistent placeholder on error
#       return(list(S = NA, value_opt = Inf))
#     })
# }
# 
# cont_jump_R <- function(Y,
#                        K,
#                        jump_penalty=1e-5,
#                        alpha=2,
#                        initial_states=NULL,
#                        max_iter=10,
#                        n_init=10,
#                        tol=NULL,
#                        mode_loss=T,
#                        #method="euclidean",
#                        grid_size=.05,
#                        prll=F,
#                        n_cores=NULL
# ) {
#   # Fit jump model using framework of Bemporad et al. (2018)
# 
#   Y=as.matrix(Y)
#   K=as.integer(K)
# 
#   TT <- nrow(Y)
#   P <- ncol(Y)
# 
#   # Initialize mu
#   mu <- colMeans(Y,na.rm = T)
# 
#   # Track missings with 0 1 matrix
#   M=ifelse(is.na(Y),T,F)
# 
# 
#   Ytil=Y
#   # Impute missing values with mean of observed states
#   for(i in 1:P){
#     Y[,i]=ifelse(M[,i],mu[i],Y[,i])
#   }
# 
# 
#     jms=list()
#     for (init in 1:n_init) {
# 
#       jms[[init]]=onerun_contJM(Y,K,
#                                 jump_penalty,alpha,grid_size,mode_loss,
#                                 max_iter,tol,initial_states,M=M)
# 
#     }
# 
# 
#   best_init=which.min(unlist(lapply(jms,function(x)x$value_opt)))
#   best_S=jms[[best_init]]$S
#   best_loss=jms[[best_init]]$value_opt
#   best_mu=jms[[best_init]]$mu
# 
# 
#   MAP=apply(best_S, 1, which.max)
# 
#   old_MAP=apply(best_S,1,which.max)
#   MAP=order_states_condMed(Y[,1],old_MAP)
# 
#   tab <- table(MAP, old_MAP)
#   new_order <- apply(tab, 1, which.max)
# 
#   # Reorder the columns of S accordingly
#   best_S <- best_S[, new_order]
# 
#   # Cluster validity indexes
#   # Partitioning entropy
#   PE=compute_entropy(best_S, base = exp(1))
# 
#   # PB
#   # E1 component
#   unc_med=apply(Y,2,median)
#   ref <- as.data.frame(t(unc_med))
#   colnames(ref) <- colnames(Y)
#   E1 <- sum(gower_dist(Y, ref))
# 
#   # DK component
#   Dmat <- gower.dist(best_mu)
#   DK=sum(Dmat[lower.tri(Dmat)])
# 
#   # Jm component
#   Jm=sum(gower.dist(Y, best_mu) * best_S)
# 
#   PB=(DK*(1/K)*E1/Jm)^2
#   PB_lambda=(DK*(1/K)*E1/best_loss)^2
# 
#   # XB
#   XB=Jm/(TT*min(Dmat[lower.tri(Dmat)]))
# 
#   return(list(best_S=best_S,
#               best_loss=best_loss,
#               best_mu=best_mu,
#               XB=XB,
#               PB=PB,
#               PB_lambda=PB_lambda,
#               PE=PE))
# }
