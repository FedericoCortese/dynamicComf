# PER IL MOMENTO LE FEATURES SONO TUTTE CONTINUE

# weight_inv_exp_dist=function(Y,s,W,zeta){
#   #(yi,yj,si,sj,W,zeta)
#   # max_w=apply(rbind(W[si,],W[sj,]),2,max)
#   # -zeta*log(exp(-abs(yi-yj)/zeta)%*%max_w)
#   pairs <- combn(TT, 2)
#   
#   results <- apply(pairs, 2, function(idx) {
#     i <- idx[1]
#     j <- idx[2]
#     
#     yi <- Y[i, ]
#     yj <- Y[j, ]
#     si <- s[i]
#     sj <- s[j]
#     
#     diff <- abs(yi - yj)
#     max_w <- pmax(W[si, ], W[sj, ])
#     val <- -zeta * log(sum(exp(-diff / zeta) * max_w))
#     c(i, j, val)
#   })
#   
#   # Convert to a full symmetric matrix
#   mat <- matrix(0, TT, TT)
#   idx <- results[1:2, ]
#   vals <- results[3, ]
#   mat[cbind(idx[1, ], idx[2, ])] <- vals
#   mat[cbind(idx[2, ], idx[1, ])] <- vals
#   return(mat)
# }

weight_inv_exp_dist_old <- function(Y, s, W, zeta) {
  TT <- nrow(Y)
  P <- ncol(Y)
  
  # Range per ogni feature (colonna)
  range_Y <- apply(Y, 2, function(col) max(col) - min(col))
  range_Y[range_Y == 0] <- 1  # evita divisione per zero
  
  # Tutte le combinazioni (i < j)
  pairs <- combn(TT, 2)
  i_idx <- pairs[1, ]
  j_idx <- pairs[2, ]
  
  Yi <- Y[i_idx, , drop = FALSE]
  Yj <- Y[j_idx, , drop = FALSE]
  si <- s[i_idx]
  sj <- s[j_idx]
  
  # Distanze normalizzate tipo Gower
  diff <- abs(Yi - Yj) / matrix(range_Y, nrow = nrow(Yi), ncol = P, byrow = TRUE)
  
  # max_w vettorializzato
  max_w <- mapply(function(a, b) pmax(W[a, ], W[b, ]), si, sj, SIMPLIFY = "array")
  max_w <- t(max_w)  # n_pairs x P
  
  # Distanza finale
  dist_vals <- -zeta * log(rowSums(exp(-diff / zeta) * max_w))
  
  # Matrice simmetrica
  mat <- matrix(0, TT, TT)
  mat[cbind(i_idx, j_idx)] <- dist_vals
  mat[cbind(j_idx, i_idx)] <- dist_vals
  return(mat)
}

weight_inv_exp_dist <- function(Y, s, W, zeta) {
  TT <- nrow(Y)
  P <- ncol(Y)
  
  # 1. Normalizzazione Gower
  range_Y <- apply(Y, 2, function(col) {
    r <- max(col) - min(col)
    if (r == 0) 1 else r
  })
  Y_scaled <- sweep(Y, 2, range_Y, FUN = "/")
  
  # 2. Genera indici delle coppie (i < j)
  pairs <- combn(TT, 2)
  i_idx <- pairs[1, ]
  j_idx <- pairs[2, ]
  n_pairs <- ncol(pairs)
  
  # 3. Estrai le righe corrispondenti
  Yi <- Y_scaled[i_idx, , drop = FALSE]
  Yj <- Y_scaled[j_idx, , drop = FALSE]
  diff <- abs(Yi - Yj)
  
  # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
  W_si <- W[s[i_idx], , drop = FALSE]
  W_sj <- W[s[j_idx], , drop = FALSE]
  max_w <- pmax(W_si, W_sj)
  
  # 5. Calcola la distanza finale
  weighted_exp <- exp(-diff / zeta) * max_w
  dist_vals <- -zeta * log(rowSums(weighted_exp))
  
  # 6. Ricostruzione della matrice simmetrica
  mat <- matrix(0, TT, TT)
  mat[cbind(i_idx, j_idx)] <- dist_vals
  mat[cbind(j_idx, i_idx)] <- dist_vals
  
  return(mat)
}


WCD=function(s,Y,K){
  Tk=table(s)
  wcd=matrix(0,nrow=K,ncol=ncol(Y))
  for(k in unique(s)){
    for(p in 1:P){
      temp=Y[s==k,p]
      dist_temp=abs(outer(temp, temp, "-"))
      dist_temp[upper.tri(dist_temp)] <- 0
      dist_temp=dist_temp/diff(range(temp))
      wcd[k,p]=sum(dist_temp)/Tk[k]
    }
  }
  return(wcd)
}

# within_cluster_deviance <- function(Y, s, K) {
#   P <- ncol(Y)
#   
#   # Funzione per devianza cluster-feature
#   cluster_feature_deviance <- function(k) {
#     rows_k <- which(s == k)
#     if (length(rows_k) <= 1) {
#       return(rep(0, P))
#     }
#     Yk <- Y[rows_k, , drop = FALSE]
#     apply(Yk, 2, function(col_vals) sum(abs(outer(col_vals, col_vals, "-"))))
#   }
#   
#   # Applica a tutti i cluster
#   dev_list <- lapply(1:K, cluster_feature_deviance)
#   do.call(rbind, dev_list)
# }

initialize_states <- function(Y, K) {
  n <- nrow(Y)
  
  ### Repeat the following few times?
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(cluster::daisy(Y, metric = "gower"))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  
  # Faster solution 
  dist_matrix <- StatMatch::gower.dist(Y, centroids)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}


sparse_robust_fit=function(Y,K,lambda,zeta0,alpha=.1,verbose=F,tol=1e-16,
                           n_init=10,n_inner=10,n_outer=10){
  
  P=ncol(Y)
  TT=nrow(Y)
  
  best_loss <- NULL
  best_s <- NULL
  best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  Gamma <- lambda * (1 - diag(K))
  
    zeta=zeta0
    
    for (outer in 1:n_outer){
      
      for( init in 1:n_init){
      s=initialize_states(Y,K)
      # Compute distances as in (33)-(34) of Friedman(2004) or (14) of Kampert (2017)

      loss_old <- 1e10
      for(inner in 1:n_inner){
        DW=weight_inv_exp_dist(Y,s,W,zeta)
        
        # JM clustering
        # M-Step
        medoids=cluster::pam(x=DW,k=K,diss=TRUE)
        Ymedoids=Y[medoids$medoids,]
        
        # E-Step
        
        #dist_Y_to_medoids <- proxy::dist(Y, Ymedoids, method = "Manhattan")
        dist_Y_to_medoids <- StatMatch::gower.dist(Y, Ymedoids)
        loss_by_state <- as.matrix(dist_Y_to_medoids)
        
        V <- loss_by_state
        for (t in (TT-1):1) {
          V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
        }
        s_old=s
        s[1] <- which.min(V[1,])
        for (t in 2:TT) {
          s[t] <- which.min(V[t,] + Gamma[s[t-1],])
        }
        loss <- min(V[1,])
        # if (length(unique(s)) < K) {
        #   break
        # }
        
        if (length(unique(s)) < K) {
          s=s_old
          break
        }
        
        if (verbose) {
          cat(sprintf('Init %d | Inner iteration %d: %.6e\n', init, inner, loss))
        }
        if (!is.null(tol)) {
          epsilon <- loss_old - loss
          if (epsilon < tol) {
            break
          }
        } else if (all(s == s_old)) {
          break
        }
        loss_old <- loss
        
      }
      
      if (is.null(best_s) || (loss_old < best_loss)) {
        best_loss <- loss_old
        best_s <- s
        best_medoids = Ymedoids
      }
      }
      wcd=exp(-WCD(s,Y,K)/zeta0)
      W=wcd/rowSums(wcd)
      
      eps_W=mean((W-W_old)^2)
      
      if (!is.null(tol)) {
        if (eps_W < tol) {
          break
        }
      }
      
      W_old=W
      zeta=zeta+alpha*zeta0
      
      if (verbose) {
        cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
      }
      
    }
    # if (is.null(best_s) || is.null(best_W) || (loss_old < best_loss)) {
    #   best_loss <- loss_old
    #   best_s <- s
    #   best_medoids = Ymedoids
    #   best_W = W
    #   best_init=init
    # }
  
  
  return(list(W=W,best_s=best_s,best_medoids=best_medoids))
}

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

apply_noise_by_cluster <- function(Y, s, feat_list) {
  Y_noised <- Y
  K <- length(feat_list)
  
  for (k in 1:K) {
    cluster_rows <- which(s == k)
    if (length(cluster_rows) <= 1) next  # nulla da mescolare se solo una riga
    all_features <- seq_len(ncol(Y))
    irrelevant_feats <- setdiff(all_features, feat_list[[k]])
    
    for (j in irrelevant_feats) {
      # mescola le osservazioni della colonna j solo tra le righe del cluster k
      Y_noised[cluster_rows, j] <- sample(Y[, j],size=length(Y_noised[cluster_rows, j]))
    }
  }
  
  return(Y_noised)
}


