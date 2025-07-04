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

# initialize_states <- function(Y, K, feat_type=NULL) {
#   
#   # centr: if TRUE, the functionoutputs initialized centroids
#   
#   TT <- nrow(Y)
#   P= ncol(Y)
#   
#   if(is.null(feat_type)){
#     feat_type=rep(0,P)
#   }
# 
#   ### Repeat the following few times?
#   centr_indx=sample(1:TT, 1)
#   centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
# 
#   #closest_dist_old <- as.matrix(cluster::daisy(Y, metric = "gower"))
#   closest_dist <- gower_dist(Y, Y,feat_type = feat_type,scale="m")
#   closest_dist <- closest_dist[centr_indx,]
# 
#   for (i in 2:K) {
#     prob <- closest_dist / sum(closest_dist)
#     next_centr_indx <- sample(1:TT, 1, prob = prob)
#     next_centroid <- Y[next_centr_indx, , drop = FALSE]
#     centroids <- rbind(centroids, next_centroid)
#   }
# 
#   # Faster solution
#   dist_matrix <- gower_dist(Y, centroids,feat_type = feat_type,scale="m")
#   init_stats <- apply(dist_matrix, 1, which.min)
#   
#   return(init_stats)
# 
# }

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

# Rcpp --------------------------------------------------------------------

library(Rcpp)
Rcpp::sourceCpp("robJM_R.cpp")


# OLD -------------------------------------------------------------------------

weight_inv_exp_dist_old <- function(Y,
                                    s, 
                                    W, zeta) {
  TT <- nrow(Y)
  P <- ncol(Y)
  
  sk=apply(Y,2,function(x)IQR(x)/1.35)
  
  # 2. Genera indici delle coppie (i < j)
  pairs <- combn(TT, 2)
  i_idx <- pairs[1, ]
  j_idx <- pairs[2, ]
  n_pairs <- ncol(pairs)
  
  # 3. Estrai le righe corrispondenti
  Yi <- Y[i_idx, , drop = FALSE]
  Yj <- Y[j_idx, , drop = FALSE]
  diff <- abs(Yi - Yj)
  #sk=apply(diff,2,function(x)sum(x)/TT^2)
  
  diff=sweep(diff, 2, sk, FUN = "/")
  #diff=sweep(diff, 2, range_Y, FUN = "/")
  
  
  # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
  W_si <- W[s[i_idx], , drop = FALSE]
  W_sj <- W[s[j_idx], , drop = FALSE]
  # W_si <- W[i_idx, , drop = FALSE]
  # W_sj <- W[j_idx, , drop = FALSE]
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

# weight_inv_exp_dist <- function(Y,
#                                 s, 
#                                 W, zeta) {
#   
#   indx_num=sapply(Y,is.numeric)
#   
#   TT <- nrow(Y)
#   P <- ncol(Y)
#   Y_orig=Y
#   
#   # continuous vars
#   Y=Y[,indx_num]
#   sk=apply(Y,2,function(x)IQR(x)/1.35)
#   
#   # 2. Genera indici delle coppie (i < j)
#   pairs <- combn(TT, 2)
#   i_idx <- pairs[1, ]
#   j_idx <- pairs[2, ]
#   n_pairs <- ncol(pairs)
#   
#   diff=matrix(0,ncol=P,nrow=n_pairs)
#   
#   # 3. Estrai le righe corrispondenti
#   Yi <- Y[i_idx, , drop = FALSE]
#   Yj <- Y[j_idx, , drop = FALSE]
#   
#   diff_cont <- abs(Yi - Yj)
#   #sk=apply(diff,2,function(x)sum(x)/TT^2)
#   
#   diff_cont=sweep(diff_cont, 2, sk, FUN = "/")
#   diff[,indx_num]=as.matrix(diff_cont[, names(indx_num)[indx_num]])
#   #diff=sweep(diff, 2, range_Y, FUN = "/")
#   
#   # Categorical vars
#   if(sum(indx_num)!=P){
#     
#     # FUNZIONA SOLO SE LE VARIABILI CATEGORIALI SONO ALMENO DUE
#     
#     Y=Y_orig[,!indx_num]
#     
#     # 3. Estrai le righe corrispondenti
#     Yi <- Y[i_idx, , drop = FALSE]
#     Yj <- Y[j_idx, , drop = FALSE]
#     diff_cat <- apply(Yi != Yj, 2, as.numeric)
#     #sk=apply(diff,2,function(x)sum(x)/TT^2)
#     
#     diff[,!indx_num]=as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
#     
#     
#   }
#   
#   
#   # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
#   W_si <- W[s[i_idx], , drop = FALSE]
#   W_sj <- W[s[j_idx], , drop = FALSE]
#   max_w <- pmax(W_si, W_sj)
#   
#   # 5. Calcola la distanza finale
#   weighted_exp <- exp(-diff / zeta) * max_w
#   dist_vals <- -zeta * log(rowSums(weighted_exp))
#   
#   # 6. Ricostruzione della matrice simmetrica
#   mat <- matrix(0, TT, TT)
#   mat[cbind(i_idx, j_idx)] <- dist_vals
#   mat[cbind(j_idx, i_idx)] <- dist_vals
#   
#   return(mat)
# }
# 
# 
# weight_inv_exp_dist_medoids <- function(Y, Ymedoids, s, W, zeta) {
#   indx_num <- sapply(Y, is.numeric)
#   
#   TT <- nrow(Y)
#   K <- nrow(Ymedoids)
#   P <- ncol(Y)
#   
#   Y_orig <- Y
#   Ymedoids_orig <- Ymedoids
#   
#   # Continuous variables
#   Y <- Y[, indx_num, drop = FALSE]
#   Ymedoids_cont <- Ymedoids[, indx_num, drop = FALSE]
#   sk <- apply(Y, 2, function(x) IQR(x) / 1.35)
#   
#   # Pre-allocate output matrix
#   mat <- matrix(0, nrow = TT, ncol = K)
#   
#   # Ciclo su tutti i medoids
#   for (k in 1:K) {
#     diff_cont <- abs(sweep(Y, 2, as.numeric(Ymedoids_cont[k,]), FUN = "-"))
#     
#     diff_cont <- sweep(diff_cont, 2, sk, FUN = "/") # normalizzazione
#     diff <- matrix(0, nrow = TT, ncol = P)
#     diff[, indx_num] <- as.matrix(diff_cont[, names(indx_num)[indx_num]])
#     
#     # Categorical variables
#     if (sum(indx_num) != P) {
#       Y_cat <- Y_orig[, !indx_num, drop = FALSE]
#       medoid_cat <- Ymedoids_orig[k, !indx_num, drop = FALSE]
#       diff_cat <- sweep(as.matrix(Y_cat), 2, as.matrix(medoid_cat), FUN = "!=") * 1
#       diff[, !indx_num] <- as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
#     }
#     
#     # Pesatura con pesi W
#     W_si <- W[s, , drop = FALSE]
#     W_sj <- matrix(W[k, ], nrow = TT, ncol = ncol(W), byrow = TRUE)
#     max_w <- pmax(W_si, W_sj)
#     
#     weighted_exp <- exp(-diff / zeta) * max_w
#     dist_vals <- -zeta * log(rowSums(weighted_exp))
#     
#     mat[, k] <- dist_vals
#   }
#   
#   return(mat)  # matrice T x K
# }

WCD_old=function(s,Y,K){
  #TT <- nrow(Y)
  P <- ncol(Y)
  
  wcd=matrix(0,nrow=K,ncol=P)
  
  sk=apply(Y,2,function(x)IQR(x)/1.35)
  
  for(i in 1:K){
    Ys=Y[s==i,]
    TTk <- nrow(Ys)
    pairs <- combn(TTk, 2)
    i_idx <- pairs[1, ]
    j_idx <- pairs[2, ]
    n_pairs <- ncol(pairs)
    
    Yi <- Ys[i_idx, , drop = FALSE]
    Yj <- Ys[j_idx, , drop = FALSE]
    diff <- abs(Yi - Yj)
    #sk=apply(diff,2,function(x)sum(x)/TTk^2)
    
    diff=sweep(diff, 2, sk, FUN = "/")
    # diff=sweep(diff, 2, range_Y, FUN = "/")
    
    for(p in 1:P){
      mat <- matrix(0, TTk, TTk)
      mat[cbind(i_idx, j_idx)] <- diff[,p]
      mat[cbind(j_idx, i_idx)] <- diff[,p]
      wcd[i,p]=mean(apply(mat,1,median))
    }
    # 
    #wcd[i,]=colSums(diff)/TTk^2
    #wcd[i,]=apply(diff,2,median)/TTk
  }
  return(wcd)
  
}

# WCD=function(s,Y,K){
#   #TT <- nrow(Y)
#   P <- ncol(Y)
#   
#   Y_orig=Y
#   
#   wcd=matrix(0,nrow=K,ncol=P)
#   
#   indx_num=sapply(Y,is.numeric)
#   
#   # 1. Normalizzazione Gower
#   # range_Y <- apply(Y, 2, function(col) {
#   #   r <- max(col) - min(col)
#   #   if (r == 0) 1 else r
#   # })
#   sk=apply(Y[,indx_num],2,function(x)IQR(x)/1.35)
#   
#   for(i in 1:K){
#     Y=Y_orig[,indx_num]
#     Ys=Y[s==i,]
#     TTk <- nrow(Ys)
#     pairs <- combn(TTk, 2)
#     i_idx <- pairs[1, ]
#     j_idx <- pairs[2, ]
#     n_pairs <- ncol(pairs)
#     
#     diff=matrix(0,ncol=P,nrow=n_pairs)
#     
#     Yi <- Ys[i_idx, , drop = FALSE]
#     Yj <- Ys[j_idx, , drop = FALSE]
#     
#     diff_cont <- abs(Yi - Yj)
#     #sk=apply(diff,2,function(x)sum(x)/TT^2)
#     
#     diff_cont=sweep(diff_cont, 2, sk, FUN = "/")
#     diff[,indx_num]=as.matrix(diff_cont[, names(indx_num)[indx_num]])
#     
#     
#     if(sum(indx_num)!=P){
#       
#       # FUNZIONA SOLO SE LE VARIABILI CATEGORIALI SONO ALMENO DUE
#       
#       Y=Y_orig[,!indx_num]
#       
#       # 3. Estrai le righe corrispondenti
#       Yi <- Y[i_idx, , drop = FALSE]
#       Yj <- Y[j_idx, , drop = FALSE]
#       diff_cat <- apply(Yi != Yj, 2, as.numeric)
#       #sk=apply(diff,2,function(x)sum(x)/TT^2)
#       
#       diff[,!indx_num]=as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
#       
#       
#     }
#     
#     
#     for(p in 1:P){
#       mat <- matrix(0, TTk, TTk)
#       mat[cbind(i_idx, j_idx)] <- diff[,p]
#       mat[cbind(j_idx, i_idx)] <- diff[,p]
#       wcd[i,p]=mean(apply(mat,1,median))
#     }
#     # 
#     #wcd[i,]=colSums(diff)/TTk^2
#     #wcd[i,]=apply(diff,2,median)/TTk
#   }
#   return(wcd)
#   
# }


# Sim ---------------------------------------------------------------------

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



# COSA=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F){
#   
#   # Simple version of COSA (no outlier detection)
#   
#   P=ncol(Y)
#   TT=nrow(Y)
#   
#   # best_loss <- NULL
#   # best_s <- NULL
#   # best_W = NULL
#   
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   s=initialize_states(Y,K)
#   # Ymedoids=cluster::pam(Y,k=K)
#   # s=Ymedoids$clustering
#   # Ymedoids=Ymedoids$medoids
#   #s=sample(1:K,TT,replace = T)
#   
#   for (outer in 1:n_outer){
#     
#     ## Clustering
#     #for(inner in 1:n_inner){
#     
#     #Compute distances
#     DW=weight_inv_exp_dist(Y,
#                            s,
#                            W,zeta)
#     medoids=cluster::pam(x=DW,k=K,diss=TRUE)
#     #Ymedoids=Y[medoids$medoids,]
#     s=medoids$clustering
#     
#     # Compute weights
#     
#     Spk=WCD(s,Y,K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     w_loss=sum(W*Spk)  
#     eps_W=mean((W-W_old)^2)
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
#     
#     # print(W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
#     }
#     
#   }
#   
#   
#   return(list(W=W,s=s,medoids=medoids,w_loss=w_loss))
# }

library(DescTools)

lof_star=function(x,knn){
  lof_x=DescTools::LOF(x, knn)
  mean_lof=mean(lof_x)
  sd_lof=sd(lof_x)
  lof_st=(lof_x-mean_lof)/sd_lof
  return(lof_st)
}

v_1_R <- function(Y, knn = 10, c = 2, M = NULL, feat_type = NULL,scale="i") {
  # Y           : data frame or matrix (T x P)
  # knn         : number of neighbors
  # c           : cutoff for modified score (default 2)
  # M           : optional fixed M; if NULL, compute per-point from neighbors
  # feat_type   : vector of 0/1/2 for cont/cat/ord (default all continuous)
  
  T <- nrow(Y)
  P <- ncol(Y)
  if (is.null(feat_type)) feat_type <- rep(0L, P)
  
  # 1) Gower dissimilarity
  D <- gower_dist(as.matrix(Y), as.matrix(Y), feat_type = feat_type,scale=scale)
  
  # 2) k-distance and neighborhoods
  d_knn <- numeric(T)
  N_knn <- vector("list", T)
  for(i in seq_len(T)) {
    # distance to all others, sorted
    dists <- sort(D[i, -i])
    d_knn[i] <- dists[knn]
    # include ties, exclude self
    N_knn[[i]] <- which(D[i, ] <= d_knn[i] & seq_len(T) != i)
  }
  
  # 3) Reachability distances
  reach_dist <- vector("list", T)
  for(i in seq_len(T)) {
    o <- N_knn[[i]]
    reach_dist[[i]] <- vapply(o,
                              function(j) max(d_knn[j], D[i, j]),
                              numeric(1)
    )
  }
  
  # 4) Local Reachability Density
  lrd <- sapply(reach_dist, function(rd) {
    if (length(rd) == 0) NA_real_ else 1/mean(rd)
  })
  
  # 5) Standard LOF
  lof <- numeric(T)
  for(i in seq_len(T)) {
    o <- N_knn[[i]]
    if (length(o) == 0 || is.na(lrd[i])) {
      lof[i] <- NA_real_
    } else {
      lof[i] <- mean(lrd[o] / lrd[i])
    }
  }
  
  # 6) Scaled LOF*
  lof_star <- numeric(T)
  for(i in seq_len(T)) {
    o <- N_knn[[i]]
    if (length(o) < 2 || is.na(lof[i])) {
      lof_star[i] <- NA_real_
    } else {
      mu_nb <- mean(lof[o])
      sd_nb <- sd(lof[o])
      if (sd_nb == 0) sd_nb <- 1
      lof_star[i] <- (lof[i] - mu_nb) / sd_nb
    }
  }
  
  # 7) Modified score v_t
  v <- numeric(T)
  for(i in seq_len(T)) {
    ls <- lof_star[i]
    o  <- N_knn[[i]]
    if (is.na(ls) || length(o) == 0) {
      v[i] <- NA_real_; next
    }
    
    # Determine M_i
    if (!is.null(M)) {
      M_i <- M
    } else {
      # use neighbors' lof*
      nb_star <- lof_star[o]
      med_nb  <- median(nb_star, na.rm = TRUE)
      mad_nb  <- median(abs(nb_star - med_nb), na.rm = TRUE)
      M_i     <- med_nb + mad_nb
    }
    
    # piecewise
    if (ls <= M_i) {
      v[i] <- 1
    } else if (ls >= c) {
      v[i] <- 0
    } else {
      t <- (ls - M_i) / (c - M_i)
      v[i] <- (1 - t^2)^2
    }
  }
  
  return(v)
}


v_1_old=function(x,knn=10,c=2,M=NULL){

  lof_st=lof_star(x,knn)

  if(is.null(M)){
    M=median(lof_st)+mad(lof_st)
  }

  v=rep(1,dim(x)[1])
  v[lof_st>=c]=0
  indx=which(M<lof_st&lof_st<c)
  v[indx]=(1-((lof_st[indx]-M)/(c-M))^2)^2
  return(v)
}

# robust_COSA=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F,knn=10,c=2,M=NULL){
#   
#   # Robust version of COSA
# 
#   
#   P=ncol(Y)
#   TT=nrow(Y)
#   
#   # best_loss <- NULL
#   # best_s <- NULL
#   # best_W = NULL
#   
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   s=initialize_states(Y,K)
#   
#   for (outer in 1:n_outer){
#     
#     ## Clustering
#     #for(inner in 1:n_inner){
#     
#     v1=v_1(W[s,]*Y,knn=knn,c=c,M=M)
#     v2=v_1(Y,knn=knn,c=c,M=M)
#     
#     v=apply(cbind(v1,v2),1,min)
#     
#     #Compute distances
#     DW=weight_inv_exp_dist(Y * v,
#                            s,
#                            W,zeta)
#     medoids=cluster::pam(x=DW,k=K,diss=TRUE)
#     #Ymedoids=Y[medoids$medoids,]
#     s=medoids$clustering
#     
#     # Compute weights
#     
#     Spk=WCD(s,Y * v,K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     w_loss=sum(W*Spk)  
#     eps_W=mean((W-W_old)^2)
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
#     
#     # print(W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
#     }
#     
#   }
#   
#   
#   return(list(W=W,s=s,medoids=medoids,w_loss=w_loss,v=v))
# }



JM_COSA=function(Y,zeta0,lambda,K,tol,n_outer=20,alpha=.1,verbose=F,Ts=NULL){
  P=ncol(Y)
  TT=nrow(Y)
  
  # if(is.null(Ts)){
  #   Ts=round(TT/2)
  # }
  
  Gamma <- lambda * (1 - diag(K))
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  # Ottimizzare questo che segue??
  s=initialize_states(Y,K)
  
  loss_old=1e10
  
  for (outer in 1:n_outer){
    
    # subsample=sample(1:TT,Ts,replace = F)
    # Ys=Y[subsample,]
    # ss=s[subsample]
    
    # Questo qui sotto lentissimo se T>1000
      
      DW=weight_inv_exp_dist(Y,
                             s,
                             W,zeta)
      
      # PAM only on subsample
      # DW_1=weight_inv_exp_dist(Ys,
      #                        ss,
      #                        W,zeta)
      
      medoids=cluster::pam(x=DW,k=K,diss=TRUE)
      #Ymedoids=Ys[medoids$medoids,]
      #Ymedoids=Y[medoids$medoids,]
      
      # Questo se lavoro su tutto il dato
      loss_by_state=DW[,medoids$id.med]
      
      # Questo se lavoro su sottocampioni per ottimizzare i tempi
      #loss_by_state=weight_inv_exp_dist_medoids(Y, Ymedoids, s, W, zeta)
      
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
      if (length(unique(s)) < K) {
        s=s_old
        break
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
    #}
    
    # Compute weights
    
    Spk=WCD(s,Y,K)
    wcd=exp(-Spk/zeta0)
    W=wcd/rowSums(wcd)
    
    #}
    
    eps_W=mean((W-W_old)^2)
    if (!is.null(tol)) {
      if (eps_W < tol) {
        break
      }
    }
    
    W_old=W
    zeta=zeta+alpha*zeta0
    
    #print(W)
    # print(zeta)
    # print(Spk)
    # print(zeta0)
    # print(range(DW))
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
      #cat(sprintf('Out iteration %d (# inn iterations %d): %.6e\n', outer, inner, eps_W))
    }
    
  }
  return(list(W=W,s=s,medoids=medoids))
}

simulate_sparse_hmm <- function(seed,
                                TT, P, K,
                                rel,
                                mu   = 2,
                                rho  = 0,
                                nu   = 100,
                                phi  = 0.8,
                                pers = 0.99,
                                sd_noise    = 400,
                                out_pct     = 0.02,
                                out_scale   = 100) {
  #' Simulate data from a Student-t HMM with selective features and outliers
  #'
  #' Arguments:
  #'   seed         Integer; random seed for reproducibility.
  #'   TT           Integer; number of time-series observations.
  #'   P            Integer; total number of features.
  #'   K            Integer; number of hidden states.
  #'   rel          List of length K; for each state k (1…K), rel[[k]] is an integer vector
  #'                of 1-based feature indices that are relevant in state k.
  #'   mu           Numeric; location parameter passed to sim_data_stud_t.
  #'   rho          Numeric; autocorrelation parameter for sim_data_stud_t.
  #'   nu           Numeric; degrees of freedom for the Student-t in sim_data_stud_t.
  #'   phi          Numeric; AR-parameter for sim_data_stud_t.
  #'   pers         Numeric; persistence parameter for sim_data_stud_t.
  #'   sd_noise     Numeric; s.d. of a Normal distribution used to generate noise
  #'                on irrelevant features.
  #'   out_pct      Numeric in [0,1]; fraction of observations to turn into outliers.
  #'   out_scale    Numeric; standard deviation of the additive normal outlier noise.
  #'
  #' Value:
  #'   A list with components:
  #'     • Y        — Numeric matrix of dimension TT × P containing the simulated data.
  #'     • truth    — Integer vector of length TT giving the true latent state (1…K),
  #'                  with 0 marking observations designated as outliers.
  #'     • rel_list — The same list rel provided as input, indicating per-state relevant features.
  #'
  
 
  set.seed(seed)
  # 1) genera la catena di stati e i dati base student-t
  simDat <- sim_data_stud_t(seed  = seed,
                            TT    = TT,
                            P     = P,
                            Pcat  = NULL,
                            Ktrue = K,
                            mu    = mu,
                            rho   = rho,
                            nu    = nu,
                            phi   = phi,
                            pers  = pers)
  Y0    <- as.matrix(simDat$SimData)
  truth <- simDat$mchain
  
  # 2) costruisci Y: per ogni cella (i,j), se j in rel[[ truth[i] ]] 
  #    tieni Y0[i,j], altrimenti campiona rumore t_student( df=noise_df )*noise_scale
  Y <- matrix(NA_real_, nrow = TT, ncol = P)
  for (i in seq_len(TT)) {
    k <- truth[i]
    # feature rilevanti per stato k
    keep <- rel[[k]]
    # tutte le altre
    drop <- setdiff(seq_len(P), keep)
    # preserva le rilevanti…
    Y[i, keep] <- Y0[i, keep]
    # …e genera rumore indipendente sulle irrilevanti
    Y[i, drop] <- rnorm(length(drop),mean=rnorm(1), sd=sd_noise) 
  }
  
  # 3) genera outlier: sovrascrive outlier_pct*TT righe con un rumore gaussiano di sd=out_scale
  if(out_pct>0){
    N_out <- ceiling(out_pct * TT)
    out_idx <- sample.int(TT, N_out)
    Y[out_idx, ] <- Y[out_idx, ] + matrix(rnorm(N_out * P, sd = out_scale),
                                          nrow = N_out, ncol = P)
    # marca gli outlier con stato 0
    truth[out_idx] <- 0
  }
  
  
  W_truth <- matrix(FALSE, nrow = K, ncol = P)
  
  for (k in seq_len(K)) {
    W_truth[k, rel[[k]]] <- TRUE
  }
  
  # restituisci
  list(
    Y         = Y,
    truth     = truth,
    rel_list  = rel,
    W_truth   = W_truth
  )
}

get_feat_type <- function(Y) {
  # Y: data.frame or matrix
  types <- apply(Y,2, function(col) {
    if (is.numeric(col)) {
      # numeric or integer
      0L
    } else if (is.ordered(col)) {
      # ordered factor
      2L
    } else if (is.factor(col) || is.character(col)) {
      # unordered factor or character
      1L
    } else {
      stop(sprintf("Unsupported column type: %s", class(col)[1]))
    }
  })
  unname(as.integer(types))
}

robust_sparse_jump <- function(Y,
                           zeta0,
                           lambda,
                           K,
                           tol     = NULL,
                           n_init  = 10,
                           n_outer = 20,
                           alpha   = 0.1,
                           verbose = FALSE,
                           knn     = 10,
                           c       = 5,
                           M       = NULL,
                           scale="i",
                           parallel=F,
                           n_cores=NULL) {
  library(Rcpp)
 
  Rcpp::sourceCpp("robJM_R.cpp")
  
  P  <- ncol(Y)
  TT <- nrow(Y)
  
  feat_type <- get_feat_type(Y)
  
  # Transform Y into a numeric matrix
  if(!is.matrix(Y)){
    Y=as.matrix(
      data.frame(
        lapply(Y, function(col) {
          if (is.factor(col)||is.character(col)) as.integer(as.character(col))
          else              as.numeric(col)
        })
      )
    )
  }
  
  Gamma <- lambda * (1 - diag(K))
  
  run_one <- function(init_id) {
    # 1) initial W, zeta, s
    W        <- matrix(1/P, nrow=K, ncol=P)
    W_old    <- W
    zeta     <- zeta0
    loss_old <- Inf
    
    s=initialize_states(Y, K,feat_type=feat_type,scale=scale)
    
    # Only this because I cannot weight categorical or ordinal variables
    v2=v_1(Y,                     knn=knn, c=c, M=M,scale=scale)
    # v=rep(1,TT)
    # VV=v %*% t(v)
    # VV=sqrt(VV)
    
    for (outer in seq_len(n_outer)) {
      # 2) local scales
      v1 <- v_1(W[s, , drop=FALSE] * Y, knn=knn, c=c, M=M,scale=scale)
      # v2 <- v_1(Y,                     knn=knn, c=c, M=M)
      v  <- pmin(v1, v2)
      
      
      # 3) weighted distances + PAM
      
      DW      <- weight_inv_exp_dist(v*Y, s, W, zeta
                                     ,
                                     feat_type=feat_type,scale=scale
                                     )
      
      pam_out <- cluster::pam(DW, k=K, diss=TRUE)
      medoids <- pam_out$id.med
      
      # 4) build loss-by-state
      loss_by_state <- DW[, medoids, drop=FALSE]  # TT x K
      
      # 5) DP forward: V[t,j] = loss[t,j] + min_i( V[t+1,i] + Gamma[i,j] )
      V <- loss_by_state
      for (t in (TT-1):1) {
        for (j in seq_len(K)) {
          # look at row t+1 of V plus column j of Gamma:
          V[t, j] <- loss_by_state[t, j] +
            min( V[t+1, ] + Gamma[, j] )
        }
      }
      
      # 6) backtrack to get s
      s_old <- s
      # first time‐point
      s[1] <- which.min(V[1, ])
      loss  <- V[1, s[1]]
      # subsequent
      for (t in 2:TT) {
        prev <- s[t-1]
        # pick state j minimizing V[t,j] + penalty from prev
        scores <- V[t, ] + Gamma[prev, ]
        s[t] <- which.min(scores)
      }
      
      # 7) must have all K states or revert
      if (length(unique(s)) < K) {
        s <- s_old
        break
      }
      
      # 8) loss‐convergence
      if (!is.null(tol) && (loss_old - loss) < tol) break
      loss_old <- loss
      
      # 9) update W via WCD + exp
      Spk <- WCD(s, v*Y , K
                 #, feat_type=feat_type,v=rep(1,TT),scale=scale
                 )
      
      wcd <- exp(-Spk / zeta0)
      W   <- wcd / rowSums(wcd)
      
      # 10) W‐convergence
      epsW <- mean((W - W_old)^2)
      if (!is.null(tol) && epsW < tol) break
      W_old <- W
      
      # 11) bump zeta
      zeta <- zeta + alpha * zeta0
      
      if (verbose) {
        cat(sprintf("init %2d, outer %2d → loss=%.4e, epsW=%.4e, zeta=%.3f\n",
                    init_id, outer, loss, epsW, zeta))
      }
    }
    
    list(W      = W,
         s      = s,
         medoids= medoids,
         v      = v,
         loss   = loss)
  }
  
  # run n_init times, pick the one with smallest loss
  if(parallel){
    if(is.null(n_cores)){
      n_cores=parallel::detectCores()-1
    }
    res_list=mclapply(
      seq_len(n_init),
      run_one,
      mc.cores = n_cores
    )
  }
  else{
    res_list <- lapply(seq_len(n_init), run_one)
    
  }
  losses   <- vapply(res_list, `[[`, numeric(1), "loss")
  best_run <- res_list[[ which.min(losses) ]]
  
  best_s   <- best_run$s
  best_loss<- best_run$loss
  best_W = best_run$W
  best_medoids  <- Y[best_run$medoids,]
  best_v <- best_run$v
  
  mif=which.max(apply(best_W,2,sum))
  
  # Re‐order states based on most important feature state-conditional median
  new_best_s <- order_states_condMed(Y[, mif], best_s)
  
  tab <- table(best_s, new_best_s)
  new_order <- apply(tab, 1, which.max)
  
  best_W <- best_W[new_order,]
  
  ret_list=list(W = best_W,
                s = new_best_s,
                medoids = best_medoids,
                v = best_v,
                loss = best_loss,
                zeta0 = zeta0,
                lambda = lambda,
                c = c,
                knn=knn,
                M = M)
  
  return(ret_list)
  
}

cv_robust_sparse_jump <- function(
    Y,
    true_states,
    K_grid=NULL,
    zeta0=NULL,
    lambda_grid=NULL,
    n_folds = 5,
    parallel=F,
    n_cores=NULL,
    cv_method="blocked-cv",
    knn=10,
    c=5,
    M=NULL
) {
  
  # cv_sparse_jump: Cross-validate Sparse Jump Model parameters (K and lambda)
  
  # Arguments:
  #   Y           - data matrix (N × P)
  #   true_states - vector of true states for ARI computation
  #   K_grid      - vector of candidate numbers of states
  #   zeta0       - sparsity hyperparameter (if NULL it is set to 0.2)
  #   lambda_grid - vector of candidate lambdas
  #   n_folds     - number of folds for cross-validation (default: 5)
  #   parallel    - logical; TRUE for parallel execution (default: FALSE)
  #   n_cores     - number of cores to use for parallel execution (default:  NULL)
  #   cv_method   - method for cross-validation: "blocked-cv" or "forward-chain"
  #   knn         - number of nearest neighbors for LOF (default: 10)
  #   c           - lower threshold for LOF (default: 5)
  #   M           - upper threshold for LOF (default: NULL, uses median + mad)
  
  
  # Value:
  #   A data.frame with one row per (K, zeta0, lambda) combination, containing:
  #     K      - number of states tested
  #     zeta0  - sparsity hyperparameter value
  #     lambda - jump penalty value
  #     ARI    - mean Adjusted Rand Index across folds
  
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  
  if(is.null(lambda_grid)) {
    lambda_grid <- seq(0,1,.1)  # Default range for lambda
  }
  
  if(is.null(zeta0)) {
    zeta0 <- 0.2  # Default sparsity hyperparameter
  }
  
  # Libreria per ARI
  library(mclust)
  
  N <- nrow(Y)
  P <- ncol(Y)
  
  # Suddivido gli N campioni in n_folds blocchi contigui
  if(cv_method=="blocked-cv"){
    fold_indices <- split(
      1:N,
      rep(1:n_folds, each = ceiling(N / n_folds), length.out = N)
    )
    
  }
  else if(cv_method=="forward-chain"){
    fold_indices <- lapply(seq_len(n_folds), function(k) {
      idx_end <- N - (k - 1)
      1:(idx_end-1)
    })
    names(fold_indices) <- as.character(seq_len(n_folds))
  }
  else{
    stop("cv_method must be either 'blocked-cv' or 'forward-chain'")
  }
  
  # Funzione che, per una tripla (K, kappa, lambda) e un fold (train_idx, val_idx),
  # calcola l’ARI sui punti di validazione
  fold_ari <- function(K, kappa, lambda, train_idx, val_idx) {
    # 1) Fit del modello sparse_jump su soli dati di TRAIN
    res <- robust_sparse_jump(Y=Y[train_idx, , drop = FALSE],
                          zeta0=zeta0,
                          lambda=lambda,
                          K=K,
                          tol        = 1e-16,
                          n_init     = 5,
                          n_outer    = 10,
                          alpha      = 0.1,
                          verbose    = F,
                          knn        = knn,
                          c          = c,
                          M          = M)
    states_train <- res$s
    feat_idx     <- which(colSums(res$W) > 0.025)
    
    # Se non vengono selezionate feature, restituisco ARI = 0
    if (length(feat_idx) == 0) {
      return(0)
    }
    
    # 2) Extract medoids
    medoids=res$medoids[,feat_idx]
    
    # 3) Assegno ciascun punto in VAL a uno stato: 
    #    lo stato k che minimizza la distanza gower su feature selezionate
    Y_val_feats     <- as.matrix(Y[val_idx, feat_idx, drop = FALSE])
    pred_val_states <- integer(length(val_idx))
    
    dists<- gower_dist(Y_val_feats,medoids)
    
    pred_val_states<- apply(dists,1,which.min)
    
    
    # 4) Calcolo ARI tra etichette vere e quelle predette sul blocco VAL
    return(
      adjustedRandIndex(true_states[val_idx], pred_val_states)
    )
  }
  
  # Costruisco la griglia di tutte le combinazioni di (K, kappa, lambda)
  grid <- expand.grid(
    K      = K_grid,
    # For kappa we take a representative value, we will select kappa later based on GAP stat
    zeta0  = zeta0,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Data.frame in cui raccogliere (K, kappa, lambda, media_ARI)
  results <- data.frame(
    K      = integer(0),
    zeta0  = integer(0),
    lambda = numeric(0),
    ARI    = numeric(0)
  )
  
  # Loop su ciascuna riga della griglia
  if (!parallel) {
    # applichiamo una funzione su ogni riga di 'grid'
    results_list <- lapply(seq_len(nrow(grid)), function(row) {
      K_val     <- grid$K[row]
      zeta0_val <- grid$zeta0[row]
      lambda_val<- grid$lambda[row]
      
      # calcolo ARI su ciascun fold
      ari_vals <- numeric(n_folds)
      for (f in seq_len(n_folds)) {
        val_idx   <- fold_indices[[f]]
        train_idx <- setdiff(seq_len(N), val_idx)
        ari_vals[f] <- fold_ari(K_val, zeta0_val, lambda_val, train_idx, val_idx)
      }
      mean_ari <- mean(ari_vals)
      
      # restituisco un data.frame di una sola riga
      data.frame(
        K      = K_val,
        zeta0  = zeta0,
        lambda = lambda_val,
        ARI    = mean_ari,
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
        zeta0_val <- as.integer(grid$zeta0[row])
        lambda_val<-        grid$lambda[row]
        
        # calcolo ARI su ciascun fold
        ari_vals <- numeric(n_folds)
        for (f in seq_len(n_folds)) {
          val_idx   <- fold_indices[[f]]
          train_idx <- setdiff(seq_len(N), val_idx)
          ari_vals[f] <- fold_ari(K_val, zeta0_val, lambda_val, train_idx, val_idx)
        }
        mean_ari <- mean(ari_vals)
        
        # ritorno un data.frame di una sola riga
        data.frame(
          K      = K_val,
          zeta0  = zeta0,
          lambda = lambda_val,
          ARI    = mean_ari,
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

permute_gap <- function(Y) {
  if (!is.matrix(Y)) {
    stop("Input Y must be a matrix.")
  }
  TT <- nrow(Y)
  P <- ncol(Y)
  
  Y_perm <- apply(Y, 2, function(col) {
    sample(col, size = TT, replace = FALSE)
  })
  
  if (P == 1) {
    Y_perm <- matrix(Y_perm, nrow = T, ncol = 1)
  }
  
  rownames(Y_perm) <- rownames(Y)
  colnames(Y_perm) <- colnames(Y)
  
  return(Y_perm)
}

gap_robust_sparse_jump=function(
    Y,
    K_grid=NULL,
    zeta0_grid=NULL,
    lambda=0,
    B=10,
    parallel=F,
    n_cores=NULL,
    knn=10,
    c=5,
    M=NULL
){
  
  # gap_robust_sparse_jump: Compute the Gap Statistic for Robust Sparse Jump Model
  
  # Arguments:
  #   Y           - data matrix (N × P)
  #   K_grid      - vector of candidate numbers of states (default: seq(2, 4, by = 1))
  #   zeta0_grid  - vector of candidate kappa values (default: seq(0.05, 0.4, by = 0.05))
  #   lambda      - jump penalty value (default: 0.2)
  #   B           - number of bootstrap samples (default: 10)
  #   parallel    - logical; TRUE for parallel execution (default: FALSE)
  #   n_cores     - number of cores to use for parallel execution (default: NULL)
  #   knn         - number of nearest neighbors for LOF (default: 10)
  #   c           - lower threshold for LOF (default: 5)
  #   M           - upper threshold for LOF (default: NULL, uses median + mad)
  
  # Value:
  #   A list containing:
  #     gap_stats - data.frame with Gap Statistic results
  #     plot_res  - ggplot object of Gap Statistic vs zeta0
  #     meta_df   - data.frame with detailed results for each (K, zeta0, lambda, b)
  
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  
  if(is.null(zeta0_grid)) {
    zeta0_grid <- seq(0.05,.4,.05)  # Default range for lambda
  }
  
  if(is.null(lambda)){
    lambda=0.2
  }
  
  # Libreria per ARI
  library(mclust)
  
  N <- nrow(Y)
  P <- ncol(Y)
  
  grid <- expand.grid(zeta0 = zeta0_grid, 
                      lambda = lambda, 
                      K = K_grid, 
                      b = 0:B)
  
  gap_one_run=function(zeta0,lambda,K,b){
    if (b == 0) {
      Y_input <- Y
      permuted <- FALSE
    } else {
      # Permute features for kappa
      Y_input <- permute_gap(Y)
      permuted <- TRUE
    }
    # Fit the model
    
    fit=robust_sparse_jump(Y,
                       zeta0=zeta0,
                       lambda=lambda,
                       K=K,
                       tol     = NULL,
                       n_init  = 5,
                       n_outer = 20,
                       alpha   = 0.1,
                       verbose = FALSE,
                       knn     = knn,
                       c       = c,
                       M       = M)
    
    return(list(loss=fit$loss,
                K=K,
                permuted=permuted,
                zeta0=zeta0,
                lambda=lambda
    ))
    
  }
  
  if(parallel){
    if(is.null(n_cores)){
      n_cores <- parallel::detectCores() - 1
    }
    results <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
      params <- grid[i, ]
      gap_one_run(
        zeta0 = params$zeta0,
        lambda= params$lambda,
        K     = params$K,
        b     = params$b
      )
    }, mc.cores = mc_cores)
  }
  else{
    results <- lapply(seq_len(nrow(grid)), function(i) {
      params <- grid[i, ]
      gap_one_run(
        zeta0 = params$zeta0,
        lambda= params$lambda,
        K     = params$K,
        b     = params$b
      )
    })
  }
  
  meta_df <- do.call(rbind.data.frame, c(results, make.row.names = FALSE))
  
  library(dplyr)
  gap_stats <- meta_df %>%
    group_by(K, zeta0) %>%
    summarise(
      log_O = log(loss[!permuted]),
      log_O_star_mean = mean(log(loss[permuted])),
      se_log_O_star=sd(log(loss[permuted])),
      GAP =  log_O - log_O_star_mean,
      .groups = 'drop'
    )
  
  library(ggplot2)
  
  plot_res=ggplot(gap_stats, aes(x = zeta0, y = GAP, color = factor(K))) +
    geom_line() +
    geom_point() +
    scale_color_discrete(name = "Number of clusters\n(K)") +
    labs(
      x = expression(kappa),
      y = "Gap Statistic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  return(list(gap_stats=gap_stats,
              plot_res=plot_res,
              meta_df=meta_df))
  
  
}


# robust_JM_COSA=function(Y,zeta0,lambda,K,tol,n_outer=20,alpha=.1,
#                         verbose=F,knn=10,c=2,M=NULL){
#   library(Rcpp)
#   Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
#   Rcpp::sourceCpp("wcd.cpp")
#   P=ncol(Y)
#   TT=nrow(Y)
# 
#   Gamma <- lambda * (1 - diag(K))
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   # Multiple initialization, keep best one (lower loss)
#   s=initialize_states(Y,K)
#   
#   loss_old=1e10
#   
#   for (outer in 1:n_outer){
#     
#     # subsample=sample(1:TT,Ts,replace = F)
#     # Ys=Y[subsample,]
#     # ss=s[subsample]
#     
#     v1=v_1(W[s,]*Y,knn=knn,c=c,M=M)
#     v2=v_1(Y,knn=knn,c=c,M=M)
#     
#     v=apply(cbind(v1,v2),1,min)
#     
#     #Compute distances
#     DW=weight_inv_exp_dist(as.matrix(Y * v),
#                            s,
#                            W,zeta)
#     medoids=cluster::pam(x=DW,k=K,diss=TRUE)
#     
#     # Questo se lavoro su tutto il dato
#     loss_by_state=DW[,medoids$id.med]
#     
#     # Questo se lavoro su sottocampioni per ottimizzare i tempi
#     #loss_by_state=weight_inv_exp_dist_medoids(Y, Ymedoids, s, W, zeta)
#     
#     V <- loss_by_state
#     for (t in (TT-1):1) {
#       V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
#     }
#     s_old=s
#     s[1] <- which.min(V[1,])
#     for (t in 2:TT) {
#       s[t] <- which.min(V[t,] + Gamma[s[t-1],])
#     }
#     loss <- min(V[1,])
#     if (length(unique(s)) < K) {
#       s=s_old
#       break
#     }
#     
#     epsilon <- loss_old - loss
#     if (!is.null(tol)) {
#       if (epsilon < tol) {
#         break
#       }
#     } 
#     # else if (all(s == s_old)) {
#     #   break
#     # }
#     loss_old <- loss
#     #}
#     
#     # Compute weights
#     
#     Spk=WCD(s,as.matrix(Y * v),K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     #}
#     
#     eps_W=mean((W-W_old)^2)
#     
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
# 
#     # print(W)
#     # print(epsilon)
#     # print(eps_W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Iteration %d: %.6e\n', outer, loss))
#       #cat(sprintf('Out iteration %d (# inn iterations %d): %.6e\n', outer, inner, eps_W))
#     }
#     
#   }
#   return(list(W=W,s=s,medoids=medoids,v=v,loss=loss))
# }

# RJM_COSA_gap=function(Y,
#                       zeta_grid=seq(0.1,.7,.1),
#                       lambda_grid=seq(0,1,.1),
#                       K_grid=2:6,
#                       tol=NULL,n_outer=20,alpha=.1,verbose=F,n_cores=NULL,
#                       B=10, knn=10,c=2,M=NULL){
#   
#   # B is the number of permutations
#   
#   grid <- expand.grid(zeta0 = zeta_grid, 
#                       lambda = lambda_grid, 
#                       K = K_grid, b = 0:B)
#   
#   library(foreach)
#   library(doParallel)
#   
#   if(is.null(n_cores)){
#     n_cores <- parallel::detectCores() - 1
#   } 
#   
#   # Set up cluster
#   cl <- makeCluster(n_cores)
#   registerDoParallel(cl)
#   
#   results_list <- foreach(i = 1:nrow(grid), .combine = 'list',
#                           .packages = c("cluster","Rcpp","DescTools"),
#                           .multicombine = TRUE,
#                           .export = c("Y", "robust_JM_COSA", 
#                                       #"WCD", "weight_inv_exp_dist",
#                                       "initialize_states",
#                                       "v_1","lof_star",
#                                       "grid", "tol", "n_outer", "alpha",
#                                       "knn","c","M")) %dopar% {
#                                         K_val <- grid$K[i]
#                                         zeta_val <- grid$zeta0[i]
#                                         lambda_val=grid$lambda[i]
#                                         b <- grid$b[i]
#                                         
#                                         set.seed(b + 1000 * i)
#                                         
#                                         if (b == 0) {
#                                           Y_input <- Y
#                                           permuted <- FALSE
#                                         } else {
#                                           # Permute features for zeta0
#                                           Y_input <- apply(Y, 2, sample)
#                                           # Permute rows for lambda
#                                           Y_input <- Y_input[ sample(nrow(Y_input),
#                                                                      size = nrow(Y_input),
#                                                                      replace = FALSE), ]
#                                           permuted <- TRUE
#                                         }
#                                         
#                                         res <- robust_JM_COSA(Y_input, zeta0 = zeta_val, 
#                                                               lambda = lambda_val, K = K_val, tol = tol,
#                                                               n_outer = n_outer, alpha = alpha, verbose = FALSE,
#                                                               knn=knn,c=c,M=M)
#                                         
#                                         list(
#                                           meta = data.frame(zeta0 = zeta_val, lambda=lambda_val,K = K_val, 
#                                                             loss = res$loss, permuted = permuted),
#                                           cosa = if (!permuted) list(zeta0 = zeta_val, lambda=lambda_val,K = K_val, 
#                                                                      W = res$W, s = res$s, 
#                                                                      medoids = res$medoids$medoids,
#                                                                      v=res$v) else NULL
#                                         )
#                                       }
#   
#   
#   stopCluster(cl)
#   
#   # Flatten results
#   meta_df <- do.call(rbind, lapply(results_list, `[[`, "meta"))
#   cosa_results <- Filter(Negate(is.null), lapply(results_list, `[[`, "cosa"))
#   
#   
#   # Compute GAP
#   library(dplyr)
#   gap_stats <- meta_df %>%
#     group_by(K, zeta0) %>%
#     summarise(
#       log_O = log(loss[!permuted]),
#       log_O_star_mean = mean(log(loss[permuted])),
#       se_log_O_star=sd(log(loss[permuted])),
#       GAP = log_O_star_mean - log_O,
#       .groups = 'drop'
#     )
#   
#   return(list(
#     gap_stats = gap_stats,
#     cosa_results = cosa_results
#   ))
#   
# }

plot_W=function(W){
  library(reshape)
  df <- as.data.frame(W)
  df$Cluster <- factor(paste0("Cluster_", 1:nrow(df)))
  
  # Riorganizziamo in formato lungo
  df_long <- melt(df, id.vars = "Cluster", variable.name = "Feature", value.name = "Weight")
  
  # Converti Feature in fattore per ordinare le colonne
  #df_long$Feature <- factor(df_long$Feature, levels = paste0("V", 1:ncol(fit$W)))
  
  # Bar plot
  library(ggplot2)
  p=ggplot2::ggplot(df_long, aes(x = variable, y = value, fill = Cluster)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Feature Weights by Cluster",
         x = "Feature",
         y = "Weight") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}


