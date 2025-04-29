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

weight_inv_exp_dist <- function(Y,
                                #Ymedoids,
                                #index_medoids, 
                                s, 
                                W, zeta) {
  TT <- nrow(Y)
  P <- ncol(Y)
  
  # 1. Normalizzazione Gower
  # range_Y <- apply(Y, 2, function(col) {
  #   r <- max(col) - min(col)
  #   if (r == 0) 1 else r
  # })
  
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

# WCD=function(s,Y,K){
#   Tk=table(s)
#   wcd=matrix(0,nrow=K,ncol=ncol(Y))
#   for(k in unique(s)){
#     for(p in 1:P){
#       temp=Y[s==k,p]
#       dist_temp=abs(outer(temp, temp, "-"))
#       dist_temp[upper.tri(dist_temp)] <- 0
#       dist_temp=dist_temp/diff(range(temp))
#       wcd[k,p]=sum(dist_temp)/Tk[k]^2
#     }
#   }
#   return(wcd)
# }

WCD=function(s,Y,K){
  #TT <- nrow(Y)
  P <- ncol(Y)
  
  wcd=matrix(0,nrow=K,ncol=P)
  
  # 1. Normalizzazione Gower
  # range_Y <- apply(Y, 2, function(col) {
  #   r <- max(col) - min(col)
  #   if (r == 0) 1 else r
  # })
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

# apply_noise_by_cluster <- function(Y, s, feat_list) {
#   Y_noised <- Y
#   K <- length(feat_list)
#   
#   for (k in 1:K) {
#     cluster_rows <- which(s == k)
#     if (length(cluster_rows) <= 1) next  # nulla da mescolare se solo una riga
#     all_features <- seq_len(ncol(Y))
#     irrelevant_feats <- setdiff(all_features, feat_list[[k]])
#     
#     for (j in irrelevant_feats) {
#       # mescola le osservazioni della colonna j solo tra le righe del cluster k
#       Y_noised[cluster_rows, j] <- sample(Y[, j],size=length(Y_noised[cluster_rows, j]))
#     }
#   }
#   
#   return(Y_noised)
# }

zeta0=0.2
alpha=.1
K=2
tol=1e-9
n_outer=15
verbose=T

TT=500
P=10

simDat=sim_data_stud_t(seed=123,
                       TT=TT,
                       P=P,
                       Pcat=NULL,
                       Ktrue=2,
                       mu=1,
                       rho=0,
                       nu=100,
                       phi=.8,
                       pers=0.95)

Y=simDat$SimData
true_stat=simDat$mchain

plot(Y[,2],col=true_stat,pch=19)

# feat_list=list()
# feat_list[[1]]=1:10
# feat_list[[2]]=5:15
# feat_list[[3]]=10:20
# Y_noised=apply_noise_by_cluster(Y,simDat$mchain,feat_list)
# Y=Y_noised


nu=4
# State 1, only features 1,2 and 3 are relevant
indx=which(true_stat!=1)
Sigma <- matrix(0,ncol=P-3,nrow=P-3)
diag(Sigma)=5
Y[indx,-(1:3)]=mvtnorm::rmvt(length(indx), 
                             sigma = (nu-2)*Sigma/nu, 
                             df = nu, delta = rep(0,P-3))

# State 2, only features 3,4 and 5 are relevant
indx=which(true_stat!=2)
Y[indx,-(3:5)]=mvtnorm::rmvt(length(indx), 
                             sigma = (nu-2)*Sigma/nu, 
                             df = nu, delta = rep(0,P-3))

# # State 3, only features 5,6 and 7 are relevant
# indx=which(true_stat!=3)
# Y[indx,-(5:7)]=mvtnorm::rmvt(length(indx), 
#                              sigma = (nu-2)*Sigma/nu, 
#                              df = nu, delta = rep(0,P-3))

# Other features are not relevant
# Sigma <- matrix(0,ncol=5,nrow=5)
# diag(Sigma)=5
# Y[,6:10]=mvtnorm::rmvt(TT, 
#                        sigma = (nu-2)*Sigma/nu, 
#                        df = nu, delta = rep(0,5))

Sigma <- matrix(0,ncol=P-5,nrow=P-5)
diag(Sigma)=5
Y[,6:P]=mvtnorm::rmvt(TT, 
                      sigma = (nu-2)*Sigma/nu, 
                      df = nu, delta = rep(0,P-5))

x11()
par(mfrow=c(4,3))
for (i in 1:P) {
  plot(Y[, i], col=simDat$mchain, pch=19,ylab=i)
}


COSA=function(Y,zeta0,K,tol,n_outer=20,alpha=.1,verbose=F){
  P=ncol(Y)
  TT=nrow(Y)
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  #s=initialize_states(Y,K)
  # Ymedoids=cluster::pam(Y,k=K)
  # s=Ymedoids$clustering
  # Ymedoids=Ymedoids$medoids
  s=sample(1:K,TT,replace = T)
  
  for (outer in 1:n_outer){
    
    ## Clustering
    #for(inner in 1:n_inner){
    
    #Compute distances
    DW=weight_inv_exp_dist(Y,
                           s,
                           W,zeta)
    medoids=cluster::pam(x=DW,k=K,diss=TRUE)
    Ymedoids=Y[medoids$medoids,]
    s=medoids$clustering
    
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
    
    # print(W)
    # print(zeta)
    print(Spk)
    print(zeta0)
    print(range(DW))
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
    }
    
  }
  return(list(W=W,s=s,medoids=medoids))
}

JM_COSA=function(Y,zeta0,lambda,K,tol,n_outer=20,n_inner=10,alpha=.1,verbose=F){
  P=ncol(Y)
  TT=nrow(Y)
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  Gamma <- lambda * (1 - diag(K))
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  #s=initialize_states(Y,K)
  # Ymedoids=cluster::pam(Y,k=K)
  # s=Ymedoids$clustering
  # Ymedoids=Ymedoids$medoids
  s=sample(1:K,TT,replace = T)
  
  for (outer in 1:n_outer){
    
    ## Clustering
    #for(inner in 1:n_inner){
    loss_old=1e10
    #Compute distances
    for(inner in 1:n_inner){
      DW=weight_inv_exp_dist(Y,
                             s,
                             W,zeta)
      medoids=cluster::pam(x=DW,k=K,diss=TRUE)
      Ymedoids=Y[medoids$medoids,]
      
      loss_by_state=DW[,medoids$id.med]
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
    }
    
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
    
    print(W)
    # print(zeta)
    # print(Spk)
    # print(zeta0)
    # print(range(DW))
    
    if (verbose) {
      #cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
      cat(sprintf('Out iteration %d (# inn iterations %d): %.6e\n', outer, inner, eps_W))
    }
    
  }
  return(list(W=W,s=s,medoids=medoids))
}

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


