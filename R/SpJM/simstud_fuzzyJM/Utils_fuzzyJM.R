
library(cluster)
library(StatMatch)
library(poliscidata) # for weighted mode 
library(pdfCluster) # for ARI

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


# Temporal ----------------------------------------------------------------


sim_data_mixed=function(seed=123,
                        TT,
                        P,
                        Ktrue=3,
                        mu_val=1,
                        phi=.8,
                        rho=0,
                        Pcat=NULL,
                        pers=.95,
                        pNAs=0,
                        typeNA=3){
  
  # Function to simulate mixed data with fixed parameters for the data generating process
  
  # Arguments:
  # seed: seed for the random number generator
  # TT: number of observations
  # P: number of features
  # Ktrue: number of states
  # mu: mean value for the continuous variables
  # phi: conditional probability for the categorical outcome k in state k
  # rho: correlation for the variables
  # Pcat: number of categorical variables
  # pers: self-transition probability
  # pNAs: percentage of missing values
  # typeNA is the type of missing values (0 for random, 1 for continuous, all other values will turn into no missing imputation)
  
  # value:
  # SimData: matrix of simulated data
  
  MU=mu_val
  mu_val=c(-mu_val,0,mu_val)
  
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
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
    u = MASS::mvrnorm(TT,rep(mu_val[k],P),Sigma)
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  if(Pcat!=0){
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=x,mu_val=MU,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  if(typeNA==0|typeNA==1){
    SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=typeNA)
    SimData.NA=as.data.frame(SimData.NA)
    SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
    SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  }
  else{
    SimData.NA=SimData
  }
  
  return(list(
    SimData.NA=SimData.NA,
    SimData.complete=SimData,
    mchain=x,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
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

fuzzy_jump <- function(Y, 
                       n_states, jump_penalty=1e-5, 
                       initial_states=NULL,
                       max_iter=10, n_init=10, tol=NULL, 
                       verbose=FALSE
                       # ,
                       #        time_vec=NULL
                       
) {
  # Fit jump model for mixed type data 
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # n_states: number of states
  # jump_penalty: penalty for the number of jumps
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
  
  # timeflag=FALSE
  # if(!is.null(time_vec)){
  #   timeflag=TRUE
  #   if(length(time_vec)!=nrow(Y)){
  #     stop("time_vec must have the same length of the number of observations")
  #   }
  #   else{
  #     time=sort(unique(time_vec))
  #     dtime=diff(time)
  #     dtime=dtime/as.numeric(min(dtime))
  #     dtime=as.numeric(dtime)
  #   }
  # }
  
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_S <- NULL
  
  # Which vars are categorical and which are numeric
  cat_flag=any(sapply(Y, is.factor))
  
  if(cat_flag){
    cat.indx=which(sapply(Y, is.factor))
    cont.indx=which(sapply(Y, is.numeric))
    Ycont=Y[,cont.indx]
    Ycont=apply(Ycont,2,scale)
    Y[,cont.indx]=Ycont
    Ycat=Y[,cat.indx]
    
    if(length(cat.indx)==1){
      n_levs=length(levels(Ycat))
    }
    else{
      n_levs=apply(Ycat, 2, function(x)length(unique(x[!is.na(x)])))
    }

    n_cat=length(cat.indx)
    n_cont=n_features-n_cat
    
  }
  else{
    cont.indx=1:n_features
    Ycont=Y
    Ycont=apply(Ycont,2,scale)
    Y[,cont.indx]=Ycont
    n_cont=dim(Y)[2]
    n_cat=0
  }
  
  
  
  
  for (init in 1:n_init) {
    
    # State initialization through kmeans++
    if (!is.null(initial_states)) {
      s <- initial_states
    } else {
      s=initialize_states(Y,n_states)
    }
    
    S <- matrix(0, nrow = n_obs, ncol = n_states)
    row_indices <- rep(1:n_obs)  # Row positions in SS
    # Assign 1s in a single step
    S[cbind(row_indices, s)] <- 1 
    
    mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
    if(cat_flag){
      mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    }
    
    for (i in unique(s)) {
      # substitute with medians
      mu[i,] <- apply(Ycont[s==i,], 2, median, na.rm = TRUE)
      if(cat_flag){
        if(length(cat.indx)==1){
            mo[i,]=Mode(Ycat[s==i])
        }
        else{
          mo[i,]=apply(Ycat[s==i,],2,Mode)
        }
      }
    }
    
    mu=data.frame(mu)
    mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
    
    if(cat_flag){
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        if(length(cat.indx)==1){
          mo[,i]=factor(mo[,i],levels=levels(Ycat[i]))
        }
        else{
          mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
        }
      }
      mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
      mumo[,cat.indx]=mo
    }
    
    mumo[,cont.indx]=mu
    colnames(mumo)=colnames(Y)
    
    S_old=S
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      # var.weights in gower.dist allows for weighted distance
      
      # E step
      V=gower.dist(Y,mumo)^2
      V1=1/V
      V1_lambda=1/(V+jump_penalty)
      S_til=V1_lambda/rowSums(V1_lambda)
      S=matrix(0,nrow=n_obs,ncol=n_states)
      # FUZZY
      S[n_obs,]=V1[n_obs,]/sum(V1[n_obs,])
      for(t in (n_obs-1):1){
        
        # beta=(2*(1-jump_penalty*sum(S[(t+1),]/(V[t,]+jump_penalty))))/
        #   (sum(V1_lambda[t,]))
        # S[t,]=(beta+2*jump_penalty*S[(t+1),])/(2*V[t,]+2*jump_penalty)
        
        S[t,]=S_til[t,]-
          jump_penalty*sum(S[(t+1),]/(V[t,]+jump_penalty))*S_til[t,]+
          jump_penalty*S[(t+1),]/(V[t,]+jump_penalty)
        
      }
      
      #S=t(apply(S,1,function(x) x/sum(x)))
      
      loss <- min(V[1,])
      
      # M step
      for(k in 1:n_states){
        #mu[k,]=apply(Ycont, 2, function(x) weighted_median(x, weights = S[,k]))
        mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=S[,k])})
        
        if(cat_flag){
          if(n_cat==1){
            mo[k,]=poliscidata::wtd.mode(Ycat,weights=S[,k])
            
          }
          else{
            mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=S[,k])})
            
          }
          #mo[k,]=apply(Ycat,2,function(x)weighted_mode(x,weights=S[,k]))
        }
      }
      
      if(cat_flag){
        mumo[,cat.indx]=mo
      }
      
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(Y)
      
      
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
      }
      if (!is.null(tol)) {
        epsilon <- loss_old - loss
        if (epsilon < tol) {
          break
        }
      } else if (all(S == S_old)) {
        break
      }
      loss_old <- loss
      S_old=S
    }
    if (is.null(best_S) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_S <- S
    }
    s=initialize_states(Y,n_states)
  }
  
  MAP=apply(best_S,1,which.max)
  res_Y=data.frame(Y,MAP=MAP)
  col_sort=as.integer(names(sort(tapply(res_Y[,cont.indx[1]],res_Y$MAP,mean))))
  mumo=mumo[col_sort,]
  best_S=best_S[,col_sort]
  
  return(list(best_S=best_S,
              MAP=MAP,
              Y=Y,
              condMM=mumo))
}

simstud_fuzzyJM=function(seed,lambda,TT,P,
                         K=3,mu=1,
                         phi=.8,rho=0,
                         Pcat=NULL,pers=.95,
                         pNAs=0,typeNA=2){
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                        TT=TT,
                        P=P,
                        Ktrue=K,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers)
  
  Y=simDat$SimData.complete
  # Estimate
  est=fuzzy_jump(Y, 
                 n_states=K, jump_penalty=lambda, 
                 initial_states=NULL,
                 max_iter=10, n_init=10, tol=NULL, 
                 verbose=FALSE)
  
  # MAP=apply(est$best_S,1,which.max)
  # 
  ARI=adj.rand.index(est$MAP,simDat$mchain)
  
  # Return
  return(list(
    S=est$best_S,
    MAP=est$MAP ,
    ARI=ARI
    #,
    # seed=seed,
    # lambda=lambda,
    # TT=TT,
    # P=P,
    # Ktrue=K,
    # mu=mu,
    # phi=phi,
    # rho=rho,
    # Pcat=Pcat,
    # pers=pers,
    # true_seq=simDat$mchain,
    # est_seq=MAP
  ))
  
}


# Spatio-temporal ---------------------------------------------------------


