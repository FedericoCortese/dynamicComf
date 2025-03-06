
library(cluster)
library(StatMatch)
library(poliscidata) # for weighted mode 
library(pdfCluster) # for ARI

library(parallel)
library(foreach)
library(doParallel)

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


fuzzy_jump <- function(Y, 
                       n_states, jump_penalty=1e-5, 
                       initial_states=NULL,
                       max_iter=10, n_init=10, tol=NULL, 
                       verbose=FALSE
                       # ,
                       # alpha=NULL
                       # # ,
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
      

      # E step
      V=gower.dist(Y,mumo)^2
      # if(!is.null(alpha)){
      #   V=V-alpha
      # }
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
        #mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=S[,k])})
        mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=S[,k]^2)})
        if(cat_flag){
          if(n_cat==1){
            #mo[k,]=poliscidata::wtd.mode(Ycat,weights=S[,k])
            mo[k,]=poliscidata::wtd.mode(Ycat,weights=S[,k]^2)
            
          }
          else{
            #mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=S[,k])})
            mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=S[,k]^2)})
            
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

relabel_clusters <- function(predicted, true) {
  mapping <- apply(table(predicted, true), 1, which.max)
  sapply(predicted, function(x) mapping[x])
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
  # est=fuzzy_jump(Y, 
  #                n_states=K, jump_penalty=lambda, 
  #                initial_states=NULL,
  #                max_iter=10, n_init=10, tol=NULL, 
  #                verbose=FALSE)
  success <- FALSE
  trials=1
  while (!success&trials<10) {
    est <- try(fuzzy_jump(Y, 
                          n_states = K, jump_penalty = lambda, 
                          initial_states = NULL,
                          max_iter = 10, n_init = 10, tol = NULL, 
                          verbose = FALSE), silent = TRUE)
    trials=trials+1
    
    if (!inherits(est, "try-error")) {
      success <- TRUE  # Exit the loop if no error
    } else {
      message("Retrying fuzzy_jump() due to an error...")
    }
  }
  
  est$MAP=as.factor(relabel_clusters(est$MAP,simDat$mchain))
  simDat$mchain=as.factor(simDat$mchain)
  
  BAC=caret::confusionMatrix(est$MAP,simDat$mchain)$overall[1]
  
  ARI=adj.rand.index(est$MAP,simDat$mchain)
  
  # Return
  return(list(
    S=est$best_S,
    MAP=est$MAP ,
    ARI=ARI,
    BAC=BAC,
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
    true_seq=simDat$mchain
    # est_seq=MAP
  ))
  
}


# Spatio-temporal ---------------------------------------------------------

simulate_observations <- function(mu=1,rho=0,phi=.8,n_states=3,P,Pcat,M,s,seed
                                  #,pNAs,typeNA=0
) {
  
  # This function simulates data from a multivariate normal distribution given the latent states sequence
  
  # Arguments:
  # mu: Mean values for data simulation (first state has mean = mu, last state has mean = -mu, and all intermediate states are equally spaced between them)
  # rho: Correlation between variables
  # n_states: Number of states
  # P: Number of features
  # Pcat: Number of categorical features
  # M: Number of spatial points
  # s: Latent states sequence
  # seed: Random seed
  # pNAs: Percentage of missing values
  # tpNA: Type of missing values (0 = random missing pattern, 1 = block (continuous) missing pattern)
  
  
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
  MU=mu
  mu=seq(-mu,mu,length.out=n_states)
  
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, M, P * n_states)
  SimData = matrix(0, M, P)
  
  set.seed(seed)
  for(k in 1:n_states){
    u = MASS::mvrnorm(M,rep(mu[k],P),Sigma)
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:M) {
    k = s[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
  }
  
  if(Pcat!=0){
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=s,mu=MU,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  # if(pNAs>0){
  #   SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=typeNA)
  #   SimData.NA=as.data.frame(SimData.NA)
  #   SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
  #   SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  # }
  
  #else{
  SimData.NA=SimData
  #}
  
  return(
    #list(
    SimData=SimData
    #,SimData.NA=SimData.NA
    #)
  )
  
}

# Function to create spatial points
generate_spatial_points <- function(n, max_distance = 10) {
  x <- runif(n, 0, max_distance)
  y <- runif(n, 0, max_distance)
  return(data.frame(x = x, y = y))
}

generate_spatio_temporal_data <- function(M, TT, theta, beta, K = 3,
                                          mu=1,rho=0,phi=.8,
                                          P,Pcat,seed,
                                          pGap=.2,
                                          pNAs=0) {
  
  
  # Function to generate spatio-temporal data with spatial and temporal persistence
  
  # Arguments:
  # M: Number of spatial points
  # TT: Number of time points
  # theta: Spatial correlation parameter (the lower, the higher the spatial correlation)
  # beta: Temporal correlation parameter (the higher, the higher the temporal correlation)
  # K: Number of states (only K=3 available at the moment)
  # mu: Mean values for data simulation (first state has mean = mu, last state has mean = -mu, and all intermediate states are equally spaced between them)
  # rho: Correlation between variables
  # phi: Conditional probability for the categorical outcome k in state k
  # P: Number of features
  # Pcat: Number of categorical features
  # seed: Random seed
  # pGap: Percentage of time points to be removed 
  # pNAs: Percentage of missing values (only random missing pattern is available at the moment)
  
  # Value:
  # A list with the following elements:
  # S: A matrix TTxM with the simulated states
  # Y: A data frame with the complete simulated data in long format
  # Y.NA: A data frame with the simulated data with missing values in long format
  # spatial_points: A data frame with the spatial points
  # spatial_cov: The spatial covariance matrix
  # dist_matrix: The distance matrix
  
  
  
  # Increment TT by one as the first time step will be removed
  TT=TT+1
  
  
  # Generate spatial points
  spatial_points <- generate_spatial_points(M)
  
  # Create spatial covariance matrix based on distance and theta
  dist_matrix <- as.matrix(dist(spatial_points)) # Eventually substitute with Gower distance
  
  spatial_cov <- exp(-theta * dist_matrix)
  
  # Initialize data matrix
  data <- array(0, dim = c(TT, M))
  
  S=matrix(0,TT,M)
  
  # Initial time step data (from spatial process)
  data[1, ] <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)
  
  cluster_levels <- quantile(data[1,], probs = seq(0, 1, length.out = K + 1))
  S[1,] <- cut(data[1,], breaks = cluster_levels, labels = FALSE,
               include.lowest =T)
  
  temp=simulate_observations(mu=mu,rho=rho,phi=phi,
                             n_states=K,P=P,Pcat=Pcat,M=M,
                             s=S[1,],seed=seed+seed*1000
                             #,pNAs=pNAs,typeNA=0
  )
  
  Y=temp#$SimData
  #Y.NA=temp$SimData.NA
  
  Y=data.frame(Y)
  Y$m=1:M
  Y$t=rep(0,M)
  
  S[1,]=order_states_condMean(Y[Y$t==0,dim(Y)[2]-2],S[1,])
  
  # Y.NA=data.frame(Y.NA)
  # Y.NA$m=1:M
  # Y.NA$t=rep(0,M)
  
  # Generate data for each subsequent time step
  for (t in 2:TT) {
    eta_t <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)  # Spatial noise
    data[t, ] <- beta * data[t-1, ] + eta_t  # Temporal correlation
    
    #data[t, ] <- beta * data[t-1, ] + (1-beta)*eta_t  # Temporal correlation
    
    cluster_levels <- quantile(data[t,], probs = seq(0, 1, length.out = K + 1))
    S[t,] <- cut(data[t,], breaks = cluster_levels, labels = FALSE,
                 include.lowest =T)
    
    simDat=simulate_observations(mu=mu,rho=rho,phi=phi,
                                 n_states=K,P=P,Pcat=Pcat,M=M,
                                 s=S[t,],seed=seed+seed*1000+t-1
                                 #,pNAs=pNAs,typeNA=0
    )
    
    temp=data.frame(simDat#$SimData
    )
    temp$m=1:M
    temp$t=rep(t-1,M)
    Y=rbind(Y,temp)
    
    # temp=data.frame(simDat$SimData.NA)
    # temp$m=1:M
    # temp$t=rep(t-1,M)
    #Y.NA=rbind(Y.NA,temp)
    
    S[t,]=order_states_condMean(Y[Y$t==(t-1),dim(Y)[2]-2],S[t,])
    
  }
  
  data=data[-1,]
  S=S[-1,]
  Y=Y[-which(Y$t==0),]
  #Y.NA=Y.NA[-which(Y.NA$t==0),]
  
  Y <- Y %>% relocate(t,m)
  #Y.NA <- Y.NA %>% relocate(t,m)
  
  Y.NA=apply(Y[,-(1:2)],2,punct,pNAs=pNAs,type=0)
  Y.NA=as.data.frame(Y.NA)
  Y.NA[,1:Pcat]=Y.NA[,1:Pcat]%>%mutate_all(as.factor)
  Y.NA[,-(1:Pcat)]=Y.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  Y.NA=data.frame(t=Y$t,m=Y$m,Y.NA)
  
  if(pGap>0){
    set.seed(seed)
    gaps=sort(sample(1:TT,round(TT*pGap),replace=F))
    Y=Y[-which(Y$t %in% gaps),]
    Y.NA=Y.NA[-which(Y.NA$t %in% gaps),]
  }
  
  return(list(S = S, 
              Y=Y,
              Y.NA=Y.NA,
              spatial_points = spatial_points,
              spatial_cov = spatial_cov,
              dist_matrix = dist_matrix)
  )
}

dist_fun_norm=function(D){
  
  dist_weights=D
  for(m in 1:nrow(D)){
    dist_weights[m,]=ifelse(D[m,] == 0, 0, D[m,])/sum(ifelse(D[m,] == 0, 0, D[m,])) 
  }
  return(dist_weights)
  
}

fuzzy_STJM=function(Y,
                    K,
                    dist_weights,
                    lambda,
                    gamma,
                    n_init=10,
                    max_iter=10,
                    tol=NULL,
                    verbose=F,
                    ncores_M=NULL){
  
  # dist_weights is a matrix with the weights for the distance between the spatial points
  
  Y=Y[order(Y$t,Y$m),]
  P=ncol(Y)-2
  TT=length(unique(Y$t))
  M=length(unique(Y$m))
  best_loss <- NULL
  best_S <- NULL
  
  Y.orig=Y
  
  Y=subset(Y,select=-c(t,m))
  
  # Scale numeric variables
  Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  
  Y <- data.frame(t=Y.orig$t,m=Y.orig$m,Y)
  
  # Reorder columns so that we have t,m, cat vars and cont vars
  cat.indx=which(sapply(Y, is.factor))
  cat.flag=T
  if(length(cat.indx)==0){
    cat.flag=F
  }
  cont.indx=which(sapply(Y, is.numeric))
  cont.indx=cont.indx[-(1:2)]
  Y=Y[,c(1,2,cat.indx,cont.indx)]
  
  YY=subset(Y,select=-c(t,m))
  
  cont.indx=cont.indx-2
  cat.indx=cat.indx-2
  
  Ycont=YY[,cont.indx]
  
  Ycat=YY[,cat.indx]
  levels_cat=lapply(Ycat,levels)
  names(levels_cat)=cat.indx
  
  n_cat=length(cat.indx)
  n_cont=length(cont.indx)
  
  # Missing data imputation 
  # TO DO
  # Mcont=ifelse(is.na(Ycont),T,F)
  # Mcat=ifelse(is.na(Ycat),T,F)
  # mu <- apply(Ycont, 2, median, na.rm = TRUE)
  # #mu <- colMeans(Ycont,na.rm = T)
  # mo <- apply(Ycat,2,Mode)
  # for(i in 1:n_cont){
  #   Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  # }
  # 
  # if(cat.flag){
  #   for(i in 1:n_cat){
  #     x=Ycat[,i]
  #     Ycat[which(is.na(Ycat[,i])),i]=mo[i]
  #   }
  # }
  # 
  # YY[,cont.indx]=Ycont
  # YY[,cat.indx]=Ycat
  # 
  # Y[,-(1:2)]=YY
  
  for (init in 1:n_init) {
    
    loss_old <- 1e10
    
    # State initialization
    S=matrix(0,nrow=TT,ncol=M)
    for(m in 1:M){
      S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],K)
    }
    vec_S=as.vector(t(S))
    
    SS <- matrix(0, nrow = TT * M, ncol = 2 + K)
    SS[, 1] <- rep(1:TT, each = M)  # Time indices
    SS[, 2] <- rep(1:M, times = TT) # Spatial indices
    
    # Compute row indices in SS corresponding to (t, m)
    row_indices <- rep(1:(TT * M))  # Row positions in SS
    col_indices <- vec_S + 2  # Convert S into a vector and shift for SS indexing
    
    # Assign 1s in a single step
    SS[cbind(row_indices, col_indices)] <- 1 
    
    # Convert to DataFrame and set column names
    SS <- data.frame(SS)
    colnames(SS) <- c("t", "m", 1:K)
    
    SS_old=SS
    
    mu <- matrix(0, nrow=K, ncol=length(cont.indx))
    mo <- matrix(0, nrow=K, ncol=length(cat.indx))
    
    for (i in unique(vec_S)) {
      # substitute with medians
      mu[i,] <- apply(Ycont[vec_S==i,], 2, median, na.rm = TRUE)
      mo[i,]=apply(Ycat[vec_S==i,],2,Mode)
    }
    
    mu=data.frame(mu)
    mo=data.frame(mo,stringsAsFactors=TRUE)
    for(i in 1:n_cat){
      mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
    }
    mumo=data.frame(matrix(0,nrow=K,ncol=P))
    mumo[,cat.indx]=mo
    mumo[,cont.indx]=mu
    colnames(mumo)=colnames(YY)
    
    for (it in 1:max_iter) {
      
      # E Step
      g=gower.dist(YY,mumo)
      
      if(is.null(ncores_M)){
        for (m in 1:M){
          indx=which(Y$m==m)
          omega_m=dist_weights[m,]
          Wm=sum(omega_m)
          
          g2_1=1/(g[indx,]^2+gamma*Wm+lambda)
          g2_1[1,]=1/(g[indx[1],]^2+gamma*Wm)
          
          s_til_m=g2_1/rowSums(g2_1)
          
          s_i=as.matrix(SS[SS$t==1,-(1:2)])
          
          A=sum((omega_m%*%s_i)/g2_1[1,])
          B=(omega_m%*%s_i)/g2_1[1,]
          
          SS[indx[1],-(1:2)]=s_til_m[1,]-gamma*s_til_m[1,]*A+gamma*B
          
          for(t in 2:TT){
            s_i=as.matrix(SS[SS$t==t,-(1:2)])
            A=sum((omega_m%*%s_i)/g2_1[t,])
            B=(omega_m%*%s_i)/g2_1[t,]
            SS[indx[t],-(1:2)]=
              s_til_m[t,]-
              gamma*(B-s_til_m[t,]*A)+
              lambda*(SS[indx[t-1],-(1:2)]/g2_1[t,]-
                        s_til_m[t,]*sum(SS[indx[t-1],-(1:2)]/g2_1))
            # Potrebbero non sommare a  uno, verificare
            # Intanto normalizzo
            SS[indx[t],-(1:2)]=SS[indx[t],-(1:2)]/sum(SS[indx[t],-(1:2)])
          }
        }
      }
      else{
        num_cores <- ncores_M  
        cl <- makeCluster(num_cores)
        registerDoParallel(cl)
        # SS <- foreach(m = 1:M
        #               #, .combine = '+'
        #               ) %dopar% {
        #   indx=which(Y$m==m)
        #   omega_m=dist_weights[m,]
        #   Wm=sum(omega_m)
        #   
        #   g2_1=1/(g[indx,]^2+gamma*Wm+lambda)
        #   g2_1[1,]=1/(g[indx[1],]^2+gamma*Wm)
        #   
        #   s_til_m=g2_1/rowSums(g2_1)
        #   
        #   s_i=as.matrix(SS[SS$t==1,-(1:2)])
        #   
        #   A=sum((omega_m%*%s_i)/g2_1[1,])
        #   B=(omega_m%*%s_i)/g2_1[1,]
        #   
        #   SS[indx[1],-(1:2)]=s_til_m[1,]-gamma*s_til_m[1,]*A+gamma*B
        #   
        #   for(t in 2:TT){
        #     s_i=as.matrix(SS[SS$t==t,-(1:2)])
        #     A=sum((omega_m%*%s_i)/g2_1[t,])
        #     B=(omega_m%*%s_i)/g2_1[t,]
        #     SS[indx[t],-(1:2)]=
        #       s_til_m[t,]-
        #       gamma*(B-s_til_m[t,]*A)+
        #       lambda*(SS[indx[t-1],-(1:2)]/g2_1[t,]-
        #                 s_til_m[t,]*sum(SS[indx[t-1],-(1:2)]/g2_1))
        #     # Potrebbero non sommare a  uno, verificare
        #     # Intanto normalizzo
        #     SS[indx[t],-(1:2)]=SS[indx[t],-(1:2)]/sum(SS[indx[t],-(1:2)])
        #   }
        #   return(as.matrix(SS))
        # }
        SS <- foreach(m = 1:M, .combine = rbind, 
                             .packages = "matrixStats") %dopar% {
          indx <- which(Y$m == m)
          omega_m <- dist_weights[m, ]
          Wm <- sum(omega_m)
          
          g2_1 <- 1 / (g[indx, ]^2 + gamma * Wm + lambda)
          g2_1[1, ] <- 1 / (g[indx[1], ]^2 + gamma * Wm)
          
          s_til_m <- g2_1 / rowSums(g2_1)
          
          # Create a local copy of SS for this iteration
          SS_local <- SS  # Assuming SS is pre-initialized and accessible
          
          s_i <- as.matrix(SS_local[SS_local$t == 1, -(1:2)])
          
          A <- sum((omega_m %*% s_i) / g2_1[1, ])
          B <- (omega_m %*% s_i) / g2_1[1, ]
          
          SS_local[indx[1], -(1:2)] <- s_til_m[1, ] - gamma * s_til_m[1, ] * A + gamma * B
          
          for (t in 2:TT) {
            s_i <- as.matrix(SS_local[SS_local$t == t, -(1:2)])
            A <- sum((omega_m %*% s_i) / g2_1[t, ])
            B <- (omega_m %*% s_i) / g2_1[t, ]
            
            SS_local[indx[t], -(1:2)] <-
              s_til_m[t, ] -
              gamma * (B - s_til_m[t, ] * A) +
              lambda * (SS_local[indx[t - 1], -(1:2)] / g2_1[t, ] -
                          s_til_m[t, ] * sum(SS_local[indx[t - 1], -(1:2)] / g2_1))
            
            # Normalize
            SS_local[indx[t], -(1:2)] <- SS_local[indx[t], -(1:2)] / sum(SS_local[indx[t], -(1:2)])
          }
          
          # Return the data frame for this m
          return(data.frame(t = SS_local$t[indx], m = m, 
                            SS_local[indx, -(1:2)]))
        }
        stopCluster(cl)
        colnames(SS) <- c("t", "m", 1:K)
        # sort SS by t 
        SS=SS[order(SS$t), ]
        
        }
      
      
      # M Step
      for(k in 1:K){
        
        mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=SS[,k+2]^2)})
        mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=SS[,k+2]^2)})
        # mu[k,]=apply(Ycont, 2, function(x) weighted_median(x, weights = SS[,k+2]))
        # mo[k,]=apply(Ycat,2,function(x)weighted_mode(x,weights=SS[,k+2]))
      }
      
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(YY)
      
      # Re-fill-in missings 
      # TO DO
      
      # Compute objective function
      g=gower.dist(YY,mumo)
      # sum_result <- 0
      # for (t in 1:TT) {
      #   for (m in 1:M) {
      #     s_t_m <- SS[SS$t==t&SS$m==m,-(1:2),]
      # 
      #     for (i in 1:M) {
      #       if (i != m) {
      #         s_t_i <- SS[SS$t==t&SS$m==i,-(1:2),]
      # 
      #         diff_norm_sq <- sum((s_t_i - s_t_m)^2)
      #         sum_result <- sum_result + dist_weights[m, i] * diff_norm_sq
      # 
      #       }
      #     }
      #   }
      # }
      
      value_opt=sum(g^2)+
        lambda*sum((SS[-which(SS$t==1),-(1:2)]-SS[-which(SS$t==TT),-(1:2)])^2)
      #+gamma*sum_result
      
      if (verbose) {
        cat("Init.: ", init, "// Iteration: ", it, "// Obj. Fun.: ", value_opt, "\n")
      }
      
      if (!is.null(tol)) {
        if (sum((SS-SS_old)^2)<tol) {
          break
        }
        else if( (loss_old - value_opt) < tol){
          break
        }
      } 
      
      SS_old=SS
      
      loss_old <- value_opt
    }
    if (is.null(best_S)||(loss_old < best_loss)
        ) {
      best_loss <- loss_old
      best_S <- SS
    }
    
  }
  
  return(list(S=best_S))
  
}


simstud_fuzzySTJM=function(lambda,gamma,seed,M,TT,beta, theta,
                             mu=.5,rho=0.2,
                             K=3,P=20,phi=.8,Pcat=NULL,pNAs=0,pg=0,
                           ncores_M=3){
  
  library(MASS)
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
  result <- generate_spatio_temporal_data(M, TT, theta, beta, K = K,
                                          mu=mu,rho=rho,phi=phi,
                                          P=P,Pcat=Pcat,seed=seed,pGap=pg,pNAs=pNAs)
  
  Y.compl=result$Y
  D=result$dist_matrix
  Y=result$Y.NA
  
  dist_weights=dist_fun_norm(D)
  
  # tf=I(pg>0)
  prova=fuzzy_STJM(Y,K,dist_weights,lambda,gamma,verbose=F,tol=1e-6,
                   ncores_M=ncores_M,
                   n_init = 5)
  
  best_S=prova$S
  
  MAP=data.frame(t=Y$t,m=Y$m,
    MAP=apply(best_S[,-c(1,2)],1,which.max))
  
  
  df <- as.data.frame(result$S)
  colnames(df) <- paste0("m", 1:ncol(df))
  df$t <- 1:nrow(df)
  df_long <- tidyr::pivot_longer(df, cols = starts_with("m"), 
                                 names_to = "m", 
                                 values_to = "true_state")
  df_long$m <- as.numeric(gsub("m", "", df_long$m))
  df_long=as.data.frame(df_long)
  
  # TY=unique(Y$t)
  # S_true=result$S[TY,]
  # 
  # for(t in 1:length(TY)){
  #   best_s[t,]=order_states_condMean(Y[Y$t==TY[t],dim(Y)[2]],best_s[t,])
  # }
  # 
  ARI=adj.rand.index(df_long$true_state,MAP$MAP)
  ARI
  return(list(ARI=ARI,
              S_true=df_long,
              best_S=best_S,
              lambda=lambda,
              gamma=gamma,
              seed=seed,
              M=M,TT=TT,P=P,
              lambda=lambda,
              gamma=gamma
              # ,
              # mu=mu,rho=rho,
              # K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs,pg=pg
              ))
}
