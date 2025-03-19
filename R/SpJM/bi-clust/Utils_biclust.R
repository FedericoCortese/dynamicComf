library(RcppHMM)
library(reticulate)
library(pdfCluster)
library(caret)
library(boot)
library(xtable)
library(dplyr)
library(cluster)
library(gower)
library(StatMatch)
library(SpectralClMixed)
library(multiUS)
#library(missForest)
#library(MCMCprecision)
#require("potts")
library(missMethods)
library(aricode)
#py_install("scipy")
library(mvtnorm)
library(GA)

library(parallel)
library(foreach)
library(doParallel)

library(MASS)  # For multivariate normal sampling
library(spdep)  # For spatial correlation

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

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

discretize_prob_simplex <- function(n_c, grid_size) {
  # Sample grid points on a probability simplex.
  N <- as.integer(1 / grid_size)
  
  # Generate all combinations and filter those that sum to N
  tuples <- expand.grid(rep(list(0:N), n_c))
  valid_tuples <- tuples[rowSums(tuples) == N, ]
  
  # Reverse the order and scale to get the simplex points
  lst <- as.matrix(valid_tuples[nrow(valid_tuples):1, ]) / N
  rownames(lst) <- NULL
  return(lst)
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

order_states=function(states){
  
  # This function organizes states by assigning 1 to the first observed state and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  # states is a vector of observed states
  
  N=length(states)
  states_temp=rep(0,N)
  new=1
  states_temp[1]=new
  for(i in 2:N){
    if(sum(states[i]==states[1:(i-1)])==0){
      # we enter this is if-stat. whenever a new state appeares
      states_temp[i]=new+1
      new=new+1
    }
    else{
      states_temp[i]=states_temp[which(states[1:(i-1)]==states[i])[1]]
    }
  }
  return(states_temp)
}

order_states_freq=function(states){
  
  # This function organizes states by assigning 1 to the mostly observed state and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  # states is a vector of observed states
  
  states_temp=match(states,names(sort(table(states))))
  
  return(states_temp)
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

get_cat=function(y,mc,mu,phi){
  # Function to simulate categorical data
  
  # Arguments:
  # y: continuous variable 
  # mc: Markov chain states
  # mu: numeric mean value
  # phi: conditional probability for the categorical outcome k in state k
  
  library(dplyr)
  
  mu=c(-mu,0,mu)
  #K=length(unique(mc))
  #mu=seq(-mu,mu,length.out=K)
  #phi1=(1-phi)/(K-1)
  phi1=(1-phi)/2
  
  TT=length(y)
  for(i in 1:TT){
    k=mc[i]
    switch(k,
           "1"={
             threshold=c(qnorm(phi1,mu[1]),qnorm(phi+phi1,mu[1]))
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
             threshold=c(qnorm(phi1,mu[2]),qnorm(phi+phi1,mu[2]))
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
             threshold=c(qnorm(phi1,mu[3]),qnorm(phi+phi1,mu[3]))
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
# Function to reassign labels using the Hungarian method
reassign_labels <- function(true_labels, predicted_labels) {
  unique_true_labels <- unique(true_labels)
  unique_predicted_labels <- unique(predicted_labels)
  
  cost_matrix <- compute_cost_matrix(true_labels, predicted_labels)
  
  # Ensure the cost matrix is square by adding dummy clusters if necessary
  if (nrow(cost_matrix) > ncol(cost_matrix)) {
    cost_matrix <- cbind(cost_matrix, matrix(0, nrow = nrow(cost_matrix), ncol = nrow(cost_matrix) - ncol(cost_matrix)))
  } else if (nrow(cost_matrix) < ncol(cost_matrix)) {
    cost_matrix <- rbind(cost_matrix, matrix(0, ncol = ncol(cost_matrix), nrow = ncol(cost_matrix) - nrow(cost_matrix)))
  }
  
  # Convert cost matrix to distance matrix
  distance_matrix <- max(cost_matrix) - cost_matrix
  
  # Solve the assignment problem
  assignment <- solve_LSAP(distance_matrix, maximum = FALSE)
  
  # Create a mapping based on the optimal assignment
  new_labels <- predicted_labels
  for (i in seq_along(unique_predicted_labels)) {
    new_labels[predicted_labels == unique_predicted_labels[i]] <- unique_true_labels[assignment[i]]
  }
  
  return(new_labels)
}

# Function to simulate observations given latent states sequence
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

onerun_biclust_JM=function(Y,n_states,
                          jump_penalty,
                          alpha,
                          grid_size=.05,
                          rand_search_sample=100,
                          mode_loss=T,
                          max_iter,tol,initial_states, 
                          Mcont, Mcat,cont.indx,cat.indx,levels_cat,
                          ncores_M=ncores_M){
  #tryCatch({
  
  YY=Y[,-(1:2)]
  Ycat=YY[,cat.indx]
  Ycont=YY[,cont.indx]
  n_cat=length(cat.indx)
  n_cont=length(cont.indx)
  TT=length(unique(Y$t))
  M=length(unique(Y$m))
  
  
  loss_old <- 1e10
  
  
  # State initialization
  S=matrix(0,nrow=TT,ncol=M)
  for(m in 1:M){
    S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],n_states)
  }
  # I need vectorized S many times, let's do it only once
  vec_S=as.vector(t(S))
  
  # Initialize soft clustering matrix
  # SLOW
  # SS <- matrix(0, nrow = TT * M, ncol = 2 + n_states)
  # SS[, 1] <- rep(1:TT, each = M)  # First column: t indices
  # SS[, 2] <- rep(1:M, times = TT) # Second column: m indices
  # for (t in 1:TT) {
  #   for (m in 1:M) {
  #     state <- S[t, m]
  #     SS[(t - 1) * M + m, 2 + state] <- 1  # Set the appropriate column to 1
  #   }
  # }
  # SS=data.frame(SS)
  # colnames(SS)=c("t","m",1:n_states)
  
  # FASTER
  SS <- matrix(0, nrow = TT * M, ncol = 2 + n_states)
  
  SS[, 1] <- rep(1:TT, each = M)  # Time indices
  SS[, 2] <- rep(1:M, times = TT) # Spatial indices
  
  # Compute row indices in SS corresponding to (t, m)
  row_indices <- rep(1:(TT * M))  # Row positions in SS
  col_indices <- vec_S + 2  # Convert S into a vector and shift for SS indexing
  
  # Assign 1s in a single step
  SS[cbind(row_indices, col_indices)] <- 1 
  
  # Convert to DataFrame and set column names
  SS <- data.frame(SS)
  colnames(SS) <- c("t", "m", 1:n_states)
  
  mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
  mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
  
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
  mumo=data.frame(matrix(0,nrow=n_states,ncol=P))
  mumo[,cat.indx]=mo
  mumo[,cont.indx]=mu
  colnames(mumo)=colnames(YY)
  
  # Should this be moved outside the for loop of max_iter?
  if(!is.null(grid_size)){
    # Grid search
    prob_vecs <- discretize_prob_simplex(n_states, grid_size)
  }
  else{
    # Random search
    # Consider using Dirichlet distribution to put more mass on extremes
    prob_vecs <-matrix(rep(sample_dirichlet_vector(n_states, alpha = .3),rand_search_sample), 
                       nrow = rand_search_sample)
    # prob_vecs <- matrix(runif(n_states * rand_search_sample), nrow = rand_search_sample)
    # prob_vecs=t(apply(prob_vecs,1,function(x)x/sum(x)))
  }
  
  # This should be computed for each m and t differently if we adopt recursive refinement
  pairwise_l1_dist <- as.matrix(dist(prob_vecs, method = "manhattan")) / 2
  jump_penalty_mx <- jump_penalty * (pairwise_l1_dist ^ alpha)
  
  if (mode_loss) {
    # Adding mode loss
    m_loss <- log(rowSums(exp(-jump_penalty_mx)))
    m_loss <- m_loss - m_loss[1]  # Offset a constant
    jump_penalty_mx <- jump_penalty_mx + m_loss
  }
  
  for (it in 1:max_iter) {
    
    # E Step
    
    loss_mx=gower.dist(YY,mumo)
    
    # Slighlty faster
    # loss_mx<- 0.5 * sqrt(outer(1:nrow(YY), 1:nrow(mu), 
    #                             Vectorize(function(r, k) sum((YY[r, ] - mu[k, ])^2))))
    
    
    LARGE_FLOAT=1e1000
    # Handle continuous model with probability vectors
    if (!is.null(prob_vecs)) {
      loss_mx[is.nan(loss_mx)] <- LARGE_FLOAT
      loss_mx <- loss_mx %*% t(prob_vecs) # (TxM)xN
    }
    
    N <- ncol(loss_mx)
    
    loss_mx[is.nan(loss_mx)] <- Inf
    
    values <- matrix(NA, TT*M, N) # (TxM)xN
    value_opt=rep(0,M)
    assign <- integer(TT*M) #TxM
    SS_new=SS
    #SS_new=SS[,-(1:2)]
    SS_new[,-(1:2)]=0
    
    if(is.null(ncores_M)){
      # DP iteration (bottleneck)
      for(m in 1:M){
        
        #print(m)
        indx=which(Y$m==m)
        
        # I pesi sono definiti come le distanze normalizzate tra 0 e 1
        
        
        values[indx[1], ] <- loss_mx[indx[1], ]
        
        # Bootleneck!!
        for (t in 2:TT) {
          #print(t)
          
          values[indx[t], ] <- loss_mx[indx[t], ] + 
            apply(values[indx[t-1], ] + jump_penalty_mx, 2, min)
          
        }
        
        # Find optimal path backwards
        assign[indx[TT]] <- which.min(values[indx[TT], ])
        value_opt[m] <- values[indx[TT], assign[indx[TT]]]
        
        SS_new[indx[TT],-(1:2)]=prob_vecs[assign[indx[TT]],]
        
        # Traceback
        for (t in (TT - 1):1) {
          assign[indx[t]] <- which.min(values[indx[t], ] + 
                                         jump_penalty_mx[, assign[indx[t+1]]])
          SS_new[indx[t],-(1:2)]=prob_vecs[assign[indx[t]],]
        }
        
      }
      # vector assign tells, for each m and t, which one among the N candidate vectors is the best
      value_opt=mean(value_opt)
    }
    else{
      num_cores <- ncores_M  # Use all but one core
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      # Parallel Execution using foreach
      SS_new <- foreach(m = 1:M, .combine = '+') %dopar% {
        
        indx <- which(Y$m == m)
        
        
        values[indx[1], ] <- loss_mx[indx[1], ] 
        
        
        for (t in 2:TT) {
          
          values[indx[t], ] <- loss_mx[indx[t], ] + 
            apply(values[indx[t-1], ] + jump_penalty_mx, 2, min) 
        }
        
        # Find optimal path backwards
        assign[indx[TT]] <- which.min(values[indx[TT], ])
        #value_opt[m] <- values[indx[TT], assign[indx[TT]]]
        value_opt <- values[indx[TT], assign[indx[TT]]]
        
        SS_new=matrix(0,nrow=TT*M,ncol=n_states)
        
        SS_new[indx[TT], ] <- prob_vecs[assign[indx[TT]],]
        
        # Traceback
        for (t in (TT - 1):1) {
          assign[indx[t]] <- which.min(values[indx[t], ] + 
                                         jump_penalty_mx[, assign[indx[t + 1]]])
          SS_new[indx[t],] <- prob_vecs[assign[indx[t]],]
        }
        SS_new=rbind(rep(value_opt,n_states),SS_new)
        
        return(SS_new[,])  # Avoid returning large objects inside parallel execution
      }
      
      # Stop parallel cluster after execution
      stopCluster(cl)
      
      value_opt=SS_new[1,1]/M
      SS_new=SS_new[-1,]
      SS_new=data.frame(t=Y$t,m=Y$m,SS_new)
      colnames(SS_new)=c("t","m",1:n_states)
    }
    
    # M Step
    
    for(k in 1:n_states){
      mu[k,]=apply(Ycont, 2, function(x) weighted_median(x, weights = SS_new[,k+2]))
      mo[k,]=apply(Ycat,2,function(x)weighted_mode(x,weights=SS_new[,k+2]))
    }
    
    mumo[,cat.indx]=mo
    mumo[,cont.indx]=mu
    colnames(mumo)=colnames(YY)
    
    if (!is.null(tol)) {
      epsilon <- loss_old - value_opt
      if (epsilon < tol) {
        break
      }
    } 
    
    else if (all(SS == SS_new)) {
      break
    }
    
    SS=SS_new
    
    loss_old <- value_opt
  }
  
  
  
  return(list(S=SS,value_opt=value_opt,mumo=mumo))
  #    }
  # , 
  #     error = function(e) {
  #       # Return a consistent placeholder on error
  #       return(list(S = NA, value_opt = Inf))
  #     })
}

biclust_JM <- function(Y,K,
                      jump_penalty,
                      alpha=2,grid_size=.05,
                      mode_loss=T,rand_search_sample=100,
                      n_init=10,
                      max_iter=10,tol=NULL,initial_states=NULL,
                      n_cores=NULL,prll=F,ncores_M=NULL
                      #,parallel_ga=F
) {
  
  # n_cores is the number of cores dedicated to parallel computation for different initialization
  # prll is the flag for parallel computation of that
  # parallel_ga is the number of cores dedicated to parallel computation for GA over subpopulations
  # Default is FALSE (no parallelization)
  
  Y=Y[order(Y$t,Y$m),]
  P=ncol(Y)-2
  Y.orig=Y
  
  Y=subset(Y,select=-c(t,m))
  Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  
  Y <- data.frame(t=Y.orig$t,m=Y.orig$m,Y)
  
  # Reorder columns so that we have t,m, cat vars and cont vars
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  cont.indx=cont.indx[-(1:2)]
  Y=Y[,c(1,2,cat.indx,cont.indx)]
  
  YY=subset(Y,select=-c(t,m))
  
  TT=length(unique(Y$t))
  M=length(unique(Y$m))
  
  cat.indx=which(sapply(YY, is.factor))
  cont.indx=which(sapply(YY, is.numeric))
  
  Ycont=YY[,cont.indx]
  
  Ycat=YY[,cat.indx]
  levels_cat=lapply(Ycat,levels)
  names(levels_cat)=cat.indx
  
  n_cat=length(cat.indx)
  n_cont=length(cont.indx)
  
  
  if(prll){
    if(is.null(n_cores)){
      n_cores=parallel::detectCores()-1
    }
    hp_init=expand.grid(init=1:n_init)
    jms <- parallel::mclapply(1:nrow(hp_init),
                              function(x)
                                
                                onerun_biclust_JM(Y=Y,n_states=K,
                                                 jump_penalty=jump_penalty,
                                                 alpha=alpha,grid_size=grid_size,
                                                 rand_search_sample=rand_search_sample,
                                                 mode_loss=mode_loss,
                                                 max_iter=max_iter,tol=tol,
                                                 initial_states=initial_states,
                                                 Mcont=Mcont, Mcat=Mcat,
                                                 cont.indx=cont.indx,cat.indx=cat.indx,levels_cat=levels_cat,
                                                 ncores_M = ncores_M),
                              mc.cores = n_cores)
    # Sostituire con GA_onerun_cont_STJM
  }
  else{
    jms=list()
    for (init in 1:n_init) {
      
      
      
      jms[[init]]=onerun_biclust_JM(Y=Y,n_states=K,
                                   jump_penalty=jump_penalty,
                                   alpha=alpha,grid_size=grid_size,
                                   rand_search_sample=rand_search_sample,
                                   mode_loss=mode_loss,
                                   max_iter=max_iter,tol=tol,
                                   initial_states=initial_states,
                                   Mcont=Mcont, Mcat=Mcat,
                                   cont.indx=cont.indx,cat.indx=cat.indx,levels_cat=levels_cat,
                                   ncores_M = ncores_M)
      
    }
  }
  
  best_init=which.min(unlist(lapply(jms,function(x)x$value_opt)))
  best_S=jms[[best_init]]$S
  best_mumo=jms[[best_init]]$mumo
  
  
  return(best_S)
}