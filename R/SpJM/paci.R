library(MASS)  # For multivariate normal sampling
library(spdep)  # For spatial correlation

# Function to simulate observations given latent states sequence
simulate_observations <- function(mu=1,rho=0,phi=.8,n_states=3,P,Pcat,M,s,seed,pNAs=0,tpNA=0) {
  
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
  
  if(pNAs>0){
    SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=tpNA)
    SimData.NA=as.data.frame(SimData.NA)
    SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
    SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  }
  
  else{
    SimData.NA=SimData
  }
  
  return(list(SimData=SimData,SimData.NA=SimData.NA))
  
}

# Function to create spatial points
generate_spatial_points <- function(n, max_distance = 10) {
  x <- runif(n, 0, max_distance)
  y <- runif(n, 0, max_distance)
  return(data.frame(x = x, y = y))
}

# Function to generate spatio-temporal data with spatial (theta) and temporal (rho) persistence
generate_spatio_temporal_data <- function(M, TT, theta, beta, K = 3,
                                          mu=1,rho=0,phi=.8,
                                          P,Pcat,seed) {
  
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
                        s=S[1,],seed=seed+seed*1000)
  
  Y=temp$SimData
  Y.NA=temp$SimData.NA
  
  Y=data.frame(Y)
  Y$m=1:M
  Y$t=rep(1,M)
  
  S[1,]=order_states_condMean(Y[Y$t==1,dim(Y)[2]-2],S[1,])
  
  Y.NA=data.frame(Y.NA)
  Y.NA$m=1:M
  Y.NA$t=rep(1,M)
  
  # Generate data for each subsequent time step
  for (t in 2:TT) {
    eta_t <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)  # Spatial noise
    data[t, ] <- beta * data[t-1, ] + eta_t  # Temporal correlation
    
    cluster_levels <- quantile(data[t,], probs = seq(0, 1, length.out = K + 1))
    S[t,] <- cut(data[t,], breaks = cluster_levels, labels = FALSE,
                 include.lowest =T)
    
    simDat=simulate_observations(mu=mu,rho=rho,phi=phi,
                               n_states=K,P=P,Pcat=Pcat,M=M,
                               s=S[t,],seed=seed+seed*1000+t-1)
    
    temp=data.frame(simDat$SimData)
    temp$m=1:M
    temp$t=rep(t,M)
    Y=rbind(Y,temp)
    
    temp=data.frame(simDat$SimData.NA)
    temp$m=1:M
    temp$t=rep(t,M)
    Y.NA=rbind(Y.NA,temp)
    
    S[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S[t,])
    
    #data[, t] <- beta * data[, t - 1] + (1-beta) * eta_t
    # cluster_levels <- quantile(data[,t], probs = seq(0, 1, length.out = K + 1))
    # clusters[,t] <- cut(data[,t], breaks = cluster_levels, labels = FALSE,
    #                     include.lowest =T)
    
  }
  
  # matplot(t(data[1:10,]),type='l')
  # hist(data[,4])
  # abline(v=cluster_levels,col='red')
  
  # Discretize data into clusters
  
  # cluster_levels <- quantile(data, probs = seq(0, 1, length.out = K + 1))
  # clusters <- apply(data, 2, function(col) cut(col, breaks = cluster_levels, labels = FALSE,
  #                                              include.lowest =T))
  
  return(list(data = data, 
              S = S, 
              Y=Y,
              Y.NA=Y.NA,
              spatial_points = spatial_points,
              spatial_cov = spatial_cov))
}

# Example usage
M <- 100  # Number of spatial locations
TT <- 4   # Number of time steps
theta <- .01 # Spatial persistence
beta <- .5    # Temporal persistence
K=3
P=20
Pcat=10
result <- generate_spatio_temporal_data(M, TT, theta, beta, K = 3,
                                        mu=1,rho=0,phi=.8,
                                        P=P,Pcat=Pcat,seed=1)

plot(result$spatial_points,col=result$S[1,],pch=19,cex=1.5)
plot(result$spatial_points,col=result$S[2,],pch=19,cex=1.5)
plot(result$spatial_points,col=result$S[3,],pch=19,cex=1.5)
plot(result$spatial_points,col=result$S[4,],pch=19,cex=1.5)
matplot(as.data.frame(result$data),type='l')
