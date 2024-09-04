library(RcppHMM)
library(reticulate)
library(pdfCluster)
library(boot)
library(xtable)
library(dplyr)
library(cluster)
library(gower)
library(StatMatch)
library(SpectralClMixed)
library(multiUS)
library(missForest)
library(parallel)
library(MCMCprecision)
require("potts")
#py_install("scipy")

# Import the python module
 import("scipy")
# 
# # Import python functions for SJM estimation
 source_python('SJ.py')

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

compute_cost_matrix <- function(true_labels, predicted_labels) {
  unique_true_labels <- unique(true_labels)
  unique_predicted_labels <- unique(predicted_labels)
  
  # Create an empty cost matrix
  cost_matrix <- matrix(0, nrow = length(unique_true_labels), ncol = length(unique_predicted_labels))
  
  # Fill the cost matrix
  for (i in seq_along(unique_true_labels)) {
    for (j in seq_along(unique_predicted_labels)) {
      cost_matrix[i, j] <- sum(true_labels == unique_true_labels[i] & predicted_labels == unique_predicted_labels[j])
    }
  }
  
  return(cost_matrix)
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
get_BCD <- function(Y, states) {
  
  n_states=length(unique(states))
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  
  
  Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  
  Ycont=Y[,-cat.indx]
  Ycat=Y[,cat.indx]
  n_cat=length(cat.indx)
  n_cont=n_features-n_cat
  
  # Unconditional means
  mu <- colMeans(Ycont,na.rm = T)
  mu=matrix(mu,nrow=1)
  mu=data.frame(mu)
  # Unconditional modes
  mo <- apply(Ycat,2,Mode)
  mo=matrix(mo,nrow=1)
  mo=data.frame(mo)
  for(i in 1:n_cat){
    x=Ycat[,i]
    mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
  }
  
  mumo=data.frame(matrix(0,nrow=1,ncol=n_features))
  mumo[,cat.indx]=mo
  mumo[,-cat.indx]=mu
  uncondMM=mumo
  colnames(uncondMM)=colnames(Y)
  
  # Find TD (total deviance)
  TD=rep(0,n_features)
  for(j in 1:n_features){
    TD[j]=sum(gower.dist(Y[,j],uncondMM[,j]))
  }
  
  if(length(unique(states))==1){
    WCD=TD
  }
  
  else{
    # Compute conditional means and modes
    mu=matrix(0,nrow=n_states,ncol=n_cont)
    mo=matrix(0,nrow=n_states,ncol=n_cat)
    
    for (i in unique(states)) {
      mu[i,] <- colMeans(Ycont[states==i,],na.rm=T)
      mo[i,] <- as.numeric(apply(Ycat[states==i,],2,Mode,na.rm=T))
    }
    
    mu=data.frame(mu)
    mo=data.frame(mo,stringsAsFactors=TRUE)
    for(i in 1:n_cat){
      x=Ycat[,i]
      #mo[,i]=factor(mo[,i],levels=1:n_levs[i])
      mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
    }
    
    mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
    mumo[,cat.indx]=mo
    mumo[,cont.indx]=mu
    condMM=mumo
    colnames(condMM)=colnames(Y)
    
    
    # The following can be written easier (TBD)
    WCD=matrix(0,ncol=n_states,nrow=n_features)
    
    for (i in unique(states)) {
      for(j in 1:n_features){
        temp=gower.dist(Y[,j],condMM[i,j])
        WCD[j,i]=sum(temp[states==i])
      }    
    }
    WCD=rowSums(WCD)
  }
  
  BCD=TD-WCD
  
  return(list(BCD=BCD,TD=TD,WCD=WCD))
}


# Ordinary JM -------------------------------------------------------------

sim_data=function(seed,Ktrue,N,P,cors,pers,m){
  
  # Function to simulate data from a multivariate Gaussian HMM
  # seed is the seed for the random number generator
  # Ktrue is the number of states
  # N is the number of observations
  # P is the number of features
  # cors is a vector with correlations, one for each state
  # pers is the self-transition probability 
  # m is an integer number such that m*P is the number of false features
  
  # Output:   
  # A list with the true states and the simulated data
  
  states_names = as.character(seq(1,Ktrue))
  
  Sigma <- array(0, dim =c(P,P,Ktrue))
  
  for(i in 1:Ktrue){
    Sigma[,,i] <- matrix(cors[i], ncol = P,  
                         byrow = TRUE)
    diag(Sigma[,,i])=1
  }
  
  set.seed(seed)
  Mu <- matrix(runif(P*Ktrue,min=-2,max=2), 
               nrow = P, 
               byrow = TRUE)
  
  A <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(A)=rep(pers,Ktrue)
  
  # initial probabilities
  Pi <- rep(1/Ktrue,Ktrue)
  
  HMM.cont.multi <- verifyModel(list( "Model" = "GHMM",
                                      "StateNames" = states_names,
                                      "A" = A, 
                                      "Mu" = Mu, 
                                      "Sigma" = Sigma, 
                                      "Pi" = Pi)
                                )
  set.seed(seed)
  observationSequence <- generateObservations(HMM.cont.multi, N)
  Y=t(observationSequence$Y)
  Y=apply(Y,2,scale)
  true_states=order_states(as.factor(observationSequence$X))
  
  set.seed(seed)
  if(m>=1){
    if(m==1){
      Ytil=Y[sample(N),]
    }
    else{
      Ytil=Y[sample(N),]
      for(i in 2:m){
        Ytil=cbind(Ytil,Y[sample(N),])
      } 
    }
    YY=cbind(Y,Ytil) 
  }
  else{
    YY=Y
  }
  
  return(list(true_states=true_states,
              YY=YY))
  
}

sim_est_SJM=function(seed,K,lambda,kappa,true_data,Ktrue,
                     N=1000,P=100,m=2,cors,pers,K0,alpha0){
  
  # Function to estimate the parameters of a sparse jump model
  # seed is the seed for the random number generator
  # K is the number of states
  # lambda is the penalty for the jump
  # kappa is the penalty for sparsity
  # true_data is obtained using the function "sim_data"
  # Ktrue is the true number of states
  # N is the number of observations
  # P is the number of features
  # m is an integer number such that m*P is the number of false features
  # cors is a vector with correlations, one for each state
  # pers is the self-transition probability
  # K0 is the prior number of states
  
  
  # Output: 
  # A list with the estimated weights,
  # the normalized estimated weights,
  # the estimated states,
  # the seed,
  # the number of observations,
  # the number of features,
  # the number of false features,
  # the correlations,
  # the self-transition probability,
  # the true number of states,
  # the penalty for the jump,
  # the penalty for sparsity,
  # the L1 norm of the BCSS,
  # the L2 norm of the BCSS,
  # the penalty for the number of jumps,
  # the number of selected features,
  # the adjusted Rand index for the states,
  # the adjusted Rand index for the weights,
  # the number of correctly selected features,
  # the number of incorrectly selected features
  
  
  YY=true_data[[seed]]$YY
  true_states=true_data[[seed]]$true_states
  
  res=sparse_jump(Y=as.matrix(YY), n_states=as.integer(K), 
                  max_features=kappa, 
                  jump_penalty=lambda,
                  max_iter=as.integer(10), 
                  tol=1e-4, 
                  n_init=as.integer(10), 
                  verbose=F)
  
  est_weights=res[[2]]
  norm_est_weights=est_weights/sum(est_weights)
  
  est_states=order_states(res[[1]])
  true_states=order_states(true_states)
  
  Pfalse=dim(YY)[2]-P
  true_weights_seq=c(rep(1,P),rep(0,Pfalse))
  est_weights_seq=as.numeric(est_weights!=0)
  
  ARI_states=adj.rand.index(true_states,est_states)
  ARI_weights=adj.rand.index(true_weights_seq,est_weights_seq)
  
  corr_sel_feat=sum(est_weights[1:P]!=0)
  wrong_sel_feat=sum(est_weights[-(1:P)]!=0)
  
  indx=which( est_weights!=0)
  XX=YY[,indx]
  ### First version: compute BCSS as L1 norm
  Ln1=sum(get_BCSS(as.matrix(XX),est_states))
  ### Second version: compute BCSS as L2 norm
  Ln2=sqrt(sum(get_BCSS(as.matrix(XX),est_states)^2))
  
  pen=sum(est_states[1:(N-1)]!=est_states[2:N])
  alphak=length(which(est_weights!=0))
  
  # TotalPenalty=K0*alpha0+alpha0*(K-K0)+K0*(alphak-alpha0)+pen
  
  return(list(
    #YY=YY,
    indx=indx,
    est_weights=est_weights,
    norm_est_weights=norm_est_weights,
    est_weights_seq=est_weights_seq,
    est_states=est_states,
    true_states=true_states,
    seed=seed,
    N=N,
    P=P,
    m=m,
    cors=cors,
    pers=pers,
    K=K,
    lambda=lambda,
    kappa=kappa,
    Ln1=Ln1,
    Ln2=Ln2,
    # TotalPenalty=TotalPenalty,
    pen=pen,
    alphak=alphak,
    ARI_states=ARI_states,
    ARI_weights=ARI_weights,
    corr_sel_feat=corr_sel_feat,
    wrong_sel_feat=wrong_sel_feat))
  
}


satmod_est=function(seed,true_data,Ksat=6){
  
  # Function to estimate the parameters of a saturated model
  # seed is the seed for the random number generator
  # true_data is obtained using the function "sim_data"
  # Ksat is the number of states for the saturated model
  
  # Output:
  # A list with the L1 and L2 norms of the BCSS, 
  # the estimated weights, 
  # the normalized estimated weights,
  # and the estimated states
  
  YY=true_data[[seed]]$YY
  true_states=true_data[[seed]]$true_states
  
  res_sat=sparse_jump(Y=as.matrix(YY), n_states=as.integer(Ksat), 
                      max_features=sqrt(dim(YY)[2]), 
                      jump_penalty=0,
                      max_iter=as.integer(10), 
                      tol=1e-4, 
                      n_init=as.integer(10), 
                      verbose=F)
  
  est_weights_sat=res_sat[[2]]
  norm_est_weights_sat=est_weights_sat/sum(est_weights_sat)
  
  est_states_sat=order_states(res_sat[[1]])
  
  ### First version: compute BCSS as L1 norm
  Lnsat1=sum(get_BCSS(as.matrix(YY),est_states_sat))
  ### Second version: compute weighted BCSS (some features will disappear)
  Lnsat2=sqrt(sum(get_BCSS(as.matrix(YY),est_states_sat)^2))
  
  return(list(Lnsat1=Lnsat1,
              Lnsat2=Lnsat2,
              est_weights_sat=est_weights_sat,
              norm_est_weights_sat=norm_est_weights_sat,
              est_states_sat=est_states_sat))
  
}

GIC=function(simest,satmod,Ksat=6,alpha0,K0){
  
  # Function to compute the GICs for the sparse jump model
  # simest is a list of the results of the sparse jump model estimated for varying parameters (lambda, kappa, K)
  # satmod is a list of the results of the saturated model estimated for varying parameters (lambda, kappa, K)
  # Ksat is the number of states for the saturated model
  # alpha0 is prior number of features  
  # K0 is prior number of states
  
  # Output:
  # A data frame with the average and standard deviation of the GICs for the sparse jump model
  
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  
  K=unlist(lapply(simest,function(x)x$K))
  kappa=unlist(lapply(simest,function(x)x$kappa))
  lambda=unlist(lapply(simest,function(x)x$lambda))
  #YY=unlist(lapply(simest,function(x)x$K))
  PP=unlist(lapply(simest,function(x)x$P))
  PP=PP*(1+unlist(lapply(simest,function(x)x$m)))
  N=unlist(lapply(simest,function(x)x$N))
  seed=unlist(lapply(simest,function(x)x$seed))
  ARI.weights=unlist(lapply(simest,function(x)x$ARI_weights))
  ARI.states=unlist(lapply(simest,function(x)x$ARI_states))
  pen=unlist(lapply(simest,function(x)x$pen))
  
  # if(drop_pterm){
  #   CKp=-log(K)-log(2*pi)*1/2
  #   CKp_sat=-log(Ksat)-log(2*pi)*1/2
  # }
  # else{
  CKp=-log(K)-log(2*pi)*kappa/2
  CKp_sat=-log(Ksat)-log(2*pi)*sqrt(PP)/2
  #}
  
  anFTIC=log(log(N))*log(PP)
  anAIC=rep(2,length(N))
  anBIC=log(N)
  
  alphak=unlist(lapply(simest,function(x)x$alphak))
  pers=unlist(lapply(simest,function(x)x$pers))
  pen0=(1-pers)*N*(K0-1)
  
  TotalPenalty=(alpha0+pen0)*K+K0*(alphak-alpha0+pen-pen0) 
  
  Klk=length(unique(K))*length(unique(kappa))*length(unique(lambda))
  setup=rep(1:Klk,
            each=length(unique(seed)))
  
  Ln=unlist(lapply(simest,function(x)x$Ln1))
  Lnsat=unlist(lapply(satmod,function(x)x$Lnsat1))
  Lnsat=(rep(Lnsat,Klk))
  Ln_diff=Lnsat-Ln
  CKp_diff=CKp_sat-CKp
  
  FTIC=2*CKp_diff+(Ln_diff+anFTIC*TotalPenalty)/N
  BIC=2*CKp_diff+(Ln_diff+anBIC*TotalPenalty)/N
  AIC=2*CKp_diff+(Ln_diff+anAIC*TotalPenalty)/N
  
  corr_sel_feat=unlist(lapply(simest,function(x)x$corr_sel_feat))
  wrong_sel_feat=unlist(lapply(simest,function(x)x$wrong_sel_feat))
  est_weights=data.frame(matrix(unlist(lapply(simest,function(x)x$norm_est_weights)),
                                ncol = PP[1],byrow=T))
  
  res=data.frame(seed,
                 FTIC,
                 BIC,
                 AIC,
                 Ln_diff,
                 CKp_diff,
                 pen,
                 lambda,
                 kappa,
                 K,
                 ARI.states,
                 ARI.weights,
                 setup,
                 corr_sel_feat,
                 wrong_sel_feat)
  
  res_est_weights=data.frame(setup,
                             est_weights)
  
  res=res%>%group_by(setup)%>%
    summarise(avFTIC=mean(FTIC),
              avAIC=mean(AIC),
              avBIC=mean(BIC),
              sdFTIC=sd(FTIC),
              sdAIC=sd(AIC),
              sdBIC=sd(BIC),
              avARI_states=mean(ARI.states),
              avARI_weights=mean(ARI.weights),
              avLn_diff=mean(Ln_diff),
              avCKp_diff=mean(CKp_diff),
              avPen=mean(pen),
              lambda=mean(lambda),
              kappa=mean(kappa),
              K=mean(K),
              corr_sel_feat=mean(corr_sel_feat),
              wrong_sel_feat=mean(wrong_sel_feat))
  
  res_est_weights=res_est_weights%>%group_by(setup)%>%
    summarise_all(mean)
  
  res_av=data.frame(res,res_est_weights[,-1])
  
  return(res_av)
}


# GLM-JM ------------------------------------------------------------------

MClinregSim=function(n,
                   coeff,
                   X,
                   pers=.95, 
                   init,
                   seed=123,
                   family,
                   sigma=0.5){
  # This function simulates data from a linear regression model with Markov chain underlying dynamics

  # Arguments:
  # n: number of observations
  # coeff: matrix of coefficients for the linear regression model
  # X: matrix of covariates
  # pers: self-transition probability
  # init: initial probabilities
  # seed: seed for the random number generator

  # Value:
  # SimData: vector of simulated responses
  
  X=cbind(1,X)
  
  #markov chain simulation
  reg=nrow(coeff)
  Q <- matrix(rep((1-pers)/(reg-1),reg*reg), 
              ncol = reg,
              byrow = TRUE)
  diag(Q)=rep(pers,reg)
  #reg = dim(Q)[1]
  x <- numeric(n)
  set.seed(seed)
  x[1] <- sample(1:reg, 1, prob = init)
  for(i in 2:n){
    x[i] <- sample(1:reg, 1, prob = Q[x[i - 1], ])
  }
  
  d=1
  Sim = matrix(0, n, d * reg)
  SimData = matrix(0, n, d)
  
  #pseudo-observations simulation
  # for (k in 1:reg) {
  #   #u = rCopula(n, copula::tCopula(param=P2p(R[,,k]), dim = d,df=nu[k],dispstr = "un"))
  #   u = X%*%coeff[k,]
  #   Sim[, (d * k - d + 1):(d * k)] = u
  # }
  
  switch (family,
    gaussian={
      for (k in 1:reg) {
        #u = rCopula(n, copula::tCopula(param=P2p(R[,,k]), dim = d,df=nu[k],dispstr = "un"))
        u = X%*%coeff[k,]+ rnorm(n = nrow(X), mean = 0, sd = sigma)
        Sim[, (d * k - d + 1):(d * k)] = u
      }
    },
    binomial={
      for (k in 1:reg) {
        # Compute the linear predictor
        linear_predictor <- X %*% coeff[k,]
        
        # Compute the probability using the logistic function
        probabilities <- 1 / (1 + exp(-linear_predictor))
        
        # Simulate the binomial responses
        u <- rbinom(n = nrow(X), size = 1, prob = probabilities)
        Sim[, (d * k - d + 1):(d * k)] = u
      }
      },
    poisson={
      for (k in 1:reg) {
      # Compute the linear predictor
      linear_predictor <- X %*% coeff[k,]
      
      # Compute the rate (lambda) using the exponential function
      lambda <- exp(linear_predictor)
      
      # Simulate the Poisson responses
      u <- rpois(n = nrow(X), lambda = lambda)
      Sim[, (d * k - d + 1):(d * k)] = u
      }
    }
  )
  
  for (i in 1:n) {
    k = x[i]
    SimData[i, ] = Sim[i, (d * k - d + 1):(d * k)]
  }
  
  return(list(SimData=SimData,mc=x))
  
}

jumpLR <- function(Y,X, n_states, jump_penalty=1e-5, initial_states=NULL,
                   max_iter=10, n_init=10, tol=NULL, verbose=FALSE,family="gaussian") {
  # Fit jump model using framework of Bemporad et al. (2018)
  
  # state initialization, can be provided as input or computed via init_states
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s <- init_states(Y, n_states)+1
  }
  
  n_obs <- nrow(Y)
  n_features <- ncol(X)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=1)
    #mu <- matrix(0, nrow=n_states, ncol=n_obs)
    
    cfs=matrix(0,nrow=n_states,ncol=n_features+1)
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        # Fit model by updating mean of observed states
        #if(sum(s==i)>1){
        mod <- glm(Y[s==i,] ~ X[s==i,],
                   family = family)
        #mu[i,] <- colMeans(Y[s==i,])
        temp=coefficients(mod)
        cfs[i,]=as.vector(temp)
        
        
        switch(family,
               gaussian={
                 mu[i,]=mean(cbind(1,X[s==i,])%*%cfs[i,])
               },
               poisson={
                 mu[i,]=mean(exp(cbind(1,X[s==i,])%*%cfs[i,]))
               },
               binomial={
                 tmp=cbind(1,X[s==i,])%*%cfs[i,]
                 #map mu to [0,1]
                 mu[i,]=mean(1/(1+exp(-tmp)))
               }
        )
        #}
        # else{
        #   mu[i,]=mean(Y[s==i,])
        # }
      }
      
      # Order states here?
      s=order_states(s)
      
      # Fit state sequence
      s_old <- s
      
      loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      for(k in 1:n_states){
        loss_by_state[,k]=
          apply(Y,1,function(x)dist(rbind(x,mu[k,]),method="euclidean"))^2
        
      }
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
      }
      
      s[1] <- which.min(V[1,])
      for (t in 2:n_obs) {
        s[t] <- which.min(V[t,] + Gamma[s[t-1],])
      }
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- min(V[1,])
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
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
    }
    if(is.null(initial_states)){
      s <- init_states(Y, n_states)+1
    }
    else{
      s <- initial_states
    }
  }
  
  return(list(s=best_s,coefs=cfs))
}


# Mixed JM with missings --------------------------------------------------

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]

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

jump_mixed <- function(Y, n_states, jump_penalty=1e-5, 
                       initial_states=NULL,
                       max_iter=10, n_init=10, tol=NULL, verbose=FALSE
                       # , 
                       # method="gower"
) {
  # Fit jump model using framework of Bemporad et al. (2018)
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # n_states: number of states
  # jump_penalty: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  
  # Value:
  # List with state sequence and imputed data
  
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Which vars are categorical and which are numeric
  
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  Ycont=Y[,-cat.indx]
  Ycat=Y[,cat.indx]
  
  n_levs=apply(Ycat, 2, function(x)length(unique(x)))
  # n_levs=apply(Ycat, 2, function(x)levels(x))
  
  
  n_cat=length(cat.indx)
  n_cont=n_features-n_cat
  
  # Initialize mu 
  mu <- colMeans(Ycont,na.rm = T)
  
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  
  #M=ifelse(is.na(Y),T,F)
  
  
  Ytil=Y
  # Impute missing values with mean of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    Ycat[,i]=ifelse(Mcat[,i],mo[i],Ycat[,i])
    Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  }
  
  # Ycat=Ycat%>%mutate_all(factor)
  # for(i in 1:n_cat){
  #   Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  # }
  
  Y[,-cat.indx]=Ycont
  Y[,cat.indx]=Ycat
  
  # State initialization through kmeans++
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s=initialize_states(Y,n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features-length(cat.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        # Fit model by updating mean of observed states
        #if(sum(s==i)>1){
        mu[i,] <- colMeans(Ycont[s==i,])
        mo[i,]=apply(Ycat[s==i,],2,Mode)
        # }
        # else{
        #   mu[i,]=mean(Y[s==i,])
        # }
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        mo[,i]=factor(mo[,i],levels=1:n_levs[i])
      }
      
      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:ncol(Ycont)){
        Ycont[,i]=ifelse(Mcont[,i],mu[s,i],Ycont[,i])
      }
      for(i in 1:ncol(Ycat)){
        Ycat[,i]=ifelse(Mcat[,i],mo[s,i],Ycat[,i])
        Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
      }
      
      Y[,-cat.indx]=Ycont
      Y[,cat.indx]=Ycat
      
      mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      
      # loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      # for(k in 1:n_states){
      #   loss_by_state[,k]=gower.dist(Y,mumo[k,])
      # }
      
      
      # var.weights in gower.dist allows for weighted distance
      
      loss_by_state=gower.dist(Y,mumo)
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
      }
      
      s[1] <- which.min(V[1,])
      for (t in 2:n_obs) {
        s[t] <- which.min(V[t,] + Gamma[s[t-1],])
      }
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- min(V[1,])
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
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
    }
    #s <- init_states(Y, n_states)+1
    s=initialize_states(Y,n_states)
  }
  
  return(list(best_s=best_s,
              Y=Y,
              Y.orig=Ytil,
              condMM=mumo))
}

jump_mixed2 <- function(Y, n_states, jump_penalty=1e-5, 
                  initial_states=NULL,
                  max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
                  timeflag=T
                  ) {
  # Updated version of function jump_mixed (work in progress)
  
  # Fit jump model using framework of Bemporad et al. (2018)
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors with numerical levels.
  # n_states: number of states
  # jump_penalty: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  
  # Value:
  # List with state sequence and imputed data
  
  Y=data.frame(Y)
  Ynoscaled=Y
  
  #########
  ## NEW ##
  #timeflag=F
  if(timeflag){
    time=Y[,1]
    dtime=diff(time)
    dtime=dtime/as.numeric(min(dtime))
    dtime=as.numeric(dtime)
    Y=Y[,-1]
    Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
    Ynoscaled=Ynoscaled[,-1]
  }
  else{
    Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  }
  #########
  
  n_states=as.integer(n_states)

  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Which vars are categorical and which are numeric
  
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  Ycont=Y[,-cat.indx]
  Ycat=Y[,cat.indx]
  
  # n_levs=apply(Ycat, 2, function(x)length(unique(x)))
  # n_levs=apply(Ycat, 2, function(x)levels(x))
  
  
  n_cat=length(cat.indx)
  n_cont=n_features-n_cat
  
  # Initialize mu 
  mu <- colMeans(Ycont,na.rm = T)
  
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  
  #M=ifelse(is.na(Y),T,F)
  
  
  # Ytil=Y
  # Impute missing values with mean of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    x=Ycat[,i]
    Ycat[which(is.na(Ycat[,i])),i]=mo[i]
    #Ycat[,i]=ifelse(Mcat[,i],mo[i],Ycat[,i])
    #Ycat[,i]=factor(Ycat[,i],levels=unique(x[!is.na(x)]))
  }
  
  # Ycat=Ycat%>%mutate_all(factor)
  # for(i in 1:n_cat){
  #   Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  # }
  
  Y[,-cat.indx]=Ycont
  Y[,cat.indx]=Ycat
  
  # State initialization through kmeans++
   if (!is.null(initial_states)) {
     s <- initial_states
   } else {
     s=initialize_states(Y,n_states)
   }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features-length(cat.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        mu[i,] <- colMeans(Ycont[s==i,])
        mo[i,]=apply(Ycat[s==i,],2,Mode)
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        #mo[,i]=factor(mo[,i],levels=1:n_levs[i])
        mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
        
      }

      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:ncol(Ycont)){
        Ycont[,i]=ifelse(Mcont[,i],mu[s,i],Ycont[,i])
      }
      for(i in 1:ncol(Ycat)){
        Ycat[Mcat[,i],i]=mo[s[Mcat[,i]],i]
        #Ycat[,i]=ifelse(Mcat[,i],mo[s,i],Ycat[,i])
        #Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
      }
      
      Y[,-cat.indx]=Ycont
      Y[,cat.indx]=Ycat
      
      mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(Y)
      # loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      # for(k in 1:n_states){
      #   loss_by_state[,k]=gower.dist(Y,mumo[k,])
      # }
      
      
      # var.weights in gower.dist allows for weighted distance
  
      loss_by_state=gower.dist(Y,mumo)
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        if(timeflag){
          V[t-1,] <- loss_by_state[t-1,] + apply(V[t,]/dtime[t] + Gamma, 2, min)
        }
        else{
          V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
        }
      }
      
      s[1] <- which.min(V[1,])
      for (t in 2:n_obs) {
        s[t] <- which.min(V[t,] + Gamma[s[t-1],])
      }
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- min(V[1,])
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
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
    }
    #s <- init_states(Y, n_states)+1
    s=initialize_states(Y,n_states)
  }
  
  ### NEW ###
  # Final imputation

  Ycont=Ynoscaled[,-cat.indx]
  Ycat=Ynoscaled[,cat.indx]
  
  mu=matrix(0,nrow=n_states,ncol=n_cont)
  mo=matrix(0,nrow=n_states,ncol=n_cat)
  
  for (i in unique(best_s)) {
    mu[i,] <- colMeans(Ycont[best_s==i,],na.rm=T)
    mo[i,] <- as.numeric(apply(Ycat[best_s==i,],2,Mode,na.rm=T))
  }
  
  mu=data.frame(mu)
  mo=data.frame(mo,stringsAsFactors=TRUE)
  for(i in 1:n_cat){
    x=Ycat[,i]
    #mo[,i]=factor(mo[,i],levels=1:n_levs[i])
    mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
  }
  
  for(i in 1:ncol(Ycont)){
    Ycont[,i]=ifelse(Mcont[,i],mu[best_s,i],Ycont[,i])
  }
  
  for(i in 1:ncol(Ycat)){
    Ycat[Mcat[,i],i]=mo[best_s[Mcat[,i]],i]
  }
  
  Y[,-cat.indx]=Ycont
  Y[,cat.indx]=Ycat
  mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
  mumo[,cat.indx]=mo
  mumo[,cont.indx]=mu
  colnames(mumo)=colnames(Y)
  ######
  
  best_s=factor(best_s)
  levels(best_s)=1:length(unique(best_s))
  best_s=as.numeric(best_s)
  
  
  return(list(best_s=best_s,
              Y=Y,
              Y.orig=Ynoscaled,
              condMM=mumo,
              K=n_states,
              lambda=jump_penalty))
}

GIC_mixed=function(Y,states,states_sat,K,K0=3,pers0=.90,timeflag=F,Ksat=6,lambda=NULL){
  
  # Returns error when unique(states) are less than K (TBD)
  
  if(is.null(K0)){
    K0=length(unique(states))
  }
  
  #K=length(unique(states))
  if(timeflag){
    Y=Y[,-1]
  }
  PP=ncol(Y)
  N=nrow(Y)
  kappa=sqrt(PP)
  alpha0=PP
  alphak=alpha0
  
  DD=get_BCD(Y,states)
  DDsat=get_BCD(Y,states_sat)
  
  CKp=-log(K)-log(2*pi)*kappa/2
  CKp_sat=-log(Ksat)-log(2*pi)*sqrt(PP)/2
  CKp_diff=CKp_sat-CKp
  
  anFTIC=log(log(N))*log(PP)
  anAIC=rep(2,length(N))
  anBIC=log(N)
  
  pen=sum(states[1:(N-1)]!=states[2:N])
  pen0=(1-pers0)*N*(K0-1)
  
  TotalPenalty=(alpha0+pen0)*K+K0*(alphak-alpha0+pen-pen0) 
  Ln=sum(DD$BCD,na.rm = T)
  Lnsat=sum(DDsat$BCD,na.rm = T)
  Ln_diff=Lnsat-Ln
  
  FTIC=2*CKp_diff+(Ln_diff+anFTIC*TotalPenalty)/N
  BIC=2*CKp_diff+(Ln_diff+anBIC*TotalPenalty)/N
  AIC=2*CKp_diff+(Ln_diff+anAIC*TotalPenalty)/N
  
  
  return(list(FTIC=FTIC,
              BIC=BIC,
              AIC=AIC))
}

sim_data_mixed_rand=function(TT,P,Ktrue,pers, seed=123){
  
  # Function to simulate mixed data with random parameters for the data generating process
  
  # Arguments:
  # TT: number of observations
  # P: number of features
  # Ktrue: number of states
  # pers: self-transition probability
  # seed: seed for the random number generator
  
  # Value:
  # SimData: matrix of simulated data
  # mchain: vector of simulated states
  
  
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
  Pcont=P/2
  SimCont = matrix(0, TT, Pcont * Ktrue)
  SimDataCont = matrix(0, TT, Pcont)
  
  Sigma <- array(0, dim =c(Pcont,Pcont,Ktrue))
  
  set.seed(seed)
  for(i in 1:Ktrue){
    Sigma[,,i] <- matrix(0, ncol = Pcont,  
                         byrow = TRUE)
    diag(Sigma[,,i])=runif(Pcont,.2,2)
  }
  
  set.seed(seed)
  Mu <- matrix(runif(Pcont*Ktrue,min=-5,max=5), 
               nrow = Pcont, 
               byrow = TRUE)
  
  for(k in 1:Ktrue){
    u = MASS::mvrnorm(TT,Mu[,k],Sigma[,,k])
    SimCont[, (Pcont * k - Pcont + 1):(k * Pcont)] = u
  }
  
  # Categorical variables simulation
  Pcat=P/2
  SimCat = matrix(0, TT, Pcat* Ktrue)
  SimDataCat = matrix(0, TT, Pcat)
  
  
  set.seed(seed)
  Np=sample(2:6,Pcat,repl=TRUE)
  u=rep(0,Pcat)
  for(p in 1:Pcat){
    B <- matrix(dirmult::rdirichlet(n=Ktrue, alpha=rep(1,Np[p])),
                nrow = Ktrue,
                byrow =F)
    for(k in 1:Ktrue){
      u=sample(1:Np[p],TT,replace=TRUE,prob=B[k,])
      SimCat[, Pcat*(k-1)+p]=u
    }
    
  }
  
  
  for (i in 1:TT) {
    k = x[i]
    SimDataCont[i, ] = SimCont[i, (Pcont * k - Pcont + 1):(Pcont * k)]
    SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  SimDataCont=as.data.frame(SimDataCont)
  library(dplyr)
  SimDataCat=as.data.frame(SimDataCat)
  SimDataCat=SimDataCat%>%mutate_all(as.factor)
  
  SimData=cbind(SimDataCont,SimDataCat)
  return(list(SimData=SimData,
              mchain=x,
              TT=TT,
              P=P,
              Ktrue=Ktrue,
              pers=pers, 
              seed=seed)
         )
  
}


get_cat=function(y,mc,mu,phi){
  # Function to simulate categorical data
  
  # Arguments:
  # y: continuous variable 
  # mc: Markov chain states
  # mu: numeric mean value
  # phi: conditional probability for the categorical outcome k in state k
  
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

sim_data_mixed=function(seed=123,
                        TT,
                        P,
                        Ktrue=3,
                        mu=1,
                        phi=.8,
                        rho=0,
                        Pcat=NULL,
                        pers=.95,
                        pNAs=0,
                        typeNA=2){
  
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
  
  mu=c(-mu,0,mu)
  
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
    u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  if(Pcat!=0){
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=x,mu=mu,phi=phi)
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



sim_data_mixed2=function(seed=123,
                        TT,
                        P,
                        Ktrue=3,
                        mu=1,
                        phi=.8,
                        rho=0,
                        Pcat=NULL,
                        pers=.95,
                        pNAs=0,
                        typeNA=2){
  
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
  
 # mu=c(-mu,0,mu)
  MU=mu
  mu=seq(-mu,mu,length.out=Ktrue)
  
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
    u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  if(Pcat!=0){
  SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=x,mu=MU,phi=phi)
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

simstud_JMmixed=function(seed,lambda,TT,P,
                         Ktrue=3,mu=1,
                         phi=.8,rho=0,
                         Pcat=NULL,pers=.95,
                         pNAs=0,typeNA=2){
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                          TT=TT,
                          P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  # Estimate
  est=jump_mixed(simDat$SimData.NA,
                 n_states=Ktrue,
                 jump_penalty = lambda,
                 verbose=F)

  est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  simDat$SimData.complete=simDat$SimData.complete%>%
    mutate_if(is.factor,factor,levels=c(1,2,3))
  
  imput.err=gower_dist(est$Y,simDat$SimData.complete)
  ARI=adj.rand.index(est$best_s,simDat$mchain)
  
  # Return
  return(list(
    imput.err=imput.err,
    ARI=ARI,
    seed=seed,
    lambda=lambda,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pers=pers,
    pNAs=pNAs,
    typeNA=typeNA,
    true_seq=simDat$mchain,
    est_seq=est$best_s,
    true_data=simDat$SimData.complete,
    est_data=est$Y))

}

simstud_JMmixed2=function(seed,lambda,TT,P,
                         Ktrue=3,mu=1,
                         phi=.8,rho=0,
                         Pcat=NULL,pers=.95,
                         pNAs=0,typeNA=2,
                         timeflag=T,
                         pGap=.2){
  
  # Additional timeflag for unequally spaced time series
  # pGap is the percentage of missing times 
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                        TT=round(TT*(1+pGap)),
                        P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  set.seed(seed)
  gaps=sort(sample(1:TT,round(TT*pGap),replace=F))
  
  time=1:(TT*(1+pGap))
  time=time[-gaps]
  time=as.Date(time,origin="2000-01-01")
  
  Y=simDat$SimData.NA[-gaps,]
  Y=data.frame(time,Y)
  simDat$SimData.complete=simDat$SimData.complete[-gaps,]
  simDat$mchain=simDat$mchain[-gaps]
  simDat$TT=TT

  # Estimate
  # est=jump_mixed(simDat$SimData.NA,
  #                n_states=Ktrue,
  #                jump_penalty = lambda,
  #                verbose=F)
  est=jump_mixed2(Y,
                  n_states=Ktrue,
                  jump_penalty = lambda,
                  verbose=F,timeflag=timeflag)
  
  est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  simDat$SimData.complete=simDat$SimData.complete%>%
    mutate_if(is.factor,factor,levels=c(1,2,3))
  
  imput.err=gower_dist(est$Y,simDat$SimData.complete)
  ARI=adj.rand.index(est$best_s,simDat$mchain)
  
  # Return
  return(list(
    imput.err=imput.err,
    ARI=ARI,
    seed=seed,
    lambda=lambda,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pers=pers,
    pNAs=pNAs,
    typeNA=typeNA,
    true_seq=simDat$mchain,
    est_seq=est$best_s,
    true_data=simDat$SimData.complete,
    est_data=est$Y))
  
}

simstud_speclust=function(seed,TT,P,
                          Ktrue=3,mu=1,
                          phi=.8,rho=0,
                          Pcat=NULL,pers=.95,
                          pNAs=0,typeNA=2){
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                        TT=TT,
                        P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  # Estimate
  # est=jump_mixed(simDat$SimData.NA,
  #                n_states=Ktrue,
  #                jump_penalty = lambda,
  #                verbose=F)
  est=mspec(simDat$SimData.NA, k = Ktrue,verbose=T)
  
  # est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  # simDat$SimData.complete=simDat$SimData.complete%>%
  #   mutate_if(is.factor,factor,levels=c(1,2,3))
  # 
  # imput.err=gower_dist(est$Y,simDat$SimData.complete)
  ARI=adj.rand.index(est$cluster,simDat$mchain)
  
  # Return
  return(list(
    #imput.err=imput.err,
    ARI=ARI,
    seed=seed,
    lambda=lambda,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pers=pers,
    pNAs=pNAs,
    typeNA=typeNA,
    true_seq=simDat$mchain,
    est_seq=est$cluster
    # ,
    # true_data=simDat$SimData.complete,
    # est_data=est$Y
    ))
  
}

simstud_speclust2=function(seed,TT,P,
                          Ktrue=3,mu=1,
                          phi=.8,rho=0,
                          Pcat=NULL,pers=.95,
                          pNAs=0,typeNA=2,pGap=.2){
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                        TT=round(TT*(1+pGap)),
                        P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  
  # Gaps
  set.seed(seed)
  gaps=sort(sample(1:TT,round(TT*pGap),replace=F))
  
  time=1:(TT*(1+pGap))
  time=time[-gaps]
  time=as.Date(time,origin="2000-01-01")
  
  Y=simDat$SimData.NA[-gaps,]
  Y=data.frame(time,Y)
  simDat$SimData.complete=simDat$SimData.complete[-gaps,]
  simDat$mchain=simDat$mchain[-gaps]
  simDat$TT=TT
  
  # Estimate
  # est=jump_mixed(simDat$SimData.NA,
  #                n_states=Ktrue,
  #                jump_penalty = lambda,
  #                verbose=F)
  est=mspec(simDat$SimData.complete, k = Ktrue,verbose=T)
  
  # est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  # simDat$SimData.complete=simDat$SimData.complete%>%
  #   mutate_if(is.factor,factor,levels=c(1,2,3))
  # 
  # imput.err=gower_dist(est$Y,simDat$SimData.complete)
  ARI=adj.rand.index(est$cluster,simDat$mchain)
  
  # Return
  return(list(
    #imput.err=imput.err,
    ARI=ARI,
    seed=seed,
    lambda=lambda,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pers=pers,
    pNAs=pNAs,
    typeNA=typeNA,
    true_seq=simDat$mchain,
    est_seq=est$cluster
    # ,
    # true_data=simDat$SimData.complete,
    # est_data=est$Y
  ))
  
}

simstud_kNN=function(seed,
                     TT,P,
                          Ktrue=3,mu=1,
                          phi=.8,rho=0,
                          Pcat=NULL,pers=.95,
                          pNAs=0,typeNA=2){
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                        TT=TT,
                        P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  # Estimate
  # est=jump_mixed(simDat$SimData.NA,
  #                n_states=Ktrue,
  #                jump_penalty = lambda,
  #                verbose=F)
  
  
  # Transform factors to ordinal
  simDat$SimData.NA=simDat$SimData.NA%>%mutate_if(is.factor,factor,levels=c(1,2,3),ordered=T)
  
  est=KNNimp(simDat$SimData.NA,meth="median")
  
  Yimp=data.frame(est)
  Yimp$Y=Yimp$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  
  # est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  simDat$SimData.complete=simDat$SimData.complete%>%
    mutate_if(is.factor,factor,levels=c(1,2,3))
  # 
  imput.err=gower_dist(Yimp,simDat$SimData.complete)
  #ARI=adj.rand.index(est$cluster,simDat$mchain)
  
  # Return
  return(list(
    imput.err=imput.err,
    #ARI=ARI,
    seed=seed,
    #lambda=lambda,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pers=pers,
    pNAs=pNAs,
    typeNA=typeNA,
    # true_seq=simDat$mchain,
    # est_seq=est$cluster
    #,
    true_data=simDat$SimData.complete,
    est_data=Yimp
  ))
  
}

simstud_missForest=function(seed,
                     TT,P,
                     Ktrue=3,mu=1,
                     phi=.8,rho=0,
                     Pcat=NULL,pers=.95,
                     pNAs=0,typeNA=2){
  # Simulate
  simDat=sim_data_mixed(seed=seed,
                        TT=TT,
                        P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  # Estimate
  # est=jump_mixed(simDat$SimData.NA,
  #                n_states=Ktrue,
  #                jump_penalty = lambda,
  #                verbose=F)
  est=missForest(simDat$SimData.NA, verbose = F)
  
  Yimp=est$ximp
  Yimp=Yimp%>%
    mutate_if(is.factor,factor,levels=c(1,2,3))
  # est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  simDat$SimData.complete=simDat$SimData.complete%>%
    mutate_if(is.factor,factor,levels=c(1,2,3))
  # 
  imput.err=gower_dist(Yimp,simDat$SimData.complete)
  #ARI=adj.rand.index(est$cluster,simDat$mchain)
  
  # Return
  return(list(
    imput.err=imput.err,
    #ARI=ARI,
    seed=seed,
    #lambda=lambda,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pers=pers,
    pNAs=pNAs,
    typeNA=typeNA,
    # true_seq=simDat$mchain,
    # est_seq=est$cluster
    #,
    true_data=simDat$SimData.complete,
    est_data=Yimp
  ))
  
}


# Continuous SJM with missings -------------------------------------------------------

initialize_states_jumpR <- function(Y, K,method="euclidean") {
  n <- nrow(Y)
  
  ### Repeat the following few times?
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(daisy(Y, metric = method))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  
  #dist_matrix <- gower.dist(Y, centroids)
  dist_matrix <- proxy::dist(Y,centroids,method=method)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}

jumpR <- function(Y, n_states, jump_penalty=1e-5, 
                  initial_states=NULL,
                  max_iter=10, n_init=10, tol=NULL, verbose=FALSE, method="euclidean",python=F) {
  # Fit jump model using framework of Bemporad et al. (2018)
  
  Y=as.matrix(Y)
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Initialize mu
  mu <- colMeans(Y,na.rm = T)
  
  # Track missings with 0 1 matrix
  M=ifelse(is.na(Y),T,F)
  
  
  Ytil=Y
  # Impute missing values with mean of observed states
  for(i in 1:n_features){
    Y[,i]=ifelse(M[,i],mu[i],Y[,i])
  }
  
  # State initialization through kmeans++
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    if(python){
      s <- init_states(Y, n_states)+1
    }
    else{
      s <- initialize_states_jumpR(Y, n_states,method)
    }
  }
  
  
  # Fill-in unobserved entries
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features)
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        # Fit model by updating mean of observed states
        if(sum(s==i)>1){
          mu[i,] <- colMeans(Y[s==i,])
        }
        else{
          mu[i,]=mean(Y[s==i,])
        }
      }
      
      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:n_features){
        # DA RIVEDERE mu[s,i], non sono sicuro di quell's forse ci va s==i o simili
        Y[,i]=ifelse(M[,i],mu[s,i],Y[,i])
      }
      
      loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      switch(method,
             euclidean={
               loss_by_state <- proxy::dist(Y,mu,method="euclidean")^2
               # for(k in 1:n_states){
               #   loss_by_state[,k]=apply(Y,1,function(x) 
               #     dist(rbind(x,mu[k,]),method="euclidean"))^2
               # }
             },
             manhattan={
               loss_by_state <- proxy::dist(Y,mu,method="manhattan")^2
               # for(k in 1:n_states){
               #   loss_by_state[,k]=apply(Y,1,function(x) 
               #     dist(rbind(x,mu[k,]),method="manhattan"))
               # }
             })
      
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
      }
      
      s[1] <- which.min(V[1,])
      for (t in 2:n_obs) {
        s[t] <- which.min(V[t,] + Gamma[s[t-1],])
      }
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- min(V[1,])
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
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
    }
    if(python){
      s <- init_states(Y, n_states)+1
    }
    else{
      s <- initialize_states_jumpR(Y, n_states,method)
    }
  }
  
  return(best_s)
}

sparse_jumpR <- function(Y, n_states, max_features, jump_penalty=1e-5,
                         max_iter=10, tol=1e-4, n_init=10, verbose=FALSE,method="euclidean",python=F) {
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  max_features <- pmax(1, pmin(max_features, sqrt(n_features)))
  feat_w <- rep(1 / sqrt(n_features), n_features)
  states <- NULL
  
  for (it in 1:max_iter) {
    states <- jumpR(Y * sqrt(feat_w), n_states, initial_states=states,
                    jump_penalty=jump_penalty, n_init=n_init,method=method,python=python)
    if (length(unique(states)) == 1) {
      break
    } else {
      new_w <- get_weights(Y, states, max_features, n_states)
    }
    if (sum(abs(new_w - feat_w)) / sum(abs(feat_w)) < tol) {
      break
    } else if (verbose) {
      cat('Iteration', it, ', w diff', sum(abs(new_w - feat_w)), '\n')
    }
    feat_w <- as.numeric(new_w)
  }
  
  return(list(states=states, feat_w=feat_w))
}


# Spatial JM --------------------------------------------------------------

Cmatrix=function(sp_indx){
  
  # This function creates the adjacency matrix for a spatial grid
  
  # Arguments:
  # sp_indx is a matrix with the spatial indexes of the grid
  
  # Value:
  # C is the adjacency matrix
  
  nc=ncol(sp_indx)
  nr=nrow(sp_indx)
  #indx=as.vector(t(sp_indx))
  M=nc*nr
  #indx=indx[-1]
  
  C=matrix(0,nrow=M,ncol=M)
  
  for(i in 1:M){
    # Check on top
    if(i/nc>1){
      C[i,(i%%nc+nc*(i%/%nc-1))]=1
    }
    #Check right
    if(i%%nc!=0){
      C[i,i+1]=1
    }
    
  }
  C=C+t(C)
  return(C)
}

sim_spatial_JM=function(P,C,seed,
                        #pers_fact=4,
                        rho=0,Pcat=NULL, phi=.8,
                        mu=3,pNAs=0){
  
  # This function simulates data from a 3-states spatial jump model (to be updated)
  
  # Arguments:
  # P: number of features
  # C: adjacency matrix of dimension MxM, where M is the desired number of spatial points
  # seed: seed for the random number generator
  # rho: correlation for the variables
  # Pcat: number of categorical variables
  # phi: conditional probability for the categorical outcome k in state k
  # mu: mean value for the continuous variables
  # pNAs: percentage of missing values (random missing pattern by default)
  
  # Value:
  # A list with the following elements:
  # SimData: a data frame with the simulated data. 
  # SimData.NA: a data frame with the simulated data with missing values
  # s: a vector with the simulated states
  
  # Continuous variables are simulated from a Gaussian distribution with mean specified in mumo and unitary st. dev.
  # Categorical variables are obtained censoring the continuous variables in three intervals, the central of which has probability equal to phi
  
  # The function simulates only 3-states spatial JM (to be updated).
  
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
  n_states=3
  M=dim(C)[1]
  #s=matrix(0,ncol=ncg,nrow=nrg)
  #s=rep(0,M)
  set.seed(seed)
  
  #require("potts")
  ncolor = as.integer(n_states) # transizione di fase continua per 1 <= ncolor <= 4
  nr = sqrt(M)
  nc = sqrt(M)
  init <- matrix(sample(ncolor, nr * nc, replace = TRUE), nrow = nr, ncol=nc)
  init <- packPotts(init, ncol = ncolor)
  
  beta <- log(1 + sqrt(ncolor))
  theta <- c(rep(0, ncolor), beta)
  out <- potts(init, param=theta, nbatch = 200 , blen=10, nspac=1)
  rotate <- function(x) apply(t(x), 2, rev)
  # Recover decoded matrix
  mat=unpackPotts(out$final)
  mat=rotate(mat)
  mat
  s=c(t(mat))
  
  # s[1]=sample(1:n_states,1)
  # #eff_it=1
  # for(m in 2:M){
  #   if(prod(s)==0){
  #     #n_prox=length(which(C[m,]==1))
  #     #probs=rep(1,n_states)/n_states
  #     #probs=probs+table(factor(s[which(C[m,]==1)],levels=1:n_states))+1/(pers_fact+.01)
  #     #probs=probs/sum(probs)
  #     #s[which(C[m,]==1)]=sample(1:n_states,n_prox,prob=probs,replace=T) 
  #     
  #     succ=table(factor(s[which(C[m,]==1)],levels=1:n_states))
  #     
  #     s[m]=sample(1:n_states,1,
  #                 prob=colMeans(MCMCprecision::rdirichlet(1000,rep(pers_fact,n_states)+succ)) ,
  #                 replace=T) 
  #     #eff_it=eff_it+1
  #   }
  #   else{
  #     break
  #   }
  #}
  
  # Continuous variables simulation
  mu=c(-mu,0,mu)
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
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=s,mu=mu,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  if(pNAs>0){
    SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=0)
    SimData.NA=as.data.frame(SimData.NA)
    SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
    SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  }
  else{
    SimData.NA=SimData
  }
  s=order_states_freq(s)
  
  return(list(SimData=SimData,SimData.NA=SimData.NA,s=s))
}

spatial_jump <- function(Y,C, n_states, jump_penalty=1e-5, 
                         initial_states=NULL,
                         max_iter=10, n_init=10, tol=NULL, verbose=FALSE
                         # , 
                         # method="gower"
){
  
  # This function implements the spatial jump algorithm for clustering spatial data
  
  # Arguments:
  # Y is a data frame with P columns (features) and M rows (spatial points)
  # C is a  MxM adjacency matrix
  # n_states is the number of states
  # jump_penalty is the penalty for jumping between states
  # initial_states is a vector of length M with the initial state of each point
  # max_iter is the maximum number of iterations
  # n_init is the number of initializations
  # tol is the tolerance for stopping the algorithm
  # verbose is a boolean for printing the loss at each iteration
  
  # Value:
  # best_s is the best state sequence
  # Y is the imputed data
  # Y.orig is the original data
  # condMM is the conditional mean-mode matrix
  
  
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Which vars are categorical and which are numeric
  
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  Ycont=Y[,-cat.indx]
  Ycat=Y[,cat.indx]
  
  n_levs=apply(Ycat, 2, function(x)length(unique(x)))
  # n_levs=apply(Ycat, 2, function(x)levels(x))
  
  n_cat=length(cat.indx)
  n_cont=n_features-n_cat
  
  # Initialize mu 
  mu <- colMeans(Ycont,na.rm = T)
  
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  
  #M=ifelse(is.na(Y),T,F)
  
  
  Ytil=Y
  # Impute missing values with mean of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    Ycat[,i]=ifelse(Mcat[,i],mo[i],Ycat[,i])
    Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  }
  
  Y[,-cat.indx]=Ycont
  Y[,cat.indx]=Ycat
  
  # State initialization through kmeans++
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s=initialize_states(Y,n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features-length(cat.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        # Fit model by updating mean of observed states
        #if(sum(s==i)>1){
        mu[i,] <- colMeans(Ycont[s==i,])
        mo[i,]=apply(Ycat[s==i,],2,Mode)
        # }
        # else{
        #   mu[i,]=mean(Y[s==i,])
        # }
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        mo[,i]=factor(mo[,i],levels=1:n_levs[i])
      }
      
      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:ncol(Ycont)){
        Ycont[,i]=ifelse(Mcont[,i],mu[s,i],Ycont[,i])
      }
      for(i in 1:ncol(Ycat)){
        Ycat[,i]=ifelse(Mcat[,i],mo[s,i],Ycat[,i])
        Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
      }
      
      Y[,-cat.indx]=Ycont
      Y[,cat.indx]=Ycat
      
      mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      
      # loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      # for(k in 1:n_states){
      #   loss_by_state[,k]=gower.dist(Y,mumo[k,])
      # }
      
      
      # var.weights in gower.dist allows for weighted distance
      
      loss_by_state=gower.dist(Y,mumo)
      
      for (m in 1:n_obs) {
        loss_by_state[m,] <- loss_by_state[m,] + Gamma[s[m],]*table(factor(s[which(C[m,]==1)],
                                                                           levels=1:n_states))
      }
      
      s=apply(loss_by_state,1,which.min)
      
      # for (m in 1:n_obs) {
      #   s[m] <- which.min(loss_by_state[m,])
      # }
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- mean(apply(loss_by_state,1,min))
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
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
    }
    #s <- init_states(Y, n_states)+1
    s=initialize_states(Y,n_states)
  }
  
  return(list(best_s=best_s,
              Y=Y,
              Y.orig=Ytil,
              condMM=mumo))
}

simstud_spatialJM=function(Ktrue=3,
                           seed,
                           gamma,
                           M,
                           P,
                           #Ktrue=3,
                           mu=3,
                           phi=.8,
                           rho=0,
                           Pcat=NULL,
                           pNAs=0){
  
  # Remember to assume squared areas (M=4,100,900,...)
  
  sp_indx=1:M
  sp_indx=matrix(sp_indx,ncol=sqrt(M),byrow=T)
  
  C=Cmatrix(sp_indx)
  simDat=sim_spatial_JM(P,C,seed,
                        rho=rho,Pcat=Pcat, phi=phi,
                        mu=mu,pNAs=pNAs)
  
  Y=simDat$SimData.NA
  # Estimate
  est=spatial_jump(Y,C, n_states=Ktrue, 
                   jump_penalty=gamma, 
                   initial_states=NULL,
                   max_iter=10, n_init=10, tol=NULL, 
                   verbose=F)
  
  est$Y=est$Y%>%mutate_if(is.factor,factor,levels=c(1,2,3))
  simDat$SimData=simDat$SimData%>%
    mutate_if(is.factor,factor,levels=c(1,2,3))
  
  imput.err=gower_dist(est$Y,simDat$SimData)
  ARI=adj.rand.index(est$best_s,simDat$s)
  
  # Return
  return(list(
    imput.err=imput.err,
    ARI=ARI,
    seed=seed,
    gamma=gamma,
    M=M,
    P=P,
    Ktrue=Ktrue,
    mu=mu,
    phi=phi,
    rho=rho,
    Pcat=Pcat,
    pNAs=pNAs,
    true_seq=simDat$s,
    est_seq=est$best_s,
    true_data=simDat$SimData,
    est_data=est$Y))
  
}


# Spatio-temporal JM --------------------------------------------------------

STjumpR=function(Y,n_states,C,jump_penalty=1e-5,
                 spatial_penalty=1e-5,
                 initial_states=NULL,
                  max_iter=10, n_init=10, tol=NULL, verbose=FALSE){
  # Y is a dataframe in long format with mandatory columns t and m which are time and spatial indexes
  # n_states is the number of states
  
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  library(dplyr)
  Y <- Y %>% select(t, m, everything())
  YY=subset(Y,select=-c(t,m))
  
  TT=length(unique(Y$t))
  M=length(unique(Y$m))
  
  cat.indx=which(sapply(YY, is.factor))
  cont.indx=which(sapply(YY, is.numeric))
  
  Ycont=YY[,cont.indx]
  Ycat=YY[,cat.indx]
  
  n_cat=length(cat.indx)
  n_cont=P-n_cat
  
  # Missing data imputation TBD
  
  # State initialization through kmeans++
  S=matrix(0,nrow=TT,ncol=M)
  for(m in 1:M){
    S[,m]=initialize_states(Y[which(Y$m==m),],n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(as.vector(S))) {
        mu[i,] <- colMeans(Ycont[as.vector(t(S))==i,])
        mo[i,]=apply(Ycat[as.vector(t(S))==i,],2,Mode)
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        #mo[,i]=factor(mo[,i],levels=1:n_levs[i])
        mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
        
      }
      mumo=data.frame(matrix(0,nrow=n_states,ncol=P))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(YY)
      
      # Fit state sequence
      S_old <- S
      
      loss_v=NULL
      
      # Re-fill-in missings TBD 
      
      for(m in 1:M){
        loss_by_state=gower.dist(Y[which(Y$m==m),-(1:2)],mumo)
        #loss_by_state_old=t(apply(loss_by_state,1,function(x){x/sum(x)}))

        for(k in 1:n_states){
           loss_by_state[,k]=loss_by_state[,k]-spatial_penalty*rowSums(S[,which(C[m,]==1)]==k)/length(which(C[m,]==1))
          #loss_by_state=loss_by_state+spatial_penalty*length(which(C[m,]==1))/(rowSums(S[,which(C[m,]==1)]==k)+1)
        }
        
        # Normalize
        # loss_by_state=t(apply(loss_by_state,1,function(x){x/sum(x)}))
        
        V <- loss_by_state
        for (t in (TT-1):1) {
          V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
          
          # Normalize
          # V[t-1,] <- V[t-1,]/sum(V[t-1,])
          
        }
        
        S[1,m] <- which.min(V[1,])
        for (t in 2:TT) {
          S[t,m] <- which.min(V[t,] + Gamma[S[t-1,m],])
        }
        loss_v=c(loss_v,min(V[1,]))
      }

      
      
      if (length(unique(S)) == 1) {
        break
      }
      loss <- mean(loss_v)
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
    }
    if (is.null(best_s) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_s <- S
    }
    
    for(m in 1:M){
      S[,m]=initialize_states(Y[which(Y$m==m),],n_states)
    }
  }
  
  #best_s=matrix(order_states_freq(best_s),ncol=M,byrow=T)
  
  return(list(best_s=best_s,
              Y=Y,
              K=n_states,
              lambda=jump_penalty,
              gamma=spatial_penalty))
  
}

sim_spatiotemp_JM=function(P,C,seed,
                           n_states=3,
                           rho=0,Pcat=NULL, phi=.8,
                           mu=3,pNAs=0,ST=NULL,PI=.5){
  
  # This function simulates data from a 3-states spatial jump model (to be updated)
  
  # Arguments:
  # P: number of features
  # C: adjacency matrix of dimension MxM, where M is the desired number of spatial points
  # seed: seed for the random number generator
  # rho: correlation for the variables
  # Pcat: number of categorical variables
  # phi: conditional probability for the categorical outcome k in state k
  # mu: mean value for the continuous variables
  # pNAs: percentage of missing values (random missing pattern by default)
  # ST: state sequence at time t-1
  # PI: probability of sampling the same state as the previous one
  
  # Value:
  # A list with the following elements:
  # SimData: a data frame with the simulated data. 
  # SimData.NA: a data frame with the simulated data with missing values
  # s: a vector with the simulated states
  
  # Continuous variables are simulated from a Gaussian distribution with mean specified in mumo and unitary st. dev.
  # Categorical variables are obtained censoring the continuous variables in three intervals, the central of which has probability equal to phi
  
  # The function simulates only 3-states spatial JM (to be updated).
  
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
  #n_states=3
  M=dim(C)[1]
  #s=matrix(0,ncol=ncg,nrow=nrg)
  #s=rep(0,M)
  set.seed(seed)
  
  #require("potts")
  ncolor = as.integer(n_states) # transizione di fase continua per 1 <= ncolor <= 4
  nr = sqrt(M)
  nc = sqrt(M)
  init <- matrix(sample(ncolor, nr * nc, replace = TRUE), nrow = nr, ncol=nc)
  init <- packPotts(init, ncol = ncolor)
  
  beta <- log(1 + sqrt(ncolor))
  theta <- c(rep(0, ncolor), beta)
  
  # Per dare un label sensato ai cluster occorre necessariamente simulare un tempo alla volta
  # e riordinare gli stati con order_states_condMeans.
  # Chiedere ad Antonio se nbatch e nblen cosi hanno senso
  out <- potts(init, param=theta, nbatch = 200 , blen=1, nspac=1)
  
  rotate <- function(x) apply(t(x), 2, rev)
  # Recover decoded matrix
  mat=unpackPotts(out$final)
  mat=rotate(mat)
  #mat
  s=c(t(mat))
  
  # Sample between previous state and Potts state
  if(!is.null(ST)){
    for(i in 2:M){
      s[i]=sample(c(s[i],ST[i]),1,prob=c(1-PI,PI))
    }
  }
  
  
  # Continuous variables simulation
  #mu=c(-mu,0,mu)
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
    SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=0)
    SimData.NA=as.data.frame(SimData.NA)
    SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
    SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  }
  else{
    SimData.NA=SimData
  }
  #s=order_states_freq(s)
  
  return(list(SimData=SimData,SimData.NA=SimData.NA,s=s))
}



sim_obs=function(s=s,mu=mu,rho=rho,P=P,Pcat=NULL,n_states=n_states,seed=seed,pNAs=pNAs){
  
  # This function simulates observations conditional on latent states s
  
  # Arguments:
  # s: vector of latent states
  # mu: state conditional mean
  # rho: state conditional correlation
  # P: number of variables
  # n_states: number of states
  # seed: seed
  
  mu=c(-mu,0,mu)
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
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=s,mu=mu,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  if(pNAs>0){
    SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=0)
    SimData.NA=as.data.frame(SimData.NA)
    SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
    SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  }
  else{
    SimData.NA=SimData
  }
  return(list(SimData=SimData,SimData.NA=SimData.NA))
}

simstud_STJump=function(lambda,gamma,seed,M,TT,
                        mu=1,rho=0.2,
                        K=3,P=30,phi=.8,Pcat=10,pNAs=0,PI=.9){
  
  sp_indx=1:M
  sp_indx=matrix(sp_indx,ncol=sqrt(M),byrow=T)
  S_true=matrix(0,nrow=TT,ncol=M)
  C=Cmatrix(sp_indx)
  Y=NULL
  
  t=1
  simDat=sim_spatiotemp_JM(P,C,seed=seed+seed*1000+t-1,
                           rho=rho,Pcat=Pcat, phi=phi,
                           mu=mu,pNAs=pNAs,ST=NULL,n_states=K)
  temp=data.frame(simDat$SimData)
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=simDat$s
  S_true[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S_true[t,])
  
  # Temporal persistence 
  # PI=0.7
  for(t in 2:TT){
    simDat=sim_spatiotemp_JM(P,C,seed=seed+seed*1000+t-1,
                             rho=rho,Pcat=Pcat, phi=phi,
                             mu=mu,pNAs=pNAs,ST=S_true[t-1,],PI=PI,n_states=K)
    temp=data.frame(simDat$SimData)
    temp$m=1:M
    temp$t=rep(t,M)
    Y=rbind(Y,temp)
    S_true[t,]=simDat$s
    S_true[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S_true[t,])
    
  }
  
  # Put t and m in front of all the others with dplyr
  Y <- Y %>% select(t, m, everything())
  
  st=Sys.time()
  fit <- STjumpR(Y, n_states = K, C=C, jump_penalty=lambda,spatial_penalty = gamma, verbose=F)
  end=Sys.time()
  elapsed=end-st
  
  best_s=fit$best_s
  for(t in 1:TT){
    best_s[t,]=order_states_condMean(Y$V20[Y$t==t],best_s[t,])
  }
  
  ARI=adj.rand.index(S_true,best_s)
  
  
  return(list(lambda=lambda,gamma=gamma,ARI=ARI,
              seed=seed,
              M=M,TT=TT,
              mu=mu,rho=rho,
              K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs,
              elapsed=elapsed,
              S_true=S_true,best_s=best_s))
}