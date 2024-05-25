library(RcppHMM)
library(reticulate)
library(pdfCluster)
library(boot)
library(xtable)
library(dplyr)
#py_install("scipy")

# Import the python module
import("scipy")

# Import python functions for SJM estimation
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
                                      "Pi" = Pi))
  set.seed(seed)
  observationSequence <- generateObservations(HMM.cont.multi, N)
  Y=t(observationSequence$Y)
  Y=apply(Y,2,scale)
  true_states=order_states(as.factor(observationSequence$X))
  
  set.seed(seed)
  Ytil=Y[sample(N),]
  for(i in 2:m){
    Ytil=cbind(Ytil,Y[sample(N),])
  }
  
  YY=cbind(Y,Ytil)
  
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

jumpR <- function(Y, n_states, jump_penalty=1e-5, 
                  #initial_states=NULL,
                  max_iter=10, n_init=10, tol=NULL, verbose=FALSE, method="euclidean") {
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
    Y[,i]=ifelse(M[,i],mu,Y[,i])
  }
  
  # State initialization through kmeans++
  # if (!is.null(initial_states)) {
  #   s <- initial_states
  # } else {
  s <- init_states(Y, n_states)+1
  # }
  
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
        Y[,i]=ifelse(M[,i],mu[s,i],Y[,i])
      }
      
      loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      switch(method,
             euclidean={
               for(k in 1:n_states){
                 loss_by_state[,k]=apply(Y,1,function(x) dist(rbind(x,mu[k,]),method="euclidean"))^2
               }
             },
             manhattan={
               for(k in 1:n_states){
                 loss_by_state[,k]=apply(Y,1,function(x) dist(rbind(x,mu[k,]),method="manhattan"))
               }
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
    s <- init_states(Y, n_states)+1
  }
  
  return(best_s)
}

sparse_jumpR <- function(Y, n_states, max_features, jump_penalty=1e-5,
                        max_iter=10, tol=1e-4, n_init=10, verbose=FALSE) {
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  max_features <- pmax(1, pmin(max_features, sqrt(n_features)))
  feat_w <- rep(1 / sqrt(n_features), n_features)
  states <- NULL

for (it in 1:max_iter) {
  states <- jumpR(Y * sqrt(feat_w), n_states, initial_states=states,
                  jump_penalty=jump_penalty, n_init=n_init)
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
  feat_w <- new_w
}

return(list(states=states, feat_w=feat_w))
}

jump_mixed <- function(Y, n_states, jump_penalty=1e-5, 
                  #initial_states=NULL,
                  max_iter=10, n_init=10, tol=NULL, verbose=FALSE, method="euclidean") {
  # Fit jump model using framework of Bemporad et al. (2018)
  
  Y=as.matrix(Y)
  n_states=as.integer(n_states)

  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # P-dim vector identifying continuous variables
  # cont_vars <- rep(1, n_features)
  
  
  # Initialize mu
  mu <- colMeans(Y,na.rm = T)
  
  # Track missings with 0 1 matrix
  M=ifelse(is.na(Y),T,F)
  
  
  Ytil=Y
  # Impute missing values with mean of observed states
  for(i in 1:n_features){
    Y[,i]=ifelse(M[,i],mu,Y[,i])
  }
  
  # State initialization through kmeans++
  # if (!is.null(initial_states)) {
  #   s <- initial_states
  # } else {
  s <- init_states(Y, n_states)+1
  # }
  
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
        Y[,i]=ifelse(M[,i],mu[s,i],Y[,i])
      }
      
      loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      
      for(k in 1:n_states){
        loss_by_state[,k]=apply(Y,1,function(x) dist(rbind(x,mu[k,]),method="manhattan"))
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
    s <- init_states(Y, n_states)+1
  }
  
  return(best_s)
}