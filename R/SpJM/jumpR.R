library(RcppHMM)
library(reticulate)
library(pdfCluster)
library(boot)
library(xtable)
#py_install("scipy")

# Import the python module
import("scipy")

# Import python functions for SJM estimation
source_python('SJ.py')


# Simulation --------------------------------------------------------------

source("Utils.R")

# Ktrue=2 -----------------------------------------------------------------

pers1=.99
Ns=c(300,600,1000)
seed=123
corsK2=c(.8,.4)
YYs_K2_300=lapply(seed,sim_data,Ktrue=2,N=Ns[1],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_600=lapply(seed,sim_data,Ktrue=2,N=Ns[2],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_1000=lapply(seed,sim_data,Ktrue=2,N=Ns[3],P=100,cors=corsK2,pers=pers1,m=2)

Y=YYs_K2_1000[[1]]$YY
dim(Y)
true.st=YYs_K2_1000[[1]]$true_states

jump_penalty=1e-5
initial_states=NULL
max_iter=10
n_init=10
tol=NULL
verbose=FALSE

n_states=2
Y=as.matrix(Y)
n_states=as.integer(n_states)

jumpR <- function(Y, n_states, jump_penalty=1e-5, initial_states=NULL,
                 max_iter=10, n_init=10, tol=NULL, verbose=FALSE) {
  # Fit jump model using framework of Bemporad et al. (2018)
  
  # state initialization, can be provided as input or computed via init_states
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s <- init_states(Y, n_states)+1
  }
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
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
      
      loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      for(k in 1:n_states){
        loss_by_state[,k]=apply(Y,1,function(x) dist(rbind(x,mu[k,]),method="euclidean"))^2
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

prv=jumpR(Y,n_states,jump_penalty = 100,verbose = F)
good=jump(Y,n_states,jump_penalty = 100,verbose = F)
plot(good,type='l')
plot(prv,type='l')


# Linreg ------------------------------------------------------------------


##### jump LR
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
          
          #Not sure of the following
          mu[i,]=mean(cbind(1,X[s==i,])%*%cfs[i,])
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
        loss_by_state[,k]=apply(Y,1,function(x) dist(rbind(x,mu[k,]),method="euclidean"))^2
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
  
  return(list(s=best_s,coefs=cfs))
}


TT=1000
X=matrix(0,ncol=2,nrow=TT)
X[,1]=rnorm(TT,sd=2)
X[,2]=rnorm(TT,mean=2)
X=apply(X,2,scale)
coeff=matrix(0,nrow=2,ncol=3)
coeff[1,]=c(0,2,1)
coeff[2,]=c(2,-2,.5)
pers=.95
init=c(.8,.2)
Y=MClinregSim(TT,coeff,X,pers,init)
mchain=Y$mc
Y=Y$SimData
plot(X[,1],Y)
plot(X[,2],Y)
plot(mchain,type='l')

Y=as.matrix(Y)
n_states=as.integer(2)

prvLR=jumpLR(Y,X,n_states,jump_penalty = 10,verbose = F)
round(prvLR$coefs,3)
coeff
table(prvLR$s)
table(mchain)
