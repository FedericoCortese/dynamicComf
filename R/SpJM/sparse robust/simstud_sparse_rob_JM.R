library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_sparse_robust_2.R")

zeta0=seq(0.1,0.5,0.05)
alpha=.1

# Check how to modify final evaluation to include also K>2
K=2
tol=1e-16
n_outer=10
verbose=T
lambda=seq(0,.5,.1)
nseed=50

TT=1000
P=10

c=5

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  K=K,
  lambda=lambda,
  TT=TT,
  P=P,
  c=c
)

ncores=parallel::detectCores()-1

start=Sys.time()
res_list_K2 <- mclapply(seq_len(nrow(hp)), function(i) {
  
  seed <- hp$seed[i]
  zeta0    <- hp$zeta0[i]
  K        <- hp$K[i]
  lambda   <- hp$lambda[i]
  TT       <- hp$TT[i]
  P        <- hp$P[i]
  c        <- hp$c[i]
  
  set.seed(seed)
  simDat=sim_data_stud_t(seed=123,
                         TT=TT,
                         P=P,
                         Pcat=NULL,
                         Ktrue=2,
                         mu=2,
                         rho=0,
                         nu=100,
                         phi=.8,
                         pers=0.99)
  
  Y=simDat$SimData
  true_stat=simDat$mchain
  
  nu=4
  # For State 1, only features 1,2 and 3 are relevant, the rest are noise
  indx=which(true_stat!=1)
  Sigma <- matrix(0,ncol=P-3,nrow=P-3)
  diag(Sigma)=5
  Y[indx,-(1:3)]=mvtnorm::rmvt(length(indx),
                               sigma = (nu-2)*Sigma/nu,
                               df = nu, delta = rep(0,P-3))
  
  # For State 2, only features 3,4 and 5 are relevant, the rest are noise
  indx=which(true_stat!=2)
  Y[indx,-(3:5)]=mvtnorm::rmvt(length(indx),
                               sigma = (nu-2)*Sigma/nu,
                               df = nu, delta = rep(0,P-3))
  
  
  Sigma <- matrix(0,ncol=P-5,nrow=P-5)
  diag(Sigma)=5
  
  # All other features are noise
  Y[,6:P]=mvtnorm::rmvt(TT,
                        sigma = (nu-2)*Sigma/nu,
                        df = nu, delta = rep(0,P-5))
  
  # Introduce outliers
  set.seed(123)
  out_sigma=100
  N_out=TT*0.02
  t_out=sample(1:TT,size=N_out)
  Y[t_out,]=Y[t_out,]+rnorm(N_out*P,0,out_sigma)
  
  # Set the truth for latent states sequence
  truth=simDat$mchain
  truth[t_out]=0
  
  # Set the truth for features
  W_truth=matrix(F,nrow=K,ncol=P)
  W_truth[,3]=T
  W_truth[1,1:2]=T
  W_truth[2,4:5]=T
  
  fit=robust_JM_COSA(Y=as.matrix(Y),
                     zeta0=zeta0,
                     lambda=lambda,
                     K=K,
                     tol        = 1e-16,
                     n_init     = 2,
                     n_outer    = 10,
                     alpha      = 0.1,
                     verbose    = F,
                     knn        = 10,
                     c          = c,
                     M          = NULL)
  
  est_s=fit$s
  est_s[fit$v<0.5]=0
  
  W_ind=fit$W>0.01
  
  ARI_s=mclust::adjustedRandIndex(est_s,truth)
  ARI_W=mclust::adjustedRandIndex(W_ind,W_truth)
  
  res <- list(
    seed = seed,
    zeta0 = zeta0,
    K = K,
    lambda = lambda,
    TT = TT,
    P = P,
    c = c,
    W = fit$W,
    s=fit$s,
    v=fit$v,
    truth = truth,
    ARI_s=ARI_s,
    ARI_W=ARI_W
  )
  
}, mc.cores = ncores)
end=Sys.time()

print(end-start)
