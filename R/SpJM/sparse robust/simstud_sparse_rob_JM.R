library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_sparse_robust_2.R")

zeta0=seq(0.05,0.4,by=.05)
alpha=.1

# Check how to modify final evaluation to include also K>2
K=2:3
tol=1e-16
n_outer=10
verbose=F
lambda=seq(0,1,.2)
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
  simDat=sim_data_stud_t(seed=seed,
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
  set.seed(seed)
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
                     n_init     = 5,
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

#save(res_list_K2,file='simple_simstud_rob_JM_K2.Rdata')

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_robJM/simple_simstud_rob_JM_K2.Rdata")

df_results_robJM_K2 <- do.call(rbind, lapply(res_list_K2, function(x) {
  x_sel <- x[c("seed", "zeta0", "K", "lambda", "TT", "P", "c", "ARI_s", "ARI_W")]
  as.data.frame(x_sel, stringsAsFactors = FALSE)
}))

library(dplyr)

av_results_robJM_K2 <- df_results_robJM_K2 %>%
  group_by(zeta0, K, lambda, TT, P, c) %>%
  summarise(
    mean_ARI_s  = mean(ARI_s),
    sd_ARI_s    = sd(ARI_s),
    mean_ARI_W  = mean(ARI_W),
    sd_ARI_W    = sd(ARI_W),
    n           = n(),        # numero di repliche (semi)
    .groups     = "drop"
  )


library(plotly)

# Supponendo che av_results_robJM_K2 sia il tuo tibble

# 1) superfice 3D per mean_ARI_s
fig_s <- plot_ly(
  data = av_results_robJM_K2[av_results_robJM_K2$K==2,],
  x = ~lambda, y = ~zeta0, z = ~mean_ARI_s,
  type = "mesh3d",
  intensity = ~mean_ARI_s,
  colors = colorRamp(c("blue", "yellow", "red"))
) %>%
  layout(
    scene = list(
      xaxis = list(title = "λ"),
      yaxis = list(title = "ζ₀"),
      zaxis = list(title = "Mean ARI[s]")
    ),
    title = " "
  )

# 2) superfice 3D per mean_ARI_W
fig_W <- plot_ly(
  data = av_results_robJM_K2[av_results_robJM_K2$K==2,],
  x = ~lambda, y = ~zeta0, z = ~mean_ARI_W,
  type = "mesh3d",
  intensity = ~mean_ARI_W,
  colors = colorRamp(c("blue", "yellow", "red"))
) %>%
  layout(
    scene = list(
      xaxis = list(title = "λ"),
      yaxis = list(title = "ζ₀"),
      zaxis = list(title = "Mean ARI[W]")
    ),
    title = " "
  )

# Per visualizzarli uno sotto l'altro:
fig_s
fig_W




