library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_sparse_robust_2.R")

TT=1000
P=10

zeta0=seq(0.05,0.4,by=.05)
alpha=.1
K=3
lambda=seq(0,2,.25)

tol=1e-16
verbose=F
nseed=50

c=c(5,7.5,10)

mu=3
rho=0
nu=10
pers = 0.99
K_true=3
perc_out=.05
out_sigma=100

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  K=K,
  lambda=lambda,
  TT=TT,
  P=P,
  c=c,
  K_true=K_true
)

ncores=parallel::detectCores()-1

rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

thres_out=0
thres_feat_weight=.02

start <- Sys.time()

res_list_K3 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    K <- hp$K[i]
    lambda <- hp$lambda[i]
    TT <- hp$TT[i]
    P <- hp$P[i]
    c <- hp$c[i]
    K_true <- hp$K_true[i]
    
    simDat <- sim_data_stud_t(
      seed = seed,
      TT = TT,
      P = P,
      Pcat = NULL,
      Ktrue = K_true,
      mu = mu,
      rho = rho,
      nu = nu,
      pers = pers
    )
    
    simDat_sparse <- simulate_sparse_hmm(
      Y = simDat$SimData,
      rel_ = rel_,
      true_stat = simDat$mchain,
      perc_out = perc_out,
      out_sigma = out_sigma,
      seed = seed
    )
    
    fit <- robust_sparse_jump(
      Y = as.matrix(simDat_sparse$Y),
      zeta0 = zeta0,
      lambda = lambda,
      K = K,
      tol = 1e-16,
      n_init = 3,
      n_outer = 10,
      alpha = 0.1,
      verbose = FALSE,
      knn = 10,
      c = c,
      M = NULL,
      hd=F
    )
    
    est_s <- fit$s
    est_s[fit$v <= thres_out] <- 0
    W_ind <- fit$W > thres_feat_weight
    truth <- simDat_sparse$truth
    W_truth <- simDat_sparse$W_truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    
    ARI_W <- tryCatch(
      mclust::adjustedRandIndex(W_ind, W_truth),
      error = function(e) 0
    )
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      c = c,
      W = fit$W,
      s = est_s,
      v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      ARI_W = ARI_W
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end <- Sys.time()


print(end-start)

save(res_list_K3,file='simple_simstud_rob_JM_K3.Rdata')

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_robJM/simple_simstud_rob_JM_K3.Rdata")

df_results_robJM_K3 <- do.call(rbind, lapply(res_list_K3, function(x) {
  x_sel <- x[c("seed", "zeta0", "K", "lambda", "TT", "P", "c", "ARI_s", "ARI_W")]
  as.data.frame(x_sel, stringsAsFactors = FALSE)
}))

library(dplyr)

av_results_robJM_K3 <- df_results_robJM_K3 %>%
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

# Supponendo che av_results_robJM_K3 sia il tuo tibble

# 1) superfice 3D per mean_ARI_s
fig_s <- plot_ly(
  data = av_results_robJM_K3[av_results_robJM_K3$K==3,],
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
  data = av_results_robJM_K3[av_results_robJM_K3$K==3,],
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




