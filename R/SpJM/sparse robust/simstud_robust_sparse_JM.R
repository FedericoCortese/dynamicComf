TT=500
P=50
Pcat=NULL
lambda=.1
K=3
mu=1.5
rho=0
nu=4
phi=.8
pers=.99
n_init=5
n_inner=10
n_outer=10
tol=NULL
verbose=T
alpha=.1


source("Utils_sparse_robust.R")

simDat=sim_data_stud_t(seed=123,
                       TT=TT,
                       P=P,
                       Pcat=Pcat,
                       Ktrue=K,
                       mu=mu,
                       rho=rho,
                       nu=nu,
                       phi=phi,
                       pers=pers)

Y=simDat$SimData
feat_list=list()
feat_list[[1]]=c(1,2,3)
feat_list[[2]]=c(3,4,5)
feat_list[[3]]=c(5,6,7)
Y_noised=apply_noise_by_cluster(Y,simDat$mchain,feat_list)

Y=Y_noised

zeta0=3
lambda=.1

fit=sparse_robust_fit(Y=Y,K=K,
                      lambda=lambda,
                      zeta0=zeta0,
                      alpha=.2,verbose=T,tol=1e-4,
                      n_init=5,n_outer=10)

round(fit$W,2)
mclust::adjustedRandIndex(simDat$mchain,fit$best_s)

start_time=Sys.time()
temp=weight_inv_exp_dist(Y,s,W,zeta0)
end_time=Sys.time()
end_time-start_time

start_new_time=Sys.time()
temp_new=weight_inv_exp_dist_new(Y,s,W,zeta0)
end_new_time=Sys.time()
end_new_time-start_new_time




