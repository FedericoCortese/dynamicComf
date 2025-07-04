source("Utils_sparse_robust_2.R")

alpha=.1
zeta0=0.15
lambda=0.3
verbose=T

TT=1000
P=10
K=3

rel=list()
rel[[1]]=c(1,2,3)
rel[[2]]=c(1,3,4)
rel[[3]]=c(1,4,5)

ssim=simulate_sparse_hmm(seed=sample(1:1000,1), 
                         TT=TT, 
                         P=P, 
                         K=K,
                                rel=rel,
                                mu   = 3,
                                rho  = 0.2,
                                nu   = 100,
                                phi  = 0.8,
                                pers = 0.99,
                                sd_noise = 100,
                                out_pct     = 0.05,
                                out_scale   = 200)

Y=ssim$Y
head(Y)
truth=ssim$truth
W_truth=ssim$W_truth

# Plot Y with colors as in truth
x11()
par(mfrow=c(3,3))
for(i in 1:9){
plot(Y[,i],col=truth+1,ylab = paste("Y",i,sep=""),xlab = "Time",ylim = c(-7,7))
}

startR=Sys.time()
prv=robust_sparse_jump(Y=Y,
                   zeta0=zeta0,
                   lambda=lambda,
                   K=K,
                   tol        = 1e-16,
                   n_init     = 10,
                   n_outer    = 20,
                   alpha      = 0.1,
                   verbose    = T,
                   knn        = 10,
                   c          = 5,
                   M          = NULL,
                   scale="i")
endR=Sys.time()
print(endR-startR)

table(prv$s,truth)

round(prv$W,2)

prv$s[prv$v<0.5]=0
table(prv$s,truth)


