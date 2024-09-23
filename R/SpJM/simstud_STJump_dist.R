library(parallel)
# library(snow)
# library(doSNOW)

source("Utils.R")

mu=.5
rho=.2
phi=.8
beta=.9
theta=.01
k=3
P=20
Pcat=10


# 20% gaps, no missing ----------------------------------------------------

pNAs=0
pg=.2

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
M=c(10,50)
TT=c(10,50)

hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)

start_STJsim=Sys.time()
STJsim <- parallel::mclapply(1:nrow(hp),
                             function(x)
                               simstud_STJump_dist(lambda=hp[x,]$lambda,
                                                   gamma=hp[x,]$gamma,
                                                   seed=hp[x,]$seed,
                                                   M=hp[x,]$M,
                                                   TT=hp[x,]$TT,
                                                   beta=beta, 
                                                   theta=theta,
                                                   mu=mu,
                                                   rho=rho,
                                                   K=K,P=P,
                                                   phi=phi,
                                                   Pcat=Pcat,
                                                   pNAs=pNAs,
                                                   pg=pg),
                             mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim,elapsed_STJsim,file="STJsim_dist.RData")