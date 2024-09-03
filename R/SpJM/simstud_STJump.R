source("Utils.R")

K=3
P=30
###
mu=1
rho=0.2
###
phi=.8
Pcat=10
pNAs=0

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:2
M=100
TT=10
hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)
# seed=1:100
# M=c(25,100)
# TT=c(10,50)


start_STJsim=Sys.time()
STJsim <- parallel::mclapply(1:nrow(hp),
                             function(x)
                               STJump_sim(lambda=hp[x,]$lambda,
                                          gamma=hp[x,]$gamma,
                                          seed=hp[x,]$seed,
                                          M=hp[x,]$M,
                                          TT=hp[x,]$TT,
                                          mu=mu,rho=rho,
                                          K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs),
                             mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim,elapsed_STJsim,file="STJsim.RData")