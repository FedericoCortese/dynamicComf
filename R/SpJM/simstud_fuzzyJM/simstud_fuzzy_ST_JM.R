library(parallel)
# library(snow)
# library(doSNOW)

#source("Utils.R")

lambda=seq(0,.5,by=.05)
gamma=seq(0,.5,by=.05)
TT=100
#TT=c(50,100)
#P=c(4,20)
P=12
M=c(50,4)
seeds=1:100

hp=expand.grid(lambda=lambda,gamma=gamma,TT=TT,P=P,M=M,seed=seeds)

# n_cores_total=parallel::detectCores()
# frac_ext=1/3
# n_cores_ext=n_cores_total*frac_ext
# n_cores_int=n_cores_int=n_cores_int=(n_cores_total-n_cores_ext)/n_cores_ext
# n_cores_ext=ceiling(n_cores_total*frac_ext)

# nrow_hp=nrow(hp)
# n_chunks=10
# chunk_length=nrow_hp/n_chunks

source("Utils_fuzzyJM.R")


start_=Sys.time()
fuzzy_STJM_sim=parallel::mclapply(1:nrow(hp),
                   function(x)
                   simstud_fuzzySTJM(seed=hp[x,]$seed,
                                     lambda=hp[x,]$lambda,
                                     gamma = hp[x,]$gamma,
                                     TT=hp[x,]$TT,
                                     P=hp[x,]$P,
                                     Pcat=NULL,
                                     M=hp[x,]$M,
                                     beta=.9,
                                     #theta=.01,
                                     theta=1,
                                              mu=1,rho=0,
                                              K=3,phi=.8,pNAs=0,pg=0,
                                              ncores_M=NULL),
                   mc.cores = detectCores()-1)

end_=Sys.time()
elapsed_=end_-start_
save(fuzzy_STJM_sim,elapsed_,file="fuzzy_STJM_sim.RData")
