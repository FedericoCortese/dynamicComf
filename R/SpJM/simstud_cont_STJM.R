library(parallel)
source("Utils.R")

P=10
TT=30
M=10
lambda=seq(0,1,by=.1)
gamma=seq(0,1,by=.1)
seeds=1:100
hp=expand.grid(TT=TT,M=M,P=P,
               lambda=lambda,gamma=gamma,
               seed=seeds)

hp=hp[1:30,]

n_cores_total=parallel::detectCores()
frac_ext=1/3
n_cores_ext=n_cores_total*frac_ext
n_cores_int=(n_cores_total-n_cores_ext)/n_cores_ext


# Setup 1
st=Sys.time()
cont_STJM_setup1_P10 <- parallel::mclapply(1:nrow(hp),
                                        function(x)
                                          simstud_cont_STJM(seed=hp[x,]$seed,
                                                            lambda=hp[x,]$lambda,
                                                            TT=hp[x,]$TT,P=hp[x,]$P,
                                                            M=hp[x,]$M,
                                                            Ktrue=3,mu=1,
                                                            phi=.8,rho=0,
                                                            theta=.9,beta=.01,
                                                            Pcat=NULL,pers=.95,
                                                            pNAs=0,typeNA=3,
                                                            n_cores_int=n_cores_int,
                                                            prll=T),
                                        mc.cores = n_cores_ext)
en=Sys.time()
elapsed=en-st
save(contJM_setup1_P20,elapsed,file="contJM_setup1_P20.Rdata")


