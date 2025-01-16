library(parallel)
source("Utils.R")

lambda=seq(0,1,by=.1)
TT=c(100,1000)
P=c(4,20,50)
seeds=1:100
hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

# Calcola i core per livello esterno ed interno
n_cores_total=parallel::detectCores()
n_cores_ext <- min(nrow(hp), floor(n_cores_total / 2))  # Massimo metà per livello esterno
n_cores_int <- floor((n_cores_total - n_cores_ext) / n_cores_ext)

# Setup 1
st=Sys.time()
contJM_setup1 <- parallel::mclapply(1:nrow(hp),
                                    function(x)
                                      simstud_contJM(seed=hp[x,]$seed,
                                                     lambda=hp[x,]$lambda,
                                                     TT=hp[x,]$TT,
                                                     P=hp[x,]$P,
                                                     Ktrue=3,mu=1,
                                                     phi=.8,rho=0,
                                                     Pcat=NULL,pers=.95,
                                                     pNAs=0,typeNA=3,
                                                     n_cores_int),
                                    mc.cores = n_cores_ext)
en=Sys.time()
