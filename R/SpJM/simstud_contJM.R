library(parallel)
source("Utils.R")

lambda=seq(0,1,by=.1)
TT=c(100,1000)
# Errore quando vario i P, potrei separare gli studi simulati per ogni P
# P=c(4,20,50)
P=20
seeds=1:100
hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

# Calcola i core per livello esterno ed interno
n_cores_total=parallel::detectCores()
frac_ext=1/3
n_cores_ext=n_cores_total*frac_ext
n_cores_int=(n_cores_total-n_cores_ext)/n_cores_ext

# n_cores_ext=30

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
                                                     prll=T,
                                                     n_cores_int=n_cores_int),
                                    mc.cores = n_cores_ext)
en=Sys.time()
