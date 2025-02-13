library(parallel)
# library(snow)
# library(doSNOW)

#source("Utils.R")


# P=4 ---------------------------------------------------------------------


lambda=seq(0,1,by=.05)
TT=c(50,100,500)
#P=c(4,20,50)
P=4
seeds=1:100

hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)
#hp=hp[1:10,]

source("Utils_fuzzyJM.R")

start_=Sys.time()
fuzzyJM_sim_P4 <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_fuzzyJM(seed=hp[x,]$seed,
                                        lambda=hp[x,]$lambda,
                                        TT=hp[x,]$TT,
                                        P=hp[x,]$P,
                                        K=3,mu=1,
                                        phi=.8,rho=0,
                                        Pcat=NULL,pers=.95),
                                      mc.cores = parallel::detectCores()-1)

# source("Utils_fuzzyJM.R")
# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(cluster)
#   library(StatMatch)
# })
# fuzzyJM_sim <- clusterApply(cl,
#                             1:nrow(hp),
#                             function(x)
#                               simstud_fuzzyJM(seed=hp[x,]$seed,
#                                               lambda=hp[x,]$lambda,
#                                               TT=hp[x,]$TT,
#                                               P=hp[x,]$P,
#                                               K=3,mu=1,
#                                               phi=.8,rho=0,
#                                               Pcat=NULL,pers=.95)
# )
# stopCluster(cl)


end_=Sys.time()
elapsed_=end_-start_
save(fuzzyJM_sim_P4,elapsed_,file="fuzzyJM_sim_P4.RData")


