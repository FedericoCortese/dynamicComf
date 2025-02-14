library(parallel)
# library(snow)
# library(doSNOW)

#source("Utils.R")

lambda=seq(0,1,by=.05)
TT=c(50,100,500)
P=c(4,20,50)
seeds=1:100

hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

nrow_hp=nrow(hp)
n_chunks=10
chunk_length=nrow_hp/n_chunks

source("Utils_fuzzyJM.R")

fuzzyJM_sim=list()

start_=Sys.time()

for(i in 1:n_chunks){
  fuzzyJM_sim[(1+(i-1)*chunk_length):(chunk_length*i)]=
    parallel::mclapply(1:nrow(hp),
                       function(x)
                         simstud_fuzzyJM(seed=hp[x,]$seed,
                                         lambda=hp[x,]$lambda,
                                         TT=hp[x,]$TT,
                                         P=hp[x,]$P,
                                         K=3,mu=1,
                                         phi=.8,rho=0,
                                         Pcat=NULL,pers=.95),
                       mc.cores = parallel::detectCores()-1)
  print(i)
  
}

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
save(fuzzyJM_sim,elapsed_,file="fuzzyJM_sim.RData")


