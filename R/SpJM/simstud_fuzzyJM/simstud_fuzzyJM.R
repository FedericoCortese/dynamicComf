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


# mu=1, pers=.95 ----------------------------------------------------------


fuzzyJM_sim=list()

start_=Sys.time()

for(i in 1:n_chunks){
  hp_chunk=hp[(1+(i-1)*chunk_length):(chunk_length*i),]
  fuzzyJM_sim[(1+(i-1)*chunk_length):(chunk_length*i)]=
    parallel::mclapply(1:nrow(hp_chunk),
                       function(x)
                         simstud_fuzzyJM(seed=hp_chunk[x,]$seed,
                                         lambda=hp_chunk[x,]$lambda,
                                         TT=hp_chunk[x,]$TT,
                                         P=hp_chunk[x,]$P,
                                         K=3,mu=1,
                                         phi=.8,rho=0,
                                         Pcat=NULL,
                                         pers=.95
                                         ),
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



# Results -----------------------------------------------------------------

ARI_fuzzyJM=unlist(lapply(fuzzyJM_sim, function(x) {
  if (is.list(x) && !is.null(x$ARI)) {
    return(x$ARI)
  } else {
    return(NA)  # Return NA if structure is invalid
  }
}))

BAC_fuzzyJM=unlist(lapply(fuzzyJM_sim, function(x) {
  if (is.list(x) && !is.null(x$BAC)) {
    return(x$BAC)
  } else {
    return(NA)  # Return NA if structure is invalid
  }
}))

res=data.frame(hp,ARI_fuzzyJM,BAC_fuzzyJM)

# BAC
average_bac <- aggregate(BAC_fuzzyJM ~ TT + P + lambda, data = res, FUN = mean)
average_bac

best_lambda <- average_bac[with(average_bac, ave(BAC_fuzzyJM, 
                                                 TT, P, FUN = max) == BAC_fuzzyJM), 
                           c("TT", "P", "lambda", "BAC_fuzzyJM")]
best_lambda


ggplot(average_bac, aes(x = lambda, y = BAC_fuzzyJM)) +
  geom_line() +
  facet_grid(TT ~ P) +
  labs(x = "Lambda", y = "Average BAC", title = "BAC") +
  theme_minimal()

# ARI
average_ari <- aggregate(ARI_fuzzyJM ~ TT + P + lambda, data = res, FUN = mean)
average_ari

best_lambda <- average_ari[with(average_ari, ave(ARI_fuzzyJM, TT, P, FUN = max) == ARI_fuzzyJM), c("TT", "P", "lambda", "ARI_fuzzyJM")]
best_lambda


ggplot(average_ari, aes(x = lambda, y = ARI_fuzzyJM)) +
  geom_line() +
  facet_grid(TT ~ P) +
  labs(x = "Lambda", y = "Average ARI", title = "ARI") +
  theme_minimal()


# mu=.5 pers=.99 ----------------------------------------------------------
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
fuzzyJM_sim_scen2=list()

start_=Sys.time()

for(i in 1:n_chunks){
  hp_chunk=hp[(1+(i-1)*chunk_length):(chunk_length*i),]
  fuzzyJM_sim_scen2[(1+(i-1)*chunk_length):(chunk_length*i)]=
    parallel::mclapply(1:nrow(hp_chunk),
                       function(x)
                         simstud_fuzzyJM(seed=hp_chunk[x,]$seed,
                                         lambda=hp_chunk[x,]$lambda,
                                         TT=hp_chunk[x,]$TT,
                                         P=hp_chunk[x,]$P,
                                         K=3,mu=.5,
                                         phi=.8,rho=0,
                                         Pcat=NULL,
                                         pers=.99
                         ),
                       mc.cores = parallel::detectCores()-1)
  print(i)
  
}

end_=Sys.time()
elapsed_=end_-start_
save(fuzzyJM_sim_scen2,elapsed_,file="fuzzyJM_sim_scen2.RData")



# Results -----------------------------------------------------------------

ARI_fuzzyJM_scen2=unlist(lapply(fuzzyJM_sim_scen2, function(x) {
  if (is.list(x) && !is.null(x$ARI)) {
    return(x$ARI)
  } else {
    return(NA)  # Return NA if structure is invalid
  }
}))

BAC_fuzzyJM_scen2=unlist(lapply(fuzzyJM_sim_scen2, function(x) {
  if (is.list(x) && !is.null(x$BAC)) {
    return(x$BAC)
  } else {
    return(NA)  # Return NA if structure is invalid
  }
}))

res_scen2=data.frame(hp,ARI_fuzzyJM_scen2,BAC_fuzzyJM_scen2)

# BAC
average_bac <- aggregate(BAC_fuzzyJM_scen2 ~ TT + P + lambda, 
                         data = res_scen2, FUN = mean)
average_bac

best_lambda <- average_bac[with(average_bac, ave(BAC_fuzzyJM_scen2, 
                                                 TT, P, FUN = max) == 
                                  BAC_fuzzyJM_scen2), 
                           c("TT", "P", "lambda", "BAC_fuzzyJM_scen2")]
best_lambda


ggplot(average_bac, aes(x = lambda, y = BAC_fuzzyJM_scen2)) +
  geom_line() +
  facet_grid(TT ~ P) +
  labs(x = "Lambda", y = "Average BAC", title = "BAC") +
  theme_minimal()

# ARI
average_ari <- aggregate(ARI_fuzzyJM_scen2 ~ TT + P + lambda, 
                         data = res_scen2, FUN = mean)
average_ari

best_lambda <- average_ari[with(average_ari,
                                ave(ARI_fuzzyJM_scen2, TT, P, 
                                    FUN = max) == ARI_fuzzyJM_scen2), 
                           c("TT", "P", "lambda", "ARI_fuzzyJM_scen2")]
best_lambda


ggplot(average_ari, aes(x = lambda, y = ARI_fuzzyJM_scen2)) +
  geom_line() +
  facet_grid(TT ~ P) +
  labs(x = "Lambda", y = "Average ARI", title = "ARI") +
  theme_minimal()

