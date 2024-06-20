library(doParallel)
library(snow)
library(doSNOW)
###
lambda=seq(0,1,by=.05)
TT=c(50,100,500)
P=c(25,50,75)
seeds=1:100

# lambda=seq(0,1,by=.5)
# TT=c(50)
# P=c(20)
# seeds=1:2
source("Utils.R")

hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)


# mu=1 rho=0 --------------------------------------------------------------

# start_gaps=Sys.time()
# mixedJM_gaps <- parallel::mclapply(1:nrow(hp),
#                                       function(x)
#                                         simstud_JMmixed2(
#                                           seed=hp[x,]$seed,
#                                           lambda=hp[x,]$lambda,
#                                           TT=hp[x,]$TT,
#                                           P=hp[x,]$P,
#                                           Ktrue=3,mu=1,
#                                           phi=.8,rho=0,
#                                           Pcat=NULL,pers=.95,
#                                           pNAs=0,typeNA=2,timeflag=T),
#                                       mc.cores = parallel::detectCores()-1)
# 
# 
# end_gaps=Sys.time()
# elapsed_gaps=end_gaps-start_gaps
# save(mixedJM_gaps,elapsed_gaps,file="mixedJM_gaps.RData")

start_gaps=Sys.time()
cl<-makeCluster(parallel::detectCores(),type="SOCK")
parallel::clusterExport(cl,ls())
parallel::clusterEvalQ(cl,{library(RcppHMM)
  library(reticulate)
  library(pdfCluster)
  library(boot)
  library(xtable)
  library(dplyr)
  library(cluster)
  library(gower)
  library(StatMatch)})
mixedJM_gaps <- clusterApply(cl,
                             1:nrow(hp),
                             function(x)
                               simstud_JMmixed2(
                                 seed=hp[x,]$seed,
                                 lambda=hp[x,]$lambda,
                                 TT=hp[x,]$TT,
                                 P=hp[x,]$P,
                                 Ktrue=3,mu=1,
                                 phi=.8,rho=0,
                                 Pcat=NULL,pers=.95,
                                 pNAs=0,typeNA=2,timeflag=T)
)
stopCluster(cl)
end_gaps=Sys.time()
elapsed_gaps=end_gaps-start_gaps
save(mixedJM_gaps,elapsed_gaps,file="mixedJM_gaps.RData")
rm(mixedJM_gaps,elapsed_gaps)


# mu=1 rho=0.2 ------------------------------------------------------------

# start_gaps_rho02=Sys.time()
# mixedJM_gaps_rho02 <- parallel::mclapply(1:nrow(hp),
#                                    function(x)
#                                      simstud_JMmixed2(
#                                        seed=hp[x,]$seed,
#                                        lambda=hp[x,]$lambda,
#                                        TT=hp[x,]$TT,
#                                        P=hp[x,]$P,
#                                        Ktrue=3,mu=1,
#                                        phi=.8,rho=0.2,
#                                        Pcat=NULL,pers=.95,
#                                        pNAs=0,typeNA=2,timeflag=T),
#                                    mc.cores = parallel::detectCores()-1)
# 
# 
# end_gaps_rho02=Sys.time()
# elapsed_gaps_rho02=end_gaps_rho02-start_gaps_rho02
# save(mixedJM_gaps_rho02,elapsed_gaps_rho02,file="mixedJM_gaps_rho02.RData")

start_gaps_rho02=Sys.time()
cl<-makeCluster(parallel::detectCores(),type="SOCK")
parallel::clusterExport(cl,ls())
parallel::clusterEvalQ(cl,{library(RcppHMM)
  library(reticulate)
  library(pdfCluster)
  library(boot)
  library(xtable)
  library(dplyr)
  library(cluster)
  library(gower)
  library(StatMatch)})
mixedJM_gaps_rho02 <- clusterApply(cl,
                                   1:nrow(hp),
                                   function(x)
                                     simstud_JMmixed2(
                                       seed=hp[x,]$seed,
                                       lambda=hp[x,]$lambda,
                                       TT=hp[x,]$TT,
                                       P=hp[x,]$P,
                                       Ktrue=3,mu=1,
                                       phi=.8,rho=0.2,
                                       Pcat=NULL,pers=.95,
                                       pNAs=0,typeNA=2,timeflag=T)
)
stopCluster(cl)
end_gaps_rho02=Sys.time()
elapsed_gaps_rho02=end_gaps_rho02-start_gaps_rho02
save(mixedJM_gaps_rho02,elapsed_gaps_rho02,file="mixedJM_gaps_rho02.RData")
rm(mixedJM_gaps_rho02,elapsed_gaps_rho02)


# mu=0.5 rho=0 ------------------------------------------------------------

# start_gaps_mu05=Sys.time()
# mixedJM_gaps_mu05 <- parallel::mclapply(1:nrow(hp),
#                                          function(x)
#                                            simstud_JMmixed2(
#                                              seed=hp[x,]$seed,
#                                              lambda=hp[x,]$lambda,
#                                              TT=hp[x,]$TT,
#                                              P=hp[x,]$P,
#                                              Ktrue=3,mu=.5,
#                                              phi=.8,rho=0,
#                                              Pcat=NULL,pers=.95,
#                                              pNAs=0,typeNA=2,timeflag=T),
#                                          mc.cores = parallel::detectCores()-1)
# 
# 
# end_gaps_mu05=Sys.time()
# elapsed_gaps_mu05=end_gaps_mu05-start_gaps_mu05
# save(mixedJM_gaps_mu05,elapsed_gaps_mu05,file="mixedJM_gaps_mu05.RData")
start_gaps_mu05=Sys.time()
cl<-makeCluster(parallel::detectCores(),type="SOCK")
parallel::clusterExport(cl,ls())
parallel::clusterEvalQ(cl,{library(RcppHMM)
  library(reticulate)
  library(pdfCluster)
  library(boot)
  library(xtable)
  library(dplyr)
  library(cluster)
  library(gower)
  library(StatMatch)})
mixedJM_gaps_mu05 <- clusterApply(cl,
                                  1:nrow(hp),
                                  function(x)
                                    simstud_JMmixed2(
                                      seed=hp[x,]$seed,
                                      lambda=hp[x,]$lambda,
                                      TT=hp[x,]$TT,
                                      P=hp[x,]$P,
                                      Ktrue=3,mu=.5,
                                      phi=.8,rho=0,
                                      Pcat=NULL,pers=.95,
                                      pNAs=0,typeNA=2,timeflag=T)
)
stopCluster(cl)
end_gaps_mu05=Sys.time()
elapsed_gaps_mu05=end_gaps_mu05-start_gaps_mu05
save(mixedJM_gaps_mu05,elapsed_gaps_mu05,file="mixedJM_gaps_mu05.RData")
rm(mixedJM_gaps_mu05,elapsed_gaps_mu05)
