library(parallel)
library(snow)
library(doSNOW)

# Simulation parameters --------------------------------------------------------------

source("Utils.R")

lambda=seq(0,1,by=.05)
TT=c(20,50,100)
P=c(10,20,50)
seeds=1:100

# lambda=seq(0,1,by=.5)
# TT=c(50)
# P=c(20)
# seeds=1:2

hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)
#head(hp)


# No missing --------------------------------------------------------------

start_no.miss=Sys.time()
# mixedJM_no.miss <- parallel::mclapply(1:nrow(hp),
#                                       function(x)
#                                         simstud_JMmixed(
#                                           seed=hp[x,]$seed,
#                                           lambda=hp[x,]$lambda,
#                                           TT=hp[x,]$TT,
#                                           P=hp[x,]$P,
#                                           Ktrue=3,mu=1,
#                                           phi=.8,rho=0,
#                                           Pcat=NULL,pers=.95,
#                                           pNAs=0,typeNA=2),
#                                       mc.cores = parallel::detectCores())

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
mixedJM_no.miss <- clusterApply(cl, 
                                1:nrow(hp), 
                                function(x)
                                  simstud_JMmixed(
                                    seed=hp[x,]$seed,
                                    lambda=hp[x,]$lambda,
                                    TT=hp[x,]$TT,
                                    P=hp[x,]$P,
                                    Ktrue=3,mu=1,
                                    phi=.8,rho=0,
                                    Pcat=NULL,pers=.95,
                                    pNAs=0,typeNA=2)
)
stopCluster(cl)

end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss


# Random missing ----------------------------------------------------------

pNAs=0.1
start_rand.miss10=Sys.time()
mixedJM_rand.miss10 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=0),
                                          mc.cores = parallel::detectCores())

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{library(RcppHMM)
#   library(reticulate)
#   library(pdfCluster)
#   library(boot)
#   library(xtable)
#   library(dplyr)
#   library(cluster)
#   library(gower)
#   library(StatMatch)})
# mixedJM_rand.miss10 <- clusterApply(cl, 
#                                 1:nrow(hp), 
#                                 function(x)
#                                   simstud_JMmixed(
#                                     seed=hp[x,]$seed,
#                                     lambda=hp[x,]$lambda,
#                                     TT=hp[x,]$TT,
#                                     P=hp[x,]$P,
#                                     Ktrue=3,mu=1,
#                                     phi=.8,rho=0,
#                                     Pcat=NULL,pers=.95,
#                                     pNAs=pNAs,typeNA=0)
# )
# stopCluster(cl)

end_rand.miss10=Sys.time()
elapsed_rand.miss10=end_rand.miss10-start_rand.miss10
save(mixedJM_rand.miss10,file="mixedJM_rand_miss10.RData")

pNAs=0.20
start_rand.miss20=Sys.time()
mixedJM_rand.miss20 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=0),
                                          mc.cores = parallel::detectCores())
end_rand.miss20=Sys.time()
elapsed_rand.miss20=end_rand.miss20-start_rand.miss20

pNAs=0.50
start_rand.miss50=Sys.time()
mixedJM_rand.miss50 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=0),
                                          mc.cores = parallel::detectCores())
end_rand.miss50=Sys.time()
elapsed_rand.miss50=end_rand.miss50-start_rand.miss50


# Continuous missing ------------------------------------------------------

pNAs=0.1
start_cont.miss10=Sys.time()
mixedJM_cont.miss10 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=1),
                                          mc.cores = parallel::detectCores())
end_cont.miss10=Sys.time()
elapsed_cont.miss10=end_cont.miss10-start_cont.miss10

pNAs=0.20
start_cont.miss20=Sys.time()
mixedJM_cont.miss20 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=1),
                                          mc.cores = parallel::detectCores())
end_cont.miss20=Sys.time()
elapsed_cont.miss20=end_cont.miss20-start_cont.miss20

pNAs=0.50
start_cont.miss50=Sys.time()
mixedJM_cont.miss50 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=1),
                                          mc.cores = parallel::detectCores())
end_cont.miss50=Sys.time()
elapsed_cont.miss50=end_cont.miss50-start_cont.miss50

# Save -------------------------------------------------------------------
save.image("simres_mixedJM.RData")


# Evaluation --------------------------------------------------------------

res=data.frame(hp,ARI=unlist(lapply(mixedJM_no.miss,function(x)x$ARI)))
head(res)
library(dplyr)

maxres=res%>%group_by(TT,P)%>%summarise(maxARI=max(ARI,na.rm=T))

avres=res%>%group_by(TT,P,lambda)%>%summarise(avARI=median(ARI,na.rm=T))

avres%>%group_by(TT,P)%>%summarise(maxARI=max(avARI),avlambda=mean(lambda))
