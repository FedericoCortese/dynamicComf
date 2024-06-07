library(parallel)
library(snow)
library(doSNOW)

# Simulation parameters --------------------------------------------------------------

source("Utils.R")

lambda=seq(0,1,by=.05)
TT=c(50,100,500)
P=c(25,50,75)
seeds=1:100

# lambda=seq(0,1,by=.5)
# TT=c(50)
# P=c(20)
# seeds=1:2

hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)
#head(hp)


# No missing --------------------------------------------------------------

start_no.miss=Sys.time()
mixedJM_no.miss <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_JMmixed(
                                          seed=hp[x,]$seed,
                                          lambda=hp[x,]$lambda,
                                          TT=hp[x,]$TT,
                                          P=hp[x,]$P,
                                          Ktrue=3,mu=1,
                                          phi=.8,rho=0,
                                          Pcat=NULL,pers=.95,
                                          pNAs=0,typeNA=2),
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
# mixedJM_no.miss <- clusterApply(cl, 
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
#                                     pNAs=0,typeNA=2)
# )
# stopCluster(cl)

end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss
save(mixedJM_no.miss,elapsed_no.miss,file="mixedJM_no_miss.RData")

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

# pNAs=0.50
# start_rand.miss50=Sys.time()
# mixedJM_rand.miss50 <- parallel::mclapply(1:nrow(hp),
#                                           function(x)
#                                             simstud_JMmixed(
#                                               seed=hp[x,]$seed,
#                                               lambda=hp[x,]$lambda,
#                                               TT=hp[x,]$TT,
#                                               P=hp[x,]$P,
#                                               Ktrue=3,mu=1,
#                                               phi=.8,rho=0,
#                                               Pcat=NULL,pers=.95,
#                                               pNAs=pNAs,typeNA=0),
#                                           mc.cores = parallel::detectCores())
# end_rand.miss50=Sys.time()
# elapsed_rand.miss50=end_rand.miss50-start_rand.miss50


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

# pNAs=0.50
# start_cont.miss50=Sys.time()
# mixedJM_cont.miss50 <- parallel::mclapply(1:nrow(hp),
#                                           function(x)
#                                             simstud_JMmixed(
#                                               seed=hp[x,]$seed,
#                                               lambda=hp[x,]$lambda,
#                                               TT=hp[x,]$TT,
#                                               P=hp[x,]$P,
#                                               Ktrue=3,mu=1,
#                                               phi=.8,rho=0,
#                                               Pcat=NULL,pers=.95,
#                                               pNAs=pNAs,typeNA=1),
#                                           mc.cores = parallel::detectCores())
# end_cont.miss50=Sys.time()
# elapsed_cont.miss50=end_cont.miss50-start_cont.miss50

# Save -------------------------------------------------------------------
save.image("simres_mixedJM.RData")


# Evaluation --------------------------------------------------------------

library(dplyr)

res_eval=function(res_obj,hp){
  library(dplyr)
  
  res=data.frame(hp,ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                 imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err)))
                 )
  
  # maxres=res%>%group_by(TT,P)%>%summarise(maxARI=max(ARI,na.rm=T))
  
  avres=res%>%group_by(TT,P,lambda)%>%summarise(avARI=median(ARI,na.rm=T),
                                                avErr=mean(imput.err))
  
  avres%>%group_by(TT,P)%>%summarise(maxARI=max(avARI),
                                     lambda=lambda[which.max(avARI)],
                                     avErr=avErr[which.max(avARI)])
  
}

res_eval(mixedJM_no.miss,hp)

res_eval(mixedJM_cont.miss10,hp)
res_eval(mixedJM_cont.miss20,hp)
#res_eval(mixedJM_cont.miss50,hp)

res_eval(mixedJM_rand.miss10,hp)
res_eval(mixedJM_rand.miss20,hp)
#res_eval(mixedJM_rand.miss50,hp)
