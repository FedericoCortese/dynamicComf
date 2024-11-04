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


# Spectral clustering -----------------------------------------------------

hp_comp=expand.grid(TT=TT,P=P,seed=seeds)




# Results -----------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres/mixJM gaps/mixedJM_gaps.RData")
load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres/mixJM gaps/mixedJM_gaps_mu05.RData")
load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres/mixJM gaps/mixedJM_gaps_rho02.RData")

res_eval=function(res_obj,hp,lambda0=F,ARI=T){
  library(dplyr)
  
  if(ARI){
    res=data.frame(hp,ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                   imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err)))
    )
    
    # maxres=res%>%group_by(TT,P)%>%summarise(maxARI=max(ARI,na.rm=T))
    
    if(lambda0){
      res=res[which(res$lambda==0),]
      res%>%group_by(TT,P)%>%summarise(avARI=median(ARI,na.rm=T),
                                       avErr=mean(imput.err))
    }
    
    else{
      avres=res%>%group_by(TT,P,lambda)%>%summarise(avARI=median(ARI,na.rm=T),
                                                    avErr=mean(imput.err))
      
      avres%>%group_by(TT,P)%>%summarise(maxARI=max(avARI),
                                         lambda=lambda[which.max(avARI)],
                                         avErr=avErr[which.max(avARI)])
    }
  }
  else{
    res=data.frame(hp,
                   #ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                   imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err)))
    )
    res%>%group_by(TT,P)%>%summarise(
      #avARI=median(ARI,na.rm=T),
      avErr=mean(imput.err))
  }
  
}

res_eval(mixedJM_gaps,hp)
res_eval(mixedJM_gaps,hp,lambda0=T)
res_eval(mixedJM_no.miss_specluster,hp,lambda0=T)

res_eval(mixedJM_gaps_mu05,hp,lambda0=F)
res_eval(mixedJM_gaps_mu05,hp,lambda0=T)
res_eval(mixedJM_no.miss_specluster_mu05,hp,lambda0=T)


res_eval(mixedJM_gaps_rho02,hp,lambda0=F)
res_eval(mixedJM_gaps_rho02,hp,lambda0=T)
res_eval(mixedJM_no.miss_specluster_rho02,hp,lambda0=T)

