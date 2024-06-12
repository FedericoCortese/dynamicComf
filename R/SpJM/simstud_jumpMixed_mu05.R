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

rho=0
mu=.5

# No missing --------------------------------------------------------------

start_no.miss=Sys.time()
mixedJM_mu_no.miss <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=.5,
                                              phi=.8,rho=rho,
                                              Pcat=NULL,pers=.95,
                                              pNAs=0,typeNA=2),
                                          mc.cores = parallel::detectCores()-1)

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
# mixedJM_mu_no.miss <- clusterApply(cl, 
#                                 1:nrow(hp), 
#                                 function(x)
#                                   simstud_JMmixed(
#                                     seed=hp[x,]$seed,
#                                     lambda=hp[x,]$lambda,
#                                     TT=hp[x,]$TT,
#                                     P=hp[x,]$P,
#                                     Ktrue=3,mu=.5,
#                                     phi=.8,rho=.5,
#                                     Pcat=NULL,pers=.95,
#                                     pNAs=0,typeNA=2)
# )
# stopCluster(cl)

end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss
save(mixedJM_mu_no.miss,elapsed_no.miss,file="mixedJM_mu_no_miss.RData")

# Classification performance comparison --------------------------------------------------------------

hp_comp=expand.grid(TT=TT,P=P,seed=seeds)

# KMeMo
start_no.miss_kMeMo=Sys.time()
mixedJM_mu_no.miss_kMeMo <- parallel::mclapply(1:nrow(hp_comp),
                                                function(x)
                                                  simstud_JMmixed(
                                                    seed=hp_comp[x,]$seed,
                                                    lambda=0,
                                                    TT=hp_comp[x,]$TT,
                                                    P=hp_comp[x,]$P,
                                                    Ktrue=3,
                                                    mu=.5,
                                                    phi=.8,
                                                    rho=rho,
                                                    Pcat=NULL,
                                                    pers=.95,
                                                    pNAs=0,
                                                    typeNA=2)
                                                ,
                                                mc.cores = parallel::detectCores()-1)

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
# mixedJM_mu_no.miss_kMeMo <- clusterApply(cl,
#                                       1:nrow(hp_comp),
#                                       function(x)
#                                         simstud_JMmixed(
#                                           seed=hp_comp[x,]$seed,
#                                           lambda=0,
#                                           TT=hp_comp[x,]$TT,
#                                           P=hp_comp[x,]$P,
#                                           Ktrue=3,
#                                           mu=.5,
#                                           phi=.8,
#                                           rho=.5,
#                                           Pcat=NULL,
#                                           pers=.95,
#                                           pNAs=0,
#                                           typeNA=2)
# )
# stopCluster(cl)

end_no.miss_kMeMo=Sys.time()
elapsed_no.miss_kMeMo=end_no.miss_kMeMo-start_no.miss_kMeMo
save(mixedJM_mu_no.miss_kMeMo,elapsed_no.miss_kMeMo,file="mixedJM_mu_no_miss_kMeMo.RData")
rm(mixedJM_mu_no.miss_kMeMo,elapsed_no.miss_kMeMo)

# Spectral Clustering
library(SpectralClMixed)
start_no.miss_specluster=Sys.time()
mixedJM_mu_no.miss_specluster <- parallel::mclapply(1:nrow(hp_comp),
                                                     function(x)
                                                       simstud_speclust(
                                                         seed=hp_comp[x,]$seed,
                                                         TT=hp_comp[x,]$TT,
                                                         P=hp_comp[x,]$P,
                                                         Ktrue=3,
                                                         mu=.5,
                                                         phi=.8,
                                                         rho=rho,
                                                         Pcat=NULL,
                                                         pers=.95,
                                                         pNAs=0,
                                                         typeNA=2),
                                                     mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(RcppHMM)
#   library(reticulate)
#   library(pdfCluster)
#   library(boot)
#   library(xtable)
#   library(dplyr)
#   library(cluster)
#   library(gower)
#   library(StatMatch)
#   library(SpectralClMixed)
# })
# mixedJM_mu_no.miss_specluster <- clusterApply(cl,
#                                            1:nrow(hp_comp),
#                                            function(x)
#                                              simstud_speclust(
#                                                seed=hp_comp[x,]$seed,
#                                                TT=hp_comp[x,]$TT,
#                                                P=hp_comp[x,]$P,
#                                                Ktrue=3,
#                                                mu=.5,
#                                                phi=.8,
#                                                rho=.5,
#                                                Pcat=NULL,
#                                                pers=.95,
#                                                pNAs=0,
#                                                typeNA=2)
# )
# stopCluster(cl)

end_no.miss_specluster=Sys.time()
elapsed_no.miss_specluster=end_no.miss_specluster-start_no.miss_specluster
save(mixedJM_mu_no.miss_specluster,elapsed_no.miss_specluster,file="mixedJM_mu_no_miss_specluster.RData")
rm(mixedJM_mu_no.miss_specluster,elapsed_no.miss_specluster)


# Evaluation --------------------------------------------------------------


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

res_eval(mixedJM_mu_no.miss,hp)
res_eval(mixedJM_mu_no.miss_kMeMo,hp,lambda0=T)
res_eval(mixedJM_mu_no.miss_specluster,hp,lambda0=T)


