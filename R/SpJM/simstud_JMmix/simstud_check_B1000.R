source("Utils.R")

lambda=seq(0,1,by=.05)
TT=c(100)
P=c(25)
seeds=1:1000
hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

start_no.miss=Sys.time()
mixedJM_no.miss_setup1 <- parallel::mclapply(1:nrow(hp),
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
                                      mc.cores = parallel::detectCores()-1)
end_no.miss=Sys.time()
elapsed_no.miss_setup1=end_no.miss-start_no.miss
save(mixedJM_no.miss_setup1,elapsed_no.miss_setup1,file="mixedJM_no_miss_setup1.RData")

start_no.miss=Sys.time()
mixedJM_no.miss_setup3 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_JMmixed(
                                                 seed=hp[x,]$seed,
                                                 lambda=hp[x,]$lambda,
                                                 TT=hp[x,]$TT,
                                                 P=hp[x,]$P,
                                                 Ktrue=3,mu=.5,
                                                 phi=.8,rho=0,
                                                 Pcat=NULL,pers=.95,
                                                 pNAs=0,typeNA=2),
                                             mc.cores = parallel::detectCores()-1)
end_no.miss=Sys.time()
elapsed_no.miss_setup3=end_no.miss-start_no.miss
save(mixedJM_no.miss_setup3,elapsed_no.miss_setup3,file="mixedJM_no_miss_setup3.RData")

## SpClust
library(SpectralClMixed)
start_no.miss_specluster=Sys.time()
mixedJM_specluster_setup1 <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_speclust(
                                          seed=hp[x,]$seed,
                                          TT=hp[x,]$TT,
                                          P=hp[x,]$P,
                                          Ktrue=3,
                                          mu=1,
                                          phi=.8,
                                          rho=0,
                                          Pcat=NULL,
                                          pers=.95,
                                          pNAs=0,
                                          typeNA=2),
                                      mc.cores = parallel::detectCores()-1)

end_no.miss_specluster=Sys.time()
elapsed_no.miss_specluster_setup1=end_no.miss_specluster-start_no.miss_specluster
save(mixedJM_specluster_setup1,elapsed_no.miss_specluster_setup1,file="mixedJM_no_miss_specluster_setup1.RData")

start_no.miss_specluster=Sys.time()
mixedJM_specluster_setup3 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_speclust(
                                                 seed=hp[x,]$seed,
                                                 TT=hp[x,]$TT,
                                                 P=hp[x,]$P,
                                                 Ktrue=3,
                                                 mu=.5,
                                                 phi=.8,
                                                 rho=0,
                                                 Pcat=NULL,
                                                 pers=.95,
                                                 pNAs=0,
                                                 typeNA=2),
                                             mc.cores = parallel::detectCores()-1)

end_no.miss_specluster=Sys.time()
elapsed_no.miss_specluster_setup3=end_no.miss_specluster-start_no.miss_specluster
save(mixedJM_specluster_specluster_setup3,elapsed_no.miss_specluster_setup3,file="mixedJM_no_miss_specluster_setup3.RData")

library(dplyr)
res_eval=function(res_obj,hp,lambda0=F,ARI=T){
  
  if(ARI){
    res=data.frame(hp,ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                   imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err)))
    )
    
    # maxres=res%>%group_by(TT,P)%>%summarise(maxARI=max(ARI,na.rm=T))
    
    if(lambda0){
      res=res[which(res$lambda==0),]
      res%>%group_by(TT,P)%>%summarise(avARI=median(ARI,na.rm=T),
                                       sdARI=sd(ARI,na.rm=T),
                                       avErr=mean(imput.err))
    }
    
    else{
      avres=res%>%group_by(TT,P,lambda)%>%summarise(avARI=median(ARI,na.rm=T),
                                                    sdARI=sd(ARI,na.rm=T),
                                                    avErr=mean(imput.err))
      
      avres%>%group_by(TT,P)%>%summarise(maxARI=max(avARI),
                                         sdARI=min(sdARI),
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

res_eval(mixedJM_no.miss_setup1,hp)
res_eval(mixedJM_no.miss_setup1,hp,lambda0=T)
res_eval(mixedJM_specluster_setup1,hp)

res_eval(mixedJM_no.miss_setup3,hp)
res_eval(mixedJM_no.miss_setup3,hp,lambda0=T)
res_eval(mixedJM_no.miss_specluster_setup3,hp)