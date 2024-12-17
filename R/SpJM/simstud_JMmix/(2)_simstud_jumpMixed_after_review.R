library(parallel)
library(snow)
library(doSNOW)
source("Utils.R")

# Complete data -------------------------------


lambda=seq(0,1,by=.05)
TT=c(50,100,500)
P=c(4,20,50)
seeds=1:100
hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

# Setup 1
start_no.miss=Sys.time()
mixedJM_no.miss_setup1 <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_JMmixed_missmec(
                                          seed=hp[x,]$seed,
                                          lambda=hp[x,]$lambda,
                                          TT=hp[x,]$TT,
                                          P=hp[x,]$P,
                                          Ktrue=3,mu=1,
                                          phi=.8,rho=0,
                                          Pcat=NULL,pers=.95,
                                          pNAs=0,typeNA=3),
                                      mc.cores = parallel::detectCores())
mixedJM_no.miss_setup2 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_JMmixed_missmec(
                                                 seed=hp[x,]$seed,
                                                 lambda=hp[x,]$lambda,
                                                 TT=hp[x,]$TT,
                                                 P=hp[x,]$P,
                                                 Ktrue=3,mu=1,
                                                 phi=.8,rho=0.2,
                                                 Pcat=NULL,pers=.95,
                                                 pNAs=0,typeNA=3),
                                             mc.cores = parallel::detectCores())
mixedJM_no.miss_setup3 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_JMmixed_missmec(
                                                 seed=hp[x,]$seed,
                                                 lambda=hp[x,]$lambda,
                                                 TT=hp[x,]$TT,
                                                 P=hp[x,]$P,
                                                 Ktrue=3,mu=0.5,
                                                 phi=.8,rho=0,
                                                 Pcat=NULL,pers=.95,
                                                 pNAs=0,typeNA=3),
                                             mc.cores = parallel::detectCores())
end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss
save(mixedJM_no.miss_setup1,mixedJM_no.miss_setup2,mixedJM_no.miss_setup3,
     elapsed_no.miss,file="mixedJM_no_miss.RData")

# Spectral Clustering
library(SpectralClMixed)
start_no.miss_specluster=Sys.time()
spClust_no.miss_setup1 <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_speclust(
                                          seed=hp_comp[x,]$seed,
                                          TT=hp_comp[x,]$TT,
                                          P=hp_comp[x,]$P,
                                          Ktrue=3,
                                          mu=1,
                                          phi=.8,
                                          rho=0,
                                          Pcat=NULL,
                                          pers=.95,
                                          pNAs=0,
                                          typeNA=3),
                                      mc.cores = parallel::detectCores())
spClust_no.miss_setup2 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_speclust(
                                                 seed=hp_comp[x,]$seed,
                                                 TT=hp_comp[x,]$TT,
                                                 P=hp_comp[x,]$P,
                                                 Ktrue=3,
                                                 mu=1,
                                                 phi=.8,
                                                 rho=0.2,
                                                 Pcat=NULL,
                                                 pers=.95,
                                                 pNAs=0,
                                                 typeNA=3),
                                             mc.cores = parallel::detectCores())
spClust_no.miss_setup3 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_speclust(
                                                 seed=hp_comp[x,]$seed,
                                                 TT=hp_comp[x,]$TT,
                                                 P=hp_comp[x,]$P,
                                                 Ktrue=3,
                                                 mu=.5,
                                                 phi=.8,
                                                 rho=0,
                                                 Pcat=NULL,
                                                 pers=.95,
                                                 pNAs=0,
                                                 typeNA=3),
                                             mc.cores = parallel::detectCores())
end_no.miss_specluster=Sys.time()
elapsed_no.miss_specluster=end_no.miss_specluster-start_no.miss_specluster
save(spClust_no.miss_setup1,
     spClust_no.miss_setup2,
     spClust_no.miss_setup3,
     elapsed_no.miss_specluster,
     file="spClust_no_miss.RData")


# Incomplete data -----------------------------


# MCAR
start_MCAR=Sys.time()
mixedJM_MCAR.miss10 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed_missmec(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=.1,typeNA=0),
                                          mc.cores = parallel::detectCores())

mixedJM_MCAR.miss20 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed_missmec(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=.2,typeNA=0),
                                          mc.cores = parallel::detectCores())
end_MCAR=Sys.time()
elapsed_MCAR=end_MCAR-start_MCAR
save(mixedJM_MCAR.miss10,
     mixedJM_MCAR.miss20,
     elapsed_MCAR,
     file="mixedJM_MCAR.Rdata")

# MAR

start_MAR=Sys.time()
mixedJM_MAR.miss10 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed_missmec(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=.1,typeNA=1),
                                          mc.cores = parallel::detectCores())

mixedJM_MAR.miss20 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed_missmec(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=0,
                                              Pcat=NULL,pers=.95,
                                              pNAs=.2,typeNA=1),
                                          mc.cores = parallel::detectCores())
end_MAR=Sys.time()
elapsed_MAR=end_MAR-start_MAR
save(mixedJM_MAR.miss10,
     mixedJM_MAR.miss20,
     elapsed_MAR,
     file="mixedJM_MAR.Rdata")

# MNAR

start_MNAR=Sys.time()
mixedJM_MNAR.miss10 <- parallel::mclapply(1:nrow(hp),
                                         function(x)
                                           simstud_JMmixed_missmec(
                                             seed=hp[x,]$seed,
                                             lambda=hp[x,]$lambda,
                                             TT=hp[x,]$TT,
                                             P=hp[x,]$P,
                                             Ktrue=3,mu=1,
                                             phi=.8,rho=0,
                                             Pcat=NULL,pers=.95,
                                             pNAs=.1,typeNA=2),
                                         mc.cores = parallel::detectCores())

mixedJM_MNAR.miss20 <- parallel::mclapply(1:nrow(hp),
                                         function(x)
                                           simstud_JMmixed_missmec(
                                             seed=hp[x,]$seed,
                                             lambda=hp[x,]$lambda,
                                             TT=hp[x,]$TT,
                                             P=hp[x,]$P,
                                             Ktrue=3,mu=1,
                                             phi=.8,rho=0,
                                             Pcat=NULL,pers=.95,
                                             pNAs=.2,typeNA=2),
                                         mc.cores = parallel::detectCores())
end_MNAR=Sys.time()
elapsed_MNAR=end_MNAR-start_MNAR
save(mixedJM_MNAR.miss10,
     mixedJM_MNAR.miss20,
     elapsed_MNAR,
     file="mixedJM_MNAR.Rdata")


# Results -----------------------------------------------------------------

res_eval=function(res_obj,hp,lambda0=F,ARI=T){
  
  if(ARI){
    res=data.frame(hp,
                   ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                   imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err))),
                   clust_pur=unlist(lapply(res_obj,function(x)x$clust_pur))
    )
    
    # maxres=res%>%group_by(TT,P)%>%summarise(maxARI=max(ARI,na.rm=T))
    
    if(lambda0){
      res=res[which(res$lambda==0),]
      res%>%group_by(TT,P)%>%summarise(avARI=median(ARI,na.rm=T),
                                       sdARI=sd(ARI,na.rm=T),
                                       avErr=mean(imput.err),
                                       avClustPur=mean(clust_pur))
    }
    
    else{
      avres=res%>%group_by(TT,P,lambda)%>%summarise(avARI=median(ARI,na.rm=T),
                                                    sdARI=sd(ARI,na.rm=T),
                                                    avErr=mean(imput.err),
                                                    avClustPur=mean(clust_pur))
      
      avres%>%group_by(TT,P)%>%summarise(avARI=mean(avARI),
                                         sdARI=min(sdARI),
                                         lambda=lambda[which.max(avARI)],
                                         avErr=avErr[which.max(avARI)],
                                         avClustPur=avClustPur[which.max(avClustPur)])
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

## No missings
# Setup 1

res_eval(mixedJM_no.miss_setup1,hp)
res_eval(mixedJM_no.miss_setup1,hp,lambda0=T)
res_eval(spClust_no.miss_setup1,hp,lambda0=F)

# Setup 2
res_eval(mixedJM_no.miss_setup2,hp)
res_eval(mixedJM_no.miss_setup2,hp,lambda0=T)
res_eval(spClust_no.miss_setup2,hp,lambda0=F)

# Setup 3
res_eval(mixedJM_no.miss_setup3,hp)
res_eval(mixedJM_no.miss_setup3,hp,lambda0=T)
res_eval(spClust_no.miss_setup3,hp,lambda0=F)


# Convergence analysis ----------------------------------------------------


