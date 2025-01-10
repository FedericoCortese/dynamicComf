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
lambda=seq(0,1,by=.05)
TT=c(50,100,500)
P=c(4,20,50)
seeds=1:100
hp_spCl=expand.grid(TT=TT,P=P,seed=seeds)

library(SpectralClMixed)

start_no.miss_specluster=Sys.time()

# Setup 1
cl<-makeCluster(parallel::detectCores(),type="SOCK")
parallel::clusterExport(cl,ls())
parallel::clusterEvalQ(cl,{
  library(RcppHMM)
  library(reticulate)
  library(pdfCluster)
  library(boot)
  library(xtable)
  library(dplyr)
  library(cluster)
  library(gower)
  library(StatMatch)
  library(SpectralClMixed)
})
spClust_no.miss_setup1 <- clusterApply(cl,
                                       1:nrow(hp_comp),
                                       function(x)
                                         simstud_speclust(
                                           seed=hp_spCl[x,]$seed,
                                           TT=hp_spCl[x,]$TT,
                                           P=hp_spCl[x,]$P,
                                           Ktrue=3,
                                           mu=1,
                                           phi=.8,
                                           rho=0,
                                           Pcat=NULL,
                                           pers=.95,
                                           pNAs=0,
                                           typeNA=3)
)
stopCluster(cl)
save(spClust_no.miss_setup1,
     file="spCluster_no_miss_setup1.RData")
rm(spClust_no.miss_setup1)

# Setup 2
cl<-makeCluster(parallel::detectCores(),type="SOCK")
parallel::clusterExport(cl,ls())
parallel::clusterEvalQ(cl,{
  library(RcppHMM)
  library(reticulate)
  library(pdfCluster)
  library(boot)
  library(xtable)
  library(dplyr)
  library(cluster)
  library(gower)
  library(StatMatch)
  library(SpectralClMixed)
})
spClust_no.miss_setup2 <- clusterApply(cl,
                                       1:nrow(hp_comp),
                                       function(x)
                                         simstud_speclust(
                                           seed=hp_spCl[x,]$seed,
                                           TT=hp_spCl[x,]$TT,
                                           P=hp_spCl[x,]$P,
                                           Ktrue=3,
                                           mu=1,
                                           phi=.8,
                                           rho=0.2,
                                           Pcat=NULL,
                                           pers=.95,
                                           pNAs=0,
                                           typeNA=3)
)
stopCluster(cl)
save(spClust_no.miss_setup2,
     file="spCluster_no_miss_setup2.RData")
rm(spClust_no.miss_setup2)

cl<-makeCluster(parallel::detectCores(),type="SOCK")
parallel::clusterExport(cl,ls())
parallel::clusterEvalQ(cl,{
  library(RcppHMM)
  library(reticulate)
  library(pdfCluster)
  library(boot)
  library(xtable)
  library(dplyr)
  library(cluster)
  library(gower)
  library(StatMatch)
  library(SpectralClMixed)
})
spClust_no.miss_setup3 <- clusterApply(cl,
                                       1:nrow(hp_comp),
                                       function(x)
                                         simstud_speclust(
                                           seed=hp_spCl[x,]$seed,
                                           TT=hp_spCl[x,]$TT,
                                           P=hp_spCl[x,]$P,
                                           Ktrue=3,
                                           mu=.5,
                                           phi=.8,
                                           rho=0,
                                           Pcat=NULL,
                                           pers=.95,
                                           pNAs=0,
                                           typeNA=3)
)
stopCluster(cl)
save(spClust_no.miss_setup3,
     file="spCluster_no_miss_setup3.RData")
rm(spClust_no.miss_setup3)

end_no.miss_specluster=Sys.time()
elapsed_no.miss_specluster=end_no.miss_specluster-start_no.miss_specluster


# start_no.miss_specluster=Sys.time()
# spClust_no.miss_setup1 <- parallel::mclapply(1:nrow(hp_spCl),
#                                       function(x)
#                                         simstud_speclust(
#                                           seed=hp_spCl[x,]$seed,
#                                           TT=hp_spCl[x,]$TT,
#                                           P=hp_spCl[x,]$P,
#                                           Ktrue=3,
#                                           mu=1,
#                                           phi=.8,
#                                           rho=0,
#                                           Pcat=NULL,
#                                           pers=.95,
#                                           pNAs=0,
#                                           typeNA=3),
#                                       mc.cores = parallel::detectCores())
# spClust_no.miss_setup2 <- parallel::mclapply(1:nrow(hp_spCl),
#                                              function(x)
#                                                simstud_speclust(
#                                                  seed=hp_spCl[x,]$seed,
#                                                  TT=hp_spCl[x,]$TT,
#                                                  P=hp_spCl[x,]$P,
#                                                  Ktrue=3,
#                                                  mu=1,
#                                                  phi=.8,
#                                                  rho=0.2,
#                                                  Pcat=NULL,
#                                                  pers=.95,
#                                                  pNAs=0,
#                                                  typeNA=3),
#                                              mc.cores = parallel::detectCores())
# spClust_no.miss_setup3 <- parallel::mclapply(1:nrow(hp_spCl),
#                                              function(x)
#                                                simstud_speclust(
#                                                  seed=hp_spCl[x,]$seed,
#                                                  TT=hp_spCl[x,]$TT,
#                                                  P=hp_spCl[x,]$P,
#                                                  Ktrue=3,
#                                                  mu=.5,
#                                                  phi=.8,
#                                                  rho=0,
#                                                  Pcat=NULL,
#                                                  pers=.95,
#                                                  pNAs=0,
#                                                  typeNA=3),
#                                              mc.cores = parallel::detectCores())
# end_no.miss_specluster=Sys.time()
# elapsed_no.miss_specluster=end_no.miss_specluster-start_no.miss_specluster
# save(spClust_no.miss_setup1,
#      spClust_no.miss_setup2,
#      spClust_no.miss_setup3,
#      elapsed_no.miss_specluster,
#      file="spClust_no_miss.RData")


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


## MCAR

#10%
res_eval(mixedJM_MCAR.miss10,hp)
res_eval(mixedJM_MCAR.miss10,hp,lambda0=T)

#20%
res_eval(mixedJM_MCAR.miss20,hp)
res_eval(mixedJM_MCAR.miss20,hp,lambda0=T)

## MAR

#10%
res_eval(mixedJM_MAR.miss10,hp)
res_eval(mixedJM_MAR.miss10,hp,lambda0=T)

#20%
res_eval(mixedJM_MAR.miss20,hp)
res_eval(mixedJM_MAR.miss20,hp,lambda0=T)

## MNAR

#10%
res_eval(mixedJM_MNAR.miss10,hp)
res_eval(mixedJM_MNAR.miss10,hp,lambda0=T)

#20%
res_eval(mixedJM_MNAR.miss20,hp)
res_eval(mixedJM_MNAR.miss20,hp,lambda0=T)




# l1 vs l2 ----------------------------------------------------------------

lambda=seq(0,1,by=.05)
TT=c(50,100,500)
P=c(4,20,50)
seeds=1:100
hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

# Setup 1, no missing data

mixedJM_l1vsl2 <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_JMmixed_l1vsl2(seed=hp[x,]$seed,
                                                                      TT=hp[x,]$TT,
                                                                      P=hp[x,]$P,
                                                                      lambda=hp[x,]$lambda,
                                                                      Ktrue=3,
                                                                      mu=1,
                                                                      rho=0,
                                                                      nu=4,
                                                                      pers=.95),
                                             mc.cores = parallel::detectCores())


save(mixedJM_l1vsl2,
     file="mixedJM_l1vsl2.Rdata")

res=data.frame(hp,
               ARI_l1=unlist(lapply(mixedJM_l1vsl2,function(x)x$ARI_l1)),
               ARI_l2=unlist(lapply(mixedJM_l1vsl2,function(x)x$ARI_l2))
)
    
    
avg_res <- res %>%
  group_by(TT, P, lambda) %>%
  summarise(
    avg_ARI_l1 = mean(ARI_l1,na.rm=T),
    sd_ARI_l1 = sd(ARI_l1,na.rm=T),
    avg_ARI_l2 = mean(ARI_l2,na.rm=T),
    sd_ARI_l2 = sd(ARI_l2,na.rm=T),
    .groups = "drop"
  )

max_ARI_l1 <- avg_res %>%
  group_by(TT, P) %>%
  filter(avg_ARI_l1 == max(avg_ARI_l1)) %>%
  ungroup()

# Find the maximum average ARI_l2 and corresponding lambda and standard deviation for each TT and P
max_ARI_l2 <- avg_res %>%
  group_by(TT, P) %>%
  filter(avg_ARI_l2 == max(avg_ARI_l2)) %>%
  ungroup()


print(max_ARI_l1)
print(max_ARI_l2)
  
# merged results
merged_res <- merge(max_ARI_l1[,c(1,2,4,5)], max_ARI_l2[,c(1,2,6,7)], 
                    by = c("TT", "P"))

merged_res

