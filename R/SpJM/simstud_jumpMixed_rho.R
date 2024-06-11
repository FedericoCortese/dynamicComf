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
mixedJM_rho05_no.miss <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_JMmixed(
                                          seed=hp[x,]$seed,
                                          lambda=hp[x,]$lambda,
                                          TT=hp[x,]$TT,
                                          P=hp[x,]$P,
                                          Ktrue=3,mu=1,
                                          phi=.8,rho=.5,
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
# mixedJM_rho05_no.miss <- clusterApply(cl, 
#                                 1:nrow(hp), 
#                                 function(x)
#                                   simstud_JMmixed(
#                                     seed=hp[x,]$seed,
#                                     lambda=hp[x,]$lambda,
#                                     TT=hp[x,]$TT,
#                                     P=hp[x,]$P,
#                                     Ktrue=3,mu=1,
#                                     phi=.8,rho=.5,
#                                     Pcat=NULL,pers=.95,
#                                     pNAs=0,typeNA=2)
# )
# stopCluster(cl)

end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss
save(mixedJM_rho05_no.miss,elapsed_no.miss,file="mixedJM_rho05_no_miss.RData")

# Random missing ----------------------------------------------------------

pNAs=0.1
start_rand.miss10=Sys.time()
mixedJM_rho05_rand.miss10 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=.5,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=0),
                                          mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores()-1,type="SOCK")
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
# mixedJM_rho05_rand.miss10 <- clusterApply(cl,
#                                     1:nrow(hp),
#                                     function(x)
#                                       simstud_JMmixed(
#                                         seed=hp[x,]$seed,
#                                         lambda=hp[x,]$lambda,
#                                         TT=hp[x,]$TT,
#                                         P=hp[x,]$P,
#                                         Ktrue=3,mu=1,
#                                         phi=.8,rho=.5,
#                                         Pcat=NULL,pers=.95,
#                                         pNAs=pNAs,typeNA=0)
#                                     
# )
# stopCluster(cl)

end_rand.miss10=Sys.time()
elapsed_rand.miss10=end_rand.miss10-start_rand.miss10
save(mixedJM_rho05_rand.miss10,elapsed_rand.miss10,file="mixedJM_rho05_rand_miss10.Rdata")

pNAs=0.20
start_rand.miss20=Sys.time()
mixedJM_rho05_rand.miss20 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=.5,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=0),
                                          mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores()-1,type="SOCK")
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
# mixedJM_rho05_rand.miss20 <- clusterApply(cl,
#                                     1:nrow(hp),
#                                     function(x)
#                                       simstud_JMmixed(
#                                         seed=hp[x,]$seed,
#                                         lambda=hp[x,]$lambda,
#                                         TT=hp[x,]$TT,
#                                         P=hp[x,]$P,
#                                         Ktrue=3,mu=1,
#                                         phi=.8,rho=.5,
#                                         Pcat=NULL,pers=.95,
#                                         pNAs=pNAs,typeNA=0)
#                                     
# )
# stopCluster(cl)
end_rand.miss20=Sys.time()
elapsed_rand.miss20=end_rand.miss20-start_rand.miss20
save(mixedJM_rho05_rand.miss20,elapsed_rand.miss20,file="mixedJM_rho05_rand_miss20.Rdata")

# pNAs=0.50
# start_rand.miss50=Sys.time()
# mixedJM_rho05_rand.miss50 <- parallel::mclapply(1:nrow(hp),
#                                           function(x)
#                                             simstud_JMmixed(
#                                               seed=hp[x,]$seed,
#                                               lambda=hp[x,]$lambda,
#                                               TT=hp[x,]$TT,
#                                               P=hp[x,]$P,
#                                               Ktrue=3,mu=1,
#                                               phi=.8,rho=.5,
#                                               Pcat=NULL,pers=.95,
#                                               pNAs=pNAs,typeNA=0),
#                                           mc.cores = parallel::detectCores())
# end_rand.miss50=Sys.time()
# elapsed_rand.miss50=end_rand.miss50-start_rand.miss50


# Continuous missing ------------------------------------------------------

pNAs=0.1
start_cont.miss10=Sys.time()
mixedJM_rho05_cont.miss10 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=.5,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=1),
                                          mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores()-1,type="SOCK")
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
# mixedJM_rho05_cont.miss10 <- clusterApply(cl,
#                                     1:nrow(hp),
#                                     function(x)
#                                       simstud_JMmixed(
#                                         seed=hp[x,]$seed,
#                                         lambda=hp[x,]$lambda,
#                                         TT=hp[x,]$TT,
#                                         P=hp[x,]$P,
#                                         Ktrue=3,mu=1,
#                                         phi=.8,rho=.5,
#                                         Pcat=NULL,pers=.95,
#                                         pNAs=pNAs,typeNA=1)
#                                     
# )
# stopCluster(cl)
end_cont.miss10=Sys.time()
elapsed_cont.miss10=end_cont.miss10-start_cont.miss10
save(mixedJM_rho05_cont.miss10,elapsed_cont.miss10,file="mixedJM_rho05_cont_miss10.Rdata")

rm(mixedJM_rho05_cont.miss10)

pNAs=0.20
start_cont.miss20=Sys.time()
mixedJM_rho05_cont.miss20 <- parallel::mclapply(1:nrow(hp),
                                          function(x)
                                            simstud_JMmixed(
                                              seed=hp[x,]$seed,
                                              lambda=hp[x,]$lambda,
                                              TT=hp[x,]$TT,
                                              P=hp[x,]$P,
                                              Ktrue=3,mu=1,
                                              phi=.8,rho=.5,
                                              Pcat=NULL,pers=.95,
                                              pNAs=pNAs,typeNA=1),
                                          mc.cores = parallel::detectCores()-1)
# cl<-makeCluster(parallel::detectCores()-1,type="SOCK")
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
# mixedJM_rho05_cont.miss20 <- clusterApply(cl,
#                                     1:nrow(hp),
#                                     function(x)
#                                       simstud_JMmixed(
#                                         seed=hp[x,]$seed,
#                                         lambda=hp[x,]$lambda,
#                                         TT=hp[x,]$TT,
#                                         P=hp[x,]$P,
#                                         Ktrue=3,mu=1,
#                                         phi=.8,rho=.5,
#                                         Pcat=NULL,pers=.95,
#                                         pNAs=pNAs,typeNA=1)
#                                     
# )
# stopCluster(cl)
end_cont.miss20=Sys.time()
elapsed_cont.miss20=end_cont.miss20-start_cont.miss20
save(mixedJM_rho05_cont.miss20,elapsed_cont.miss20,file="mixedJM_rho05_cont_miss20.Rdata")

rm(mixedJM_rho05_cont.miss20)

# pNAs=0.50
# start_cont.miss50=Sys.time()
# mixedJM_rho05_cont.miss50 <- parallel::mclapply(1:nrow(hp),
#                                           function(x)
#                                             simstud_JMmixed(
#                                               seed=hp[x,]$seed,
#                                               lambda=hp[x,]$lambda,
#                                               TT=hp[x,]$TT,
#                                               P=hp[x,]$P,
#                                               Ktrue=3,mu=1,
#                                               phi=.8,rho=.5,
#                                               Pcat=NULL,pers=.95,
#                                               pNAs=pNAs,typeNA=1),
#                                           mc.cores = parallel::detectCores())
# end_cont.miss50=Sys.time()
# elapsed_cont.miss50=end_cont.miss50-start_cont.miss50

# Classification performance comparison --------------------------------------------------------------

hp_comp=expand.grid(TT=TT,P=P,seed=seeds)

# KMeMo
start_no.miss_kMeMo=Sys.time()
mixedJM_rho05_no.miss <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_JMmixed(
                                          seed=hp[x,]$seed,
                                          lambda=hp[x,]$lambda,
                                          TT=hp[x,]$TT,
                                          P=hp[x,]$P,
                                          Ktrue=3,mu=1,
                                          phi=.8,rho=.5,
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
# mixedJM_rho05_no.miss_kMeMo <- clusterApply(cl,
#                                       1:nrow(hp_comp),
#                                       function(x)
#                                         simstud_JMmixed(
#                                           seed=hp_comp[x,]$seed,
#                                           lambda=0,
#                                           TT=hp_comp[x,]$TT,
#                                           P=hp_comp[x,]$P,
#                                           Ktrue=3,
#                                           mu=1,
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
save(mixedJM_rho05_no.miss_kMeMo,elapsed_no.miss_kMeMo,file="mixedJM_rho05_no_miss_kMeMo.RData")
rm(mixedJM_rho05_no.miss_kMeMo,elapsed_no.miss_kMeMo)

# Spectral Clustering
library(SpectralClMixed)
start_no.miss_specluster=Sys.time()
mixedJM_rho05_no.miss <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_JMmixed(
                                          seed=hp[x,]$seed,
                                          lambda=hp[x,]$lambda,
                                          TT=hp[x,]$TT,
                                          P=hp[x,]$P,
                                          Ktrue=3,mu=1,
                                          phi=.8,rho=.5,
                                          Pcat=NULL,pers=.95,
                                          pNAs=0,typeNA=2),
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
# mixedJM_rho05_no.miss_specluster <- clusterApply(cl,
#                                            1:nrow(hp_comp),
#                                            function(x)
#                                              simstud_speclust(
#                                                seed=hp_comp[x,]$seed,
#                                                TT=hp_comp[x,]$TT,
#                                                P=hp_comp[x,]$P,
#                                                Ktrue=3,
#                                                mu=1,
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
save(mixedJM_rho05_no.miss_specluster,elapsed_no.miss_specluster,file="mixedJM_rho05_no_miss_specluster.RData")
rm(mixedJM_rho05_no.miss_specluster,elapsed_no.miss_specluster)


