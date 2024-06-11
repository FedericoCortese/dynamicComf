source("Utils.R")

TT=c(50,100,500)
P=c(25,50,75)
seeds=1:100
hp_comp=expand.grid(TT=TT,P=P,seed=seeds)

# Random Forest
pNAs=0.1
start_rand.miss10_rf=Sys.time()
mixedJM_rand.miss10_rf <- parallel::mclapply(1:nrow(hp_comp),
                                      function(x)
                                        simstud_missForest(
                                          seed=hp_comp[x,]$seed,
                                          TT=hp_comp[x,]$TT,
                                          P=hp_comp[x,]$P,
                                          Ktrue=3,mu=1,
                                          phi=.8,rho=0,
                                          Pcat=NULL,pers=.95,
                                          pNAs=pNAs,typeNA=0),
                                      mc.cores = parallel::detectCores()-1)
end_rand.miss10_rf=Sys.time()
elapsed_rand.miss10_rf=end_rand.miss10_rf-start_rand.miss10_rf
save(mixedJM_rand.miss10_rf,elapsed_rand.miss10_rf,file="mixedJM_rand_miss10_rf.Rdata")
rm(mixedJM_rand.miss10_rf,elapsed_rand.miss10_rf)

pNAs=0.2
start_rand.miss20_rf=Sys.time()
mixedJM_rand.miss20_rf <- parallel::mclapply(1:nrow(hp_comp),
                                        function(x)
                                          simstud_missForest(
                                            seed=hp_comp[x,]$seed,
                                            TT=hp_comp[x,]$TT,
                                            P=hp_comp[x,]$P,
                                            Ktrue=3,mu=1,
                                            phi=.8,rho=0,
                                            Pcat=NULL,pers=.95,
                                            pNAs=pNAs,typeNA=0),
                                        mc.cores = parallel::detectCores()-1)

end_rand.miss20_rf=Sys.time()
elapsed_rand.miss20_rf=end_rand.miss20_rf-start_rand.miss20_rf
save(mixedJM_rand.miss20_rf,elapsed_rand.miss20_rf,file="mixedJM_rand_miss20_rf.Rdata")
rm(mixedJM_rand.miss20_rf,elapsed_rand.miss20_rf)

# Cont miss

pNAs=0.1
start_cont.miss10_rf=Sys.time()
mixedJM_cont.miss10_rf <- parallel::mclapply(1:nrow(hp_comp),
                                             function(x)
                                               simstud_missForest(
                                                 seed=hp_comp[x,]$seed,
                                                 TT=hp_comp[x,]$TT,
                                                 P=hp_comp[x,]$P,
                                                 Ktrue=3,mu=1,
                                                 phi=.8,rho=0,
                                                 Pcat=NULL,pers=.95,
                                                 pNAs=pNAs,typeNA=1),
                                             mc.cores = parallel::detectCores()-1)

end_cont.miss10_rf=Sys.time()
elapsed_cont.miss10_rf=end_cont.miss10_rf-start_cont.miss10_rf
save(mixedJM_cont.miss10_rf,elapsed_cont.miss10_rf,file="mixedJM_cont_miss10_rf.Rdata")
rm(mixedJM_cont.miss10_rf,elapsed_cont.miss10_rf)

pNAs=0.2
start_cont.miss20_rf=Sys.time()
mixedJM_cont.miss20_rf <- parallel::mclapply(1:nrow(hp_comp),
                                             function(x)
                                               simstud_missForest(
                                                 seed=hp_comp[x,]$seed,
                                                 TT=hp_comp[x,]$TT,
                                                 P=hp_comp[x,]$P,
                                                 Ktrue=3,mu=1,
                                                 phi=.8,rho=0,
                                                 Pcat=NULL,pers=.95,
                                                 pNAs=pNAs,typeNA=1),
                                             mc.cores = parallel::detectCores()-1)

end_cont.miss20_rf=Sys.time()
elapsed_cont.miss20_rf=end_cont.miss20_rf-start_cont.miss20_rf
save(mixedJM_cont.miss20_rf,elapsed_cont.miss20_rf,file="mixedJM_cont_miss20_rf.Rdata")
rm(mixedJM_cont.miss20_rf,elapsed_cont.miss20_rf)