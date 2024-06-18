source("Utils.R")

gamma=seq(0,20,by=1)
M=c(100,400,900)
P=c(25,50,75)
seeds=1:100

hp=expand.grid(M=M,P=P,gamma=gamma,seed=seeds)


# No missing --------------------------------------------------------------

start_no.miss=Sys.time()
spatialJMrho05_no.miss <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        simstud_spatialJM(Ktrue=3,
                                                          seed=hp[x,]$seed,
                                                          gamma=hp[x,]$gamma,
                                                          M=hp[x,]$M,
                                                          P=hp[x,]$P,
                                                          #Ktrue=3,
                                                          mu=3,
                                                          phi=.8,
                                                          rho=0.5,
                                                          Pcat=NULL,
                                                          pNAs=0),
                                      mc.cores = parallel::detectCores()-1)

end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss
save(spatialJMrho05_no.miss,elapsed_no.miss,file="spatialJMrho05_no_miss.RData")
