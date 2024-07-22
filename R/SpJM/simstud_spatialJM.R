source("Utils.R")

gamma=seq(0,1,by=.05)
M=c(100,400,900)
P=c(25,50,75)
seeds=1:100

hp=expand.grid(M=M,P=P,gamma=gamma,seed=seeds)


# No missing rho=0.5 --------------------------------------------------------------

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

library(dplyr)



res_eval=function(res_obj,hp,gamma0=F,ARI=T){
  library(dplyr)
  
  if(ARI){
    res=data.frame(hp,ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                   imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err)))
    )
    
    # maxres=res%>%group_by(TT,P)%>%summarise(maxARI=max(ARI,na.rm=T))
    
    if(gamma0){
      res=res[which(res$gamma==0),]
      res%>%group_by(M,P)%>%summarise(avARI=median(ARI,na.rm=T),
                                       avErr=mean(imput.err))
    }
    
    else{
      avres=res%>%group_by(M,P,gamma)%>%summarise(avARI=median(ARI,na.rm=T),
                                                    avErr=mean(imput.err))
      
      avres%>%group_by(M,P)%>%summarise(maxARI=max(avARI),
                                         gamma=gamma[which.max(avARI)],
                                         avErr=avErr[which.max(avARI)])
    }
  }
  else{
    res=data.frame(hp,
                   #ARI=unlist(lapply(res_obj,function(x)x$ARI)),
                   imput.err=unlist(lapply(res_obj,function(x)mean(x$imput.err)))
    )
    res%>%group_by(M,P)%>%summarise(
      #avARI=median(ARI,na.rm=T),
      avErr=mean(imput.err))
  }
  
}

res_spJM_rho05=res_eval(spatialJMrho05_no.miss,hp)

res_spJM_rho05_gam0=res_eval(spatialJMrho05_no.miss,hp,gamma0=T)
res_spJM_rho05_gam0

accur=function(x){
  true=order_states_freq(x$true_seq)
  est=order_states_freq(x$est_seq)
  true=factor(true)
  est=factor(est)
  
  tmp=confusionMatrix(true,est)$overall[1]
  return(tmp)
}
library(caret)
acc=unlist(lapply(spatialJMrho05_no.miss,accur))

res=data.frame(hp,acc=acc)
avres=res%>%group_by(M,P,gamma)%>%summarise(avAcc=median(acc,na.rm=T))
avres%>%group_by(M,P)%>%summarise(maxAcc=max(avAcc),gamma=gamma[which.max(avAcc)])

res0=res[which(res$gamma==0),]
res0%>%group_by(M,P)%>%summarise(avAcc=median(acc,na.rm=T))



# No missing rho=0 --------------------------------------------------------------

start_no.miss=Sys.time()
spatialJMrho0_no.miss <- parallel::mclapply(1:nrow(hp),
                                             function(x)
                                               simstud_spatialJM(Ktrue=3,
                                                                 seed=hp[x,]$seed,
                                                                 gamma=hp[x,]$gamma,
                                                                 M=hp[x,]$M,
                                                                 P=hp[x,]$P,
                                                                 #Ktrue=3,
                                                                 mu=3,
                                                                 phi=.8,
                                                                 rho=0,
                                                                 Pcat=NULL,
                                                                 pNAs=0),
                                             mc.cores = parallel::detectCores()-1)

end_no.miss=Sys.time()
elapsed_no.miss=end_no.miss-start_no.miss
save(spatialJMrho0_no.miss,elapsed_no.miss,file="spatialJMrho0_no_miss.RData")