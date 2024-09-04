source("Utils.R")


# mu=3 --------------------------------------------------------------------

K=3
P=30
###
mu=3
rho=0.2
###
phi=.8
Pcat=10
pNAs=0

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
# M=c(25,100)
# TT=c(10,50)

# lambda=.05
# gamma=.05
# seed=1:2
M=100
TT=10
hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)



start_STJsim=Sys.time()
STJsim <- parallel::mclapply(1:nrow(hp),
                             function(x)
                               simstud_STJump(lambda=hp[x,]$lambda,
                                          gamma=hp[x,]$gamma,
                                          seed=hp[x,]$seed,
                                          M=hp[x,]$M,
                                          TT=hp[x,]$TT,
                                          mu=mu,rho=rho,
                                          K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs),
                             mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim,elapsed_STJsim,file="STJsim.RData")


library(dplyr)
res=data.frame(hp[,c("lambda","gamma","seed")],ARI=unlist(lapply(STJsim,function(x)x$ARI)))

res_av=res%>%group_by(lambda,gamma)%>%summarise(avARI=mean(ARI,na.rm=T))
