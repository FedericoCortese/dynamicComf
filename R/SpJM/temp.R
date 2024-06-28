source("Utils.R")
YY=sim_data(seed=1,Ktrue=3,N=100,P=5,cors=c(0,.5,.9),pers=.95,m=2)

sqrt(5)
true_st=YY$true_states
YY=as.matrix(YY$YY)
n_states=as.integer(3)

est.nystrup=sparse_jump(YY,n_states,jump_penalty = 1,max_features = sqrt(5)/1.16 ,verbose = T)
est.mine=sparse_jumpR(YY,3,jump_penalty = 1,max_features = sqrt(5)/1.16 ,verbose = T)

round(est.nystrup[[2]],2)
round(est.mine$feat_w,2)

adj.rand.index(est.nystrup[[1]],true_st)  
adj.rand.index(est.mine$states,true_st)  


YY1=sim_data(seed=1,Ktrue=3,N=100,P=2,cors=c(0,.5,.9),pers=.95,m=0)
YY2=sim_data(seed=2,Ktrue=3,N=100,P=2,cors=c(0,.5,.9),pers=.95,m=0)
YY3=sim_data(seed=3,Ktrue=3,N=100,P=2,cors=c(0,.5,.9),pers=.95,m=0)

Z=array(0,dim=c(100,3,2))
Z[,1,]=YY1$YY
Z[,2,]=YY2$YY
Z[,3,]=YY3$YY

initial_states=cbind(YY1$true_states,YY2$true_states,YY3$true_states);initial_states

#####


## ELBOW FOR SELECTION OF LAMBDA
source("Utils.R")
YY=sim_data_mixed(seed=123,
                 TT=1000,
                 P=10,
                 Ktrue=3,
                 mu=1,
                 phi=.8,
                 rho=0,
                 Pcat=NULL,
                 pers=.95,
                 pNAs=0,
                 typeNA=2)

Y=YY$SimData.NA
#prv1=jump_mixed2(Y,3,jump_penalty = .1,timeflag = F)

lambda=seq(0,2,by=.05)
res=lapply(lambda,function(x){
  est=jump_mixed2(Y=Y,4,jump_penalty = x,timeflag = F)
  dev=get_BCD(Y,est$best_s)
  WCD=dev$WCD
  BCD=dev$BCD
  true_seq=order_states_freq(YY$mchain)
  est_seq=order_states_freq(est$best_s)
  true_seq=factor(true_seq,levels=c(1,2,3,4))
  est_seq=factor(est_seq,levels=c(1,2,3,4))
  BAC=bacc(true_seq,est_seq)
 
  return(list(WCD=sum(WCD),
              BCD=sum(BCD),
              K=length(unique(est$best_s)),
              BAC=BAC))
})

wcd_lam=data.frame(WCD=unlist(lapply(res,function(x){x$WCD})),
                   BCD=unlist(lapply(res,function(x){x$BCD})),
                   K=unlist(lapply(res,function(x){x$K})),
                   lambda,
                   BAC=unlist(lapply(res,function(x){x$BAC})))

wcd_lam

plot(wcd_lam$lambda,wcd_lam$K,type="l")
plot(wcd_lam$lambda,wcd_lam$BCD,type="l",ylab="BCD",xlab="lambda")

prv=jump_mixed2(Y=Y,3,jump_penalty = 1.25,timeflag = F)
adj.rand.index(prv$best_s,YY$mchain)

prv1=jump_mixed2(Y=Y,3,jump_penalty = 1,timeflag = F)
adj.rand.index(prv1$best_s,YY$mchain)

prv15=jump_mixed2(Y=Y,3,jump_penalty = 1.5,timeflag = F)
adj.rand.index(prv15$best_s,YY$mchain)

plot(wcd_lam$lam,wcd_lam$WCD,type="l")


states=prv$best_s
condMM=prv$condMM

prvemp=jump_mixed(enth_tab4,2,jump_penalty = .1)
prvemp$best_s


### spatial JM
source("utils.R")
Mtrue=100
sqMtrue=sqrt(Mtrue)
sp_indx=matrix(1:Mtrue,ncol=sqMtrue,byrow=T)
C=Cmatrix(sp_indx)
P=10
seed=1
spDat=sim_spatial_JM(P,C,seed=sample(1:10,1),rho=0.5,Pcat=NULL, phi=.8,mu=3,pNAs=.1)
Strue=matrix(order_states_freq(spDat$s),ncol=sqMtrue,byrow=T)
im=packPotts(Strue,ncolor=3)
windows()
image(im,main="True states")
Y=spDat$SimData
Y=spDat$SimData.NA

lambdas=seq(11,20,by=1)
ARIs=rep(0,length(lambdas))
Ss=list()

for(i in 1:length(lambdas)){
  prv=spatial_jump(Y,C,3,jump_penalty = lambdas[i])
  Ss[[i]]=prv$best_s
  ARIs[i]=adj.rand.index(prv$best_s,spDat$s)
  print(i)
}

windows()
plot(lambdas,ARIs,type="l")
ARIs

est_states=matrix(order_states_freq(Ss[[1]]),ncol=sqMtrue,byrow=T)
im_est=packPotts(est_states,ncolor=3)
windows()
image(im_est,main="Estimated states")
