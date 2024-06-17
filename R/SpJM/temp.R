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

### spatial JM
source("utils.R")
Mtrue=900
sqMtrue=sqrt(Mtrue)
sp_indx=matrix(1:Mtrue,ncol=sqMtrue,byrow=T)
C=Cmatrix(sp_indx)
P=50
seed=1
spDat=sim_spatial_JM(P,C,seed=sample(1:10,1),pers_fact=0.05,rho=0.5,Pcat=NULL, phi=.8,mu=3)
Strue=matrix(order_states_freq(spDat$s),ncol=sqMtrue,byrow=T)
im=packPotts(Strue,ncolor=3)
windows()
image(im,main="True states")
Y=spDat$Y

lambdas=seq(0,10,by=1)
ARIs=rep(0,length(lambdas))
Ss=list()

for(i in 1:length(lambdas)){
  prv=spatial_jump(Y,C,3,jump_penalty = lambdas[i])
  Ss[[i]]=prv$best_s
  ARIs[i]=adj.rand.index(prv$best_s,spDat$s)
  print(i)
}

plot(lambdas,ARIs,type="l")
ARIs

est_states=matrix(order_states_freq(Ss[[11]]),ncol=sqMtrue,byrow=T)
im_est=packPotts(est_states,ncolor=3)
windows()
image(im_est,main="Estimated states")
