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
