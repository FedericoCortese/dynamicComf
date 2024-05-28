
# Simulation --------------------------------------------------------------

source("Utils.R")

# Ktrue=2 -----------------------------------------------------------------

pers1=.99
Ns=c(300,600,100)
seed=123
corsK2=c(.8,.4)
# YYs_K2_300=lapply(seed,sim_data,Ktrue=2,N=Ns[1],P=100,cors=corsK2,pers=pers1,m=2)
# YYs_K2_600=lapply(seed,sim_data,Ktrue=2,N=Ns[2],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_1000=lapply(seed,sim_data,Ktrue=2,N=Ns[3],P=2,cors=corsK2,pers=pers1,m=1)

Ytrue=YYs_K2_1000[[1]]$YY
dim(Ytrue)
true.st=YYs_K2_1000[[1]]$true_states
# 
# initial_states=NULL
# max_iter=10
# n_init=10
# tol=NULL
# verbose=FALSE

n_states=2

# Add cat vars
Ytrue[,1]=sample(1:2,Ns[3],replace = T)
Ytrue[,2]=sample(1:4,Ns[3],replace = T)

# 20% NAs
Y=Ytrue
set.seed(12345)
Y[sample(1:Ns[3],Ns[3]*.1),]=NA
head(Y)
Y=data.frame(Y)
Y$X1=factor(Y$X1)
Y$X2=factor(Y$X2)


str(Y)
Amelia::missmap(data.frame(Y))

jump_penalty=0.1

n_obs <- nrow(Y)

set.seed(12345)
#initial_states=sample(1:n_states,n_obs,replace = T)
prv=jump_mixed(Y,n_states,jump_penalty = jump_penalty,verbose = T)

table(prv$best_s)
prv$best_s

dfYs=data.frame(prv$Y,s=prv$best_s)
tapply(dfYs$X3,dfYs$s,mean)
tapply(dfYs$X1,dfYs$s,Mode)
tapply(dfYs$X2,dfYs$s,Mode)



# j_eucl=jumpR(Y,n_states,jump_penalty = 50,verbose = F)
# adj.rand.index(true.st,j_eucl)
# j_manh=jumpR(Y,n_states,jump_penalty = 50,verbose = F,method="manhattan")
# adj.rand.index(true.st,j_manh)

table(true.st)
table(prv)
#prv=recode(prv,"2"=1,"1"=2)

adj.rand.index(true.st,prv)

plot(true.st,type='l')
lines(prv,col='red')

# Linreg ------------------------------------------------------------------


Y=MClinregSim(TT,coeff,X,pers,init,family="gaussian")
mchainbin=Y$mc
Y=Y$SimData
#prvLRbin=jumpLR(Ybin,X,n_states,jump_penalty = 5,verbose = F,family="binomial")
round(prvLRbin$coefs,3)
coeff
table(prvLRbin$s)
table(mchainbin)

##### jump LR

source("Utils.R")
TT=1000
X=matrix(0,ncol=2,nrow=TT)
X[,1]=rnorm(TT,sd=2)
X[,2]=rnorm(TT,mean=2)
X=apply(X,2,scale)
coeff=matrix(0,nrow=2,ncol=3)
coeff[1,]=c(0,2,2)
coeff[2,]=c(2,-2,0.5)
pers=.95
init=c(.8,.2)


Y=MClinregSim(TT,coeff,X,pers,init,family="gaussian")
mchain=Y$mc
Y=Y$SimData
plot(X[,1],Y)
plot(X[,2],Y)
plot(mchain,type='l')

Y=as.matrix(Y)
n_states=as.integer(2)

prvLR=jumpLR(Y,X,n_states,jump_penalty = 1,verbose = F,family="gaussian")
round(prvLR$coefs,3)
coeff
table(prvLR$s)
table(mchain)

Ypois=MClinregSim(TT,coeff,X,pers,init,family="poisson")
mchain=Ypois$mc
Y=Ypois$SimData
prvLRpois=jumpLR(Y,X,n_states,jump_penalty = 5,verbose = F,family="poisson")
round(prvLRpois$coefs,3)
coeff
table(prvLRpois$s)
table(mchainpois)

Ybin=MClinregSim(TT,coeff,X,pers,init,family="binomial")
mchainbin=Ybin$mc
Y=Ybin$SimData
prvLRbin=jumpLR(Y,X,n_states,jump_penalty = .75,verbose = F,family="binomial")
round(prvLRbin$coefs,3)
coeff
table(prvLRbin$s)
table(mchainbin)
