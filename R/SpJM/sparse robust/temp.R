source("Utils_sparse_robust_2.R")

zeta0=0.15
alpha=.1
K=2
tol=1e-16
n_outer=15
verbose=T
lambda=.25

TT=1000
P=10

simDat=sim_data_stud_t(seed=123,
                       TT=TT,
                       P=P,
                       Pcat=NULL,
                       Ktrue=2,
                       mu=1,
                       rho=0,
                       nu=100,
                       phi=.8,
                       pers=0.95)

Y=simDat$SimData
true_stat=simDat$mchain

plot(Y[,2],col=true_stat,pch=19)




nu=4
# State 1, only features 1,2 and 3 are relevant
indx=which(true_stat!=1)
Sigma <- matrix(0,ncol=P-3,nrow=P-3)
diag(Sigma)=5
Y[indx,-(1:3)]=mvtnorm::rmvt(length(indx),
                             sigma = (nu-2)*Sigma/nu,
                             df = nu, delta = rep(0,P-3))

# # State 2, only features 3,4 and 5 are relevant
indx=which(true_stat!=2)
Y[indx,-(3:5)]=mvtnorm::rmvt(length(indx),
                             sigma = (nu-2)*Sigma/nu,
                             df = nu, delta = rep(0,P-3))


Sigma <- matrix(0,ncol=P-5,nrow=P-5)
diag(Sigma)=5
Y[,6:P]=mvtnorm::rmvt(TT,
                      sigma = (nu-2)*Sigma/nu,
                      df = nu, delta = rep(0,P-5))

# Introduce outliers
set.seed(1)
N_out=TT*0.02
t_out=sample(1:TT,size=N_out)
Y[t_out,]=Y[t_out,]+rnorm(N_out*P,0,10)

truth=simDat$mchain
truth[t_out]=0

x11()
par(mfrow=c(4,3))
for (i in 1:P) {
  plot(Y[, i], col=truth+1, pch=19,ylab=i)
}

# Y$X1=factor(round(Y$X1))
# Y$X2=factor(round(Y$X2))


zeta_grid=seq(0.1,.3,.1)
lambda_grid=seq(0,.3,.1)
K_grid=2:3
B=2

prv=RJM_COSA_gap(Y,
                 zeta_grid,
                 lambda_grid,
                 K_grid,
                 tol=NULL,n_outer=15,alpha=.1,verbose=F,n_cores=NULL,
                 B, knn=10,c=2,M=NULL)

library(ggplot2)

x11()
ggplot(prv$gap_stats, 
       aes(x = zeta0, 
           y = GAP, 
           color = factor(K), 
           group = factor(K))) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set1", name = "K") +
  theme_minimal() +
  labs(
    x = expression(zeta[0]),
    y = "GAP",
    title = "Gap Statistic vs. Zeta0, by Number of Clusters"
  )