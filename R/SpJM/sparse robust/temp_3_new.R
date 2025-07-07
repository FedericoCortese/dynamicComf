source("Utils_sparse_robust_2.R")

zeta0=0.15
alpha=.1
K=3
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
                       Ktrue=3,
                       mu=3,
                       rho=0,
                       nu=10,
                       phi=.8,
                       pers=0.99)

Y=simDat$SimData
true_stat=simDat$mchain

# State 1: 1,2,3
# State 2: 1,3,4
# State 3: 1,4,5
rel_=list()
rel_[[1]]=c(1,2,3)
rel_[[2]]=c(1,3,4)
rel_[[3]]=c(1,4,5)

inv_rel_=invert_rel(rel_,P)

irrelevant_features=which(sapply(inv_rel_,function(x)length(x)==0))
subsetY <- Y[,irrelevant_features]
Y[,irrelevant_features]=subsetY[sample(nrow(subsetY)), ]

relevant_features=which(sapply(inv_rel_,function(x)length(x)!=0))

for(p in 1:relevant_features){
  relevant_states=inv_rel_[[p]]
  indx=which(true_stat%in%relevant_states)
  Ytemp=Y[-indx,p]
  if(length(Ytemp)!=0){
    Y[-indx,p] = Ytemp[sample(length(Ytemp))]
  }
}

# Introduce outliers
set.seed(1)
out_sigma=100
p_out=0.02
N_out=TT*p_out
t_out=sample(1:TT,size=N_out)
Y[t_out,]=Y[t_out,]+rnorm(N_out*P,0,out_sigma)

truth=simDat$mchain
truth[t_out]=0

x11()
par(mfrow=c(4,3))
for (i in 1:P) {
  plot(Y[, i], col=truth+1, pch=19,ylab=i
       ,ylim=c(-10,10)
  )
}


# State 1, 
indx=which(true_stat==1)
# original subset
subsetY <- Y[indx, -c(1,2,3)]
# permute its rows
Y[indx, -c(1,2,3)] <- subsetY[sample(nrow(subsetY)), ]



# # State 2, 
indx=which(true_stat==2)
# original subset
subsetY <- Y[indx, -c(1,2,3)]
# permute its rows
Y[indx, -c(1,3,4)] <- subsetY[sample(nrow(subsetY)), ]

## State 3 
indx=which(true_stat==3)
# original subset
subsetY <- Y[indx, -c(1,4,5)]
# permute its rows
Y[indx, -c(1,4,5)] <- subsetY[sample(nrow(subsetY)), ]


subsetY <- Y[,6:P]
Y[,6:P]=subsetY[sample(nrow(subsetY)), ]

# Introduce outliers
set.seed(1)
out_sigma=100
p_out=0.02
N_out=TT*p_out
t_out=sample(1:TT,size=N_out)
Y[t_out,]=Y[t_out,]+rnorm(N_out*P,0,out_sigma)

truth=simDat$mchain
truth[t_out]=0

x11()
par(mfrow=c(4,3))
for (i in 1:P) {
  plot(Y[, i], col=truth+1, pch=19,ylab=i
       ,ylim=c(-10,10)
       )
}

# Y$X1=factor(round(Y$X1))
# Y$X2=factor(round(Y$X2))


# zeta_grid=seq(0.1,.3,.1)
# lambda_grid=seq(0,.3,.1)
# K_grid=2:3
# B=2

# prv=RJM_COSA_gap(Y,
#                  zeta_grid,
#                  lambda_grid,
#                  K_grid,
#                  tol=NULL,n_outer=15,alpha=.1,verbose=F,n_cores=NULL,
#                  B, knn=10,c=2,M=NULL)
source("Utils_sparse_robust_2.R")
startR=Sys.time()
prv=robust_sparse_jump(Y=as.matrix(Y),
                       zeta0=.25,
                       lambda=.25,
                       K=3,
                       tol        = 1e-16,
                       n_init     = 3,
                       n_outer    = 20,
                       alpha      = 0.1,
                       verbose    = T,
                       knn        = 10,
                       c          = 10,
                       M          = NULL)
endR=Sys.time()
print(endR-startR)

round(prv$W,2)
prv$s[prv$v<0.5]=0
table(prv$s,truth)

# cpp NOT RELIABLE
Rcpp::sourceCpp("rob_JM.cpp")
library(cluster)
startCpp=Sys.time()
prv_cpp=robust_JM_COSA(as.matrix(Y),
                       zeta0=.2,
                       lambda=.2,
                       K=2,
                       tol=1e-16,
                       n_init =10,
                       n_outer   = 20,
                       alpha  = 0.1,
                       verbose  = TRUE,
                       knn       = 10,
                       c      = 2.0
) 
endCpp=Sys.time()
print(endCpp-startCpp)

round(prv_cpp$W,2)

prv_cpp$s[prv_cpp$v<0.5]=0

table(prv_cpp$s,truth)

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