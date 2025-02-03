
source("Utils.R")

Sim=sim_data_mixed(seed=123,
                   TT=1000,
                   P=5,
                   Ktrue=3,
                   mu=1,
                   phi=.8,
                   rho=0,
                   Pcat=0,
                   pers=.95,
                   pNAs=0,
                   typeNA=3)

Y=Sim$SimData.complete
K=3
jump_penalty = .5
grid_size =.05
verbose=T
tol=NULL
n_init=3
max_iter=5
initial_states = NULL
alpha=2
mode_loss=T

fit=cont_jumpR(Y,
               K=3,
               jump_penalty = .5,
               grid_size = .05,
               verbose=T,
               tol=NULL,
               n_init=3,
               max_iter=5,
               initial_states = NULL)

adj.rand.index(apply(fit,1,which.max),Sim$mchain)


### spatiotemp

source("Utils.R")


M=10
TT=50
theta=.01
beta=.9
K=3
mu=1
rho=0
phi=0.8
P=6
Pcat=2
seed=123
pg=0
pNAs=0

result <- generate_spatio_temporal_data(M, TT, theta, beta, K = K,
                                        mu=mu,rho=rho,phi=phi,
                                        P=P,Pcat=Pcat,seed=seed,pGap=pg,pNAs=pNAs)

Y.compl=result$Y
D=result$dist_matrix
Y=result$Y.NA
#Y=Y[,-(3:4)]
head(Y)

jump_penalty = .1
grid_size =NULL
verbose=F
tol=NULL
spatial_penalty = .1
alpha=2
n_states=3

prova=cont_STJM(Y,K,D,
                jump_penalty,
                spatial_penalty,
                alpha=2,grid_size=NULL,mode_loss=T,
                rand_search_sample = 100,
                n_init=10,
                max_iter=10,tol=NULL,initial_states=NULL,
                n_cores=10,prll=T
)

