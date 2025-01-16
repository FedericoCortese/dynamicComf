
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

