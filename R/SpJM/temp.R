
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
lambda=jump_penalty
grid_size =NULL
verbose=F
tol=NULL
spatial_penalty = 0.05
gamma=spatial_penalty
alpha=2
n_states=3

# n_cores=10
# parallel_ga = 3

st=Sys.time()
prova=cont_STJM(Y,K,D,
                jump_penalty,
                spatial_penalty,
                n_init=10,
                max_iter=10,tol=NULL,initial_states=NULL,
                n_cores=30,prll=T,parallel_ga = F
)
en=Sys.time()
en-st

true_SS=data.frame(result$S)
colnames(true_SS)=paste0(1:M)
res_SS_long=true_SS %>%
  mutate(t = 1:n()) %>%
  tidyr::pivot_longer(cols = -t, names_to = "m", values_to = "true") %>%
  mutate(m=as.integer(m))%>%
  arrange(t, m)

est_SS=apply(prova[,-c(1:2)],1,which.max)
res_SS_long$est=est_SS

adj.rand.index(res_SS_long$true,res_SS_long$est)

table(res_SS_long$true,res_SS_long$est)
