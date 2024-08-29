source("Utils.R")

TT=4
Ktrue=4
seed=1
M=100
P=20
mu=3
phi=.8
rho=0.2
Pcat=0
pNAs=0
sp_indx=1:M
sp_indx=matrix(sp_indx,ncol=sqrt(M),byrow=T)

S_true=matrix(0,nrow=TT,ncol=M)

C=Cmatrix(sp_indx)
# simDat=sim_spatial_JM(P,C,seed,
#                       rho=rho,Pcat=Pcat, phi=phi,
#                       mu=mu,pNAs=pNAs)

Y=NULL



# Mixed Potts-Autoregressive ---------------------------------------------------------------------

# States at time 1
t=1
simDat=sim_spatiotemp_JM(P,C,t,
                      rho=rho,Pcat=Pcat, phi=phi,
                      mu=mu,pNAs=pNAs,ST=NULL,n_states=Ktrue)
temp=simDat$SimData
temp$m=1:M
temp$t=rep(t,M)
Y=rbind(Y,temp)
S_true[t,]=simDat$s
S_true[t,]=order_states_condMean(Y$V20[Y$t==t],S_true[t,])

# Temporal persistence 
PI=0.9
for(t in 2:TT){
  simDat=sim_spatiotemp_JM(P,C,t,
                        rho=rho,Pcat=Pcat, phi=phi,
                        mu=mu,pNAs=pNAs,ST=S_true[t-1,],PI=PI,n_states=Ktrue)
  temp=simDat$SimData
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=simDat$s
  S_true[t,]=order_states_condMean(Y$V20[Y$t==t],S_true[t,])
  
}

# Put t and m in front of all the others with dplyr
Y <- Y %>% select(t, m, everything())

data_matrix <- matrix(S_true[4,], nrow = sqrt(M), ncol = sqrt(M),byrow = T)
data_df <- as.data.frame(as.table(data_matrix))
colnames(data_df) <- c("X", "Y", "Value")
ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()
        # ,
        # legend.position = "none"
        ) +
  coord_fixed() 

lambda=0.8
gamma=0.8
initial_states=NULL
max_iter=10
n_init=10
tol=NULL
verbose=T

fit <- STjumpR(Y, n_states = 3, C, jump_penalty=lambda,spatial_penalty = gamma, verbose=T)

best_s=fit$best_s
for(t in 1:TT){
  best_s[t,]=order_states_condMean(Y$V20[Y$t==t],best_s[t,])
}

tt=2

data_matrix <- matrix(S_true[tt,], nrow = sqrt(M), ncol = sqrt(M),byrow = T)
data_df <- as.data.frame(as.table(data_matrix))
colnames(data_df) <- c("X", "Y", "Value")

Ptrue=ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()
          ,
          legend.position = "none"
    ) +
    coord_fixed() 
  
data_matrix <- matrix(best_s[tt,], nrow = sqrt(M), ncol = sqrt(M),byrow = T)
data_df <- as.data.frame(as.table(data_matrix))
colnames(data_df) <- c("X", "Y", "Value")
  
Pest=ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()
          ,
          legend.position = "none"
    ) +
    coord_fixed() 
  
library(ggpubr)
arr_plot=ggarrange(Ptrue, Pest)
annotate_figure(arr_plot, top = text_grob(paste0("t=",tt), 
                                            color = "black", face = "bold", size = 14))


# Potts only --------------------------------------------------------------

require("potts")

Y=NULL
M=dim(C)[1]

set.seed(2024)
ncolor = as.integer(3) # transizione di fase continua per 1 <= ncolor <= 4
nr = sqrt(M)
nc = sqrt(M)
init <- matrix(sample(ncolor, nr * nc, replace = TRUE), nrow = nr, ncol=nc)
init <- packPotts(init, ncol = ncolor)

beta <- log(1 + sqrt(ncolor))
theta <- c(rep(0, ncolor), beta)
out <- potts(init, param=theta, nbatch = 200 , blen=1, nspac=1)

#nit=10
sim= vector("list",TT)
sim[[1]] = out$final
rotate <- function(x) apply(t(x), 2, rev)
# Recover decoded matrix
mat=unpackPotts(sim[[1]])
mat=rotate(mat)
s=c(t(mat))

t=1
temp=sim_obs(s=s,mu=mu,rho=rho,P=P,Pcat=P/2,n_states=3,seed=seed,pNAs=pNAs)
temp=temp$SimData
temp$m=1:M
temp$t=rep(t,M)
Y=rbind(Y,temp)
S_true[t,]=s
S_true[t,]=order_states_condMean(Y$V20[Y$t==t],S_true[t,])

for(t in 2:TT){
  out = potts(sim[[t-1]], param=theta, nbatch = 1 , blen=1, nspac=1)
  sim[[t]] = out$final
  mat=unpackPotts(sim[[t]])
  mat=rotate(mat)
  s=c(t(mat))
  
  temp=sim_obs(s=s,mu=mu,rho=rho,P=P,Pcat=P/2,n_states=3,seed=seed,pNAs=pNAs)
  temp=temp$SimData
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=s
  S_true[t,]=order_states_condMean(Y$V20[Y$t==t],S_true[t,])
  
}

Y <- Y %>% select(t, m, everything())

data_matrix <- matrix(S_true[4,], nrow = sqrt(M), ncol = sqrt(M),byrow = T)
data_df <- as.data.frame(as.table(data_matrix))
colnames(data_df) <- c("X", "Y", "Value")
ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()
        # ,
        # legend.position = "none"
  ) +
  coord_fixed() 