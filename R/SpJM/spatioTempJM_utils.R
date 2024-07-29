source("Utils.R")

TT=4
Ktrue=3
seed=1
M=100
P=20
mu=3
phi=.8
rho=0.2
Pcat=NULL
pNAs=0
sp_indx=1:M
sp_indx=matrix(sp_indx,ncol=sqrt(M),byrow=T)

S_true=matrix(0,nrow=TT,ncol=M)

C=Cmatrix(sp_indx)
# simDat=sim_spatial_JM(P,C,seed,
#                       rho=rho,Pcat=Pcat, phi=phi,
#                       mu=mu,pNAs=pNAs)

Y=NULL


# States at time 1
t=1
simDat=sim_spatiotemp_JM(P,C,t,
                      rho=rho,Pcat=Pcat, phi=phi,
                      mu=mu,pNAs=pNAs,ST=NULL)
temp=simDat$SimData
temp$m=1:M
temp$t=rep(t,M)
Y=rbind(Y,temp)
S_true[t,]=simDat$s
S_true[t,]=order_states_condMean(Y$V20[Y$t==t],S_true[t,])

# Temporal persistence 
PI=.7
for(t in 2:TT){
  simDat=sim_spatiotemp_JM(P,C,t,
                        rho=rho,Pcat=Pcat, phi=phi,
                        mu=mu,pNAs=pNAs,ST=S_true[t-1,],PI=PI)
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

lambda=.1
gamma=.05
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

tt=4

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

