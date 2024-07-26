source("Utils.R")

TT=6
Ktrue=3
seed=1
M=225
P=25
mu=3
phi=.7
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
for(t in 1:TT){
  simDat=sim_spatial_JM(P,C,t,
                        rho=rho,Pcat=Pcat, phi=phi,
                        mu=mu,pNAs=pNAs)
  temp=simDat$SimData
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=order_states_freq(simDat$s)
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

lambda=.05
gamma=.075
initial_states=NULL
max_iter=10
n_init=10
tol=NULL
verbose=T

fit <- STjumpR(Y, n_states = 3, C, jump_penalty=lambda,spatial_penalty = gamma, verbose=T)

best_s=fit$best_s

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
ggarrange(Ptrue, Pest)