TT=500
P=50
Pcat=NULL
lambda=.1
K=3
mu=1.5
rho=0
nu=4
phi=.8
pers=.99
n_init=5
n_inner=10
n_outer=10
tol=NULL
verbose=T
alpha=.1


source("Utils_sparse_robust.R")

simDat=sim_data_stud_t(seed=123,
                       TT=TT,
                       P=P,
                       Pcat=Pcat,
                       Ktrue=K,
                       mu=mu,
                       rho=rho,
                       nu=nu,
                       phi=phi,
                       pers=pers)

Y=simDat$SimData
feat_list=list()
feat_list[[1]]=1:10
feat_list[[2]]=5:15
feat_list[[3]]=15:25
Y_noised=apply_noise_by_cluster(Y,simDat$mchain,feat_list)

Y=Y_noised

zeta0=.25
lambda=.25

fit=sparse_robust_fit(Y=Y,K=K,
                      lambda=lambda,
                      zeta0=zeta0,
                      alpha=.1,verbose=T,tol=1e-4,
                      n_init=5,n_outer=20)

round(fit$W,2)
mclust::adjustedRandIndex(simDat$mchain,fit$best_s)

library(reshape)
df <- as.data.frame(fit$W)
df$Cluster <- factor(paste0("Cluster_", 1:nrow(df)))

# Riorganizziamo in formato lungo
df_long <- melt(df, id.vars = "Cluster", variable.name = "Feature", value.name = "Weight")

# Converti Feature in fattore per ordinare le colonne
#df_long$Feature <- factor(df_long$Feature, levels = paste0("V", 1:ncol(fit$W)))

# Bar plot
library(ggplot2)
ggplot2::ggplot(df_long, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Feature Weights by Cluster",
       x = "Feature",
       y = "Weight") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
