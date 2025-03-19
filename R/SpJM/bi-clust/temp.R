lambda=.1
gamma=.05
seed=123
M=50
TT=100
mu=.5
rho=.2
phi=.8
beta=.99
theta=.001
K=3
P=6
phi=.8
Pcat=2
pNAs=0
pg=0
pers=.95
m=1.01
source("Utils_biclust.R")

result <- generate_spatio_temporal_data(M, TT, theta, beta, K = K,
                                        mu=mu,rho=rho,phi=phi,
                                        P=P,Pcat=Pcat,seed=seed,pGap=pg,pNAs=pNAs)
Y=result$Y.NA

tmp=biclust_JM(Y=Y,K=3,
                           jump_penalty=.1,
                           alpha=2,grid_size=.05,
                           mode_loss=T,rand_search_sample=100,
                           n_init=10,
                           max_iter=10,tol=NULL,initial_states=NULL,
                           n_cores=NULL,prll=F,ncores_M=3
                           #,parallel_ga=F
)

tmp$MAP=as.factor(apply(tmp[,-(1:2)],1,which.max))

head(tmp)

library(ggplot2)
library(reshape2)
data_melted <- dcast(tmp, t ~ m, value.var = "MAP")

# Convertiamo in formato matrice (escludendo la colonna t)
rownames(data_melted) <- data_melted$t
data_melted$t <- NULL
data_matrix <- as.matrix(data_melted)

# Creazione della heatmap con ggplot2
heatmap_plot <- ggplot(melt(data_matrix), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  #scale_fill_gradient(low = "white", high = "red") + 
  labs(title = "Heatmap di MAP", x = "m", y = "t", fill = "MAP") +
  theme_minimal()

# Visualizziamo il grafico
print(heatmap_plot)
