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

map_counts <- tmp %>%
  group_by(m) %>%
  summarise(
    count_1 = sum(MAP == 1),
    count_2 = sum(MAP == 2),
    count_3 = sum(MAP == 3)
  )

# Ordinamento delle colonne: prima MAP = 1, poi MAP = 2, poi MAP = 3
map_counts <- map_counts %>%
  arrange(desc(count_1), desc(count_2), desc(count_3))


# Creiamo una tabella pivot con t sulle righe e m sulle colonne
data_matrix <- dcast(data, t ~ m, value.var = "MAP")

# Riordinamento delle colonne secondo il criterio stabilito
ordered_cols <- c("t", as.character(map_counts$m))  # Manteniamo la colonna t
data_matrix <- data_matrix[, ordered_cols]

# Convertiamo il dataframe in formato lungo per ggplot
data_melted <- melt(data_matrix, id.vars = "t")

# Convertiamo in formato matrice (escludendo la colonna t)
rownames(data_melted) <- data_melted$t
data_melted$t <- NULL
data_matrix <- as.matrix(data_melted)

# Creazione della heatmap con ggplot2
heatmap_plot <- ggplot(data_melted, aes(x = variable, y = factor(t), fill = value)) +
  geom_tile() +
 # scale_fill_gradient(low = "white", high = "red") + 
  labs(title = "Heatmap di MAP con colonne riordinate", x = "m (riordinato)", y = "t", fill = "MAP") +
  theme_minimal()

# Visualizziamo il grafico
print(heatmap_plot)
