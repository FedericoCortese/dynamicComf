library(parallel)
source("Utils.R")

lambda=seq(0,1,by=.1)
TT=c(100,1000)
# Errore quando vario i P, potrei separare gli studi simulati per ogni P
# P=c(4,20,50)
P=20
seeds=1:100
hp=expand.grid(TT=TT,P=P,lambda=lambda,seed=seeds)

# Calcola i core per livello esterno ed interno
n_cores_total=parallel::detectCores()
frac_ext=1/3
n_cores_ext=n_cores_total*frac_ext
n_cores_int=(n_cores_total-n_cores_ext)/n_cores_ext

# n_cores_ext=30

# Setup 1
st=Sys.time()
contJM_setup1_P20 <- parallel::mclapply(1:nrow(hp),
                                    function(x)
                                      simstud_contJM(seed=hp[x,]$seed,
                                                     lambda=hp[x,]$lambda,
                                                     TT=hp[x,]$TT,
                                                     P=hp[x,]$P,
                                                     Ktrue=3,mu=1,
                                                     phi=.8,rho=0,
                                                     Pcat=NULL,pers=.95,
                                                     pNAs=0,typeNA=3,
                                                     prll=T,
                                                     n_cores_int=n_cores_int),
                                    mc.cores = n_cores_ext)
en=Sys.time()
elapsed=en-st
save(contJM_setup1_P20,elapsed,file="contJM_setup1_P20.Rdata")



# Results -----------------------------------------------------------------

hp_res=hp
hp_res$ARI=NA
for (i in 1:nrow(hp_res)){
  simDat=sim_data_mixed_missmec(seed=hp_res[i,]$seed,
                                TT=hp_res[i,]$TT,
                                P=hp_res[i,]$P,
                                Ktrue=3,
                                mu=1,
                                phi=.8,
                                rho=0,
                                Pcat=NULL,
                                pers=.95,
                                pNAs=0,
                                typeNA=3)
  est_seq=apply(contJM_setup1_P20[[i]],1,which.max)
  hp_res$ARI[i]=mclust::adjustedRandIndex(simDat$mchain,est_seq)
  
}

# Calculate average ARI and standard deviation for each combination of P, lambda, and TT
library(dplyr)

# Summarize the data
res_summary <- hp_res %>%
  group_by(P, lambda, TT) %>%
  summarise(
    avg_ARI = mean(ARI),
    sd_ARI = sd(ARI),
    .groups = 'drop'
  )

# Find the lambda with the maximum average ARI for each combination of P and TT
result <- res_summary %>%
  group_by(P, TT) %>%
  filter(avg_ARI == max(avg_ARI)) %>%
  summarise(
    optimal_lambda = lambda[1], # Lambda with max avg_ARI
    max_avg_ARI = avg_ARI[1],   # Maximum average ARI
    corresponding_sd = sd_ARI[1], # Standard deviation for this lambda
    .groups = 'drop'
  )

ggplot(res_summary, aes(x = lambda, y = avg_ARI, color = interaction(P, TT))) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = avg_ARI - sd_ARI, ymax = avg_ARI + sd_ARI, fill = interaction(P, TT)), alpha = 0.2) +
  labs(
    title = "Average ARI vs Lambda with Confidence Bands",
    x = "Lambda",
    y = "Average ARI",
    color = "P, TT",
    fill = "P, TT"
  ) +
  theme_minimal()
