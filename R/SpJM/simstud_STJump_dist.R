library(parallel)
library(plotly)
library(dplyr)

source("Utils.R")

mu=.5
rho=.2
phi=.8
beta=.9
theta=.01
K=3
P=20
Pcat=10

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
M=c(10,50)
TT=c(10,50)

hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)


# 20% gaps, no missing ----------------------------------------------------

pNAs=0
pg=.2


start_STJsim=Sys.time()
STJsim_pg <- parallel::mclapply(1:nrow(hp),
                             function(x)
                               simstud_STJump_dist(lambda=hp[x,]$lambda,
                                                   gamma=hp[x,]$gamma,
                                                   seed=hp[x,]$seed,
                                                   M=hp[x,]$M,
                                                   TT=hp[x,]$TT,
                                                   beta=beta, 
                                                   theta=theta,
                                                   mu=mu,
                                                   rho=rho,
                                                   K=K,P=P,
                                                   phi=phi,
                                                   Pcat=Pcat,
                                                   pNAs=pNAs,
                                                   pg=pg),
                             mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim_pg,elapsed_STJsim,file="STJsim_dist_pg.RData")


load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/STJsim_dist_pg.RData")
# res=data.frame(hp,ARI=unlist(lapply(STJsim,function(x)x$ARI)))
# 
# res_av=res%>%group_by(M,TT,lambda,gamma)%>%summarise(avARI=mean(ARI,na.rm=T),
#                                                      sdARI=sd(ARI,na.rm=T))
# 
# res_max=res_av%>%group_by(M,TT)%>%summarise(maxARI=max(avARI),
#                                             sd_maxARI=sdARI[which.max(avARI)],
#                                             lambdaARI=lambda[which.max(avARI)],
#                                             gammaARI=gamma[which.max(avARI)])

BAC_gap20=as.vector(unlist(lapply(STJsim_pg,BAC)))
res_BAC_gap20=data.frame(hp,BAC=BAC_gap20)
res_BAC_gap20_av=res_BAC_gap20%>%group_by(M,TT,lambda,gamma)%>%summarise(avBAC=mean(BAC,na.rm=T),
                                                                       sdBAC=sd(BAC,na.rm=T))
res_BAC_gap20_max=res_BAC_gap20_av%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                              sd_maxBAC=sdBAC[which.max(avBAC)],
                                                              lambdaBAC=lambda[which.max(avBAC)],
                                                              gammaBAC=gamma[which.max(avBAC)])


# # Extract unique combinations of T and M
# unique_combinations <- unique(res_av[, c("TT", "M")])
# 
# # Initialize an empty list to store plots
# plots <- list()
# 
# # Loop over each combination of T and M to create surface plots
# for (i in 1:nrow(unique_combinations)) {
#   T_value <- as.numeric(unique_combinations[i, "TT"])
#   M_value <- as.numeric(unique_combinations[i, "M"])
#   
#   # Filter data for the specific T and M combination
#   subset_data <- res_av[which(res_av$M==M_value&res_av$TT==T_value),]
#   #subset(res_av, TT == T_value & M == M_value)
#   
#   # Create a matrix or grid for the surface plot
#   lambda_values <- unique(subset_data$lambda)
#   gamma_values <- unique(subset_data$gamma)
#   
#   # Create a matrix of ARI values
#   z_matrix <- matrix(NA, nrow = length(lambda_values), ncol = length(gamma_values))
#   for (j in 1:length(lambda_values)) {
#     for (k in 1:length(gamma_values)) {
#       z_matrix[j, k] <- subset_data$avARI[subset_data$lambda == lambda_values[j] & subset_data$gamma == gamma_values[k]]
#     }
#   }
#   
#   # Create the surface plot
#   # p <- plot_ly(x = lambda_values, y = gamma_values, z = z_matrix, type = "surface", showscale = FALSE) %>%
#   #   layout(
#   #     title = list(text = paste("T =", T_value,"<br>", "M =", M_value), y = 0.9, x = 0.2, xanchor = 'center', yanchor = 'top'),
#   #     scene = list(
#   #       xaxis = list(title = "x=lambda"),
#   #       yaxis = list(title = "y=gamma"),
#   #       zaxis = list(title = "ARI")
#   #     ),
#   #     showlegend = FALSE
#   #   )
#   p <- plot_ly(x = lambda_values, y = gamma_values, z = z_matrix, type = "surface", showscale = FALSE,
#                hovertemplate = paste(
#                  "<b>λ</b>: %{x}<br>",
#                  "<b>γ</b>: %{y}<br>",
#                  "<b>ARI</b>: %{z}<extra></extra>"
#                )) %>%
#     layout(
#       title = list(text = paste("T =", T_value, "<br>", "M =", M_value), y = 0.9, x = 0.2, xanchor = 'center', yanchor = 'top'),
#       scene = list(
#         xaxis = list(title = "λ"),
#         yaxis = list(title = "γ"),
#         zaxis = list(title = "ARI")
#       ),
#       showlegend = FALSE
#     )
#   
#   
#   
#   # Store the plot in the list
#   plots[[i]] <- p
# }
# 
# 
# plots[[1]]
# plots[[2]]
# plots[[3]]
# plots[[4]]


# 5 % NA ------------------------------------------------------------------

pNAs=0.05
pg=0

start_STJsim=Sys.time()
STJsim_NA5 <- parallel::mclapply(1:nrow(hp),
                                function(x)
                                  simstud_STJump_dist(lambda=hp[x,]$lambda,
                                                      gamma=hp[x,]$gamma,
                                                      seed=hp[x,]$seed,
                                                      M=hp[x,]$M,
                                                      TT=hp[x,]$TT,
                                                      beta=beta, 
                                                      theta=theta,
                                                      mu=mu,
                                                      rho=rho,
                                                      K=K,P=P,
                                                      phi=phi,
                                                      Pcat=Pcat,
                                                      pNAs=pNAs,
                                                      pg=pg),
                                mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim_NA5,elapsed_STJsim,file="STJsim_dist_NA5.RData")

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/STJsim_dist_NA5.RData")
# res_NA5=data.frame(hp,ARI=unlist(lapply(STJsim_NA5,function(x)x$ARI)))
# 
# res_NA5_av=res_NA5%>%group_by(M,TT,lambda,gamma)%>%summarise(avARI=mean(ARI,na.rm=T),
#                                                              sdARI=sd(ARI,na.rm=T))
# 
# res_NA5_max=res_NA5_av%>%group_by(M,TT)%>%summarise(maxARI=max(avARI),
#                                                     sd_maxARI=sdARI[which.max(avARI)],
#                                                     lambdaARI=lambda[which.max(avARI)],
#                                                     gammaARI=gamma[which.max(avARI)])
# 
# BAC
BAC_NA5=as.vector(unlist(lapply(STJsim_NA5,BAC)))
res_BAC_NA5=data.frame(hp,BAC=BAC_NA5)
res_BAC_NA5_av=res_BAC_NA5%>%group_by(M,TT,lambda,gamma)%>%summarise(avBAC=mean(BAC,na.rm=T),
                                                                     sdBAC=sd(BAC,na.rm=T))
res_BAC_NA5_max=res_BAC_NA5_av%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                            sd_maxBAC=sdBAC[which.max(avBAC)],
                                                            lambdaBAC=lambda[which.max(avBAC)],
                                                            gammaBAC=gamma[which.max(avBAC)])
# 
# 

# 20 % NA -----------------------------------------------------------------

pNAs=0.2
pg=0

start_STJsim=Sys.time()
STJsim_NA <- parallel::mclapply(1:nrow(hp),
                             function(x)
                               simstud_STJump_dist(lambda=hp[x,]$lambda,
                                                   gamma=hp[x,]$gamma,
                                                   seed=hp[x,]$seed,
                                                   M=hp[x,]$M,
                                                   TT=hp[x,]$TT,
                                                   beta=beta, 
                                                   theta=theta,
                                                   mu=mu,
                                                   rho=rho,
                                                   K=K,P=P,
                                                   phi=phi,
                                                   Pcat=Pcat,
                                                   pNAs=pNAs,
                                                   pg=pg),
                             mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim_NA,elapsed_STJsim,file="STJsim_dist_NA20.RData")

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/STJsim_dist_NA20.RData")
# res_NA20=data.frame(hp,ARI=unlist(lapply(STJsim_NA,function(x)x$ARI)))
# 
# res_NA20_av=res_NA%>%group_by(M,TT,lambda,gamma)%>%summarise(avARI=mean(ARI,na.rm=T),
#                                                              sdARI=sd(ARI,na.rm=T))
# 
# res_NA20_max=res_NA_av%>%group_by(M,TT)%>%summarise(maxARI=max(avARI),
#                                             lambdaARI=lambda[which.max(avARI)],
#                                             gammaARI=gamma[which.max(avARI)])
# 
# # BAC 
BAC_NA20=as.vector(unlist(lapply(STJsim_NA,BAC)))
res_BAC_NA20=data.frame(hp,BAC=BAC_NA20)
res_BAC_NA20_av=res_BAC_NA20%>%group_by(M,TT,lambda,gamma)%>%summarise(avBAC=mean(BAC,na.rm=T),
                                                                     sdBAC=sd(BAC,na.rm=T))
res_BAC_NA20_max=res_BAC_NA20_av%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                            sd_maxBAC=sdBAC[which.max(avBAC)],
                                                            lambdaBAC=lambda[which.max(avBAC)],
                                                            gammaBAC=gamma[which.max(avBAC)])
