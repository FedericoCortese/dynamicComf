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

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
M=c(10,50)
TT=c(10,50)

hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)


# P=20 --------------------------------------------------------------------
P=20
Pcat=10
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

custom_labeller <- labeller(
  TT = function(value) paste("T =", value),
  M = function(value) paste("M =", value)
)
p_p20_gap20_lambda <- ggplot(res_BAC_gap20_av[res_BAC_gap20_av$gamma <= 0.14 & 
                                                    res_BAC_gap20_av$lambda <= 0.20, ], 
                             aes(x = lambda, y = avBAC, color = factor(gamma))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "Average BAC",
    color = expression(gamma)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )


#png(width = 700, height = 700,filename="p_p10_gap20_lambda.png")
p_p20_gap20_lambda
#dev.off()


p_p20_gap20_gamma=ggplot(res_BAC_gap20_av
                         [res_BAC_gap20_av$gamma <= 0.20 & res_BAC_gap20_av$lambda <= 0.14, ]
                         , 
                         aes(x = gamma, y = avBAC, color = factor(lambda))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(gamma),
    y = "Average BAC",
    color = expression(lambda)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal()+
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )

#jpeg(width = 700, height = 400,filename="p_p10_gap20_gamma.jpeg",quality=100)
p_p20_gap20_gamma


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

custom_labeller <- labeller(
  TT = function(value) paste("T =", value),
  M = function(value) paste("M =", value)
)
p_p20_NA5_lambda <- ggplot(res_BAC_NA5_av[res_BAC_NA5_av$gamma <= 0.14 & 
                                                res_BAC_NA5_av$lambda <= 0.20, ], 
                           aes(x = lambda, y = avBAC, color = factor(gamma))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "Average BAC",
    color = expression(gamma)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )


#png(width = 700, height = 700,filename="p_p10_gap20_lambda.png")
p_p20_NA5_lambda
#dev.off()


p_p20_NA5_gamma=ggplot(res_BAC_NA5_av
                       [res_BAC_NA5_av$gamma <= 0.20 & res_BAC_NA5_av$lambda <= 0.14, ]
                       , 
                       aes(x = gamma, y = avBAC, color = factor(lambda))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(gamma),
    y = "Average BAC",
    color = expression(lambda)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal()+
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )

#jpeg(width = 700, height = 400,filename="p_p10_gap20_gamma.jpeg",quality=100)
p_p20_NA5_gamma

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

custom_labeller <- labeller(
  TT = function(value) paste("T =", value),
  M = function(value) paste("M =", value)
)
p_p20_NA20_lambda <- ggplot(res_BAC_NA20_av[res_BAC_NA20_av$gamma <= 0.14 & 
                                                  res_BAC_NA20_av$lambda <= 0.20, ], 
                            aes(x = lambda, y = avBAC, color = factor(gamma))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "Average BAC",
    color = expression(gamma)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )


#png(width = 700, height = 700,filename="p_p10_gap20_lambda.png")
p_p20_NA20_lambda
#dev.off()


p_p20_NA20_gamma=ggplot(res_BAC_NA20_av[res_BAC_NA20_av$gamma <= 0.20 & 
                                              res_BAC_NA20_av$lambda <= 0.14, ]
                        , 
                        aes(x = gamma, y = avBAC, color = factor(lambda))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(gamma),
    y = "Average BAC",
    color = expression(lambda)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal()+
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )

#jpeg(width = 700, height = 400,filename="p_p10_gap20_gamma.jpeg",quality=100)
p_p20_NA20_gamma

res_BAC_NA20_max=res_BAC_NA20_av%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                            sd_maxBAC=sdBAC[which.max(avBAC)],
                                                            lambdaBAC=lambda[which.max(avBAC)],
                                                            gammaBAC=gamma[which.max(avBAC)])


# P=10 --------------------------------------------------------------------

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

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
M=c(10,50)
TT=c(10,50)

hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)

P=10
Pcat=5


# 20% gaps ----------------------------------------------------------------


pNAs=0
pg=.2

start_STJsim=Sys.time()
STJsim_pg_P10 <- parallel::mclapply(1:nrow(hp),
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
save(STJsim_pg_P10,elapsed_STJsim,file="STJsim_dist_pg_P10.RData")

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/STJsim_dist_pg_P10.RData")

BAC_gap20_P10=as.vector(unlist(lapply(STJsim_pg_P10,BAC)))
res_BAC_gap20_P10=data.frame(hp,BAC=BAC_gap20_P10)
res_BAC_gap20_av_P10=res_BAC_gap20_P10%>%group_by(M,TT,lambda,gamma)%>%summarise(avBAC=mean(BAC,na.rm=T),
                                                                         sdBAC=sd(BAC,na.rm=T))

custom_labeller <- labeller(
  TT = function(value) paste("T =", value),
  M = function(value) paste("M =", value)
)
p_p10_gap20_lambda <- ggplot(res_BAC_gap20_av_P10[res_BAC_gap20_av_P10$gamma <= 0.14 & 
                                             res_BAC_gap20_av_P10$lambda <= 0.20, ], 
                      aes(x = lambda, y = avBAC, color = factor(gamma))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "Average BAC",
    color = expression(gamma)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )


#png(width = 700, height = 700,filename="p_p10_gap20_lambda.png")
p_p10_gap20_lambda
#dev.off()


p_p10_gap20_gamma=ggplot(res_BAC_gap20_av_P10
       [res_BAC_gap20_av_P10$gamma <= 0.20 & res_BAC_gap20_av_P10$lambda <= 0.14, ]
       , 
       aes(x = gamma, y = avBAC, color = factor(lambda))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(gamma),
    y = "Average BAC",
    color = expression(lambda)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal()+
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )

#jpeg(width = 700, height = 400,filename="p_p10_gap20_gamma.jpeg",quality=100)
p_p10_gap20_gamma
#dev.off()

res_BAC_gap20_max_P10=res_BAC_gap20_av_P10%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                                sd_maxBAC=sdBAC[which.max(avBAC)],
                                                                lambdaBAC=lambda[which.max(avBAC)],
                                                                gammaBAC=gamma[which.max(avBAC)])




# 5% NA -------------------------------------------------------------------


pNAs=0.05
pg=0

start_STJsim=Sys.time()
STJsim_NA5_P10 <- parallel::mclapply(1:nrow(hp),
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
save(STJsim_NA5_P10,elapsed_STJsim,file="STJsim_dist_NA5_P10.RData")

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/STJsim_dist_NA5_P10.RData")

BAC_NA5_P10=as.vector(unlist(lapply(STJsim_NA5_P10,BAC)))
res_BAC_NA5_P10=data.frame(hp,BAC=BAC_NA5_P10)
res_BAC_NA5_av_P10=res_BAC_NA5_P10%>%group_by(M,TT,lambda,gamma)%>%summarise(avBAC=mean(BAC,na.rm=T),
                                                                                 sdBAC=sd(BAC,na.rm=T))

custom_labeller <- labeller(
  TT = function(value) paste("T =", value),
  M = function(value) paste("M =", value)
)
p_p10_NA5_lambda <- ggplot(res_BAC_NA5_av_P10[res_BAC_NA5_av_P10$gamma <= 0.14 & 
                                                  res_BAC_NA5_av_P10$lambda <= 0.20, ], 
                             aes(x = lambda, y = avBAC, color = factor(gamma))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "Average BAC",
    color = expression(gamma)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )


#png(width = 700, height = 700,filename="p_p10_gap20_lambda.png")
p_p10_NA5_lambda
#dev.off()


p_p10_NA5_gamma=ggplot(res_BAC_NA5_av_P10
                         [res_BAC_NA5_av_P10$gamma <= 0.20 & res_BAC_NA5_av_P10$lambda <= 0.14, ]
                         , 
                         aes(x = gamma, y = avBAC, color = factor(lambda))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(gamma),
    y = "Average BAC",
    color = expression(lambda)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal()+
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )

#jpeg(width = 700, height = 400,filename="p_p10_gap20_gamma.jpeg",quality=100)
p_p10_NA5_gamma
#dev.off()

res_BAC_NA5_max_P10=res_BAC_NA5_av_P10%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                                        sd_maxBAC=sdBAC[which.max(avBAC)],
                                                                        lambdaBAC=lambda[which.max(avBAC)],
                                                                        gammaBAC=gamma[which.max(avBAC)])


# 20% NA ------------------------------------------------------------------

pNAs=0.2
pg=0

start_STJsim=Sys.time()
STJsim_NA20_P10 <- parallel::mclapply(1:nrow(hp),
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
save(STJsim_NA20_P10,elapsed_STJsim,file="STJsim_dist_NA20_P10.RData")

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/STJsim_dist_NA20_P10.RData")

BAC_NA20_P10=as.vector(unlist(lapply(STJsim_NA20_P10,BAC)))
res_BAC_NA20_P10=data.frame(hp,BAC=BAC_NA20_P10)
res_BAC_NA20_av_P10=res_BAC_NA20_P10%>%group_by(M,TT,lambda,gamma)%>%summarise(avBAC=mean(BAC,na.rm=T),
                                                                             sdBAC=sd(BAC,na.rm=T))

custom_labeller <- labeller(
  TT = function(value) paste("T =", value),
  M = function(value) paste("M =", value)
)
p_p10_NA20_lambda <- ggplot(res_BAC_NA20_av_P10[res_BAC_NA20_av_P10$gamma <= 0.14 & 
                                                  res_BAC_NA20_av_P10$lambda <= 0.20, ], 
                           aes(x = lambda, y = avBAC, color = factor(gamma))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "Average BAC",
    color = expression(gamma)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )


#png(width = 700, height = 700,filename="p_p10_gap20_lambda.png")
p_p10_NA20_lambda
#dev.off()


p_p10_NA20_gamma=ggplot(res_BAC_NA20_av_P10[res_BAC_NA20_av_P10$gamma <= 0.20 & 
                                              res_BAC_NA20_av_P10$lambda <= 0.14, ]
                       , 
                       aes(x = gamma, y = avBAC, color = factor(lambda))) +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  labs(
    x = expression(gamma),
    y = "Average BAC",
    color = expression(lambda)
  ) +
  facet_grid(TT ~ M, labeller = custom_labeller) +
  theme_minimal()+
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(size = 14),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Grey border around each subplot
  )

#jpeg(width = 700, height = 400,filename="p_p10_gap20_gamma.jpeg",quality=100)
p_p10_NA20_gamma
#dev.off()

res_BAC_NA20_max_P10=res_BAC_NA20_av_P10%>%group_by(M,TT)%>%summarise(maxBAC=max(avBAC),
                                                                    sd_maxBAC=sdBAC[which.max(avBAC)],
                                                                    lambdaBAC=lambda[which.max(avBAC)],
                                                                    gammaBAC=gamma[which.max(avBAC)])