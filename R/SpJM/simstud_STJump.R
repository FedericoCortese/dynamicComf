source("Utils.R")


# mu=1 --------------------------------------------------------------------

K=3
P=20
###
mu=1
rho=0
###
phi=.8
Pcat=10
pNAs=0

PI=.9

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
M=c(25,400)
TT=c(5,50)

hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)



start_STJsim=Sys.time()
STJsim <- parallel::mclapply(1:nrow(hp),
                             function(x)
                               simstud_STJump(lambda=hp[x,]$lambda,
                                          gamma=hp[x,]$gamma,
                                          seed=hp[x,]$seed,
                                          M=hp[x,]$M,
                                          TT=hp[x,]$TT,
                                          mu=mu,rho=rho,
                                          K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs,PI=PI),
                             mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim,elapsed_STJsim,file="STJsim.RData")



# Plot results ------------------------------------------------------------
K=3
P=20
mu=1
rho=0
phi=.8
Pcat=10
pNAs=0
PI=.9
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:100
M=c(25,400)
TT=c(5,50)
hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)
load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_SpatioTemporal JM/STJsim.RData")

library(dplyr)
res=data.frame(hp,ARI=unlist(lapply(STJsim,function(x)x$ARI)))

res_av=res%>%group_by(M,TT,lambda,gamma)%>%summarise(avARI=mean(ARI,na.rm=T))

res_max=res_av%>%group_by(M,TT)%>%summarise(maxARI=max(avARI),
                                            lambdaARI=lambda[which.max(avARI)],
                                            gammaARI=gamma[which.max(avARI)])

# Load necessary library
library(plotly)

# Extract unique combinations of T and M
unique_combinations <- unique(res_av[, c("TT", "M")])

# Initialize an empty list to store plots
plots <- list()

# Loop over each combination of T and M to create surface plots
for (i in 1:nrow(unique_combinations)) {
  T_value <- as.numeric(unique_combinations[i, "TT"])
  M_value <- as.numeric(unique_combinations[i, "M"])
  
  # Filter data for the specific T and M combination
  subset_data <- res_av[which(res_av$M==M_value&res_av$TT==T_value),]
  #subset(res_av, TT == T_value & M == M_value)
  
  # Create a matrix or grid for the surface plot
  lambda_values <- unique(subset_data$lambda)
  gamma_values <- unique(subset_data$gamma)
  
  # Create a matrix of ARI values
  z_matrix <- matrix(NA, nrow = length(lambda_values), ncol = length(gamma_values))
  for (j in 1:length(lambda_values)) {
    for (k in 1:length(gamma_values)) {
      z_matrix[j, k] <- subset_data$avARI[subset_data$lambda == lambda_values[j] & subset_data$gamma == gamma_values[k]]
    }
  }
  
  # Create the surface plot
  # p <- plot_ly(x = lambda_values, y = gamma_values, z = z_matrix, type = "surface", showscale = FALSE) %>%
  #   layout(
  #     title = list(text = paste("T =", T_value,"<br>", "M =", M_value), y = 0.9, x = 0.2, xanchor = 'center', yanchor = 'top'),
  #     scene = list(
  #       xaxis = list(title = "x=lambda"),
  #       yaxis = list(title = "y=gamma"),
  #       zaxis = list(title = "ARI")
  #     ),
  #     showlegend = FALSE
  #   )
  p <- plot_ly(x = lambda_values, y = gamma_values, z = z_matrix, type = "surface", showscale = FALSE,
               hovertemplate = paste(
                 "<b>λ</b>: %{x}<br>",
                 "<b>γ</b>: %{y}<br>",
                 "<b>ARI</b>: %{z}<extra></extra>"
               )) %>%
    layout(
      title = list(text = paste("T =", T_value, "<br>", "M =", M_value), y = 0.9, x = 0.2, xanchor = 'center', yanchor = 'top'),
      scene = list(
        xaxis = list(title = "λ"),
        yaxis = list(title = "γ"),
        zaxis = list(title = "ARI")
      ),
      showlegend = FALSE
    )
  
  
  
  # Store the plot in the list
  plots[[i]] <- p
}


plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]

