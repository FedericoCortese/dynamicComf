source("Utils.R")

# Close to 52 weeks
#TT=50
TT=8
Ktrue=3
seed=1

# Close to 21, number of italian regions 
#M=25
M=25
P=20

# MU=1 REEEEALLY BAD RESULTS, MU=3 GOOD RESULTS
mu=2

phi=.8
rho=0.2
Pcat=10
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
simDat=sim_spatiotemp_JM(P,C,seed=t+10,
                         rho=rho,Pcat=Pcat, phi=phi,
                         mu=mu,pNAs=pNAs,ST=NULL,n_states=Ktrue)
temp=data.frame(simDat$SimData)
temp$m=1:M
temp$t=rep(t,M)
Y=rbind(Y,temp)
S_true[t,]=simDat$s
S_true[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S_true[t,])

# Temporal persistence 
PI=0.9
for(t in 2:TT){
  simDat=sim_spatiotemp_JM(P,C,seed=t+10,
                           rho=rho,Pcat=Pcat, phi=phi,
                           mu=mu,pNAs=pNAs,ST=S_true[t-1,],PI=PI,n_states=Ktrue)
  temp=data.frame(simDat$SimData)
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=simDat$s
  S_true[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S_true[t,])
  
}

# Put t and m in front of all the others with dplyr
Y <- Y %>% select(t, m, everything())

library(shiny)
library(ggplot2)
library(ggpubr)

# Define the UI
ui <- fluidPage(
  titlePanel("States evolution"),
  sidebarLayout(
    sidebarPanel(
      # Slider to select time point (tt)
      sliderInput("time", "Select Time Point:",
                  min = 1, max = TT, value = 1)  # Adjust max as needed
    ),
    mainPanel(
      # Plot output
      plotOutput("statePlot")
    )
  )
)
server <- function(input, output) {
  
  output$statePlot <- renderPlot({
    tt <- input$time
    
    # True state plot
    data_matrix <- matrix(S_true[tt,], nrow = sqrt(M), ncol = sqrt(M), byrow = TRUE)
    data_df <- as.data.frame(as.table(data_matrix))
    colnames(data_df) <- c("X", "Y", "Value")
    
    # Ensure all levels are present
    data_df$Value <- factor(data_df$Value, levels = c("1", "2", "3"))
    
    ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
      geom_tile(color = "white") +
      ggtitle(paste0("t=",tt))+
      scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      coord_fixed() 
    
    
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

lambda=0.1
gamma=0.05
initial_states=NULL
max_iter=10
n_init=10
tol=NULL
verbose=T

st=Sys.time()
fit <- STjumpR(Y, n_states = 3, C, jump_penalty=lambda,spatial_penalty = gamma, verbose=T)
end=Sys.time()
end-st

best_s=fit$best_s
for(t in 1:TT){
  best_s[t,]=order_states_condMean(Y$V20[Y$t==t],best_s[t,])
}

adj.rand.index(S_true,best_s)

# Define the UI
ui <- fluidPage(
  titlePanel("Comparison of True and Estimated States"),
  sidebarLayout(
    sidebarPanel(
      # Slider to select time point (tt)
      sliderInput("time", "Select Time Point:",
                  min = 1, max = TT, value = 1)  # Adjust max as needed
    ),
    mainPanel(
      # Plot output
      plotOutput("statePlot")
    )
  )
)

# Define the server
server <- function(input, output) {
  
  output$statePlot <- renderPlot({
    tt <- input$time
    
    # True state plot
    data_matrix <- matrix(S_true[tt,], nrow = sqrt(M), ncol = sqrt(M), byrow = TRUE)
    data_df <- as.data.frame(as.table(data_matrix))
    colnames(data_df) <- c("X", "Y", "Value")
    
    Ptrue <- ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      coord_fixed() 
    
    # Estimated state plot
    data_matrix <- matrix(best_s[tt,], nrow = sqrt(M), ncol = sqrt(M), byrow = TRUE)
    data_df <- as.data.frame(as.table(data_matrix))
    colnames(data_df) <- c("X", "Y", "Value")
    
    Pest <- ggplot(data_df, aes(x = X, y = Y, fill = factor(Value))) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("1" = "pink", "2" = "orange", "3" = "violet")) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      coord_fixed() 
    
    # Arrange the plots
    arr_plot <- ggarrange(Ptrue, Pest)
    
    # Annotate with time point
    annotate_figure(arr_plot, top = text_grob(paste0("t=", tt), 
                                              color = "black", face = "bold", size = 14))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)


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
out <- potts(init, param=theta, nbatch = 1 , blen=1, nspac=1)

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

data_matrix <- matrix(S_true[3,], nrow = sqrt(M), ncol = sqrt(M),byrow = T)
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
