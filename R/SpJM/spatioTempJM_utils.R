source("Utils.R")

K=3
P=30
###
mu=1
rho=0.2
###
phi=.8
Pcat=10
pNAs=0

###
lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
seed=1:2
M=100
TT=10
hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed,M=M,TT=TT)
# seed=1:100
# M=c(25,100)
# TT=c(10,50)

STJump_sim=function(lambda,gamma,seed,M,TT,
                    mu=1,rho=0.2,
                    K=3,P=30,phi=.8,Pcat=10,pNAs=0){
  
  sp_indx=1:M
  sp_indx=matrix(sp_indx,ncol=sqrt(M),byrow=T)
  S_true=matrix(0,nrow=TT,ncol=M)
  C=Cmatrix(sp_indx)
  Y=NULL
  
  t=1
  simDat=sim_spatiotemp_JM(P,C,seed=seed,
                           rho=rho,Pcat=Pcat, phi=phi,
                           mu=mu,pNAs=pNAs,ST=NULL,n_states=K)
  temp=data.frame(simDat$SimData)
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=simDat$s
  S_true[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S_true[t,])
  
  # Temporal persistence 
  PI=0.7
  for(t in 2:TT){
    simDat=sim_spatiotemp_JM(P,C,seed=seed,
                             rho=rho,Pcat=Pcat, phi=phi,
                             mu=mu,pNAs=pNAs,ST=S_true[t-1,],PI=PI,n_states=K)
    temp=data.frame(simDat$SimData)
    temp$m=1:M
    temp$t=rep(t,M)
    Y=rbind(Y,temp)
    S_true[t,]=simDat$s
    S_true[t,]=order_states_condMean(Y[Y$t==t,dim(Y)[2]-2],S_true[t,])
    
  }
  
  # Put t and m in front of all the others with dplyr
  Y <- Y %>% select(t, m, everything())
  
  st=Sys.time()
  fit <- STjumpR(Y, n_states = K, C=C, jump_penalty=lambda,spatial_penalty = gamma, verbose=F)
  end=Sys.time()
  elapsed=end-st
  
  best_s=fit$best_s
  for(t in 1:TT){
    best_s[t,]=order_states_condMean(Y$V20[Y$t==t],best_s[t,])
  }
  
  ARI=adj.rand.index(S_true,best_s)
  
  
  return(list(lambda=lambda,gamma=gamma,ARI=ARI,
              seed=seed,
              M=M,TT=TT,
              mu=mu,rho=rho,
              K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs,
              elapsed=elapsed,
              S_true=S_true,best_s=best_s))
}


start_STJsim=Sys.time()
STJsim <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        STJump_sim(lambda=hp[x,]$lambda,
                                                   gamma=hp[x,]$gamma,
                                                   seed=hp[x,]$seed,
                                                   M=hp[x,]$M,
                                                   TT=hp[x,]$TT,
                                                   mu=mu,rho=rho,
                                                   K=K,P=P,phi=phi,Pcat=Pcat,pNAs=pNAs),
                                      mc.cores = parallel::detectCores())
end_STJsim=Sys.time()
elapsed_STJsim=end_STJsim-start_STJsim
save(STJsim,elapsed_STJsim,file="STJsim.RData")

# Graph -------------------------------------------------------------------
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


