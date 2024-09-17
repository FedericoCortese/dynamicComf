source("Utils.R")

# Example usage
M <- 20# Number of spatial locations
TT <- 5   # Number of time steps
theta <- .1 # Spatial persistence
beta <- .7    # Temporal persistence
K=3
P=20
Pcat=10
result <- generate_spatio_temporal_data(M, TT, theta, beta, K = 3,
                                        mu=1,rho=0.2,phi=.8,
                                        P=P,Pcat=Pcat,seed=1)

Y=result$Y
D=result$dist_matrix

library(shiny)
# Define UI
ui <- fluidPage(
  titlePanel("Dynamic Plot for Varying t"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("t", "Select t:", 
                  min = 1, max = TT, 
                  value = 1, step = 1)
    ),
    
    mainPanel(
      plotOutput("dynamicPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  output$dynamicPlot <- renderPlot({
    t <- input$t
    
    # Define color mapping based on selected t
    plot(result$spatial_points, 
         col = result$S[t, ], 
         pch = 19, cex = 1.5, 
         main = paste("Plot for t =", t))
  })
}

# Run the app
shinyApp(ui = ui, server = server)

lambda=0.05
gamma=0.05
prova=STjumpDist(Y,3,D,jump_penalty = lambda,spatial_penalty = gamma,verbose = T,timeflag = F)
best_s=prova$best_s
for(t in 1:TT){
  best_s[t,]=order_states_condMean(Y$V20[Y$t==t],best_s[t,])
}
adj.rand.index(best_s,result$S)

server <- function(input, output) {
  
  output$dynamicPlot <- renderPlot({
    t <- input$t
    
    # Define color mapping based on selected t
    par(mfrow=c(1,2))
    plot(result$spatial_points, 
         col = result$S[t, ], 
         pch = 19, cex = 1.5, 
         main = paste("True"))
    plot(result$spatial_points, 
         col = best_s[t, ], 
         pch = 19, cex = 1.5, 
         main = paste("Estimated with lambda = ", lambda," gamma = ", gamma))
  })
}

# Run the app
shinyApp(ui = ui, server = server)
