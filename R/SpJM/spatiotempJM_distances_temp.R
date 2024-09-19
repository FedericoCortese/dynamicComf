source("Utils.R")

# Example usage
M <- 50 # Number of spatial locations
theta <- .001 # Spatial persistence
beta <- 0.7   # Temporal persistence
K=3
P=15
Pcat=5

pg=0.2
TT <- 10  # Number of time steps
TT=TT+round(TT*pg)
pNAs=0

result <- generate_spatio_temporal_data(M, TT, theta, beta, K = 3,
                                        mu=.5,rho=0,phi=.8,
                                        P=P,Pcat=Pcat,seed=1,pGap=pg,pNAs=pNAs)

Y.compl=result$Y
D=result$dist_matrix
Y=result$Y.NA

# Count % NAs for each column of Y.NA
colMeans(is.na(Y))

# maxD=max(D,na.rm=T)
# DD=D/maxD

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

tf=I(pg>0)
lambda=0.05
gamma=0.1
prova=STjumpDist(Y,3,D,jump_penalty = lambda,spatial_penalty = gamma,verbose = T,timeflag = tf)
best_s=prova$best_s

TY=unique(Y$t)
S_true=result$S[TY,]

for(t in 1:length(TY)){
  best_s[t,]=order_states_condMean(Y[Y$t==TY[t],dim(Y)[2]],best_s[t,])
}
adj.rand.index(best_s,S_true)

ui <- fluidPage(
  titlePanel("Dynamic Plot for Varying t"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("t", "Select t:", 
                  min = min(TY), 
                  max = max(TY), 
                  value = TY[1], 
                  step = NULL)  # step is NULL to allow discrete values
    ),
    
    mainPanel(
      plotOutput("dynamicPlot")
    )
  )
)

server <- function(input, output, session) {
  # Update slider to only allow values from TY
  observe({
    updateSliderInput(session, "t", 
                      min = min(TY), 
                      max = max(TY), 
                      value = TY[1], 
                      step = NULL)
  })
  
  output$dynamicPlot <- renderPlot({
    # Ensure t is one of the values from TY
    t_value <- input$t
    
    # Define color mapping based on selected t_value
    par(mfrow=c(1,2))
    plot(result$spatial_points, 
         col = S_true[t_value, ], 
         pch = 19, cex = 1.5, 
         main = paste("True"))
    plot(result$spatial_points, 
         col = best_s[t_value, ], 
         pch = 19, cex = 1.5, 
         main = bquote("ST-JM with " ~ lambda == .(lambda) ~ " and " ~ gamma == .(gamma)))
  })
}

# Run the app
shinyApp(ui = ui, server = server)
