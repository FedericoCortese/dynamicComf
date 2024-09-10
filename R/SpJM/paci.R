library(MASS)  # For multivariate normal sampling
library(spdep)  # For spatial correlation

# Function to create spatial points
generate_spatial_points <- function(n, max_distance = 10) {
  x <- runif(n, 0, max_distance)
  y <- runif(n, 0, max_distance)
  return(data.frame(x = x, y = y))
}

# Function to generate spatio-temporal data with spatial (theta) and temporal (rho) persistence
generate_spatio_temporal_data <- function(M, TT, theta, rho, K = 4) {
  
  # Generate spatial points
  spatial_points <- generate_spatial_points(M)
  
  # Create spatial covariance matrix based on distance and theta
  dist_matrix <- as.matrix(dist(spatial_points)) # Eventually substitute with Gower distance
  
  spatial_cov <- exp(-theta * dist_matrix)
  
  # Initialize data matrix
  data <- array(0, dim = c(M, TT))
  clusters <- array(0, dim = c(M, TT))
  
  # Initial time step data (from spatial process)
  data[, 1] <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)
  
  cluster_levels <- quantile(data[,1], probs = seq(0, 1, length.out = K + 1))
  clusters[,1] <- cut(data[,1], breaks = cluster_levels, labels = FALSE,
                                               include.lowest =T)
  
  # Generate data for each subsequent time step
  for (t in 2:TT) {
    eta_t <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)  # Spatial noise
    #data[, t] <- rho * data[, t - 1] + eta_t  # Temporal correlation
    data[, t] <- rho * data[, t - 1] + (1-rho) * eta_t
    cluster_levels <- quantile(data[,t], probs = seq(0, 1, length.out = K + 1))
    clusters[,t] <- cut(data[,t], breaks = cluster_levels, labels = FALSE,
                        include.lowest =T)
    
  }
  
  # Discretize data into clusters
  
  # cluster_levels <- quantile(data, probs = seq(0, 1, length.out = K + 1))
  # clusters <- apply(data, 2, function(col) cut(col, breaks = cluster_levels, labels = FALSE,
  #                                              include.lowest =T))
  
  return(list(data = data, 
              clusters = clusters, 
              spatial_points = spatial_points,
              spatial_cov = spatial_cov))
}

# Example usage
M <- 100  # Number of spatial locations
TT <- 4   # Number of time steps
theta <- 0.01 # Spatial persistence (the lower the more persistent)
rho <- .9    # Temporal persistence
K=4
result <- generate_spatio_temporal_data(M, TT, theta, rho,K)

plot(result$spatial_points,col=result$clusters[,1],pch=19,cex=1.5)
plot(result$spatial_points,col=result$clusters[,2],pch=19,cex=1.5)
plot(result$spatial_points,col=result$clusters[,3],pch=19,cex=1.5)
plot(result$spatial_points,col=result$clusters[,4],pch=19,cex=1.5)

