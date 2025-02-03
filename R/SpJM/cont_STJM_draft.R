# Load necessary libraries
library(philentropy)  # For Hellinger distance

# Function to compute Hellinger distance between two probability vectors
hellinger_distance <- function(p, q) {
  return(sqrt(sum((sqrt(p) - sqrt(q))^2) / 2))
}

# Function to compute Manhattan distance
manhattan_distance <- function(vec1, vec2) {
  return(sum(abs(vec1 - vec2)))
}


# Compute the objective function
compute_objective <- function(SS_m, SS_all, loss_mx_m, 
                              W, gamma, lambda,m) {
  
  # SS contains only state probabilities
  # It cannot contain columns t and m as it is the argument of the objective function
  # Also, it refers to a single point m (see option II)
  
  # Need for an object like a tm_indx (first two columns of former SS)
  
  # Get unique time steps and unique clusters
  TT=dim(loss_mx_m)[1]
  M_vals <- unique(SS_all$m)

  # Initialize objective value
  objective_value <- 0
  
  # Loop over time steps and clusters
  for (t in 1:TT) {

      # Extract s_{t,m} (current soft assignment)
      s_tm <- SS_m[t,]
      
      # Retrieve precomputed Gower distance from loss_mx
      g_term <- sum(s_tm * loss_mx_m[t,])
      
      # Compute Hellinger term h() with other clusters
      h_term <- 0
      for (i in M_vals[M_vals != m]) {
        s_prev <- as.numeric(SS_all[SS_all$t == t & SS_all$m == i, -c(1,2)])  
        # s_{t,i} from previous iteration
        h_term <- h_term + hellinger_distance(s_prev, s_tm) * W[m, i]
      }
      
      # Compute temporal smoothness term (only if t+1 exists)
      if ((t + 1) %in% 1:TT) {
        s_next <- as.numeric(SS_all[SS_all$t == (t + 1) & SS_all$m == m, 
                                    -c(1,2)])
        time_smoothness <- manhattan_distance(s_next, s_tm)
      } else {
        time_smoothness <- 0
      }
      
      # Aggregate terms
      objective_value <- objective_value + g_term + gamma * h_term + 
        (lambda / 4) * time_smoothness
    
  }
  
  return(objective_value)
}

W=ifelse(D==0,0,1/D)

tm_indx=SS[,1:2]
loss_mx_tm=cbind(tm_indx,loss_mx)
m=1
loss_mx_m=loss_mx_tm[which(loss_mx_tm$m==m),-(1:2)]
SS_m=SS[which(SS$m==m),-c(1:2)]
compute_objective(SS_m, SS_all, loss_mx_m, 
                              W, gamma, lambda,m)

# Load necessary libraries
library(ParBayesianOptimization)

# Softmax function to transform unconstrained variables into probability simplex
softmax <- function(v) {
  exp_v <- exp(v - max(v))  # Subtract max for numerical stability
  return(exp_v / sum(exp_v))
}

# Wrapper function for Bayesian Optimization
objective_wrapper <- function(par, SS_all, loss_mx_m, W, gamma, lambda, m) {
  
  # Transform parameters from unconstrained space to probability simplex
  SS_trans <- t(apply(par,1,softmax))
  
  # Compute the objective function for this candidate solution
  ret=compute_objective(SS_m = matrix(SS_trans, nrow = dim(loss_mx_m)[1], byrow = TRUE), 
                         SS_all = SS_all, 
                         loss_mx_m = loss_mx_m, 
                         W = W, 
                         gamma = gamma, 
                         lambda = lambda, 
                         m = m)
  return(list(Score=-ret))
}

objective_wrapper(SS_m, SS_all, loss_mx_m, W, gamma, lambda, m)

# Bayesian Optimization Setup
bayes_optimize <- function(SS_all, loss_mx_m, W, gamma, lambda, m, K, init_points = 10, n_iter = 30) {
  
  # Define bounds for unconstrained optimization (before Softmax transformation)
  # bounds <- list()
  # for (t in 1:TT) {
  #   for (k in 1:K) {
  #     bounds[[paste0("p_", t, "_", k)]] <- c(0, 1)  # Probabilities must be in [0,1]
  #   }
  # }
  
  # Run Bayesian Optimization
  results <- bayesOpt(
    FUN = function(...) objective_wrapper(par = c(...), SS_all = SS_all, loss_mx_m = loss_mx_m, 
                                          W = W, gamma = gamma, lambda = lambda, m = m),
    bounds = NULL,
    initPoints = init_points,  # Number of random initial points
    iters.n = n_iter,          # Number of optimization iterations
    acq = "ei",                # Expected Improvement acquisition function
    parallel = FALSE
  )
  
  # Extract best solution
  best_params <- getBestPars(results)
  
  # Convert back to valid probability vector using Softmax
  best_s_tm <- softmax(unlist(best_params))
  
  return(best_s_tm)
}

# Example Usage
gamma <- 1.0
lambda <- 0.5
K <- 3  # Number of probability components
m <- 1  # Cluster to optimize

set.seed(42)
best_probability_vector <- bayes_optimize(SS_all, loss_mx_m, W, gamma, lambda, m, K)
cat("Best Probability Vector:", best_probability_vector, "\n")


# Example usage
gamma <- 1.0
lambda <- 0.5
objective_result <- compute_objective(SS, loss_mx, W, gamma, lambda)
cat("Objective Function Value:", objective_result, "\n")

adj.rand.index(apply(fit,1,which.max),Sim$mchain)


### spatiotemp

source("Utils.R")


M=10
TT=50
theta=.01
beta=.9
K=3
mu=1
rho=0
phi=0.8
P=6
Pcat=2
seed=123
pg=0
pNAs=0

result <- generate_spatio_temporal_data(M, TT, theta, beta, K = K,
                                        mu=mu,rho=rho,phi=phi,
                                        P=P,Pcat=Pcat,seed=seed,pGap=pg,pNAs=pNAs)

Y.compl=result$Y
D=result$dist_matrix
Y=result$Y.NA
#Y=Y[,-(3:4)]
head(Y)

jump_penalty = .1
grid_size =.05
verbose=F
tol=NULL
spatial_penalty = .1
alpha=2
n_states=3

prova=cont_STJM(Y,K,D,
                jump_penalty,
                spatial_penalty,
                alpha=2,grid_size=NULL,mode_loss=T,
                rand_search_sample = 5,
                n_init=10,
                max_iter=10,tol=NULL,initial_states=NULL,
                n_cores=NULL,prll=F
)


