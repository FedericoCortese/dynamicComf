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

W=ifelse(D==0,0,1/D)

# Compute the objective function
compute_objective <- function(SS, loss_mx, W, gamma, lambda) {
  
  # Get unique time steps and unique clusters
  T_vals <- unique(SS$t)
  M_vals <- unique(SS$m)
  
  # Initialize objective value
  objective_value <- 0
  
  # Loop over time steps and clusters
  for (t in T_vals) {
    for (m in M_vals) {
      
      # Extract s_{t,m} (current soft assignment)
      s_tm <- as.numeric(SS[SS$t == t & SS$m == m, -c(1,2)])
      
      # Retrieve precomputed Gower distance from loss_mx
      g_term <- sum(s_tm * loss_mx[m, ])
      
      # Compute Hellinger term h() with other clusters
      h_term <- 0
      for (i in M_vals[M_vals != m]) {
        s_prev <- as.numeric(SS[SS$t == t & SS$m == i, -c(1,2)])  
        # s_{t,i} from previous iteration
        h_term <- h_term + hellinger_distance(s_prev, s_tm) * W[m, i]
      }
      
      # Compute temporal smoothness term (only if t+1 exists)
      if ((t + 1) %in% T_vals) {
        s_next <- as.numeric(SS[SS$t == (t + 1) & SS$m == m, -c(1,2)])
        time_smoothness <- manhattan_distance(s_next, s_tm)
      } else {
        time_smoothness <- 0
      }
      
      # Aggregate terms
      objective_value <- objective_value + g_term + gamma * h_term + (lambda / 4) * time_smoothness
    }
  }
  
  return(objective_value)
}

<<<<<<< HEAD
# Example usage
gamma <- 1.0
lambda <- 0.5
objective_result <- compute_objective(SS, loss_mx, W, gamma, lambda)
cat("Objective Function Value:", objective_result, "\n")
=======
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

>>>>>>> 943d4952e25e83d7005c5b6f1beeae96a164c2a9
