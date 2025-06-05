library(Rcpp)
sourceCpp("cont_jump.cpp")


source("Utils_fuzzyJM_2.R")

TT = 1000
P = 10
mu = 1
Sigma_rho = 0
ar_rho = 0.99
tau = .2

soft_scen=simulate_fuzzy_mixture_mv(
  TT = TT,
  P = P,
  mu = mu,
  Sigma_rho = Sigma_rho,
  ar_rho = ar_rho,
  tau = tau,
  seed = NULL
)

plot(soft_scen$pi_1,type='l')

Y=soft_scen[,(1:P)+1]

res_cont <- cont_jump(
  Y_in         = as.matrix(Y),
  K            = 2,
  jump_penalty = 1,
  alpha        = 2.0,
  max_iter     = 10,
  n_init       = 5,
  tol          = 1e-8,
  mode_loss    = F,
  grid_size    = 0.01
)

res_fuzzy=fuzzy_jump_cpp(Y, 
                                     K=2, 
                                     lambda=1, 
                                     m=1.25,
                                     max_iter=10, 
                                     n_init=10, tol=1e-8, 
                                     verbose=T
                                     
)

# Comparison

true_S=cbind(soft_scen$pi_1, 1-soft_scen$pi_1)

hellinger_distance_matrix(res_fuzzy$best_S,true_S)
hellinger_distance_matrix(res_cont$best_S,true_S)

plot(soft_scen$pi_1,type='l')
lines(res_fuzzy$best_S[,1],col='red')
lines(1-res_cont$best_S[,1],col='blue')
