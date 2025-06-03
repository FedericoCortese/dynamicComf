source("Utils_fuzzyJM_2.R")

# prv=sim_data_stud_t(seed=123,
#                              TT=1000,
#                              P=10,
#                              Pcat=NULL,
#                              Ktrue=2,
#                              mu=1,
#                              rho=0,
#                              nu=4,
#                              phi=.8,
#                              pers=.99)

prv=simulate_fuzzy_mixture_mv(
    TT = 1000,
    P = 5,
    mu = 1,
    Sigma_rho = 0,
    ar_rho = 0.9,
    tau = 0.5,
    seed = NULL
)

plot(prv$pi_1,type='l')

Y=prv$SimData

temp=fuzzy_jump_cpp(Y, 
                    K=2, 
                    lambda = 0.05, 
                    m      = 1.01,
                    max_iter = 20, 
                    n_init   = 10, 
                    tol      = 1e-8, 
                    verbose  = FALSE
                                 
)

temp_par=fuzzy_jump_cpp_parallel(Y, 
                                    K=2, 
                                    lambda = 0.05, 
                                    m      = 1.01,
                                    max_iter = 20, 
                                    n_init   = 10, 
                                    tol      = 1e-8, 
                                    verbose  = FALSE
                                 )

table(temp$MAP,prv$mchain)
table(temp_par$MAP,prv$mchain)

compute_entropy(temp_par$best_S)

temp_par_null=fuzzy_jump_cpp_parallel(Y, 
                                 K=2, 
                                 lambda = 0, 
                                 m      = 1.01,
                                 max_iter = 20, 
                                 n_init   = 10, 
                                 tol      = 1e-8, 
                                 verbose  = FALSE
)

table(temp_par_null$MAP,prv$mchain)

compute_entropy(temp_par_null$best_S)

temp_par_l1=fuzzy_jump_cpp_parallel(Y, 
                                      K=2, 
                                      lambda = 0.2, 
                                      m      = 1.01,
                                      max_iter = 20, 
                                      n_init   = 10, 
                                      tol      = 1e-8, 
                                      verbose  = FALSE
)

table(temp_par_l1$MAP,prv$mchain)

compute_entropy(temp_par_l1$best_S)

temp_gap=fuzzy_jump_gap(Y,
                           K_grid = 2:4,
                           lambda_grid = seq(0, .2, 0.1),
                           B = 50, # numero permutazioni
                           m = 1.01,
                           max_iter = 6,
                           n_init = 10,
                           tol = 1e-10,
                           verbose = FALSE,
                           n_cores = NULL)


# hard vs soft comparison -------------------------------------------------

# hard clustering scenario
source("Utils_fuzzyJM_2.R")

TT = 100
P = 2
mu = 1
Sigma_rho = 0
ar_rho = 0.99
tau = 5

hard_scen=simulate_fuzzy_mixture_mv(
  TT = TT,
  P = P,
  mu = mu,
  Sigma_rho = Sigma_rho,
  ar_rho = ar_rho,
  tau = tau,
  seed = NULL
)

plot(hard_scen$pi_1,type='l')

Y=hard_scen[,(1:P)+1]

fjm_par_m1=fuzzy_jump_cpp_parallel(Y, 
                                   K=2, 
                                   lambda = 0.1, 
                                   m      = 1.01,
                                   max_iter = 15, 
                                   n_init   = 5, 
                                   tol      = 1e-8, 
                                   verbose  = FALSE
)

table(fjm_par_m1$MAP,hard_scen$MAP)
compute_entropy(fjm_par_m1$best_S)

true_distr=data.frame(p1=hard_scen$pi_1,p2=1-hard_scen$pi_1)

plot(fjm_par_m1$best_S[,1],type='l')
lines(hard_scen$pi_1,col='red')

hellinger_ts <- apply(cbind(fjm_par_m1$best_S, true_distr), 1, function(row) {
  p <- row[1:2]
  q <- row[3:4]
  hellinger_distance_vec(p, q)
})

mean(hellinger_ts)

fjm_par_m2=fuzzy_jump_cpp_parallel(Y, 
                                   K=2, 
                                   lambda = 0.1, 
                                   m      = 2,
                                   max_iter = 15, 
                                   n_init   = 5, 
                                   tol      = 1e-8, 
                                   verbose  = FALSE
)

table(fjm_par_m2$MAP,hard_scen$MAP)
compute_entropy(fjm_par_m2$best_S)

true_distr=data.frame(p1=hard_scen$pi_1,p2=1-hard_scen$pi_1)

plot(fjm_par_m2$best_S[,1],type='l')
lines(hard_scen$pi_1,col='red')

hellinger_ts <- apply(cbind(fjm_par_m2$best_S, true_distr), 1, function(row) {
  p <- row[1:2]
  q <- row[3:4]
  hellinger_distance_vec(p, q)
})
mean(hellinger_ts)


# soft clustering scenario

source("Utils_fuzzyJM_2.R")

TT = 500
P = 10
mu = 1
Sigma_rho = 0
ar_rho = 0.9
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

x11()
plot(soft_scen$pi_1,type='l',col='red')

Y=soft_scen[,(1:P)+1]

fjm_par_m1_soft=fuzzy_jump_cpp_parallel(Y, 
                                        K=2, 
                                        lambda = 0.1, 
                                        m      = 1.5,
                                        max_iter = 10, 
                                        n_init   = 10, 
                                        tol      = 1e-8,
                                        ncores   = NULL
                                      
)

table(fjm_par_m1_soft$MAP,soft_scen$MAP)
compute_entropy(fjm_par_m1_soft$best_S)



x11()
plot(fjm_par_m1_soft$best_S[,1],type='l')
lines(soft_scen$pi_1,col='red')

true_distr=data.frame(p1=soft_scen$pi_1,p2=1-soft_scen$pi_1)


hellinger_ts <- apply(cbind(fjm_par_m1_soft$best_S, true_distr), 1, function(row) {
  p <- row[1:2]
  q <- row[3:4]
  hellinger_distance_vec(p, q)
})
mean(hellinger_ts)


fjm_par_m2_soft=fuzzy_jump_cpp_parallel(Y, 
                                        K=2, 
                                        lambda = 0.05, 
                                        m      = 1.5,
                                        max_iter = 15, 
                                        n_init   = 5, 
                                        tol      = 1e-8, 
                                        verbose  = FALSE
)

table(fjm_par_m2_soft$MAP,soft_scen$MAP)
compute_entropy(fjm_par_m2_soft$best_S)

plot(fjm_par_m2_soft$best_S[,1],type='l')
lines(soft_scen$pi_1,col='red')

true_distr=data.frame(p1=soft_scen$pi_1,p2=1-soft_scen$pi_1)


hellinger_ts <- apply(cbind(fjm_par_m2_soft$best_S, true_distr), 1, function(row) {
  p <- row[1:2]
  q <- row[3:4]
  hellinger_distance_vec(p, q)
})
mean(hellinger_ts)


# gap ---------------------------------------------------------------------

source("Utils_fuzzyJM.R")
prv=fuzzyJM_gap(Y,
                        K_grid    = 2:3,
                        lambda_grid = seq(0.01, .5, length.out=2),
                        m_grid     = seq(1.01,2,length.out=2),
                        tol       = NULL,
                        max_iter   = 10,
                        verbose   = FALSE,
                        n_cores   = NULL,
                        B         = 10,
                        n_init=10)
