library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)

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

K_grid=2:3
lambda_grid=seq(0,1,length.out=2)
m_grid=seq(1.01,2,length.out=2)
B          = 10
max_iter   = 5
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 19

# 0) Compila il C++ una volta sola
Rcpp::sourceCpp("simplex_pgd.cpp")

# 1) Definisci la griglia

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  b = 0:B)

# 2) Parallel apply
res_list <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  mi    <- hp$m[i]
  bi    <- hp$b[i]
  
  Yinput <- if (bi == 0) Y else apply(Y, 2, sample)
  
  attempt <- 1
  last_loss <- NA_real_
  
  repeat {
    result <- tryCatch({
      fit <- fuzzy_jump_cpp(
        Yinput,
        K        = Ki,
        lambda   = li,
        m        = mi,
        max_iter = max_iter,
        n_init   = n_init,
        tol      = tol,
        verbose  = FALSE
      )
      list(success = TRUE, loss = fit$loss)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      last_loss <- result$loss
      break
    }
    if (attempt >= max_retries) {
      warning(sprintf("Row %d: reached max attempts (%d), setting loss=NA", i, max_retries))
      break
    }
    attempt <- attempt + 1
    Sys.sleep(0.1)
  }
  
  list(
    K         = Ki,
    lambda    = li,
    m         = mi,
    permuted  = as.logical(bi),
    loss      = last_loss,
    attempts  = attempt
  )
}, mc.cores = ncores)
# 3) Trasforma in data.frame
results <- do.call(rbind, lapply(res_list, as.data.frame))
results <- as.data.frame(results)