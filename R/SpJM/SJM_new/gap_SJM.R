# Install and load Rcpp if not already installed
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  install.packages("Rcpp")
}
library(Rcpp)
Rcpp::sourceCpp("Utils_cpp.cpp")

source("Utils.R")
# Sim data

TT=1000
P=10
Ktrue=2
mu=1
rho=0
nu=6
pers=.95

simDat=sim_data_stud_t(TT=TT, P=P, Ktrue=Ktrue, mu=mu, rho=rho, nu=nu, Pcat=NULL,pers=pers)

Y=simDat$SimData
Y=as.matrix(Y)
ground_truth=simDat$mchain

K_grid=2:4
lambda_grid=seq(0,1,length.out=11)
kappa_grid=seq(1,sqrt(P),length.out=5)
B=10

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  kappa= kappa_grid,
                  b = 0:B)

# 2) Parallel apply
start=Sys.time()
res_list <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  kappai <- hp$kappa[i]
  bi  <- hp$b[i]
  
  Yinput <- if (bi == 0) Y else apply(Y, 2, sample)
  
  attempt <- 1
  last_loss <- NA_real_
  max_retries <- 10
  
  repeat {
    result <- tryCatch({
      fit  <- sparse_jump(
        Y_in        = Yinput,
        n_states    = Ki,
        max_features= kappai,
        jump_penalty= li,
        max_iter    = 10,
        tol         = 1e-8,
        n_init      = 10,
        verbose     = FALSE
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
    kappa         = kappai,
    permuted  = as.logical(bi),
    loss      = last_loss,
    attempts  = attempt
  )
}, mc.cores = ncores)
end=Sys.time()
end-start


