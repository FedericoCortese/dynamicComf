library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)


# hard clustering ---------------------------------------------------------

source("Utils_fuzzyJM_2.R")

TT = 1000
P = 10
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

K_grid=2:4
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,2,length.out=3)
B          = 10
max_iter   = 5
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

# 0) Compila il C++ una volta sola
Rcpp::sourceCpp("simplex_pgd.cpp")

# 1) Definisci la griglia

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  b = 0:B)

# 2) Parallel apply
start=Sys.time()
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
end=Sys.time()
start-end
# 3) Trasforma in data.frame
results <- do.call(rbind, lapply(res_list, as.data.frame))
results <- as.data.frame(results)

gap_df <- results %>%
  # Escludiamo eventuali NA di loss
  filter(!is.na(loss)) %>%
  # Raggruppiamo per ogni combinazione di K, lambda e m
  group_by(K, lambda, m) %>%
  summarise(
    # log-loss sul dato non permutato (dovrebbe essere un singolo valore)
    logloss      = log(loss[permuted == FALSE]),
    # aspettativa del log-loss sui dati permutati
    Eg_logloss   = mean(log(loss[permuted == TRUE]), na.rm = TRUE),
    # Gap = E[log(loss_perm)] − log(loss_orig)
    gap          = Eg_logloss - logloss,
    .groups = "drop"
  )

# Dai un’occhiata
print(gap_df)

uniqueKs <- sort(unique(gap_df$K))
plots <- vector("list", length(uniqueKs))  # pre-allocate list
names(plots) <- paste0("K", uniqueKs)

for (i in seq_along(uniqueKs)) {
  k <- uniqueKs[i]
  dfk <- filter(gap_df, K == k)
  
  p_k <- ggplot(dfk, aes(x = lambda, y = gap, group = factor(m), color = factor(m))) +
    geom_line() +
    geom_point() +
    labs(
      title = paste0("Gap vs Lambda (K=", k, ")"),
      x     = expression(lambda),
      y     = "Gap statistic",
      color = expression(m)
    ) +
    theme_minimal()
  
  plots[[i]] <- p_k
}

# To display them one after the other:
library(patchwork)
combined <- wrap_plots(plots, ncol = 2, nrow=2,guides = "collect") &
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

print(combined)

# m=1.01

df_m <- filter(gap_df, m == 1.01)

# Single plot: gap vs lambda, one curve per K
ggplot(df_m, aes(x = lambda, y = gap, group = factor(K), color = factor(K))) +
  geom_line() +
  geom_point() +
  labs(
    title = expression(paste("Gap vs ", lambda, " at ", m, " = 1.01")),
    x     = expression(lambda),
    y     = "Gap statistic",
    color = "K"
  ) +
  theme_minimal()


# soft clustering --------------------------------------------------------------------

source("Utils_fuzzyJM_2.R")

TT = 1000
P = 10
mu = 1
Sigma_rho = 0
ar_rho = 0.8
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

K_grid=2:4
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,3,length.out=6)
B          = 10
max_iter   = 5
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

# 0) Compila il C++ una volta sola
Rcpp::sourceCpp("simplex_pgd.cpp")

# 1) Definisci la griglia

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  b = 0:B)

# 2) Parallel apply
start=Sys.time()
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
end=Sys.time()
start-end
# 3) Trasforma in data.frame
results <- do.call(rbind, lapply(res_list, as.data.frame))
results <- as.data.frame(results)

gap_df <- results %>%
  # Escludiamo eventuali NA di loss
  filter(!is.na(loss)) %>%
  # Raggruppiamo per ogni combinazione di K, lambda e m
  group_by(K, lambda, m) %>%
  summarise(
    # log-loss sul dato non permutato (dovrebbe essere un singolo valore)
    logloss      = log(loss[permuted == FALSE]),
    # aspettativa del log-loss sui dati permutati
    Eg_logloss   = mean(log(loss[permuted == TRUE]), na.rm = TRUE),
    # Gap = E[log(loss_perm)] − log(loss_orig)
    gap          = Eg_logloss - logloss,
    .groups = "drop"
  )

# Dai un’occhiata
print(gap_df)

uniqueKs <- sort(unique(gap_df$K))
plots <- vector("list", length(uniqueKs))  # pre-allocate list
names(plots) <- paste0("K", uniqueKs)

for (i in seq_along(uniqueKs)) {
  k <- uniqueKs[i]
  dfk <- filter(gap_df, K == k)
  
  p_k <- ggplot(dfk, aes(x = lambda, y = gap, group = factor(m), color = factor(m))) +
    geom_line() +
    geom_point() +
    labs(
      title = paste0("Gap vs Lambda (K=", k, ")"),
      x     = expression(lambda),
      y     = "Gap statistic",
      color = expression(m)
    ) +
    theme_minimal()
  
  plots[[i]] <- p_k
}

# To display them one after the other:
library(patchwork)
combined <- wrap_plots(plots, ncol = 2, nrow=2,guides = "collect") &
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

print(combined)

# m=2.6

df_m <- filter(gap_df, m == m_grid[5])

# Single plot: gap vs lambda, one curve per K
ggplot(df_m, aes(x = lambda, y = gap, group = factor(K), color = factor(K))) +
  geom_line() +
  geom_point() +
  labs(
    title = expression(paste("Gap vs ", lambda, " at ", m, " = 2.6")),
    x     = expression(lambda),
    y     = "Gap statistic",
    color = "K"
  ) +
  theme_minimal()

K=2
m=1.5
lambda=.05
fit <- fuzzy_jump_cpp_parallel(
  Y,
  K        = K,
  lambda   = lambda,
  m        = m,
  max_iter = 10,
  n_init   = 5,
  tol      = 1e-8,
  ncores=29
)

true_S=cbind(soft_scen$pi_1,1-soft_scen$pi_1)
hellinger_distance_matrix(fit$best_S,true_S)

table(fit$MAP,soft_scen$MAP)

plot(soft_scen$pi_1,type='l',col='red',ylim=c(0,1))
lines(fit$best_S[,1],type='l',col='black')



