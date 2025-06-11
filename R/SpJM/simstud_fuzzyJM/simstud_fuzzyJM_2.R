library(parallel)
library(Rcpp)
library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_fuzzyJM_2.R")

# K=2 soft ---------------------------------------------------------------------

max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

K_grid=2:3
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,2,length.out=5)
seed=1:50

TT=c(1000,2000)
P=c(5,10)

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

mu=1
Sigma_rho=0
ar_rho = 0.9
tau = .2
ncores=detectCores()-1

start=Sys.time()
res_list_soft_K2 <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  mi    <- hp$m[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  soft_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=cbind(soft_scen$pi_1,1-soft_scen$pi_1)
  
  Yinput=soft_scen[,(1:Pi)+1]
  
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
      list(success = TRUE, loss = fit$loss,
           PE=fit$PE,PB=fit$PB,PB_lambda=fit$PB_lambda,XB=fit$XB,
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
      last_loss <- result$loss
      PE <- result$PE
      PB <- result$PB
      PB_lambda <- result$PB_lambda
      XB <- result$XB
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
    S=S,
    ground_truth=ground_truth,
    K         = Ki,
    lambda    = li,
    m         = mi,
    loss      = last_loss,
    attempts  = attempt,
    PE=PE,
    PB=PB,
    PB_lambda=PB_lambda,
    XB=XB
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

lista_risultati <- lapply(res_list_soft_K2, function(el) {
  # estraggo S e ground_truth come vettori
  cS  <- el$S
  cGT <- el$ground_truth
  
  # calcolo la distanza di Hellinger
  hd <- hellinger_distance_matrix(cS, cGT)
  
  # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
  data.frame(
    jhellinger_dist = hd,
    lambda          = el$lambda,
    K               = el$K,
    m               = el$m,
    loss            = el$loss,
    attempts        = el$attempts,
    PE              = el$PE,
    PB              = el$PB,
    PB_lambda       = el$PB_lambda,
    XB              = el$XB,
    stringsAsFactors = FALSE
  )
})

# 3. Combino tutti i data.frame in uno solo
results_df_soft_K2_fuzzy <- do.call(rbind, lista_risultati)
results_df_soft_K2_fuzzy$TT=hp$TT
results_df_soft_K2_fuzzy$P=hp$P
results_df_soft_K2_fuzzy$seed=hp$seed

#save(results_df_soft_K2_fuzzy,file='hellinger_df_soft_K2_fuzzy.Rdata')

# load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_soft_K2_fuzzy.Rdata")

head(results_df_soft_K2_fuzzy)

avg_hd_soft_K2_fuzzy <- results_df_soft_K2_fuzzy %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

# Best hellinger:

best_hd_soft_K2_fuzzy <- avg_hd_soft_K2_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

# k-prot
avg_lambda0_byTP <- avg_hd_soft_K2_fuzzy %>%
  filter(lambda == 0, m==1.01, K==2)


# Plot varying lambda
plot_data <- avg_hd_soft_K2_fuzzy %>%
  filter(K == 2, m %in% unique(hp$m)) %>%
  mutate(m_label = case_when(
    m == 1.01   ~ "m = 1.01",
    m == 1.2575 ~ "m = 1.25",
    m == 1.505  ~ "m = 1.50",
    m == 1.7525 ~ "m = 1.75",
    m == 2.00   ~ "m = 2.00"
  ))

# Create custom labeller for TT and P
custom_labeller <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data, aes(x = lambda, y = mean_hellinger,
                      color = m_label, group = m_label)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller) +
  scale_color_discrete(name = NULL) +
  labs(
    x = expression(lambda),
    y = "Mean Hellinger Distance"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Plot varying m
plot_data <- avg_hd_soft_K2_fuzzy %>%
  filter(lambda == 0.20) %>%
  mutate(K = factor(K))  # ensure K is treated as categorical for color/legend

# Custom facet labels
custom_labeller <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data, aes(x = m, y = mean_hellinger,
                      color = K, group = K)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller) +
  scale_color_discrete(name = "K") +
  labs(
    x = "m",
    y = "Mean Hellinger Distance"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# K=2 hard ----------------------------------------------------------------

max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

K_grid=2:3
lambda_grid=seq(0,.5,length.out=11)
m_grid=seq(1.01,2,length.out=5)
seed=1:50

TT=c(1000,2000)
P=c(5,10)

hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  m = m_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

mu=1
Sigma_rho=0
ar_rho = 0.99
tau = 5
ncores=detectCores()-1

start=Sys.time()
res_list_hard_K2 <- mclapply(seq_len(nrow(hp)), function(i) {
  Ki    <- hp$K[i]
  li    <- hp$lambda[i]
  mi    <- hp$m[i]
  TTi    <- hp$TT[i]
  Pi    <- hp$P[i]
  seedi <- hp$seed[i]
  
  hard_scen=simulate_fuzzy_mixture_mv(
    TT = TTi,
    P = Pi,
    mu = mu,
    Sigma_rho = Sigma_rho,
    ar_rho = ar_rho,
    tau = tau,
    seed = seedi
  )
  
  ground_truth=cbind(hard_scen$pi_1,1-hard_scen$pi_1)
  
  Yinput=hard_scen[,(1:Pi)+1]
  
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
      list(success = TRUE, loss = fit$loss,
           PE=fit$PE,PB=fit$PB,PB_lambda=fit$PB_lambda,XB=fit$XB,
           S=fit$best_S)
    }, error = function(e) {
      message(sprintf("Row %d, attempt %d/%d failed: %s",
                      i, attempt, max_retries, e$message))
      list(success = FALSE, loss = NA_real_)
    })
    
    if (result$success) {
      S= result$S
      last_loss <- result$loss
      PE <- result$PE
      PB <- result$PB
      PB_lambda <- result$PB_lambda
      XB <- result$XB
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
    S=S,
    ground_truth=ground_truth,
    K         = Ki,
    lambda    = li,
    m         = mi,
    TT = TTi,
    P = Pi,
    seed= seedi,
    mu=mu, 
    Sigma_rho=Sigma_rho,
    ar_rho=ar_rho,
    tau=tau,
    phi=phi,
    loss      = last_loss,
    attempts  = attempt,
    PE=PE,
    PB=PB,
    PB_lambda=PB_lambda,
    XB=XB
  )
}, mc.cores = ncores)
end=Sys.time()
end-start

lista_risultati <- lapply(res_list_hard_K2, function(el) {
  # estraggo S e ground_truth come vettori
  cS  <- el$S
  cGT <- el$ground_truth
  
  # calcolo la distanza di Hellinger
  hd <- hellinger_distance_matrix(cS, cGT)
  
  # costruisco un piccolo data.frame con jhellinger_dist e tutti gli altri campi
  data.frame(
    jhellinger_dist = hd,
    lambda          = el$lambda,
    K               = el$K,
    m               = el$m,
    loss            = el$loss,
    attempts        = el$attempts,
    PE              = el$PE,
    PB              = el$PB,
    PB_lambda       = el$PB_lambda,
    XB              = el$XB,
    stringsAsFactors = FALSE
  )
})

# 3. Combino tutti i data.frame in uno solo
results_df_hard_K2_fuzzy <- do.call(rbind, lista_risultati)

results_df_hard_K2_fuzzy$TT=hp$TT
results_df_hard_K2_fuzzy$P=hp$P
results_df_hard_K2_fuzzy$seed=hp$seed

save(results_df_hard_K2_fuzzy,file='hellinger_df_hard_K2_fuzzy.Rdata')

head(results_df_hard_K2_fuzzy)

avg_hd_hard_K2_fuzzy <- results_df_hard_K2_fuzzy %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

# Best hellinger:

best_hd_hard_K2_fuzzy <- avg_hd_hard_K2_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

# k-prot
avg_lambda0_byTP <- avg_hd_hard_K2_fuzzy %>%
  filter(lambda == 0, m==1.01, K==2) 


# Plot varying lambda
plot_data <- avg_hd_hard_K2_fuzzy %>%
  filter(K == 2, m %in% unique(hp$m)) %>%
  mutate(m_label = case_when(
    m == 1.01   ~ "m = 1.01",
    m == 1.2575 ~ "m = 1.25",
    m == 1.505  ~ "m = 1.50",
    m == 1.7525 ~ "m = 1.75",
    m == 2.00   ~ "m = 2.00"
  ))

# Create custom labeller for TT and P
custom_labeller <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data, aes(x = lambda, y = mean_hellinger,
                      color = m_label, group = m_label)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller) +
  scale_color_discrete(name = NULL) +
  labs(
    x = expression(lambda),
    y = "Mean Hellinger Distance"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Plot varying m
plot_data <- avg_hd_hard_K2_fuzzy %>%
  filter(lambda == 0.05) %>%
  mutate(K = factor(K))  # ensure K is treated as categorical for color/legend

# Custom facet labels
custom_labeller <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data, aes(x = m, y = mean_hellinger,
                      color = K, group = K)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller) +
  scale_color_discrete(name = "K") +
  labs(
    x = "m",
    y = "Mean Hellinger Distance"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
