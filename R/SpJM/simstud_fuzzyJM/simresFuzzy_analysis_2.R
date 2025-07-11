library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

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

# K=3 soft ----------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_list_soft_K3.Rdata")

res_summary_K3_soft <- do.call(rbind, lapply(res_list_soft_K3, function(el) {
  # 1) Coerciamo S in matrice (anche se fosse un vettore 1D)
  S_mat <- as.matrix(el$S)
  
  # 2) Estraiamo ground_truth come matrice
  gt_mat <- as.matrix(el$ground_truth)
  
  # 3) Troviamo il numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # 4) Sottoselezioniamo le prime nc colonne (sempre drop=FALSE)
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  
  # 5) Calcoliamo l'MSE su tutti gli elementi
  diffs <- (S_sub - gt_sub)^2
  MSE   <- mean(diffs)
  
  # 6) Ritorniamo la riga di risultato
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = el$m,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_soft_K3_fuzzy <- res_summary_K3_soft %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_soft_K3_fuzzy <- avg_mse_soft_K3_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()

best_mse_soft_K3_fuzzy

# Competitive models

# k-prot

avg_lambda0_byTP <- avg_mse_soft_K3_fuzzy %>%
  filter(lambda == 0, m==1.01, K==3)

avg_lambda0_byTP

# Cont JM
load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_cont_soft_K3.Rdata")

res_summary_K3_soft_cont <- do.call(rbind, lapply(res_list_soft_K3_cont, function(el) {
  # Coercizione in matrice
  S_mat  <- as.matrix(el$S)
  gt_mat <- as.matrix(el$ground_truth)
  
  # Numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # Sottoselezione e calcolo MSE
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  MSE    <- mean((S_sub - gt_sub)^2)
  
  # Estrai m se presente, altrimenti NA
  m_val <- if (!is.null(el$m)) el$m else NA_real_
  
  # Ritorna la riga di output
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = m_val,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_soft_K3_cont <- res_summary_K3_soft_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_soft_K3_cont <- avg_mse_soft_K3_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()

best_mse_soft_K3_cont

# Some Diagnostics Plots

plot_data <- avg_mse_soft_K3_fuzzy %>%
  filter(K == 3, m %in% unique(hp$m)) %>%
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
ggplot(plot_data, aes(x = lambda, y = av_MSE,
                      color = m_label, group = m_label)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller) +
  scale_color_discrete(name = NULL) +
  labs(
    x = expression(lambda),
    y = "av. MSE"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Plot varying m
plot_data <- avg_hd_soft_K3_fuzzy %>%
  filter(lambda == 0.40) %>%
  filter(lambda == 0.1) %>%
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



# K=2 soft ----------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_list_soft_K2.Rdata")
res_summary_K2_soft <- do.call(rbind, lapply(res_list_soft_K2, function(el) {
  # 1) Coerciamo S in matrice (anche se fosse un vettore 1D)
  S_mat <- as.matrix(el$S)
  
  # 2) Estraiamo ground_truth come matrice
  gt_mat <- as.matrix(el$ground_truth)
  
  # 3) Troviamo il numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # 4) Sottoselezioniamo le prime nc colonne (sempre drop=FALSE)
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  
  # 5) Calcoliamo l'MSE su tutti gli elementi
  diffs <- (S_sub - gt_sub)^2
  MSE   <- mean(diffs)
  
  # 6) Ritorniamo la riga di risultato
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = el$m,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_soft_K2_fuzzy <- res_summary_K2_soft %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_soft_K2_fuzzy <- avg_mse_soft_K2_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()  
best_mse_soft_K2_fuzzy

# Competitive models
# k-prot
avg_lambda0_byTP_K2 <- avg_mse_soft_K2_fuzzy %>%
  filter(lambda == 0, m==1.01, K==2)

avg_lambda0_byTP_K2


# Cont JM
load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_cont_soft_K2.Rdata")
res_summary_K2_soft_cont <- do.call(rbind, lapply(res_list_soft_K2_cont, function(el) {
  # Coercizione in matrice
  S_mat  <- as.matrix(el$S)
  gt_mat <- as.matrix(el$ground_truth)
  
  # Numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # Sottoselezione e calcolo MSE
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  MSE    <- mean((S_sub - gt_sub)^2)
  
  # Estrai m se presente, altrimenti NA
  m_val <- if (!is.null(el$m)) el$m else NA_real_
  
  # Ritorna la riga di output
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = m_val,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))
avg_mse_soft_K2_cont <- res_summary_K2_soft_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )
best_mse_soft_K2_cont <- avg_mse_soft_K2_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()

best_mse_soft_K2_cont

# Some Diagnostics Plots

plot_data_K2 <- avg_mse_soft_K2_fuzzy %>%
  filter(K == 2, m %in% unique(hp$m)) %>%
  mutate(m_label = case_when(
    m == 1.01   ~ "m = 1.01",
    m == 1.2575 ~ "m = 1.25",
    m == 1.505  ~ "m = 1.50",
    m == 1.7525 ~ "m = 1.75",
    m == 2.00   ~ "m = 2.00"
  ))

# Create custom labeller for TT and P
custom_labeller_K2 <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data_K2, aes(x = lambda, y = av_MSE,
                          color = m_label, group = m_label)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller_K2) +
  scale_color_discrete(name = NULL) +
  labs(
    x = expression(lambda),
    y = "av. MSE"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
# Plot varying m

plot_data_K2 <- avg_mse_soft_K2_fuzzy %>%
  filter(lambda == 0.40) %>%
  mutate(K = factor(K))  # ensure K is treated as categorical for color/legend

# Custom facet labels
custom_labeller_K2 <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data_K2, aes(x = m, y = av_MSE,
                          color = K, group = K)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller_K2) +
  scale_color_discrete(name = "K") +
  labs(
    x = "m",
    y = "Mean Squared Error"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )



# K=3 hard ----------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_list_hard_K3.Rdata")

res_summary_K3_hard <- do.call(rbind, lapply(res_list_hard_K3, function(el) {
  # 1) Coerciamo S in matrice (anche se fosse un vettore 1D)
  S_mat <- as.matrix(el$S)
  
  # 2) Estraiamo ground_truth come matrice
  gt_mat <- as.matrix(el$ground_truth)
  
  # 3) Troviamo il numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # 4) Sottoselezioniamo le prime nc colonne (sempre drop=FALSE)
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  
  # 5) Calcoliamo l'MSE su tutti gli elementi
  diffs <- (S_sub - gt_sub)^2
  MSE   <- mean(diffs)
  
  # 6) Ritorniamo la riga di risultato
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = el$m,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_hard_K3_fuzzy <- res_summary_K3_hard %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_hard_K3_fuzzy <- avg_mse_hard_K3_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()

best_mse_hard_K3_fuzzy

# Competitive models
# k-prot

avg_lambda0_byTP_hard <- avg_mse_hard_K3_fuzzy %>%
  filter(lambda == 0, m==1.01, K==3)
avg_lambda0_byTP_hard

# Cont JM
load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_cont_hard_K3.Rdata")

res_summary_K3_hard_cont <- do.call(rbind, lapply(res_list_hard_K3_cont, function(el) {
  # Coercizione in matrice
  S_mat  <- as.matrix(el$S)
  gt_mat <- as.matrix(el$ground_truth)
  
  # Numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # Sottoselezione e calcolo MSE
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  MSE    <- mean((S_sub - gt_sub)^2)
  
  # Estrai m se presente, altrimenti NA
  m_val <- if (!is.null(el$m)) el$m else NA_real_
  
  # Ritorna la riga di output
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = m_val,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_hard_K3_cont <- res_summary_K3_hard_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_hard_K3_cont <- avg_mse_hard_K3_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()

best_mse_hard_K3_cont

# Some Diagnostics Plots
plot_data_hard <- avg_mse_hard_K3_fuzzy %>%
  filter(K == 3, m %in% unique(hp$m)) %>%
  mutate(m_label = case_when(
    m == 1.01   ~ "m = 1.01",
    m == 1.2575 ~ "m = 1.25",
    m == 1.505  ~ "m = 1.50",
    m == 1.7525 ~ "m = 1.75",
    m == 2.00   ~ "m = 2.00"
  ))

# Create custom labeller for TT and P
custom_labeller_hard <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data_hard, aes(x = lambda, y = av_MSE,
                            color = m_label, group = m_label)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller_hard) +
  scale_color_discrete(name = NULL) +
  labs(
    x = expression(lambda),
    y = "av. MSE"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
# Plot varying m

plot_data_hard <- avg_mse_hard_K3_fuzzy %>%
  filter(lambda == 0.1) %>%
  mutate(K = factor(K))  # ensure K is treated as categorical for color/legend

# Custom facet labels
custom_labeller_hard <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data_hard, aes(x = m, y = av_MSE,
                            color = K, group = K)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller_hard) +
  scale_color_discrete(name = "K") +
  labs(
    x = "m",
    y = "Mean Squared Error"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
# K=2 hard ----------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_list_hard_K2.Rdata")

res_summary_K2_hard <- do.call(rbind, lapply(res_list_hard_K2, function(el) {
  # 1) Coerciamo S in matrice (anche se fosse un vettore 1D)
  S_mat <- as.matrix(el$S)
  
  # 2) Estraiamo ground_truth come matrice
  gt_mat <- as.matrix(el$ground_truth)
  
  # 3) Troviamo il numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # 4) Sottoselezioniamo le prime nc colonne (sempre drop=FALSE)
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  
  # 5) Calcoliamo l'MSE su tutti gli elementi
  diffs <- (S_sub - gt_sub)^2
  MSE   <- mean(diffs)
  
  # 6) Ritorniamo la riga di risultato
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = el$m,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_hard_K2_fuzzy <- res_summary_K2_hard %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_hard_K2_fuzzy <- avg_mse_hard_K2_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()

best_mse_hard_K2_fuzzy


# Competitive models
# k-prot
avg_lambda0_byTP_hard_K2 <- avg_mse_hard_K2_fuzzy %>%
  filter(lambda == 0, m==1.01, K==2)
avg_lambda0_byTP_hard_K2

# Cont JM
load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/res_cont_hard_K2.Rdata")

res_summary_K2_hard_cont <- do.call(rbind, lapply(res_list_hard_K2_cont, function(el) {
  # Coercizione in matrice
  S_mat  <- as.matrix(el$S)
  gt_mat <- as.matrix(el$ground_truth)
  
  # Numero di colonne comuni
  nc <- min(ncol(S_mat), ncol(gt_mat))
  
  # Sottoselezione e calcolo MSE
  S_sub  <- S_mat[, 1:nc, drop = FALSE]
  gt_sub <- gt_mat[, 1:nc, drop = FALSE]
  MSE    <- mean((S_sub - gt_sub)^2)
  
  # Estrai m se presente, altrimenti NA
  m_val <- if (!is.null(el$m)) el$m else NA_real_
  
  # Ritorna la riga di output
  data.frame(
    K      = el$K,
    lambda = el$lambda,
    m      = m_val,
    TT     = el$TT,
    P      = el$P,
    seed   = el$seed,
    MSE    = MSE
  )
}))

avg_mse_hard_K2_cont <- res_summary_K2_hard_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    av_MSE = mean(MSE),
    sd_MSE = sd(MSE),
    .groups = "drop"
  )

best_mse_hard_K2_cont <- avg_mse_hard_K2_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = av_MSE, n = 1, with_ties = FALSE) %>%
  ungroup()
best_mse_hard_K2_cont

# Some Diagnostics Plots
plot_data_hard_K2 <- avg_mse_hard_K2_fuzzy %>%
  filter(K == 2, m %in% unique(hp$m)) %>%
  mutate(m_label = case_when(
    m == 1.01   ~ "m = 1.01",
    m == 1.2575 ~ "m = 1.25",
    m == 1.505  ~ "m = 1.50",
    m == 1.7525 ~ "m = 1.75",
    m == 2.00   ~ "m = 2.00"
  ))

# Create custom labeller for TT and P
custom_labeller_hard_K2 <- labeller(
  TT = function(x) paste("T =", x),
  P  = function(x) paste("P =", x)
)

# Plot
ggplot(plot_data_hard_K2, aes(x = lambda, y = av_MSE,
                               color = m_label, group = m_label)) +
  geom_line(size = .7) +
  facet_grid(TT ~ P, labeller = custom_labeller_hard_K2) +
  scale_color_discrete(name = NULL) +
  labs(
    x = expression(lambda),
    y = "av. MSE"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
