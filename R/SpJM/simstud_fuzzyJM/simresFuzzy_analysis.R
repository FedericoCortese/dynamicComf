library(dplyr)
library(tidyr)
library(ggplot2)

<<<<<<< HEAD
=======
max_iter   = 10
n_init     = 10
tol        = 1e-8
verbose    = FALSE
ncores   = 25
max_retries=50

>>>>>>> 705082e2ca4af0e77243b8cfc65f4d4157cb14b6
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

# K=2 soft ----------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_soft_K2_fuzzy.Rdata")

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

best_hd_soft_K2_fuzzy

# k-prot
avg_lambda0_byTP <- avg_hd_soft_K2_fuzzy %>%
  filter(lambda == 0, m==1.01, K==2)

avg_lambda0_byTP

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

# Cont JM

head(results_df_soft_K2_cont)

avg_hd_soft_K2_cont <- results_df_soft_K2_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

best_hd_soft_K2_fuzzy <- avg_hd_soft_K2_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

best_hd_soft_K2_fuzzy

# K=2 hard ----------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_hard_K2_fuzzy.Rdata")


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

best_hd_hard_K2_fuzzy

# k-prot
avg_lambda0_byTP <- avg_hd_hard_K2_fuzzy %>%
  filter(lambda == 0, m==1.01, K==2) 

avg_lambda0_byTP

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

# Cont JM

head(results_df_hard_K2_cont)

avg_hd_hard_K2_cont <- results_df_hard_K2_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

best_hd_hard_K2_cont <- avg_hd_hard_K2_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

best_hd_hard_K2_cont


# K=3 soft ----------------------------------------------------------------

<<<<<<< HEAD
#load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_soft_K3_fuzzy.Rdata")
=======
load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_soft_K3_fuzzy.Rdata")

>>>>>>> 705082e2ca4af0e77243b8cfc65f4d4157cb14b6
head(results_df_soft_K3_fuzzy)

avg_hd_soft_K3_fuzzy <- results_df_soft_K3_fuzzy %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

# Best hellinger:

best_hd_soft_K3_fuzzy <- avg_hd_soft_K3_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

best_hd_soft_K3_fuzzy

# k-prot
avg_lambda0_byTP <- avg_hd_soft_K3_fuzzy %>%
  filter(lambda == 0, m==1.01, K==3)

avg_lambda0_byTP

# Plot varying lambda
plot_data <- avg_hd_soft_K3_fuzzy %>%
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
plot_data <- avg_hd_soft_K3_fuzzy %>%
<<<<<<< HEAD
  filter(lambda == 0.40) %>%
=======
  filter(lambda == 0.1) %>%
>>>>>>> 705082e2ca4af0e77243b8cfc65f4d4157cb14b6
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

# Cont JM

# load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_soft_K3_cont.Rdata")

head(results_df_soft_K3_cont)
K_grid=2:3
lambda_grid=c(0,.5,1,5,10,50,100,1000)
seed=1:50

TT=c(1000,2000)
P=c(5,10)


hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

results_df_soft_K3_cont$TT=hp$TT
results_df_soft_K3_cont$P=hp$P

avg_hd_soft_K3_cont <- results_df_soft_K3_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

best_hd_soft_K3_fuzzy <- avg_hd_soft_K3_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

best_hd_soft_K3_fuzzy


# K=3 hard ----------------------------------------------------------------

<<<<<<< HEAD
#load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_hard_K3_fuzzy.Rdata")
=======
load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_hard_K3_fuzzy.Rdata")

>>>>>>> 705082e2ca4af0e77243b8cfc65f4d4157cb14b6
head(results_df_hard_K3_fuzzy)

avg_hd_hard_K3_fuzzy <- results_df_hard_K3_fuzzy %>%
  group_by(TT, P, K, lambda, m) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

# Best hellinger:

best_hd_hard_K3_fuzzy <- avg_hd_hard_K3_fuzzy %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

best_hd_hard_K3_fuzzy

# k-prot
avg_lambda0_byTP <- avg_hd_hard_K3_fuzzy %>%
  filter(lambda == 0, m==1.01, K==3)

avg_lambda0_byTP

# Plot varying lambda
plot_data <- avg_hd_hard_K3_fuzzy %>%
<<<<<<< HEAD
  filter(K == 2, m %in% unique(results_df_hard_K3_fuzzy$m)) %>%
=======
  filter(K == 3, m %in% unique(hp$m)) %>%
>>>>>>> 705082e2ca4af0e77243b8cfc65f4d4157cb14b6
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
plot_data <- avg_hd_hard_K3_fuzzy %>%
  filter(lambda == 0.10) %>%
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

# Cont JM

# load("D:/CNR/OneDrive - CNR/Comfort - HMM/simres_fuzzyJM_fuzzySTJM/hellinger_df_hard_K3_cont.Rdata")

head(results_df_hard_K3_cont)
K_grid=2:3
lambda_grid=c(0,.5,1,5,10,50,100,1000)
seed=1:50

TT=c(1000,2000)
P=c(5,10)


hp <- expand.grid(K = K_grid,
                  lambda = lambda_grid,
                  TT=TT,
                  P=P,
                  seed=seed)

results_df_hard_K3_cont$TT=hp$TT
results_df_hard_K3_cont$P=hp$P

avg_hd_hard_K3_cont <- results_df_hard_K3_cont %>%
  group_by(TT, P, K, lambda) %>%
  summarize(
    mean_hellinger = mean(jhellinger_dist),
    sd_hellinger = sd(jhellinger_dist),
    .groups = "drop"
  )

best_hd_hard_K3_fuzzy <- avg_hd_hard_K3_cont %>%
  group_by(TT, P) %>%
  slice_min(order_by = mean_hellinger, n = 1, with_ties = FALSE) %>%
  ungroup()

best_hd_hard_K3_fuzzy
