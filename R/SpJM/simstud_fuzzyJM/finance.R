
# Prices data --------------------------------------------------------------

# install.packages("quantmod")
library(quantmod)

start_date <- "2019-01-01"

# 1. Stocks & bonds
getSymbols(c("SPY", "AGG"), src = "yahoo", from = start_date)
spy <- Cl(SPY)
agg <- Cl(AGG)

# 2. Gold (proxy: GLD ETF)
getSymbols("GLD", src = "yahoo", from = start_date)
gold <- Cl(GLD)

# 3. Bitcoin (USD)
btc  <- Cl(getSymbols("BTC-USD",
                      src        = "yahoo",
                      from       = start_date,
                      auto.assign=FALSE))

# 4. EUR–USD FX rate
eur_usd <- Cl(getSymbols("EURUSD=X",
                         src        = "yahoo",
                         from       = start_date,
                         auto.assign=FALSE))

# 5. Merge and drop NAs
all_data <- na.omit(
  merge(spy, agg, gold, btc, eur_usd)
)

# Plot each in a different panel

library(ggplot2)
library(tidyr)
all_data_long <- pivot_longer(
  data = data.frame(date = index(all_data), coredata(all_data)),
  cols = -date,
  names_to = "asset",
  values_to = "price"
)

ggplot(all_data_long, aes(x = date, y = price)) +
  geom_line() +
  facet_wrap(~ asset, scales = "free_y") +
  labs(title = "Asset Prices Over Time",
       x = "Date",
       y = "Price") +
  theme_minimal()

# 6. Inspect
head(all_data)
cor(all_data)

# Corrplot of prices
library(corrplot)
cor_matrix_prices <- cor(all_data, use = "pairwise.complete.obs")
corrplot(cor_matrix_prices, method = "circle", type = "upper",
         title = "Correlation of Asset Prices",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)


# Log-returns -------------------------------------------------------------


# 1) compute daily log‐returns on the merged price series
returns_xts <- na.omit(diff(log(all_data)))

# 2) turn the xts into a data.frame with a date column
returns_df <- data.frame(
  date = index(returns_xts),
  coredata(returns_xts)
)

# 3) (optional) rename columns to indicate returns
colnames(returns_df)[-1] <- paste0(colnames(all_data), "_ret")

# 4) inspect
head(returns_df)

# Plot each in a different panel
library(ggplot2)
library(tidyr)
returns_long <- pivot_longer(
  returns_df,
  cols = -date,
  names_to = "asset",
  values_to = "return"
)

ggplot(returns_long, aes(x = date, y = return)) +
  geom_line() +
  facet_wrap(~ asset, scales = "free_y") +
  labs(title = "Daily Returns of Assets",
       x = "Date",
       y = "Return") +
  theme_minimal()

# 5) Correlation plot of the returns
library(corrplot)
cor_matrix <- cor(returns_xts, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "circle", type = "upper",
         title = "Correlation of Daily Returns",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)


# Additional features -----------------------------------------------------

# 7-days Moving average of returns
library(dplyr)
Y=returns_df

# Moving standard deviations

Y$sd_SPY.Close_ret=zoo::rollapply(Y$SPY.Close_ret, 7, sd, fill=NA)
Y$sd_AGG.Close_ret=zoo::rollapply(Y$AGG.Close_ret, 7, sd, fill=NA)
Y$sd_GLD.Close_ret=zoo::rollapply(Y$GLD.Close_ret, 7, sd, fill=NA)
Y$sd_BTC.USD.Close_ret=zoo::rollapply(Y$BTC.USD.Close_ret, 7, sd, fill=NA)
Y$sd_EURUSD.X.Close_ret=zoo::rollapply(Y$EURUSD.X.Close_ret, 7, sd, fill=NA)

# Plot the moving standard deviations

ggplot(Y_long, aes(x = date, y = return)) +
  geom_line() +
  facet_wrap(~ asset, scales = "free_y") +
  labs(title = "7-Day Moving Standard Deviations of Daily Returns",
       x = "Date",
       y = "Standard Deviation") +
  theme_minimal()


Y=Y[complete.cases(Y),]

dim(Y)


# Comparison --------------------------------------------------------------

# 1) Define the sequence of lambda values
lambdas <- seq(0, 1, by = 0.1)

# 2) Fit the model for each lambda via lapply
fits <- lapply(lambdas, function(i) {
  fuzzy_jump_cpp(
    Y         = Y[, -1],
    K         = 2,
    lambda    = i,
    m         = 1.1,
    max_iter  = 10,
    n_init    = 5,
    tol       = 1e-8,
    verbose   = FALSE,
    parallel  = FALSE,
    n_cores   = NULL
  )
})
names(fits) <- paste0("λ=", lambdas)

# 3) Compare each consecutive pair via Map (no explicit SIMPLIFY here)
comparison_list <- Map(
  function(f1, f2, λ1, λ2) {
    MAP1    <- f1$MAP
    MAP2    <- f2$MAP
    S1      <- f1$best_S
    S2      <- f2$best_S
    data.frame(
      lambda1 = λ1,
      lambda2 = λ2,
      ARI      = adjustedRandIndex(MAP1, MAP2),
      avgMSE   = mean((S1 - S2)^2)
    )
  },
  fits[-length(fits)],    # all but last fit
  fits[-1],               # all but first fit
  lambdas[-length(lambdas)],
  lambdas[-1]
)

# 4) Combine into one data.frame
results <- do.call(rbind, comparison_list)
print(results)

ggplot(results, aes(x = lambda1, y = avgMSE)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = expression(lambda),
    y = "avg. MSE"
  ) +
  theme_minimal(base_size = 16) +                # increases all text
  theme(
    axis.title   = element_text(size = 18),      # axis labels
    axis.text    = element_text(size = 14),      # tick labels
    plot.title   = element_text(size = 20),      # if you add a title
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Final Fit ---------------------------------------------------------------------



source("Utils_fuzzyJM_2.R")
fit_finance = fuzzy_jump_cpp(Y=Y[,-1], 
                               K=2, 
                               lambda = .5, 
                               m = 1.1,
                               max_iter = 10, 
                               n_init = 5, 
                               tol = 1e-8, 
                               verbose = T,
                               parallel = FALSE,
                               n_cores = NULL
)

results_finance=data.frame(all_data[-(1:7)],Y[,-1],MAP=fit_finance$MAP,
                           p1=fit_finance$best_S[,1])

tapply(results_finance$AGG.Close_ret, results_finance$MAP, mean)*100
tapply(results_finance$AGG.Close_ret, results_finance$MAP, sd)*100

tapply(results_finance$BTC.USD.Close_ret, results_finance$MAP, mean)*100
tapply(results_finance$BTC.USD.Close_ret, results_finance$MAP, sd)*100

tapply(results_finance$EURUSD.X.Close_ret, results_finance$MAP, mean)*100
tapply(results_finance$EURUSD.X.Close_ret, results_finance$MAP, sd)*100

tapply(results_finance$GLD.Close_ret, results_finance$MAP, mean)*100
tapply(results_finance$GLD.Close_ret, results_finance$MAP, sd)*100

tapply(results_finance$SPY.Close_ret, results_finance$MAP, mean)*100
tapply(results_finance$SPY.Close_ret, results_finance$MAP, sd)*100



head(results_finance)

library(tidyverse)

# 1. Convert rownames to a date column
df <- results_finance %>%
  rownames_to_column(var = "date") %>%
  mutate(date = as.Date(date))

# 2. Pivot the five price columns into long format, keeping p1
df_long <- df %>%
  select(date, p1, SPY.Close:EURUSD.X.Close) %>%
  pivot_longer(
    cols      = -c(date, p1),
    names_to  = "asset",
    values_to = "price"
  )

# 3. Plot with lines colored by p1
ggplot(df_long, aes(x = date, y = price, color = p1, group = 1)) +
  geom_line(size = 1) +
  facet_wrap(~ asset, ncol = 3, scales = "free_y") +
  scale_color_gradientn(
    colors = c("green","yellow", "red"),
    limits = c(0, 1),
    name   = "P(Bear)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.title.x     = element_blank()
  ) +
  labs(
    y     = "Closing Price"
    # ,
    # title = "Daily Prices"
  )



library(tidyverse)

# prepare df_long as before…

ggplot(df_long, aes(date, price, color = p1, group = 1)) +
  geom_line(size = 1) +
  facet_wrap(~ asset, ncol = 3, scales = "free_y") +
  scale_color_gradientn(
    colors   = c("green","yellow","red"),
    limits   = c(0, 1),
    name     = "P(Bear)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text        = element_text(face = "bold"),
    axis.title.x      = element_blank(),
    legend.position   = c(0.9, 0.1),      # x,y in [0,1], bottom‑right
    legend.justification = c("right","bottom"),
    legend.background = element_rect(fill = "white", color = "grey50")
  ) +
  labs(y = "Closing Price")
