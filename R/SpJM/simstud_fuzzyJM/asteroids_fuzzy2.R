trans_theta=function(theta){
  #theta_trans <- data$theta
  theta[which(theta > pi)] <- theta[which(theta > pi)] - 2 * pi
  #data$theta <- theta_trans
  return(theta)
}

max_min_feat=function(data,tt_thres_maxmin=2.5,l=5){
  library(dplyr)
  maxs <- as.numeric(splus2R::peaks(data$theta, span = l))
  mins <- as.numeric(splus2R::peaks(-data$theta, span = l))
  
  # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
  maxs[which(data$theta>(tt_thres_maxmin))] <- 0
  mins[which(data$theta<(-tt_thres_maxmin))] <- 0
  
  #
  data$maxs=maxs
  data$mins=mins
  
  find_closest <- function(t_val, t_events, theta_events) {
    if (length(t_events) == 0) return(NA)  # No events found
    diffs <- abs(t_events - t_val)  # Compute absolute time differences
    closest_index <- which.min(diffs)  # Find index of closest event
    return(theta_events[closest_index])  # Return corresponding theta value
  }
  
  # Extract time and values of local maxima and minima
  t_maxs <- data$t[data$maxs == 1]
  theta_maxs <- data$theta[data$maxs == 1]
  
  t_mins <- data$t[data$mins == 1]
  theta_mins <- data$theta[data$mins == 1]
  
  # Compute the closest local max and min for each row
  data <- data %>%
    mutate(
      value_max = sapply(t, find_closest, t_maxs, theta_maxs),
      value_min = sapply(t, find_closest, t_mins, theta_mins)
    )
  
  # TP indicator
  data$I_TP=as.numeric(data$value_max*data$value_min>0 & data$value_max>data$value_min)
  
  # HS indicator
  data$I_HS=as.numeric(data$value_max*data$value_min<0&data$value_max<data$value_min)
  
  # QS indicator
  data$I_QS=as.numeric(data$value_max*data$value_min<0 & data$value_max>data$value_min)
  
  # CP indicator
  data$I_CP=as.numeric(data$value_max*data$value_min>0&data$value_max<=data$value_min)
  
  # Check if distance between max and min is small to distinguish between NR and QS
  # data$I_QS=data$I_QS*data$I_diffmaxmin
  
  
  
  data$mean_osc=data$I_HS*(data$value_min+data$value_max+2*pi)/2+
    data$I_QS*(data$value_min+data$value_max)/2+
    data$I_TP*(data$value_min+data$value_max)/2
  
  data=data[,c("t","maxs","mins",
               "value_min","value_max",
               #"I_TP",
               "I_HS","I_QS",
               #"I_CP",
               "mean_osc")]
}


load("propagation_164207_new_v2.Rdata")



df164207$theta=trans_theta(df164207$theta)

df164207$dtheta=c(NA,diff(df164207$theta))
df164207$de=c(NA,diff(df164207$e))
df164207$domega=c(NA,diff(df164207$omega))

df164207_maxmin_theta=max_min_feat(df164207,
                                   tt_thres_maxmin=2.5,l=5)

df164207=merge(df164207_maxmin_theta,
               df164207,by="t")

df164207$I_domega=(df164207$domega<=-0.0025)


df164207_select=df164207[,c("t","value_min","value_max",
                            "I_HS","I_QS",
                            "mean_osc",
                            "I_domega","type")]

df164207_select=df164207[,c("t",
                            "theta",
                            "omega",
                            #"domega",
                            "value_min","value_max",
                            #"I_HS","I_QS",
                            #"mean_osc",
                            "I_domega",
                            "type")]

str(df164207_select)
# df164207_select$I_HS=as.factor(df164207_select$I_HS)
# df164207_select$I_QS=as.factor(df164207_select$I_QS)
df164207_select$I_domega=as.factor(as.numeric(df164207_select$I_domega))

str(df164207_select)
df164207_select=df164207_select[complete.cases(df164207_select),]

Y=subset(df164207_select,select=-c(t,type))

str(Y)

source("Utils_fuzzyJM_2.R")
lambda=.5
m=1.5

start=Sys.time()
fit_164207=fuzzy_jump_cpp(Y=Y, 
                          K=2, 
                          lambda = lambda, 
                          m = m,
                          max_iter = 5, 
                          n_init = 10, 
                          tol = 1e-8, 
                          verbose = T
)
end=Sys.time()
end-start

results_164207=data.frame(df164207_select,MAP=fit_164207$MAP,
                           p1=fit_164207$best_S[,1])

tapply(results_164207$theta, results_164207$MAP, mean)
tapply(results_164207$theta, results_164207$MAP, sd)

tapply(results_164207$omega, results_164207$MAP, mean)
tapply(results_164207$omega, results_164207$MAP, sd)


library(ggplot2)
library(tidyverse)

# 1) Pivot θ and ω into long format
df_long2 <- results_164207 %>%
  select(time, theta, omega, p1) %>%
  pivot_longer(
    cols      = c(theta, omega),
    names_to  = "variable",
    values_to = "value"
  )

ggplot(df_long2, aes(x = time, y = value, color = p1, group = 1)) +
  geom_point(size = 1) +
  facet_wrap(
    ~ variable,
    ncol           = 1,
    scales         = "free_y",
    labeller       = label_parsed,
    strip.position = "left"
  ) +
  scale_color_gradientn(
    colors = c("cyan", "yellow", "magenta"),
    limits = c(0, 1),
    name   = expression(P(QS))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement      = "outside",
    strip.text.y.left    = element_text(angle = 90, face = "italic", size = 16),
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    #legend.position      = c(0.85, 0.20),               # inside plot
    legend.background    = element_rect(fill = "white", # grey box
                                        color = "grey50"),
    legend.key           = element_rect(fill = "white", color = NA)
  ) +
  scale_x_continuous()
