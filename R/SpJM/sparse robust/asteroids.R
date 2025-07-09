inv_trans_theta=function(data){
  theta_inv=data$theta
  theta_inv[which(theta_inv<0)]=theta_inv[which(theta_inv<0)]+2*pi
  data$theta <- theta_inv
  return(data)
}
trans_theta=function(data){
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta <- theta_trans
  return(data)
}

max_min_feat=function(x,tt_thres_maxmin=2.5,l=5){
  library(dplyr)
  maxs <- as.numeric(splus2R::peaks(x$theta, span = l))
  mins <- as.numeric(splus2R::peaks(-x$theta, span = l))
  
  # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
  maxs[which(x$theta>(tt_thres_maxmin))] <- 0
  mins[which(x$theta<(-tt_thres_maxmin))] <- 0
  
  #
  x$maxs=maxs
  x$mins=mins
  
  find_closest <- function(t_val, t_events, theta_events) {
    if (length(t_events) == 0) return(NA)  # No events found
    diffs <- abs(t_events - t_val)  # Compute absolute time differences
    closest_index <- which.min(diffs)  # Find index of closest event
    return(theta_events[closest_index])  # Return corresponding theta value
  }
  
  # Extract time and values of local maxima and minima
  t_maxs <- x$t[x$maxs == 1]
  theta_maxs <- x$theta[x$maxs == 1]
  
  t_mins <- x$t[x$mins == 1]
  theta_mins <- x$theta[x$mins == 1]
  
  # Compute the closest local max and min for each row
  x <- x %>%
    mutate(
      value_max = sapply(t, find_closest, t_maxs, theta_maxs),
      value_min = sapply(t, find_closest, t_mins, theta_mins)
    )
  
  # TP indicator
  x$I_TP=as.numeric(x$value_max*x$value_min>0 & x$value_max>x$value_min)
  
  # HS indicator
  x$I_HS=as.numeric(x$value_max*x$value_min<0&x$value_max<x$value_min)
  
  # QS indicator
  x$I_QS=as.numeric(x$value_max*x$value_min<0 & x$value_max>x$value_min)
  
  # CP indicator
  x$I_CP=as.numeric(x$value_max*x$value_min>0&x$value_max<=x$value_min)
  
  
  x$mean_osc=x$I_HS*(x$value_min+x$value_max+2*pi)/2+
    x$I_QS*(x$value_min+x$value_max)/2+
    x$I_TP*(x$value_min+x$value_max)/2
  
  x=subset(x,select=-c(maxs,mins))
  
  return(x)
}

# a_feat=function(data,l){
#   cust_fun=function(x){
#     all(x==T)
#   }
#   
#   ind_a=I(data$a<1)
#   ind_a2=I(data$a>1)
#   
#   # Short
#   ind_a_mov <- runner::runner(ind_a, k = l, f = cust_fun, na_pad = TRUE)
#   ind_a2_mov <- runner::runner(ind_a2, k = l, f = cust_fun, na_pad = TRUE)
#   data$I_a=as.numeric(ind_a_mov+ind_a2_mov)
#   
#   return(data)
#   
# }

compute_feat=function(x,l_short,l_long,tt_thres_maxmin=2.5){
  
  # Transform theta
  x=trans_theta(x)
  
  # First diff
  x$dtheta=c(NA,diff(x$theta))
  
  # Moving sd of theta
  x$sd_theta=zoo::rollapply(x$theta, l_long, sd, fill=NA)
  x$sd_dtheta=zoo::rollapply(x$dtheta, l_long, sd, fill=NA)
  
  # Moving av and sd of a
  x$ma_a=zoo::rollapply(x$a, l_long, mean, fill=NA)
  x$sd_a=zoo::rollapply(x$a, l_long, sd, fill=NA)
  
  x=max_min_feat(x,tt_thres_maxmin=tt_thres_maxmin,l=l_short)
  
  return(x)
  
}

# propagation_2001GO2_new_v2

df_2001GO2=read.table("data_asteroids/propagation_2001GO2_new_v2.txt",header = T)
unique(df_2001GO2$type)

# propagation_2002AA29_new_v2
df_2002AA29=read.table("data_asteroids/propagation_2002AA29_new_v2.txt",header = T)
df_2002AA29=trans_theta(df_2002AA29)
unique(df_2002AA29$type)
plot(df_2002AA29$theta,type='p',col=df_2002AA29$type+1)


# 2014OL339 ---------------------------------------------------------------

# propagation_2014OL339_new_v2

df_2014OL339=read.table("data_asteroids/propagation_2014OL339_new_v2.txt",header = T)

# Transform theta
df_2014OL339=trans_theta(df_2014OL339)

# Select relevant variables
df_2014OL339=df_2014OL339[,c("t","a","theta","type")]
df_2014OL339$type[which(df_2014OL339$type==100)]=10

# Compute features
features_2014OL339=compute_feat(df_2014OL339,l_short=5,l_long=100,
                          tt_thres_maxmin=2.5)


# Remove NAs
features_2014OL339=features_2014OL339[complete.cases(features_2014OL339),]

# Extract ground truth and time
gt_2014OL339=features_2014OL339$type
time=features_2014OL339$t

# Remove ground truth 
features_2014OL339=subset(features_2014OL339,select=-c(type,t))

head(features_2014OL339)

plot(features_2014OL339$dtheta,col=gt_2014OL339+2)
plot(features_2014OL339$sd_dtheta,col=gt_2014OL339+2)
plot(features_2014OL339$mean_osc,col=gt_2014OL339+2)
plot(features_2014OL339$I_QS,col=gt_2014OL339+2)
plot(features_2014OL339$sd_a,col=gt_2014OL339+2)

source("Utils_sparse_robust_2.R")
Y_2014OL339=scale(features_2014OL339)
cv_2014OL339=cv_robust_sparse_jump(
    Y=Y_2014OL339,
    true_states=gt_2014OL339,
    K_grid=2:4,
    zeta0=seq(0.05, 0.5, by=0.05),
    lambda_grid=seq(0,1,.1),
    n_folds = 5,
    parallel=F,
    n_cores=NULL,
    cv_method="blocked-cv",
    knn=10,
    c=10,
    M=NULL
)


# propagation_2015SO2_new_v2
df_2015SO2=read.table("data_asteroids/propagation_2015SO2_new_v2.txt",header = T)
df_2015SO2=trans_theta(df_2015SO2)
unique(df_2015SO2$type)
plot(df_2015SO2$theta,type='p',col=df_2015SO2$type+2)

# propagation_2015XX169_new_v2
df_2015XX169=read.table("data_asteroids/propagation_2015XX169_new_v2.txt",header = T)
df_2015XX169=trans_theta(df_2015XX169)
unique(df_2015XX169$type)
plot(df_2015XX169$theta,type='p',col=df_2015XX169$type+2)

# propagation_2016CA138_new_v2
df_2016CA138=read.table("data_asteroids/propagation_2016CA138_new_v2.txt",header = T)
df_2016CA138=trans_theta(df_2016CA138)
unique(df_2016CA138$type)
df_2016CA138$type[which(df_2016CA138$type==100)]=10
plot(df_2016CA138$theta,type='p',col=df_2016CA138$type+2)

# propagation_2016CO246_new_v2 <--- START WITH THIS
df_2016CO246=read.table("data_asteroids/propagation_2016CO246_new_v2.txt",header = T)
df_2016CO246=trans_theta(df_2016CO246)
unique(df_2016CO246$type)
df_2016CO246$type[which(df_2016CO246$type==100)]=10

plot(df_2016CO246$theta,type='p',col=df_2016CO246$type+2)

# propagation_2016HO3_new_v2
df_2016HO3=read.table("data_asteroids/propagation_2016HO3_new_v2.txt",header = T)
df_2016HO3=trans_theta(df_2016HO3)
plot(df_2016HO3$theta,type='p',col=df_2016HO3$type+2)

# propagation_2019GM1_new_v2
df_2019GM1=read.table("data_asteroids/propagation_2019GM1_new_v2.txt",header = T)
df_2019GM1=trans_theta(df_2019GM1)
unique(df_2019GM1$type)
plot(df_2019GM1$theta,type='p',col=df_2019GM1$type+2)


# propagation_2020PN1_new_v2 <--- THIS IS ALSO GOOD STARTING POINT
df_2020PN1=read.table("data_asteroids/propagation_2020PN1_new_v2.txt",header = T)
df_2020PN1=trans_theta(df_2020PN1)
unique(df_2020PN1$type)
df_2020PN1$type[which(df_2020PN1$type==100)]=10

plot(df_2020PN1$theta,type='p',col=df_2020PN1$type+2)

# propagation_2020PP1_new_v2
df_2020PP1=read.table("data_asteroids/propagation_2020PP1_new_v2.txt",header = T)
df_2020PP1=trans_theta(df_2020PP1)
plot(df_2020PP1$theta,type='p',col=df_2020PP1$type+2)


# propagation_164207_new_v2
df_164207=read.table("data_asteroids/propagation_164207_new_v2.txt",header = T)
df_164207=trans_theta(df_164207)
plot(df_164207$theta,type='p',col=df_164207$type+2)


