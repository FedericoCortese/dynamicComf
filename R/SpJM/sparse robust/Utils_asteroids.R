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