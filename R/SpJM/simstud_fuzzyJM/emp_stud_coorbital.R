trans_theta=function(theta){
  #theta_trans <- data$theta
  theta[which(theta > pi)] <- theta[which(theta > pi)] - 2 * pi
  #data$theta <- theta_trans
  return(theta)
}

max_min_feat=function(data,tt_thres_maxmin=2.5,l=5){
  
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

# max_min_feat=function(data,tt_thres_maxmin=3,
#                       #tt_thres_diffmaxmin=pi/4,
#                       l=5
# ){
#   
#   last_max_value <- function(is_max, values) {
#     # Initialize a vector to store the result
#     result <- numeric(length(values))
#     
#     # Initialize the last observed max
#     last_max <- NA
#     
#     # Iterate through the values
#     for (i in seq_along(values)) {
#       if (is_max[i] == 1) {
#         # If the current position is a max, store the value
#         last_max <- values[i]
#       }
#       # Store the last observed max
#       result[i] <- last_max
#     }
#     
#     return(result)
#   }
#   
#   # Define the function for minima
#   last_min_value <- function(is_min, values) {
#     # Initialize a vector to store the result
#     result <- numeric(length(values))
#     
#     # Initialize the last observed min
#     last_min <- NA
#     
#     # Iterate through the values
#     for (i in seq_along(values)) {
#       if (is_min[i] == 1) {
#         # If the current position is a min, store the value
#         last_min <- values[i]
#       }
#       # Store the last observed min
#       result[i] <- last_min
#     }
#     
#     return(result)
#   }
#   
#   # Theta transformation
#   # theta_trans <- data$theta
#   # theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
#   # data$theta <- theta_trans
#   
#   maxs <- as.numeric(splus2R::peaks(data$theta, span = l))
#   mins <- as.numeric(splus2R::peaks(-data$theta, span = l))
#   
#   # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
#   maxs[which(data$theta>(tt_thres_maxmin))] <- 0
#   mins[which(data$theta<(-tt_thres_maxmin))] <- 0
#   
#   #
#   data$maxs=maxs
#   data$mins=mins
#   
#   data_rev=data.frame(apply(data,2,rev))
#   
#   # Add value_min, value_max and diffmaxmin
#   data$value_min <- last_min_value(mins, data$theta)
#   data$value_max <- last_max_value(maxs, data$theta)
#   
#   
#   # TP indicator
#   data$I_TP=as.numeric(data$value_max*data$value_min>0 & data$value_max>data$value_min)
#   
#   # HS indicator
#   data$I_HS=as.numeric(data$value_max*data$value_min<0&data$value_max<data$value_min)
#   
#   # QS indicator
#   data$I_QS=as.numeric(data$value_max*data$value_min<0 & data$value_max>data$value_min)
#   
#   # CP indicator
#   data$I_CP=as.numeric(data$value_max*data$value_min>0&data$value_max<=data$value_min)
#   
#   # Check if distance between max and min is small to distinguish between NR and QS
#   # data$I_QS=data$I_QS*data$I_diffmaxmin
#   
#   
#   
#   data$mean_osc=data$I_HS*(data$value_min+data$value_max+2*pi)/2+
#     data$I_QS*(data$value_min+data$value_max)/2+
#     data$I_TP*(data$value_min+data$value_max)/2
#   
#   data=data[,c("t","maxs","mins",
#                "value_min","value_max",
#                #"I_TP",
#                "I_HS","I_QS",
#                #"I_CP",
#                "mean_osc")]
#   
#   return(data)
# }

# 164207 ------------------------------------------------------------------

# df164207=read.table("propagation_164207_new_v2.txt")
# str(df164207)
# var_names=df164207[1,]
# df164207=df164207[-1,]
# colnames(df164207)=var_names
# df164207=apply(df164207,2,as.numeric)
# df164207=as.data.frame(df164207)
# save(df164207,file="propagation_164207_new_v2.Rdata")

load("propagation_164207_new_v2.Rdata")

df164207$theta=trans_theta(df164207$theta)

df164207$dtheta=c(NA,diff(df164207$theta))
df164207$de=c(NA,diff(df164207$e))
df164207$domega=c(NA,diff(df164207$omega))

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~omega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~domega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

df164207$I_domega=(df164207$domega<=-0.0025)

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~as.numeric(I_domega),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

df164207_maxmin_theta=max_min_feat(df164207,
                                   tt_thres_maxmin=2.5,l=5)

# df164207_rev=data.frame(apply(df164207,2,rev))
# str(df164207_rev)
# 
# df164207_maxmin_theta_rev=max_min_feat(df164207_rev,
#                                    tt_thres_maxmin=2.5,l=5)
# 
# df164207_maxmin_theta_2=merge(df164207_maxmin_theta,
#                               df164207_maxmin_theta_rev,
#                               by="t")
# 
# df164207_maxmin_theta_2$I_HS=
#   df164207_maxmin_theta_2$I_HS.x*df164207_maxmin_theta_2$I_HS.y
# 
# df164207_maxmin_theta_2$I_QS=
#   df164207_maxmin_theta_2$I_QS.x*df164207_maxmin_theta_2$I_QS.y
# 
# 
# plot(df164207_maxmin_theta_2$I_HS,pch=20,col='pink')
# points(df164207_maxmin_theta_2$I_HS.x,pch=20,col='red')
# points(df164207_maxmin_theta_2$I_HS.y,pch=20,col='orange')
# points(df164207_maxmin_theta_2$I_HS,pch=20,col='cyan')
# 
# plot(df164207$theta,pch=20,col=df164207_maxmin_theta_2$I_QS+1)
# lines(df164207_maxmin_theta_2$I_QS,col=2)
# 
# plot(df164207$theta,pch=20,col=df164207_maxmin_theta_2$I_HS+1)
# lines(df164207_maxmin_theta_2$I_HS)
# 
# plot(df164207_maxmin_theta_2$I_QS,type='l')
# lines(df164207_maxmin_theta_2$I_HS,col=2)
# 
# 
# df164207_maxmin_theta_2$I_QS=
#   df164207_maxmin_theta_2$I_QS.x*df164207_maxmin_theta_2$I_QS.y
# 
# true_val_mins=rep(NA,dim(df164207_maxmin_theta_2)[1])
# indx_val_min=which(df164207_maxmin_theta_2$value_min.x==
#   df164207_maxmin_theta_2$value_min.y)
# true_val_mins[indx_val_min]=df164207_maxmin_theta_2$value_min.x[indx_val_min]
# 
# plot(df164207_maxmin_theta_2$value_min.x,pch=20,col='pink')
# points(true_val_mins,pch=20,col='red')
# 
# df164207_maxmin_theta_2$value_max=last_max_value(df164207_maxmin_theta_2$maxs,
#                                                  df164207$theta)
# 
# plot(df164207_maxmin_theta_2$maxs.x,pch=20,col='pink')
# points(df164207_maxmin_theta_2$maxs,pch=20,col='red')
# 
# df164207=merge(df164207_maxmin_theta_2[,c("t","I_HS","I_QS","mins","maxs","value_min","value_max")],
#                df164207,by="t")

df164207=merge(df164207_maxmin_theta,
               df164207,by="t")

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~mean_osc,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~I_QS,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df164207,
  x = ~t,
  y = ~value_max,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

# df164207$msd_theta=zoo::rollapply(df164207$theta, width = 30, 
#                     FUN = sd, 
#                     by = 1, 
#                     align = "right", 
#                     fill = NA)
# 
# df164207$msd_domega=zoo::rollapply(df164207$domega, width = 30, 
#                     FUN = sd, 
#                     by = 1, 
#                     align = "right", 
#                     fill = NA)
# 
# plotly::plot_ly(
#   data = df164207,
#   x = ~t,
#   y = ~msd_theta,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )
# 
# plotly::plot_ly(
#   data = df164207,
#   x = ~t,
#   y = ~msd_domega,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )

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

plotly::plot_ly(
  data = df164207_select,
  x = ~t,
  y = ~as.numeric(I_QS),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

source("Utils_fuzzyJM.R")
#lambda=0.75
lambda=.2

fit_164207=fuzzy_jump(Y,2,lambda,verbose=TRUE,tol=1e-6)

summary(fit_164207$best_S)

table(fit_164207$MAP,df164207_select$type)

mclust::adjustedRandIndex(fit_164207$MAP,df164207_select$type)

res_164207=data.frame(
  t=df164207_select$t,
  # theta=tail(df164207$theta,dim(Y)[1]),
  omega=tail(df164207$omega,dim(Y)[1]),
  Y,
  MAP=fit_164207$MAP,
  prob_state_1=fit_164207$best_S[,1],
  prob_state_2=fit_164207$best_S[,2]
)

tapply(res_164207$theta,res_164207$MAP,mean)
tapply(res_164207$theta,res_164207$MAP,sd)

tapply(res_164207$omega,res_164207$MAP,mean)
tapply(res_164207$omega,res_164207$MAP,sd)


plotly::plot_ly(
  data = res_164207,
  x = ~t,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~prob_state_1,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta, 
                "<br>State 1:", prob_state_1,
                "<br>State 2:", prob_state_2),
  hoverinfo = "text"
) 

plotly::plot_ly(
  data = res_164207,
  x = ~t,
  y = ~omega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~prob_state_1,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta, 
                "<br>State 1:", prob_state_1,
                "<br>State 2:", prob_state_2),
  hoverinfo = "text"
) 

# Static plot
# Load necessary libraries
library(ggplot2)
library(patchwork)  # For arranging plots
sz=15

# Define the base plot function with custom text sizes
plot_scatter <- function(y_var,x_label, y_label,range=1000:4000) {
  ggplot(res_164207[range,], aes(x = t, y = !!sym(y_var), color = prob_state_1)) +
    #geom_line(aes(y=!!sym(y_var)),color="grey80")+
    geom_point(size = 1.5) +
    scale_color_gradient(low = "lightgreen", high = "lightcoral",limits = c(0, 1),
                         labels = scales::number_format(accuracy = 0.1),
                         guide = guide_colorbar(title.position = "top", 
                                                title.hjust = 0.5)) +  
    labs(x=x_label, y = y_label, color = expression(s[QS])) +
    theme_minimal() +
    theme(
      #legend.title.align = 0.5,
      legend.position = "top",                 # Legend on top
      legend.text = element_text(size = sz-3),   # Legend text size
      legend.title = element_text(size = sz-5),  # Legend title size
      axis.title = element_text(size = sz-2),    # Axis labels size
      axis.text = element_text(size = sz-3)      # Axis tick labels size
    )
}

# Create the two plots
plot_theta <- plot_scatter("theta"," ", expression(theta))
plot_omega <- plot_scatter("omega","Time", expression(omega))

# Arrange plots vertically with a shared legend at the top
final_plot <- (plot_theta / plot_omega) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Print the final combined plot
print(final_plot)


# 2002AA29 ----------------------------------------------------------------


# df2002AA29=read.table("propagation_2002AA29_new_v2.txt")
# var_names=df2002AA29[1,]
# df2002AA29=df2002AA29[-1,]
# colnames(df2002AA29)=var_names
# df2002AA29=apply(df2002AA29,2,as.numeric)
# df2002AA29=as.data.frame(df2002AA29)
# save(df2002AA29,file="propagation_2002AA29_new_v2.Rdata")

load("propagation_2002AA29_new_v2.Rdata")

df2002AA29$msd_theta=zoo::rollapply(df2002AA29$theta, width = 10,
                    FUN = sd,
                    by = 1,
                    align = "right",
                    fill = NA)
# 
# df2002AA29$msd_domega=zoo::rollapply(df2002AA29$domega, width = 10,
#                     FUN = sd,
#                     by = 1,
#                     align = "right",
#                     fill = NA)

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~msd_theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)
# 
# plotly::plot_ly(
#   data = df2002AA29,
#   x = ~t,
#   y = ~msd_domega,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )

df2002AA29$theta=trans_theta(df2002AA29$theta)

df2002AA29$dtheta=c(NA,diff(df2002AA29$theta))
df2002AA29$de=c(NA,diff(df2002AA29$e))
df2002AA29$domega=c(NA,diff(df2002AA29$omega))

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~omega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~domega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

df2002AA29$I_domega=(df2002AA29$domega<=-0.05)

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~as.numeric(I_domega),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

df2002AA29_maxmin_theta=max_min_feat(df2002AA29,
                                     tt_thres_maxmin=2.5,l=3)

df2002AA29=merge(df2002AA29_maxmin_theta,df2002AA29,by="t")

# plotly::plot_ly(
#   data = df2002AA29,
#   x = ~t,
#   y = ~theta,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~maxs*2+mins,         # Colore basato su 'truth'
#   #colors = c("red", "blue"),  # Palette colori personalizzabile
#   text = ~paste("Time:", t, "<br>Theta:", theta),
#   hoverinfo = "text"
# )

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~mean_osc,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~I_QS,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2002AA29,
  x = ~t,
  y = ~value_max,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>value_max:", value_max),
  hoverinfo = "text"
)

# df2002AA29$msd_theta=zoo::rollapply(df2002AA29$theta, width = 10,
#                     FUN = sd,
#                     by = 1,
#                     align = "right",
#                     fill = NA)
# 
# df2002AA29$msd_domega=zoo::rollapply(df2002AA29$domega, width = 10,
#                     FUN = sd,
#                     by = 1,
#                     align = "right",
#                     fill = NA)

# plotly::plot_ly(
#   data = df2002AA29,
#   x = ~t,
#   y = ~msd_theta,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )
# 
# plotly::plot_ly(
#   data = df2002AA29,
#   x = ~t,
#   y = ~msd_domega,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )


df2002AA29_select=df2002AA29[,c("t",
                                "theta",
                                #"domega",
                                "value_min","value_max",
                                "I_HS","I_QS","mean_osc",
                                "msd_theta",
                                #"I_domega",
                                "type")]

str(df2002AA29_select)
df2002AA29_select$I_HS=as.factor(df2002AA29_select$I_HS)
df2002AA29_select$I_QS=as.factor(df2002AA29_select$I_QS)
#df2002AA29_select$I_domega=as.factor(df2002AA29_select$I_domega)

str(df2002AA29_select)
df2002AA29_select=df2002AA29_select[complete.cases(df2002AA29_select),]

Y=subset(df2002AA29_select,select=-c(t,type))

str(Y)

source("Utils_fuzzyJM.R")

lambda=0

fit_2002AA29=fuzzy_jump(Y,2,lambda,verbose=TRUE,tol=1e-6)

summary(fit_2002AA29$best_S)

table(fit_2002AA29$MAP,df2002AA29_select$type)

res_2002AA29=data.frame(
  t=df2002AA29_select$t,
  # theta=tail(df2002AA29$theta,dim(Y)[1]),
  omega=tail(df2002AA29$omega,dim(Y)[1]),
  Y,
  MAP=fit_2002AA29$MAP,
  prob_state_1=fit_2002AA29$best_S[,1],
  prob_state_2=fit_2002AA29$best_S[,2]
)

plotly::plot_ly(
  data = res_2002AA29,
  x = ~t,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~prob_state_1,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta, 
                "<br>State 1:", prob_state_1,
                "<br>State 2:", prob_state_2),
  hoverinfo = "text"
)



# 2016HO3 --------------------------------------------------------------------

# df2016HO3=read.table("propagation_2016HO3_new_v2.txt")
# var_names=df2016HO3[1,]
# df2016HO3=df2016HO3[-1,]
# colnames(df2016HO3)=var_names
# df2016HO3=apply(df2016HO3,2,as.numeric)
# df2016HO3=as.data.frame(df2016HO3)
# save(df2016HO3,file="propagation_2016HO3_new_v2.Rdata")
# 
load("propagation_2016HO3_new_v2.Rdata")

df2016HO3$theta=trans_theta(df2016HO3$theta)

df2016HO3$dtheta=c(NA,diff(df2016HO3$theta))
df2016HO3$de=c(NA,diff(df2016HO3$e))
df2016HO3$domega=c(NA,diff(df2016HO3$omega))

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~omega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~domega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

df2016HO3$I_domega=(df2016HO3$domega<=-0.01)

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~as.numeric(I_domega),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

df2016HO3_maxmin_theta=max_min_feat(df2016HO3,
                                    tt_thres_maxmin=2.5,l=3)

df2016HO3=merge(df2016HO3_maxmin_theta,df2016HO3,by="t")

# plotly::plot_ly(
#   data = df2016HO3,
#   x = ~t,
#   y = ~theta,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~maxs*2+mins,         # Colore basato su 'truth'
#   #colors = c("red", "blue"),  # Palette colori personalizzabile
#   text = ~paste("Time:", t, "<br>Theta:", theta),
#   hoverinfo = "text"
# )

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~mean_osc,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~I_QS,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plotly::plot_ly(
  data = df2016HO3,
  x = ~t,
  y = ~value_max,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~type,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>value_max:", value_max),
  hoverinfo = "text"
)

# df2016HO3$msd_theta=zoo::rollapply(df2016HO3$theta, width = 10,
#                     FUN = sd,
#                     by = 1,
#                     align = "right",
#                     fill = NA)
# 
# df2016HO3$msd_domega=zoo::rollapply(df2016HO3$domega, width = 10,
#                     FUN = sd,
#                     by = 1,
#                     align = "right",
#                     fill = NA)

# plotly::plot_ly(
#   data = df2016HO3,
#   x = ~t,
#   y = ~msd_theta,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )
# 
# plotly::plot_ly(
#   data = df2016HO3,
#   x = ~t,
#   y = ~msd_domega,
#   type = 'scatter',
#   mode = 'markers',
#   marker = list(size = 8),
#   color = ~type,         # Colore basato su 'truth'
#   colors = c("red", "blue"),  # Palette colori personalizzabile
#   #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
#   hoverinfo = "text"
# )


df2016HO3_select=df2016HO3[,c("t",
                              "theta",
                              "domega",
                              "value_min","value_max",
                              "I_HS","I_QS","mean_osc",
                              "I_domega",
                              "type")]

str(df2016HO3_select)
df2016HO3_select$I_HS=as.factor(df2016HO3_select$I_HS)
df2016HO3_select$I_QS=as.factor(df2016HO3_select$I_QS)
df2016HO3_select$I_domega=as.factor(df2016HO3_select$I_domega)

str(df2016HO3_select)
df2016HO3_select=df2016HO3_select[complete.cases(df2016HO3_select),]

plot(as.numeric(df2016HO3_select$value_max),pch=20,col=df2016HO3_select$type+4)

Y=subset(df2016HO3_select,select=-c(t,type))

str(Y)

source("Utils_fuzzyJM.R")

lambda=5

fit_2016HO3=fuzzy_jump(Y,2,lambda,verbose=TRUE,tol=1e-6)

summary(fit_2016HO3$best_S)

table(fit_2016HO3$MAP,df2016HO3_select$type)

mclust::adjustedRandIndex(fit_2016HO3$MAP,df2016HO3_select$type)

res_2016HO3=data.frame(
  t=df2016HO3_select$t,
  # theta=tail(df2016HO3$theta,dim(Y)[1]),
  omega=tail(df2016HO3$omega,dim(Y)[1]),
  Y,
  MAP=fit_2016HO3$MAP,
  prob_state_1=fit_2016HO3$best_S[,1],
  prob_state_2=fit_2016HO3$best_S[,2]
)

plotly::plot_ly(
  data = res_2016HO3,
  x = ~t,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~prob_state_1,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta, 
                "<br>State 1:", prob_state_1,
                "<br>State 2:", prob_state_2),
  hoverinfo = "text"
) 

plotly::plot_ly(
  data = res_2016HO3,
  x = ~t,
  y = ~omega,
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~prob_state_1,         # Colore basato su 'truth'
  colors = c("red", "blue"),  # Palette colori personalizzabile
  text = ~paste("Time:", t, "<br>Theta:", theta, 
                "<br>State 1:", prob_state_1,
                "<br>State 2:", prob_state_2),
  hoverinfo = "text"
) 

# 100006174 ---------------------------------------------------------------

load("100006174_cleaned.Rdata")

library(plotly)
plot_ly(data = df100006174) %>%
  add_trace(x = ~t, 
            y = ~theta, 
            type = 'scatter', 
            mode = 'lines', 
            line = list(color = 'lightgray', width = 1),
            name = 'Trend') 

source("Utils_fuzzyJM.R")
lambda=.1

str(Y)

# Convert indicators into factors
Y[,5]=as.factor(Y[,5])
Y[,6]=as.factor(Y[,6])
Y[,7]=as.factor(Y[,7])
Y[,8]=as.factor(Y[,8])
Y[,9]=as.factor(Y[,9])

str(Y)

fit_100006174=fuzzy_jump(Y,4,lambda,verbose=TRUE,tol=1e-6)

summary(fit_100006174$best_S)
table(fit_100006174$MAP)

res_1000016174=data.frame(Date=timesY,
                          tail(df100006174,dim(Y)[1]),
                          marg_prob_1=fit_100006174$best_S[,1],
                          marg_prob_2=fit_100006174$best_S[,2],
                          marg_prob_3=fit_100006174$best_S[,3],
                          marg_prob_4=fit_100006174$best_S[,4],
                          MAP=as.factor(fit_100006174$MAP))


res_1000016174$theta=trans_theta(res_1000016174$theta)

library(plotly)

# Definisci una palette di base per ciascuna marg_prob
base_colors <- list(
  marg_prob_1 = col2rgb("black"),
  marg_prob_3 = col2rgb("cyan"),
  marg_prob_4 = col2rgb("yellow"),
  marg_prob_5 = col2rgb("magenta")
)

# Funzione per mescolare colori secondo le marg_prob
blend_colors <- function(row) {
  # Somma ponderata degli RGB in base alle marg_prob
  r <- sum(row[c("marg_prob_1", "marg_prob_2", "marg_prob_3", "marg_prob_4")] * 
             sapply(base_colors, function(c) c[1])) / sum(row[c("marg_prob_1", "marg_prob_2", 
                                                                "marg_prob_3", "marg_prob_4")])
  g <- sum(row[c("marg_prob_1", "marg_prob_2", "marg_prob_3", "marg_prob_4")] * 
             sapply(base_colors, function(c) c[2])) / sum(row[c("marg_prob_1", "marg_prob_2", 
                                                                "marg_prob_3", "marg_prob_4")])
  b <- sum(row[c("marg_prob_1", "marg_prob_2", "marg_prob_3", "marg_prob_4")] * 
             sapply(base_colors, function(c) c[3])) / sum(row[c("marg_prob_1", "marg_prob_2", 
                                                                "marg_prob_3", "marg_prob_4")])
  
  rgb(
    round(r, 0), 
    round(g, 0), 
    round(b, 0), 
    maxColorValue = 255
  )
}


for(i in 1:nrow(res_1000016174)){
  res_1000016174$color[i] <- blend_colors(res_1000016174[i,])
}

# Plot con linea grigia + punti con colori sfumati compositi
fig <- plot_ly(
  data = res_1000016174
) %>%
  # Linea continua grigia
  add_trace(
    x = ~Date,
    y = ~theta,
    type = 'scatter',
    mode = 'lines',
    line = list(color = 'gray', width = 1),
    name = "Linea continua"
  ) %>%
  # Punti sfumati con colore derivato dalle marg_prob
  add_trace(
    x = ~Date,
    y = ~theta,
    type = 'scatter',
    mode = 'markers',
    marker = list(
      size = 7,
      color = ~color
    ),
    text = ~paste(
      "Date:", Date, "<br>",
      "Theta:", round(theta, 4), "<br>",
      "MAP:", MAP, "<br>",
      "Marg Prob 1 :", round(marg_prob_1, 5), "<br>",
      "Marg Prob 2 :", round(marg_prob_2, 5), "<br>",
      "Marg Prob 3 :", round(marg_prob_3, 5), "<br>",
      "Marg Prob 4 :", round(marg_prob_4, 5)
    ),
    hoverinfo = "text",
    name = "Punti (Colori Sfumati)"
  ) %>%
  layout(
    title = "Oscillazione di Theta con Colori Compositi",
    xaxis = list(title = "Tempo (Date)"),
    yaxis = list(title = "Theta"),
    hovermode = "closest"
  )

# Visualizza il grafico
fig




# 100011836 ---------------------------------------------------------------


load("100011836_cleaned.Rdata")

library(plotly)
plot_ly(data = df100011836) %>%
  add_trace(x = ~t, 
            y = ~theta, 
            type = 'scatter', 
            mode = 'lines', 
            line = list(color = 'lightgray', width = 1),
            name = 'Trend') 

source("Utils_fuzzyJM.R")
lambda=.25

str(Y)

# Convert indicators into factors
Y[,5]=as.factor(Y[,5])
Y[,6]=as.factor(Y[,6])
Y[,7]=as.factor(Y[,7])
Y[,8]=as.factor(Y[,8])
Y[,9]=as.factor(Y[,9])

str(Y)

fit_100011836=fuzzy_jump(Y,5,lambda,verbose=TRUE,tol=1e-6)

summary(fit_100011836$best_S)
table(fit_100011836$MAP)

res_100011836=data.frame(Date=timesY,
                          tail(df100011836,dim(Y)[1]),
                          marg_prob_1=fit_100011836$best_S[,1],
                          marg_prob_2=fit_100011836$best_S[,2],
                          marg_prob_3=fit_100011836$best_S[,3],
                          marg_prob_4=fit_100011836$best_S[,4],
                          marg_prob_5=fit_100011836$best_S[,5],
                          MAP=as.factor(fit_100011836$MAP))

res_100011836$theta=trans_theta(res_100011836$theta)


library(plotly)

# Definisci una palette di base per ciascuna marg_prob
base_colors <- list(
  marg_prob_1 = col2rgb("black"),
  marg_prob_2 = col2rgb("green"),
  marg_prob_3 = col2rgb("cyan"),
  marg_prob_4 = col2rgb("yellow"),
  marg_prob_5 = col2rgb("magenta")
)

# Funzione per mescolare colori secondo le marg_prob
blend_colors <- function(row) {
  # Somma ponderata degli RGB in base alle marg_prob
  r <- sum(row[c("marg_prob_1", "marg_prob_2", "marg_prob_3", "marg_prob_4", "marg_prob_5")] * 
             sapply(base_colors, function(c) c[1])) / sum(row[c("marg_prob_1", "marg_prob_2", 
                                                                "marg_prob_3", "marg_prob_4", "marg_prob_5")])
  g <- sum(row[c("marg_prob_1", "marg_prob_2", "marg_prob_3", "marg_prob_4", "marg_prob_5")] * 
             sapply(base_colors, function(c) c[2])) / sum(row[c("marg_prob_1", "marg_prob_2", 
                                                                "marg_prob_3", "marg_prob_4", "marg_prob_5")])
  b <- sum(row[c("marg_prob_1", "marg_prob_2", "marg_prob_3", "marg_prob_4", "marg_prob_5")] * 
             sapply(base_colors, function(c) c[3])) / sum(row[c("marg_prob_1", "marg_prob_2", 
                                                                "marg_prob_3", "marg_prob_4", "marg_prob_5")])
  
  rgb(
    round(r, 0), 
    round(g, 0), 
    round(b, 0), 
    maxColorValue = 255
  )
}


for(i in 1:nrow(res_100011836)){
  res_100011836$color[i] <- blend_colors(res_100011836[i,])
}

# Plot con linea grigia + punti con colori sfumati compositi
fig <- plot_ly(
  data = res_100011836
) %>%
  # Linea continua grigia
  add_trace(
    x = ~Date,
    y = ~theta,
    type = 'scatter',
    mode = 'lines',
    line = list(color = 'gray', width = 1),
    name = "Linea continua"
  ) %>%
  # Punti sfumati con colore derivato dalle marg_prob
  add_trace(
    x = ~Date,
    y = ~theta,
    type = 'scatter',
    mode = 'markers',
    marker = list(
      size = 7,
      color = ~color
    ),
    text = ~paste(
      "Date:", Date, "<br>",
      "Theta:", round(theta, 4), "<br>",
      "MAP:", MAP, "<br>",
      "Marg Prob 1 :", round(marg_prob_1, 5), "<br>",
      "Marg Prob 2 :", round(marg_prob_2, 5), "<br>",
      "Marg Prob 3 :", round(marg_prob_3, 5), "<br>",
      "Marg Prob 4 :", round(marg_prob_4, 5), "<br>",
      "Marg Prob 5 :", round(marg_prob_5, 5)
    ),
    hoverinfo = "text",
    name = "Punti (Colori Sfumati)"
  ) %>%
  layout(
    title = "Oscillazione di Theta con Colori Compositi",
    xaxis = list(title = "Tempo (Date)"),
    yaxis = list(title = "Theta"),
    hovermode = "closest"
  )

# Visualizza il grafico
fig


