# functions --------------------------------------------------------------------

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


# analysis --------------------------------------------------------------------



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



source("Utils_fuzzyJM.R")
#lambda=0.75
lambda=.2
m=1.5

#fit_164207=fuzzy_jump(Y,2,lambda,verbose=TRUE,tol=1e-6)
start=Sys.time()
fit_164207=fuzzy_jump_coord_par(Y,K=2,lambda=lambda,m=m)
end=Sys.time()
end-start

matplot(fit_164207$best_S,type='l')

summary(fit_164207$best_S)

matplot(fit_164207$best_S,col='red')

table(fit_164207$MAP,df164207_select$type)

df164207_select$type<- ifelse(df164207_select$type == -1, 2, 1)

res_compare=data.frame(MAP=fit_164207$MAP,
                       true=df164207_select$type)

res_compare$MAP=as.factor(res_compare$MAP)
res_compare$true=as.factor(res_compare$true)
res_compare$true=recode_factor(res_compare$true,"1"="2","2"="1")

caret::confusionMatrix(res_compare$MAP,res_compare$true)

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
##########
library(ggplot2)
library(patchwork)

sz = 18  # Text size for labels

# Function to normalize values to [0,1]
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Function to plot theta/omega with prob_state_1 on a secondary axis
plot_dual_axis <- function(y_var, y_label, show_x_axis = TRUE, range = 1800:4500) {
  data_subset <- res_164207[range,]
  
  # Normalize both y_var (theta/omega) and prob_state_1
  y_scaled <- normalize(data_subset[[y_var]])
  prob_scaled <- normalize(data_subset$prob_state_1)
  
  ggplot(data_subset, aes(x = t)) +
    # Plot normalized theta/omega (black line, left axis)
    geom_line(aes(y = y_scaled), color = "black", size = 1) +
    # Plot normalized prob_state_1 (red dashed line, right axis)
    geom_line(aes(y = prob_scaled), color = "red", size = 1, alpha = 0.5, linetype = "longdash") +
    # Scale axes back to original range
    scale_y_continuous(
      name = y_label,  # Left y-axis label for theta/omega
      sec.axis = sec_axis(~ ., name = expression(s[QS]))  # Right y-axis for prob_state_1
    ) +
    scale_x_continuous(labels = scales::scientific_format()) +  # Set x-axis to scientific notation
    labs(x = ifelse(show_x_axis, "Time", ""), y = y_label) +  # Hide x-axis label for top plot
    theme_minimal() +
    theme(
      legend.position = "none",  # No legend needed
      axis.title.y.left = element_text(color = "black", size = sz, margin = margin(r = 10)),  # Avoid overlap
      axis.title.y.right = element_text(color = "red", size = sz, margin = margin(l = 10)), 
      axis.text.y.left = element_text(size = sz-8),
      axis.text.y.right = element_text(size = sz-8),
      axis.text.x = element_text(size = ifelse(show_x_axis, sz-8, 0)),  # Hide x-axis labels for top plot
      axis.title.x = element_text(size = ifelse(show_x_axis, sz-3, 0)),  # Hide x-axis title for top plot
      axis.ticks.x = element_blank()  # Remove x-ticks for the top plot
    )
}

# Create the two plots (hide x-axis on the first one)
plot_theta <- plot_dual_axis("theta", expression(theta), show_x_axis = FALSE)
plot_omega <- plot_dual_axis("omega", expression(omega), show_x_axis = TRUE)

# Arrange plots vertically with a single shared x-axis
final_plot <- plot_theta / plot_omega + plot_layout(heights = c(1, 1.2))  # Adjust height for clarity

# Print the final combined plot
print(final_plot)
