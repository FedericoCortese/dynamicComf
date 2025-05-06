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
               "I_TP",
               "I_HS","I_QS",
               "I_CP",
               "mean_osc")]
}

a_feat=function(data,l=5){
  cust_fun=function(x){
    all(x==T)
  }
  
  ind_a=I(data$a<1)
  ind_a2=I(data$a>1)
  
  ind_a_mov <- runner::runner(ind_a, k = l, f = cust_fun, na_pad = TRUE)
  ind_a2_mov <- runner::runner(ind_a2, k = l, f = cust_fun, na_pad = TRUE)
  #data$ind_a_short=as.factor(as.numeric(ind_a_mov+ind_a2_mov))
  data$ind_a_l=as.numeric(ind_a_mov+ind_a2_mov)
  return(data[,c("t","ind_a_l")])
}

# package loading ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(lubridate)

source("Utils_sparse_robust_2.R")


# 164207 ------------------------------------------------------------------

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
df164207$I_domega=(df164207$domega<=-0.0025)

df164207$sd_theta=zoo::rollapply(df164207$theta, 5, sd, fill=NA)
df164207$sd_dtheta=zoo::rollapply(df164207$dtheta, 5, sd, fill=NA)


df164207_maxmin_theta=max_min_feat(df164207,
                                   tt_thres_maxmin=2.5,l=5)
df164207=merge(df164207_maxmin_theta,
               df164207,by="t")
# df164207_select=df164207[,c("t","value_min","value_max",
#                             "I_HS","I_QS",
#                             "mean_osc",
#                             "I_domega","type")]

df164207_select=df164207[,c("t",
                            "theta",
                            "omega",
                            "dtheta",
                            "domega",
                            "sd_theta",
                            "sd_dtheta",
                            #"value_min","value_max",
                            #"I_HS","I_QS",
                            #"mean_osc",
                            #"I_domega",
                            "type")]

str(df164207_select)
# df164207_select$I_HS=as.factor(df164207_select$I_HS)
# df164207_select$I_QS=as.factor(df164207_select$I_QS)
#df164207_select$I_domega=as.factor(as.numeric(df164207_select$I_domega))

str(df164207_select)
df164207_select=df164207_select[complete.cases(df164207_select),]

Y=subset(df164207_select,select=-c(t,type))

str(Y)



source("Utils_sparse_robust_2.R")
#lambda=0.75
#fit_164207=fuzzy_jump(Y,2,lambda,verbose=TRUE,tol=1e-6)
start=Sys.time()
fit_164207_rob=
  JM_COSA(Y,zeta0=.2,lambda=.25,
          K=2,tol=NULL,n_outer=10,alpha=.1,verbose=T)
end=Sys.time()
end-start

plot(Y$theta,col=fit_164207_rob$s)

# 100011836 ---------------------------------------------------------------


load("C:/Users/federico/OneDrive - CNR/SJM-for_ACOMC/cleaned data for analysis in Nonlinear Dynamics/100011836_cleaned.Rdata")

df100011836$theta=trans_theta(df100011836$theta)

df100011836$da=c(NA,diff(df100011836$a))
df100011836$dtheta=c(NA,diff(df100011836$theta))
df100011836$de=c(NA,diff(df100011836$e))
df100011836$domega=c(NA,diff(df100011836$omega))

df100011836$sd_theta=zoo::rollapply(df100011836$theta, 5, sd, fill=NA)
df100011836$sd_dtheta=zoo::rollapply(df100011836$dtheta, 5, sd, fill=NA)

df100011836_maxmin_theta=max_min_feat(df100011836,
                                   tt_thres_maxmin=2.5,l=5)

df100011836=merge(df100011836,df100011836_maxmin_theta,by='t')

df100011836$maxmins=ifelse(df100011836$maxs==1,2,0)
df100011836$maxmins=df100011836$maxmins+ifelse(df100011836$mins==1,1,0)

df100011836_ind_a=a_feat(df100011836,l=10)

df100011836=merge(df100011836,df100011836_ind_a,by='t')


plotly::plot_ly(
  data = df100011836,
  x = ~t,
  y = ~as.numeric(theta),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~maxmins+1,         # Colore basato su 'truth'
  # colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plot(df100011836$theta[1:100],type='l')

plotly::plot_ly(
  data = df100011836,
  x = ~t,
  y = ~as.numeric(theta),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~I_CP+1,         # Colore basato su 'truth'
  # colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)


plotly::plot_ly(
  data = df100011836,
  x = ~t,
  y = ~as.numeric(a),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~ind_a_l,         # Colore basato su 'truth'
  # colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

# QS feature is not well defined

Y=df100011836[,c("t",
                 "a",
                 "da",
                 "theta",
                 "dtheta",
                 "e",
                 "de",
                 "omega",
                 "domega",
                 "sd_theta",
                 "sd_dtheta",
                 #"I_TP","I_HS",
                 #"I_QS",
                # "I_CP",
                 #"ind_a_l",
                 "mean_osc")]

Y=Y[complete.cases(Y),]

timesY=Y$t
Y=Y[,-1]

str(Y)

# Y$I_CP=factor(Y$I_CP)
# Y$I_TP=factor(Y$I_TP)
# Y$I_HS=factor(Y$I_HS)
# Y$I_QS=factor(Y$I_QS)
# Y$ind_a_l=factor(Y$ind_a_l)

str(Y)
dim(Y)

TT=dim(Y)[1]
zoom=400:(TT/3)
Y=Y[zoom,]

fit_100011836=JM_COSA(Y,zeta0=.2,lambda=.25,
                      K=5,tol=1e-16,n_outer=20,alpha=.1,verbose=T)

W_100011836=fit_100011836$W
colnames(W_100011836)=colnames(Y)

round(W_100011836,3)

s_100011836=fit_100011836$s
medoids_100011836=fit_100011836$medoids

df_res_100011836=data.frame(t=timesY[zoom],
  Y,
                            State=s_100011836
)


tapply(df_res_100011836$theta,df_res_100011836$State,mean)
tapply(df_res_100011836$theta,df_res_100011836$State,sd)

tapply(df_res_100011836$a,df_res_100011836$State,mean)


df_res_100011836 <- df_res_100011836 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a))

#df_segments_a$zoom_group <- ifelse(df_segments_a$t < df_segments_a$t[dim(Y)[1] / 2], "First Half", "Second Half")
label_size=18
p_a_100011836=ggplot(data = df_segments_a[zoom,]) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = a, color = as.factor(State))) +
  scale_color_manual(values = c(4, 2, 3, 1,6),
                     labels = c("NR", "QS", "CP", "HS","TP"),
                     name = "Orbital regime") +
  # scale_x_continuous(breaks=seq(range(df_segments_a$t[zoom])[1],
  #                               range(df_segments_a$t[zoom])[2],
  #                               length.out=4),
  #                    labels = label_scientific())+
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(a ~ "(AU)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
  # +
#   facet_zoom(x = t < df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1)

#ggplotly(p_a_res)

p_a_100011836

df_segments_theta <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

#df_segments_theta$zoom_group <- ifelse(df_segments_theta$t < df_segments_theta$t[dim(Y)[1] / 2], "First Half", "Second Half")

label_size=18
# Plot con facet_zoom
p_theta_100011836=ggplot(data = df_segments_theta[zoom,]) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = theta, color = as.factor(State))) +
  scale_color_manual(values = c(4, 2, 3, 1,6),
                     labels = c("NR", "QS", "CP", "HS","TP"),
                     name = "Orbital regime") +
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(theta ~ "(rad)")) +
  scale_y_continuous(
    breaks = c(-pi, 0, pi),  # Specify where to place the labels
    labels = c(expression(-pi), expression(0), expression(pi))  # Use LaTeX-style labels
  )+
  # scale_x_continuous(breaks=seq(range(df_segments_theta$t[zoom])[1],
  #                               range(df_segments_theta$t[zoom])[2],
  #                               length.out=4),
  #                    labels = label_scientific())+
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  )+
  guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1)

p_theta_100011836

ggarrange(p_a_100011836,p_theta_100011836,nrow=2,common.legend = T)



# 100006174 -------------------------------------------------------------------------


load("C:/Users/federico/OneDrive - CNR/SJM-for_ACOMC/cleaned data for analysis in Nonlinear Dynamics/100006174_cleaned.Rdata")


df100006174$theta=trans_theta(df100006174$theta)

df100006174$da=c(NA,diff(df100006174$a))
df100006174$dtheta=c(NA,diff(df100006174$theta))
df100006174$de=c(NA,diff(df100006174$e))
df100006174$domega=c(NA,diff(df100006174$omega))

#l=5
l=10

df100006174$sd_theta=zoo::rollapply(df100006174$theta, l, sd, fill=NA)
df100006174$sd_dtheta=zoo::rollapply(df100006174$dtheta, l, sd, fill=NA)

df100006174$sd_a=zoo::rollapply(df100006174$a, l, sd, fill=NA)
df100006174$sd_da=zoo::rollapply(df100006174$da, l, sd, fill=NA)

df100006174$sd_e=zoo::rollapply(df100006174$e, l, sd, fill=NA)
df100006174$sd_de=zoo::rollapply(df100006174$de, l, sd, fill=NA)

df100006174$sd_omega=zoo::rollapply(df100006174$omega, l, sd, fill=NA)
df100006174$sd_domega=zoo::rollapply(df100006174$domega, l, sd, fill=NA)

df100006174_maxmin_theta=max_min_feat(df100006174,
                                      tt_thres_maxmin=2.5,l=l)

df100006174=merge(df100006174,df100006174_maxmin_theta,by='t')

df100006174$maxmins=ifelse(df100006174$maxs==1,2,0)
df100006174$maxmins=df100006174$maxmins+ifelse(df100006174$mins==1,1,0)

df100006174_ind_a=a_feat(df100006174,l=l)

df100006174=merge(df100006174,df100006174_ind_a,by='t')


plotly::plot_ly(
  data = df100006174,
  x = ~t,
  y = ~as.numeric(theta),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~maxmins+1,         # Colore basato su 'truth'
  # colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

plot(df100006174$theta[1:100],type='l')

plotly::plot_ly(
  data = df100006174,
  x = ~t,
  y = ~as.numeric(theta),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~I_CP+1,         # Colore basato su 'truth'
  # colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)


plotly::plot_ly(
  data = df100006174,
  x = ~t,
  y = ~as.numeric(a),
  type = 'scatter',
  mode = 'markers',
  marker = list(size = 8),
  color = ~ind_a_l,         # Colore basato su 'truth'
  # colors = c("red", "blue"),  # Palette colori personalizzabile
  #text = ~paste("Time:", time, "<br>Theta:", theta, "<br>State:", MAP),
  hoverinfo = "text"
)

# QS feature is not well defined

Y=df100006174[,c("t",
                 "a",
                 "da",
                 "theta",
                 "dtheta",
                 "e",
                 "de",
                 "omega",
                 "domega",
                 "sd_theta",
                 "sd_dtheta",
                 "sd_a",
                 "sd_da",
                 "sd_e",
                 "sd_de",
                 "sd_omega",
                 "sd_domega"
                 #"I_TP","I_HS",
                 #"I_QS",
                 # "I_CP",
                 #"ind_a_l",
                 #"mean_osc"
                 )]

Y=Y[complete.cases(Y),]

timesY=Y$t
Y=Y[,-1]

str(Y)

# Y$I_CP=factor(Y$I_CP)
# Y$I_TP=factor(Y$I_TP)
# Y$I_HS=factor(Y$I_HS)
# Y$I_QS=factor(Y$I_QS)
# Y$ind_a_l=factor(Y$ind_a_l)

dim(Y)

TT=dim(Y)[1]
zoom=400:(TT/3)
Y=Y[zoom,]

fit_100006174=JM_COSA(Y,zeta0=.2,lambda=.25,
                      K=4,tol=1e-16,n_outer=20,alpha=.1,verbose=T)

W_100006174=fit_100006174$W
colnames(W_100006174)=colnames(Y)

round(W_100006174,3)

s_100006174=fit_100006174$s
medoids_100006174=fit_100006174$medoids

df_res_100006174=data.frame(t=timesY[zoom],
                            Y,
                            State=s_100006174
)

label_size=18

df_res_100006174 <- df_res_100006174 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a))

p_a_100006174=ggplot(data = df_segments_a[zoom,]) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = a, color = as.factor(State))) +
  scale_color_manual(values = c(4, 1, 2, 3),
                     labels = c("NR", "HS", "QS", "CP"),
                     name = "Orbital regime") +
  #scale_x_continuous(labels = label_scientific())+
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(a ~ "(AU)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1)

#ggplotly(p_a_res)
p_a_100006174

df_segments_theta <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

#df_segments_theta$zoom_group <- ifelse(df_segments_theta$t < df_segments_theta$t[dim(Y)[1] / 2], "First Half", "Second Half")

# Plot con facet_zoom
p_theta_100006174=ggplot(data = df_segments_theta[zoom,]) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = theta, color = as.factor(State))) +
  scale_color_manual(values = c(4, 1, 2, 3),
                     labels = c("NR", "HS", "QS", "CP"),
                     name = "Orbital regime") +
  #scale_x_continuous(labels = label_scientific())+
  scale_y_continuous(
    breaks = c(-pi, 0, pi),  # Specify where to place the labels
    labels = c(expression(-pi), expression(0), expression(pi))  # Use LaTeX-style labels
  )+
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(theta ~ "(rad)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  ) + guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1)

X11()
ggarrange(p_a_100006174,p_theta_100006174,nrow=2,common.legend = T)

