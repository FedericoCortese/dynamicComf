###############################################
### Empirical study of the STJM model paper ###
###############################################


# Packages loading --------------------------------------------------------
library(dplyr)
library(lubridate) 
library(tidyr)
library(tidyverse)
library(shiny)
library(geosphere)
library(caret)

source("Utils.R")

# Raw data loading --------------------------------------------------------
# Rich dataset of weather data
load("C:/Users/federico/OneDrive - CNR/General - Comfort/Kaggle/cleaned data/cleaned_data.RData")

# Dataset pulito da Antonio
load("C:/Users/federico/OneDrive - CNR/General - Comfort/Kaggle/singapore/q_outdoor.RData")

cozie_q.wd <- cozie_q.wd %>%
  mutate(date = as.Date(time))  # Extract the date (year-month-day) from the time column
wdn="1 hour"
cozie_q.wd$time=round_date(cozie_q.wd$time,wdn)
cozie_q.wd$ws_timestamp_start=round_date(cozie_q.wd$time,wdn)
cozie_q.wd$ws_timestamp_location=round_date(cozie_q.wd$ws_timestamp_location,wdn)

# Retrieve data for future comparison with STJM estimates
time_range=c(as.Date('2023-04-18'),as.Date('2023-04-25'))
cozie_compare=cozie_q.wd %>%
  group_by(date) %>%
  filter(date >= time_range[1] & date <= time_range[2])


# # Count the number of observations for each specific day
# daily_counts <- cozie_q.wd %>%
#   group_by(date) %>%
#   summarise(observation_count = n())
# 
# # Print the result
# print(daily_counts)
# plot(daily_counts$date, daily_counts$observation_count, type = "l", xlab = "Date", ylab = "Number of observations", main = "Number of observations per day")
# 
# max_date <- daily_counts$date[
# which.max(daily_counts$observation_count)]  # Find the date with the most observations
# 
# # Select only the data for the date with the most observations
# max_date_data <- cozie_q.wd %>%
#   filter(date == max_date)

stats=intersect(colnames(air_temp),colnames(rainfall))
stats=intersect(stats,colnames(RH))
stats=intersect(stats,colnames(wind_speed))
stats=intersect(stats,colnames(wind_dir))

loc_weath=locations %>%
  filter(id %in% stats)
plot(y=loc_weath$latitude,x=loc_weath$longitude, xlab = "Long", 
     ylab = "Lat", main = " ", pch = 19, col = "blue")
points(
  y=cozie_compare$ws_latitude, x=cozie_compare$ws_longitude,
  col='red',pch=19)


# Weather data ------------------------------------------------------------

# For air_temp, retrieve the data for the same time range as the comfort data
data_air_temp <- air_temp%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
apply(data_air_temp, 2, function(x) 100*sum(is.na(x))/length(x))

# Drop S102
data_air_temp=subset(data_air_temp, select = -c(S102))

# Summarize by hourly mean
wdn="1 hour"
data_air_temp=data_air_temp%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE) 

data_air_temp=data.frame(data_air_temp)
max(diff(data_air_temp$time))

# For RH, retrieve the data for the same time range as the comfort data

data_rh <- RH%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
apply(data_rh, 2, function(x) 100*sum(is.na(x))/length(x))

# Drop S102
data_rh=subset(data_rh, select = -c(S102))

# Summarize by hourly mean
wdn="1 hour"
data_rh=data_rh%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_rh$time))

# For rainfall, retrieve the data for the same time range as the comfort data

data_rainfall <- rainfall%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
sts=apply(data_rainfall, 2, function(x) 100*sum(is.na(x))/length(x))

# Drop stations with 100% NAs
drp=which(sts==100)
data_rainfall=subset(data_rainfall, select = -drp)

# Summarize by hourly mean
wdn="1 hour"
data_rainfall=data_rainfall%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_rainfall$time))

# For wind_speed, retrieve the data for the same time range as the comfort data

data_wind_speed <- wind_speed%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
sts=apply(data_wind_speed, 2, function(x) 100*sum(is.na(x))/length(x))
sts

# Drop stations with 100% NAs
drp=which(sts==100)
data_wind_speed=subset(data_wind_speed, select = -drp)

# Summarize by hourly mean
wdn="1 hour"
data_wind_speed=data_wind_speed%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_wind_speed$time))

# For wind_dir, retrieve the data for the same time range as the comfort data

data_wind_dir <- wind_dir%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
sts=apply(data_wind_dir, 2, function(x) 100*sum(is.na(x))/length(x))
sts

# Drop stations with 100% NAs
drp=which(sts==100)
data_wind_dir=subset(data_wind_dir, select = -drp)

# Summarize by hourly mean
wdn="1 hour"
data_wind_dir=data_wind_dir%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_wind_dir$time))

TT=nrow(data_air_temp)
data_air_temp$t=1:TT


# Introduce temporal gaps
pGap=.05
time_set=data_air_temp$time
time_set=time_set[-which(time_set%in%unique(cozie_compare$time))]
set.seed(1234)
gaps=sort(sample(time_set,round(TT*pGap),replace=F))
gaps=which(data_air_temp$time%in%gaps)

times_complete=data_air_temp$time

data_air_temp=data.frame(data_air_temp)
data_air_temp$t=1:TT
data_air_temp=data_air_temp[-gaps,]

data_rh=data.frame(data_rh)
data_rh$t=1:TT
data_rh=data_rh[-gaps,]

data_rainfall=data.frame(data_rainfall)
data_rainfall$t=1:TT
data_rainfall=data_rainfall[-gaps,]

data_wind_speed=data.frame(data_wind_speed)
data_wind_speed$t=1:TT
data_wind_speed=data_wind_speed[-gaps,]

data_wind_dir=data.frame(data_wind_dir)
data_wind_dir$t=1:TT
data_wind_dir=data_wind_dir[-gaps,]


# Keep only relevant stations
stats=intersect(colnames(data_air_temp),colnames(data_rainfall))
stats=intersect(stats,colnames(data_rh))
stats=intersect(stats,colnames(data_wind_speed))
stats=intersect(stats,colnames(data_wind_dir))

data_air_temp=data_air_temp[,stats]
data_rainfall=data_rainfall[,stats]
data_rh=data_rh[,stats]
data_wind_speed=data_wind_speed[,stats]
data_wind_dir=data_wind_dir[,stats]

# Long format
data_air_temp_long <- data_air_temp %>%
  gather(key = "station", value = "air_temp", -c(t,time))
data_rh_long <- data_rh %>%
  gather(key = "station", value = "rh", -c(t,time))
data_rainfall_long <- data_rainfall %>%
  gather(key = "station", value = "rainfall", -c(t,time))
data_wind_speed_long <- data_wind_speed %>%
  gather(key = "station", value = "wind_speed", -c(t,time))
data_wind_dir_long <- data_wind_dir %>%
  gather(key = "station", value = "wind_dir", -c(t,time))

# Merge all the weather data
data_weather_long <- data_air_temp_long %>%
  left_join(data_rh_long, by = c("t","time", "station")) %>%
  left_join(data_rainfall_long, by = c("t","time", "station")) %>%
  left_join(data_wind_speed_long, by = c("t","time", "station")) %>%
  left_join(data_wind_dir_long, by = c("t","time", "station"))

# Retrieve relevant locations
loc_weath=locations %>%
  filter(id %in% stats)


# Save the data
save(cozie_compare,loc_weath,
     data_air_temp, data_rh, data_rainfall, data_wind_speed, data_wind_dir, 
     data_weather_long, 
     file = "data_weather.RData")


# Load data  --------------------------------------------------------------

# Load weather data
load("data_weather.RData")


# ggplot each variable, for each station, in a grid
data_weather_long %>%
  ggplot(aes(x = time, y = air_temp)) +
  geom_line(col='red',linewidth=.6) +
  facet_wrap(~station, scales = "free_y") +
  labs(title = "Air temperature (°C)", x = "Time", y = " ") +
  theme_minimal()

data_weather_long %>%
  ggplot(aes(x = time, y = rh)) +
  geom_line(col='orange',linewidth=.6) +
  facet_wrap(~station, scales = "free_y") +
  labs(title = "Relative humidity (%)", x = "Time", y = " ") +
  theme_minimal()

data_weather_long %>%
  ggplot(aes(x = time, y = rainfall)) +
  geom_line(col='green',linewidth=.6) +
  facet_wrap(~station, scales = "free_y") +
  labs(title = "Rainfall (mm)", x = "Time", y = " ") +
  theme_minimal()

data_weather_long %>%
  ggplot(aes(x = time, y = wind_speed)) +
  geom_line(col='blue',linewidth=.6) +
  facet_wrap(~station, scales = "free_y") +
  labs(title = "Wind speed (m/s)", x = "Time", y = " ") +
  theme_minimal()

# data_weather_long %>%
#   ggplot(aes(x = time, y = wind_dir)) +
#   geom_line(col='purple',linewidth=.6) +
#   facet_wrap(~station, scales = "free_y") +
#   labs(title = "Wind direction (°)", x = "Time", y = " ") +
#   theme_minimal()

# Summary of the data
summary(data_weather_long)

# Missing values
apply(data_air_temp, 2, function(x) 100*sum(is.na(x))/length(x))
apply(data_rh, 2, function(x) 100*sum(is.na(x))/length(x))
apply(data_rainfall, 2, function(x) 100*sum(is.na(x))/length(x))
apply(data_wind_speed, 2, function(x) 100*sum(is.na(x))/length(x))
apply(data_wind_dir, 2, function(x) 100*sum(is.na(x))/length(x))

# Spatial distribution of stations with respect to feedback locations
plot(y=loc_weath$latitude,x=loc_weath$longitude, xlab = "Long", 
     ylab = "Lat", main = " ", pch = 19, col = "blue")
points(
  y=cozie_compare$ws_latitude, x=cozie_compare$ws_longitude,
  col='red',pch=19)


# Match geo and weather data ---------------------------------------------

# Load geo data
load("geo_feat.Rdata")

D_geo_weath=distm(loc_weath[,c("longitude","latitude")], 
        geo_feat[,c("longitude","latitude")], 
        fun = distGeo)/1000
indx=apply(D_geo_weath, 1,which.min)

# 21 is teh closest with no 0 entries for the indexes
indx8=sort(D_geo_weath[8,])[21]
i8=which(D_geo_weath[8,]==indx8)
geo_feat[i8,]

indx[8]=i8

# Maximum distance between locations of SVF and GVF and weather station is less than 800 m
max(diag(D_geo_weath[,indx]))

data_geo_weath=geo_feat[indx,]
data_geo_weath$station=loc_weath$id

# Does the order of the stations match? YES!
plot(data_geo_weath$longitude,data_geo_weath$latitude, xlab = "Long", 
     ylab = "Lat", main = " ", pch = 19, col = as.numeric(data_geo_weath$station))
points(loc_weath$longitude,loc_weath$latitude, col = as.numeric(loc_weath$id), pch = 2)



# Features construction ---------------------------------------------------

M=dim(loc_weath)[1]
data_stat_number=data.frame(station=loc_weath$id, m=1:M,
                            longitude=loc_weath$longitude,
                            latitude=loc_weath$latitude)
data_stat_number=data_stat_number[order(data_stat_number$m),]

# data_stat_number=merge(data_stat_number,data_geo_weath[,c("station",                                                           "green_view_mean",                                                           "sky_view_mean",
#                                                           "building_view_mean")],
#                        by="station")
# 
# Y=data_weather_long %>%
#   left_join(data_stat_number[,c("station","m","green_view_mean", "sky_view_mean")], by = c("station"))

Y=data_weather_long %>%
  left_join(data_stat_number[,c("station","m")], by = c("station"))



# Order based on t and m
Y=Y[order(Y$t,Y$m),]
rownames(Y)=NULL

# Add categorical features
#Y$rainy=ifelse(Y$rainfall>0,1,0)+1
#Y$windy=ifelse(Y$wind_speed>3,1,0)+1

# Construct rainy feature according to Table 1 of Marsico 2021
# Y$rainy=factor(
#   cut(Y$rainfall, breaks = c(-Inf, 2.5, 10, 50, Inf), labels = c(1, 2, 3, 4)),
#   levels = c(1, 2, 3, 4)
# )
# Y$rainy=droplevels(Y$rainy)
# 
# # Same for Beaufort scale, remember to start from light air
Y$windy=factor(
  cut(Y$wind_speed, breaks = c(-Inf, 1.5,
                               3.3,5.4,7.9,10.7,
                               13.8,Inf), labels = 1:7,
  levels = 1:7
))
Y$windy=droplevels(Y$windy)

Y$hour=hour(Y$time)

# Change NaN in NA
for(i in 2:ncol(Y)){
  Y[,i]=ifelse(is.nan(Y[,i]),NA,Y[,i])
}

Y$hour=as.factor(Y$hour)
#Y$rainy=as.factor(Y$rainy)
Y$windy=as.factor(Y$windy)

# Drop time and station columns
# times=Y[,1]
# Y=Y[,-c(1,3)]
# head(Y)
# str(Y)

# Drop wind direction
Y=subset(Y,select=-wind_dir)

summary(Y)
cor(Y[complete.cases(Y),4:7])


# Merge with solar power data ---------------------------------------------------

source("Utils.R")
S100=read.csv("S100.csv",header = T)
S100=pw_time(S100,var_name="S100")
S104=read.csv("S104.csv",header = T)
S104=pw_time(S104,var_name="S104")
S107=read.csv("S107.csv",header = T)
S107=pw_time(S107,var_name="S107")
S108=read.csv("S108.csv",header = T)
S108=pw_time(S108,var_name="S108")
S109=read.csv("S109.csv",header = T)
S109=pw_time(S109,var_name="S109")
S111=read.csv("S111.csv",header = T)
S111=pw_time(S111,var_name="S111")
S115=read.csv("S115.csv",header = T)
S115=pw_time(S115,var_name="S115")
S116=read.csv("S116.csv",header = T)
S116=pw_time(S116,var_name="S116")
S121=read.csv("S121.csv",header = T)
S121=pw_time(S121,var_name="S121")
S115=read.csv("S115.csv",header = T)
S115=pw_time(S115,var_name="S115")
S116=read.csv("S116.csv",header = T)
S116=pw_time(S116,var_name="S116")
S121=read.csv("S121.csv",header = T)
S121=pw_time(S121,var_name="S121")
S24=read.csv("S24.csv",header = T)
S24=pw_time(S24,var_name="S24")
S43=read.csv("S43.csv",header = T)
S43=pw_time(S43,var_name="S43")
S44=read.csv("S44.csv",header = T)
S44=pw_time(S44,var_name="S44")
S50=read.csv("S50.csv",header = T)
S50=pw_time(S50,var_name="S50")
S60=read.csv("S60.csv",header = T)
S60=pw_time(S60,var_name="S60")

# Merge
data_power=rbind(S100,S104,S107,
                 S108,S109,S111,
                 S115,S116,S121,
                 S24,S43,S44,
                 S50,S60)

Y_complete=merge(Y,data_power,by=c("time","station"))
Y_complete=Y_complete[order(Y_complete$time),]

Y=Y_complete[,-(1:2)]
times=Y_complete[,1]


# Summ stats --------------------------------------------------------------

Y <- Y %>% dplyr::select(t, m, everything())
summary(Y)
cor(Y[complete.cases(Y),-c(1:2,7:9)])

#save(Y,data_stat_number,times,cozie_compare,file="Y.Rdata")


# Load Y ------------------------------------------------------------------
#load("Y.Rdata")
# Remove air pressure
Y=subset(Y,select=-c(PS,ALLSKY_SFC_SW_DWN))
cor(Y[complete.cases(Y),-c(1:2,7:8)])
summary(Y)


# STJM fit ----------------------------------------------------------------


source("Utils.R")
D=distm(data_stat_number[,c("longitude","latitude")], 
        data_stat_number[,c("longitude","latitude")], 
        fun = distGeo)/1000

lambda=.05
gamma=.05

fit=STjumpDist(Y,3,D,
           jump_penalty=lambda,
           spatial_penalty=gamma,
           initial_states=NULL,
           max_iter=10, n_init=10, 
           tol=1e-4, 
           verbose=T,timeflag=T)

# State characterization

State=c(t(fit$best_s))
State=order_states_condMean(Y$air_temp,State)

S_est=matrix(State,ncol=M,byrow = T)

Y_res=data.frame(Y,State,times)
table(Y_res$State)

tapply(Y_res$air_temp,Y_res$State,mean,na.rm=T)
tapply(Y_res$rh,Y_res$State,mean,na.rm=T)
tapply(Y_res$rainfall,Y_res$State,mean,na.rm=T)
tapply(Y_res$wind_speed,Y_res$State,mean,na.rm=T)
#tapply(Y_res$wind_dir,Y_res$State,mean,na.rm=T)
#tapply(Y_res$rainy,Y_res$State,Mode)
tapply(Y_res$windy,Y_res$State,Mode)
tapply(Y_res$hour,Y_res$State,Mode)
# tapply(Y_res$green_view_mean,Y_res$State,mean,na.rm=T)
# tapply(Y_res$sky_view_mean,Y_res$State,mean,na.rm=T)
# tapply(Y_res$ALLSKY_SFC_SW_DWN,Y_res$State,mean,na.rm=T)
tapply(Y_res$ALLSKY_SFC_LW_DWN,Y_res$State,mean,na.rm=T)

# Graphic representation of the results
TY=unique(Y$t)
timesY=unique(times)
TT=length(TY)

# Strano che ws_timestamp_start ha SOLO orari dalle 2 di notte alle 11 di mattina, forse meglio prendere time
q_data=cozie_compare[,c("q_thermal_preference","time","ws_longitude","ws_latitude")]
#q_data=cozie_compare[,c("q_thermal_preference","ws_timestamp_start","ws_longitude","ws_latitude")]
colnames(q_data)[2:4]=c("time","longitude","latitude")
# wdn="1 hour"
# q_data$time=round_date(q_data$time,wdn)
q_data$State=q_data$q_thermal_preference
q_data$State=recode(q_data$State, "Warmer" = 1, "No change" = 2, "Cooler" = 3)

# Hourly average of the data
data_hour_av=Y_res%>%
  group_by(times)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Define UI
ui <- fluidPage(
  titlePanel("Dynamic Plot for Varying t"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("t", "Select t:", 
                  min = min(TY), 
                  max = max(TY), 
                  value = TY[1], 
                  step = NULL)  # step is NULL to allow discrete values
    ),
    
    mainPanel(
      plotOutput("dynamicPlot")
    )
  )
)

server <- function(input, output, session) {
  # Update slider to only allow values from TY
  observe({
    updateSliderInput(session, "t", 
                      min = min(TY), 
                      max = max(TY), 
                      value = TY[1], 
                      step = 1)
  })
  
  output$dynamicPlot <- renderPlot({
    # Ensure t is one of the values from TY
    t_value <- input$t
    
    # Define color mapping based on selected t_value
    #par(mfrow=c(1,2))
    # plot(result$spatial_points, 
    #      col = S_true[t_value, ], 
    #      pch = 19, cex = 1.5, 
    #      main = paste("True"))
    plot(data_stat_number[,c("longitude","latitude")], 
         col = S_est[t_value, ], 
         pch = 17, cex = 1.5, 
         #main = bquote("ST-JM " ~ lambda == .(lambda) ~ " and " ~ gamma == .(gamma))
         main = paste("t = ", timesY[t_value])
         )
    if(any(q_data$time==timesY[t_value])){
      points(q_data$longitude[which(q_data$time==timesY[t_value])], 
             q_data$latitude[which(q_data$time==timesY[t_value])], 
             col = q_data$State[which(q_data$time==timesY[t_value])], pch = 19)
    }
    #legend("topright", legend = 1:3, fill = 1:3)
    legend("topright", legend = c("Cool","Neutral","Hot"), fill = c(1,2,3))
  })
}

# Run the app
shinyApp(ui = ui, server = server)

# Comparison with feedback ------------------------------------------------

# Distances (in km) between weather stations and feedback locations
cozie_dist=distm(q_data[,c("longitude","latitude")],
                 data_stat_number[,c("longitude","latitude")], 
      fun = distGeo)/1000

summary(cozie_dist)

# Radius for the neighbourhood
R=max(apply(cozie_dist, 1, min ))

# Which stations are in the neighbourhood of each feedback location
neigh=cozie_dist<=R
colnames(neigh)=loc_weath$id
neigh=data.frame(time=q_data$time,neigh)
colnames(neigh)[2:ncol(neigh)]=data_stat_number$m
sum(apply(neigh[,-1],1,sum))

S_est=data.frame(time=timesY,fit$best_s)
colnames(S_est)[2:ncol(S_est)]=paste("S", data_stat_number$m)
indx=left_join(neigh,S_est,by="time")
indx=indx[,2:15]*indx[,16:29]
indx[indx==0]=NA

modes=factor(apply(indx,1,Mode,na.rm=T))
comparison=factor(q_data$State)
confusionMatrix(modes,comparison)

adj.rand.index(modes,comparison)



# SJM fit with K=2 --------------------------------------------------------
source("Utils.R")
D=distm(data_stat_number[,c("longitude","latitude")], 
        data_stat_number[,c("longitude","latitude")], 
        fun = distGeo)/1000

lambda=.05
gamma=.05
fit2=STjumpDist(Y,2,D,
           jump_penalty=lambda,
           spatial_penalty=gamma,
           initial_states=NULL,
           max_iter=10, n_init=5, 
           tol=NULL, 
           verbose=T,timeflag=T)

State2=c(t(fit2$best_s))
State2=order_states_condMean(Y$air_temp,State2)

Y_res2=data.frame(Y,State2,times)

tapply(Y_res2$air_temp,Y_res2$State,mean,na.rm=T)
tapply(Y_res2$rh,Y_res2$State,mean,na.rm=T)
tapply(Y_res2$rainfall,Y_res2$State,mean,na.rm=T)
tapply(Y_res2$wind_speed,Y_res2$State,mean,na.rm=T)
#tapply(Y_res2$wind_dir,Y_res2$State,mean,na.rm=T)
tapply(Y_res2$rainy,Y_res2$State,Mode)
tapply(Y_res2$windy,Y_res2$State,Mode)
tapply(Y_res2$hour,Y_res2$State,Mode)
# tapply(Y_res2$green_view_mean,Y_res2$State,mean,na.rm=T)
# tapply(Y_res2$sky_view_mean,Y_res2$State,mean,na.rm=T)
tapply(Y_res2$ALLSKY_SFC_SW_DWN,Y_res2$State,mean,na.rm=T)
tapply(Y_res2$ALLSKY_SFC_LW_DWN,Y_res2$State,mean,na.rm=T)

table(Y_res2$State)

# Consider only two response categories, warmer-no change or cooler
q_data$State2=recode(q_data$q_thermal_preference, "Warmer" = 1, "No change" = 1, "Cooler" = 2)

S_est2=data.frame(time=timesY,fit2$best_s)
colnames(S_est2)[2:ncol(S_est2)]=paste("S", data_stat_number$m)
indx2=left_join(neigh,S_est2,by="time")
indx2=indx2[,2:15]*indx2[,16:29]
indx2[indx2==0]=NA

modes=factor(apply(indx2,1,Mode,na.rm=T))
comparison=factor(q_data$State2)
confusionMatrix(modes,comparison)

# Accuracy around 40%, not good

adj.rand.index(modes,comparison)