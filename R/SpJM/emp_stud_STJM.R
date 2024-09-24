###############################################
### Empirical study of the STJM model paper ###
###############################################

# Rich dataset of weather data
load("C:/Users/federico/OneDrive - CNR/General - Comfort/Kaggle/cleaned data/cleaned_data.RData")

# Dataset pulito da Antonio
load("C:/Users/federico/OneDrive - CNR/General - Comfort/Kaggle/singapore/q_outdoor.RData")

library(dplyr)
library(lubridate) 
library(tidyr)
library(tidyverse)

cozie_q.wd <- cozie_q.wd %>%
  mutate(date = as.Date(time))  # Extract the date (year-month-day) from the time column

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
  select(-date)

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
  select(-date)

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
  select(-date)

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
  select(-date)

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
  select(-date)

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
TT=nrow(data_air_temp)
# Introduce temporal gaps
pGap=.05
set.seed(123)
gaps=sort(sample(1:TT,round(TT*pGap),replace=F))
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

data_weather_long %>%
  ggplot(aes(x = time, y = wind_dir)) +
  geom_line(col='purple',linewidth=.6) +
  facet_wrap(~station, scales = "free_y") +
  labs(title = "Wind direction (°)", x = "Time", y = " ") +
  theme_minimal()

# Retrieve relevant locations
loc_weath=locations %>%
  filter(id %in% stats)


# Save the data
save(cozie_compare,loc_weath,
     data_air_temp, data_rh, data_rainfall, data_wind_speed, data_wind_dir, 
     data_weather_long, 
     file = "data_weather.RData")

# Load the data
load("data_weather.RData")

# Summary of the data
summary(data_weather_long)

# Missinf values
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


library(geosphere)
# Distances (in km) between weather stations and feedback locations
cozie_dist=distm(cozie_compare[,c("ws_longitude","ws_latitude")],loc_weath[,c("longitude","latitude")], 
      fun = distGeo)/1000
summary(cozie_dist)

# Radius for the neighbourhood
R=max(apply(cozie_dist, 1, min ))

# Which stations are in the neighbourhood of each feedback location
neigh=cozie_dist<R
colnames(neigh)=loc_weath$id
# Distances between stations
D=distm(loc_weath[,c("longitude","latitude")], loc_weath[,c("longitude","latitude")], 
      fun = distGeo)/1000

# Round times in cozie_compare to hours
wdn="1 hour"
cozie_compare$time=round_date(cozie_compare$time,wdn)

tmp=list()
for(i in loc_weath$id){
  tmp[[i]]=cozie_compare$time[which(neigh[,i])]
  cozie_compare$time[i]
}
