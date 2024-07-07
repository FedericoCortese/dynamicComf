library(ARPALData)
#!!!!! Frequency cannot be hourly as municipal estimates are provided on a daily basis !!!!!


# AQ data -----------------------------------------------------------------

# 100451 ----> MILANO
dat=get_ARPA_Lombardia_AQ_municipal_data(
  ID_station=100451,
  Date_begin = "2023-12-01",
  Date_end = "2024-06-30",
  Frequency = "daily",
  parallel = T
)

dat_AQ=subset(dat,select = -c(IDStation,NameStation))

# dat$NameStation[grep("Milano",dat$NameStation)]
# dat$IDStation[grep("Milano",dat$NameStation)]

windows()
par(mfrow=c(2,3))
plot(y=dat_AQ$NO2_mean,x=dat_AQ$Date,type='l',main = "NO2 mean",xlab = " ",ylab = " ")
plot(y=dat_AQ$NO2_max_day,x=dat_AQ$Date,type='l',main = "NO2 max",xlab = " ",ylab = " ")
plot(y=dat_AQ$Ozone_max_day,x=dat_AQ$Date,type='l',main = "Ozone max day",xlab = " ",ylab = " ")
plot(y=dat_AQ$Ozone_max_8h,x=dat_AQ$Date,type='l',main = "Ozone max 8h",xlab = " ",ylab = " ")
plot(y=dat_AQ$PM10_mean,x=dat_AQ$Date,type='l',main = "PM10 mean",xlab = " ",ylab = " ")
plot(y=dat_AQ$PM2.5_mean,x=dat_AQ$Date,type='l',main = "PM2.5 mean",xlab = " ",ylab = " ")




# Weather data ------------------------------------------------------------
milano_stat=c(100,501,502,503,620,869,1327,1425)

dat_w=get_ARPA_Lombardia_W_data(
  ID_station = milano_stat,
  Date_begin = "2023-12-01",
  Date_end = "2024-06-30",
  Frequency = "daily",
  Var_vec = NULL,
  Fns_vec = NULL,
  by_sensor = FALSE,
  verbose = TRUE,
  parallel = T,
  parworkers = NULL,
  parfuturetype = "multisession"
)

library(dplyr)
dat_av_w=dat_w%>%
  group_by(Date)%>%
  summarise_if(is.numeric,mean,na.rm = T)

dat_av_w=subset(dat_av_w,select = -c(IDStation,Water_height))

# Substitute NaN with NA
dat_av_w=dat_av_w%>%
  mutate_all(~replace(.,is.nan(.),NA))

windows()
par(mfrow=c(2,3))
plot(y=dat_av_w$Temperature,x=dat_av_w$Date,type='l',main = "Temperature",xlab = " ",ylab = " ")
plot(y=dat_av_w$Relative_humidity,x=dat_av_w$Date,type='l',main = "Humidity",xlab = " ",ylab = " ")
plot(y=dat_av_w$Wind_speed,x=dat_av_w$Date,type='l',main = "Wind Speed",xlab = " ",ylab = " ")
plot(y=dat_av_w$Wind_direction,x=dat_av_w$Date,type='l',main = "Wind Direction",xlab = " ",ylab = " ")
plot(y=dat_av_w$Global_radiation,x=dat_av_w$Date,type='l',main = "Global Radiation",xlab = " ",ylab = " ")
plot(y=dat_av_w$Rainfall,x=dat_av_w$Date,type='l',main = "Rainfall",xlab = " ",ylab = " ")

# Drop whetever rainfall<0
dat_av_w$Rainfall[dat_av_w$Rainfall<0]=NA

# Drop temperature outliers
wdn=171:175
dat_av_w$Temperature[wdn]=NA

# Drop wind speed outlier
dat_av_w$Wind_speed[dat_av_w$Wind_speed>15]=NA

windows()
par(mfrow=c(2,3))
plot(y=dat_av_w$Temperature,x=dat_av_w$Date,type='l',main = "Temperature",xlab = " ",ylab = " ")
plot(y=dat_av_w$Relative_humidity,x=dat_av_w$Date,type='l',main = "Humidity",xlab = " ",ylab = " ")
plot(y=dat_av_w$Wind_speed,x=dat_av_w$Date,type='l',main = "Wind Speed",xlab = " ",ylab = " ")
plot(y=dat_av_w$Wind_direction,x=dat_av_w$Date,type='l',main = "Wind Direction",xlab = " ",ylab = " ")
plot(y=dat_av_w$Global_radiation,x=dat_av_w$Date,type='l',main = "Global Radiation",xlab = " ",ylab = " ")
plot(y=dat_av_w$Rainfall,x=dat_av_w$Date,type='l',main = "Rainfall",xlab = " ",ylab = " ")


# Merge data --------------------------------------------------------------

dat=merge(dat_AQ,dat_av_w,by="Date")

save(dat,file="dat.RData")

# Feature construction ----------------------------------------------------

# Day of the week
weekday=c(3:7,1:2)
res=dim(dat)[1]%%7
weekday=rep(weekday,dim(dat)[1]/7)

weekday=c(weekday,tail(weekday,1)+(1:res))

dat$weekday=weekday

# Is weekend?
dat$weekend=ifelse(dat$weekday%in%c(6,7),1,0)

# Rainy day?
dat$rainy=ifelse(dat$Rainfall>0,1,0)

par(mfrow=c(1,2))
hist(dat$Wind_speed,xlab = " ",main = "Wind Speed")
hist(dat$Rainfall,xlab = " ",main = "Rainfall")


# First differences -------------------------------------------------------

# dat$NO2_mean_diff=c(NA,diff(dat$NO2_mean))
# dat$O3_max_day_diff=c(NA,diff(dat$Ozone_max_day))
# dat$PM10_mean_diff=c(NA,diff(dat$PM10_mean))
# dat$PM25_mean_diff=c(NA,diff(dat$PM2.5_mean))
# dat$Temperature_diff=c(NA,diff(dat$Temperature))
# dat$Relative_humidity_diff=c(NA,diff(dat$Relative_humidity))
# dat$Wind_speed_diff=c(NA,diff(dat$Wind_speed))
# dat$Wind_direction_diff=c(NA,diff(dat$Wind_direction))
# dat$Global_radiation_diff=c(NA,diff(dat$Global_radiation))
# dat$Rainfall_diff=c(NA,diff(dat$Rainfall))
# 
# # Plot
# windows()
# par(mfrow=c(2,5))
# plot(y=dat$NO2_mean_diff,x=dat$Date,type='l',main = "NO2 mean diff",xlab = " ",ylab = " ")
# plot(y=dat$O3_max_day_diff,x=dat$Date,type='l',main = "O3 max day diff",xlab = " ",ylab = " ")
# plot(y=dat$PM10_mean_diff,x=dat$Date,type='l',main = "PM10 mean diff",xlab = " ",ylab = " ")
# plot(y=dat$PM25_mean_diff,x=dat$Date,type='l',main = "PM2.5 mean diff",xlab = " ",ylab = " ")
# plot(y=dat$Temperature_diff,x=dat$Date,type='l',main = "Temperature diff",xlab = " ",ylab = " ")
# plot(y=dat$Relative_humidity_diff,x=dat$Date,type='l',main = "Humidity diff",xlab = " ",ylab = " ")
# plot(y=dat$Wind_speed_diff,x=dat$Date,type='l',main = "Wind Speed diff",xlab = " ",ylab = " ")
# plot(y=dat$Wind_direction_diff,x=dat$Date,type='l',main = "Wind Direction diff",xlab = " ",ylab = " ")
# plot(y=dat$Global_radiation_diff,x=dat$Date,type='l',main = "Global Radiation diff",xlab = " ",ylab = " ")
# plot(y=dat$Rainfall_diff,x=dat$Date,type='l',main = "Rainfall diff",xlab = " ",ylab = " ")


# Other features ----------------------------------------------------------

library(zoo)
wdn=7
dat$corr_NO2_temp=c(rep(NA,wdn-1),rollapply(apply(dat[,c("NO2_mean","Temperature")],2,scale),
                                                       width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                       by.column=FALSE))
dat$corr_NO2_humidity=c(rep(NA,wdn-1),rollapply(apply(dat[,c("NO2_mean","Relative_humidity")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
dat$corr_NO2_windspeed=c(rep(NA,wdn-1),rollapply(apply(dat[,c("NO2_mean","Wind_speed")],2,scale),
                                        width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                        by.column=FALSE))
dat$corr_NO2_radiation=c(rep(NA,wdn-1),rollapply(apply(dat[,c("NO2_mean","Global_radiation")],2,scale),
                                         width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                         by.column=FALSE))
dat$corr_NO2_rainfall=c(rep(NA,wdn-1),rollapply(apply(dat[,c("NO2_mean","Rainfall")],2,scale),
                                        width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                        by.column=FALSE))
dat$corr_NO2_winddir=c(rep(NA,wdn-1),rollapply(apply(dat[,c("NO2_mean","Wind_direction")],2,scale),
                                        width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                        by.column=FALSE))

dat$corr_O3_temp=c(rep(NA,wdn-1),rollapply(apply(dat[,c("Ozone_max_day","Temperature")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
dat$corr_O3_humidity=c(rep(NA,wdn-1),rollapply(apply(dat[,c("Ozone_max_day","Relative_humidity")],2,scale),
                                                width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                by.column=FALSE))
dat$corr_O3_windspeed=c(rep(NA,wdn-1),rollapply(apply(dat[,c("Ozone_max_day","Wind_speed")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
dat$corr_O3_radiation=c(rep(NA,wdn-1),rollapply(apply(dat[,c("Ozone_max_day","Global_radiation")],2,scale),
                                             width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                             by.column=FALSE))
dat$corr_O3_rainfall=c(rep(NA,wdn-1),rollapply(apply(dat[,c("Ozone_max_day","Rainfall")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
dat$corr_O3_winddir=c(rep(NA,wdn-1),rollapply(apply(dat[,c("Ozone_max_day","Wind_direction")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))

dat$corr_PM10_temp=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM10_mean","Temperature")],2,scale),
                                              width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                              by.column=FALSE))
dat$corr_PM10_humidity=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM10_mean","Relative_humidity")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))
dat$corr_PM10_windspeed=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM10_mean","Wind_speed")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))
dat$corr_PM10_radiation=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM10_mean","Global_radiation")],2,scale),
                                                   width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                   by.column=FALSE))
dat$corr_PM10_rainfall=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM10_mean","Rainfall")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))
dat$corr_PM10_winddir=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM10_mean","Wind_direction")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))

dat$corr_PM25_temp=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM2.5_mean","Temperature")],2,scale),
                                              width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                              by.column=FALSE))
dat$corr_PM25_humidity=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM2.5_mean","Relative_humidity")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))
dat$corr_PM25_windspeed=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM2.5_mean","Wind_speed")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))
dat$corr_PM25_radiation=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM2.5_mean","Global_radiation")],2,scale),
                                                   width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                   by.column=FALSE))
dat$corr_PM25_rainfall=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM2.5_mean","Rainfall")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))
dat$corr_PM25_winddir=c(rep(NA,wdn-1),rollapply(apply(dat[,c("PM2.5_mean","Wind_direction")],2,scale),
                                                  width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                                  by.column=FALSE))


# Final check -------------------------------------------------------------

# Transform weekday, weekend and rainy into factors
dat$weekday=factor(dat$weekday)
dat$weekend=factor(dat$weekend)
dat$rainy=factor(dat$rainy)

str(dat)

save(dat,file="dat.RData")

# AQ Index ----------------------------------------------------------------

library(con2aqi)
con2aqi("no2",con=dat$NO2_mean)
con2aqi("pm25",con=dat$PM2.5_mean)
con2aqi("pm10",con=dat$PM10_mean)

# Ozone: 1 ppb=1.96 mugr/m3

o3_aqi=as.vector(con2aqi("o3",con=dat$Ozone_max_8h*0.00051,type="8h"))
no2_aqi=as.vector(con2aqi("no2",con=dat$NO2_mean))
pm25_aqi=as.vector(con2aqi("pm25",con=dat$PM2.5_mean))
pm10_aqi=as.vector(con2aqi("pm10",con=dat$PM10_mean))

all_aqi=cbind(o3_aqi,no2_aqi,pm25_aqi,pm10_aqi)
AQI=apply(all_aqi,1,max)
AQI_fact=rep(0,length(AQI))
for(i in 1:length(AQI)){
  if(AQI[i]<=50){
    AQI_fact[i]=1
  }else if(AQI[i]>50&AQI[i]<=100){
    AQI_fact[i]=2
  }else if(AQI[i]>100&AQI[i]<=150){
    AQI_fact[i]=3
  }else if(AQI[i]>150&AQI[i]<=200){
    AQI_fact[i]=4
  }else if(AQI[i]>200&AQI[i]<=300){
    AQI_fact[i]=5
  }else{
    AQI_fact[i]=6
  }
}

AQI_fact=factor(AQI_fact,levels=1:6)
save(AQI_fact,file="AQI_fact.RData")
