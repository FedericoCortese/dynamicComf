library(lubridate)
pm25=read.csv("ARPAL_2/pm25.csv")
colnames(pm25)=c("Date","pm25")
pm25$Date=ymd(pm25$Date)
pm25$pm25=as.numeric(pm25$pm25)
str(pm25)

pm10=read.csv("ARPAL_2/pm10.csv")
colnames(pm10)=c("Date","pm10")
pm10$Date=ymd(pm10$Date)
pm10$pm10=as.numeric(pm10$pm10)
str(pm10)

o3=read.csv("ARPAL_2/o3.csv")
colnames(o3)=c("Date","o3")
o3$Date=ymd(o3$Date)
o3$o3=as.numeric(o3$o3)
str(o3)

no2=read.csv("ARPAL_2/no2.csv")
colnames(no2)=c("Date","no2")
no2$Date=ymd(no2$Date)
no2$no2=as.numeric(no2$no2)
str(no2)

no2max=read.csv("ARPAL_2/no2max.csv")
colnames(no2max)=c("Date","no2max")
no2max$Date=ymd(no2max$Date)
no2max$no2max=as.numeric(no2max$no2max)
str(no2max)

o3max=read.csv("ARPAL_2/o3max.csv")
colnames(o3max)=c("Date","o3max")
o3max$Date=ymd(o3max$Date)
o3max$o3max=as.numeric(o3max$o3max)
str(o3max)

temp=read.csv("ARPAL_2/temp.csv")
temp=subset(temp,select=-Id.Sensore)
colnames(temp)=c("Date","temp","temp_min","temp_max")
temp$Date=ymd(temp$Date)
temp[temp==-999]=NA
str(temp)

rel_hum=read.csv("ARPAL_2/rel_hum.csv")
rel_hum=subset(rel_hum,select=-Id.Sensore)
colnames(rel_hum)=c("Date","rel_hum","rel_hum_min","rel_hum_max")
rel_hum$Date=ymd(rel_hum$Date)
rel_hum[rel_hum==-999]=NA
str(rel_hum)

wind_speed=read.csv("ARPAL_2/windspeed.csv")
wind_speed=subset(wind_speed,select=-Id.Sensore)
colnames(wind_speed)=c("Date","wind_speed","wind_speed_min","wind_speed_max")
wind_speed$Date=ymd(wind_speed$Date)
wind_speed[wind_speed==-999]=NA
str(wind_speed)

rainfall=read.csv("ARPAL_2/rainfall.csv")
rainfall=subset(rainfall,select=-Id.Sensore)
colnames(rainfall)=c("Date","rainfall","rainfall_min","rainfall_max")
rainfall$Date=ymd(rainfall$Date)
rainfall[rainfall==-999]=NA
str(rainfall)

globrad=read.csv("ARPAL_2/glob_rad.csv")
globrad=subset(globrad,select=-Id.Sensore)
colnames(globrad)=c("Date","globrad","globrad_min","globrad_max")
globrad$Date=ymd(globrad$Date)
globrad[globrad==-999]=NA
str(globrad)

# Merge all data
data=merge(pm25,pm10,by="Date",all=TRUE)
data=merge(data,o3,by="Date",all=TRUE)
data=merge(data,no2,by="Date",all=TRUE)
data=merge(data,no2max,by="Date",all=TRUE)
data=merge(data,o3max,by="Date",all=TRUE)
data=merge(data,temp,by="Date",all=TRUE)
data=merge(data,rel_hum,by="Date",all=TRUE)
data=merge(data,wind_speed,by="Date",all=TRUE)
data=merge(data,rainfall,by="Date",all=TRUE)
data=merge(data,globrad,by="Date",all=TRUE)
str(data)

# Plot data
windows()
par(mfrow=c(2,2))
plot(data$Date,data$pm25,type="l",xlab="Date",ylab="pm25")
plot(data$Date,data$pm10,type="l",xlab="Date",ylab="pm10")
plot(data$Date,data$o3,type="l",xlab="Date",ylab="o3")
plot(data$Date,data$no2,type="l",xlab="Date",ylab="no2")

windows()
par(mfrow=c(2,2))
plot(data$Date,data$temp,type="l",xlab="Date",ylab="temp")
plot(data$Date,data$rel_hum,type="l",xlab="Date",ylab="rel_hum")
plot(data$Date,data$wind_speed,type="l",xlab="Date",ylab="wind_speed")
plot(data$Date,data$rainfall,type="l",xlab="Date",ylab="rainfall")
plot(data$Date,data$globrad,type="l",xlab="Date",ylab="glob rad")


# Save data
save(data,file="ARPAL_2/data.Rdata")



# Features ----------------------------------------------------------------

# Weekday
data$weekday=weekdays(data$Date)
data$weekday=as.factor(data$weekday)
data$weekday=recode_factor(data$weekday,"lunedì"="1","martedì"="2","mercoledì"="3","giovedì"="4",
                           "venerdì"="5","sabato"="6","domenica"="7")
str(data)

# Month 
data$month=month(data$Date)
data$month=as.factor(data$month)
str(data)

# Weekend?
data$weekend=ifelse(data$weekday=="6" | data$weekday=="7",1,0)
data$weekend=as.factor(data$weekend)
str(data)

# Holiday
holidays=c("01-01",
           "01-06",
           "04-25",
           "05-01",
           "06-02",
           "08-15",
           "11-01",
           "12-07",
           "12-08",
           "12-25",
           "12-26",
           "2021-04-04",
           "2022-04-17",
           "2023-04-09",
           "2024-03-31")
data$holiday=0
for(i in 1:length(data$Date)){
  for(j in 1:length(holidays)){
    if(grepl(holidays[j], as.character(data$Date[i]))){
      data$holiday[i]=1
    }
  }
}
data$holiday=factor(data$holiday)
str(data)

# Rainy day?
#data$rainy=ifelse(data$rainfall>0,1,0)
data$rainy=ifelse(data$rainfall>mean(data$rainfall,na.rm = T),1,0)
data$rainy=as.factor(data$rainy)
str(data)

# Windy day?
data$windy=ifelse(data$wind_speed>median(data$wind_speed,na.rm = T),1,0)
data$windy=as.factor(data$windy)
str(data)

data=subset(data,select=-weekday)

# 7-days Correlations
library(zoo)
wdn=7
data$corr_no2_temp=c(rep(NA,wdn-1),rollapply(apply(data[,c("no2","temp")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_no2_relhum=c(rep(NA,wdn-1),rollapply(apply(data[,c("no2","rel_hum")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_no2_windspeed=c(rep(NA,wdn-1),rollapply(apply(data[,c("no2","wind_speed")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_no2_rainfall=c(rep(NA,wdn-1),rollapply(apply(data[,c("no2","rainfall")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_no2_globrad=c(rep(NA,wdn-1),rollapply(apply(data[,c("no2","globrad")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))

data$corr_pm25_temp=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm25","temp")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm25_relhum=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm25","rel_hum")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm25_windspeed=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm25","wind_speed")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm25_rainfall=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm25","rainfall")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm25_globrad=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm25","globrad")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))

data$corr_pm10_temp=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm10","temp")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm10_relhum=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm10","rel_hum")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm10_windspeed=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm10","wind_speed")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm10_rainfall=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm10","rainfall")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_pm10_globrad=c(rep(NA,wdn-1),rollapply(apply(data[,c("pm10","globrad")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))

data$corr_o3_temp=c(rep(NA,wdn-1),rollapply(apply(data[,c("o3","temp")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_o3_relhum=c(rep(NA,wdn-1),rollapply(apply(data[,c("o3","rel_hum")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_o3_windspeed=c(rep(NA,wdn-1),rollapply(apply(data[,c("o3","wind_speed")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_o3_rainfall=c(rep(NA,wdn-1),rollapply(apply(data[,c("o3","rainfall")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))
data$corr_o3_globrad=c(rep(NA,wdn-1),rollapply(apply(data[,c("o3","globrad")],2,scale),
                                            width=wdn, function(x) cor(x[,1],x[,2],use="na.or.complete"), 
                                            by.column=FALSE))


# 7-days MA ---------------------------------------------------------------

wdn=7
data$ma_no2=rollapply(data$no2,width=wdn,mean,align="right",fill=NA)
data$ma_pm25=rollapply(data$pm25,width=wdn,mean,align="right",fill=NA)
data$ma_pm10=rollapply(data$pm10,width=wdn,mean,align="right",fill=NA)
data$ma_o3=rollapply(data$o3,width=wdn,mean,align="right",fill=NA)
data$ma_temp=rollapply(data$temp,width=wdn,mean,align="right",fill=NA)
data$ma_rel_hum=rollapply(data$rel_hum,width=wdn,mean,align="right",fill=NA)
data$ma_wind_speed=rollapply(data$wind_speed,width=wdn,mean,align="right",fill=NA)
data$ma_rainfall=rollapply(data$rainfall,width=wdn,mean,align="right",fill=NA)
data$ma_globrad=rollapply(data$globrad,width=wdn,mean,align="right",fill=NA)
dat=data
save(dat,file="ARPAL_2/data.Rdata")

Amelia::missmap(data)

# AQI ---------------------------------------------------------------------

library(con2aqi)

o3_aqi=as.vector(con2aqi("o3",con=data$o3*0.00051,type="8h"))
no2_aqi=as.vector(con2aqi("no2",con=data$no2))
pm25_aqi=as.vector(con2aqi("pm25",con=data$pm25))
pm10_aqi=as.vector(con2aqi("pm10",con=data$pm10))

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


# JM mix ------------------------------------------------------------------

load("ARPAL_2/data.Rdata")
load("ARPAL_2/AQI_fact.Rdata")

str(data)

# Percentage of missing values
round(colMeans(is.na(dat)) * 100, 2)

Amelia::missmap(dat)

# Summary statistics
summary(dat)

# Correlation plot
windows()
corrplot::corrplot(cor(dat[,c(2:5,8,11,14,17,20)],use="complete.obs"),method="number")

tapply(dat$pm25,dat$windy,mean)
tapply(dat$pm25,dat$rainy,mean)
tapply(dat$pm25,dat$weekend,mean)
tapply(dat$pm25,dat$weekday,mean)
tapply(dat$pm25,dat$holiday,mean)

tapply(dat$pm10,dat$windy,mean)
tapply(dat$pm10,dat$rainy,mean)
tapply(dat$pm10,dat$weekend,mean)
tapply(dat$pm10,dat$weekday,mean)
tapply(dat$pm10,dat$holiday,mean)

tapply(dat$o3,dat$windy,mean)
tapply(dat$o3,dat$rainy,mean)
tapply(dat$o3,dat$weekend,mean)
tapply(dat$o3,dat$weekday,mean)
tapply(dat$o3,dat$holiday,mean)

tapply(dat$no2,dat$windy,mean)
tapply(dat$no2,dat$rainy,mean)
tapply(dat$no2,dat$weekend,mean)
tapply(dat$no2,dat$weekday,mean)
tapply(dat$no2,dat$holiday,mean)

source("Utils.R")

dat_notime=dat[,-1]

lambda=seq(0,1,by=.05)
K=2:6
hp=expand.grid(K=K,lambda=lambda)

start_=Sys.time()
est <- parallel::mclapply(1:nrow(hp),
                          function(x)
                            jump_mixed2(dat_notime, 
                                        n_states=hp[x,]$K, 
                                        jump_penalty=hp[x,]$lambda, 
                                        initial_states=NULL,
                                        max_iter=10, n_init=10, tol=NULL, 
                                        verbose=FALSE,
                                        timeflag=F
                            ),
                          mc.cores = parallel::detectCores()-1)

end_=Sys.time()
elapsed=end_-start_

# est=lapply(lambda,function(l){
#   jump_mixed2(dat_notime, 4, jump_penalty=l, 
#               initial_states=NULL,
#               max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
#               timeflag=F
#   )
# })


ARI_res=unlist(lapply(est,function(e){adj.rand.index(e$best_s,AQI_fact)}))

source("Utils.R")
GIC=unlist(lapply(est,function(e){GIC_mixed(e$Y,e$best_s,est[[1]]$best_s,K=e$K)$FTIC}))
res=data.frame(ARI_res,GIC,lambda,K)

plot(res$lambda,res$ARI_res,type="l",xlab="lambda",ylab="ARI",main="ARI vs lambda")

best_est=est[[16]]

table(best_est$best_s)

best_est$condMM

states=factor(best_est$best_s,levels=1:4,labels=c("Good","Moderate","US","Unhealthy"))
true_states=factor(AQI_fact,levels=1:4,labels=c("Good","Moderate","US","Unhealthy"))


library(caret)
confusionMatrix(states,true_states)