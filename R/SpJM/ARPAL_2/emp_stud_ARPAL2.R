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

# # Lagged rainfall
# data$rainfall_lag1=lag(data$rainfall,1)
# data$rainfallmax_lag1=lag(data$rainfall_max,1)
# data$rainfallmin_lag1=lag(data$rainfall_min,1)

# Plot data
windows()
par(mfrow=c(2,2))
plot(data$Date,data$pm25,type="l",xlab="Date",ylab="pm25")
plot(data$Date,data$pm10,type="l",xlab="Date",ylab="pm10")
plot(data$Date,data$o3,type="l",xlab="Date",ylab="o3")
plot(data$Date,data$no2,type="l",xlab="Date",ylab="no2")

# GGplot
library(ggplot2)
ppm25=ggplot(data,aes(x=Date,y=pm25))+geom_line()+xlab("Date")+ylab("PM2.5")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))
ppm10=ggplot(data,aes(x=Date,y=pm10))+geom_line()+xlab("Date")+ylab("PM10")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))
po3=ggplot(data,aes(x=Date,y=o3))+geom_line()+xlab("Date")+ylab("O3")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))
pno2=ggplot(data,aes(x=Date,y=no2))+geom_line()+xlab("Date")+ylab("NO2")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))

library(gridExtra)
windows()
grid.arrange(ppm25,ppm10,po3,pno2,ncol=2)

windows()
par(mfrow=c(3,2))
plot(data$Date,data$temp,type="l",xlab="Date",ylab="temp")
plot(data$Date,data$rel_hum,type="l",xlab="Date",ylab="rel_hum")
plot(data$Date,data$wind_speed,type="l",xlab="Date",ylab="wind_speed")
plot(data$Date,data$rainfall,type="l",xlab="Date",ylab="rainfall")
plot(data$Date,data$globrad,type="l",xlab="Date",ylab="glob rad")

#ggplot
ptemp=ggplot(data,aes(x=Date,y=temp))+geom_line()+xlab("Date")+ylab("Temp")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))
prel_hum=ggplot(data,aes(x=Date,y=rel_hum))+geom_line()+xlab("Date")+ylab("RH")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))
pwind_speed=ggplot(data,aes(x=Date,y=wind_speed))+geom_line()+xlab("Date")+ylab("WS")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))
prainfall=ggplot(data,aes(x=Date,y=rainfall))+geom_line()+xlab("Date")+ylab("RF")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))

pglobrad=ggplot(data,aes(x=Date,y=globrad))+geom_line()+xlab("Date")+ylab("GR")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))


library(gridExtra)
windows()
grid.arrange(ptemp,prel_hum,pwind_speed,prainfall,pglobrad,ncol=2)


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
data$rainy=ifelse(data$rainfall>0,1,0)
#data$rainy=ifelse(data$rainfall>mean(data$rainfall,na.rm = T),1,0)
data$rainy=as.factor(data$rainy)
str(data)

# Windy day?
#data$windy=ifelse(data$wind_speed>0.5,1,0)
data$windy=ifelse(data$wind_speed>median(data$wind_speed,na.rm = T),1,0)
#data$windy=ifelse(data$wind_speed>1.7,1,0)

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



# Remove rainfall and windspeed -------------------------------------------

#data=subset(data,select=-c(rainfall,rainfall_min,rainfall_max,wind_speed,wind_speed_min,wind_speed_max))

save(data,file="ARPAL_2/data.Rdata")

dat=data
save(dat,file="ARPAL_2/data.Rdata")


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
save(AQI_fact,file="ARPAL_2/AQI_fact.RData")

plot(x=dat$Date,as.numeric(AQI_fact),type='l',xlab="Date",ylab="AQI",main="AQI")

#with ggplot2
library(ggplot2)
ggplot(dat,aes(x=Date,y=as.numeric(AQI_fact)))+geom_line()+xlab("Date")+ylab("AQI")+ggtitle("AQI")+scale_x_date(date_breaks = "3 month",date_labels = "%b %Y")

# JM mix ------------------------------------------------------------------

load("ARPAL_2/data.Rdata")
load("ARPAL_2/AQI_fact.Rdata")

# Percentage of missing values

round(colMeans(is.na(dat)) * 100, 2)
windows()
Amelia::missmap(dat)

round(colMeans(is.na(dat)) * 100, 2)

Amelia::missmap(dat)

# Summary statistics
dat%>%summarise_if(is.numeric,mean,na.rm=T)
dat%>%summarise_if(is.numeric,sd,na.rm=T)
dat%>%summarise_if(is.numeric,min,na.rm=T)
dat%>%summarise_if(is.numeric,max,na.rm=T)

apply(dat[,c("windy","rainy","weekend","holiday","month")],2,table)
#NAs % of categorical variables
sum(which(is.na(dat$windy)))/nrow(dat)
sum(which(is.na(dat$rainy)))/nrow(dat)
sum(which(is.na(dat$weekend)))/nrow(dat)


  # Correlation plot
windows()

#corrplot::corrplot(cor(data[,c(2:5,8,11,14,53:54)],use="complete.obs"),method="number")

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
#K=2:6
K=4
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

save(est,file="est.Rdata")

# est=lapply(lambda,function(l){
#   jump_mixed2(dat_notime, 4, jump_penalty=l, 
#               initial_states=NULL,
#               max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
#               timeflag=F
#   )
# })

library(pdfCluster)
ARI_res=unlist(lapply(est,function(e){adj.rand.index(e$best_s,AQI_fact)}))

source("Utils.R")
GIC=unlist(lapply(est,function(e){GIC_mixed(e$Y,e$best_s,est[[1]]$best_s,K=e$K,K0=3,pers0 = .9)$FTIC}))

res=data.frame(ARI_res,GIC,lambda,K)

plot(res$lambda,res$ARI_res,type="l",xlab="lambda",ylab="ARI",main="ARI vs lambda")

best_est=est[[7]]

table(best_est$best_s)

best_est$condMM

states=recode(best_est$best_s, "3" = "Good", "2" = "Moderate","1"="US","4"="Unhealthy")
states=factor(states,levels=c("Good","Moderate","US","Unhealthy"))
#states=factor(best_est$best_s,levels=1:4,labels=c("Good","US","Unhealthy","Moderate"))
true_states=factor(AQI_fact,levels=1:4,labels=c("Good","Moderate","US","Unhealthy"))

adj.rand.index(states,true_states)
sum(states==true_states)/length(states)


# AQI
P1=ggplot(dat,aes(x=Date,y=as.numeric(true_states)))+
  geom_line(size=1)+xlab("Date")+ylab(" ")+ggtitle("AQI")+
  scale_x_date(date_breaks = "3 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))

P2=ggplot(dat,aes(x=Date,y=as.numeric(states)))+
  geom_line(size=1)+xlab("Date")+ylab(" ")+ggtitle("JM-mix")+
  scale_x_date(date_breaks = "3 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))

library(gridExtra)
windows()
grid.arrange(P1,P2,ncol=1)

library(caret)
confusionMatrix(states,true_states)

#AQI plot##

data_state=data.frame(dat,State=states,AQI=true_states)

ppm25=ggplot(data_state,aes(x=Date,y=pm25))+geom_line()+xlab("Date")+ylab("PM2.5")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

ppm10=ggplot(data_state,aes(x=Date,y=pm10))+geom_line()+xlab("Date")+ylab("PM10")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

po3=ggplot(data_state,aes(x=Date,y=o3))+geom_line()+xlab("Date")+ylab("O3")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

pno2=ggplot(data_state,aes(x=Date,y=no2))+geom_line()+xlab("Date")+ylab("NO2")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))


library(ggpubr)
windows()
#grid.arrange(ppm25,ppm10,po3,pno2,ncol=2)
ggarrange(ppm25,ppm10,po3,pno2, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

#ggplot
ptemp=ggplot(data_state,aes(x=Date,y=temp))+geom_line()+xlab("Date")+ylab("Temp")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

prel_hum=ggplot(data_state,aes(x=Date,y=rel_hum))+geom_line()+xlab("Date")+ylab("RH")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

pwind_speed=ggplot(data_state,aes(x=Date,y=wind_speed))+geom_line()+xlab("Date")+ylab("WS")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

prainfall=ggplot(data_state,aes(x=Date,y=rainfall))+geom_line()+xlab("Date")+ylab("RF")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))

pglobrad=ggplot(data_state,aes(x=Date,y=globrad))+geom_line()+xlab("Date")+ylab("GR")+
  scale_x_date(date_breaks = "6 month",date_labels = "%b %Y")+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=22))+
  geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), 
                ymin = -Inf, ymax = Inf, fill = State), alpha = .2) +
  scale_fill_manual(values = alpha(c("green","yellow", "#FF6600", "#FF0033")))


library(gridExtra)
windows()
#grid.arrange(ptemp,prel_hum,pwind_speed,prainfall,pglobrad,ncol=2)
ggarrange(ptemp,prel_hum,pwind_speed,prainfall,pglobrad, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")


# State conditional correlations
# State conditional correlation
dat_good=data_state[data_state$State=="Good",]
dat_moderate=data_state[data_state$State=="Moderate",]
dat_us=data_state[data_state$State=="US",]
dat_unhealthy=data_state[data_state$State=="Unhealthy",]

# Extract numeric variables
cont_feat=c(2:5,8,11,14,17,20)

dat_good=dat_good[cont_feat]
dat_moderate=dat_moderate[cont_feat]
dat_us=dat_us[cont_feat]
dat_unhealthy=dat_unhealthy[cont_feat]

cor_good=cor(dat_good,use="complete.obs")
cor_moderate=cor(dat_moderate,use="complete.obs")
cor_us=cor(dat_us,use="complete.obs")
cor_unhealthy=cor(dat_unhealthy,use="complete.obs")

round(cor_good,2)
round(cor_moderate,2)
round(cor_us,2)
round(cor_unhealthy,2)
