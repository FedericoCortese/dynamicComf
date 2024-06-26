library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)

# Data cleaning -----------------------------------------------------------

enth_surv=read.csv("enth_surveys_calc.csv" )
enth_tab=read.csv("enth_tabular_merged.csv")

str(enth_surv)
str(enth_tab)

# space_id, indoor_latitude and indoor_longitude can be used for parsing some environmental features
# Remove unnecessary columns
enth_tab=subset(enth_tab,select=-c(space_id,
                                   building_name,
                                   response_speed,
                                   indoor_floor,
                                   body_presence,
                                   user_id,
                                   indoor_latitude,
                                   indoor_longitude,
                                   indoor_floor,
                                   change
))

enth_tab$time=as.POSIXct(enth_tab$time,format="%Y-%m-%d %H:%M:%S")

str(enth_tab)

#wdn="15 mins"
wdn="5 mins"
enth_tab2=enth_tab%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE) 

count_consecutive=function(x,tol=4){
  d_t=diff(x$time)/60
  counter=d_t<tol
  av=ave(counter, cumsum(!counter), FUN = cumsum)
  return(list(av=av,
              max.t=max(av),
              t=which.max(av)))
}

# "Average" person
enth_tab_av=enth_tab2%>%group_by(time)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

cc=count_consecutive(enth_tab_av)
#cbind(cc$av,enth_tab_av$time)
sort(cc$av,decreasing = TRUE)
cc$max.t
cc$t
# Waaaaay better
#enth_tab_av$time[(156-56):155]
#enth_tab_av$time[(428-45+1):(428+1)] # 15 mins
enth_tab_av$time[(194-80):194-1] # 5 mins
#wdn2=(cc$t-1-cc$max.t):cc$t
#wdn2=(156-56):155
wdn2=(194-80):194-1 #5 mins

#wdn2=(428-45+1):(428+1) # 15 mins

enth_tab3=enth_tab_av[wdn2,]
Amelia::missmap(enth_tab3)

str(enth_tab3)

# NAs percentage for each variable
apply(enth_tab3,2,function(x) sum(is.na(x))/length(x))*100

# Remove vars with 100% NAs
enth_tab3=subset(enth_tab3,select=-c(co2_indoor,voc_indoor,pm25_indoor,noise_indoor))

Amelia::missmap(enth_tab3)
summary(enth_tab3)

colnames(enth_tab3)

# Indoor or outdoor?
old_io=enth_tab3$indoor.outdoor
dst=enth_tab3$indoor.outdoor-rep(9,length(enth_tab3$indoor.outdoor))

for(i in 1:length(dst)){
  if(dst[i]==1){
    enth_tab3$indoor.outdoor[i]=NA
  }else if(dst[i]<1){
    enth_tab3$indoor.outdoor[i]=9
  }
  else{
    enth_tab3$indoor.outdoor[i]=11
  }
}

# Thermal preferences
old_th=enth_tab3$thermal
dst_th=enth_tab3$thermal-rep(10,length(enth_tab3$thermal))

for(i in 1:length(enth_tab3$thermal)){
  if(dst_th[i]==.5){
    enth_tab3$thermal[i]=NA
  }else if(dst_th[i]<.5){
    enth_tab3$thermal[i]=10
  }
  else{
    enth_tab3$thermal[i]=11
  }
}

# Create temp, humidity, pressure variables

enth_tab3[,-1]=enth_tab3[,-1]%>% mutate_all(~ifelse(is.nan(.), NA, .))
enth_tab3$temp=NA
enth_tab3$humidity=NA
enth_tab3$pressure=NA

for(i in 1:length(enth_tab3$indoor.outdoor)){
  if(is.na(enth_tab3$indoor.outdoor[i])){
    
    enth_tab3$temp[i]=mean(enth_tab3$temp_zone[i],enth_tab3$temp_outdoor[i],na.rm=TRUE,trim=.5)
    # i=39
    # Errore in if (trim > 0 && n) { : 
    # valore mancante dove è richiesto TRUE/FALSE
    
    enth_tab3$humidity[i]=mean(enth_tab3$humidity_zone[i],enth_tab3$humidity_outdoor[i],na.rm=TRUE,trim=.5)
    #enth_tab3$pressure[i]=NA
  }
  else if(enth_tab3$indoor.outdoor[i]==9){
    enth_tab3$temp[i]=enth_tab3$temp_zone[i]
    enth_tab3$humidity[i]=enth_tab3$humidity_zone[i]
    enth_tab3$pressure[i]=enth_tab3$pressure_outdoor[i]
  }
  else{
    enth_tab3$temp[i]=enth_tab3$temp_outdoor[i]
    enth_tab3$humidity[i]=enth_tab3$humidity_outdoor[i]
    enth_tab3$pressure[i]=enth_tab3$pressure_outdoor[i]
  }
}

enth_tab4=select(enth_tab3,subset=-c(temp_zone,temp_outdoor,humidity_zone,humidity_outdoor,pressure_outdoor,
                                     X0.3um_count_outdoor, X0.5um_count_outdoor, 
                                     X1.0um_count_outdoor, X10.0um_count_outdoor, X5.0um_count_outdoor,
                                     X2.5um_count_outdoor))

enth_tab4$met=round(enth_tab4$met)
enth_tab4$comfort=round(enth_tab4$comfort)
enth_tab4[,c("comfort","thermal","indoor.outdoor","met")]=enth_tab4[,c("comfort","thermal","indoor.outdoor","met")]%>%
  mutate_all(as.factor)
head(enth_tab4)
str(enth_tab4)
save(enth_tab4,file="enth_tab4.RData")

Amelia::missmap(enth_tab4)
# source("Utils.R")
# lambda=.2
# n_states=2
# enth_tab5=enth_tab4
# enth_tab5=enth_tab5%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
# 
# est=jump_mixed(enth_tab5[,-1], n_states, jump_penalty=lambda, 
#                initial_states=NULL,
#                max_iter=10, n_init=10, tol=NULL, verbose=T) 
# 
# res=data.frame(enth_tab4,state=est$best_s)
# tapply(res$temp,res$state,mean,na.rm=T)
# tapply(res$humidity,res$state,mean,na.rm=T)
# tapply(res$pressure,res$state,mean,na.rm=T)
# 
# # Included? Maybe not
# table(res$thermal,res$state)
# 
# # No comfort, only thermal
# table(res$comfort,res$state)
# 
# tapply(res$met,res$state,Mode,na.rm=T)
# tapply(res$indoor.outdoor,res$state,Mode,na.rm=T)
# tapply(res$heartrate,res$state,mean,na.rm=T)
# tapply(enth_tab4$resting_heartrate,res$state,mean,na.rm=T)
# tapply(enth_tab4$skin_temp,res$state,mean,na.rm=T)
# tapply(enth_tab4$nb_temp,res$state,mean,na.rm=T)
# 
# # Try categorical
# tapply(enth_tab4$air_vel,res$state,mean,na.rm=T)

enth_tab4$air_vel=round(enth_tab4$air_vel)
enth_tab4$air_vel=as.factor(enth_tab4$air_vel) # 10=slightly perceived, 11=perceived
enth_tab4$clothing=floor(enth_tab4$clothing)
enth_tab4$clothing=as.factor(enth_tab4$clothing)

# Ground truth "thermal"
gt_thermal=enth_tab4$thermal
enth_tab4=select(enth_tab4,subset=-c(thermal,comfort,
                                     pm1.0_outdoor,pm10.0_outdoor,
                                     pm2.5_outdoor))

str(enth_tab4)

enth_tab5=enth_tab4
enth_tab5[,-1]=enth_tab5[,-1]%>% mutate_if(is.numeric,~ifelse(is.nan(.), NA, .))

# Missing %
apply(enth_tab5,2,function(x) sum(is.na(x))/length(x))*100

# Data augmentation -------------------------------------------------------

env_vars=enth_tab5[,c("time","temp","humidity","pressure")]
tmp=enth_tab5[,c("time","temp","humidity","pressure")]

dt=as.numeric(diff(enth_tab5$time))/min(as.numeric(diff(enth_tab5$time)))
env_vars$dt=c(NA,dt)

str(env_vars)

for(i in 2:(nrow(env_vars)-1)){
  if(env_vars$dt[i]==1|env_vars$dt[i+1]==1){
    if(is.na(env_vars$temp[i])){
      env_vars$temp[i]=mean(env_vars$temp[i-1],env_vars$temp[i+1],na.rm=TRUE,trim = .5)
    }
    if(is.na(env_vars$humidity[i])){
      env_vars$humidity[i]=mean(env_vars$humidity[i-1],env_vars$humidity[i+1],na.rm=TRUE,trim = .5)
    }
    if(is.na(env_vars$pressure[i])){
      env_vars$pressure[i]=mean(env_vars$pressure[i-1],env_vars$pressure[i+1],na.rm=TRUE,trim = .5)
    }
  }
}

# % of missing values
apply(env_vars,2,function(x) sum(is.na(x))/length(x))*100

# compare with pre data augmentation
apply(tmp,2,
      function(x) sum(is.na(x))/length(x))*100

enth_tab5[,c("temp","humidity","pressure")]=env_vars[,c("temp","humidity","pressure")]

# Interaction terms
# enth_tab5$temp_hum=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$humidity))
# enth_tab5$temp_press=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$pressure))
# enth_tab5$hum_press=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$pressure))
# enth_tab5$temp_heart=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$heartrate))
# enth_tab5$hum_heart=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$heartrate))
# enth_tab5$press_heart=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$heartrate))
# enth_tab5$temp_resthr=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$resting_heartrate))
# enth_tab5$hum_resthr=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$resting_heartrate))
# enth_tab5$press_resthr=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$resting_heartrate))
# enth_tab5$temp_nbtemp=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$nb_temp))
# enth_tab5$hum_nbtemp=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$nb_temp))
# enth_tab5$press_nbtemp=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$nb_temp))
# enth_tab5$temp_skintemp=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$skin_temp))
# enth_tab5$hum_skintemp=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$skin_temp))
# enth_tab5$press_skintemp=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$skin_temp))

# Correlations
library(zoo)
wdn=10
enth_tab5$corr_temp_humidity=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("temp","humidity")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_temp_pressure=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("temp","pressure")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_humidity_pressure=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("humidity","pressure")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_temp_heart=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("temp","heartrate")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_humidity_heart=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("humidity","heartrate")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_pressure_heart=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("pressure","heartrate")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_temp_resthr=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("temp","resting_heartrate")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_humidity_resthr=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("humidity","resting_heartrate")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_pressure_resthr=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("pressure","resting_heartrate")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_temp_nbtemp=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("temp","nb_temp")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_humidity_nbtemp=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("humidity","nb_temp")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_pressure_nbtemp=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("pressure","nb_temp")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_temp_skintemp=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("temp","skin_temp")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_humidity_skintemp=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("humidity","skin_temp")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))
enth_tab5$corr_pressure_skintemp=c(rep(NA,wdn-1),rollapply(apply(enth_tab5[,c("pressure","skin_temp")],2,scale),
                                                         width=wdn, function(x) cor(x[,1],x[,2],use="complete.obs"), 
                                                         by.column=FALSE))

# Interaction terms
# Interaction terms
enth_tab5$temp_hum=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$humidity))
enth_tab5$temp_press=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$pressure))
enth_tab5$hum_press=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$pressure))
enth_tab5$temp_heart=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$heartrate))
enth_tab5$hum_heart=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$heartrate))
enth_tab5$press_heart=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$heartrate))
enth_tab5$temp_resthr=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$resting_heartrate))
enth_tab5$hum_resthr=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$resting_heartrate))
enth_tab5$press_resthr=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$resting_heartrate))
enth_tab5$temp_nbtemp=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$nb_temp))
enth_tab5$hum_nbtemp=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$nb_temp))
enth_tab5$press_nbtemp=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$nb_temp))
enth_tab5$temp_skintemp=as.numeric(scale(enth_tab5$temp))*as.numeric(scale(enth_tab5$skin_temp))
enth_tab5$hum_skintemp=as.numeric(scale(enth_tab5$humidity))*as.numeric(scale(enth_tab5$skin_temp))
enth_tab5$press_skintemp=as.numeric(scale(enth_tab5$pressure))*as.numeric(scale(enth_tab5$skin_temp))

# replace NaN with NA
enth_tab5[,-1]=enth_tab5[,-1]%>% mutate_if(is.numeric,~ifelse(is.nan(.), NA, .))

# Save 

save(enth_tab5,gt_thermal,file="enth_tab5.RData")
