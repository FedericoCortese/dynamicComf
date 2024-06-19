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

wdn="10 mins"
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
cbind(cc$av,enth_tab_av$time)
sort(cc$av)
cc$max.t
cc$t
# Waaaaay better
enth_tab_av$time[(156-56):155]

#wdn2=(cc$t-1-cc$max.t):cc$t
wdn2=(156-56):155

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
source("Utils.R")
lambda=.2
n_states=2
enth_tab5=enth_tab4
enth_tab5=enth_tab5%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))

est=jump_mixed(enth_tab5[,-1], n_states, jump_penalty=lambda, 
               initial_states=NULL,
               max_iter=10, n_init=10, tol=NULL, verbose=T) 

res=data.frame(enth_tab4,state=est$best_s)
tapply(res$temp,res$state,mean,na.rm=T)
tapply(res$humidity,res$state,mean,na.rm=T)
tapply(res$pressure,res$state,mean,na.rm=T)

# Included? Maybe not
table(res$thermal,res$state)

# No comfort, only thermal
table(res$comfort,res$state)

tapply(res$met,res$state,Mode,na.rm=T)
tapply(res$indoor.outdoor,res$state,Mode,na.rm=T)
tapply(res$heartrate,res$state,mean,na.rm=T)
tapply(enth_tab4$resting_heartrate,res$state,mean,na.rm=T)
tapply(enth_tab4$skin_temp,res$state,mean,na.rm=T)
tapply(enth_tab4$nb_temp,res$state,mean,na.rm=T)

# Try categorical
tapply(enth_tab4$air_vel,res$state,mean,na.rm=T)
