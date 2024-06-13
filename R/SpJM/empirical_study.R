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

wdn2=(cc$t-1-cc$max.t):cc$t

enth_tab3=enth_tab_av[wdn2,]
Amelia::missmap(enth_tab3)

str(enth_tab3)

# NAs percentage for each variable
apply(enth_tab3,2,function(x) sum(is.na(x))/length(x))*100

# Remove vars with 100% NAs
enth_tab3=subset(enth_tab3,select=-c(co2_indoor,voc_indoor,pm25_indoor,noise_indoor))

enth_tab

plot(enth_tab3$indoor.outdoor)

Amelia::missmap(enth_tab3)
