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
enth_tab=subset(enth_tab,select=-c(building_name,
                                   response_speed,
                                   indoor_floor
                                   ))

enth_tab$time=as.POSIXct(enth_tab$time,format="%Y-%m-%d %H:%M:%S")

wdn="15 mins"
enth_tab2=enth_tab%>%
  group_by(user_id,time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE) 

enth01_all=enth_tab[which(enth_tab$user_id=="enth01"),]
enth02_all=enth_tab[which(enth_tab$user_id=="enth02"),]
enth03_all=enth_tab[which(enth_tab$user_id=="enth03"),]
enth04_all=enth_tab[which(enth_tab$user_id=="enth04"),]
enth05_all=enth_tab[which(enth_tab$user_id=="enth05"),]
enth07_all=enth_tab[which(enth_tab$user_id=="enth07"),]
enth09_all=enth_tab[which(enth_tab$user_id=="enth09"),]
enth10_all=enth_tab[which(enth_tab$user_id=="enth10"),]
enth11_all=enth_tab[which(enth_tab$user_id=="enth11"),]
enth13_all=enth_tab[which(enth_tab$user_id=="enth13"),]
enth15_all=enth_tab[which(enth_tab$user_id=="enth15"),]
enth16_all=enth_tab[which(enth_tab$user_id=="enth16"),]
enth17_all=enth_tab[which(enth_tab$user_id=="enth17"),]
enth20_all=enth_tab[which(enth_tab$user_id=="enth20"),]
enth22_all=enth_tab[which(enth_tab$user_id=="enth22"),]
enth25_all=enth_tab[which(enth_tab$user_id=="enth25"),]
enth28_all=enth_tab[which(enth_tab$user_id=="enth28"),]

count_consecutive=function(x,tol=4){
  d_t=diff(x$time)/60
  counter=d_t<tol
  av=ave(counter, cumsum(!counter), FUN = cumsum)
  return(max(av))
}

count_consecutive(enth01_all)
count_consecutive(enth02_all)
count_consecutive(enth03_all)
count_consecutive(enth04_all)
count_consecutive(enth05_all)
count_consecutive(enth07_all)
count_consecutive(enth09_all)
count_consecutive(enth10_all)
count_consecutive(enth11_all)
count_consecutive(enth13_all)
count_consecutive(enth15_all)
count_consecutive(enth16_all)
count_consecutive(enth17_all)
count_consecutive(enth20_all)
count_consecutive(enth22_all)
count_consecutive(enth25_all)
count_consecutive(enth28_all)

# "Average" person
enth_tab_av=enth_tab%>%group_by(time)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

count_consecutive(enth_tab_av)
# Waaaaay better
