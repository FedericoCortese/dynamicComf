load("cleaned_data.RData")
source("Utils.R")

# Summarize on hourly basis -----------------------------------------------

summarize_weather=function(data,wdn){
  data2=data %>% 
    group_by(time=floor_date(time,wdn))%>%
    summarise_if(is.numeric, mean, na.rm = TRUE) 
  data2[,-1]=data2[,-1]%>% mutate_all(~ifelse(is.nan(.), NA, .))
  eqdist_times=seq(from=data2$time[1],
                   to=tail(data2$time,1),
                   by=wdn)
  data3=data2%>% 
    right_join(data.frame(time=eqdist_times),by="time")%>%
    arrange(time)
  return(data3)
}

wdn="1 hour"

# Apply to all weather variables
air_temp=summarize_weather(air_temp,wdn)
rainfall=summarize_weather(rainfall,wdn)
RH=summarize_weather(RH,wdn)
wsp=summarize_weather(wind_speed,wdn)
wdir=summarize_weather(wind_dir,wdn)

# Air temp
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(air_temp)){
  plot(x=air_temp$time,y=as.vector(unlist(air_temp[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(air_temp[,i]))
  title(main=colnames(air_temp)[i])
}
mtext("Air temperatures", side = 3, line = - 2, outer = TRUE)

# Relative humidity
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(RH)){
  plot(x=RH$time,y=as.vector(unlist(RH[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(RH[,i]))
  title(main=colnames(RH)[i])
}
mtext("Relative humidity", side = 3, line = - 2, outer = TRUE)


# Rainfall
windows()
par(mfrow=c(4,9),mar=c(2,2,6,2))
for(i in 2:36){
  plot(x=rainfall$time,y=as.vector(unlist(rainfall[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(rainfall[,i]))
  title(main=colnames(rainfall)[i])
}
mtext("Rainfall (1)", side = 3, line = - 2, outer = TRUE)

# Plot for the next 36 stations
windows()
par(mfrow=c(4,9),mar=c(2,2,6,2))
for(i in 37:72){
  plot(x=rainfall$time,y=as.vector(unlist(rainfall[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(rainfall[,i]))
  title(main=colnames(rainfall)[i])
}
mtext("Rainfall (2)", side = 3, line = - 2, outer = TRUE)

# Wind speed
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(wsp)){
  plot(x=wsp$time,y=as.vector(unlist(wsp[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(wsp[,i]))
  title(main=colnames(wsp)[i])
}
mtext("Wind speed", side = 3, line = - 2, outer = TRUE)

# Wind direction
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(wdir)){
  plot(x=wdir$time,y=as.vector(unlist(wdir[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(wdir[,i]))
  title(main=colnames(wdir)[i])
}
mtext("Wind direction", side = 3, line = - 2, outer = TRUE)

#  NA counts --------------------------------------------------------------

na_count=function(x){
  x %>% summarise_all(funs(sum(is.na(.))))
}

# Drop stations with more than 40%  of missings
# Also, impute as NA values that have been linearly interpolated

air_temp_na=na_count(air_temp); air_temp_na=air_temp_na/dim(air_temp)[1]*100;air_temp_na
air_temp=air_temp[,-which(air_temp_na>40)]

RH_na=na_count(RH); RH_na=RH_na/dim(RH)[1]*100
RH=RH[,-which(RH_na>40)]

rain_na=na_count(rainfall); rain_na=rain_na/dim(rainfall)[1]*100
rainfall=rainfall[,-which(rain_na>40)]

wsp_na=na_count(wind_speed); wsp_na=wsp_na/dim(wind_speed)[1]*100
wsp=wsp[,-which(wsp_na>40)]

wdir_na=na_count(wind_dir); wdir_na=wdir_na/dim(wind_dir)[1]*100
wdir=wdir[,-which(wdir_na>40)]

# New plot
# Air temp
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(air_temp)){
  plot(x=air_temp$time,y=as.vector(unlist(air_temp[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(air_temp[,i]))
  title(main=colnames(air_temp)[i])
}
mtext("Air temperatures", side = 3, line = - 2, outer = TRUE)

# Relative humidity
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(RH)){
  plot(x=RH$time,y=as.vector(unlist(RH[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(RH[,i]))
  title(main=colnames(RH)[i])
}
mtext("Relative humidity", side = 3, line = - 2, outer = TRUE)


# Rainfall
windows()
par(mfrow=c(4,9),mar=c(2,2,6,2))
for(i in 2:36){
  plot(x=rainfall$time,y=as.vector(unlist(rainfall[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(rainfall[,i]))
  title(main=colnames(rainfall)[i])
}
mtext("Rainfall (1)", side = 3, line = - 2, outer = TRUE)

# Plot for the next 36 stations
windows()
par(mfrow=c(4,9),mar=c(2,2,6,2))
for(i in 37:72){
  plot(x=rainfall$time,y=as.vector(unlist(rainfall[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(rainfall[,i]))
  title(main=colnames(rainfall)[i])
}
mtext("Rainfall (2)", side = 3, line = - 2, outer = TRUE)


# Wind speed
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(wsp)){
  plot(x=wsp$time,y=as.vector(unlist(wsp[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(wsp[,i]))
  title(main=colnames(wsp)[i])
}
mtext("Wind speed", side = 3, line = - 2, outer = TRUE)

# Wind direction
windows()
par(mfrow=c(4,5),mar=c(2,2,6,2))
for(i in 2:ncol(wdir)){
  plot(x=wdir$time,y=as.vector(unlist(wdir[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(wdir[,i]))
  title(main=colnames(wdir)[i])
}
mtext("Wind direction", side = 3, line = - 2, outer = TRUE)

count_gaps=function(x){
  gaps=list()
  gaps[[1]]=x[,1]
  for(i in 2:ncol(x)){
    d_na <- as.numeric(is.na(x[,i]))
    Nal=rle(d_na)
    tosave=Nal$values==1
    Nals=Nal$lengths[tosave]
    gaps[[i]]=Nals
  }
  names(gaps)=colnames(x)
  return(gaps)
}

# Gap lengths for air_data_wide as % of total length
NA_count_air=count_gaps(air_temp);NA_count_air
#lapply(NA_count[-1],function(x) round(x/dim(air_data_wide)[1]*100,1))



# Spatially interpolate missing values ------------------------------------

locations2=locations[,c("id","longitude","latitude")]

air_temp1=weightdist_imp(x_data=air_temp,locations2=locations2)
RH1=weightdist_imp(x_data=RH,locations2=locations2)
rainfall1=weightdist_imp(x_data=rainfall,locations2=locations2)
wsp1=weightdist_imp(x_data=wsp,locations2=locations2)
wdir1=weightdist_imp(x_data=wdir,locations2=locations2)

NA_count_air=count_gaps(air_temp1);NA_count_air


# New plot
# Air temp
windows()
par(mfrow=c(4,4),mar=c(2,2,6,2))
for(i in 2:ncol(air_temp1)){
  plot(x=air_temp1$time,y=as.vector(unlist(air_temp1[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(air_temp1[,i]))
  title(main=colnames(air_temp1)[i])
}
mtext("Air temperatures", side = 3, line = - 2, outer = TRUE)

# Relative humidity
windows()
par(mfrow=c(4,4),mar=c(2,2,6,2))
for(i in 2:ncol(RH1)){
  plot(x=RH1$time,y=as.vector(unlist(RH1[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(RH1[,i]))
  title(main=colnames(RH1)[i])
}
mtext("Relative humidity", side = 3, line = - 2, outer = TRUE)


# Rainfall
windows()
par(mfrow=c(4,9),mar=c(2,2,6,2))
for(i in 2:36){
  plot(x=rainfall1$time,y=as.vector(unlist(rainfall1[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(rainfall1[,i]))
  title(main=colnames(rainfall1)[i])
}
mtext("Rainfall (1)", side = 3, line = - 2, outer = TRUE)

# Plot for the next 36 stations
windows()
par(mfrow=c(4,9),mar=c(2,2,6,2))
for(i in 37:68){
  plot(x=rainfall1$time,y=as.vector(unlist(rainfall1[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(rainfall1[,i]))
  title(main=colnames(rainfall1)[i])
}
mtext("Rainfall (2)", side = 3, line = - 2, outer = TRUE)


# Wind speed
windows()
par(mfrow=c(4,4),mar=c(2,2,6,2))
for(i in 2:ncol(wsp)){
  plot(x=wsp1$time,y=as.vector(unlist(wsp1[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(wsp1[,i]))
  title(main=colnames(wsp1)[i])
}
mtext("Wind speed", side = 3, line = - 2, outer = TRUE)

# Wind direction
windows()
par(mfrow=c(4,4),mar=c(2,2,6,2))
for(i in 2:ncol(wdir1)){
  plot(x=wdir1$time,y=as.vector(unlist(wdir1[,i])),type="l",col="blue",
       xlab=" ",ylab=" ",
       main=colnames(wdir1[,i]))
  title(main=colnames(wdir1)[i])
}
mtext("Wind direction", side = 3, line = - 2, outer = TRUE)


# Retrieve relevant times, long and lat to be predicted -------------------

load("df_cozie.Rdata")
data_input=df_cozie[,c("time","ws_longitude","ws_latitude")]
colnames(data_input)=c("time","longitude","latitude")

# Naive prediction --------------------------------------------------------

# Air temp
air_temp_pred=weightdist_pred(x_input=data_input,locations2=locations2,
                              x_data=air_temp1,name="air_temp")

# Relative humidity
RH_pred=weightdist_pred(x_input=data_input,locations2=locations2,
                        x_data=RH1,name="RH")

# Rainfall
rainfall_pred=weightdist_pred(x_input=data_input,locations2=locations2,
                              x_data=rainfall1,name="rainfall")

# Wind speed
wsp_pred=weightdist_pred(x_input=data_input,locations2=locations2,
                         x_data=wsp1,name="wsp")

# Wind direction
wdir_pred=weightdist_pred(x_input=data_input,locations2=locations2,
                          x_data=wdir1,name="wdir")

# # Remove NAs
# air_temp_pred2=air_temp_pred[complete.cases(air_temp_pred),]
# RH_pred2=RH_pred[complete.cases(RH_pred),]
# rainfall_pred2=rainfall_pred[complete.cases(rainfall_pred),]
# wsp_pred2=wsp_pred[complete.cases(wsp_pred),]
# wdir_pred2=wdir_pred[complete.cases(wdir_pred),]

# Merge by time, longitude, latiitude
data_pred=merge(air_temp_pred,RH_pred,by=c("time","longitude","latitude"))
data_pred=merge(data_pred,rainfall_pred,by=c("time","longitude","latitude"))
data_pred=merge(data_pred,wsp_pred,by=c("time","longitude","latitude"))
data_pred=merge(data_pred,wdir_pred,by=c("time","longitude","latitude"))

# Remove NAs
data_pred=data_pred[complete.cases(data_pred),]

#id_input=df_cozie[,c("time","ws_longitude","ws_latitude","id_participant")]
id_input=subset(df_cozie,select=-c(ws_survey_count,
                                   ws_timestamp_start,
                                   ws_timestamp_location,
                                   dT,
                                   id_unique))
changename=which(colnames(id_input)=="ws_longitude"|colnames(id_input)=="ws_latitude")
colnames(id_input)[changename]=c("latitude","longitude")

data_pred=merge(data_pred,id_input,by=c("time","longitude","latitude"))
table(data_pred$id)


# Only outdoor ------------------------------------------------------------

data_model=data_pred[data_pred$q_location=="Outdoor",]
data_model=subset(data_model,select=-c(q_location,q_location_office,
                                       q_location_transport))

# Prop odds ratio --------------------------------------------------------------------
data_model$q_thermal_preference=ordered(data_model$q_thermal_preference,
                                 levels=c("Warmer","No change", "Cooler"))
colnames(data_model)
dat=data_model[,-c(1:3,9)];colnames(dat)

summary(dat$q_thermal_preference)
barplot(table(dat$q_thermal_preference)/length(dat$q_thermal_preference),
        col=rainbow(length(unique(dat$q_thermal_preference))))

summary(dat[c("air_temp","rainfall","RH","wsp",
              "Green.View.Mean","Sky.View.Mean","Building.View.Mean")])

library(doBy)
#Air Temperature
summaryBy(air_temp ~ q_thermal_preference, data=dat, FUN=c(mean, sd, min, max))

#Rainfall
summaryBy(rainfall ~ q_thermal_preference, data=dat, FUN=c(mean, sd, min, max))

#Relative humidity
summaryBy(RH ~ q_thermal_preference, data=dat, FUN=c(mean, sd, min, max))

#Wind speed
summaryBy(wsp ~ q_thermal_preference, data=dat, FUN=c(mean, sd, min, max))

# Add combined effect of air temp and humidity
dat2=data.frame(dat,
                air_RH=dat$air_temp*dat$RH)
# dat_scaled=apply(dat2[,-c(11,16)],2,scale)
# dat_scaled=data.frame(dat_scaled,
#                       q_noise_nearby=dat$q_noise_nearby,
#                       q_thermal_preference=dat$q_thermal_preference)

# Scale all numeric variables with dplyr
dat2=dat2 %>% 
  mutate_if(is.numeric, scale)

fit_sc=MASS::polr(q_thermal_preference~air_temp+RH+air_RH,data=dat2)
summary(fit_sc)

fit_=MASS::polr(q_thermal_preference~.,data=dat)
summary(fit_)

# Wide format - SJM - dynamic map of comfort ------------------------------------



# LMM ---------------------------------------------------------------------
library(LMest)
?est_lm_cov_manifest

