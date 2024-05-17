library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

# Data loading ------------------------------------------------------------
enth_surv=read.csv("enth_surveys_calc.csv" )
enth_tab=read.csv("enth_tabular_merged.csv")
# enth_tab=subset(enth_tab, select=-c(response_speed,space_id,building_name,body_presence,change))

# Data cleaning -----------------------------------------------------------
str(enth_tab)

enth_tab$time=as.POSIXct(enth_tab$time, format="%Y-%m-%d %H:%M:%S")

# Order by user_id
enth_tab=enth_tab[order(enth_tab$user_id, enth_tab$time),]

# # Replace all NA with 999
# enth_tab[is.na(enth_tab)] <- 999
# 
# enth_tab2=unique( enth_tab )
# 
# which((diff(enth_tab$time)==0))
# 
# enth_tab2[140,]==enth_tab2[141,]
# 
# enth_tab[139:141,c("time","user_id")]
# enth_tab[207:209,]

# Round time to hour

enth_tab$air_vel=as.factor(enth_tab$air_vel)
# "11"= "Perceived" = 1, "10"= "Not Perceived" = 0
enth_tab$air_vel=recode(enth_tab$air_vel, "11"= 1, "10"= 0)

enth_tab$body_presence=as.factor(enth_tab$body_presence)
enth_tab$change=as.factor(enth_tab$change)

enth_tab$clothing=as.ordered(enth_tab$clothing)
# "8"="Very Light" = 0,"9" = "Light"=1,  "10"= "Medium"=2, "11"= "Heavy"=3
enth_tab$clothing=recode(enth_tab$clothing, "8"=0,"9" = 1, 
                         "10"= 2, "11"= 3)

enth_tab$comfort=as.ordered(enth_tab$comfort)
# 9=  "Not Comfy" =1, "10"= "Comfy"=0
enth_tab$comfort=recode(enth_tab$comfort, "9" = 1, "10"= 0)

enth_tab$indoor.outdoor=as.factor(enth_tab$indoor.outdoor)

enth_tab$thermal=as.ordered(enth_tab$thermal)
# "9" = "Warmer" = -1, "10"= "No Change" = 0, "11"= "Cooler" =1
enth_tab$thermal=recode(enth_tab$thermal, "9" = -1, "10"= 0, "11"= 1)

enth_tab$body_presence=as.factor(enth_tab$body_presence)
enth_tab$user_id=as.factor(enth_tab$user_id)

enth_tab$met=as.factor(enth_tab$met)
# "8" = "Resting" = 0, "9"= "Sitting" =1, "10"= "Standing" =2, "11"= "Exercising" =3
enth_tab$met=recode(enth_tab$met, "8" = 0, "9"= 1, "10"= 2, "11"= 3)


str(enth_surv)
#enth_surv$sex=as.factor(enth_surv$sex)
enth_surv$sex=recode(enth_surv$sex, "Female" = 0, "Male"= 1)

enth_surv$used_weather=as.factor(enth_surv$used_weather)

#enth_surv$satisfaction_weather=as.factor(enth_surv$satisfaction_weather)
enth_surv$satisfaction_weather=recode(enth_surv$satisfaction_weather, "Extremely Dissatisfied"=-3,
                                      "Moderately Dissatisfied" =-2  ,        
                                      "Moderately Satisfied"=2, "Neither Satisfied nor Dissatisfied"=0,
                                     "Slightly Dissatisfied"=-1 ,"Slightly Satisfied"=1)

#enth_surv$sweating=as.ordered(enth_surv$sweating)
enth_surv$sweating=recode(enth_surv$sweating, "3"=0, "4"=1, "5"=2)

#enth_surv$enjoy_ourdoor=as.factor(enth_surv$enjoy_ourdoor)
enth_surv$enjoy_ourdoor=recode(enth_surv$enjoy_ourdoor, "Yes"=1,"No"=0)

enth_surv$user_id=as.factor(enth_surv$user_id)


library(Amelia)
missmap(enth_tab, main = "Missing values vs observed",margins = c(10,2))


# Indoor Only -------------------------------------------------------------



# Outdoor only ------------------------------------------------------------
# Extract outdoor observations (9 means outdoor)
enth_out=enth_tab[enth_tab$indoor.outdoor==9,]

# Select only relevant features
enth_out=enth_out[,c("user_id",
                     "time",
                     "air_vel",
                     "clothing",
                     "comfort",
                     "heartrate",
                     "met",
                     "resting_heartrate",
                     "thermal",
                     "nb_temp",
                     "skin_temp",
                     "humidity_outdoor",
                     "pm1.0_outdoor",
                     "pm10.0_outdoor",
                     "pm2.5_outdoor",
                     "pressure_outdoor",
                     "temp_outdoor")]

summary(enth_out$thermal)

#which((diff(enth_out$time)==0))

enth_out=enth_out[-which((diff(enth_out$time)==0)),]

which((diff(enth_out$time)==0))
which(is.na(enth_out$thermal))
enth_out=enth_out[-which(is.na(enth_out$thermal)),]


library(Amelia)
missmap(enth_out, main = "Missing values vs observed",margins = c(10,2))
# missmap(enth_surv, main = "Missing values vs observed",margins = c(10,2))


# Exploratory Data Analysis -----------------------------------------------
library(doBy)

enth_out_EDA=enth_out[complete.cases(enth_out),]

# Impact of weather variables on thermal and general comfort
summaryBy(temp_outdoor ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(temp_outdoor ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(humidity_outdoor ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(humidity_outdoor ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pressure_outdoor ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pressure_outdoor ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)

con.comf.airvel<-table(enth_out_EDA$air_vel,enth_out_EDA$comfort)
con.comf.airvel/rowSums(con.comf.airvel)*100
con.therm.airvel<-table(enth_out_EDA$air_vel,enth_out_EDA$thermal)
con.therm.airvel/rowSums(con.therm.airvel)*100

# Impact of air quality on thermal and general comfort
summaryBy(pm1.0_outdoor ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pm1.0_outdoor ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pm2.5_outdoor ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pm2.5_outdoor ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pm10.0_outdoor ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(pm10.0_outdoor ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)

# Impact of physiological variables on thermal and general comfort
con.therm.met<-table(enth_out_EDA$met,enth_out_EDA$thermal)
con.therm.met/rowSums(con.therm.met)*100

con.comf.met<-table(enth_out_EDA$met,enth_out_EDA$comfort)
con.comf.met/rowSums(con.comf.met)*100

summary(dat$clothing)
con.comf.clot<-table(enth_out_EDA$clothing,enth_out_EDA$comfort)
con.comf.clot/rowSums(con.comf.clot)*100
con.comf.therm<-table(enth_out_EDA$clothing,enth_out_EDA$thermal)
con.comf.therm/rowSums(con.comf.therm)*100

summaryBy(heartrate ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(heartrate ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)

summaryBy(resting_heartrate ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(resting_heartrate ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)

summaryBy(nb_temp ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(nb_temp ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)

summaryBy(skin_temp ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
summaryBy(skin_temp ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)

# Impact of demographic variables on thermal and general comfort
# enth_out_EDA$BMI=enth_out_EDA$weight/(enth_out_EDA$height/100)^2
# summaryBy(BMI ~ thermal, data=enth_out_EDA, FUN=c(median, sd, min, max),na.rm=T)
# summaryBy(BMI ~ comfort, data=enth_out_EDA, FUN=c(median, sd, min, max),na.rm=T)
# 
# summaryBy(years_here ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
# 
# con.comf.sweat<-table(dat$sweating,dat$comfort)
# con.comf.sweat/rowSums(con.comf.sweat)*100
# con.therm.sweat<-table(dat$sweating,dat$thermal)
# con.therm.sweat/rowSums(con.therm.sweat)*100
# 
# summaryBy(outdoor_hr_weekday ~ thermal, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
# summaryBy(outdoor_hr_weekday ~ comfort, data=enth_out_EDA, FUN=c(mean, sd, min, max),na.rm=T)
# 
# con.comf.enjoy<-table(dat$enjoy_ourdoor,dat$comfort)
# con.comf.enjoy/rowSums(con.comf.enjoy)*100
# con.therm.enjoy<-table(dat$enjoy_ourdoor,dat$thermal)
# con.therm.enjoy/rowSums(con.therm.enjoy)*100

# Prop Odds Ratio  --------------------------------------------------------

# model <- glm(thermal ~ 
#                temp_outdoor+skin_temp+
#                BMI+humidity_outdoor+
#                clothing+met+
#                enjoy_ourdoor,family=binomial(link='logit'),data=dat)
# summary(model)


# # HMM for continuous vars with miss  --------------------------------------
# 
# require(LMest)
# require(mvtnorm)
# require(MultiLCIRT)
# 
# source('HMMcont_miss/RCode/lmbasic.cont.MISS.R')
# source('HMMcont_miss/RCode/lk_comp_cont_MISS.R')
# source('HMMcont_miss/RCode/lmcovlatent.cont.MISS.R')
# source('HMMcont_miss/RCode/lk_comp_latent_cont_MISS.R')
# source('HMMcont_miss/RCode/prob_post_cov_cont.R')
# source('HMMcont_miss/RCode/est_multilogit.R')
# source('HMMcont_miss/RCode/prob_multilogit.R')
# source("HMMcont_miss/RCode/bootstrap.MISS.R")
# 
# # Construct objects
# ui=unique(enth_out$user_id)
# ut=sort(unique(enth_out$time))
# 
# n=length(unique(ui))
# TT=length(unique(ut))
# 
# temp=data.frame(time=rep(ut,each=n),user_id=rep(unique(ui),length(ut)))
# 
# # Join temp and enth_out
# temp=merge(temp,enth_out,by=c("time","user_id"),all=T)
# 
# # Construct response vars object
# Yweather=temp[,c("user_id","time","temp_outdoor","humidity_outdoor","pressure_outdoor")]
# Yphysio=temp[,c("user_id","time","heartrate","resting_heartrate","nb_temp","skin_temp")]
# Yaq=temp[,c("user_id","time","pm1.0_outdoor","pm2.5_outdoor","pm10.0_outdoor")]
# 
# YYweather=array(0,dim=c(n,TT,3))
# YYphysio=array(0,dim=c(n,TT,4))
# YYaq=array(0,dim=c(n,TT,3))
# 
# for(i in 1:n){
#   YYweather[i,,1]=temp[temp$user_id==ui[i],"temp_outdoor"]
#   YYweather[i,,2]=temp[temp$user_id==ui[i],"humidity_outdoor"]
#   YYweather[i,,3]=temp[temp$user_id==ui[i],"pressure_outdoor"]
#   YYphysio[i,,1]=temp[temp$user_id==ui[i],"heartrate"]
#   YYphysio[i,,2]=temp[temp$user_id==ui[i],"resting_heartrate"]
#   YYphysio[i,,3]=temp[temp$user_id==ui[i],"nb_temp"]
#   YYphysio[i,,4]=temp[temp$user_id==ui[i],"skin_temp"]
#   YYaq[i,,1]=temp[temp$user_id==ui[i],"pm1.0_outdoor"]
#   YYaq[i,,2]=temp[temp$user_id==ui[i],"pm2.5_outdoor"]
#   YYaq[i,,3]=temp[temp$user_id==ui[i],"pm10.0_outdoor"]
# }

#YY1=log(YY1)


# LMM with covariates on the latent model ---------------------------------------------------------------------

enth_surv$BMI=enth_surv$weight/(enth_surv$height/100)^2
enth_surv$outdoor_hr_day=apply(cbind(enth_surv$outdoor_hr_weekday,enth_surv$outdoor_hr_weekend),1,mean)

XX1=array(0,dim=c(n,TT,15))
for(i in 1:n){
  XX1[i,,1]=enth_surv$sex[which(enth_surv$user_id==ui[i])]
  XX1[i,,2]=enth_surv$sweating[which(enth_surv$user_id==ui[i])]
  XX1[i,,3]=enth_surv$BMI[which(enth_surv$user_id==ui[i])]
  XX1[i,,4]=enth_surv$enjoy_ourdoor[which(enth_surv$user_id==ui[i])]
  XX1[i,,5]=enth_surv$satisfaction_weather[which(enth_surv$user_id==ui[i])]
  XX1[i,,6]=enth_surv$years_here[which(enth_surv$user_id==ui[i])]
  XX1[i,,7]=enth_surv$outdoor_hr_day[which(enth_surv$user_id==ui[i])]
  XX1[i,,8]=enth_surv$shoulder_circumference[which(enth_surv$user_id==ui[i])]
  XX1[i,,9]=enth_surv$hsps[which(enth_surv$user_id==ui[i])]
  XX1[i,,10]=enth_surv$swls[which(enth_surv$user_id==ui[i])]
  XX1[i,,11]=enth_surv$extraversion[which(enth_surv$user_id==ui[i])]
  XX1[i,,12]=enth_surv$agreeableness[which(enth_surv$user_id==ui[i])]
  XX1[i,,13]=enth_surv$conscientiousness[which(enth_surv$user_id==ui[i])]
  XX1[i,,14]=enth_surv$emotional_stability[which(enth_surv$user_id==ui[i])]
  XX1[i,,15]=enth_surv$openness_to_experiences[which(enth_surv$user_id==ui[i])]
}

prv=matrices2long(YY1[,,8],XX1)
colnames(prv)[3:18]=c("thermal","sex","sweating","BMI","enjoy_ourdoor",
                      "satisfaction_weather","years_here","outdoor_hr_day",
                      "shoulder_circumference","hsps","swls","extraversion",
                      "agreeableness","conscientiousness","emotional_stability",
                      "openness_to_experiences")

prv$thermal=1+prv$thermal

out <- lmest(responsesFormula = thermal ~ NULL,
             latentFormula = ~ sex+BMI+enjoy_ourdoor+
               satisfaction_weather+years_here+outdoor_hr_day+
               shoulder_circumference+hsps+swls+extraversion+
               agreeableness+conscientiousness+emotional_stability+
               openness_to_experiences,
             index = c("id","time"),
             data = prv,
             k = 3,
             start = 1,
             modBasic = 1,
             seed = 123)

summary(out)
out$Psi

outDec <- lmestDecoding(out)

dim(outDec$Ug)

prv$thermal[which(prv$id==4)]
outDec$Ug[4,]




dim(YY1)
n <- dim(YY1)[1]                   # number of individuals
TT <- dim(YY1)[2]                  # number of time occasions
r <- dim(YY1)[3]                   # number of continuous outcomes            
head(YY1[,1,])
dim(XX1)
head(XX1[,1,])

nrep <- 10 
Kmax <- 4
tol <- 10^-4

modva = vector("list",Kmax)
for(k in 1:Kmax){
  print(k)
  if(k>=7) nrep<-2
  modva[[k]] <- lmbasic.cont.MISS(YY1,k=k,modBasic=1,start=0,tol=tol)
  if(k>1){
    for(k1 in 1:(nrep*(k-1))){
      print(c(k,k1))
      tmp <- lmbasic.cont.MISS(YY,k=k,modBasic=1, start=1,tol=tol)
      if(tmp$lk>modva[[k]]$lk){
        modva[[k]] = tmp
      }
    }
  }
}

#round(mp$Mu,3)
exp(round(mp2$Mu,3))
mp2$Si
round(mp2$Pi[,,2],3)
round(mp2$piv,3)

# Too many NAS leads to errors on the HM estimation
# na_ids=apply(YY1,1,function(x) sum(is.na(x)))/(4*TT)
# na_ids=which(na_ids>.99)

# YY1_drop=YY1[-na_ids,,]
# XX1_drop=XX1[-na_ids,,]


# Individual covariates affecting the initial probabilities
X1 <- XX1[,1,]
# Individual covariates affecting the transition probabilities
X2 <- XX1[,2:TT,]
mpcov <- lmcovlatent.cont.MISS(Y=YY1,k=2,X1=X1,
                               X2=X2,start=1,tol=tol,output=TRUE)



#
modva = vector("list",Kmax)
# for(k in 1:Kmax){
#   print(k)
#   if(k>=7) nrep<-2
#   modva[[k]] <- lmbasic.cont.MISS(YY,k=k,modBasic=1,start=0,tol=tol)
#   if(k>1){
    for(k1 in 1:(nrep*(k-1))){
      print(c(k,k1))
      tmp <- lmbasic.cont.MISS(YY1,k=3,modBasic=1, start=1,tol=tol)
      if(tmp$lk>modva[[k]]$lk){
        modva[[k]] = tmp
      }
     }
#   }
# }
