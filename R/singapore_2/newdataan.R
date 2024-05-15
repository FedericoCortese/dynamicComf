
# Data loading ------------------------------------------------------------
enth_surv=read.csv("enth_surveys_calc.csv" )
enth_tab=read.csv("enth_tabular_merged.csv")

# Data cleaning -----------------------------------------------------------
str(enth_tab)
enth_tab$air_vel=as.factor(enth_tab$air_vel)
enth_tab$air_vel=recode(enth_tab$air_vel, "11"= "Perceived", "10"= "Not Perceived")

enth_tab$body_presence=as.factor(enth_tab$body_presence)
enth_tab$change=as.factor(enth_tab$change)

enth_tab$clothing=as.ordered(enth_tab$clothing)
enth_tab$clothing=recode(enth_tab$clothing, "8"="Very Light","9" = "Light", 
                         "10"= "Medium", "11"= "Heavy")

enth_tab$comfort=as.ordered(enth_tab$comfort)
enth_tab$comfort=recode(enth_tab$comfort, "9" = "Not Comfy", "10"= "Comfy")

enth_tab$indoor.outdoor=as.factor(enth_tab$indoor.outdoor)

enth_tab$thermal=as.ordered(enth_tab$thermal)
enth_tab$thermal=recode(enth_tab$thermal, "9" = "Warmer", "10"= "No Change", "11"= "Cooler")

enth_tab$body_presence=as.factor(enth_tab$body_presence)
enth_tab$user_id=as.factor(enth_tab$user_id)

enth_tab$met=as.factor(enth_tab$met)
enth_tab$met=recode(enth_tab$met, "8" = "Resting", "9"= "Sitting", "10"= "Standing", "11"= "Exercising")


str(enth_surv)
enth_surv$sex=as.factor(enth_surv$sex)
enth_surv$used_weather=as.factor(enth_surv$used_weather)
enth_surv$satisfaction_weather=as.factor(enth_surv$satisfaction_weather)
enth_surv$sweating=as.ordered(enth_surv$sweating)
enth_surv$enjoy_ourdoor=as.factor(enth_surv$enjoy_ourdoor)

library(Amelia)
missmap(enth_tab, main = "Missing values vs observed",margins = c(10,2))

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
                     "X0.3um_count_outdoor",
                     "X0.5um_count_outdoor",
                     "X1.0um_count_outdoor",
                     "X10.0um_count_outdoor",
                     "X2.5um_count_outdoor",
                     "X5.0um_count_outdoor",
                     "humidity_outdoor",
                     "pm1.0_outdoor",
                     "pm10.0_outdoor",
                     "pm2.5_outdoor",
                     "pressure_outdoor",
                     "temp_outdoor")]

summary(enth_out$thermal)

library(Amelia)
missmap(enth_out, main = "Missing values vs observed",margins = c(10,2))
missmap(enth_surv, main = "Missing values vs observed",margins = c(10,2))

# Remove NAs
enth_out_complete=enth_out[complete.cases(enth_out),]

# Merge with survey data
dat=merge(enth_out_complete, enth_surv, by="user_id")


# Exploratory Data Analysis -----------------------------------------------
library(doBy)

# Impact of weather variables on thermal and general comfort
summaryBy(temp_outdoor ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(temp_outdoor ~ comfort, data=dat, FUN=c(mean, sd, min, max))
summaryBy(humidity_outdoor ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(humidity_outdoor ~ comfort, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pressure_outdoor ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pressure_outdoor ~ comfort, data=dat, FUN=c(mean, sd, min, max))

con.comf.airvel<-table(dat$air_vel,dat$comfort)
con.comf.airvel/rowSums(con.comf.airvel)*100
con.therm.airvel<-table(dat$air_vel,dat$thermal)
con.therm.airvel/rowSums(con.therm.airvel)*100

# Impact of air quality on thermal and general comfort
summaryBy(pm1.0_outdoor ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pm1.0_outdoor ~ comfort, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pm2.5_outdoor ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pm2.5_outdoor ~ comfort, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pm10.0_outdoor ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(pm10.0_outdoor ~ comfort, data=dat, FUN=c(mean, sd, min, max))

# Impact of physiological variables on thermal and general comfort
con.therm.met<-table(dat$met,dat$thermal)
con.therm.met[2:3,2:3]/rowSums(con.therm.met[2:3,2:3])*100

con.comf.met<-table(dat$met,dat$comfort)
con.comf.met[2:3,]/rowSums(con.comf.met[2:3,])*100

summary(dat$clothing)
con.comf.clot<-table(dat$clothing,dat$comfort)
con.comf.clot/rowSums(con.comf.clot)*100
con.comf.therm<-table(dat$clothing,dat$thermal)
con.comf.therm[,2:3]/rowSums(con.comf.therm[,2:3])*100

summaryBy(heartrate ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(heartrate ~ comfort, data=dat, FUN=c(mean, sd, min, max))

summaryBy(resting_heartrate ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(resting_heartrate ~ comfort, data=dat, FUN=c(mean, sd, min, max))

summaryBy(nb_temp ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(nb_temp ~ comfort, data=dat, FUN=c(mean, sd, min, max))

summaryBy(skin_temp ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(skin_temp ~ comfort, data=dat, FUN=c(mean, sd, min, max))

# Impact of demographic variables on thermal and general comfort
dat$BMI=dat$weight/(dat$height/100)^2
summaryBy(BMI ~ thermal, data=dat, FUN=c(median, sd, min, max))
summaryBy(BMI ~ comfort, data=dat, FUN=c(median, sd, min, max))

summaryBy(years_here ~ thermal, data=dat, FUN=c(mean, sd, min, max))

con.comf.sweat<-table(dat$sweating,dat$comfort)
con.comf.sweat/rowSums(con.comf.sweat)*100
con.therm.sweat<-table(dat$sweating,dat$thermal)
con.therm.sweat/rowSums(con.therm.sweat)*100

summaryBy(outdoor_hr_weekday ~ thermal, data=dat, FUN=c(mean, sd, min, max))
summaryBy(outdoor_hr_weekday ~ comfort, data=dat, FUN=c(mean, sd, min, max))

con.comf.enjoy<-table(dat$enjoy_ourdoor,dat$comfort)
con.comf.enjoy/rowSums(con.comf.enjoy)*100
con.therm.enjoy<-table(dat$enjoy_ourdoor,dat$thermal)
con.therm.enjoy/rowSums(con.therm.enjoy)*100

# Prop Odds Ratio  --------------------------------------------------------

model <- glm(thermal ~ 
               temp_outdoor+skin_temp+
               BMI+humidity_outdoor+
               clothing+met+
               enjoy_ourdoor,family=binomial(link='logit'),data=dat)
summary(model)
