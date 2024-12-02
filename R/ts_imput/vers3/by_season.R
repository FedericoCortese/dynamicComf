# Set time zone
Sys.setenv(TZ='UTC')

# Only 2023 ---------------------------------------------------------------
load("prov_genova_tempANDrh.RData")
dim(dat_temp_wide3)
summary(dat_temp_wide3)

range(dat_temp_wide3$time)

time_start=which(dat_temp_wide3$time=="2022-12-21 00:00:00")
time_end=which(dat_temp_wide3$time=="2023-12-20 23:00:00")

dat_air_2023=dat_temp_wide3[time_start:time_end,]
dat_rh_2023=dat_rh_wide3[time_start:time_end,]
save(dat_air_2023,dat_rh_2023,locations3,file="dat_2023.RData")


# Load 2023 ---------------------------------------------------------------

load("dat_2023.RData")
source("Utils2.R")

summary(dat_air_2023)

locations3$stat_numb=(1:nrow(locations3))+1

map3=leaflet() %>% 
  addTiles() %>%
  addCircleMarkers(data=locations3[locations3$`NOME STAZIONE` %in% colnames(dat_air_2023),],
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",Latitude, "<br>Lon:", Longitude
                                  , "<br>Height (m):", Height, "<br>Station:", `NOME STAZIONE`,"<br>Number:",stat_numb
                   )
  )
map3

n_stations=ncol(dat_air_2023)-1

# I Winter I --------------------------------------------------------

winter_start=which(dat_air_2023$time=="2022-12-21 00:00:00 UTC")
winter_end=which(dat_air_2023$time=="2023-03-20 23:00:00 UTC")

# Select season
dat_air_winter=dat_air_2023[winter_start:winter_end,]

rh_winter=dat_rh_2023[winter_start:winter_end,]


# Fill with NA ------------------------------------------------------------


TT_winter=dim(dat_air_winter)[1]
na_len=TT_winter*c(.05,.1,.2)
#na_start=rep(TT/2-7,3)
na_start=floor(TT_winter-na_len)

air_5_winter=dat_air_winter
air_10_winter=dat_air_winter
air_20_winter=dat_air_winter

# For cycle, randomly fill with NAs 
for(i in 2:ncol(dat_air_winter)){
  set.seed(i)
  na_st=sample(1:na_start[1],1)
  air_5_winter[na_st:(na_st+na_len[1]),i]=NA
  na_st=sample(1:na_start[2],1)
  air_10_winter[na_st:(na_st+na_len[2]),i]=NA
  na_st=sample(1:na_start[3],1)
  air_20_winter[na_st:(na_st+na_len[3]),i]=NA
}

# windows()
# par(mfrow=c(6,6),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5_winter)){
#   plot(x=air_5_winter$time,y=as.vector(unlist(air_5_winter[,i])),type="l",col="black",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5_winter[,i]))
#  # title(main=colnames(air_5_winter)[i])
# }
# mtext("Air temperatures - 5% NAs", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(6,6),mar=c(2,2,6,2))
# for(i in 2:ncol(air_10_winter)){
#   plot(x=air_10_winter$time,y=as.vector(unlist(air_10_winter[,i])),type="l",col="black",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5_winter[,i]))
#   # title(main=colnames(air_5_winter)[i])
# }
# mtext("Air temperatures - 10% NAs", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(6,6),mar=c(2,2,6,2))
# for(i in 2:ncol(air_20_winter)){
#   plot(x=air_20_winter$time,y=as.vector(unlist(air_20_winter[,i])),type="l",col="black",
#        xlab=" ",ylab=" ",
#        main=colnames(air_20_winter[,i]))
#   # title(main=colnames(air_5_winter)[i])
# }
# mtext("Air temperatures - 20% NAs", side = 3, line = - 2, outer = TRUE)

# Imputation --------------------------------------------------------------


# X-SARIMA

air_5_winter_sarima=fill_sarima(air_5_winter,rh_winter,period = 24)
air_10_winter_sarima=fill_sarima(air_10_winter,rh_winter,period = 24)
air_20_winter_sarima=fill_sarima(air_20_winter,rh_winter,period = 24)

# windows()
# par(mfrow=c(6,6),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5_winter_sarima)){
#   plot(x=air_5_winter$time,y=as.vector(unlist(air_5_winter_sarima[,i])),type="l",col="red",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5_winter[,i]))
#   lines(x=air_5_winter$time,y=as.vector(unlist(air_5_winter[,i])),col="black")
#   # title(main=colnames(air_5_winter)[i])
# }
# mtext("Air temperatures - 5% NAs", side = 3, line = - 2, outer = TRUE)

# NAIVE

air_5_winter_naive=MA_imp(air_5_winter)
air_5_winter_naive=air_5_winter_naive$naive
#
air_10_winter_naive=MA_imp(air_10_winter)
air_10_winter_naive=air_10_winter_naive$naive
#
air_20_winter_naive=MA_imp(air_20_winter)
air_20_winter_naive=air_20_winter_naive$naive

# windows()
# par(mfrow=c(6,6),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5_winter_sarima)){
#   plot(x=air_5_winter$time,y=as.vector(unlist(air_5_winter_naive[,i])),type="l",col="red",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5_winter[,i]))
#   lines(x=air_5_winter$time,y=as.vector(unlist(air_5_winter[,i])),col="black")
#   # title(main=colnames(air_5_winter)[i])
# }
# mtext("Air temperatures - 5% NAs", side = 3, line = - 2, outer = TRUE)

# Linear Regression

air_5_winter_lr=lin_reg_imp(air_5_winter,rh_winter)
air_10_winter_lr=lin_reg_imp(air_10_winter,rh_winter)
air_20_winter_lr=lin_reg_imp(air_20_winter,rh_winter)

# windows()
# par(mfrow=c(6,6),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5_winter_sarima)){
#   plot(x=air_5_winter$time,y=as.vector(unlist(air_5_winter_lr[,i])),type="l",col="red",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5_winter[,i]))
#   lines(x=air_5_winter$time,y=as.vector(unlist(air_5_winter[,i])),col="black")
#   # title(main=colnames(air_5_winter)[i])
# }
# mtext("Air temperatures - 5% NAs", side = 3, line = - 2, outer = TRUE)

## Decomposition ##

# X-SARIMA 

air_5_winter_loess_sarima=LOESS.df(air_5_winter_sarima)
air_10_winter_loess_sarima=LOESS.df(air_10_winter_sarima)
air_20_winter_loess_sarima=LOESS.df(air_20_winter_sarima)

# NAIVE

air_5_winter_loess_naive=LOESS.df(air_5_winter_naive)
air_10_winter_loess_naive=LOESS.df(air_10_winter_naive)
air_20_winter_loess_naive=LOESS.df(air_20_winter_naive)

# Linear Regression

air_5_winter_loess_lr=LOESS.df(air_5_winter_lr)
air_10_winter_loess_lr=LOESS.df(air_10_winter_lr)
air_20_winter_loess_lr=LOESS.df(air_20_winter_lr)


# CV Kriging --------------------------------------------------------------


indx=2:(n_stations+1)

#start = Sys.time()

# X-SARIMA
krig_air_5_winter_sarima <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air_5_winter_loess_sarima$residuals,
                                                            locations3,
                                                            relevant_times=which(is.na(air_5_winter[,x]))
                                                            )
                                         ,
                                         mc.cores = parallel::detectCores()-1)
# end = Sys.time()
# end-start

save(krig_air_5_winter_sarima,file="krig_air_5_winter_sarima.RData")
#rm(krig_air_5_winter_sarima)

krig_air_10_winter_sarima <- parallel::mclapply(indx,
                                               function(x)CV_STkr(x,
                                                                  air_10_winter_loess_sarima$residuals,
                                                                  locations3,
                                                                  relevant_times=which(is.na(air_10_winter[,x]))
                                                                  ),
                                               mc.cores = parallel::detectCores()-1)
save(krig_air_10_winter_sarima,file="krig_air_10_winter_sarima.RData")
#rm(krig_air_10_winter_sarima)

krig_air_20_winter_sarima <- parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_20_winter_loess_sarima$residuals,
                                                                   locations3
                                                                   ,
                                                                   relevant_times=which(is.na(air_20_winter[,x]))
                                                                   ),
                                                mc.cores = parallel::detectCores()-1)
save(krig_air_20_winter_sarima,file="krig_air_20_winter_sarima.RData")
#rm(krig_air_20_winter_sarima)
# NAIVE

krig_air_5_winter_naive <- parallel::mclapply(indx,
                                               function(x)CV_STkr(x,
                                                                  air_5_winter_loess_naive$residuals,
                                                                  locations3,
                                                                  relevant_times=which(is.na(air_5_winter[,x]))
                                                                  ),
                                               mc.cores = parallel::detectCores()-1)
save(krig_air_5_winter_naive,file="krig_air_5_winter_naive.RData")

krig_air_10_winter_naive <- parallel::mclapply(indx,
                                              function(x)CV_STkr(x,
                                                                 air_10_winter_loess_naive$residuals,
                                                                 locations3,
                                                                 relevant_times=which(is.na(air_10_winter[,x]))
                                                                 ),
                                              mc.cores = parallel::detectCores()-1)
save(krig_air_10_winter_naive,file="krig_air_10_winter_naive.RData")

krig_air_20_winter_naive <- parallel::mclapply(indx,
                                              function(x)CV_STkr(x,
                                                                 air_20_winter_loess_naive$residuals,
                                                                 locations3,
                                                                 relevant_times=which(is.na(air_20_winter[,x]))
                                                                 ),
                                              mc.cores = parallel::detectCores()-1)

save(krig_air_20_winter_naive,file="krig_air_20_winter_naive.RData")

# Lin Regression

krig_air_5_winter_lr <- parallel::mclapply(indx,
                                              function(x)CV_STkr(x,
                                                                 air_5_winter_loess_lr$residuals,
                                                                 locations3,
                                                                 relevant_times=which(is.na(air_5_winter[,x]))
                                                                 ),
                                              mc.cores = parallel::detectCores()-1)

save(krig_air_5_winter_lr,file="krig_air_5_winter_lr.RData")
krig_air_10_winter_lr <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air_10_winter_loess_lr$residuals,
                                                              locations3,
                                                              relevant_times=which(is.na(air_10_winter[,x]))
                                                              ),
                                           mc.cores = parallel::detectCores()-1)

save(krig_air_10_winter_lr,file="krig_air_10_winter_lr.RData")

krig_air_20_winter_lr <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air_20_winter_loess_lr$residuals,
                                                              locations3,
                                                              relevant_times=which(is.na(air_20_winter[,x]))
                                                              ),
                                           mc.cores = parallel::detectCores()-1)

save(krig_air_20_winter_lr,file="krig_air_20_winter_lr.RData")

end = Sys.time()
elapsed=end-start


# Recover trend and seas --------------------------------------------------

# X-SARIMA
air_5_sarima_winter_recover=df_recover(krig_air_5_winter_sarima,
                                air_5_winter_loess_sarima, 
                                loess=T,
                                locations3)

air_10_sarima_winter_recover=df_recover(krig_air_10_winter_sarima,
                                       air_10_winter_loess_sarima, 
                                       loess=T,
                                       locations3)

air_20_sarima_winter_recover=df_recover(krig_air_20_winter_sarima,
                                       air_20_winter_loess_sarima, 
                                       loess=T,
                                       locations3)

# NAIVE

air_5_naive_winter_recover=df_recover(krig_air_5_winter_naive,
                                      air_5_winter_loess_naive, 
                                      loess=T,
                                      locations3)

air_10_naive_winter_recover=df_recover(krig_air_10_winter_naive,
                                       air_10_winter_loess_naive, 
                                       loess=T,
                                       locations3)

air_20_naive_winter_recover=df_recover(krig_air_20_winter_naive,
                                       air_20_winter_loess_naive, 
                                       loess=T,
                                       locations3)

# Lin Regression

air_5_lr_winter_recover=df_recover(krig_air_5_winter_lr,
                                   air_5_winter_loess_lr, 
                                   loess=T,
                                   locations3)

air_10_lr_winter_recover=df_recover(krig_air_10_winter_lr,
                                    air_10_winter_loess_lr, 
                                    loess=T,
                                    locations3)

air_20_lr_winter_recover=df_recover(krig_air_20_winter_lr,
                                    air_20_winter_loess_lr, 
                                    loess=T,
                                    locations3)

## Recover original dimension df (maybe later)
# air_5_winter_recover=air_5_winter
# air_5_winter_recover[is.na(air_5_winter_recover)]=0
# air_5_sarima_winter_recover_2=air_5_sarima_winter_recover[,-1]+air_5_winter_recover[,-1]
# air_5_sarima_winter_recover_2=data.frame(time=air_5_winter$time,
#                                          air_5_sarima_winter_recover_2)

# Lumped Kriging ----------------------------------------------------------

indx=2:(n_stations+1)
locations4=locations3
colnames(locations4)[1:3]=c("id","longitude","latitude")

start_lump=Sys.time()

krig_air_5_winter_lump <- parallel::mclapply(indx,
                                function(x)CV_lump.kr(x,air_5_winter,locations4),
                                mc.cores = parallel::detectCores()-1)
end_lump=Sys.time()
end_lump-start_lump

save(krig_air_5_winter_lump,file="krig_air_5_winter_lump.RData")

krig_air_10_winter_lump <- parallel::mclapply(indx,
                                             function(x)CV_lump.kr(x,air_10_winter,
                                                                   locations4),
                                             mc.cores = parallel::detectCores()-1)
save(krig_air_10_winter_lump,file="krig_air_10_winter_lump.RData")

krig_air_20_winter_lump <- parallel::mclapply(indx,
                                             function(x)CV_lump.kr(x,air_20_winter,
                                                                   locations4),
                                             mc.cores = parallel::detectCores()-1)
save(krig_air_20_winter_lump,file="krig_air_20_winter_lump.RData")


end_lump=Sys.time()

# II Spring II ------------------------------------------------------------------

spring_start=which(dat_air_2023$time=="2023-03-21 00:00:00 UTC")
spring_end=which(dat_air_2023$time=="2023-06-20 23:00:00 UTC")


# III Summer III ------------------------------------------------------------------

summer_start=which(dat_air_2023$time=="2023-06-21 00:00:00 UTC")
summer_end=which(dat_air_2023$time=="2023-09-22 23:00:00 UTC")


# IV Autumn IV ------------------------------------------------------------------

autumn_start=which(dat_air_2023$time=="2023-09-23 00:00:00 UTC")
autumn_end=which(dat_air_2023$time=="2023-12-20 23:00:00 UTC")