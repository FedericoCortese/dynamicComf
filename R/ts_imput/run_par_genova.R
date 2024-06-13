source("vers2/Utils2.R")
load("prov_genova_tempANDrh.RData")
map3=leaflet() %>% 
  addTiles() %>%
  addCircleMarkers(data=locations3,
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",Latitude, "<br>Lon:", Longitude
                                  , "<br>Height (m):", Height, "<br>Station:", `NOME STAZIONE`
                   )
  )
map3

# Keep 17 of these 34 stations
dat_temp_wide3_old=dat_temp_wide3
dat_rh_wide3_old=dat_rh_wide3
kstat=c("time","BARGAGLI","CROCETTADIORERO","DAVAGNA","FALLAROSA",
        "GENOVACENTROFUNZIONALE","GENOVAQUEZZI","GENOVAS.ILARIO",
        "ISOVERDE","MADONNADELLEGRAZIE","MONTECAPPELLINO",
        "MONTEPENNELLO","MONTOGGIO","SELLAGIASSINA","TANADORSO",
        "TORRIGLIAGARAVENTA","VALBREVENNAGORRA","VIGANEGO")
dat_temp_wide3=dat_temp_wide3[,kstat]
dat_rh_wide3=dat_rh_wide3[,kstat]
locations3_old=locations3
locations3=locations3[locations3$`NOME STAZIONE` %in% kstat,]

map3=leaflet() %>% 
  addTiles() %>%
  addCircleMarkers(data=locations3[locations3$`NOME STAZIONE` %in% kstat,],
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",Latitude, "<br>Lon:", Longitude
                                  , "<br>Height (m):", Height, "<br>Station:", `NOME STAZIONE`
                   )
  )
map3

# Descriptive analysis --------------------------------------------------

summary(dat_temp_wide3)
summary(dat_rh_wide3)

# Air temp
windows()
par(mfrow=c(5,7),mar=c(2,2,6,2))
for(i in 2:ncol(dat_temp_wide3)){
  plot(x=dat_temp_wide3$time,y=as.vector(unlist(dat_temp_wide3[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(dat_temp_wide3[,i]))
  title(main=colnames(dat_temp_wide3)[i])
}
mtext("Air temperatures", side = 3, line = - 2, outer = TRUE)

windows()
par(mfrow=c(1,1),mar=c(2,2,6,2))
i=7
plot(x=dat_temp_wide3$time,y=as.vector(unlist(dat_temp_wide3[,i])),type="l",col="black",
     xlab=" ",ylab=" ",
     main=colnames(dat_temp_wide3[,i]))
title(main=colnames(dat_temp_wide3)[i])

# Relative humidity
windows()
par(mfrow=c(5,7),mar=c(2,2,6,2))
for(i in 2:ncol(dat_rh_wide3)){
  plot(x=dat_rh_wide3$time,y=as.vector(unlist(dat_rh_wide3[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(dat_rh_wide3[,i]))
  title(main=colnames(dat_rh_wide3)[i])
}
mtext("Relative humidity", side = 3, line = - 2, outer = TRUE)

windows()
par(mfrow=c(1,1),mar=c(2,2,6,2))
i=7
plot(x=dat_rh_wide3$time,y=as.vector(unlist(dat_rh_wide3[,i])),type="l",col="black",
     xlab=" ",ylab=" ",
     main=colnames(dat_rh_wide3[,i]))
title(main=colnames(dat_rh_wide3)[i])

# Scatterplot of temperature and rel humidity for each station
windows()
par(mfrow=c(5,7),mar=c(2,2,6,2))
for(i in 2:ncol(dat_temp_wide3)){
  plot(x=as.vector(unlist(dat_temp_wide3[,i])),y=as.vector(unlist(dat_rh_wide3[,i])),
       #xlab="Temperature (°C)",ylab="Relative humidity (%)",
       main=colnames(dat_temp_wide3[,i]))
  title(main=colnames(dat_temp_wide3)[i])
}
title("Scatterplot of temperature and relative humidity for each station")

# Correlation between temperature and relative humidity 
diag(round(cor(dat_temp_wide3[,-1],dat_rh_wide3[,-1]),2))


# CV ----------------------------------------------------------------------

one_year=1:8760
air_year=dat_temp_wide3[one_year,]
rh_year=dat_rh_wide3[one_year,]

# First window  ----------------------------------------------------------

# first_seas=1:2190
# L=floor(dim(dat_rh_wide3)[1]/20)
L=dim(air_year)[1]
first_wdn=1:(L/6)
# few_stat=10+1

air_short=air_year[first_wdn,]
rh_short=rh_year[first_wdn,]

# air_short=dat_temp_wide3[one_year,]
# rh_short=dat_rh_wide3[one_year,]
# air_short=dat_temp_wide3[first_seas,]
# rh_short=dat_rh_wide3[first_seas,]

# air_short=dat_temp_wide3[1:1000,]
# rh_short=dat_rh_wide3[1:1000,]

source("vers2/Utils2.R")
TT=dim(air_short)[1]
#TT/2-7
na_len=TT*c(.05,.1,.2)
#na_start=rep(TT/2-7,3)
na_start=floor(TT-na_len)

air_5=air_short
air_10=air_short
air_20=air_short

# For cycle, randomly fill with NAs 
for(i in 2:ncol(air_short)){
  set.seed(i)
  na_st=sample(1:na_start[1],1)
  air_5[na_st:(na_st+na_len[1]),i]=NA
  na_st=sample(1:na_start[2],1)
  air_10[na_st:(na_st+na_len[2]),i]=NA
  na_st=sample(1:na_start[3],1)
  air_20[na_st:(na_st+na_len[3]),i]=NA
}

# air_5[na_start[1]:(na_start[1]+na_len[1]),2:ncol(air_5)]=NA
# air_10[na_start[2]:(na_start[2]+na_len[2]),2:ncol(air_10)]=NA
# air_20[na_start[3]:(na_start[3]+na_len[3]),2:ncol(air_20)]=NA

# Plot
windows()
par(mfrow=c(5,7),mar=c(2,2,6,2))
for(i in 2:ncol(air_5)){
  plot(x=air_5$time,y=as.vector(unlist(air_5[,i])),type="l",col="grey",
       xlab=" ",ylab=" ",
       main=colnames(air_5[,i]))
  title(main=colnames(air_5)[i])
}
mtext("Air temperatures - 5% NAs", side = 3, line = - 2, outer = TRUE)

windows()
par(mfrow=c(5,7),mar=c(2,2,6,2))
for(i in 2:ncol(air_10)){
  plot(x=air_10$time,y=as.vector(unlist(air_10[,i])),type="l",col="grey40",
       xlab=" ",ylab=" ",
       main=colnames(air_10[,i]))
  title(main=colnames(air_10)[i])
}
mtext("Air temperatures - 10% NAs", side = 3, line = - 2, outer = TRUE)

windows()
par(mfrow=c(5,7),mar=c(2,2,6,2))
for(i in 2:ncol(air_20)){
  plot(x=air_20$time,y=as.vector(unlist(air_20[,i])),type="l",col="grey20",
       xlab=" ",ylab=" ",
       main=colnames(air_20[,i]))
  title(main=colnames(air_20)[i])
}
mtext("Air temperatures - 20% NAs", side = 3, line = - 2, outer = TRUE)

save(air_5,air_10,air_20,file="air_tempGENOVA_NAs.RData")

# Load NA data ------------------------------------------------------------
load("air_tempGENOVA_NAs.RData")

# 2) Imputation without preliminary detrend-deseas ------------------------

# 2.1) SARIMA ------------------------------------------------------------------

sarima_times=rep(0,3)
st=Sys.time()
air_5_sarima=fill_sarima(air_5,rh_short,period = 24)
en=Sys.time()
sarima_times[1]=en-st

st=Sys.time()
air_10_sarima=fill_sarima(air_10,rh_short,period = 24)
en=Sys.time()
sarima_times[2]=en-st

st=Sys.time()
air_20_sarima=fill_sarima(air_20,rh_short,period = 24)
en=Sys.time()
sarima_times[3]=en-st

# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5)){
#   plot(x=air_5$time,y=as.vector(unlist(air_5_sarima[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5[,i]))
#   lines(x=air_5$time,y=as.vector(unlist(air_5[,i])),col="black")
#   title(main=colnames(air_5)[i])
# }
# mtext("Air temperatures - 5% NAs - SARIMA", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_10)){
#   plot(x=air_10$time,y=as.vector(unlist(air_10_sarima[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_10[,i]))
#   lines(x=air_10$time,y=as.vector(unlist(air_10[,i])),col="black")
#   title(main=colnames(air_10)[i])
# }
# mtext("Air temperatures - 10% NAs - SARIMA", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_20)){
#   plot(x=air_20$time,y=as.vector(unlist(air_20_sarima[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_20[,i]))
#   lines(x=air_20$time,y=as.vector(unlist(air_20[,i])),col="black")
#   title(main=colnames(air_20)[i])
# }
# mtext("Air temperatures - 20% NAs - SARIMA", side = 3, line = - 2, outer = TRUE)
# 


# 2.2) Temporal kriging --------------------------------------------------------

tkgr_times=rep(0,3)

st=Sys.time()
air_5_tkr=temp_krige(air_5,rh_short)
en=Sys.time()
tkgr_times[1]=en-st

st=Sys.time()
air_10_tkr=temp_krige(air_10,rh_short)
en=Sys.time()
tkgr_times[2]=en-st

st=Sys.time()
air_20_tkr=temp_krige(air_20,rh_short)
en=Sys.time()
tkgr_times[3]=en-st

# Plot
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5)){
#   plot(x=air_5$time,y=as.vector(unlist(air_5_tkr[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5[,i]))
#   lines(x=air_5$time,y=as.vector(unlist(air_5[,i])),col="black")
#   title(main=colnames(air_5)[i])
# }
# mtext("Air temperatures - 5% NAs - temporal kriging", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_10)){
#   plot(x=air_10$time,y=as.vector(unlist(air_10_tkr[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_10[,i]))
#   lines(x=air_10$time,y=as.vector(unlist(air_10[,i])),col="black")
#   title(main=colnames(air_10)[i])
# }
# mtext("Air temperatures - 10% NAs - temporal kriging", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_20)){
#   plot(x=air_20$time,y=as.vector(unlist(air_20_tkr[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_20[,i]))
#   lines(x=air_20$time,y=as.vector(unlist(air_20[,i])),col="black")
#   title(main=colnames(air_20)[i])
# }
# mtext("Air temperatures - 20% NAs - temporal kriging", side = 3, line = - 2, outer = TRUE)
# 


# 2.3) SDEM and NAIVE ----------------------------------------------------------

ma_times=rep(0,3)

st=Sys.time()
air_sdem5=MA_imp(air_5)
en=Sys.time()
ma_times[1]=en-st

st=Sys.time()
air_sdem10=MA_imp(air_10)
st=Sys.time()
ma_times[2]=en-st

st=Sys.time()
air_sdem20=MA_imp(air_20)
en=Sys.time()
ma_times[3]=en-st

air5_SDEM=air_sdem5$SDEM
air10_SDEM=air_sdem10$SDEM
air20_SDEM=air_sdem20$SDEM

# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5)){
#   plot(x=air_5$time,y=as.vector(unlist(air5_SDEM[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5[,i]))
#   lines(x=air_5$time,y=as.vector(unlist(air_5[,i])),col="black")
#   title(main=colnames(air_5)[i])
# }
# mtext("Air temperatures - 5% NAs - SDEM", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_10)){
#   plot(x=air_10$time,y=as.vector(unlist(air10_SDEM[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_10[,i]))
#   lines(x=air_10$time,y=as.vector(unlist(air_10[,i])),col="black")
#   title(main=colnames(air_10)[i])
# }
# mtext("Air temperatures - 10% NAs - SDEM", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_20)){
#   plot(x=air_20$time,y=as.vector(unlist(air20_SDEM[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_20[,i]))
#   lines(x=air_20$time,y=as.vector(unlist(air_20[,i])),col="black")
#   title(main=colnames(air_20)[i])
# }
# mtext("Air temperatures - 20% NAs - SDEM", side = 3, line = - 2, outer = TRUE)


air5_naive=air_sdem5$naive
air10_naive=air_sdem10$naive
air20_naive=air_sdem20$naive


# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5)){
#   plot(x=air_5$time,y=as.vector(unlist(air5_naive[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5[,i]))
#   lines(x=air_5$time,y=as.vector(unlist(air_5[,i])),col="black")
#   title(main=colnames(air_5)[i])
# }
# mtext("Air temperatures - 5% NAs - Naive", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_10)){
#   plot(x=air_10$time,y=as.vector(unlist(air10_naive[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_10[,i]))
#   lines(x=air_10$time,y=as.vector(unlist(air_10[,i])),col="black")
#   title(main=colnames(air_10)[i])
# }
# mtext("Air temperatures - 10% NAs - Naive", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_20)){
#   plot(x=air_20$time,y=as.vector(unlist(air20_naive[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_20[,i]))
#   lines(x=air_20$time,y=as.vector(unlist(air_20[,i])),col="black")
#   title(main=colnames(air_20)[i])
# }
# mtext("Air temperatures - 20% NAs - Naive", side = 3, line = - 2, outer = TRUE)
# 

# 2.4) Linear regression -------------------------------------------------------

air5_lr=lin_reg_imp(air_5,rh_short)
air10_lr=lin_reg_imp(air_10,rh_short)
air20_lr=lin_reg_imp(air_20,rh_short)

save.image("dat_genova_imput.Rdata")


#Plot
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_5)){
#   plot(x=air_5$time,y=as.vector(unlist(air5_lr[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_5[,i]))
#   lines(x=air_5$time,y=as.vector(unlist(air_5[,i])),col="black")
#   
#   title(main=colnames(air_5)[i])
# }
# mtext("Air temperatures - 5% NAs - Linear regression", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_10)){
#   plot(x=air_10$time,y=as.vector(unlist(air10_lr[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_10[,i]))
#   lines(x=air_10$time,y=as.vector(unlist(air_10[,i])),col="black")
#   title(main=colnames(air_10)[i])
# }
# mtext("Air temperatures - 10% NAs - Linear regression", side = 3, line = - 2, outer = TRUE)
# 
# windows()
# par(mfrow=c(4,3),mar=c(2,2,6,2))
# for(i in 2:ncol(air_20)){
#   plot(x=air_20$time,y=as.vector(unlist(air20_lr[,i])),col="red"
#        ,type="l",
#        xlab=" ",ylab=" ",
#        main=colnames(air_20[,i]))
#   lines(x=air_20$time,y=as.vector(unlist(air_20[,i])),col="black")
#   title(main=colnames(air_20)[i])
# }
# mtext("Air temperatures - 20% NAs - Linear regression", side = 3, line = - 2, outer = TRUE)

# Comparison
mean(unlist(as.vector(rmse(air_short,air_5_sarima))))
mean(unlist(as.vector(rmse(air_short,air_5_tkr))))
#mean(unlist(as.vector(rmse(air_short,air5_SDEM))))
mean(unlist(as.vector(rmse(air_short,air5_lr))))
mean(unlist(as.vector(rmse(air_short,air5_naive))))

mean(unlist(as.vector(rmse(air_short,air_10_sarima))))
mean(unlist(as.vector(rmse(air_short,air_10_tkr))))
#mean(unlist(as.vector(rmse(air_short,air10_SDEM))))
mean(unlist(as.vector(rmse(air_short,air10_lr))))
mean(unlist(as.vector(rmse(air_short,air10_naive))))


mean(unlist(as.vector(rmse(air_short,air_20_sarima))))
mean(unlist(as.vector(rmse(air_short,air_20_tkr))),na.rm = T)
#mean(unlist(as.vector(rmse(air_short,air20_SDEM))))
mean(unlist(as.vector(rmse(air_short,air20_lr))))
mean(unlist(as.vector(rmse(air_short,air20_naive))))

# 3) Universal kriging (parallel) -----------------------------------------

locations2=locations3

indx=2:ncol(air_short)
#x=2


start = Sys.time()
air5_sarima_full <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air_5_sarima,
                                                          locations2,
                                                          ordinary=F),
                                       mc.cores = parallel::detectCores()-1)
end = Sys.time()
elapsed_air5_sarima_full=end-start

save(elapsed_air5_sarima_full,air5_sarima_full,file="elapsed_air5_sarima_full.RData")

start = Sys.time()
air5_tkr_full <- parallel::mclapply(indx,
                                    function(x)CV_STkr(x,
                                                       air_5_tkr,
                                                       locations2,
                                                       ordinary=F),
                                    mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_tkr_full=end-start

save(elapsed_air5_tkr_full,air5_tkr_full,file="elapsed_air5_tkr_full.RData")

start = Sys.time()
air5_SDEM_full <- parallel::mclapply(indx,
                                     function(x)CV_STkr(x,
                                                        air5_SDEM,
                                                        locations2,
                                                        ordinary=F),
                                     mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_SDEM_full=end-start

save(elapsed_air5_SDEM_full,air5_SDEM_full,file="elapsed_air5_SDEM_full.RData")

start = Sys.time()
air5_naive_full <- parallel::mclapply(indx,
                                      function(x)CV_STkr(x,
                                                         air5_naive,
                                                         locations2,
                                                         ordinary=F),
                                      mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_naive_full=end-start

save(elapsed_air5_naive_full,air5_naive_full,file="elapsed_air5_naive_full.RData")

start = Sys.time()
air10_sarima_full <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air_10_sarima,
                                                           locations2,
                                                           ordinary=F),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_sarima_full=end-start

save(elapsed_air10_sarima_full,air10_sarima_full,file="elapsed_air10_sarima_full.RData")

start = Sys.time()
air10_tkr_full <- parallel::mclapply(indx,
                                     function(x)CV_STkr(x,
                                                        air_10_tkr,
                                                        locations2,
                                                        ordinary=F),
                                     mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_tkr_full=end-start

save(elapsed_air10_tkr_full,air10_tkr_full,file="elapsed_air10_tkr_full.RData")

start = Sys.time()
air10_SDEM_full <- parallel::mclapply(indx,
                                      function(x)CV_STkr(x,
                                                         air10_SDEM,
                                                         locations2,
                                                         ordinary=F),
                                      mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_SDEM_full=end-start

save(elapsed_air10_SDEM_full,air10_SDEM_full,file="elapsed_air10_SDEM_full.RData")

start = Sys.time()
air10_naive_full <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air10_naive,
                                                          locations2,
                                                          ordinary=F),
                                       mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_naive_full=end-start

save(elapsed_air10_naive_full,air10_naive_full,file="elapsed_air10_naive_full.RData")

start = Sys.time()
air20_sarima_full <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air_20_sarima,
                                                           locations2,
                                                           ordinary=F),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_sarima_full=end-start

save(elapsed_air20_sarima_full,air20_sarima_full,file="elapsed_air20_sarima_full.RData")

start = Sys.time()
air20_tkr_full <- parallel::mclapply(indx,
                                     function(x)CV_STkr(x,
                                                        air_20_tkr,
                                                        locations2,
                                                        ordinary=F),
                                     mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_tkr_full=end-start

save(elapsed_air20_tkr_full,air20_tkr_full,file="elapsed_air20_tkr_full.RData")

start = Sys.time()
air20_SDEM_full <- parallel::mclapply(indx,
                                      function(x)CV_STkr(x,
                                                         air20_SDEM,
                                                         locations2,
                                                         ordinary=F),
                                      mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_SDEM_full=end-start

save(elapsed_air20_SDEM_full,air20_SDEM_full,file="elapsed_air20_SDEM_full.RData")

start = Sys.time()
air20_naive_full <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air20_naive,
                                                          locations2,
                                                          ordinary=F),
                                       mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_naive_full=end-start

save(elapsed_air20_naive_full,air20_naive_full,file="elapsed_air20_naive_full.RData")

start = Sys.time()
air5_linreg_full <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air5_lr,
                                                          locations2,
                                                          ordinary=F),
                                       mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_lr_full=end-start

save(elapsed_air5_lr_full,air5_linreg_full,file="elapsed_air5_lr_full.RData")

start = Sys.time()
air10_linreg_full <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air10_lr,
                                                           locations2,
                                                           ordinary=F),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_lr_full=end-start

save(elapsed_air10_lr_full,air10_linreg_full,file="elapsed_air10_lr_full.RData")

start = Sys.time()
air20_linreg_full <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air20_lr,
                                                           locations2,
                                                           ordinary=F),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_lr_full=end-start

save(elapsed_air20_lr_full,air20_linreg_full,file="elapsed_air20_lr_full.RData")

# 4) Preliminary detrend-deseas -------------------------------------------


## LOESS
# air 5%
air5_loess_sarima=LOESS.df(air_5_sarima)

#air_5_tkr=select(air_5_tkr,subset=-c(indx,knot))
air5_loess_tkr=LOESS.df(air_5_tkr)
air5_loess_SDEM=LOESS.df(air5_SDEM)
air5_loess_naive=LOESS.df(air5_naive)

# air 10%
air10_loess_sarima=LOESS.df(air_10_sarima)
#air_10_tkr=select(air_10_tkr,subset=-c(indx,knot))
air10_loess_tkr=LOESS.df(air_10_tkr)
air10_loess_SDEM=LOESS.df(air10_SDEM)
air10_loess_naive=LOESS.df(air10_naive)

# air 20%
air20_loess_sarima=LOESS.df(air_20_sarima)
#air_20_tkr=select(air_20_tkr,subset=-c(indx,knot))
air20_loess_tkr=LOESS.df(air_20_tkr)
air20_loess_SDEM=LOESS.df(air20_SDEM)
air20_loess_naive=LOESS.df(air20_naive)

# NEW linreg
air5_loess_lr=LOESS.df(air5_lr)
air10_loess_lr=LOESS.df(air10_lr)
air20_loess_lr=LOESS.df(air20_lr)

## HW
# air 5%
air5_hw_sarima=HoltWint.df(air_5_sarima,24)
air5_hw_tkr=HoltWint.df(air_5_tkr,24)
air5_hw_SDEM=HoltWint.df(air5_SDEM,24)
air5_hw_naive=HoltWint.df(air5_naive,24)
# air 10%
air10_hw_sarima=HoltWint.df(air_10_sarima,24)
air10_hw_tkr=HoltWint.df(air_10_tkr,24)
air10_hw_SDEM=HoltWint.df(air10_SDEM,24)
air10_hw_naive=HoltWint.df(air10_naive,24)
# air 20%
air20_hw_sarima=HoltWint.df(air_20_sarima,24)
air20_hw_tkr=HoltWint.df(air_20_tkr,24)
air20_hw_SDEM=HoltWint.df(air20_SDEM,24)
air20_hw_naive=HoltWint.df(air20_naive,24)

#NEW linreg
air5_hw_lr=HoltWint.df(air5_lr,24)
air10_hw_lr=HoltWint.df(air10_lr,24)
air20_hw_lr=HoltWint.df(air20_lr,24)

# 5) Ordinary kriging --------------------------------------------------------

#hw
indx=2:ncol(air_short)

start = Sys.time()
air5_sarima_hw_res <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air5_hw_sarima$residuals,
                                                            locations),
                                         mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_sarima_hw_res=end-start

start = Sys.time()
air5_tkr_hw_res <- parallel::mclapply(indx,
                                      function(x)CV_STkr(x,
                                                         air5_hw_tkr$residuals,
                                                         locations),
                                      mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_tkr_hw_res=end-start

start = Sys.time()
air5_SDEM_hw_res <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air5_hw_SDEM$residuals,
                                                          locations),
                                       mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_SDEM_hw_res=end-start

start = Sys.time()
air5_naive_hw_res <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air5_hw_naive$residuals,
                                                           locations),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_naive_hw_res=end-start


start = Sys.time()
air10_sarima_hw_res <- parallel::mclapply(indx,
                                          function(x)CV_STkr(x,
                                                             air10_hw_sarima$residuals,
                                                             locations),
                                          mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_sarima_hw_res=end-start

start = Sys.time()
air10_tkr_hw_res <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air10_hw_tkr$residuals,
                                                          locations),
                                       mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_tkr_hw_res=end-start

start = Sys.time()
air10_SDEM_hw_res <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air10_hw_SDEM$residuals,
                                                           locations),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_SDEM_hw_res=end-start

start = Sys.time()
air10_naive_hw_res <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air10_hw_naive$residuals,
                                                            locations),
                                         mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_naive_hw_res=end-start


start = Sys.time()
air20_sarima_hw_res <- parallel::mclapply(indx,
                                          function(x)CV_STkr(x,
                                                             air20_hw_sarima$residuals,
                                                             locations),
                                          mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_sarima_hw_res=end-start

start = Sys.time()
air20_tkr_hw_res <- parallel::mclapply(indx,
                                       function(x)CV_STkr(x,
                                                          air20_hw_tkr$residuals,
                                                          locations),
                                       mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_tkr_hw_res=end-start

start = Sys.time()
air20_SDEM_hw_res <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air20_hw_SDEM$residuals,
                                                           locations),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_SDEM_hw_res=end-start

start = Sys.time()
air20_naive_hw_res <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air20_hw_naive$residuals,
                                                            locations),
                                         mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_naive_hw_res=end-start

#loess
start = Sys.time()
air5_sarima_loess_res <- parallel::mclapply(indx,
                                            function(x)CV_STkr(x,
                                                               air5_loess_sarima$residuals,
                                                               locations),
                                            mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_sarima_loess_res=end-start

start = Sys.time()
air5_tkr_loess_res <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air5_loess_tkr$residuals,
                                                            locations),
                                         mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_tkr_loess_res=end-start

start = Sys.time()
air5_SDEM_loess_res <- parallel::mclapply(indx,
                                          function(x)CV_STkr(x,
                                                             air5_loess_SDEM$residuals,
                                                             locations),
                                          mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_SDEM_loess_res=end-start

start = Sys.time()
air5_naive_loess_res <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air5_loess_naive$residuals,
                                                              locations),
                                           mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_naive_loess_res=end-start


start = Sys.time()
air10_sarima_loess_res <- parallel::mclapply(indx,
                                             function(x)CV_STkr(x,
                                                                air10_loess_sarima$residuals,
                                                                locations),
                                             mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_sarima_loess_res=end-start

start = Sys.time()
air10_tkr_loess_res <- parallel::mclapply(indx,
                                          function(x)CV_STkr(x,
                                                             air10_loess_tkr$residuals,
                                                             locations),
                                          mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_tkr_loess_res=end-start

start = Sys.time()
air10_SDEM_loess_res <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air10_loess_SDEM$residuals,
                                                              locations),
                                           mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_SDEM_loess_res=end-start

start = Sys.time()
air10_naive_loess_res <- parallel::mclapply(indx,
                                            function(x)CV_STkr(x,
                                                               air10_loess_naive$residuals,
                                                               locations),
                                            mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_naive_loess_res=end-start


start = Sys.time()
air20_sarima_loess_res <- parallel::mclapply(indx,
                                             function(x)CV_STkr(x,
                                                                air20_loess_sarima$residuals,
                                                                locations),
                                             mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_sarima_loess_res=end-start

start = Sys.time()
air20_tkr_loess_res <- parallel::mclapply(indx,
                                          function(x)CV_STkr(x,
                                                             air20_loess_tkr$residuals,
                                                             locations),
                                          mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_tkr_loess_res=end-start

start = Sys.time()
air20_SDEM_loess_res <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air20_loess_SDEM$residuals,
                                                              locations),
                                           mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_SDEM_loess_res=end-start

start = Sys.time()
air20_naive_loess_res <- parallel::mclapply(indx,
                                            function(x)CV_STkr(x,
                                                               air20_loess_naive$residuals,
                                                               locations),
                                            mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_naive_loess_res=end-start

# HW
start = Sys.time()
air5_lr_hw_res <- parallel::mclapply(indx,
                                     function(x)CV_STkr(x,
                                                        air5_hw_lr$residuals,
                                                        locations),
                                     mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_lr_hw_res=end-start

start = Sys.time()
air10_lr_hw_res <- parallel::mclapply(indx,
                                      function(x)CV_STkr(x,
                                                         air10_hw_lr$residuals,
                                                         locations),
                                      mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_lr_hw_res=end-start

start = Sys.time()
air20_lr_hw_res <- parallel::mclapply(indx,
                                      function(x)CV_STkr(x,
                                                         air20_hw_lr$residuals,
                                                         locations),
                                      mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_lr_hw_res=end-start

#LOESS
start = Sys.time()
air5_lr_loess_res <- parallel::mclapply(indx,
                                        function(x)CV_STkr(x,
                                                           air5_loess_lr$residuals,
                                                           locations),
                                        mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air5_lr_loess_res=end-start

start = Sys.time()
air10_lr_loess_res <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air10_loess_lr$residuals,
                                                            locations),
                                         mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air10_lr_loess_res=end-start

start = Sys.time()
air20_lr_loess_res <- parallel::mclapply(indx,
                                         function(x)CV_STkr(x,
                                                            air20_loess_lr$residuals,
                                                            locations),
                                         mc.cores = parallel::detectCores())
end = Sys.time()
elapsed_air20_lr_loess_res=end-start

save.image("run_parallel_Genova2.RData")

