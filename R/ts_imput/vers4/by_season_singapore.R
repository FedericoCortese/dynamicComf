# Set time zone
Sys.setenv(TZ='Asia/Singapore')
source("Utils3.R")

# Load data ---------------------------------------------------------------

load("air_short.Rdata")
load("rh_short.Rdata")
load("locations.Rdata")
locations3=locations

locations3$stat_numb=(1:nrow(locations3))+1

map3=leaflet() %>% 
  addTiles() %>%
  addCircleMarkers(data=locations3[locations3$id %in% colnames(air_short),],
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",latitude, "<br>Lon:", longitude
                                  , "<br>Station:", id,"<br>Number:",stat_numb
                   )
  )
map3

n_stations=ncol(air_short)-1


# Fill with NA ------------------------------------------------------------
locations4=locations3[,c("id","latitude","longitude")]
dat_air_winter=air_short
rh_winter=rh_short

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

windows()
par(mfrow=c(3,4),mar=c(2,2,6,2))
for(i in 2:ncol(air_5_winter)){
  plot(x=air_5_winter$time,y=as.vector(unlist(air_5_winter[,i])),type="l",col="black",
       xlab=" ",ylab=" ",
       main=colnames(air_5_winter[,i]))
 # title(main=colnames(air_5_winter)[i])
}
mtext("Air temperatures - 5% NAs", side = 3, line = - 2, outer = TRUE)
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
                                                            locations4,
                                                            relevant_times=which(is.na(air_5_winter[,x]))
                                                            )
                                         ,
                                         mc.cores = parallel::detectCores()-1)

krig_air_5_winter_sarima_full <- parallel::mclapply(indx,
                                               function(x)CV_STkr(x,
                                                                  air_5_winter_sarima,
                                                                  #air_5_winter_loess_sarima$residuals,
                                                                  locations4,
                                                                  relevant_times=which(is.na(air_5_winter[,x]))
                                               ,ordinary=F)
                                               ,
                                               mc.cores = parallel::detectCores()-1)

# end = Sys.time()
# end-start


krig_air_10_winter_sarima <- parallel::mclapply(indx,
                                               function(x)CV_STkr(x,
                                                                  air_10_winter_loess_sarima$residuals,
                                                                  locations4,
                                                                  relevant_times=which(is.na(air_10_winter[,x]))
                                                                  ),
                                               mc.cores = parallel::detectCores()-1)

krig_air_10_winter_sarima_full=parallel::mclapply(indx,
                                                  function(x)CV_STkr(x,
                                                                     air_10_winter_sarima,
                                                                     locations4,
                                                                     relevant_times=which(is.na(air_10_winter[,x]))
                                                                     ,ordinary=F),
                                                  mc.cores = parallel::detectCores()-1)


krig_air_20_winter_sarima <- parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_20_winter_loess_sarima$residuals,
                                                                   locations4
                                                                   ,
                                                                   relevant_times=which(is.na(air_20_winter[,x]))
                                                                   ),
                                                mc.cores = parallel::detectCores()-1)

krig_air_20_winter_sarima_full=parallel::mclapply(indx,
                                                  function(x)CV_STkr(x,
                                                                     air_20_winter_sarima,
                                                                     locations4,
                                                                     relevant_times=which(is.na(air_20_winter[,x]))
                                                                     ,ordinary=F),
                                                  mc.cores = parallel::detectCores()-1)

# NAIVE

krig_air_5_winter_naive <- parallel::mclapply(indx,
                                               function(x)CV_STkr(x,
                                                                  air_5_winter_loess_naive$residuals,
                                                                  locations4,
                                                                  relevant_times=which(is.na(air_5_winter[,x]))
                                                                  ),
                                               mc.cores = parallel::detectCores()-1)

krig_air_5_winter_naive_full=parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_5_winter_naive,
                                                                   locations4,
                                                                   relevant_times=which(is.na(air_5_winter[,x]))
                                                                   ,ordinary=F),
                                                mc.cores = parallel::detectCores()-1)

krig_air_10_winter_naive <- parallel::mclapply(indx,
                                              function(x)CV_STkr(x,
                                                                 air_10_winter_loess_naive$residuals,
                                                                 locations4,
                                                                 relevant_times=which(is.na(air_10_winter[,x]))
                                                                 ),
                                              mc.cores = parallel::detectCores()-1)

krig_air_10_winter_naive_full=parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_10_winter_naive,
                                                                   locations4,
                                                                   relevant_times=which(is.na(air_10_winter[,x]))
                                                                   ,ordinary=F),
                                                mc.cores = parallel::detectCores()-1)

krig_air_20_winter_naive <- parallel::mclapply(indx,
                                              function(x)CV_STkr(x,
                                                                 air_20_winter_loess_naive$residuals,
                                                                 locations4,
                                                                 relevant_times=which(is.na(air_20_winter[,x]))
                                                                 ),
                                              mc.cores = parallel::detectCores()-1)

krig_air_20_winter_naive_full=parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_20_winter_naive,
                                                                   locations4,
                                                                   relevant_times=which(is.na(air_20_winter[,x]))
                                                                   ,ordinary=F),
                                                mc.cores = parallel::detectCores()-1)

# Lin Regression

krig_air_5_winter_lr <- parallel::mclapply(indx,
                                              function(x)CV_STkr(x,
                                                                 air_5_winter_loess_lr$residuals,
                                                                 locations4,
                                                                 relevant_times=which(is.na(air_5_winter[,x]))
                                                                 ),
                                              mc.cores = parallel::detectCores()-1)

krig_air_5_winter_lr_full=parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_5_winter_lr,
                                                                   locations4,
                                                                   relevant_times=which(is.na(air_5_winter[,x]))
                                                                   ,ordinary=F),
                                                mc.cores = parallel::detectCores()-1)

krig_air_10_winter_lr <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air_10_winter_loess_lr$residuals,
                                                              locations4,
                                                              relevant_times=which(is.na(air_10_winter[,x]))
                                                              ),
                                           mc.cores = parallel::detectCores()-1)

krig_air_10_winter_lr_full=parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_10_winter_lr,
                                                                   locations4,
                                                                   relevant_times=which(is.na(air_10_winter[,x]))
                                                                   ,ordinary=F),
                                                mc.cores = parallel::detectCores()-1)

krig_air_20_winter_lr <- parallel::mclapply(indx,
                                           function(x)CV_STkr(x,
                                                              air_20_winter_loess_lr$residuals,
                                                              locations4,
                                                              relevant_times=which(is.na(air_20_winter[,x]))
                                                              ),
                                           mc.cores = parallel::detectCores()-1)

krig_air_20_winter_lr_full=parallel::mclapply(indx,
                                                function(x)CV_STkr(x,
                                                                   air_20_winter_lr,
                                                                   locations4,
                                                                   relevant_times=which(is.na(air_20_winter[,x]))
                                                                   ,ordinary=F),
                                                mc.cores = parallel::detectCores()-1)

end = Sys.time()
elapsed=end-start


# Recover trend and seas --------------------------------------------------

locations5=locations4[locations4$id%in%colnames(dat_air_winter),]
locations5=locations5[,c("id","longitude","latitude")]

# X-SARIMA
air_5_sarima_winter_recover=df_recover(krig_air_5_winter_sarima,
                                air_5_winter_loess_sarima, 
                                loess=T,
                                locations2=locations5)

air_10_sarima_winter_recover=df_recover(krig_air_10_winter_sarima,
                                       air_10_winter_loess_sarima, 
                                       loess=T,
                                       locations5)

air_20_sarima_winter_recover=df_recover(krig_air_20_winter_sarima,
                                       air_20_winter_loess_sarima, 
                                       loess=T,
                                       locations5)

# NAIVE

air_5_naive_winter_recover=df_recover(krig_air_5_winter_naive,
                                      air_5_winter_loess_naive, 
                                      loess=T,
                                      locations5)

air_10_naive_winter_recover=df_recover(krig_air_10_winter_naive,
                                       air_10_winter_loess_naive, 
                                       loess=T,
                                       locations5)

air_20_naive_winter_recover=df_recover(krig_air_20_winter_naive,
                                       air_20_winter_loess_naive, 
                                       loess=T,
                                       locations5)

# Lin Regression

air_5_lr_winter_recover=df_recover(krig_air_5_winter_lr,
                                   air_5_winter_loess_lr, 
                                   loess=T,
                                   locations5)

air_10_lr_winter_recover=df_recover(krig_air_10_winter_lr,
                                    air_10_winter_loess_lr, 
                                    loess=T,
                                    locations5)

air_20_lr_winter_recover=df_recover(krig_air_20_winter_lr,
                                    air_20_winter_loess_lr, 
                                    loess=T,
                                    locations5)

## Recover original dimension df (maybe later)
# air_5_winter_recover=air_5_winter
# air_5_winter_recover[is.na(air_5_winter_recover)]=0
# air_5_sarima_winter_recover_2=air_5_sarima_winter_recover[,-1]+air_5_winter_recover[,-1]
# air_5_sarima_winter_recover_2=data.frame(time=air_5_winter$time,
#                                          air_5_sarima_winter_recover_2)

# Lumped Kriging ----------------------------------------------------------

indx=2:(n_stations+1)
#locations4=locations3
locations4=locations4[,c("id","longitude","latitude")]

start_lump=Sys.time()

krig_air_5_winter_lump <- parallel::mclapply(indx,
                                function(x)CV_lump.kr(x,air_5_winter,locations4),
                                mc.cores = parallel::detectCores()-1)
end_lump=Sys.time()
end_lump-start_lump

save(krig_air_5_winter_lump,file="krig_air_5_winter_lump.RData")

air_5_lump_winter=data.frame(time=dat_air_winter$time)
for(i in 1:length(krig_air_5_winter_lump)){
  pred=krig_air_5_winter_lump[[i]]$y[,1]
  #vars=krig_air_5_winter_lump[[i]]$y[,2]
  temp_times=dat_air_winter$time[which(is.na(air_5_winter[,(i+1)]))]
  temp=data.frame(temp_times,pred)
  colnames(temp)=c("time",krig_air_5_winter_lump[[i]]$stat)
  air_5_lump_winter=merge(air_5_lump_winter,temp,by="time",all.x=TRUE)
}

krig_air_10_winter_lump <- parallel::mclapply(indx,
                                             function(x)CV_lump.kr(x,air_10_winter,
                                                                   locations4),
                                             mc.cores = parallel::detectCores()-1)
save(krig_air_10_winter_lump,file="krig_air_10_winter_lump.RData")

air_10_lump_winter=data.frame(time=dat_air_winter$time)
for(i in 1:length(krig_air_10_winter_lump)){
  pred=krig_air_10_winter_lump[[i]]$y[,1]
  #vars=krig_air_10_winter_lump[[i]]$y[,2]
  temp_times=dat_air_winter$time[which(is.na(air_10_winter[,(i+1)]))]
  temp=data.frame(temp_times,pred)
  colnames(temp)=c("time",krig_air_10_winter_lump[[i]]$stat)
  air_10_lump_winter=merge(air_10_lump_winter,temp,by="time",all.x=TRUE)
}

krig_air_20_winter_lump <- parallel::mclapply(indx,
                                             function(x)CV_lump.kr(x,air_20_winter,
                                                                   locations4),
                                             mc.cores = parallel::detectCores()-1)
save(krig_air_20_winter_lump,file="krig_air_20_winter_lump.RData")

air_20_lump_winter=data.frame(time=dat_air_winter$time)
for(i in 1:length(krig_air_20_winter_lump)){
  pred=krig_air_20_winter_lump[[i]]$y[,1]
  #vars=krig_air_20_winter_lump[[i]]$y[,2]
  temp_times=dat_air_winter$time[which(is.na(air_20_winter[,(i+1)]))]
  temp=data.frame(temp_times,pred)
  colnames(temp)=c("time",krig_air_20_winter_lump[[i]]$stat)
  air_20_lump_winter=merge(air_20_lump_winter,temp,by="time",all.x=TRUE)
}

end_lump=Sys.time()

# Line plots --------------------------------------------------------------


# 5%
# X-SARIMA
pdf("5NAs-X-SARIMA_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_5_winter, 
               data_pred = air_5_sarima_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="5% NAs - SARIMA")
dev.off()



# NAIVE
pdf("5NAs-NAIVE_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_5_winter, 
               data_pred = air_5_naive_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="5% NAs - NAIVE")
dev.off()

# LR
pdf("5NAs-LR_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_5_winter, 
               data_pred = air_5_lr_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="5% NAs - LR")
dev.off()

# LUMP
pdf("5NAs-LUMP_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_5_winter, 
               data_pred = air_5_lump_winter,
               ncol=4,
               x_lab_size=8,
               title="5% NAs - LUMP")
dev.off()

# 10%
# X-SARIMA
pdf("10NAs-X-SARIMA_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_10_winter, 
               data_pred = air_10_sarima_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="10% NAs - SARIMA")
dev.off()

# NAIVE
pdf("10NAs-NAIVE_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_10_winter, 
               data_pred = air_10_naive_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="10% NAs - NAIVE")
dev.off()

# LR
pdf("10NAs-LR_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_10_winter, 
               data_pred = air_10_lr_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="10% NAs - LR")
dev.off()

# LUMP
pdf("10NAs-LUMP_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_10_winter, 
               data_pred = air_10_lump_winter,
               ncol=4,
               x_lab_size=8,
               title="10% NAs - LUMP")
dev.off()


# 20%
# X-SARIMA
pdf("20NAs-X-SARIMA_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_20_winter, 
               data_pred = air_20_sarima_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="20% NAs - SARIMA")
dev.off()

# NAIVE
pdf("20NAs-NAIVE_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_20_winter, 
               data_pred = air_20_naive_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="20% NAs - NAIVE")
dev.off()

# LR
pdf("20NAs-LR_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_20_winter, 
               data_pred = air_20_lr_winter_recover,
               ncol=4,
               x_lab_size=8,
               title="20% NAs - LR")
dev.off()

# LUMP
pdf("20NAs-LUMP_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
generate_plots(data_true = dat_air_winter, 
               data_NA = air_20_winter, 
               data_pred = air_20_lump_winter,
               ncol=4,
               x_lab_size=8,
               title="20% NAs - LUMP")
dev.off()


# RMSE --------------------------------------------------------------------

# 5%
MSE_air_5_lump_winter=(dat_air_winter[,-1]-air_5_lump_winter[,-1])^2
RMSE_air_5_lump_bystat=sqrt(apply(MSE_air_5_lump_winter,2,mean,na.rm=T))
mean(RMSE_air_5_lump_bystat)

MSE_air_5_sarima_winter=(dat_air_winter[,-1]-air_5_sarima_winter_recover[,-1])^2
RMSE_air_5_sarima_bystat=sqrt(apply(MSE_air_5_sarima_winter,2,mean,na.rm=T))
mean(RMSE_air_5_sarima_bystat)

MSE_air_5_naive_winter=(dat_air_winter[,-1]-air_5_naive_winter_recover[,-1])^2
RMSE_air_5_naive_bystat=sqrt(apply(MSE_air_5_naive_winter,2,mean,na.rm=T))
mean(RMSE_air_5_naive_bystat)

MSE_air_5_lr_winter=(dat_air_winter[,-1]-air_5_lr_winter_recover[,-1])^2
RMSE_air_5_lr_bystat=sqrt(apply(MSE_air_5_lr_winter,2,mean,na.rm=T))
mean(RMSE_air_5_lr_bystat)

# 10%
MSE_air_10_lump_winter=(dat_air_winter[,-1]-air_10_lump_winter[,-1])^2
RMSE_air_10_lump_bystat=sqrt(apply(MSE_air_10_lump_winter,2,mean,na.rm=T))
mean(RMSE_air_10_lump_bystat)

MSE_air_10_sarima_winter=(dat_air_winter[,-1]-air_10_sarima_winter_recover[,-1])^2
RMSE_air_10_sarima_bystat=sqrt(apply(MSE_air_10_sarima_winter,2,mean,na.rm=T))
mean(RMSE_air_10_sarima_bystat)

MSE_air_10_naive_winter=(dat_air_winter[,-1]-air_10_naive_winter_recover[,-1])^2
RMSE_air_10_naive_bystat=sqrt(apply(MSE_air_10_naive_winter,2,mean,na.rm=T))
mean(RMSE_air_10_naive_bystat)

MSE_air_10_lr_winter=(dat_air_winter[,-1]-air_10_lr_winter_recover[,-1])^2
RMSE_air_10_lr_bystat=sqrt(apply(MSE_air_10_lr_winter,2,mean,na.rm=T))
mean(RMSE_air_10_lr_bystat)

# 20%
MSE_air_20_lump_winter=(dat_air_winter[,-1]-air_20_lump_winter[,-1])^2
RMSE_air_20_lump_bystat=sqrt(apply(MSE_air_20_lump_winter,2,mean,na.rm=T))
mean(RMSE_air_20_lump_bystat)

MSE_air_20_sarima_winter=(dat_air_winter[,-1]-air_20_sarima_winter_recover[,-1])^2
RMSE_air_20_sarima_bystat=sqrt(apply(MSE_air_20_sarima_winter,2,mean,na.rm=T))
mean(RMSE_air_20_sarima_bystat)

MSE_air_20_naive_winter=(dat_air_winter[,-1]-air_20_naive_winter_recover[,-1])^2
RMSE_air_20_naive_bystat=sqrt(apply(MSE_air_20_naive_winter,2,mean,na.rm=T))
mean(RMSE_air_20_naive_bystat)

MSE_air_20_lr_winter=(dat_air_winter[,-1]-air_20_lr_winter_recover[,-1])^2
RMSE_air_20_lr_bystat=sqrt(apply(MSE_air_20_lr_winter,2,mean,na.rm=T))
mean(RMSE_air_20_lr_bystat)

# Group all RMSEs in one dataframe with columns: method, NA%, RMSE

RMSE_air_5=data.frame(method=rep(c("LUMP","SARIMA","NAIVE","LR"),each=n_stations),
                      RMSE=c(RMSE_air_5_lump_bystat,
                             RMSE_air_5_sarima_bystat,
                             RMSE_air_5_naive_bystat,
                             RMSE_air_5_lr_bystat),
                      NAs=rep(5,n_stations*4),
                      station=rep(colnames(dat_air_winter)[-1],4))

RMSE_air_10=data.frame(method=rep(c("LUMP","SARIMA","NAIVE","LR"),each=n_stations),
                       RMSE=c(RMSE_air_10_lump_bystat,
                              RMSE_air_10_sarima_bystat,
                              RMSE_air_10_naive_bystat,
                              RMSE_air_10_lr_bystat),
                       NAs=rep(10,n_stations*4),
                       station=rep(colnames(dat_air_winter)[-1],4))

RMSE_air_20=data.frame(method=rep(c("LUMP","SARIMA","NAIVE","LR"),each=n_stations),
                       RMSE=c(RMSE_air_20_lump_bystat,
                              RMSE_air_20_sarima_bystat,
                              RMSE_air_20_naive_bystat,
                              RMSE_air_20_lr_bystat),
                       NAs=rep(20,n_stations*4),
                       station=rep(colnames(dat_air_winter)[-1],4))

RMSE_air=rbind(RMSE_air_5,RMSE_air_10,RMSE_air_20)

locations6=locations5[,1:3]
colnames(locations5)=c("station","longitude","latitude")
RMSE_air_full=merge(RMSE_air,locations5,by="station")

# Boxplots by method

bxplt_5=ggplot(RMSE_air_5,aes(x=method,y=RMSE,fill=method))+
  geom_boxplot()+
  ggtitle("5% NAs")+
  theme_classic()

bxplt_10=ggplot(RMSE_air_10,aes(x=method,y=RMSE,fill=method))+
  geom_boxplot()+
  ggtitle("10% NAs")+
  theme_classic()

bxplt_20=ggplot(RMSE_air_20,aes(x=method,y=RMSE,fill=method))+
  geom_boxplot()+
  ggtitle("20% NAs")+
  theme_classic()

pdf("bxplt_RMSE_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
ggarrange(bxplt_5,bxplt_10,bxplt_20,ncol=3,common.legend = T)
dev.off()


# Map of Singapore ----------------------------------------------------------

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# Get the map of Italy
sing <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin == "Singapore")

# Determine bounds for the region of interest
long_min <- min(locations5$longitude) - 0.1
long_max <- max(locations5$longitude) + 0.1
lat_min <- min(locations5$latitude) - 0.1
lat_max <- max(locations5$latitude) + 0.1

# Plot with zoomed-in bounds
ggplot(data = sing) +
  geom_sf(fill = "gray90", color = "black") +
  geom_point(data = locations5, aes(x = longitude, y = latitude), color = "blue", size = 3) +
  geom_text(data = locations5, aes(x = longitude, y = latitude, label = id),
            hjust = 1.5, vjust = 0.5, size = 2.3) +
  coord_sf(xlim = c(long_min, long_max), ylim = c(lat_min, lat_max), expand = FALSE) + # Zoom to bounds
  labs(title = "Map of locations3 (Zoomed)", x = "longitude", y = "latitude") +
  theme_minimal()

# Map RMSE by method or NA% 
plot_by_NA <- ggplot(data = sing) +
  geom_sf(fill = "gray90", color = "black") + # Background map of Italy
  geom_point(data = RMSE_air_full, aes(x = longitude, y = latitude, color = RMSE), 
             size = 3, alpha = 0.7) +
  scale_color_gradient(low = "green", high = "red", name = "RMSE") +
  facet_wrap(~NAs, labeller = labeller(NAs = label_both)) +
  coord_sf(xlim = c(long_min, long_max), ylim = c(lat_min, lat_max), expand = FALSE) + # Zoom to Liguria
  labs(
    #title = "Station Map Grouped by NA% (Color by RMSE)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Frame for subplots
    strip.text = element_text(face = "bold") # Bold text for facet labels
  )

# Plot grouped by Method with color representing RMSE
plot_by_method <- ggplot(data = sing) +
  geom_sf(fill = "gray90", color = "black") + # Background map of Italy
  geom_point(data = RMSE_air_full, aes(x = longitude, y = latitude, color = RMSE), 
             size = 3, alpha = 0.7) +
  scale_color_gradient(low = "green", high = "red", name = "RMSE") +
  facet_wrap(~method, labeller = labeller(method = label_both)) +
  coord_sf(xlim = c(long_min, long_max), ylim = c(lat_min, lat_max), expand = FALSE) + # Zoom to Liguria
  labs(
    #title = "Station Map Grouped by Method (Color by RMSE)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Frame for subplots
    strip.text = element_text(face = "bold") # Bold text for facet labels
  )

# Print the plots
pdf("RMSE_by_NA_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
plot_by_NA
dev.off()

pdf("RMSE_by_method_winter.pdf", width = 11.69, height = 8.27) # A4 landscape dimensions in inches
plot_by_method
dev.off()

