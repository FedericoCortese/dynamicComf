library(dplyr)
library(lubridate)
library(gstat)
library(sp)

load("cleaned_data.RData")
wd = air_temp
wname = "air_temp"

#source("weather_summ.R") # default media ogni 15 minuti
source("Utils2.R")
tmp = weather_summ(wd,window="1 hour")
wd_summ = tmp$mean_x
n_summ = data.frame(tmp$n_x) # numero di valori validi per la media in ciascun intervallo di tempo
rm(tmp)

init=which(wd_summ$time=="2023-01-22 00:00:00")
fin=which(wd_summ$time=="2023-03-22 00:00:00")

subs=init:fin
# riduzione al sottoinsieme
wd_summ_sub = wd_summ[subs,]
n_summ_sub = n_summ[subs,]

cols=which(apply(wd_summ_sub,2,function(x) sum(is.na(x)))==0)

wd_summ_sub = wd_summ_sub[,cols]

TT=dim(wd_summ_sub)[1]
n_stats=dim(wd_summ_sub)[2]-1

# Fill with NA
air_short=wd_summ_sub
TT=dim(air_short)[1]
na_len=TT*c(.05,.1,.2)
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


windows()
par(mfrow=c(2,3),mar=c(2,2,6,2))
for(i in 2:ncol(air_20)){
  plot(x=air_20$time,y=as.vector(unlist(air_20[,i])),type="l",col="grey20",
       xlab=" ",ylab=" ",
       main=colnames(air_20[,i]))
  title(main=colnames(air_20)[i])
}
mtext("Air temperatures - 20% NAs", side = 3, line = - 2, outer = TRUE)

data=air_20

# Keep out one station
id_stat=2

CV_lump.kr=function(id_stat,data,locations){
  
  start=Sys.time()
  
  test=data[,c(1,id_stat)]
  # indexes on which I have to predict
  ttindex=which(is.na(test[,id_stat]))
  nsubs = length(ttindex)
  subs=1:nsubs
  
  train=data[,-id_stat]
  
  y = matrix(NA,length(ttindex),2)
  for(ind in unique(ttindex)){
    #print(ind)
    # finestra mobile appena possibile
    if(ind>nsubs) subs = (ind-nsubs+1):ind
    # riduzione al sottoinsieme
    train_sub = train[subs,]
    #n_summ_sub = n_summ[subs,]
    
    
    locations_wd = locations
    coordinates(locations_wd) = ~ longitude + latitude
    locations_wd = locations_wd[match(names(train)[-1], locations_wd$id),]
    proj4string(locations_wd) = "+proj=longlat +datum=WGS84"
    row.names(locations_wd) = locations_wd$id
    # creo oggetto ST
    space = list(value = names(train_sub)[-1])
    train_ST = stConstruct (train_sub[,-1], 
                            space, 
                            train_sub$time, 
                            SpatialObj= locations_wd)
    
    spdf.lst <- lapply(1:nsubs,
                       function(i) {
                         x = train_ST[,i]
                         x$ti = i # ti = time index
                         x$id = train_ST@sp$id
                         return(x)}
    )
    
    spdf <- do.call(rbind, spdf.lst)
    # summary(spdf) # Object of class SpatialPointsDataFrame
    
    if(ind==ttindex[1] & ind<=nsubs){
      # cutoff basato su esame di dists = spDists(locations_wd,longlat=T)
      # e visualizzazione di alcuni variogrammi
      vl <- variogram(values ~ ti, spdf[!is.na(spdf$values),], cutoff=25, dX=0)
      isill = quantile(vl$gamma,0.95)
      ihrange = vl$dist[match(T,vl$gamma>=isill)]
      vlm <- fit.variogram(vl, model=vgm(isill, "Exp", ihrange/3),fit.method=6) # non assegnando nugget nel modello viene stimato zero
      # plot(vl, plot.numbers=T, model=vlm, main=paste(wname, " ", nsub, " 15-min intervals lumped: ", len, "days"))
    }
    if(ind>nsubs){ # devo rifare la stima quando la finestra inizia a muoversi
      vl <- variogram(values ~ ti, spdf[!is.na(spdf$values),], cutoff=25, dX=0)
      isill = quantile(vl$gamma,0.95)
      ihrange = vl$dist[match(T,vl$gamma>=isill)]
      vlm <- fit.variogram(vl, model=vgm(isill, "Exp", ihrange/3),fit.method=6) # non assegnando nugget nel modello viene stimato zero
    }
    
    # Leave-one(station)-out
    slocs=locations[which(locations$id==colnames(data)[id_stat]),c("latitude","longitude")]
    slocs <- SpatialPoints(slocs,
                           proj4string=CRS(proj4string(spdf)))

    # kriging
    if(ind <=nsubs) 
      k.lump = krige(values ~ 1,spdf[!is.na(spdf$values) & spdf$ti==ind,], newdata=slocs,model=vlm)
    else
      k.lump = krige(values ~ 1,
                     spdf[!is.na(spdf$values) & spdf$ti==nsubs,], 
                     newdata=slocs,model=vlm)
    
    y[which(ttindex==ind),] = as.matrix(k.lump@data)
  }
  end=Sys.time()
  elapsed=end-start
  
  return(list(y=y,
              elapsed=elapsed,
              id_stat=id_stat,
              stat=colnames(data)[id_stat]))
}

plot_indx=which(is.na(data$S43))
res=data.frame(time=air_short$time[plot_indx],
               true=air_short$S43[plot_indx],
               pred=y[,1])

plot(res$true,type='l')
lines(res$pred,col='red')

# Run in parallel
indx=2:ncol(air_short)
air5_lump <- parallel::mclapply(indx,
                                function(x)CV_lump.kr(x,air_5,locations),
                                mc.cores = parallel::detectCores())
air10_lump <- parallel::mclapply(indx,
                                function(x)CV_lump.kr(x,air_10,locations),
                                mc.cores = parallel::detectCores())
air20_lump <- parallel::mclapply(indx,
                                function(x)CV_lump.kr(x,air_20,locations),
                                mc.cores = parallel::detectCores())

