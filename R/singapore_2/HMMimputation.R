library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
require(LMest)
require(mvtnorm)
require(MultiLCIRT)

source('HMMcont_miss/RCode/lmbasic.cont.MISS.R')
source('HMMcont_miss/RCode/lk_comp_cont_MISS.R')
source('HMMcont_miss/RCode/lmcovlatent.cont.MISS.R')
source('HMMcont_miss/RCode/lk_comp_latent_cont_MISS.R')
source('HMMcont_miss/RCode/prob_post_cov_cont.R')
source('HMMcont_miss/RCode/est_multilogit.R')
source('HMMcont_miss/RCode/prob_multilogit.R')
source("HMMcont_miss/RCode/bootstrap.MISS.R")


load("enths.Rdata")

# Construct objects
ui=unique(enth_out$user_id)
ut=sort(unique(enth_out$time))

n=length(unique(ui))
TT=length(unique(ut))

temp=data.frame(time=rep(ut,each=n),user_id=rep(unique(ui),length(ut)))

# Join temp and enth_out
temp=merge(temp,enth_out,by=c("time","user_id"),all=T)

# Construct response vars object
Yweather=temp[,c("user_id","time","temp_outdoor","humidity_outdoor","pressure_outdoor")]
Yphysio=temp[,c("user_id","time","heartrate","resting_heartrate","nb_temp","skin_temp")]
Yaq=temp[,c("user_id","time","pm1.0_outdoor","pm2.5_outdoor","pm10.0_outdoor")]

YYweather=array(0,dim=c(n,TT,3))
YYphysio=array(0,dim=c(n,TT,4))
YYaq=array(0,dim=c(n,TT,3))

for(i in 1:n){
  YYweather[i,,1]=log(temp[temp$user_id==ui[i],"temp_outdoor"])
  YYweather[i,,2]=log(temp[temp$user_id==ui[i],"humidity_outdoor"])
  YYweather[i,,3]=log(temp[temp$user_id==ui[i],"pressure_outdoor"])
  YYphysio[i,,1]=log(temp[temp$user_id==ui[i],"heartrate"])
  YYphysio[i,,2]=log(temp[temp$user_id==ui[i],"resting_heartrate"])
  YYphysio[i,,3]=log(temp[temp$user_id==ui[i],"nb_temp"])
  YYphysio[i,,4]=log(temp[temp$user_id==ui[i],"skin_temp"])
  YYaq[i,,1]=temp[temp$user_id==ui[i],"pm1.0_outdoor"]
  YYaq[i,,2]=temp[temp$user_id==ui[i],"pm2.5_outdoor"]
  YYaq[i,,3]=temp[temp$user_id==ui[i],"pm10.0_outdoor"]
}


nrep <- 10 
Kmax <- 4
tol <- 10^-4

modva = vector("list",Kmax)
modvb = vector("list",Kmax)
modvc = vector("list",Kmax)

for(k in 1:Kmax){
  print(k)
  if(k>=7) nrep<-2
  modva[[k]] <- lmbasic.cont.MISS(YYweather,k=k,modBasic=1,start=0,tol=tol)
  modvb[[k]] <- lmbasic.cont.MISS(YYphysio,k=k,modBasic=1,start=0,tol=tol)
  modvc[[k]] <- lmbasic.cont.MISS(YYaq,k=k,modBasic=1,start=0,tol=tol)
  
  if(k>1){
    for(k1 in 1:(nrep*(k-1))){
      print(c(k,k1))
      tmpa <- lmbasic.cont.MISS(YYweather,k=k,modBasic=1, start=1,tol=tol)
      if(tmpa$lk>modva[[k]]$lk){
        modva[[k]] = tmpa
      }
      tmpb <- lmbasic.cont.MISS(YYphysio,k=k,modBasic=1, start=1,tol=tol)
      if(tmpb$lk>modvb[[k]]$lk){
        modvb[[k]] = tmpb
      }
      tmpc <- lmbasic.cont.MISS(YYaq,k=k,modBasic=1, start=1,tol=tol)
      if(tmpc$lk>modvb[[k]]$lk){
        modvc[[k]] = tmpc
      }
    }
  }
}

lapply(modva, function(x) x$bic)
lapply(modvb, function(x) x$bic)
lapply(modvc, function(x) x$bic)

dim(modva[[4]]$Yimp)
plot(exp(modva[[4]]$Yimp[16,,3]),type="l")
  plot(exp(modvb[[4]]$Yimp[16,,4]),type="l")
