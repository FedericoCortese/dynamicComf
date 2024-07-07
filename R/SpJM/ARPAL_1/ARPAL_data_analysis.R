load("dat.Rdata")
load("AQI_fact.Rdata")

# EDA ---------------------------------------------------------------------

str(dat)

# Percentage of missing values
round(colMeans(is.na(dat)) * 100, 2)

Amelia::missmap(dat)

# Summary statistics
summary(dat)

# Correlation plot
windows()
corrplot::corrplot(cor(dat[,2:13],use="complete.obs"),method="number")

round(cor(dat[,2:13],use="complete.obs"),2)

tapply(dat$PM2.5_mean, dat$rainy, mean)
tapply(dat$PM2.5_mean, dat$weekend,mean)
tapply(dat$PM2.5_mean, dat$weekday,mean)

tapply(dat$NO2_mean, dat$rainy, mean)
tapply(dat$NO2_mean, dat$weekend,mean)
tapply(dat$NO2_mean, dat$weekday,mean)

tapply(dat$Ozone_max_day, dat$rainy, mean)
tapply(dat$Ozone_max_day, dat$weekend,mean)
tapply(dat$Ozone_max_day, dat$weekday,mean)



# JM mix ------------------------------------------------------------------

# K=4 as per observed AQI levels
source("Utils.R")
lambda=seq(0,1,.01)
dat_notime=dat[,-1]
est=lapply(lambda,function(l){
  jump_mixed2(dat_notime, 4, jump_penalty=l, 
              initial_states=NULL,
              max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
              timeflag=F
  )
})

GIC_mixed(dat_notime,est[[2]]$best_s,est[[1]]$best_s,K=4)

ARI_res=unlist(lapply(est,function(e){adj.rand.index(e$best_s,AQI_fact)}))
GIC=unlist(lapply(est,function(e){GIC_mixed(e$Y,e$best_s,est[[1]]$best_s,K=4)$FTIC}))
res=data.frame(ARI_res,GIC,lambda)

plot(res$lambda,res$ARI_res,type="l",xlab="lambda",ylab="ARI",main="ARI vs lambda")

best_est=est[[29]]

table(best_est$best_s)

best_est$condMM

states=factor(best_est$best_s,levels=1:4,labels=c("Good","Moderate","US","Unhealthy"))
true_states=factor(AQI_fact,levels=1:4,labels=c("Good","Moderate","US","Unhealthy"))


library(caret)
confusionMatrix(states,true_states)
