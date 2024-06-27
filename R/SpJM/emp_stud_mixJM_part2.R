load("enth_tab5.Rdata")

# Percentage of NAs for each var
sapply(enth_tab5,function(x){sum(is.na(x))/length(x)})*100

# Summ stat ---------------------------------------------------------------

enth_tab_complete=data.frame(enth_tab5,thermal=gt_thermal)

str(enth_tab_complete)

summary(enth_tab_complete)

Amelia::missmap(enth_tab_complete)

# Compute sd for numeric vars:
enth_tab_complete%>%summarise_if(is.numeric,sd,na.rm=T)
enth_tab_complete%>%summarise_if(is.numeric,mean,na.rm=T)


table(enth_tab_complete$thermal,enth_tab_complete$air_vel)
table(enth_tab_complete$thermal,enth_tab_complete$clothing)
table(enth_tab_complete$thermal,enth_tab_complete$met)
table(enth_tab_complete$thermal,enth_tab_complete$indoor.outdoor)

tapply(enth_tab_complete$heartrate,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$temp,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$skin_temp,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$pressure,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$humidity,enth_tab_complete$thermal,mean,na.rm=T)

tapply(enth_tab_complete$corr_temp_humidity,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$corr_temp_pressure, enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$corr_humidity_nbtemp ,enth_tab_complete$thermal,mean,na.rm=T)


source("Utils.R")

lambda=seq(0,1,.01)
est=lapply(lambda,function(l){
  jump_mixed2(enth_tab5, 2, jump_penalty=l, 
              initial_states=NULL,
              max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
              timeflag=T
  )
})

ARI_res=unlist(lapply(est,function(e){adj.rand.index(e$best_s,gt_thermal)}))
res=data.frame(ARI_res,lambda)
plot(res$lambda,res$ARI_res,type="l",xlab="lambda",ylab="ARI",main="ARI vs lambda")

s_est=matrix(unlist(lapply(est,function(e)e$best_s)),nrow=81,byrow = F)
s_est=apply(s_est,2,order_states)

indx=which(is.na(gt_thermal))
gt_thermal2=gt_thermal[-indx]
gt_thermal2=order_states(gt_thermal2)
s_est=s_est[-indx,]

library(caret)
acc=apply(s_est,2,function(x){confusionMatrix(factor(x,levels=c(1,2)),
                                          factor(gt_thermal2,levels=c(1,2)))$overall[1]}
      )

data.frame(accuracy=acc,lambda)

best_mod=est[[20]]

best_mod$condMM

res_complete=data.frame(enth_tab_complete,state=best_mod$best_s)
tapply(res_complete$temp,res_complete$state,mean,na.rm=T)
table(res_complete$state,res_complete$thermal)/rowSums(table(res_complete$state,res_complete$thermal))

# Kruskal wallis test for numerical variables
sapply(enth_tab5,function(x){kruskal.test(x~best_mod$best_s)$p.value})

# Chisq test for categorical vars
sapply(enth_tab5,function(x){chisq.test(table(x,best_mod$best_s))$p.value})


chisq.test(table(res_complete$air_vel,best_mod$best_s))$p.value
chisq.test(table(res_complete$clothing,best_mod$best_s))$p.value
