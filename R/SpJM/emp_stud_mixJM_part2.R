load("enth_tab5.Rdata")


# Summ stat ---------------------------------------------------------------

enth_tab_complete=data.frame(enth_tab5,thermal=gt_thermal)

str(enth_tab_complete)

table(enth_tab_complete$thermal,enth_tab_complete$air_vel)
table(enth_tab_complete$thermal,enth_tab_complete$clothing)
table(enth_tab_complete$thermal,enth_tab_complete$met)
table(enth_tab_complete$thermal,enth_tab_complete$indoor.outdoor)

tapply(enth_tab_complete$heartrate,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$temp,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$skin_temp,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$pressure,enth_tab_complete$thermal,mean,na.rm=T)
tapply(enth_tab_complete$humidity,enth_tab_complete$thermal,mean,na.rm=T)

source("Utils.R")

lambda=seq(0,1,.01)
est=lapply(lambda,function(l){
  jump_mixed2(enth_tab_complete, 2, jump_penalty=l, 
              initial_states=NULL,
              max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
              timeflag=T
  )
})

# ARI_res=unlist(lapply(est,function(e){adj.rand.index(e$best_s,gt_thermal)}))
# res=data.frame(ARI_res,lambda)
# plot(res$lambda,res$ARI_res,type="l",xlab="lambda",ylab="ARI",main="ARI vs lambda")

res[which.max(res$ARI_res),]
best_est=est[[which.max(res$ARI_res)]]

best_est$best_s
best_est$condMM
