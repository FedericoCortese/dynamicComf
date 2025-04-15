load("best_est_aq.Rdata")

best_est$K
Y=best_est$Y
Y=Y[,-20]
str(Y)

# fuzzy JM

source("Utils_fuzzyJM.R")



fit=fuzzy_jump(Y,n_states=4,jump_penalty = .1,verbose=T)
fit$MAP
x11()
plot(Y$pm10,col=fit$MAP)
legend("topright",legend=1:4,col=1:4,lty=1)
tapply(Y$pm25,fit$MAP,mean)
tapply(Y$o3,fit$MAP,mean)
tapply(Y$pm10,fit$MAP,mean)
tapply(Y$temp,fit$MAP,mean)
tapply(Y$rel_hum,fit$MAP,mean)
tapply(Y$rainfall,fit$MAP,mean)
tapply(Y$wind_speed,fit$MAP,mean)

x11()
par(mfrow = c(2,2))

for (i in 1:4) {
  par(new = FALSE)
  plot(Y$pm10, type = "l", col = "grey", 
       ylab = "PM10", xlab = "Time", main = paste("Component", i))
  par(new = TRUE)
  plot(fit$best_S[, i], type = "l", col = i, axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Probability", side = 4, line = 3)
}



