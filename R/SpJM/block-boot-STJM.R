library(tidyr)
library(dplyr)
source("Utils.R")


# block-boot --------------------------------------------------------------

# Converting to wide format
Y_6_wide <- Y_6 %>%
  pivot_wider(names_from = m, 
              values_from = c(air_temp, rh, rainfall, wind_speed, windy, hour, 
                              air_temp_rollmean, rh_rollmean, rainfall_rollmean, 
                              wind_speed_rollmean, air_temp_rollsd, rh_rollsd, 
                              rainfall_rollsd, wind_speed_rollsd, UTCI))

Y=Y_6_wide 

STJM_blockboot=function(Y,K=3,D,lambda,gamma){
  
  Y$t=1:nrow(Y)
  # Y is in wide format
  Y <- Y %>%
    mutate(across(-t, as.numeric))
  
  # Convert back to long format
  Y <- Y %>%
    pivot_longer(cols = -t, 
                 names_to = c(".value", "m"), 
                 names_pattern = "(.+)_(\\d+)") %>%
    mutate(m = as.integer(m))
  
  Y$windy=as.factor(Y$windy)
  Y$hour=as.factor(Y$hour)
  
  fit=STjumpDist(Y,n_states=K,D,
                 jump_penalty=lambda,
                 spatial_penalty=gamma,
                 initial_states=NULL,
                 max_iter=10, n_init=10, 
                 tol=1e-4, 
                 verbose=F,timeflag=F)
  
  M=length(unique(fit$Y$m))
  State=c(t(fit$best_s))
  State=order_states_condMean(fit$Y$air_temp,State)
  
  S_est=matrix(State,ncol=M,byrow = T)
  
  # Y is in wide format
  Y <- Y %>%
    mutate(across(-t, as.numeric))
  
  
  Y$windy=as.factor(Y$windy)
  Y$hour=as.factor(Y$hour)
  #
  Y=Y[order(Y$t,Y$m),]
  
  
  Y_res=data.frame(Y,State=State)
  mumo=matrix(0,3,7)
  
  mumo=data.frame(
    air_temp= tapply(Y_res$air_temp, Y_res$State, median, na.rm = TRUE),
    rh = tapply(Y_res$rh, Y_res$State, median, na.rm = TRUE),
    rainfall = tapply(Y_res$rainfall, Y_res$State, median, na.rm = TRUE),
    wind = tapply(Y_res$wind_speed, Y_res$State, median, na.rm = TRUE),
    windy = tapply(Y_res$windy, Y_res$State, Mode),
    hour = tapply(Y_res$hour, Y_res$State, Mode) - 1,
    UTCI = tapply(Y_res$UTCI, Y_res$State, median, na.rm = TRUE)
  )
  
  return(unlist(mumo))
}

lambda=.05
gamma=.05
D=distm(data_stat_number[,c("longitude","latitude")], 
        data_stat_number[,c("longitude","latitude")], 
        fun = distGeo)/1000

#temp=STJM_blockboot(Y,K=3,D,lambda,gamma)

TT=length(unique(Y_6$t))
l1=round(TT^(2/3))
nboot=1000
#nboot=10
K=3
library(boot)
blockboot_singapore=tsboot(Y_6_wide, STJM_blockboot, R = nboot, l = l1, sim = "fixed",
             parallel = "multicore", ncpus = detectCores()-1, n.sim=TT,
             lambda=lambda,gamma=gamma,D=D,K=K)
save(blockboot_singapore,file="blockboot_singapore.RData")


# BAC ---------------------------------------------------------------------

fit_ordered=function(Y,K,D,gamma,lambda){
  fit=STjumpDist(Y,3,D,
                 jump_penalty=lambda,
                 spatial_penalty=gamma,
                 initial_states=NULL,
                 max_iter=10, n_init=10, 
                 tol=1e-4, 
                 verbose=F,timeflag=T)
  State=c(t(fit$best_s))
  State=order_states_condMean(fit$Y$air_temp,State)
  State=factor(State,levels=1:K)
  return(list(State=State,gamma=gamma,lambda=lambda))
}

lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)
hp=expand.grid(lambda=lambda,gamma=gamma)

start_BAC=Sys.time()
STJsim_BAC <- parallel::mclapply(1:nrow(hp),
                                function(x)
                                  fit_ordered(Y=Y_6,K=3,D=D,
                                              lambda=hp[x,]$lambda,
                                              gamma=hp[x,]$gamma),
                                mc.cores = parallel::detectCores())
end_BAC=Sys.time()
elapsed_BAC=end_BAC-start_BAC
save(STJsim_BAC,elapsed_BAC,file="STJsim_BAC.RData")
