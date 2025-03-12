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
#save(blockboot_singapore,file="blockboot_singapore.RData")
load("blockboot_singapore.RData")

apply(blockboot_singapore$t,2,median)
apply(blockboot_singapore$t,2,sd)

## BAC ---------------------------------------------------------------------

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
#save(STJsim_BAC,elapsed_BAC,file="STJsim_BAC.RData")
load("STJsim_BAC.RData")

lambda=seq(0,0.25,by=.05)
gamma=seq(0,0.25,by=.05)

lambdas=unlist(lapply(STJsim_BAC,function(x) x$lambda))
gammas=unlist(lapply(STJsim_BAC,function(x) x$gamma))



BAC_lambdas=matrix(0,nrow=length(unique(lambda))*(length(unique(gamma))-1),ncol=4)

lambda_bench=0
for(i in 1:length(unique(lambdas))){
  indx=which(lambdas==lambda_bench)
  for(j in 1:(length(indx)-1)){
    BAC_lambdas[indx[j],1]=lambda_bench
    BAC_lambdas[indx[j],2]=STJsim_BAC[[indx[j+1]]]$gamma
    BAC_lambdas[indx[j],3]=STJsim_BAC[[indx[j]]]$gamma
    BAC_lambdas[indx[j],4]=
      caret::confusionMatrix(STJsim_BAC[[indx[j]]]$State,STJsim_BAC[[indx[j+1]]]$State)$overall[1]
  }
  lambda_bench=lambda_bench+.05
}

BAC_lambdas=as.data.frame(BAC_lambdas)
colnames(BAC_lambdas)=c("lambda","gamma","gamma_prec","BAC")


BAC_gammas=matrix(NA,nrow=length(unique(gamma))*length(unique(lambda)),ncol=4)

gamma_bench=0
for(i in 1:length(unique(gammas))){
  indx=which(gammas==gamma_bench)
  for(j in 1:(length(indx)-1)){
    BAC_gammas[indx[j],1]=gamma_bench
    BAC_gammas[indx[j],2]=STJsim_BAC[[indx[j+1]]]$lambda
    BAC_gammas[indx[j],3]=STJsim_BAC[[indx[j]]]$lambda
    BAC_gammas[indx[j],4]=
      caret::confusionMatrix(STJsim_BAC[[indx[j]]]$State,STJsim_BAC[[indx[j+1]]]$State)$overall[1]
  }
  gamma_bench=gamma_bench+.05
}

BAC_gammas=BAC_gammas[complete.cases(BAC_gammas),]

BAC_gammas=as.data.frame(BAC_gammas)
colnames(BAC_gammas)=c("gamma","lambda","lambda_prec","BAC")


bac_lam=ggplot(BAC_lambdas, aes(x = gamma, y = BAC, color = as.factor(lambda), group = lambda)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    #title = expression("BAC vs " * gamma * " for Different " * lambda * " Values"),
    x = expression(gamma),
    y = expression(BAC(gamma, gamma[prec])),
    color = expression(lambda)
  ) +
  theme_minimal()+
  theme(legend.position = "bottom")

bac_lam

bac_gamma=ggplot(BAC_gammas, aes(x = lambda, y = BAC, color = as.factor(gamma), group = gamma)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    #title = expression("BAC vs " * gamma * " for Different " * lambda * " Values"),
    x = expression(lambda),
    y = expression(BAC(lambda, lambda[prec])),
    color = expression(gamma)
  ) +
  theme_minimal()+
  theme(legend.position = "bottom")

bac_gamma

ggpubr::ggarrange(bac_lam, bac_gamma, ncol = 2,common.legend = F)
