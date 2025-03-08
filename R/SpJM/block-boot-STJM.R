STjumpDist_wide=function(Y,n_states,
                    D,
                    jump_penalty=1e-5,
                    spatial_penalty=1e-5,
                    initial_states=NULL,
                    max_iter=10, n_init=10, tol=1e-4, verbose=FALSE,timeflag=T){
  
  
  # Y is in wide format
  Y <- Y %>%
    mutate(across(-t, as.numeric))
  
  # Convert back to long format
  Y_long <- Y %>%
    pivot_longer(cols = -t, 
                 names_to = c(".value", "m"), 
                 names_pattern = "(.+)_(\\d+)") %>%
    mutate(m = as.integer(m))
  
  Y=Y_long
  
  Y$windy=as.factor(Y_long$windy)
  Y$hour=as.factor(Y_long$hour)
  #
  Y=Y[order(Y$t,Y$m),]
  
  P=ncol(Y)-2
  # Time differences
  Y.orig=Y
  
  if(timeflag){
    time=sort(unique(Y.orig$t))
    dtime=diff(time)
    dtime=dtime/as.numeric(min(dtime))
    dtime=as.numeric(dtime)
    Y=subset(Y,select=-c(t,m))
    Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  }
  
  else{
    Y=subset(Y,select=-c(t,m))
    Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  }
  
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  library(dplyr)
  Y <- data.frame(t=Y.orig$t,m=Y.orig$m,Y)
  
  # Reorder columns so that we have t,m, cat vars and cont vars
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  cont.indx=cont.indx[-(1:2)]
  Y=Y[,c(1,2,cat.indx,cont.indx)]
  
  YY=subset(Y,select=-c(t,m))
  
  TT=length(unique(Y$t))
  M=length(unique(Y$m))
  
  cat.indx=which(sapply(YY, is.factor))
  cont.indx=which(sapply(YY, is.numeric))
  
  Ycont=YY[,cont.indx]
  
  Ycat=YY[,cat.indx]
  levels_cat=lapply(Ycat,levels)
  names(levels_cat)=cat.indx
  
  n_cat=length(cat.indx)
  n_cont=length(cont.indx)
  
  ###
  # Missing data imputation TBD
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  # Initialize mu 
  #mu <- colMeans(Ycont,na.rm = T)
  mu=apply(Ycont,2,median,na.rm=T)
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  # Impute missing values with mean of observed states
  
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    x=Ycat[,i]
    Ycat[which(is.na(Ycat[,i])),i]=mo[i]
  }
  YY[,-cat.indx]=Ycont
  YY[,cat.indx]=Ycat
  
  Y[,-(1:2)]=YY
  
  ###
  
  # State initialization through kmeans++
  S=matrix(0,nrow=TT,ncol=M)
  for(m in 1:M){
    S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    
    for (it in 1:max_iter) {
      for (i in unique(as.vector(S))) {
        #mu[i,] <- colMeans(Ycont[as.vector(t(S))==i,])
        mu[i,]=apply(Ycont[as.vector(t(S))==i,],2,median)
        mo[i,]=apply(Ycat[as.vector(t(S))==i,],2,Mode)
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
      }
      mumo=data.frame(matrix(0,nrow=n_states,ncol=P))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(YY)
      
      
      # Fit state sequence
      S_old <- S
      
      loss_v=NULL
      
      # Re-fill-in missings 
      
      Ycont <- Mcont * mu[as.vector(t(S)), ] + (!Mcont) * Ycont
      
      mo_2=apply(mo,2,as.numeric)
      Ycat <- Mcat * mo_2[as.vector(t(S)), ] + (!Mcat) * apply(Ycat,2,as.numeric)
      Ycat=data.frame(Ycat)
      for(cc in 1:n_cat){
        Ycat[,cc]=factor(Ycat[,cc],levels=levels_cat[[cc]])
      }
      
      # for(i in 1:nrow(Mcont)){
      #   Ycont[i,]=unlist(ifelse(Mcont[i,],mu[as.vector(t(S))[i],],Ycont[i,]))
      # }
      # for(i in 1:nrow(Mcat)){
      #   Ycat[i,]=unlist(ifelse(Mcat[i,],mo[as.vector(t(S))[i],],Ycat[i,]))
      # }
      
      YY[,-cat.indx]=Ycont
      YY[,cat.indx]=Ycat
      
      Y[,-(1:2)]=YY
      
      ###
      for(m in 1:M){
        loss_by_state=gower.dist(Y[which(Y$m==m),-(1:2)],mumo)
        
        for(k in 1:n_states){
          #temp <- t(t((S[,-m]==k))/D[m,-m])
          temp <- t(t((S[,-m]==k))*exp(-D[m,-m]))
          loss_by_state[,k]=loss_by_state[,k]-spatial_penalty*rowSums(temp)
        }
        
        V <- loss_by_state
        for (t in (TT-1):1) {
          #V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
          if(timeflag){
            V[t-1,] <- loss_by_state[t-1,] + apply(V[t,]/dtime[t] + Gamma, 2, min)
          }
          else{
            V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
          }
        }
        
        S[1,m] <- which.min(V[1,])
        for (t in 2:TT) {
          S[t,m] <- which.min(V[t,] + Gamma[S[t-1,m],])
        }
        loss_v=c(loss_v,min(V[1,]))
      }
      if (length(unique(S)) == 1) {
        break
      }
      loss <- mean(loss_v)
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
      }
      if (!is.null(tol)) {
        epsilon <- loss_old - loss
        if (epsilon < tol) {
          break
        }
      } else if (all(S == S_old)) {
        break
      }
      loss_old <- loss
      
      ###
      
    }
    if (is.null(best_s) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_s <- S
    }
    
    for(m in 1:M){
      S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],n_states)
    }
  }
  return(list(best_s=best_s,
              Y=Y,
              K=n_states,
              lambda=jump_penalty,
              gamma=spatial_penalty,
              loss=best_loss))
  
}

# Y_6 from emp_stud_STJM_3
library(tidyr)
library(dplyr)

# Converting to wide format
Y_6_wide <- Y_6 %>%
  pivot_wider(names_from = m, values_from = c(air_temp, rh, rainfall, wind_speed, windy, hour, 
                                              air_temp_rollmean, rh_rollmean, rainfall_rollmean, 
                                              wind_speed_rollmean, air_temp_rollsd, rh_rollsd, 
                                              rainfall_rollsd, wind_speed_rollsd, UTCI))

Y=Y_6_wide 
STJM_blockboot=function(Y,K=3,D,lambda,gamma){
  
  fit=STjumpDist_wide(Y,n_states=K,D,
                 jump_penalty=lambda,
                 spatial_penalty=gamma,
                 initial_states=NULL,
                 max_iter=10, n_init=10, 
                 tol=1e-4, 
                 verbose=F,timeflag=T)
  
  M=length(unique(Y_long$m))
  State=c(t(fit$best_s))
  State=order_states_condMean(Y_long$air_temp,State)
  
  S_est=matrix(State,ncol=M,byrow = T)
  
  Y_res=data.frame(Y_long,State=State)
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
  
  return(mumo)
}


TT=length(unique(Y_6$t))
l1=round(TT^(2/3))
#nboot=1000
nboot=2

library(boot)
prova=tsboot(Y_6_wide, STJM_blockboot, R = nboot, l = l1, sim = "geom",
       parallel = "multicore", ncpus = detectCores()-1, n.sim=TT,
       lambda=lambda,gamma=gamma,D=D,K=K)

