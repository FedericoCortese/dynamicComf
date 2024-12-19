source("Utils.R")

conv_analysis_jump_mixed <- function(Y, n_states, 
                                     jump_penalty=1e-5, 
                       initial_states=NULL,
                       max_iter=50, n_init=10, tol=1e-5, verbose=FALSE) {
  
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Which vars are categorical and which are numeric
  
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  Ycont=Y[,-cat.indx]
  Ycat=Y[,cat.indx]
  
  n_levs=apply(Ycat, 2, function(x)length(unique(x)))
  # n_levs=apply(Ycat, 2, function(x)levels(x))
  
  
  n_cat=length(cat.indx)
  n_cont=n_features-n_cat
  
  # Initialize mu 
  # mu <- colMeans(Ycont,na.rm = T)
  mu <- apply(Ycont, 2, median, na.rm = TRUE)
  
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  
  #M=ifelse(is.na(Y),T,F)
  
  
  Ytil=Y
  # Impute missing values with median of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    Ycat[,i]=ifelse(Mcat[,i],mo[i],Ycat[,i])
    Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  }
  
  # Ycat=Ycat%>%mutate_all(factor)
  # for(i in 1:n_cat){
  #   Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  # }
  
  Y[,-cat.indx]=Ycont
  Y[,cat.indx]=Ycat
  
  # State initialization through kmeans++
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s=initialize_states(Y,n_states)
  }
  
  data_convergence=NULL
  id_best_iter=NULL
  id_best_init=NULL
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features-length(cat.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        # Fit model by updating prototypes of observed states
        #if(sum(s==i)>1){
        #mu[i,] <- colMeans(Ycont[s==i,])
        mu[i,] <- apply(Ycont[s==i,], 2, median, na.rm = TRUE)
        mo[i,]=apply(Ycat[s==i,],2,Mode)
        # }
        # else{
        #   mu[i,]=mean(Y[s==i,])
        # }
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        mo[,i]=factor(mo[,i],levels=1:n_levs[i])
      }
      
      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:ncol(Ycont)){
        Ycont[,i]=ifelse(Mcont[,i],mu[s,i],Ycont[,i])
      }
      for(i in 1:ncol(Ycat)){
        Ycat[,i]=ifelse(Mcat[,i],mo[s,i],Ycat[,i])
        Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
      }
      
      Y[,-cat.indx]=Ycont
      Y[,cat.indx]=Ycat
      
      mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      
      # loss_by_state=matrix(0,nrow=n_obs,ncol=n_states)
      # for(k in 1:n_states){
      #   loss_by_state[,k]=gower.dist(Y,mumo[k,])
      # }
      
      
      # var.weights in gower.dist allows for weighted distance
      
      loss_by_state=gower.dist(Y,mumo)
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
      }
      
      s[1] <- which.min(V[1,])
      for (t in 2:n_obs) {
        s[t] <- which.min(V[t,] + Gamma[s[t-1],])
      }
      
      target1=sum(loss_by_state[cbind(1:nrow(loss_by_state), 
                                      s[1:nrow(loss_by_state)])])
      
      n_jumps=sum(s[1:(n_obs-1)]!=s[2:n_obs])
      
      obs_mov=sum(s_old!=s)
      
      target2=target1+jump_penalty*n_jumps
      data_convergence=rbind(data_convergence,
                             data.frame(n_init=init,
                                        max_iter=it,
                                        targetJM=target2,
                                        target_kprot=target1,
                                        obs_mov=obs_mov))
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- min(V[1,])
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
      }
      if (!is.null(tol)) {
        epsilon <- loss_old - loss
        if (epsilon < tol) {
          break
        }
      } 
      # Possibly comment this to see evolution of target fun
      else if (all(s == s_old)) {
        break
      }
      loss_old <- loss
    }
    
    if (is.null(best_s) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_s <- s
      id_best_iter=it
      id_best_init=init
    }
    #s <- init_states(Y, n_states)+1
    s=initialize_states(Y,n_states)
  }
  
  colnames(mumo)=colnames(Y)
  
  # Add target function value
  
  return(list(best_s=best_s,
              Y=Y,
              Y.orig=Ytil,
              condMM=mumo,
              target1=target1,
              target2=target2,
              data_convergence=data_convergence,
              id_best_init=id_best_init,
              id_best_iter=id_best_iter))
}

SimData=sim_data_mixed(TT=200,P=20)
Y=SimData$SimData.complete

temp=conv_analysis_jump_mixed(Y,n_states=3, jump_penalty=.1,
                              max_iter=50, n_init=30,tol=NULL,verbose = T)
temp$data_convergence
temp$id_best_init
temp$id_best_iter

res=temp$data_convergence

sz=22
pconv1 = res %>%
  group_by(n_init) %>%
  ggplot(aes(x = max_iter, y = targetJM, color = factor(n_init))) +
  geom_line(linewidth = .8) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = sz), # Size of x-axis label
    axis.title.y = element_text(size = sz), # Size of y-axis label
    axis.text.x = element_text(size = sz-2),  # Size of x-axis tick text
    axis.text.y = element_text(size = sz-2)   # Size of y-axis tick text
  ) +
  labs(x = "Iteration", y = "Objective function")

pconv2 = res %>%
  group_by(n_init) %>%
  ggplot(aes(x = max_iter, y = obs_mov, color = factor(n_init))) +
  geom_line(linewidth = .8) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = sz), # Size of x-axis label
    axis.title.y = element_text(size = sz), # Size of y-axis label
    axis.text.x = element_text(size = sz-2),  # Size of x-axis tick text
    axis.text.y = element_text(size = sz-2)   # Size of y-axis tick text
  ) +
  labs(x = "Iteration", y = "# Observations moved")


# Merge
library(gridExtra)
grid.arrange(pconv1,pconv2,ncol=2)


res%>%group_by(n_init)%>%
  ggplot(aes(x=max_iter,y=target_kprot,color=factor(n_init)))+
  geom_line(linewidth=.8)+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x = "Iteration", y = "Objective function")+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = sz), # Size of x-axis label
    axis.title.y = element_text(size = sz), # Size of y-axis label
    axis.text.x = element_text(size = sz-2),  # Size of x-axis tick text
    axis.text.y = element_text(size = sz-2)   # Size of y-axis tick text
  ) +
  labs(x = "Iteration", y = "Objective function")



# Run 100 times for ninit=1,2,5,10,50, save mean and SD of results only

n_init=c(1,2,5,10,20)
seeds=1:100

conv_fun_sd=function(seed,n_init,
                     TT=200,
                     P=20,
                     Ktrue=3,
                     mu=1,
                     phi=.8,
                     rho=0,
                     Pcat=10,
                     pers=.95,
                     pNAs=0,
                     typeNA=3){
  simDat=sim_data_mixed(seed=seed,
                        TT=TT,
                        P=P,
                        Ktrue=Ktrue,
                        mu=mu,
                        phi=phi,
                        rho=rho,
                        Pcat=Pcat,
                        pers=pers,
                        pNAs=pNAs,
                        typeNA=typeNA)
  # Estimate
  est=jump_mixed(simDat$SimData.NA,
                 n_states=Ktrue,
                 jump_penalty = .1,
                 n_init=n_init,
                 verbose=F)
  target1=est$target1
  target2=est$target2
  return(list(target1=target1,target2=target2))
}

hp_conv=expand.grid(n_init=n_init,seed=seeds)

temp2 <- apply(hp_conv, 1, function(x) conv_fun_sd(as.numeric(x["seed"]), as.numeric(x["n_init"])))

res_sd=data.frame(hp_conv,
                  target1=unlist(lapply(temp2,function(x)x$target1)),
                  target2=unlist(lapply(temp2,function(x)x$target2)))


head(res_sd)

library(ggplot2)

# Melt the data if you want separate boxplots for `target1` and `target2`
library(reshape2)
res_sd_melted <- melt(res_sd, id.vars = c("n_init", "seed"), measure.vars = c("target1", "target2"))

# Create the boxplot
ggplot(res_sd_melted, aes(x = factor(n_init), y = value, fill = variable)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "n_init", y = "Target Value", fill = "Target Type") +
  theme(legend.position = "top")

library(dplyr)
library(ggplot2)

# Calculate standard deviation of target1 for each n_init
sd_summary <- res_sd %>%
  group_by(n_init) %>%
  summarise(sd_target1 = sd(target1, na.rm = TRUE))

# Create the line plot
ggplot(sd_summary, aes(x = n_init, y = sd_target1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "n_init", y = "SD of Target1", title = "Standard Deviation of Target1 by n_init")

