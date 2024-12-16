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
                                        target=target2,
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


temp=conv_analysis_jump_mixed(SimData,n_states=3, jump_penalty=.1,
                              max_iter=50, n_init=30,tol=NULL,verbose = T)
temp$data_convergence
temp$id_best_init
temp$id_best_iter
