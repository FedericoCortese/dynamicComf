

onerun_cont_STJM=function(Y,K,
                       jump_penalty,
                       spatial_penalty,
                       alpha,grid_size,mode_loss=T,
                       max_iter,tol,initial_states){
  tryCatch({
    
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
    # Missing data imputation 
    Mcont=ifelse(is.na(Ycont),T,F)
    Mcat=ifelse(is.na(Ycat),T,F)
    mu <- colMeans(Ycont,na.rm = T)
    mo <- apply(Ycat,2,Mode)
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
    
    #TT=nrow(Y)
    TT=length(unique(Y$t))
    loss_old <- 1e10
    

    # State initialization
    S=matrix(0,nrow=TT,ncol=M)
    for(m in 1:M){
      S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],n_states)
    }
    
    # Initialize soft clustering matrix
    SS <- matrix(0, nrow = TT * M, ncol = 2 + n_states)
    SS[, 1] <- rep(1:TT, each = M)  # First column: t indices
    SS[, 2] <- rep(1:M, times = TT) # Second column: m indices
    for (t in 1:TT) {
      for (m in 1:M) {
        state <- S[t, m]
        SS[(t - 1) * M + m, 2 + state] <- 1  # Set the appropriate column to 1
      }
    }
    SS=data.frame(SS)
    colnames(SS)=c("t","m",1:n_states)
    
    #mu=matrix(NA,nrow=K,ncol=P)
    mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    
    for (i in unique(as.vector(S))) {
      mu[i,] <- colMeans(Ycont[as.vector(t(S))==i,])
      mo[i,]=apply(Ycat[as.vector(t(S))==i,],2,Mode)
    }
    
    for (it in 1:max_iter) {
      
      # E Step
      
      # Compute loss matrix
      #st=Sys.time()
      loss_mx <- matrix(NA, nrow(Y), nrow(mu)) # (TxM)xK
      # Bottleneck, need to vectorize
      for (r in 1:nrow(YY)) {
        for (k in 1:nrow(mu)) {
          loss_mx[r, k] <- .5*sqrt(sum((YY[r, ] - mu[k, ])^2))
        }
      }
      # en=Sys.time()
      # en-st
      
      #st2=Sys.time()
      loss_mx_new <- 0.5 * sqrt(outer(1:nrow(YY), 1:nrow(mu), 
                                  Vectorize(function(r, k) sum((YY[r, ] - mu[k, ])^2))))
      # en2=Sys.time()
      # en2-st2
      # Continuous model
      
      prob_vecs <- discretize_prob_simplex(K, grid_size)
      pairwise_l1_dist <- as.matrix(dist(prob_vecs, method = "manhattan")) / 2
      jump_penalty_mx <- jump_penalty * (pairwise_l1_dist ^ alpha)
      
      if (mode_loss) {
        # Adding mode loss
        m_loss <- log(rowSums(exp(-jump_penalty_mx)))
        m_loss <- m_loss - m_loss[1]  # Offset a constant
        jump_penalty_mx <- jump_penalty_mx + m_loss
      }
      
      LARGE_FLOAT=1e1000
      # Handle continuous model with probability vectors
      if (!is.null(prob_vecs)) {
        loss_mx[is.nan(loss_mx)] <- LARGE_FLOAT
        loss_mx <- loss_mx %*% t(prob_vecs) # (TxM)xN
      }
      
      N <- ncol(loss_mx)
      
      loss_mx[is.nan(loss_mx)] <- Inf
      
      # DP iteration (bottleneck)
      for(m in 1:M){
        
        values <- matrix(NA, TT*M, N) # (TxM)xN
        assign <- integer(TT*M) #TxM
        
        # Initial step
        indx=which(Y$m==m)
        values[1, ] <- loss_mx[indx[1], ]
        
        # Add spatial penalty
        SSt=SS[SS$t==1,]
        
        # Per ogni m, devo penalizzare le proposte in prob_vecs in base alla distanza di Hellinger
        # pesata per la distanza spaziale
        # Forse conviene ragionare da subito in termini di adiacenze altrimenti e' un macello
        
        dist_vec=apply(prob_vecs,1,function(x)hellinger_distance(SSt[m,-(1:2)],x))
        
        dist_weights <- ifelse(D[m,] == 0, 0, 1 / D[m,])  
        #dist_weights <- dist_weights / sum(dist_weights)
        
        dist_vec=dist_vec*dist_weights
        loss_mx[indx[1], ]=loss_mx[indx[1], ]+spatial_penalty*dist_vec
        
        for (t in 2:TT) {
          
          values[t, ] <- loss_mx[t, ] + 
            apply(values[t - 1, ] + jump_penalty_mx, 2, min)
        }
        
        S=matrix(NA,nrow=TT,ncol=K)
        
        # Find optimal path backwards
        assign[TT] <- which.min(values[TT, ])
        value_opt <- values[TT, assign[TT]]
        
        S[TT,]=prob_vecs[assign[TT],]
        
        # Traceback
        for (t in (TT - 1):1) {
          assign[t] <- which.min(values[t, ] + jump_penalty_mx[, assign[t + 1]])
          S[t,]=prob_vecs[assign[t],]
        }
        
      }
      
      
      # M Step
      
      for(k in 1:K){
        # What if the denominator is 0??
        mu[k,]=apply(Y, 2, function(x) weighted.mean(x, w = S[,k]))
      }
      
      # Re-fill-in missings
      for(i in 1:P){
        Y[,i]=ifelse(M[,i],mu[apply(S,1,which.max),i],Y[,i])
      }
      
      if (!is.null(tol)) {
        epsilon <- loss_old - value_opt
        if (epsilon < tol) {
          break
        }
      } 
      
      else if (all(S == S_old)) {
        break
      }
      
      S_old=S
      
      loss_old <- value_opt
    }
    
    
    
    return(list(S=S,value_opt=value_opt))}, error = function(e) {
      # Return a consistent placeholder on error
      return(list(S = NA, value_opt = Inf))
    })
}