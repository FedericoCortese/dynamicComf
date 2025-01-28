# 1) Provare a sostituire qualche ciclo for con apply
# 2) Sostituire media euclidea con distanza di Gower

onerun_cont_STJM=function(Y,n_states,D,
                       jump_penalty,
                       spatial_penalty,
                       alpha,grid_size,mode_loss=T,
                       max_iter,tol,initial_states){
  tryCatch({
    
    TT=length(unique(Y$t))
    M=length(unique(Y$m))
    
    
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
    
    mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    
    for (i in unique(as.vector(S))) {
      # substitute with medians
      mu[i,] <- apply(Ycont[as.vector(t(S))==i,], 2, median, na.rm = TRUE)
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
    
    for (it in 1:max_iter) {
      
      # E Step
      
      loss_mx=gower.dist(YY,mumo)
      
      # Slighlty faster
      # loss_mx<- 0.5 * sqrt(outer(1:nrow(YY), 1:nrow(mu), 
      #                             Vectorize(function(r, k) sum((YY[r, ] - mu[k, ])^2))))
      
      
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
      
      values <- matrix(NA, TT*M, N) # (TxM)xN
      value_opt=rep(0,M)
      assign <- integer(TT*M) #TxM
      SS_new=SS
      SS_new[,-(1:2)]=0
      
      # DP iteration (bottleneck)
      for(m in 1:M){
        #print(m)
        # Verificare se i pesi spaziali vanno bene
        dist_weights <- ifelse(D[m,] == 0, 0, 1 / D[m,]) 
        
        spat_weigh_Ndim=rep(0,N)
        
        for(n in 1:N){
          spat_weigh_Ndim[n]=
            sum(apply(SS[SS$t==1,-(1:2)],1,
                      function(x)hellinger_distance(prob_vecs[n,],x))*dist_weights)
        }
        
        #dist_weights <- dist_weights / sum(dist_weights)
        indx=which(Y$m==m)
      
        values[indx[1], ] <- loss_mx[indx[1], ]+spatial_penalty*spat_weigh_Ndim

        # Bootleneck!!
        for (t in 2:TT) {
          #print(t)
          spat_weigh_Ndim=rep(0,N)
          
          for(n in 1:N){
            spat_weigh_Ndim[n]=
              sum(apply(SS[SS$t==t,-(1:2)],1,function(x)hellinger_distance(prob_vecs[n,],x))*dist_weights)
          }
          
          values[indx[t], ] <- loss_mx[indx[t], ] + 
            apply(values[indx[t-1], ] + jump_penalty_mx, 2, min)+
            spatial_penalty*spat_weigh_Ndim
        }
        
        # Find optimal path backwards
        assign[indx[TT]] <- which.min(values[indx[TT], ])
        value_opt[m] <- values[indx[TT], assign[indx[TT]]]
        
        SS_new[indx[TT],-(1:2)]=prob_vecs[assign[indx[TT]],]
        
        # Traceback
        for (t in (TT - 1):1) {
          assign[indx[t]] <- which.min(values[indx[t], ] + 
                                          jump_penalty_mx[, assign[indx[t+1]]])
          SS_new[indx[t],-(1:2)]=prob_vecs[assign[indx[t]],]
        }
        
      }
      value_opt=mean(value_opt)
      
      # M Step
      
      for(k in 1:K){
        mu[k,]=apply(Ycont, 2, function(x) weighted_median(x, weights = SS_new[,k+2]))
        mo[k,]=apply(Ycat,2,function(x)weighted_mode(x,weights=SS_new[,k+2]))
      }
      
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(YY)
      
      # Re-fill-in missings 
      Ycont <- Mcont * mu[as.vector(t(S)), ] + (!Mcont) * Ycont
      
      mo_2=apply(mo,2,as.numeric)
      Ycat <- Mcat * mo_2[as.vector(t(S)), ] + (!Mcat) * apply(Ycat,2,as.numeric)
      Ycat=data.frame(Ycat)
      for(cc in 1:n_cat){
        Ycat[,cc]=factor(Ycat[,cc],levels=levels_cat[[cc]])
      }
      
      YY[,-cat.indx]=Ycont
      YY[,cat.indx]=Ycat
      
      Y[,-(1:2)]=YY
      
      
      if (!is.null(tol)) {
        epsilon <- loss_old - value_opt
        if (epsilon < tol) {
          break
        }
      } 
      
      else if (all(SS == SS_new)) {
        break
      }
      
      SS=SS_new
      
      loss_old <- value_opt
    }
    
    
    
    return(list(S=SS,value_opt=value_opt,mumo=mumo))}, 
    error = function(e) {
      # Return a consistent placeholder on error
      return(list(S = NA, value_opt = Inf))
    })
}


prova=onerun_cont_STJM(Y,K=3,D,
                       jump_penalty = .1,
                       spatial_penalty = .1,
                       alpha=2,
                       grid_size = .05,
                       mode_loss=T,
                       max_iter=5,
                       tol=NULL,
                       initial_states = NULL)

apply(prova$S[prova$S$m==1,-(1:2)],1,which.max)
result$S[,1]
apply(prova$S[prova$S$m==2,-(1:2)],1,which.max)
result$S[,2]

