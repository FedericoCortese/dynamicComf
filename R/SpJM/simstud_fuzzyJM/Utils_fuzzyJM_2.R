library(cluster)
library(StatMatch)
library(poliscidata) # for weighted mode 
library(pdfCluster) # for ARI
library(parallel)
library(foreach)
library(doParallel)

weighted_median <- function(x, weights) {
  # Ensure x and weights are of the same length
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length.")
  }
  
  # Sort x and weights by x
  sorted_indices <- order(x)
  x <- x[sorted_indices]
  weights <- weights[sorted_indices]
  
  # Compute cumulative weights
  cumulative_weights <- cumsum(weights)
  total_weight <- sum(weights)
  
  # Find the smallest x such that the cumulative weight is >= 50% of the total weight
  weighted_median <- x[which(cumulative_weights >= total_weight / 2)[1]]
  
  return(weighted_median)
}


weighted_mode <- function(x, weights) {
  # Ensure x and weights are of the same length
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length.")
  }
  
  # Aggregate weights for each unique value of x
  unique_x <- unique(x)
  aggregated_weights <- sapply(unique_x, function(val) sum(weights[x == val]))
  
  # Find the value of x with the maximum aggregated weight
  weighted_mode <- unique_x[which.max(aggregated_weights)]
  
  return(weighted_mode)
}


Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

get_cat=function(y,mc,mu_val,phi){
  # Function to simulate categorical data
  
  # Arguments:
  # y: continuous variable 
  # mc: Markov chain states
  # mu: numeric mean value
  # phi: conditional probability for the categorical outcome k in state k
  
  library(dplyr)
  
  mu_val=c(-mu_val,0,mu_val)
  #K=length(unique(mc))
  #mu=seq(-mu,mu,length.out=K)
  #phi1=(1-phi)/(K-1)
  phi1=(1-phi)/2
  
  TT=length(y)
  for(i in 1:TT){
    k=mc[i]
    switch(k,
           "1"={
             threshold=c(qnorm(phi1,mu_val[1]),qnorm(phi+phi1,mu_val[1]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=1
             }
             else if(y[i]<threshold[1]){
               y[i]=2
             }
             else{
               y[i]=3
             }
           },
           "2"={
             threshold=c(qnorm(phi1,mu_val[2]),qnorm(phi+phi1,mu_val[2]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=2
             }
             else if(y[i]<threshold[1]){
               y[i]=3
             }
             else{
               y[i]=1
             }
           },
           "3"={
             threshold=c(qnorm(phi1,mu_val[3]),qnorm(phi+phi1,mu_val[3]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=3
             }
             else if(y[i]<threshold[1]){
               y[i]=1
             }
             else{
               y[i]=2
             }
           }
    )
  }
  return(y)
  
}

order_states_condMean=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional mean for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMeans=sort(tapply(y,s,mean,na.rm=T))
  
  states_temp=match(s,names(condMeans))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

order_states_condMed=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional median for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMed=sort(tapply(y,s,median,na.rm=T))
  
  states_temp=match(s,names(condMed))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

punct=function(x,pNAs,typeNA){
  
  # x is a vector (column of the dataset)
  # pNAs is the percentage of missing values
  # typeNA is the type of missing values (0 for random, 1 for continuous, all other values will turn into no missing imputation)
  
  TT=length(x)
  pTT=round(TT*pNAs)
  if(typeNA==0){
    NAindx=sample(1:TT,pTT,replace = F)
    x[NAindx]=NA
  }
  else if(typeNA==1){
    NAindx=sample(1:(TT-pTT),1,replace = F)
    NAindx=seq(NAindx,NAindx+pTT)
    x[NAindx]=NA
  }
  
  return(x)
  
}

initialize_states <- function(Y, K) {
  n <- nrow(Y)
  
  ### Repeat the following few times?
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(daisy(Y, metric = "gower"))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  ###
  
  # init_stats=rep(0,n)
  # For cycle
  # for(i in 1:n){
  #   init_stats[i]=which.min(gower.dist(Y[i,],centroids))
  # }
  
  # Using sapply and vapply
  # init_stats2 <- sapply(1:n, function(i) which.min(gower.dist(Y[i,], centroids)))
  # init_stats3 <- vapply(1:n, function(i) which.min(gower.dist(Y[i,], centroids)), integer(1))
  
  # Faster solution 
  dist_matrix <- gower.dist(Y, centroids)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}
######
objective_function <- function(s, g_values, lambda, s_t_prec,s_t_succ,m) {
  #sum(s^m * g_values) + lambda * sum((s_t_prec - s)^2+(s_t_succ - s)^2) 
  sum(s^m * g_values) + lambda *( sum(abs(s_t_prec - s))^2+sum(abs(s_t_succ - s))^2) 
}

objective_function_1T <- function(s, g_values, lambda, s_t_1,m) {
  #sum(s^m * g_values) + lambda * sum((s_t_1 - s)^2) 
  sum(s^m * g_values) + lambda * sum(abs(s_t_1 - s))^2 
}

fuzzy_jump_coord <- function(Y, 
                             K, 
                             lambda=1e-5, 
                             m=1,
                             max_iter=5, 
                             n_init=10, tol=1e-16, 
                             verbose=FALSE
                             
) {
  # Fit jump model for mixed type data 
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # K: number of states
  # lambda: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  # time_vec is a vector of time points, needed if times are not equally sampled
  
  # Value:
  # best_s: estimated state sequence
  # Y: imputed data
  # Y.orig: original data
  # condMM: state-conditional medians and modes
  
  K=as.integer(K)
  
  TT <- nrow(Y)
  P <- ncol(Y)
  best_loss <- NULL
  best_S <- NULL
  
  # Which vars are categorical and which are numeric
  cat_flag=any(sapply(Y, is.factor))
  
  if(cat_flag){
    cat.indx=which(sapply(Y, is.factor))
    cont.indx=which(sapply(Y, is.numeric))
    Ycont=Y[,cont.indx]
    Ycont=apply(Ycont,2,scale)
    Y[,cont.indx]=Ycont
    Ycat=Y[,cat.indx]
    
    if(length(cat.indx)==1){
      n_levs=length(levels(Ycat))
    }
    else{
      n_levs=apply(Ycat, 2, function(x)length(unique(x[!is.na(x)])))
    }
    
    n_cat=length(cat.indx)
    n_cont=P-n_cat
    
  }
  else{
    cont.indx=1:P
    Ycont=Y
    Ycont=apply(Ycont,2,scale)
    Y[,cont.indx]=Ycont
    n_cont=dim(Y)[2]
    n_cat=0
  }
  
  
  
  
  for (init in 1:n_init) {
    
    # k-prot++
    s=initialize_states(Y,K)
    S <- matrix(0, nrow = TT, ncol = K)
    row_indices <- rep(1:TT)  
    S[cbind(row_indices, s)] <- 1 
    
    
    # initialize mumo
    mu <- matrix(NA, nrow=K, ncol=length(cont.indx))
    if(cat_flag){
      mo <- matrix(NA, nrow=K, ncol=length(cat.indx))
    }
    
    for(k in 1:K){
      mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=S[,k]^m)})
      if(cat_flag){
        if(n_cat==1){
          mo[k,]=poliscidata::wtd.mode(Ycat,weights=S[,k]^m)
          
        }
        else{
          mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=S[,k]^m)})
        }
      }
    }
    
    mumo=data.frame(matrix(0,nrow=K,ncol=P))
    if(cat_flag){
      mumo[,cat.indx]=mo
    }
    
    mumo[,cont.indx]=mu
    colnames(mumo)=colnames(Y)
    
    if(cat_flag){
      for(p in 1:n_cat){
        if(n_cat==1){
          mumo[,cat.indx[p]]=factor(mumo[,cat.indx[p]],levels=levels(Ycat[p]))
        }
        else{
          mumo[,cat.indx[p]]=factor(mumo[,cat.indx[p]],levels=levels(Ycat[,p]))
        }
      }
    }
    
    S_old=S
    V=gower.dist(Y,mumo)
    loss_old=sum(V*S^m)+lambda*sum(abs(S[1:(TT-1),]-S[2:TT,]))^2
    
    for (it in 1:max_iter) {
      
      # S(1)
      
      result <- Rsolnp::solnp(#pars = S[1,],
        pars = rep(1/K,K),
        fun = function(s) objective_function_1T(s, 
                                                g_values=V[1,], 
                                                lambda=lambda, 
                                                s_t_1=S[2,], 
                                                m=m),
        eqfun = function(s) sum(s),
        eqB = 1,
        LB = rep(0, K),
        control = list(trace = 0))
      
      S[1,] <- result$pars
      
      # S(2) to S(T-1)
      
      for(t in 2:(TT-1)){
        result=Rsolnp::solnp(
          #pars = S[t,],
          pars = rep(1/K,K),
          fun = function(s) objective_function(s, 
                                               g_values=V[t,], 
                                               lambda=lambda, 
                                               s_t_prec=S[t-1,],
                                               s_t_succ=S[t+1,],
                                               m=m),
          eqfun = function(s) sum(s),
          eqB = 1,
          LB = rep(0, K),
          control = list(trace = 0))
        S[t,]=result$pars
        
      }
      
      # S(T)
      result <- Rsolnp::solnp(
        #pars = S[TT,],
        pars = rep(1/K,K),
        fun = function(s) objective_function_1T(s, 
                                                g_values=V[TT,], 
                                                lambda=lambda, 
                                                s_t_1=S[TT-1,], 
                                                m=m),
        eqfun = function(s) sum(s),
        eqB = 1,
        LB = rep(0, K),
        control = list(trace = 0))
      
      S[TT,] <- result$pars
      
      
      for(k in 1:K){
        #mu[k,]=apply(Ycont, 2, function(x) weighted_median(x, weights = S[,k]))
        #mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=S[,k])})
        mu[k,]=apply(Ycont,2,function(x){poliscidata::wtd.median(x,weights=S[,k]^m)})
        if(cat_flag){
          if(n_cat==1){
            #mo[k,]=poliscidata::wtd.mode(Ycat,weights=S[,k])
            mo[k,]=poliscidata::wtd.mode(Ycat,weights=S[,k]^m)
            
          }
          else{
            #mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=S[,k])})
            mo[k,]=apply(Ycat,2,function(x){poliscidata::wtd.mode(x,weights=S[,k]^m)})
            
          }
          #mo[k,]=apply(Ycat,2,function(x)weighted_mode(x,weights=S[,k]))
        }
      }
      
      if(cat_flag){
        mumo[,cat.indx]=mo
      }
      
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(Y)
      if(cat_flag){
        for(p in 1:n_cat){
          if(n_cat==1){
            mumo[,cat.indx[p]]=factor(mumo[,cat.indx[p]],levels=levels(Ycat[p]))
          }
          else{
            mumo[,cat.indx[p]]=factor(mumo[,cat.indx[p]],levels=levels(Ycat[,p]))
          }
        }
      }
      
      V=gower.dist(Y,mumo)
      loss=sum(V*S^m)+lambda*sum(abs(S[1:(TT-1),]-S[2:TT,]))^2
      
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
      S_old=S
    }
    if (is.null(best_S) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_S <- S
    }
  }
  
  old_MAP=apply(best_S,1,which.max)
  MAP=order_states_condMed(Y[,cont.indx[1]],old_MAP)
  
  tab <- table(MAP, old_MAP)
  new_order <- apply(tab, 1, which.max)
  
  # Reorder the columns of S accordingly
  best_S <- best_S[, new_order]
  
  # res_Y=data.frame(Y,MAP=MAP)
  # col_sort=as.integer(names(sort(tapply(res_Y[,cont.indx[1]],
  #                                       res_Y$MAP,mean))))
  #mumo=mumo[col_sort,]
  # best_S=best_S[,col_sort]
  
  return(list(best_S=best_S,
              MAP=MAP,
              Y=Y
              # ,
              # condMM=mumo
  ))
}

fuzzy_jump_coord_par <- function(Y, 
                                 K, 
                                 lambda=1e-5, 
                                 m=1.01,
                                 max_iter=5, 
                                 n_init=10, tol=1e-16, 
                                 verbose=FALSE,
                                 n_cores=NULL
                                 
) {
  # Fit jump model for mixed type data 
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # K: number of states
  # lambda: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  # time_vec is a vector of time points, needed if times are not equally sampled
  
  # Value:
  # best_s: estimated state sequence
  # Y: imputed data
  # Y.orig: original data
  # condMM: state-conditional medians and modes
  
  K=as.integer(K)
  
  TT <- nrow(Y)
  P <- ncol(Y)
  best_loss <- NULL
  best_S <- NULL
  
  # Which vars are categorical and which are numeric
  cat_flag=any(sapply(Y, is.factor))
  
  if(cat_flag){
    cat.indx=which(sapply(Y, is.factor))
    cont.indx=which(sapply(Y, is.numeric))
    Ycont=Y[,cont.indx]
    Ycont=apply(Ycont,2,scale)
    Y[,cont.indx]=Ycont
    Ycat=Y[,cat.indx]
    
    if(length(cat.indx)==1){
      n_levs=length(levels(Ycat))
    }
    else{
      n_levs=apply(Ycat, 2, function(x)length(unique(x[!is.na(x)])))
    }
    
    n_cat=length(cat.indx)
    n_cont=P-n_cat
    
  }
  else{
    cont.indx=1:P
    Ycont=Y
    Ycont=apply(Ycont,2,scale)
    Y[,cont.indx]=Ycont
    n_cont=dim(Y)[2]
    n_cat=0
  }
  
  
  # Imposta i core da usare
  if(is.null(n_cores)){
    n_cores <- parallel::detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Esegui le inizializzazioni in parallelo
  results <- foreach(init = 1:n_init, .packages = c("poliscidata", "Rsolnp",
                                                    "cluster","StatMatch"),
                     .export = c("initialize_states",
                                 "objective_function_1T",
                                 "objective_function")) %dopar% {
                                   
                                   s = initialize_states(Y, K)
                                   S <- matrix(0, nrow = TT, ncol = K)
                                   row_indices <- rep(1:TT)
                                   S[cbind(row_indices, s)] <- 1
                                   
                                   mu <- matrix(NA, nrow=K, ncol=length(cont.indx))
                                   if(cat_flag){
                                     mo <- matrix(NA, nrow=K, ncol=length(cat.indx))
                                   }
                                   
                                   for(k in 1:K){
                                     mu[k,] = apply(Ycont, 2, function(x) poliscidata::wtd.median(x, weights = S[,k]^m))
                                     if(cat_flag){
                                       if(n_cat == 1){
                                         mo[k,] = poliscidata::wtd.mode(Ycat, weights = S[,k]^m)
                                       } else {
                                         mo[k,] = apply(Ycat, 2, function(x) poliscidata::wtd.mode(x, weights = S[,k]^m))
                                       }
                                     }
                                   }
                                   
                                   mumo <- data.frame(matrix(0, nrow=K, ncol=P))
                                   if(cat_flag){
                                     mumo[,cat.indx] <- mo
                                   }
                                   mumo[,cont.indx] <- mu
                                   colnames(mumo) <- colnames(Y)
                                   if(cat_flag){
                                     for(p in 1:n_cat){
                                       if(n_cat==1){
                                         mumo[,cat.indx[p]] <- factor(mumo[,cat.indx[p]], levels=levels(Ycat[p]))
                                       }
                                       else{
                                         mumo[,cat.indx[p]] <- factor(mumo[,cat.indx[p]], levels=levels(Ycat[,p]))
                                       }
                                     }
                                   }
                                   
                                   S_old <- S
                                   V <- gower.dist(Y, mumo)
                                   loss_old <- sum(V * S^m) + lambda * sum(abs(S[1:(TT-1), ] - S[2:TT, ]))^2
                                   
                                   for (it in 1:max_iter) {
                                     result <- Rsolnp::solnp(
                                       pars = rep(1/K, K),
                                       fun = function(s) objective_function_1T(s, g_values=V[1,], lambda=lambda, s_t_1=S[2,], m=m),
                                       eqfun = function(s) sum(s),
                                       eqB = 1,
                                       LB = rep(0, K),
                                       control = list(trace = 0)
                                     )
                                     S[1,] <- result$pars
                                     
                                     for(t in 2:(TT-1)){
                                       result <- Rsolnp::solnp(
                                         pars = rep(1/K, K),
                                         fun = function(s) objective_function(s, g_values=V[t,], lambda=lambda, s_t_prec=S[t-1,], s_t_succ=S[t+1,], m=m),
                                         eqfun = function(s) sum(s),
                                         eqB = 1,
                                         LB = rep(0, K),
                                         control = list(trace = 0)
                                       )
                                       S[t,] <- result$pars
                                     }
                                     
                                     result <- Rsolnp::solnp(
                                       pars = rep(1/K, K),
                                       fun = function(s) objective_function_1T(s, g_values=V[TT,], lambda=lambda, s_t_1=S[TT-1,], m=m),
                                       eqfun = function(s) sum(s),
                                       eqB = 1,
                                       LB = rep(0, K),
                                       control = list(trace = 0)
                                     )
                                     S[TT,] <- result$pars
                                     
                                     
                                     for(k in 1:K){
                                       mu[k,] = apply(Ycont, 2, function(x) poliscidata::wtd.median(x, weights = S[,k]^m))
                                       if(cat_flag){
                                         if(n_cat == 1){
                                           mo[k,] = poliscidata::wtd.mode(Ycat, weights = S[,k]^m)
                                         } else {
                                           mo[k,] = apply(Ycat, 2, function(x) poliscidata::wtd.mode(x, weights = S[,k]^m))
                                         }
                                       }
                                     }
                                     
                                     if(cat_flag){
                                       mumo[,cat.indx] <- mo
                                     }
                                     mumo[,cont.indx] <- mu
                                     colnames(mumo) <- colnames(Y)
                                     if(cat_flag){
                                       for(p in 1:n_cat){
                                         if(n_cat==1){
                                           mumo[,cat.indx[p]] <- factor(mumo[,cat.indx[p]], levels=levels(Ycat[p]))
                                         }
                                         else{
                                           mumo[,cat.indx[p]] <- factor(mumo[,cat.indx[p]], levels=levels(Ycat[,p]))
                                         }
                                       }
                                     }
                                     
                                     V <- gower.dist(Y, mumo)
                                     
                                     loss <- sum(V * S^m) + lambda * sum(abs(S[1:(TT-1), ] - S[2:TT, ]))^2
                                     
                                     if(!is.null(tol)){
                                       if((loss_old - loss) < tol) break
                                     } else if(all(S == S_old)) {
                                       break
                                     }
                                     loss_old <- loss
                                     S_old <- S
                                   }
                                   
                                   list(S = S, loss = loss_old)
                                 }
  
  # Trova il migliore
  losses <- sapply(results, function(res) res$loss)
  best_index <- which.min(losses)
  best_S <- results[[best_index]]$S
  best_loss <- losses[best_index]
  
  # Ferma il cluster
  stopCluster(cl)
  
  
  
  old_MAP=apply(best_S,1,which.max)
  MAP=order_states_condMed(Y[,cont.indx[1]],old_MAP)
  
  tab <- table(MAP, old_MAP)
  new_order <- apply(tab, 1, which.max)
  
  # Reorder the columns of S accordingly
  best_S <- best_S[, new_order]
  
  return(list(best_S=best_S,
              MAP=MAP,
              Y=Y
  ))
}

get_cat_t <- function(y, mc, mu_val, phi, df = 5) {
  # y: vettore continuo (una variabile)
  # mc: stati latenti (valori da 1 a K)
  # mu_val: vettore di K valori medi (es. c(-mu, 0, mu))
  # phi: probabilità condizionata
  # df: gradi di libertà della t-Student
  
  TT <- length(y)
  K <- length(mu_val)
  phi1 <- (1 - phi) / (K - 1)  # non serve direttamente, solo se vuoi usare probabilità
  
  for (i in 1:TT) {
    k <- mc[i]
    mu_k <- mu_val[k]
    
    # Calcola l'intervallo centrato in mu_k che contiene prob = phi sotto la t-Student
    half_width <- qt((1 + phi) / 2, df = df)
    lower <- mu_k - half_width
    upper <- mu_k + half_width
    
    # Verifica se y[i] cade nell'intervallo
    if (y[i] >= lower && y[i] <= upper) {
      y[i] <- k  # corretto: valore coerente con lo stato
    } else {
      # Scegli casualmente un altro k tra gli altri K-1
      y[i] <- sample(setdiff(1:K, k), 1)
    }
  }
  return(y)
}

sim_data_stud_t=function(seed=123,
                         TT,
                         P,
                         Pcat,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=4,
                         phi=.8,
                         pers=.95){
  
  
  MU=seq(mu, -mu, length.out=Ktrue)
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    # u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    u = mvtnorm::rmvt(TT, sigma = (nu-2)*Sigma/nu, df = nu, delta = rep(MU[k],P))
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  SimData=data.frame(SimData)
  
  if(!is.null(Pcat)){
    for (j in 1:Pcat) {
      SimData[, j] <- get_cat_t(SimData[, j], x, MU, phi=phi, df = nu)
      SimData[, j]=factor(SimData[, j],levels=1:Ktrue)
    }  
  }
  
  
  
  return(list(
    SimData=SimData,
    mchain=x,
    TT=TT,
    P=P,
    K=Ktrue,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}


simstud_fuzzyJM_coord=function(seed,
                               lambda,
                               TT,
                               P,
                               K,
                               mu=1.5,
                               phi=.8,
                               rho=0,
                               Pcat=NULL,
                               pers=.9,
                               m=1.5,
                               nu=4){
  # Simulate
  simDat=sim_data_stud_t(seed=123,
                         TT=TT,
                         P=P,
                         Pcat=Pcat,
                         Ktrue=K,
                         mu=mu,
                         rho=rho,
                         nu=nu,
                         phi=phi,
                         pers=pers)
  
  Y=simDat$SimData
  
  success <- FALSE
  trials=1
  while (!success&trials<10) {
    est <- try(fuzzy_jump_coord(Y=Y, 
                                K=K, 
                                lambda=lambda, 
                                m=m,
                                max_iter=10, 
                                n_init=10, tol=1e-8, 
                                verbose=FALSE
                                
    ), silent = TRUE)
    trials=trials+1
    
    if (!inherits(est, "try-error")) {
      success <- TRUE  # Exit the loop if no error
    } else {
      message("Retrying fuzzy_jump() due to an error...")
    }
  }
  
  
  
  old_MAP=factor(est$MAP,levels=1:K)
  MAP=factor(relabel_clusters(est$MAP,simDat$mchain),levels=1:K)
  #est$MAP=factor(relabel_clusters(est$MAP,simDat$mchain),levels=1:K)
  simDat$mchain=factor(simDat$mchain,levels=1:K)
  
  BAC=caret::confusionMatrix(MAP,simDat$mchain)$overall[1]
  ARI=adj.rand.index(MAP,simDat$mchain)
  
  tab <- table(MAP, old_MAP)
  new_order <- apply(tab, 1, which.max)
  best_S=est$best_S[, new_order]
  
  fRI=fclust::ARI.F(simDat$mchain,best_S)
  
  # Return
  return(list(
    S=best_S,
    MAP=MAP ,
    ARI=ARI,
    BAC=BAC,
    fRI=fRI,
    seed=seed,
    lambda=lambda,
    TT=TT,
    P=P,
    K=K,
    Pcat=Pcat,
    true_seq=simDat$mchain
  ))
  
}

######