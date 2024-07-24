source("Utils.R")

TT=4
Ktrue=3
seed=1
M=100
P=25
mu=3
phi=.8
rho=0
Pcat=NULL
pNAs=0
sp_indx=1:M
sp_indx=matrix(sp_indx,ncol=sqrt(M),byrow=T)

Y=list()
S_true=matrix(0,nrow=TT,ncol=M)

C=Cmatrix(sp_indx)
# simDat=sim_spatial_JM(P,C,seed,
#                       rho=rho,Pcat=Pcat, phi=phi,
#                       mu=mu,pNAs=pNAs)

Y=NULL
for(t in 1:TT){
  simDat=sim_spatial_JM(P,C,t,
                        rho=rho,Pcat=Pcat, phi=phi,
                        mu=mu,pNAs=pNAs)
  temp=simDat$SimData
  temp$m=1:M
  temp$t=rep(t,M)
  Y=rbind(Y,temp)
  S_true[t,]=simDat$s
}

### Function ###
n_states=3
cat.indx=which(sapply(Y, is.factor))
cont.indx=which(sapply(Y, is.numeric))

Ycont=Y[,cont.indx]
Ycat=Y[,cat.indx]

n_cat=length(cat.indx)
n_cont=P-n_cat

# Missing data imputation TBD

# State initialization through kmeans++

S=matrix(0,nrow=TT,ncol=M)

for(m in 1:M){
  S[,m]=initialize_states(Y[which(Y$m==m),],n_states)
}

for (init in 1:n_init) {
  mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
  mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
  loss_old <- 1e10
  for (it in 1:max_iter) {
    
    for (i in unique(as.vector(S))) {
      mu[i,] <- colMeans(Ycont[as.vector(t(S))==i,])
      mo[i,]=apply(Ycat[as.vector(t(S))==i,],2,Mode)
    }
    
    mu=data.frame(mu)
    mo=data.frame(mo,stringsAsFactors=TRUE)
    for(i in 1:n_cat){
      #mo[,i]=factor(mo[,i],levels=1:n_levs[i])
      mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
      
    }
    mumo=data.frame(matrix(0,nrow=n_states,ncol=P))
    mumo[,cat.indx]=mo
    mumo[,cont.indx]=mu
    colnames(mumo)=colnames(Y)
    
    # Fit state sequence
    S_old <- S
    
    # Re-fill-in missings TBD 
    
    loss_by_state=gower.dist(Y,mumo)
    
    V <- loss_by_state
    for (t in (n_obs-1):1) {
      
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
      
    }
    
    s[1] <- which.min(V[1,])
    for (t in 2:n_obs) {
      s[t] <- which.min(V[t,] + Gamma[s[t-1],])
    }
    
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
    } else if (all(s == s_old)) {
      break
    }
    loss_old <- loss
  }
  if (is.null(best_s) || (loss_old < best_loss)) {
    best_loss <- loss_old
    best_s <- s
  }
  #s <- init_states(Y, n_states)+1
  s=initialize_states(Y,n_states)
}
