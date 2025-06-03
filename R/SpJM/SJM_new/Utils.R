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