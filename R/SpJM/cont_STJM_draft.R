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
  
   
  
}