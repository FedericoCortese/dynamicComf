inv_trans_theta=function(data){
  theta_inv=data$theta
  theta_inv[which(theta_inv<0)]=theta_inv[which(theta_inv<0)]+2*pi
  data$theta <- theta_inv
  return(data)
}
trans_theta=function(data){
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta <- theta_trans
  return(data)
}
# propagation_2001GO2_new_v2

df_2001GO2=read.table("data_asteroids/propagation_2001GO2_new_v2.txt",header = T)
unique(df_2001GO2$type)

# propagation_2002AA29_new_v2
df_2002AA29=read.table("data_asteroids/propagation_2002AA29_new_v2.txt",header = T)
df_2002AA29=trans_theta(df_2002AA29)
plot(df_2002AA29$theta,type='p',col=df_2002AA29$type+2)

# propagation_2014OL339_new_v2

df_2014OL339=read.table("data_asteroids/propagation_2014OL339_new_v2.txt",header = T)
df_2014OL339=trans_theta(df_2014OL339)
unique(df_2014OL339$type)
plot(df_2014OL339$theta,type='p',col=df_2014OL339$type+2)

# propagation_2015SO2_new_v2
df_2015SO2=read.table("data_asteroids/propagation_2015SO2_new_v2.txt",header = T)
df_2015SO2=trans_theta(df_2015SO2)
plot(df_2015SO2$theta,type='p',col=df_2015SO2$type+2)

# propagation_2015XX169_new_v2
df_2015XX169=read.table("data_asteroids/propagation_2015XX169_new_v2.txt",header = T)
df_2015XX169=trans_theta(df_2015XX169)
plot(df_2015XX169$theta,type='p',col=df_2015XX169$type+2)

# propagation_2016CA138_new_v2
df_2016CA138=read.table("data_asteroids/propagation_2016CA138_new_v2.txt",header = T)
df_2016CA138=trans_theta(df_2016CA138)
plot(df_2016CA138$theta,type='p',col=df_2016CA138$type+2)

# propagation_2016CO246_new_v2
df_2016CO246=read.table("data_asteroids/propagation_2016CO246_new_v2.txt",header = T)
df_2016CO246=trans_theta(df_2016CO246)
plot(df_2016CO246$theta,type='p',col=df_2016CO246$type+2)

# propagation_2016HO3_new_v2
df_2016HO3=read.table("data_asteroids/propagation_2016HO3_new_v2.txt",header = T)
df_2016HO3=trans_theta(df_2016HO3)
plot(df_2016HO3$theta,type='p',col=df_2016HO3$type+2)

# propagation_2019GM1_new_v2
df_2019GM1=read.table("data_asteroids/propagation_2019GM1_new_v2.txt",header = T)
df_2019GM1=trans_theta(df_2019GM1)
plot(df_2019GM1$theta,type='p',col=df_2019GM1$type+2)

# propagation_2020PN1_new_v2
df_2020PN1=read.table("data_asteroids/propagation_2020PN1_new_v2.txt",header = T)
df_2020PN1=trans_theta(df_2020PN1)
plot(df_2020PN1$theta,type='p',col=df_2020PN1$type+2)

# propagation_2020PP1_new_v2
df_2020PP1=read.table("data_asteroids/propagation_2020PP1_new_v2.txt",header = T)
df_2020PP1=trans_theta(df_2020PP1)
plot(df_2020PP1$theta,type='p',col=df_2020PP1$type+2)


# propagation_164207_new_v2
df_164207=read.table("data_asteroids/propagation_164207_new_v2.txt",header = T)
df_164207=trans_theta(df_164207)
plot(df_164207$theta,type='p',col=df_164207$type+2)


