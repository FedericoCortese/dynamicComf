
source("Utils_asteroids.R")

# propagation_2001GO2_new_v2

df_2001GO2=read.table("data_asteroids/propagation_2001GO2_new_v2.txt",header = T)
unique(df_2001GO2$type)

# propagation_2002AA29_new_v2
df_2002AA29=read.table("data_asteroids/propagation_2002AA29_new_v2.txt",header = T)
df_2002AA29=trans_theta(df_2002AA29)
unique(df_2002AA29$type)
plot(df_2002AA29$theta,type='p',col=df_2002AA29$type+1)


# 2014OL339 ---------------------------------------------------------------

# propagation_2014OL339_new_v2

df_2014OL339=read.table("data_asteroids/propagation_2014OL339_new_v2.txt",header = T)

# Recode states
old_vals <- sort(unique(df_2014OL339$type))
df_2014OL339$type <- match(df_2014OL339$type, old_vals)


# 2015XX169 ---------------------------------------------------------------

df_2015XX169=read.table("data_asteroids/propagation_2015XX169_new_v2.txt",header = T)
df_2015XX169=trans_theta(df_2015XX169)

# Select relevant variables
df_2015XX169=df_2015XX169[,c("t","a","theta","type")]
#df_2015XX169$type[which(df_2015XX169$type==100)]=10

# Recode states
old_vals <- sort(unique(df_2015XX169$type))
df_2015XX169$type <- match(df_2015XX169$type, old_vals)

plot(df_2015XX169$theta,col=df_2015XX169$type)

# Compute features

df_2015XX169_maxmins_theta=max_min_feat(df_2015XX169,
                                        "theta",tt_thres_maxmin=2.5,
                                        l=5)

df_2015XX169_stat_theta=compute_feat(df_2015XX169,"theta",l=50,ma_flag   = F)

df_2015XX169_stat_a=compute_feat(df_2015XX169,"a",l=50,sd_flag   = F)
df_2015XX169_stat_a=df_2015XX169_stat_a[,c("t","ma_a")]


# Merge by t
features_2015XX169=merge(df_2015XX169,
                         df_2015XX169_maxmins_theta,by="t")
features_2015XX169=merge(features_2015XX169,
                         df_2015XX169_stat_theta,by="t")
features_2015XX169=merge(features_2015XX169,
                         df_2015XX169_stat_a,by="t")

features_2015XX169=features_2015XX169[complete.cases(features_2015XX169), ]

# Extract ground truth and time
gt_2015XX169=df_2015XX169$type
time_2015XX169=df_2015XX169$t
a_2015XX169=df_2015XX169$a
theta_2015XX169=df_2015XX169$theta

summary(features_2015XX169)

sel_features_2015XX169=features_2015XX169[,c("value_max_theta", "value_min_theta",
                                             "dtheta","sd_theta","sd_dtheta",
                                             #"I_TP",
                                             "I_HS",
                                             "I_QS","I_CP",
                                             "ma_a")]

source("Utils_sparse_robust_2.R")

# Select lambda, K, and c
cv_2015XX169=cv_robust_sparse_jump(
  Y=sel_features_2015XX169,
  true_states=gt_2015XX169,
  K_grid=2:4,
  zeta0=.4,
  lambda_grid=seq(0,1,.1),
  n_folds = 5,
  parallel=T,
  n_cores=detectCores()-1,
  cv_method="blocked-cv",
  knn=10,
  c_grid=c(7.5,10),
  M=NULL
)

# Select zeta0 based on the best lambda, K and c
gap_2015XX169=gap_robust_sparse_jump(
  Y=sel_features_2015XX169,
  K_grid=NULL,
  zeta0_grid=seq(0.05,.5,0.05),
  lambda=0,
  B=10,
  parallel=F,
  n_cores=NULL,
  knn=10,
  c=10,
  M=NULL
)

# Final fit

fit_2015XX169=robust_sparse_jump(
  Y=sel_features_2015XX169,
  K=3,
  zeta0=.2,
  lambda=.7,
  c=10,
  knn=10,
  M=NULL,
  n_init=5,
  verbose=T,
  tol=1e-8
)

est_s_2015XX169=fit_2015XX169$s
est_s_2015XX169[fit_2015XX169$v==0]=0
plot(features_2015XX169$theta,col=est_s_2015XX169+1)
plot(features_2015XX169$a,col=est_s_2015XX169+1)

est_W_2015XX169= data.frame(round(fit_2015XX169$W,2))
colnames(est_W_2015XX169)=names(sel_features_2015XX169)

est_W_2015XX169



# 2016CA138 ---------------------------------------------------------------

df_2016CA138=read.table("data_asteroids/propagation_2016CA138_new_v2.txt",header = T)
df_2016CA138=trans_theta(df_2016CA138)
unique(df_2016CA138$type)

# Select relevant variables
df_2016CA138=df_2016CA138[,c("t","a","theta","type")]
df_2016CA138$type[which(df_2016CA138$type==100)]=10

# Compute features

df_2016CA138_maxmins_theta=max_min_feat(df_2016CA138,"theta",tt_thres_maxmin=2.5,
                                        l=5)

df_2016CA138_stat_theta=compute_feat(df_2016CA138,"theta",l=100,ma_flag   = F)

df_2016CA138_stat_a=compute_feat(df_2016CA138,"a",l=100,sd_flag   = F)
df_2016CA138_stat_a=df_2016CA138_stat_a[,c("t","ma_a")]


# Merge by t
features_2016CA138=merge(df_2016CA138,
                         df_2016CA138_maxmins_theta,by="t")
features_2016CA138=merge(features_2016CA138,
                         df_2016CA138_stat_theta,by="t")
features_2016CA138=merge(features_2016CA138,
                         df_2016CA138_stat_a,by="t")

features_2016CA138=features_2016CA138[complete.cases(features_2016CA138), ]

# Extract ground truth and time
gt_2016CA138=df_2016CA138$type
time_2016CA138=df_2016CA138$t
a_2016CA138=df_2016CA138$a
theta_2016CA138=df_2016CA138$theta

summary(features_2016CA138)

sel_features_2016CA138=features_2016CA138[,c("value_max_theta", "value_min_theta",
                                             "dtheta","sd_theta","sd_dtheta",
                                             "I_TP","I_HS",
                                             "I_QS","I_CP",
                                             "ma_a")]

source("Utils_sparse_robust_2.R")

# Select lambda, K, and c
cv_2016CA138=cv_robust_sparse_jump(
  Y=sel_features_2016CA138,
  true_states=gt_2016CA138,
  K_grid=2:4,
  zeta0=.4,
  lambda_grid=seq(0,1,.1),
  n_folds = 5,
  parallel=F,
  n_cores=NULL,
  cv_method="blocked-cv",
  knn=10,
  c_grid=c(7.5,10),
  M=NULL
)

# Select zeta0 based on the best lambda, K and c
gap_2016CA138=gap_robust_sparse_jump(
  Y=sel_features_2016CA138,
  K_grid=NULL,
  zeta0_grid=seq(0.05,.5,0.05),
  lambda=0,
  B=10,
  parallel=F,
  n_cores=NULL,
  knn=10,
  c=10,
  M=NULL
)

# Final fit

fit_2016CA138=robust_sparse_jump(
  Y=sel_features_2016CA138,
  K=3,
  zeta0=.2,
  lambda=.7,
  c=10,
  knn=10,
  M=NULL,
  n_init=5,
  verbose=T,
  tol=1e-8
)

est_s_2016CA138=fit_2016CA138$s
est_s_2016CA138[fit_2016CA138$v==0]=0
plot(features_2016CA138$theta,col=est_s_2016CA138+1)
plot(features_2016CA138$a,col=est_s_2016CA138+1)

est_W_2016CA138= data.frame(round(fit_2016CA138$W,2))
colnames(est_W_2016CA138)=names(sel_features_2016CA138)

est_W_2016CA138


# propagation_2016CO246_new_v2 <--- START WITH THIS
df_2016CO246=read.table("data_asteroids/propagation_2016CO246_new_v2.txt",header = T)
df_2016CO246=trans_theta(df_2016CO246)
unique(df_2016CO246$type)
df_2016CO246$type[which(df_2016CO246$type==100)]=10

plot(df_2016CO246$theta,type='p',col=df_2016CO246$type+2)

# propagation_2016HO3_new_v2
df_2016HO3=read.table("data_asteroids/propagation_2016HO3_new_v2.txt",header = T)
df_2016HO3=trans_theta(df_2016HO3)
plot(df_2016HO3$theta,type='p',col=df_2016HO3$type+2)

# propagation_2019GM1_new_v2
df_2019GM1=read.table("data_asteroids/propagation_2019GM1_new_v2.txt",header = T)
df_2019GM1=trans_theta(df_2019GM1)
unique(df_2019GM1$type)
plot(df_2019GM1$theta,type='p',col=df_2019GM1$type+2)


# propagation_2020PN1_new_v2 <--- THIS IS ALSO GOOD STARTING POINT
df_2020PN1=read.table("data_asteroids/propagation_2020PN1_new_v2.txt",header = T)
df_2020PN1=trans_theta(df_2020PN1)
unique(df_2020PN1$type)
df_2020PN1$type[which(df_2020PN1$type==100)]=10

plot(df_2020PN1$theta,type='p',col=df_2020PN1$type+2)

# propagation_2020PP1_new_v2
df_2020PP1=read.table("data_asteroids/propagation_2020PP1_new_v2.txt",header = T)
df_2020PP1=trans_theta(df_2020PP1)
plot(df_2020PP1$theta,type='p',col=df_2020PP1$type+2)


# propagation_164207_new_v2
df_164207=read.table("data_asteroids/propagation_164207_new_v2.txt",header = T)
df_164207=trans_theta(df_164207)
plot(df_164207$theta,type='p',col=df_164207$type+2)


