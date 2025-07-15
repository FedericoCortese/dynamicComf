
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

# Transform theta
df_2014OL339=trans_theta(df_2014OL339)

# Select relevant variables
df_2014OL339=df_2014OL339[,c("t","a","theta","type")]
df_2014OL339$type[which(df_2014OL339$type==100)]=10

# Compute features
features_2014OL339=compute_feat(df_2014OL339,l_short=5,l_long=100,
                          tt_thres_maxmin=2.5)


# Remove NAs
features_2014OL339=features_2014OL339[complete.cases(features_2014OL339),]

# Extract ground truth and time
gt_2014OL339=features_2014OL339$type
time=features_2014OL339$t

# Remove ground truth 
features_2014OL339=subset(features_2014OL339,select=-c(type,t))

apply(features_2014OL339,2,summary)

# Remove I_TP and I_HS (zero variability)
#features_2014OL339=subset(features_2014OL339,select=-c(I_TP,I_HS))
sel_features_2014OL339=subset(features_2014OL339,
                              select=-c(
                                 theta,a
                                         ,dtheta,
                                        I_TP, I_HS,
                                        I_QS,I_CP, mean_osc
                                        ))

head(sel_features_2014OL339)

# plot(features_2014OL339$dtheta,col=gt_2014OL339+2)
# plot(features_2014OL339$sd_dtheta,col=gt_2014OL339+2)
# plot(features_2014OL339$mean_osc,col=gt_2014OL339+2)
# plot(features_2014OL339$I_QS,col=gt_2014OL339+2)
# plot(features_2014OL339$sd_a,col=gt_2014OL339+2)

source("Utils_sparse_robust_2.R")

# Select lambda, K, and c
cv_2014OL339=cv_robust_sparse_jump(
    Y=sel_features_2014OL339,
    true_states=gt_2014OL339,
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
gap_2014OL339=gap_robust_sparse_jump(
    Y=sel_features_2014OL339,
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

fit_2014OL339=robust_sparse_jump(
    Y=sel_features_2014OL339,
    K=3,
    zeta0=.3,
    lambda=.7,
    c=10,
    knn=10,
    M=NULL,
    n_init=3,
    verbose=T,
    tol=1e-16
)

est_s_2014OL339=fit_2014OL339$s
est_s_2014OL339[fit_2014OL339$v==0]=0
plot(features_2014OL339$theta,col=est_s_2014OL339+1)
plot(features_2014OL339$a,col=est_s_2014OL339+1)

est_W_2014OL339= data.frame(round(fit_2014OL339$W,2))
colnames(est_W_2014OL339)=names(sel_features_2014OL339)

est_W_2014OL339

# propagation_2015SO2_new_v2
df_2015SO2=read.table("data_asteroids/propagation_2015SO2_new_v2.txt",header = T)
df_2015SO2=trans_theta(df_2015SO2)
unique(df_2015SO2$type)
plot(df_2015SO2$theta,type='p',col=df_2015SO2$type+2)

# propagation_2015XX169_new_v2
df_2015XX169=read.table("data_asteroids/propagation_2015XX169_new_v2.txt",header = T)
df_2015XX169=trans_theta(df_2015XX169)
unique(df_2015XX169$type)
plot(df_2015XX169$theta,type='p',col=df_2015XX169$type+2)

# propagation_2016CA138_new_v2
df_2016CA138=read.table("data_asteroids/propagation_2016CA138_new_v2.txt",header = T)
df_2016CA138=trans_theta(df_2016CA138)
unique(df_2016CA138$type)
df_2016CA138$type[which(df_2016CA138$type==100)]=10
plot(df_2016CA138$theta,type='p',col=df_2016CA138$type+2)

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


