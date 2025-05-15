source("Utils_fuzzyJM_2.R")

prv=sim_data_stud_t(seed=123,
                             TT=1000,
                             P=10,
                             Pcat=2,
                             Ktrue=3,
                             mu=1.5,
                             rho=0,
                             nu=4,
                             phi=.8,
                             pers=.95)

Y=prv$SimData

fuzzy_jump_coord_par(Y, 
                                 K=3, 
                                 lambda=.1, 
                                 m=1.01,
                                 max_iter=8, 
                                 n_init=10, tol=1e-16, 
                                 verbose=FALSE,
                                 n_cores=3
                                 
)
