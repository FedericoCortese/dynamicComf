source("Utils_fuzzyJM_2.R")

prv=sim_data_stud_t(seed=123,
                             TT=1000,
                             P=10,
                             Pcat=NULL,
                             Ktrue=2,
                             mu=1.5,
                             rho=0,
                             nu=4,
                             phi=.8,
                             pers=.99)

Y=prv$SimData

temp=fuzzy_jump_cpp(Y, 
                    K=2, 
                    lambda = 0.05, 
                    m      = 1.01,
                    max_iter = 20, 
                    n_init   = 10, 
                    tol      = 1e-8, 
                    verbose  = FALSE
                                 
)

temp_par=fuzzy_jump_cpp_parallel(Y, 
                                    K=2, 
                                    lambda = 0.05, 
                                    m      = 1.01,
                                    max_iter = 20, 
                                    n_init   = 10, 
                                    tol      = 1e-8, 
                                    verbose  = FALSE
                                 )

table(temp$MAP,prv$mchain)
table(temp_par$MAP,prv$mchain)
