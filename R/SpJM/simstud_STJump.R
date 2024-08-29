# lambda in (0,1) step 0.1
# gamma in (0,0.5) step 0.05
lambda=seq(0,1,.1)
gamma=seq(0,0.5,.05)
seed=1:100
hp=expand.grid(lambda=lambda,gamma=gamma,seed=seed)

# First decide empirical study and then choose the simulation study parameters (T and M)

# First sim stud ----------------------------------------------------------
# TT=4, M=100
# 1) mu=2, sigma=0
# 2) mu=.5, sigma=.2


# Supplementary sim stud --------------------------------------------------
# TT=8, M=400, mu=2, sigma=0
# TT=2, M=49, mu=2, sigma=0