
# 1. Generate data from a specific model in the whole class of NPSEM-IE

library(dplyr)

set.seed(1234)

# A

piA = 0.7
nA_range = 1000:1200

E_df = data.frame(E0 = c(), E1 = c(), E2 = c(), E3 = c(), E4 = c())

for(nA in nA_range){
    A = rbinom(nA, 1, piA) # Generate A
    
    epsZ = rnorm(n = nA, mean = 0, sd = 1); alp0 = 0.2; alp1 = 0.3
    Z = alp0 + alp1*A + epsZ # Generate Z

    epsM = rnorm(n = nA, mean = 0, sd = 1); beta0 = 1; beta1 = 2; beta2 = 3
    M = beta0 + beta1*A + beta2*Z + epsM # Generate M

    epsY = rnorm(n = nA, mean = 0, sd = 1)
    gamma0 = 1; gamma1 = 2; gamma2 = 3; gamma3 = 4
    Y = gamma0 + gamma1*A + gamma2*Z + gamma3*M + epsY # Generate Y

    E0 = gamma0 + gamma1 + gamma2*(alp0 + alp1) + gamma3*(beta0 + beta1) + gamma3*beta2*(alp0 + alp1)
    E1 = gamma0 + gamma2*(alp0 + alp1) + gamma3*(beta0 + beta1) + gamma3*beta2*(alp0 + alp1)

    df = data.frame(A, Z, M, Y) %>%
        mutate(E2 = (gamma0 + gamma2*Z + gamma3*M) * dnorm(Z, alp0, 1) * dnorm(M, beta0 + beta1 + beta2*(Z+alp1), 1))
    E2 = sum(df$E2)

    E3 = gamma0 + gamma2*alp0 + gamma3*(beta0 + beta1) + gamma3*beta2*alp0
    E4 = gamma0 + gamma2*alp0 + gamma3*(beta0 + beta2*alp0)

    E_df = E_df %>%
        bind_rows(., data.frame(E0, E1, E2, E3, E4))

    # Product of coefficient
}

E_df = bind_rows(E_df, data.frame(E0 = 0, E1 = 1, E2=2,E3=3,E4=4))










