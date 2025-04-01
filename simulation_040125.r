# 1. Generate data from a specific model in the whole class of NPSEM-IE

library(dplyr)

set.seed(1234)

# A

piA = 0.7
nA_range = seq(from = 2000, to = 3500, by = 50)

final_df = data.frame(n = c(), AY_npsem = c(), AY_poe = c(),
                        AZY_npsem = c(), AZY_poe = c(),
                        AZMY_npsem = c(), AZMY_poe = c(),
                        AMY_npsem = c(), AMY_poe = c())

for(nA in nA_range){
    A = rbinom(nA, 1, piA) # Generate A
    
    epsZ = rnorm(n = nA, mean = 0, sd = 1); alp0 = 0.2; alp1 = 0.3
    Z = alp0 + alp1*A + epsZ # Generate Z

    alp0_est = lm(Z~A)$coefficient[1]
    alp1_est = lm(Z~A)$coefficient[2]

    epsM = rnorm(n = nA, mean = 0, sd = 1); beta0 = 1; beta1 = 2; beta2 = 3
    M = beta0 + beta1*A + beta2*Z + epsM # Generate M

    beta0_est = lm(M ~ A + Z)$coefficient[1]
    beta1_est = lm(M ~ A + Z)$coefficient[2]
    beta2_est = lm(M ~ A + Z)$coefficient[3]

    epsY = rnorm(n = nA, mean = 0, sd = 1)
    gamma0 = 1; gamma1 = 2; gamma2 = 3; gamma3 = 4
    Y = gamma0 + gamma1*A + gamma2*Z + gamma3*M + epsY # Generate Y

    gamma0_est = lm(Y ~ A + Z + M)$coefficient[1]
    gamma1_est = lm(Y ~ A + Z + M)$coefficient[2]
    gamma2_est = lm(Y ~ A + Z + M)$coefficient[3]
    gamma3_est = lm(Y ~ A + Z + M)$coefficient[4]

    E0 = gamma0 + gamma1 + gamma2*(alp0 + alp1) + 
        gamma3*(beta0 + beta1) + gamma3*beta2*(alp0 + alp1)

    E1 = gamma0 + gamma2*(alp0 + alp1) + 
        gamma3*(beta0 + beta1) + gamma3*beta2*(alp0 + alp1)

    fun.s2 = function(z, m) (gamma0 + gamma2*z + gamma3*m)*dnorm(z,alp0,1)*dnorm(m,beta0 + beta1 + beta2*(z+alp1),1)
    llim = -Inf
    ulim = Inf
    E2 = integrate(function(m) { 
                sapply(m, function(m) {
                        integrate(function(z) fun.s2(z,m), llim, ulim)$value})
        }, llim, ulim)$value

    E3 = gamma0 + gamma2*alp0 + gamma3*(beta0 + beta1 + beta2*alp0)
    E4 = gamma0 + gamma2*alp0 + gamma3*(beta0 + beta2*alp0)

    result = data.frame(n = nA,
                        AY_true = gamma1,
                        AY_poe = gamma1_est,
                        AZY_true = E1 - E2,
                        AZY_poe = alp1_est * gamma2_est,
                        AZMY_true = E2 - E3,
                        AZMY_poe = alp1_est*beta2_est*gamma3_est,
                        AMY_true = beta1 * gamma3,
                        AMY_poe = beta1_est * gamma3_est
                        )

    final_df = final_df %>% bind_rows(., result)

}


final_df %>% View()


