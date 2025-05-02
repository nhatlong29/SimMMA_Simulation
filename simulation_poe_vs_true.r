
library(dplyr)
set.seed(1234)

piA = 0.7
nA_range = seq(from = 2000, to = 3500, by = 50)

final_df = data.frame(n = c(), AY_true = c(), AY_poe = c(),
                        AZY_true = c(), AZY_poe = c(),
                        AZMY_true = c(), AZMY_poe = c(),
                        AMY_true = c(), AMY_poe = c())

for(nA in nA_range){
    A = rbinom(nA, 1, piA) # Generate A
    
    epsZ = rnorm(n = nA, mean = 0, sd = 1); alp0 = 0.2; alp1 = 0.3
    Z = alp0 + alp1*A + epsZ # Generate Z

    epsM = rnorm(n = nA, mean = 0, sd = 1); beta0 = 1; beta1 = 2; beta2 = 3
    M = beta0 + beta1*A + beta2*Z + epsM # Generate M

    epsY = rnorm(n = nA, mean = 0, sd = 1)
    gamma0 = 1; gamma1 = 2; gamma2 = 3; gamma3 = 4; gamma4 = 3.5
    Y = gamma0 + gamma1*A + gamma2*Z + gamma3*M + gamma4*Z*M + epsY # Generate Y

    E0 = gamma0 + gamma1 + gamma2*(alp0 + alp1) + 
        gamma3*(beta0 + beta1) + gamma3*beta2*(alp0 + alp1) + gamma4*(alp0+alp1)*(beta0 +beta1+beta2*(alp0+alp1))

    E1 = gamma0 + gamma2*(alp0 + alp1) + 
        gamma3*(beta0 + beta1) + gamma3*beta2*(alp0 + alp1) + gamma4*(alp0+alp1)*(beta0 +beta1 + beta2*(alp0+alp1))

    fun.s2 = function(z, m) (gamma0 + gamma2*z + gamma3*m + gamma4 *z*m)*dnorm(z,alp0,1)*dnorm(m,beta0 + beta1 + beta2*(z+alp1),1)
    llim = -Inf
    ulim = Inf
    E2 = integrate(function(m) { 
                sapply(m, function(m) {
                        integrate(function(z) fun.s2(z,m), llim, ulim)$value})
        }, llim, ulim)$value

    E3 = gamma0 + gamma2*alp0 + gamma3*(beta0 + beta1 + beta2*alp0) + gamma4*alp0*(beta0+beta1+beta2*alp0)
    E4 = gamma0 + gamma2*alp0 + gamma3*(beta0 + beta2*alp0) + gamma4*alp0*(beta0+beta2*alp0)

    result = data.frame(n = nA,
                        AY_true = gamma1,
                        AY_poe = gamma1,
                        AZY_true = E1 - E2,
                        AZY_poe = alp1 * gamma2,
                        AZMY_true = E2 - E3,
                        AZMY_poe = alp1*beta2*gamma3,
                        AMY_true = beta1 * gamma3,
                        AMY_poe = beta1 * gamma3
                        )

    final_df = final_df %>% bind_rows(., result)

}


final_df %>% View()