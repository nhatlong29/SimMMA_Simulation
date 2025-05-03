
# Change gamma[7] + gamma[8]A (interaction between Z&M) when Y(s2') = gamma0 + gamma1Z + gamma2M + gamma3ZM
# Change rho (cor between Z(0) and Z(1) in the binary setting)
# Change AZ (by proposing correlation between Z(1) and Z(0) in continuous setting)
# Change ZM: beta[4] (interaction between A&Z)
# Change ZY: gamma[5]A + gamma[7]M + gamma[8]AM

# Second setting: A binary; Z&M&Y continuous

library(dplyr)

generation = function(n, alpha0, alpha1, beta, gamma, eps) {
  
  #Z = expit(alpha[1] + alpha[2]*A)
  #M = beta[1] + beta[2]A + beta[3]*Z + beta[4]*A*Z + epsM
  #Y = gamma[1] + gamma[2]*A + gamma[3]*Z + gamma[4]*M + gamma[5]*A*Z + gamma[6]*A*M + gamma[7]*Z*M + gamma[8]*A*Z*M + epsY
  
  # Generate (Z0,Z1)
  Z0 = rnorm(n, alpha0, 1)
  Z1 = Z0 + alpha1
  # Generate (T0,T1)
  T0 = rnorm(n, alpha0, 1)
  T1 = rnorm(n, alpha0 + mean(alpha1), 1)

  # Generate M(1), M(0), M(1,Z(0)), M(1,Z(1)), M(1,T(1))
  epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
  M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
  M0 = beta[1] + beta[3]*Z0 + epsM
  M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
  M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM

  # Generate Ys1, Ys2, Ys3
  epsY = rnorm(n=n, mean=eps[[3]][1], sd=eps[[3]][2])
  Ys1 = gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1 + epsY
  Ys2 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1 + epsY
  Ys3 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0 + epsY
  
  # Generate Ys1'; Ys2'; Ys2"; Ys3"
  Ys1.prime = gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1 + epsY
  Ys2.prime = gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1 + epsY
  Ys2.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1 + epsY
  Ys3.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0 + epsY

  return(colMeans(data.frame(Z0,Z1,M1,M0,Ys1,Ys2,Ys3,Ys1.prime,Ys2.prime,Ys2.dprime,Ys3.dprime)))
}

# Example
alpha0 = 0.5
set.seed(1234)
alpha1 = 2
beta = 1*c(0.5, 0.5, 1.5, 2)
gamma = 1*c(1.5, 1.5, 1.5, 1.5, 1, 1, 0, 1)
eps = list(c(NaN, NaN), c(0, 1), c(0, 1))

## 1. Changing gamma[7]

gamma7 = c(-4, -3, -0.5, 1, 1.5)
result = data.frame(value = c(), AZY_ne_re = c(), AZMY_ne_re = c())

for(gamma7 in gamma7){
  gamma[7] = gamma7
  res = generation(n=1e7, alpha0=alpha0, alpha1=alpha1, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(value = gamma7, AZY_ne_re = (res[6] - res[5]) - (res[9] - res[8]),
                          AZMY_ne_re = (res[7] - res[6]) - (res[11] - res[10]))
              )
}

# 2. Z(1) - Z(0) is not constant anymore

mu1 = c(-3, -1, 0, 4, 5) # mean of alpha1
result = data.frame(AZY_ne_re = c(), AZMY_ne_re = c())

for(mu1 in mu1){
  alpha1 = rnorm(n, mu1, 10) # large variance
  res = generation(n=1e7, alpha0=alpha0, alpha1=alpha1, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
            data.frame(AZY_ne_re = (res[6] - res[5]) - (res[9] - res[8]),
                        AZMY_ne_re = (res[7] - res[6]) - (res[11] - res[10]))
            )
}

# 3. Changing beta[4]

beta4 = c(-3, -1, 0, 5, 7)
result = data.frame(value = c(), AZY_ne_re = c(), AZMY_ne_re = c())

for(beta4 in beta4){
  beta[4] = beta4
  res = generation(n=1e7, alpha0=alpha0, alpha1=alpha1, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(value = beta4, AZY_ne_re = (res[6] - res[5]) - (res[9] - res[8]),
                          AZMY_ne_re = (res[7] - res[6]) - (res[11] - res[10]))
              )
}

