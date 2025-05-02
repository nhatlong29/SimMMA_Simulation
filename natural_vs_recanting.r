library(dplyr)
expit = function(x) {exp(x)/(1+exp(x))}

generation = function(n, rho, alpha, beta, gamma, eps, mode.y) {
    
    #Z = expit(alpha[1] + alpha[2]*A)
    #M = beta[1] + beta[2]A + beta[3]*Z + beta[4]*A*Z + epsM
    #Y = gamma[1] + gamma[2]*A + gamma[3]*Z + gamma[4]*M + gamma[5]*A*Z + gamma[6]*A*M + gamma[7]*Z*M + gamma[8]*A*Z*M + epsY 
    
    id = seq(1:n)
    mu0 = expit(alpha[1])
    mu1 = expit(alpha[1] + alpha[2])
    
    # Generate (Z0,Z1)
    Z0 = rbinom(n, 1, mu0)
    p.z0.1 = mu1 + rho*sqrt((1-mu0)*mu1*(1-mu1)/mu0)
    p.z0.0 = mu1 - rho*sqrt(mu0*mu1*(1-mu1)/(1-mu0))
    Z1 = rbinom(n,1,Z0*p.z0.1 + (1-Z0)*p.z0.0)
    
    # Generate (T0,T1)
    T0 = rbinom(n, 1, mu0)
    T1 = rbinom(n, 1, mu1)
    
    # Generate M(1), M(0), M(1,Z(0)), M(1,Z(1)), M(1,T(1))
    epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
    M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
    M0 = beta[1] + beta[3]*Z0 + epsM
    M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
    M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
    
    if (mode.y=='linear') {
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
    } else if (mode.y=='binary') {
      # Generate binary Y
      Ys1 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1))
      Ys2 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1))
      Ys3 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0))
      
      Ys1.prime = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1))
      Ys2.prime = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1))
      Ys2.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1))
      Ys3.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0))
    } else {stop('ERROR: mode Y is not available')}
    return(colMeans(data.frame(Z0,Z1,M1,M0,Ys1,Ys2,Ys3,Ys1.prime,Ys2.prime,Ys2.dprime,Ys3.dprime)))}

# Example
alpha = 1*c(0.5, 0.5)
beta = 1*c(0.5, 0.5, 1.5, 2)
gamma = 1*c(1.5, 1.5, 1.5, 1.5, 1, 1, 0, 1)
eps = list(c(NaN, NaN), c(0, 1), c(0, 1))
(res = generation(n=1e7, rho=0.75, alpha=alpha, beta=beta, gamma=gamma, eps=eps, mode.y='binary'))

# Natural effects E(phi2), E(phi3)
(ne = c(res[6] - res[5], res[7] - res[6]))

# Recanting twin effects E(phi2), E(phi3)
(re = c(res[9] - res[8], res[11] - res[10]))

