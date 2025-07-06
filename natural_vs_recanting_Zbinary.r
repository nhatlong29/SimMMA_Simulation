# Change gamma[7] + gamma[8]A (interaction between Z&M) when Y(s2') = gamma0 + gamma1Z + gamma2M + gamma3ZM
# Change rho (cor between Z(0) and Z(1) in the binary setting)
# Change AZ (by proposing correlation between Z(1) and Z(0) in continuous setting)
# Change ZM: beta[4] (interaction between A&Z)
# Change ZY: gamma[5]A + gamma[7]M + gamma[8]AM

# First setting: A&Z binary; M&Y continuous

library(dplyr)
expit = function(x) {exp(x)/(1+exp(x))}

generation = function(n=1e7, rho, alpha, beta, gamma, eps, mode.y, mode.m) {
  
  #Z = expit(alpha[1] + alpha[2]*A)
  #M = beta[1] + beta[2]A + beta[3]*Z + beta[4]*A*Z + epsM
  #Y = gamma[1] + gamma[2]*A + gamma[3]*Z + gamma[4]*M + gamma[5]*A*Z + gamma[6]*A*M + gamma[7]*Z*M + gamma[8]*A*Z*M + epsY 

  mu0 = expit(alpha[1])
  mu1 = expit(alpha[1] + alpha[2]) # mu1 should be < mu0 to make sure p.z0.1 <1, equivalent to alpha[2] <0
  
  # Generate (Z0,Z1)
  Z0 = rbinom(n, 1, mu0)
  p.z0.1 = mu1 + rho*sqrt((1-mu0)*mu1*(1-mu1)/mu0) # 0<= p.z0.1 <1, equivalent to alpha[1] >= -alpha[2]/2 and alpha[2] <0
  p.z0.0 = mu1 - rho*sqrt(mu0*mu1*(1-mu1)/(1-mu0)) # -sqrt((1-mu0)*(1-mu1)/(mu0*mu1)) < rho <= sqrt(mu1*(1-mu0)/(mu0*(1-mu1)))
  Z1 = rbinom(n,1,Z0*p.z0.1 + (1-Z0)*p.z0.0) # Z0*P[Z(1) =1|Z(0)=1] + (1-Z0)*P[Z(1) =1|Z(0)=0]

  # Generate (T0,T1)
  T0 = rbinom(n, 1, mu0)
  T1 = rbinom(n, 1, mu1)

  # Generate M(1), M(0), M(1,Z(0)), M(1,Z(1)), M(1,T(1))
  if (mode.m == 'linear') {
    epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
    M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
    M0 = beta[1] + beta[3]*Z0 + epsM
    M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
    M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
  } else if (mode.m == 'binary') {
    M1 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1))
    M0 = rbinom(n, 1, expit(beta[1] + beta[3]*Z0))
    M1Z0 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0))
    M1T1 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1))
  } else (stop('ERROR: mode M is not available'))
  
  if (mode.y=='linear') {
    # Generate Ys1, Ys2, Ys3
    epsY = rnorm(n=n, mean=eps[[3]][1], sd=eps[[3]][2])
    Ys0 = gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1 + gamma[5]*Z1 + gamma[6]*M1 + gamma[7]*Z1*M1 + gamma[8]*1*Z1*M1 + epsY
    Ys1 = gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1 + epsY
    Ys2 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1 + epsY
    Ys3 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0 + epsY
    Ys4 = gamma[1] + gamma[3]*Z0 + gamma[4]*M0 + gamma[7]*Z0*M0 + epsY

    # Generate Ys1'; Ys2'; Ys2"; Ys3"
    Ys1.prime = gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1 + epsY
    Ys2.prime = gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1 + epsY
    Ys2.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1 + epsY
    Ys3.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0 + epsY
  } else if (mode.y=='binary') {
    Ys0 = rbinom(n, 1, expit(gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1 + gamma[5]*Z1 + gamma[6]*M1 + gamma[7]*Z1*M1 + gamma[8]*1*Z1*M1))
    Ys1 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1))
    Ys2 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1))
    Ys3 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0))
    Ys4 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M0 + gamma[7]*Z0*M0))

    # Generate Ys1'; Ys2'; Ys2"; Ys3"
    Ys1.prime = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1))
    Ys2.prime = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1))
    Ys2.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1))
    Ys3.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0))
  } else (stop('ERROR: mode Y is not available'))
  # ositive AZY; AZMY -> Pr(AZY>0)...
  AZY_ne_pos = ifelse(Ys2 - Ys1 > 0, 1, 0)
  AZY_re_pos = ifelse(Ys2.prime - Ys1.prime > 0, 1, 0)
  AZMY_ne_pos = ifelse(Ys3 - Ys2 > 0, 1, 0)
  AZMY_re_pos = ifelse(Ys3.dprime - Ys2.dprime > 0, 1, 0)

  return(colMeans(data.frame(Z0,Z1,M1,M0,Ys1,Ys2,Ys3,Ys1.prime,Ys2.prime,Ys2.dprime,Ys3.dprime,Ys0,Ys4)))
  }

ranking = function(data) {
  row.names(data) = NULL
  data[,'p2xp3'] = data[,'AZY_re']*data[,'AZMY_re']
  ne = sapply(1:nrow(data), function(i) rank(data[i,c('AY','AMY','AZY_ne','AZMY_ne')],na.last = "keep"))
  ne = t(ne)
  re = sapply(1:nrow(data), function(i) rank(data[i,c('AY','AMY','AZY_re','AZMY_re')], na.last = "keep"))
  re = t(re)
  data[,c('AY_ne','AMY_ne','AZY_ne','AZMY_ne')] = ne[,c('AY','AMY','AZY_ne','AZMY_ne')]
  data[,c('AY_re','AMY_re','AZY_re','AZMY_re')] = re[,c('AY','AMY','AZY_re','AZMY_re')]
  return(data[,c('a1','g6','rho','beta2','Z_e','p2xp3','AY_ne','AMY_ne','AZY_ne','AZMY_ne','AY_re','AMY_re','AZY_re','AZMY_re')])
}


alpha2 = c(-1.5, -0.5, 0.5, 1.5)
gamma7 = c(-1.5, -1, 0, 1, 1.5)
rho = c(-0.75, -0.2, 0.2, 0.75)
beta3 = c(-1.5, -1, 1, 1.5)
mode.y = 'linear'
mode.m = 'binary'

result = data.frame(a1 = c(),
                 g6 = c(),
                 r = c(),
                 beta2 = c(),
                 AY = c(),
                 AMY = c(),
                 AZY_ne = c(), 
                 AZY_re = c(), 
                 AZMY_ne = c(), 
                 AZMY_re = c(),
                 Z_e = c())
for (a2 in alpha2) {
  for (g7 in gamma7) {
    for (r in rho) {
      for (b3 in beta3) {
        set.seed(1234)
        res = generation(n = 1e7, 
                        rho = r, 
                        alpha = c(0.5, a2), 
                        beta = c(0.5, 0.5, b3, 2), 
                        gamma = c(1.5, 1.5, 1.5, 1.5, 1, 1, g7, 1), 
                        eps = list(c(NaN, NaN), c(0, 1), c(0, 1)),
                        mode.y = mode.y,
                        mode.m = mode.m)
        result = result %>%
                bind_rows(.,
              data.frame(a1 = a2,
                        g6 = g7,
                        rho = r,
                        beta2 = b3,
                        AY = (res[5] - res[12]),
                        AMY = (res[13] - res[7]),
                        AZY_ne = (res[6] - res[5]),
                        AZY_re = (res[9] - res[8]),
                        AZMY_ne = (res[7] - res[6]),
                        AZMY_re = (res[11] - res[10]),
                        Z_e = (res[5] - res[8] + res[9] - res[10] + res[11] - res[7])
                        )
              )
      }
    }
  }
}
write.csv(result, paste0('SimMMA_Simulation/res_Zbi_M',mode.m,'_Y',mode.y,'.csv'))
r = ranking(result)
write.csv(r, paste0('SimMMA_Simulation/r_Zbi_M',mode.m,'_Y',mode.y,'.csv'))

if(FALSE) {
# parameters original settings
alpha = 1*c(0.5, -1.5)
beta = 1*c(0.5, 0.5, 1.5, 2)
gamma = 1*c(1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1)
rho = 0.75
eps = list(c(NaN, NaN), c(0, 1), c(0, 1))

result = data.frame(type =c(),
                    value = c(),
                    AZY_ne = c(), 
                    AZY_re = c(), 
                    AZY_ne_re = c(), 
                    AZMY_ne = c(), 
                    AZMY_re = c(),
                    AZMY_ne_re = c(),
                    AZY_ne_pos = c(),
                    AZY_re_pos = c(),
                    AZMY_ne_pos = c(),
                    AZMY_re_pos = c())

# 1. alpha[2] simulation
alpha1.l = c(-1.5 ,-1, -0.5, 0, 0.5, 1, 1.5)
alpha.sim = alpha
for(alpha1.sim in alpha1.l){
  set.seed(1234)
  alpha.sim[2] = alpha1.sim
  res = generation(n=1e7, rho=rho, alpha=alpha.sim, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = 'alpha1',
                         value = alpha1.sim,
                         AZY_ne = (res[6] - res[5]), 
                         AZY_re = (res[9] - res[8]), 
                         AZY_ne_re = ((res[6] - res[5]) - (res[9] - res[8])),
                         AZMY_ne = (res[7] - res[6]),
                         AZMY_re = (res[11] - res[10]),
                         AZMY_ne_re = ((res[7] - res[6]) - (res[11] - res[10])),
                         AZY_ne_pos = res[12],
                         AZY_re_pos = res[13],
                         AZMY_ne_pos = res[14],
                         AZMY_re_pos = res[15])
    )
}


# 1. gamma[7] simulation
gamma7.l = c(-2.5, -2.0, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
gamma7.l = c(-4, -3, -2, -1.5 ,-1, -0.5, 0, 0.5, 1, 1.5)
gamma.sim = gamma
for(gamma7.sim in gamma7.l){
  set.seed(1234)
  gamma.sim[7] = gamma7.sim
  res = generation(n=1e7, rho=rho, alpha=alpha, beta=beta, gamma=gamma.sim, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = 'gamma6',
                        value = gamma7.sim,
                        AZY_ne = (res[6] - res[5]), 
                        AZY_re = (res[9] - res[8]), 
                        AZY_ne_re = ((res[6] - res[5]) - (res[9] - res[8])),
                        AZMY_ne = (res[7] - res[6]),
                        AZMY_re = (res[11] - res[10]),
                        AZMY_ne_re = ((res[7] - res[6]) - (res[11] - res[10])),
                        AZY_ne_pos = res[12],
                        AZY_re_pos = res[13],
                        AZMY_ne_pos = res[14],
                        AZMY_re_pos = res[15])
              )
}

# 2. rho simulation
rho.l = seq(-0.5, 0.75, 0.25)
for(rho.sim in rho.l){
  set.seed(1234)
  res = generation(n=1e7, rho=rho.sim, alpha=alpha, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = 'rho',
                        value = rho.sim,
                        AZY_ne = (res[6] - res[5]), 
                        AZY_re = (res[9] - res[8]), 
                        AZY_ne_re = ((res[6] - res[5]) - (res[9] - res[8])),
                        AZMY_ne = (res[7] - res[6]),
                        AZMY_re = (res[11] - res[10]),
                        AZMY_ne_re = ((res[7] - res[6]) - (res[11] - res[10])),
                        AZY_ne_pos = res[12],
                        AZY_re_pos = res[13],
                        AZMY_ne_pos = res[14],
                        AZMY_re_pos = res[15])
              )
}

# 3. beta[4] simulation
beta4.l = c(-1, -0.5, 0, 0.5, 1)
beta.sim = beta
for(beta4.sim in beta4.l){
  set.seed(1234)
  beta.sim[4] = beta4.sim
  res = generation(n=1e7, rho=rho, alpha=alpha, beta=beta.sim, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = 'beta3',
                        value = beta4.sim,
                        AZY_ne = (res[6] - res[5]),
                        AZY_re = (res[9] - res[8]), 
                        AZY_ne_re = (res[6] - res[5]) - (res[9] - res[8]),
                        AZMY_ne = (res[7] - res[6]), 
                        AZMY_re = (res[11] - res[10]),
                        AZMY_ne_re = (res[7] - res[6]) - (res[11] - res[10]),
                        AZY_ne_pos = res[12],
                        AZY_re_pos = res[13],
                        AZMY_ne_pos = res[14],
                        AZMY_re_pos = res[15])
              )
}


# hard write
sink("SimMMA_Simulation/res_Zbinary.txt")
print("Simulation results for Z binary")
print("Original settings")
print(paste("alpha = ", paste(alpha, collapse = ",")))
print(paste("beta = ", paste(beta, collapse = ",")))
print(paste("gamma = ", paste(gamma, collapse = ",")))
print(paste("rho = ", rho))
print(paste("eps = ", paste(eps, collapse = ",")))
print("Simulation results")
print(result)
sink()
}