# Change gamma[7] + gamma[8]A (interaction between Z&M) when Y(s2') = gamma0 + gamma1Z + gamma2M + gamma3ZM
# Change rho (cor between Z(0) and Z(1) in the binary setting)
# Change AZ (by proposing correlation between Z(1) and Z(0) in continuous setting)
# Change ZM: beta[4] (interaction between A&Z)
# Change ZY: gamma[5]A + gamma[7]M + gamma[8]AM

# First setting: A&Z binary; M&Y continuous

library(dplyr, warn.conflicts = F)
library(kableExtra, warn.conflicts = F)
library(ggplot2)
library(tidyr)

expit = function(x) {exp(x)/(1+exp(x))}

generation = function(n, rho, alpha, beta, gamma, eps) {
  
  #Z = expit(alpha[1] + alpha[2]*A)
  #M = beta[1] + beta[2]A + beta[3]*Z + beta[4]*A*Z + epsM
  #Y = gamma[1] + gamma[2]*A + gamma[3]*Z + gamma[4]*M + gamma[5]*A*Z + gamma[6]*A*M + gamma[7]*Z*M + gamma[8]*A*Z*M + epsY 

  mu0 = expit(alpha[1])
  mu1 = expit(alpha[1] + alpha[2])

  # Generate (Z0,Z1)
  Z0 = rbinom(n, 1, mu0)
  p.z0.1 = mu1 + rho*sqrt((1-mu0)*mu1*(1-mu1)/mu0) # -sqrt((mu0*mu1)/((1-mu0)*(1-mu1))) < rho < sqrt(mu0*(1-mu1)/(mu1*(1-mu0)))
  p.z0.0 = mu1 - rho*sqrt(mu0*mu1*(1-mu1)/(1-mu0)) # -sqrt((1-mu0)*(1-mu1)/(mu0*mu1)) < rho < sqrt(mu1*(1-mu0)/(mu0*(1-mu1)))
  Z1 = rbinom(n,1,Z0*p.z0.1 + (1-Z0)*p.z0.0) # Z0*P[Z(1) =1|Z(0)=1] + (1-Z0)*P[Z(1) =1|Z(0)=0]

  # Generate (T0,T1)
  T0 = rbinom(n, 1, mu0)
  T1 = rbinom(n, 1, mu1)

  # Generate M(1), M(0), M(1,Z(0)), M(1,Z(1)), M(1,T(1))
  epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
  M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
  M0 = beta[1] + beta[3]*Z0 + epsM
  M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
  M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
  
  # Generate Ys1, Ys2, Ys3
  epsY = rnorm(n=n, mean=eps[[3]][1], sd=eps[[3]][2])
  Ys0 = gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1 + gamma[5]*Z1 +
        gamma[6]*M1 + gamma[7]*Z1*M1 + gamma[8]*1*Z1*M1 + epsY
  Ys1 = gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1 + epsY
  Ys2 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1 + epsY
  Ys3 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0 + epsY
  Ys4 = gamma[1] + gamma[3]*Z0 + gamma[4]*M0 + gamma[7]*Z0*M0 + epsY
  
  # Generate Ys1'; Ys2'; Ys2"; Ys3"
  Ys1.prime = gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1 + epsY
  Ys2.prime = gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1 + epsY
  Ys2.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1 + epsY
  Ys3.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0 + epsY

  # Calculate natural and recanting-twin effects
  # AZY_ne = Ys2 - Ys1 = gamma[3]*(Z0-Z1) + gamma[7]*M1*(Z0-Z1)
  # AZY_re = Ys2.prime - Ys1.prime = gamma[3]*(Z0-Z1) + gamma[7]*M1T1*(Z0-Z1)
  # AZMY_ne = Ys3 - Ys2 = gamma[4]*(M1Z0-M1) + gamma[7]*Z0*(M1Z0-M1)
  # AZMY_re = Ys3.dprime - Ys2.dprime = gamma[4]*(M1Z0-M1) + gamma[7]*T0*(M1Z0-M1)
  # AMY = Ys4 - Ys3 = gamma[4]*(M0-M1Z0) + gamma[7]*Z0*(M0-M1Z0)
  # AY = Ys1 - Ys0 = - gamma[2] - gamma[5]*Z1 - gamma[6]*M1 - gamma[8]*1*Z1*M1

  # Positive AZY; AZMY -> Pr(AZY>0)...
  AZY_ne_pos = ifelse(Ys2 - Ys1 > 0, 1, 0)
  AZY_re_pos = ifelse(Ys2.prime - Ys1.prime > 0, 1, 0)
  AZMY_ne_pos = ifelse(Ys3 - Ys2 > 0, 1, 0)
  AZMY_re_pos = ifelse(Ys3.dprime - Ys2.dprime > 0, 1, 0)


  return(colMeans(data.frame(Ys0,Ys1,Ys2,Ys3,Ys4,Ys1.prime,Ys2.prime,Ys2.dprime,Ys3.dprime)))
  }

# parameters original settings
alpha = 1*c(0.5, -0.5)
beta = 1*c(0.5, 0.5, 1.5, 2)
gamma = 1*c(1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1)
rho = 0.75
eps = list(c(NaN, NaN), c(0, 1), c(0, 1))

# Changing alpha[2]

alpha2.1 = c(-1.5, -0.5, 0.5, 1)
alpha.sim = alpha
result = data.frame(alpha1 =c(),
                    rho = c(),
                    rho.lower = c(),
                    rho.upper = c(),
                    AZY_ne = c(), 
                    AZY_re = c(),
                    AZMY_ne = c(), 
                    AZMY_re = c(),
                    AY = c(),
                    AMY = c())

set.seed(1234)

for(alpha2.sim in alpha2.1){

  alpha.sim[2] = alpha2.sim
  mu0 = expit(alpha.sim[1])
  mu1 = expit(alpha.sim[1] + alpha.sim[2])
  rho1.upper = min(sqrt(mu0*(1-mu1)/(mu1*(1-mu0))), sqrt(mu1*(1-mu0)/(mu0*(1-mu1)))) #rho changes because alpha changes
  rho1.lower = max(-sqrt((mu0*mu1)/((1-mu0)*(1-mu1))), -sqrt((1-mu0)*(1-mu1)/(mu0*mu1)))
  rho1 = 0.4
  res = generation(n=1e7, rho=rho1, alpha=alpha.sim, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(alpha1 = alpha2.sim,
                        rho = rho1,
                        rho.lower = rho1.lower,
                        rho.upper = rho1.upper,
                        AZY_ne = res["Ys2"] - res["Ys1"], 
                        AZY_re = res["Ys2.prime"] - res["Ys1.prime"],
                        AZMY_ne = res["Ys3"] - res["Ys2"],
                        AZMY_re = res["Ys3.dprime"] - res["Ys2.dprime"],
                        AY = res["Ys1"] - res["Ys0"],
                        AMY = res["Ys4"] - res["Ys3"])
              )
}

result %>% kable()

result %>%
  pivot_longer(., cols = c("AZY_ne", "AZY_re", "AZMY_ne", "AZMY_re", "AY", "AMY"),
              names_to = "path",
              values_to = "path_value") %>%
  ggplot(data = ., aes(x = alpha1, y = path_value, colour = path)) +
    geom_point(stat = "identity") +
    geom_line() +
    labs(y = "Value of specific path", x = expression(alpha[1])) +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())


# After choosing


# 1. Changing gamma[7]

gamma7.1 = c(-1, -0.5, 0.5, 1)
gamma.sim = gamma
result = data.frame(type =c(),
                    value = c(),
                    AZY_ne = c(), 
                    AZY_re = c(), 
                    AZMY_ne = c(), 
                    AZMY_re = c(),
                    AY = c(),
                    AMY = c())
set.seed(1234)

for(gamma7.sim in gamma7.1){

  gamma.sim[7] = gamma7.sim
  res = generation(n=1e7, rho=0.75, alpha=alpha, beta=beta, gamma=gamma.sim, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = 'gamma6',
                        value = gamma7.sim,
                        AZY_ne = res["Ys2"] - res["Ys1"], 
                        AZY_re = res["Ys2.prime"] - res["Ys1.prime"],
                        AZMY_ne = res["Ys3"] - res["Ys2"],
                        AZMY_re = res["Ys3.dprime"] - res["Ys2.dprime"],
                        AY = res["Ys1"] - res["Ys0"],
                        AMY = res["Ys4"] - res["Ys3"])
              )
}

result %>% kable(digits = 3)
result %>% kable(format = "latex", digits = 3)

result %>%
  pivot_longer(., cols = c("AZY_ne", "AZY_re", "AZMY_ne", "AZMY_re", "AY", "AMY"),
              names_to = "path",
              values_to = "path_value") %>%
  ggplot(data = ., aes(x = value, y = path_value, colour = path)) +
    geom_point(stat = "identity") +
    geom_line() +
    labs(y = "Value of specific path", x = expression(gamma[6])) +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())


# 2. Changing rho

rho.sim = c(-0.6, 0.3, 0.5, 0.75) # -sqrt((1-mu0)*(1-mu1)/(mu0*mu1)) = - 0.7788 < rho <= sqrt(mu1*(1-mu0)/(mu0*(1-mu1))) = 0.7788

result = data.frame(type = c(),
                    value = c(), 
                    AZY_ne = c(),
                    AZY_re = c(), 
                    AZMY_ne = c(), 
                    AZMY_re = c(),
                    AY = c(),
                    AMY = c())
set.seed(1234)

for(rho in rho.sim){
  res = generation(n=1e7, rho=rho, alpha=alpha, beta=beta, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = "rho",
                        value = rho, 
                        AZY_ne = res["Ys2"] - res["Ys1"], 
                        AZY_re = res["Ys2.prime"] - res["Ys1.prime"],
                        AZMY_ne = res["Ys3"] - res["Ys2"],
                        AZMY_re = res["Ys3.dprime"] - res["Ys2.dprime"],
                        AY = res["Ys1"] - res["Ys0"],
                        AMY = res["Ys4"] - res["Ys3"])
              )
}

result %>%
  pivot_longer(., cols = c("AZY_ne", "AZY_re", "AZMY_ne", "AZMY_re", "AY", "AMY"),
              names_to = "path",
              values_to = "path_value") %>%
  ggplot(data = ., aes(x = value, y = path_value, colour = path)) +
    geom_point(stat = "identity") +
    geom_line() +
    labs(y = "Value of specific path", x = expression(rho)) +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

result %>% kable()
result %>% kableExtra::kable(format = "latex", digits = 3)

# 3. Changing beta[3]


beta3.1 = c(-1.5, -1, 0, 1, 2)
beta.sim = c(0.5, 0.5, 1.5, 2)

result = data.frame(type = c(),
                    value = c(), 
                    AZY_ne = c(), 
                    AZY_re = c(), 
                    AZMY_ne = c(), 
                    AZMY_re = c(),
                    AY = c(),
                    AMY = c()
                    )
set.seed(1234)

for(beta3 in beta3.1){
  beta.sim[3] = beta3
  res = generation(n=1e7, rho=0.75, alpha=alpha, beta=beta.sim, gamma=gamma, eps=eps)
  result = result %>% 
    bind_rows(., 
              data.frame(type = "beta2",
                        value = beta3, 
                        AZY_ne = res["Ys2"] - res["Ys1"], 
                        AZY_re = res["Ys2.prime"] - res["Ys1.prime"],
                        AZMY_ne = res["Ys3"] - res["Ys2"],
                        AZMY_re = res["Ys3.dprime"] - res["Ys2.dprime"],
                        AY = res["Ys1"] - res["Ys0"],
                        AMY = res["Ys4"] - res["Ys3"])
              )
}

result %>%
  pivot_longer(., cols = c("AZY_ne", "AZY_re", "AZMY_ne", "AZMY_re", "AY", "AMY"),
              names_to = "path",
              values_to = "path_value") %>%
  ggplot(data = ., aes(x = value, y = path_value, colour = path)) +
    geom_point(stat = "identity") +
    geom_line() +
    labs(y = "Value of specific path", x = expression(beta[2])) +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())



result %>% kableExtra::kable()

