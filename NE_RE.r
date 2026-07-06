library(parallel)

expit = function(x) {exp(x)/(1+exp(x))}

generation = function(n, rho, alpha, beta, gamma, eps, mode.z, mode.y, mode.m) {
  
  if (mode.z == 'con') {
    # Generate (Z0,Z1)
    muZ = c(alpha[1],alpha[1]+alpha[2])
    covZ = matrix(c(eps[[1]][2], rho, rho, eps[[1]][2]), nrow=2)
    matZ = mvrnorm(n=n, mu=muZ, Sigma=covZ)
    Z0 = matZ[,1]
    Z1 = matZ[,2]
    # Generate (T0,T1)
    T0 = rnorm(n, muZ[1], 1)
    T1 = rnorm(n, muZ[2], 1)
  } else if (mode.z == 'bi') {
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
  } else (stop('ERROR: mode Z is not available'))
  

  # Generate M(1), M(0), M(1,Z(0)), M(1,Z(1)), M(1,T(1))
  if (mode.m == 'con') {
    epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
    M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
    M0 = beta[1] + beta[3]*Z0 + epsM
    M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
    M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
    # RT 3.0
    M1.3 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
    M0.3 = beta[1] + beta[3]*T0 + epsM
    M1Z0.3 = beta[1] + beta[2] + beta[3]*T0 + beta[4]*T0 + epsM
    M1T1.3 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
  } else if (mode.m == 'bi') {
    M1 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1))
    M0 = rbinom(n, 1, expit(beta[1] + beta[3]*Z0))
    M1Z0 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0))
    M1T1 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1))
    # RT 3.0
    M1.3 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1))
    M0.3 = rbinom(n, 1, expit(beta[1] + beta[3]*T0))
    M1Z0.3 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*T0 + beta[4]*T0))
    M1T1.3 = rbinom(n, 1, expit(beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1))
  } else (stop('ERROR: mode M is not available'))

  if (mode.y=='con') {
    # Generate Ys0, Ys1, Ys2, Ys3, Ys4 (NE)
    epsY = rnorm(n=n, mean=eps[[3]][1], sd=eps[[3]][2])
    Ys0 = gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1 + gamma[5]*Z1 + gamma[6]*M1 + gamma[7]*Z1*M1 + gamma[8]*1*Z1*M1 + epsY
    Ys1 = gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1 + epsY
    Ys2 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1 + epsY
    Ys3 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0 + epsY
    Ys4 = gamma[1] + gamma[3]*Z0 + gamma[4]*M0 + gamma[7]*Z0*M0 + epsY
    # Generate Ys0'; Ys1'; Ys2'; Ys2"; Ys3"; Ys4" (RE.1.0 & 2.0)
    Ys0.prime = gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[5]*Z1 + gamma[6]*M1T1 + gamma[7]*Z1*M1T1 + gamma[8]*1*Z1*M1T1 + epsY
    Ys1.prime = gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1 + epsY
    Ys2.prime = gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1 + epsY
    Ys2.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1 + epsY
    Ys3.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0 + epsY
    Ys4.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M0 + gamma[7]*T0*M0 + epsY
    # Generate Ys0.star, Ys1.star, Ys2.star, Ys3.star, Ys4.star (RE.3.0)
    Ys0.star = gamma[1] + gamma[2] + gamma[3]*T1 + gamma[4]*M1.3 + gamma[5]*T1 + gamma[6]*M1.3 + gamma[7]*T1*M1.3 + gamma[8]*1*T1*M1.3 + epsY
    Ys1.star = gamma[1] + gamma[3]*T1 + gamma[4]*M1.3 + gamma[7]*T1*M1.3 + epsY
    Ys2.star = gamma[1] + gamma[3]*T0 + gamma[4]*M1.3 + gamma[7]*T0*M1.3 + epsY
    Ys3.star = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0.3 + gamma[7]*T0*M1Z0.3 + epsY
    Ys4.star = gamma[1] + gamma[3]*T0 + gamma[4]*M0.3 + gamma[7]*T0*M0.3 + epsY
  } else if (mode.y=='bi') {
    # Generate Ys0, Ys1, Ys2, Ys3, Ys4 (NE)
    Ys0 = rbinom(n, 1, expit(gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1 + gamma[5]*Z1 + gamma[6]*M1 + gamma[7]*Z1*M1 + gamma[8]*1*Z1*M1))
    Ys1 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1))
    Ys2 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1))
    Ys3 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0))
    Ys4 = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M0 + gamma[7]*Z0*M0))
    # Generate Ys0'; Ys1'; Ys2'; Ys2"; Ys3"; Ys4" (RE.1.0 & 2.0)
    Ys0.prime = rbinom(n, 1, expit(gamma[1] + gamma[2] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[5]*Z1 + gamma[6]*M1T1 + gamma[7]*Z1*M1T1))
    Ys1.prime = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1))
    Ys2.prime = rbinom(n, 1, expit(gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1))
    Ys2.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1))
    Ys3.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0))
    Ys4.dprime = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M0 + gamma[7]*T0*M0))
    # Generate Ys0.star, Ys1.star, Ys2.star, Ys3.star, Ys4.star (RE.3.0)
    Ys0.star = rbinom(n, 1, expit(gamma[1] + gamma[2] + gamma[3]*T1 + gamma[4]*M1.3 + gamma[5]*T1 + gamma[6]*M1.3 + gamma[7]*T1*M1.3 + gamma[8]*1*T1*M1.3))
    Ys1.star = rbinom(n, 1, expit(gamma[1] + gamma[3]*T1 + gamma[4]*M1.3 + gamma[7]*T1*M1.3))
    Ys2.star = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1.3 + gamma[7]*T0*M1.3))
    Ys3.star = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0.3 + gamma[7]*T0*M1Z0.3))
    Ys4.star = rbinom(n, 1, expit(gamma[1] + gamma[3]*T0 + gamma[4]*M0.3 + gamma[7]*T0*M0.3))
  } else (stop('ERROR: mode Y is not available'))

  return(colMeans(data.frame(Z0, Z1, M1, M0, 
                            Ys0, Ys1, Ys2, Ys3, Ys4,
                            Ys0.prime, Ys1.prime, Ys2.prime, Ys2.dprime, Ys3.dprime, Ys4.dprime,
                            Ys0.star, Ys1.star, Ys2.star, Ys3.star, Ys4.star)))
}

ranking = function(data) {
  row.names(data) = NULL
  ne = sapply(1:nrow(data), function(i) rank(data[i,c('P1','P2','P3','P4')],na.last = "keep"))
  ne = t(ne)
  re_1.0 = sapply(1:nrow(data), function(i) rank(data[i,c('P1_1.0','P2_1.0','P3_1.0','P4_1.0')], na.last = "keep"))
  re_1.0 = t(re_1.0)
  re_2.0 = sapply(1:nrow(data), function(i) rank(data[i,c('P1_2.0','P2_2.0','P3_2.0','P4_2.0')], na.last = "keep"))
  re_2.0 = t(re_2.0)
  re_3.0 = sapply(1:nrow(data), function(i) rank(data[i,c('P1_3.0','P2_3.0','P3_3.0','P4_3.0')], na.last = "keep"))
  re_3.0 = t(re_3.0)

  data[,c('P1','P2','P3','P4')] = ne[,c('P1','P2','P3','P4')]
  data[,c('P1_1.0','P2_1.0','P3_1.0','P4_1.0')] = re_1.0[,c('P1_1.0','P2_1.0','P3_1.0','P4_1.0')]
  data[,c('P1_2.0','P2_2.0','P3_2.0','P4_2.0')] = re_2.0[,c('P1_2.0','P2_2.0','P3_2.0','P4_2.0')]
  data[,c('P1_3.0','P2_3.0','P3_3.0','P4_3.0')] = re_3.0[,c('P1_3.0','P2_3.0','P3_3.0','P4_3.0')]

  return(data[,c('mode.z','mode.m','mode.y',
                'rho', 'alpha_1', 'beta_2', 'gamma_2', 'gamma_3', 'gamma_6',
                'Ze_1.0','Ze_2.0','Ze_3.0',
                'P1','P2','P3','P4',
                'P1_1.0','P2_1.0','P3_1.0','P4_1.0',
                'P1_2.0','P2_2.0','P3_2.0','P4_2.0',
                'P1_3.0','P2_3.0','P3_3.0','P4_3.0')])
}
n = 1e7
rho = c(-0.2)
alpha1 = c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5)
beta2 = c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5)
gamma2 = c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5)
gamma3 = c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5)
gamma6 = c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5)
mode = c('con', 'bi')
params = expand.grid(mode.z = mode,
                    mode.m = mode,
                    mode.y = 'bi',
                    rho = rho,
                    alpha_1 = alpha1,
                    beta_2 = beta2,
                    gamma_2 = gamma2,
                    gamma_3 = gamma3,
                    gamma_6 = gamma6)

nworkers <- detectCores()
cl <- makeCluster(nworkers)
tasks <- c(1:nworkers)

worker <- function(task) {
    result <- c()
    indices <- split(c(1:nrow(params)), cut(c(1:nrow(params)), breaks = nworkers, labels = FALSE))[[task]]

    for (i in indices) {
        set.seed(1234)
        r = params[i,"rho"]
        a1 = params[i,"alpha_1"]
        b2 = params[i,"beta_2"]
        g2 = params[i,"gamma_2"]
        g3 = params[i,"gamma_3"]
        g6 = params[i,"gamma_6"]
        mode.z = params[i,"mode.z"]
        mode.y = params[i,"mode.y"]
        mode.m = params[i,"mode.m"]
        res = generation(n = n, 
                        rho = r, 
                        alpha = c(0.5, a1), 
                        beta = c(0.5, 0.5, b2, 2), 
                        gamma = c(1.5, 1.5, g2, g3, 1, 1, g6, 1), 
                        eps = list(c(0, 1), c(0, 1), c(0, 1)),
                        mode.z = mode.z,
                        mode.y = mode.y,
                        mode.m = mode.m)
        result <- do.call(rbind, list(result,
            data.frame(mode.z = mode.z, mode.y = mode.y, mode.m = mode.m,
                        rho = r, alpha_1 = a1, beta_2 = b2, gamma_2 = g2, gamma_3 = g3, gamma_6 = g6, 
                        Ys0 = res[5], Ys1 = res[6], Ys2 = res[7], Ys3 = res[8], Ys4 = res[9],
                        Ys0.prime = res[10], Ys1.prime = res[11], Ys2.prime = res[12], Ys2.dprime = res[13], Ys3.dprime = res[14], Ys4.dprime = res[15],
                        Ys0.star = res[16], Ys1.star = res[17], Ys2.star = res[18], Ys3.star = res[19], Ys4.star = res[20],
                        P1 = res[5] - res[6], P2 = res[6] - res[7], P3 = res[7] - res[8], P4 = res[8] - res[9],
                        P1_1.0 = res[5] - res[6], P2_1.0 = res[11] - res[12], P3_1.0 = res[13] - res[14], P4_1.0 = res[8] - res[9],
                        Ze_1.0 = res[6] - res[11] + res[12] - res[13] + res[14] - res[8],

                        diff_P2_1.0 = res[6] - res[7] - res[11] + res[12],
                        diff_PZ_1.0 = abs( res[6] - res[8] - res[11] + res[14]),
                        uppb_P2_z_1.0 = max(c(abs( max(c(-1,res[6]-res[8]-1)) - (res[11]-res[12])), abs( min(c(1,res[6]-res[8]+1)) - (res[11]-res[12])))),
                        uppb_P2_1_1.0 = max(c(abs((res[6]-1) - (res[11]-res[12])),abs(res[6] - (res[11]-res[12])))),
                        uppb_P2_bo_1.0 = 2 + min(c(res[6],1-res[11],res[12])) - max(c(res[6],1-res[11],res[12])) - abs(res[12]-(1-res[11])),
                        uppb_P2_bn_1.0 = 2 + min(c(res[6],1-res[11],res[12])) - max(c(res[6],1-res[11],res[12])) - max( abs(res[12]-(1-res[11])), abs(res[6]-(1-res[11])) ),

                        uppb_P3_z_1.0 = max(c(abs( max(c(-1,res[6]-res[8]-1)) - (res[13]-res[14])), abs( min(c(1,res[6]-res[8]+1)) - (res[13]-res[14])))),
                        uppb_P3_1_1.0 = max(c(abs((1-res[8]) - (res[13]-res[14])),abs(- res[8] - (res[13]-res[14])))),
                        uppb_P3_bo_1.0 = 2 + min(c(res[8],res[13],1-res[14])) - max(c(res[8],res[13],1-res[14])) - abs(res[13]-(1-res[14])),
                        uppb_P3_bn_1.0 = 2 + min(c(res[8],res[13],1-res[14])) - max(c(res[8],res[13],1-res[14])) - max( abs(res[13]-(1-res[14])), abs(res[8]-(1-res[14])) ),

                        P1_2.0 = res[10] - res[11], P2_2.0 = res[11] - res[12], P3_2.0 = res[13] - res[14], P4_2.0 = res[14] - res[15],
                        Ze_2.0 = res[5] - res[10] + res[15] - res[9],
                        diff_P2_2.0 = res[6] - res[7] - res[11] + res[12],
                        
                        P1_3.0 = res[16] - res[17], P2_3.0 = res[17] - res[18], P3_3.0 = res[18] - res[19], P4_3.0 = res[19] - res[20],
                        Ze_3.0 = res[5] - res[16] + res[20] - res[9],
                        diff_P2_3.0 = res[6] - res[7] - res[17] + res[18]
                        )))
    }
    list(res = result)
}
clusterExport(cl, c("params", "n", "expit", "generation", "tasks", "worker", "nworkers"))

clusterEvalQ(cl, {
    library(MASS)
    library(dplyr)
    NULL
})
results <- parLapplyLB(cl, tasks, worker)
stopCluster(cl)

# Combine results from all workers
final_result <- do.call(rbind, lapply(results, function(x) x$res))

write.csv(final_result, paste0("res_rneg0.2_ybi.csv"))
#write.csv(ranking(final_result), paste0("r_rneg0.2_ybi.csv"))