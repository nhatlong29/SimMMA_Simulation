library(dplyr)
expit = function(x) {exp(x)/(1+exp(x))}
#Y = gamma[1] + gamma[2]*A + gamma[3]*Z + gamma[4]*M + gamma[5]*A*Z + gamma[6]*A*M + gamma[7]*Z*M + gamma[8]*A*Z*M + epsY 
integration = function(rho, alpha, beta, gamma, eps) {
    z.val = c(0,1)
    m.range = c(-Inf,Inf)
    mu0 = expit(alpha[1])
    mu1 = expit(alpha[1] + alpha[2])
    #YS1
    integrate.s1 = function(z) {
        integrate(function(m) {
            (gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m)*dnorm(m,beta[1] + beta[2] + beta[3]*z + beta[4]*z,eps[[2]][2])*ifelse(z==0,1-expit(alpha[1]+alpha[2]),expit(alpha[1]+alpha[2]))
        }, lower = m.range[1], upper = m.range[2])$value
    }
    #YS2
    Y_a0.z.m = function(z,m) {gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m}
    fm_a1.zprime = function(zprime) {beta[1] + beta[2] + beta[3]*zprime + beta[4]*zprime}
    pz_a0 = function(z) {z*expit(alpha[1])+(1-z)*(1-expit(alpha[1]))}
    pz_z = function(z,zprime) {
        given_z1 = (rho * mu0^(1/2) * (1-mu0)^(1/2) * mu1^(1/2)*(1-mu1)^(1/2) + mu0*mu1)/mu0
        given_z0 = (mu1 - rho * mu0^(1/2) * (1-mu0)^(1/2) * mu1^(1/2)*(1-mu1)^(1/2) - mu0*mu1)/(1 - mu0)

        return(given_z1*z*zprime + (1-given_z1)*z*(1-zprime) + given_z0*(1-z)*zprime + (1-given_z0)*(1-z)*(1-zprime))
    }
    integrate.s2 = function(z,zprime) {
        integrate(function(m) {
            Y_a0.z.m(z,m) * dnorm(m,fm_a1.zprime(zprime),eps[[2]][2]) * pz_z(z,zprime) * pz_a0(z)
        }, lower = m.range[1], upper = m.range[2])$value
        }
    #YS3
    integrate.s3 = function(z) {
        integrate(function(m) {
            (gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m)*dnorm(m,beta[1] + beta[2] + beta[3]*z + beta[4]*z,eps[[2]][2])*ifelse(z==0,1-expit(alpha[1]),expit(alpha[1]))
        }, lower = m.range[1], upper = m.range[2])$value
    }
    #YS1.prime
    integrate.s1.prime = function(z,z.add) {
        integrate(function(m) {
            (gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m)*ifelse(z==0,1-expit(alpha[1]+alpha[2]),expit(alpha[1]+alpha[2]))*dnorm(m,beta[1] + beta[2] + beta[3]*z.add + beta[4]*z.add,eps[[2]][2])*ifelse(z.add==0,1-expit(alpha[1]+alpha[2]),expit(alpha[1]+alpha[2]))
        }, lower = m.range[1], upper = m.range[2])$value
    }
    #YS2.prime
    integrate.s2.prime = function(z,z.add) {
        integrate(function(m) {
            (gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m)*ifelse(z==0,1-expit(alpha[1]),expit(alpha[1]))*dnorm(m,beta[1] + beta[2] + beta[3]*z.add + beta[4]*z.add,eps[[2]][2])*ifelse(z.add==0,1-expit(alpha[1]+alpha[2]),expit(alpha[1]+alpha[2]))
        }, lower = m.range[1], upper = m.range[2])$value
    }
    #YS3.dprime
    integrate.s3.dprime = function(z,z.prime) {
        integrate(function(m) {
            (gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m)*ifelse(z==0,1-expit(alpha[1]),expit(alpha[1]))*dnorm(m,beta[1] + beta[2] + beta[3]*z.prime + beta[4]*z.prime,eps[[2]][2])*ifelse(z.prime==0,1-expit(alpha[1]),expit(alpha[1]))
        }, lower = m.range[1], upper = m.range[2])$value
    }
    s1 = sum(sapply(z.val, function(z) integrate.s1(z)))
    s2 = sum(sapply(z.val, function(z) sum(sapply(z.val, function(zprime) integrate.s2(z,zprime)))))
    s3 = sum(sapply(z.val, function(z) integrate.s3(z)))
    s1.prime = sum(sapply(z.val, function(z) sum(sapply(z.val, function(z.add) integrate.s1.prime(z,z.add)))))
    s2.prime = sum(sapply(z.val, function(z) sum(sapply(z.val, function(z.add) integrate.s2.prime(z,z.add)))))
    s3.dprime = sum(sapply(z.val, function(z) sum(sapply(z.val, function(z.prime) integrate.s3.dprime(z,z.prime)))))
return(c(s1,s2,s3,s1.prime,s2.prime,s3.dprime))
}

sigma = function(rho, alpha, beta, gamma, eps) {
    z.val = c(0,1)
    m.val = c(0,1)
    mu0 = expit(alpha[1])
    mu1 = expit(alpha[1] + alpha[2])
    Y_a0.z.m = function(z,m) {gamma[1] + gamma[3]*z + gamma[4]*m + gamma[7]*z*m}
    pz_a1 = function(z) {z*expit(alpha[1] + alpha[2])+(1-z)*(1-expit(alpha[1] + alpha[2]))}
    pz_a0 = function(z) {z*expit(alpha[1])+(1-z)*(1-expit(alpha[1]))}
    pm_a1.z = function(z,m) {m*expit(beta[1]+beta[2]+beta[3]*z+beta[4]*z)+(1-m)*(1-expit(beta[1]+beta[2]+beta[3]*z+beta[4]*z))}
    pm_a1.zprime = function(zprime,m) {m*expit(beta[1]+beta[2]+beta[3]*zprime+beta[4]*zprime)+(1-m)*(1-expit(beta[1]+beta[2]+beta[3]*zprime+beta[4]*zprime))}
    pz_z = function(z,zprime) {
        given_z1 = (rho * mu0^(1/2) * (1-mu0)^(1/2) * mu1^(1/2)*(1-mu1)^(1/2) + mu0*mu1)/mu0
        given_z0 = (mu1 - rho * mu0^(1/2) * (1-mu0)^(1/2) * mu1^(1/2)*(1-mu1)^(1/2) - mu0*mu1)/(1 - mu0)
        return(given_z1*z*zprime + (1-given_z1)*z*(1-zprime) + given_z0*(1-z)*zprime + (1-given_z0)*(1-z)*(1-zprime))
    }

    #YS1
    sigma.s1 = function(z,m) {Y_a0.z.m(z,m) * pm_a1.z(z,m) * pz_a1(z)}
    #YS2
    sigma.s2 = function(z,zprime,m) {Y_a0.z.m(z,m) * pm_a1.zprime(zprime,m) * pz_z(z,zprime) * pz_a0(z)}
    #YS3
    sigma.s3 = function(z,m) {Y_a0.z.m(z,m) * pm_a1.z(z,m) * pz_a0(z)}
    #YS1.prime
    sigma.s1.prime = function(z,zprime,m) {Y_a0.z.m(z,m) * pz_a1(z) * pm_a1.z(zprime,m) * pz_a1(zprime)}
    #YS2.prime
    sigma.s2.prime = function(z,zprime,m) {Y_a0.z.m(z,m) * pz_a0(z) * pm_a1.z(zprime,m) * pz_a1(zprime)}
    #YS3.dprime
    sigma.s3.dprime = function(z,zprime,m) {Y_a0.z.m(z,m) * pz_a0(z) * pm_a1.z(zprime,m) * pz_a0(zprime)}

    s1 = sum(sapply(z.val, function(z) sum(sapply(m.val, function(m) sigma.s1(z,m)))))
    s2 = sum(sapply(z.val, function(z) sum(sapply(z.val, function(zprime) sum(sapply(m.val, function(m) sigma.s2(z,zprime,m)))))))
    s3 = sum(sapply(z.val, function(z) sum(sapply(m.val, function(m) sigma.s3(z,m)))))
    s1.prime = sum(sapply(z.val, function(z) sum(sapply(z.val, function(zprime) sum(sapply(m.val, function(m) sigma.s1.prime(z,zprime,m)))))))
    s2.prime = sum(sapply(z.val, function(z) sum(sapply(z.val, function(zprime) sum(sapply(m.val, function(m) sigma.s2.prime(z,zprime,m)))))))
    s3.dprime = sum(sapply(z.val, function(z) sum(sapply(z.val, function(zprime) sum(sapply(m.val, function(m) sigma.s3.dprime(z,zprime,m)))))))
return(c(s1,s2,s3,s1.prime,s2.prime,s3.dprime))
}


generation = function(n, A, rho, alpha, beta, gamma, eps) {
    id = seq(1:n)
    Z0 = rbinom(n, 1, expit(alpha[1]))
    T0 = rbinom(n, 1, expit(alpha[1]))
    mu0 = expit(alpha[1])
    mu1 = expit(alpha[1] + alpha[2])
    Z1 = sapply(id, function(i) ifelse(Z0[i]==1,rbinom(1,1,(rho * mu0^(1/2) * (1-mu0)^(1/2) * mu1^(1/2)*(1-mu1)^(1/2) + mu0*mu1)/mu0),
                                                rbinom(1,1,(mu1 - rho * mu0^(1/2) * (1-mu0)^(1/2) * mu1^(1/2)*(1-mu1)^(1/2) - mu0*mu1)/(1 - mu0))))
    T1 = rbinom(n, 1, expit(alpha[1] + alpha[2]))
    epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
    M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
    M0 = beta[1] + beta[3]*Z0 + epsM
    M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
    M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
    corellation = cor(Z0,Z1,method="pearson")
    epsY = rnorm(n=n, mean=eps[[3]][1], sd=eps[[3]][2])
    Ys1 = gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1 + epsY
    Ys2 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1 + epsY
    Ys3 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0 + epsY

    Ys1.prime = gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1 + epsY
    Ys2.prime = gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1 + epsY
    Ys2.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1 + epsY
    Ys3.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0 + epsY
    return(data.frame(Z0,Z1,M1,M0,corellation,Ys1,Ys2,Ys3,Ys1.prime,Ys2.prime,Ys2.dprime,Ys3.dprime))
}

generation.1 = function(n, A, rho, alpha, beta, gamma, eps) {
    id = seq(1:n)
    z <- function(rho, p, q) {
        z0 <- rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
        if (any(c(z0, 1-q-z0, 1-p-z0, z0+p+q-1) < 0)) {
            warning("(rho, p, q) = (", rho, ", ", p, ", ", q, ") is an impossibility!")
            z0 <- NA
        }
        z0
    }
    p = expit(alpha[1])
    q = expit(alpha[1] + alpha[2])
    z0 = z(rho, p, q)
    prob = c(`(0,0)`=z0, `(1,0)`=1-q-z0, `(0,1)`=1-p-z0, `(1,1)`=z0+p+q-1)
    pick = sample.int(4, n, replace=TRUE, prob=prob)
    Z1 = floor((pick-1)/2)
    Z0 = 1 - pick %% 2

    T0 = rbinom(n, 1, expit(alpha[1]))
    T1 = rbinom(n, 1, expit(alpha[1] + alpha[2]))
    epsM = rnorm(n=n, mean=eps[[2]][1], sd=eps[[2]][2])
    M1 = beta[1] + beta[2] + beta[3]*Z1 + beta[4]*Z1 + epsM
    M0 = beta[1] + beta[3]*Z0 + epsM
    M1Z0 = beta[1] + beta[2] + beta[3]*Z0 + beta[4]*Z0 + epsM
    M1T1 = beta[1] + beta[2] + beta[3]*T1 + beta[4]*T1 + epsM
    corellation = cor(Z0,Z1,method="pearson")
    epsY = rnorm(n=n, mean=eps[[3]][1], sd=eps[[3]][2])

    Ys1 = gamma[1] + gamma[3]*Z1 + gamma[4]*M1 + gamma[7]*Z1*M1 + epsY
    Ys2 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1 + gamma[7]*Z0*M1 + epsY
    Ys3 = gamma[1] + gamma[3]*Z0 + gamma[4]*M1Z0 + gamma[7]*Z0*M1Z0 + epsY

    Ys1.prime = gamma[1] + gamma[3]*Z1 + gamma[4]*M1T1 + gamma[7]*Z1*M1T1 + epsY
    Ys2.prime = gamma[1] + gamma[3]*Z0 + gamma[4]*M1T1 + gamma[7]*Z0*M1T1 + epsY
    Ys2.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1 + gamma[7]*T0*M1 + epsY
    Ys3.dprime = gamma[1] + gamma[3]*T0 + gamma[4]*M1Z0 + gamma[7]*T0*M1Z0 + epsY
    return(data.frame(Z0,Z1,M1,M0,corellation,Ys1,Ys2,Ys3,Ys1.prime,Ys2.prime,Ys2.dprime,Ys3.dprime))
}


gen = c()
for (rho in seq(-0.5,1,0.5)) {
    result = data.frame(matrix(NA, nrow = 100, ncol = 8))
    for (i in 1:100) {
        set.seed(i)
        rho = rho
        n = 1e5
        alpha = c(0.5, 0)
        beta = c(0.5, 0.5, 0.5, 0.25)
        gamma = c(0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25)
        eps = list(c(NaN, NaN), c(0, 1), c(0, 1))
        dat = generation(n, 0.5, rho, alpha, beta, gamma, eps)
        result[i,] = colMeans(dat[,5:12])
    }
    colnames(result) = c("Cor","Ys1", "Ys2", "Ys3", "Ys1.prime", "Ys2.prime", "Ys2.dprime", "Ys3.dprime")
    gen = do.call(rbind,list(gen,colMeans(result)))
}


gen.1 = c()
for (rho in seq(-0.5,1,0.5)) {
    result = data.frame(matrix(NA, nrow = 100, ncol = 8))
    for (i in 1:100) {
        set.seed(i)
        rho = rho
        n = 1e5
        alpha = c(0.5, 0)
        beta = c(0.5, 0.5, 0.5, 0.25)
        gamma = c(0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25)
        eps = list(c(NaN, NaN), c(0, 1), c(0, 1))
        dat = generation.1(n, 0.5, rho, alpha, beta, gamma, eps)
        result[i,] = colMeans(dat[,5:12])
    }
    colnames(result) = c("Cor","Ys1", "Ys2", "Ys3", "Ys1.prime", "Ys2.prime", "Ys2.dprime", "Ys3.dprime")
    gen.1 = do.call(rbind,list(gen.1,colMeans(result)))
}


int = c()
for (rho in seq(-0.5,1,0.5)) {
    rho = rho
    n = 1e5
    alpha = c(0.5, 0)
    beta = c(0.5, 0.5, 0.5, 0.25)
    gamma = c(0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25)
    vec.s = integration(rho, alpha, beta, gamma, eps)
    int = do.call(rbind,list(int,vec.s))
}
colnames(int) = c("Ys1", "Ys2", "Ys3", "Ys1.prime", "Ys2.prime", "Ys3.dprime")

