library(invgamma)
library(mvtnorm)
library(copula)
library(LaplacesDemon)
library(Metrics)
library(parallel)
library(pbapply)
# library(dirichlet)


# sample size and true coefficients
N <- 1000
coeffs <- c(2, -4, 6)
sig2 <- 5

rho_ze <- .7
rho_xe <- .3


# Simulate sample
testdata <- as.data.frame(matrix(NA, nrow = N, ncol = 0))

sigma_mat1 <- matrix(c(1, rho_ze, 0,
                       rho_ze, 1, rho_xe,
                       0, rho_xe, 1),
                     ncol = 3, nrow = 3, byrow = TRUE)

eps <- rmvnorm(n = N, mean = c(0, 0, 0), sigma = sigma_mat1, method = "chol")

testdata$const <- rep(1, N)
testdata$z <- (qlnorm(p = pnorm(eps[, 2])) - exp(.5)) / sqrt((exp(1) - 1)*exp(1))
testdata$x <- qnorm(pnorm(eps[, 3]))
testdata$e <- qnorm(pnorm(eps[, 1]), sd = sqrt(sig2))


# dependent variable
testdata$y <- cbind(testdata$const, testdata$z, testdata$x) %*% coeffs + testdata$e


### Metropolis-Hastings and Gibbs

# auxiliary functions
aux1 <- function(h, param, W1, data) {
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  W01 <- W1
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  
  M01 <- matrix(ncol = 1, nrow = 3,
                data = c(qnorm(data$z1)[h], qnorm(data$x1)[h], 
                         qnorm(pnorm(data$y - a - b1*data$z - b2*data$x, sd = sqrt(param[4])))[h]))
  
  P01 <- t(M01)%*%(solve(W01) - diag(3))%*%M01
  return(P01)
  
}

# (log) posterior with normal priors on regression coefficients (and IG hyperpriors),
# IG for all variances (residual variance, variance of prior for regression
# coefficients), inverse Wishart for Copula correlation matrix, categorial 
# distribution for regressors, Dirichlet priors for probability masses 
# formalising PDF of regressors' distribution
post2 <- function(param, W, data) {
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  s <- param[4]
  p1 <- param[5]
  p2 <- param[6]
  p3 <- param[7]
  
  sa <- param[(7 + 2*N) + 1]
  sb1 <- param[(7 + 2*N) + 2]
  sb2 <- param[(7 + 2*N) + 3]
  
  e <- data$y - a - b1*data$z - b2*data$x
  
  X_hat <- matrix(0, nrow = N, ncol = 3)
  X_hat[, 1] <- data$z1
  X_hat[, 2] <- data$x1
  X_hat[, 3] <- pnorm(q = e, mean = 0, sd = sqrt(s))
  
  Phi1 <- matrix(nrow = 3, ncol = 3, byrow = TRUE,
                 data = c(1, p1, p2,
                          p1, 1, p3,
                          p2, p3, 1))
  
  
  t1 <- sum(dCopula(copula = normalCopula(param = P2p(Phi1), dim = 3, 
                                           dispstr = "un"), u = X_hat, 
                     log = TRUE)) + 
           sum(dnorm(mean = 0, sd = sqrt(s), x = e, log = TRUE)) +
    sum(LaplacesDemon::dcat(x = c(table(data$z)), p = data$z1g, log = TRUE)) +
    sum(LaplacesDemon::dcat(x = c(table(data$x)), p = data$x1g, log = TRUE)) +
    dinvwishart(Sigma = W, nu = 3, S = diag(3), log = TRUE) + 
    invgamma::dinvgamma(x = s, shape = .001, rate = .001, log = TRUE) + 
    LaplacesDemon::ddirichlet(x = t(data$z1g), alpha = table(data$z) + rep(1, N), log = TRUE) +
    LaplacesDemon::ddirichlet(x = t(data$x1g), alpha = table(data$x) + rep(1, N), log = TRUE) +
    dnorm(x = a, mean = 0, sd = sqrt(sa), log = TRUE) + 
    dnorm(x = b1, mean = 0, sd = sqrt(sb1), log = TRUE) +
    dnorm(x = b2, mean = 0, sd = sqrt(sb2), log = TRUE) +
    invgamma::dinvgamma(x = sa, shape = .001, rate = .001, log = TRUE) +
    invgamma::dinvgamma(x = sb1, shape = .001, rate = .001, log = TRUE) +
    invgamma::dinvgamma(x = sb2, shape = .001, rate = .001, log = TRUE)
  return(t1)
  
}

# proposals
proposals1 <- function(param, W1, data) {
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]

  P <- param[4]*solve(t(cbind(data$const, qnorm(data$z1), qnorm(data$x1))) %*% cbind(data$const, qnorm(data$z1), qnorm(data$x1)))
  
  # Proposals for beta
  p1 <- rmvnorm(n = 1, mean = c(a, b1, b2), sigma = P)
  # p1 <- rmvt(n = 1, delta = c(a, b1, b2), sigma = P, df = 1)

  return(p1)
  
}
proposals2 <- function(param, W1, data) {
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  
  # Proposals for sigma2
  P01 <- sapply(X = 1:N, FUN = aux1, param = param, W1 = W1, data = data)
  
  p2 <- invgamma::rinvgamma(n = 1, shape = N/2, 
                            rate = abs(sum((data$y - a - b1*data$z - b2*data$x)^2)/2 - sum(P01)/2))
  
  return(p2)
  
}
d_proposals <- function(param, param1, W1, data) {
  
  data$z1g <- param[8:(N+7)]
  data <- data[order(data$z), ]
  data$z1 <- cumsum(data$z1g)*(length(data$z1g) / (length(data$z1g) + .01))
  
  data$x1g <- param[(N+8):(7 + 2*N)]
  data <- data[order(data$x), ]
  data$x1 <- cumsum(data$x1g)*(length(data$x1g) / (length(data$x1g) + .01))
  
  a <- param[1]
  b1 <- param[2]
  b2 <- param[3]

  a0 <- param1[1]
  b01 <- param1[2]
  b02 <- param1[3]

  P01 <- sapply(X = 1:N, FUN = aux1, param = param1, W1 = W1, data = data)
  
  p2 <- invgamma::dinvgamma(x = param[4], shape = N/2, log = TRUE,
                            rate = abs(sum((data$y - a - b1*data$z - b2*data$x)^2)/2 - sum(P01)/2))
  
  return(p2)
  
}

# MCMC algorithm
metropolis_Gibbs_MCMC1 <- function(startvalue, iterations, data) {
  
  dat1 <- data
  
  chain <- matrix(data = NA, nrow = iterations + 1, ncol = (7 + 2*N) + 3)
  
  chain[1, ] <- startvalue
  
  W01 <- diag(3)
  W1 <- diag(3)
  
  
  for (i in 1:iterations) {
    
    
    ################################# Proposals ################################
    
    proposal <- chain[i, ]
    proposal[1:3] <- proposals1(param = chain[i, ], W1 = W01, data = dat1)
    
    
    # acceptance probability
    probab <- min(exp(post2(param = proposal, W = W1, data = dat1) - 
                        post2(param = chain[i, ], W = W1, data = dat1)), 1)
    if (is.nan(probab)) {probab <- 0}
    
    # Metropolis step betas
    if (runif(1) < probab) {
      
      chain[i + 1, ] <- proposal
      
      dat1$z1g <- proposal[8:(N+7)]
      dat1 <- dat1[order(dat1$z), ]
      dat1$z1 <- cumsum(dat1$z1g)*(length(dat1$z1g) / (length(dat1$z1g) + .01))
      
      dat1$x1g <- proposal[(N+8):(7 + 2*N)]
      dat1 <- dat1[order(dat1$x), ]
      dat1$x1 <- cumsum(dat1$x1g)*(length(dat1$x1g) / (length(dat1$x1g) + .01))
      
    } else {
      
      chain[i + 1, ] <- chain[i, ]
      
      dat1$z1g <- proposal[8:(N+7)]
      dat1 <- dat1[order(dat1$z), ]
      dat1$z1 <- cumsum(dat1$z1g)*(length(dat1$z1g) / (length(dat1$z1g) + .01))
      
      dat1$x1g <- proposal[(N+8):(7 + 2*N)]
      dat1 <- dat1[order(dat1$x), ]
      dat1$x1 <- cumsum(dat1$x1g)*(length(dat1$x1g) / (length(dat1$x1g) + .01))
      
    }
    
    
    # Metropolis step sigma2
    proposal <- chain[i + 1, ]
    proposal[4] <- proposals2(param = chain[i + 1, ], W1 = W01, data = dat1)
    
    probab <- min(exp(post2(param = proposal, W = W1, data = dat1) - 
                        post2(param = chain[i + 1, ], W = W1, data = dat1) + 
                        d_proposals(param = chain[i + 1, ], param1 = proposal, W1 = W01, data = dat1) - 
                        d_proposals(param = proposal, param1 = chain[i + 1, ], W1 = W01, data = dat1)), 1)
    if (is.nan(probab)) {probab <- 0}
    
    if (runif(1) < probab) {
      
      chain[i + 1, ] <- proposal
      
    } 
    
    
    ############################ Gibbs step Wishart ############################
    
    X01 <- matrix(ncol = 3, nrow = N, 
                  data = c(qnorm(dat1$z1),
                           qnorm(dat1$x1),
                           qnorm(pnorm(dat1$y - chain[i+1, 1] - 
                                         chain[i+1, 2]*dat1$z - chain[i+1, 3]*dat1$x, 
                                       sd = sqrt(chain[i+1, 4])))))
    
    W1 <- rinvwishart(nu = N + 3, S = diag(3) + t(X01)%*%X01)
    
    while (1 > 0) {
      
      W1 <- rinvwishart(nu = N + 3, S = diag(3) + t(X01)%*%X01) 
      W1[2, 3] <- 0
      W1[3, 2] <- 0
      
      if (min(eigen(W1)$values) > 0) { break }
      
    }
    
    W01 <- solve(sqrt(diag(W1))*diag(3))%*%W1%*%solve(sqrt(diag(W1))*diag(3))
    
    chain[i+1, 5:7] <- P2p(W01)
    # chain[i+1, 5] <- cor(qnorm(pobs(X[, 2])), qnorm(pobs(X[, 3])))
    # chain[i+1, 7] <- 0
    # W1[2, 3] <- 0
    # W1[3, 2] <- 0
    
    
    ########################### Gibbs step Dirichlet ###########################
    
    eps_n <- mvtnorm::rmvnorm(n = N, mean = rep(0, 3), sigma = W01, method = "eigen")
    
    pt_z <- qgamma(p = pnorm(eps_n[, 1]), shape = table(dat1$z) + rep(1, length(dat1$z)), rate = 1)
    chain[i+1, c(8:(N + 7))] <- pt_z/sum(pt_z)
    
    pt_x <- qgamma(p = pnorm(eps_n[, 2]), shape = table(dat1$x) + rep(1, length(dat1$x)), rate = 1)
    chain[i+1, c((N + 8):(2*N + 7))] <- pt_x/sum(pt_x)
    
    
    dat1$z1g <- NA
    dat1$z1 <- NA
    dat1$x1g <- NA
    dat1$x1 <- NA
    
    ### Hyperpriors
    
    chain[i + 1, (7 + 2*N) + 1] <- invgamma::rinvgamma(n = 1, shape = .001 + .5,
                                                       rate = .5*(chain[i + 1, 1]^2 + 2*.001))
    chain[i + 1, (7 + 2*N) + 2] <- invgamma::rinvgamma(n = 1, shape = .001 + .5,
                                                       rate = .5*(chain[i + 1, 2]^2 + 2*.001))
    chain[i + 1, (7 + 2*N) + 3] <- invgamma::rinvgamma(n = 1, shape = .001 + .5,
                                                       rate = .5*(chain[i + 1, 3]^2 + 2*.001))
    
  }
  
  return(chain)
  
}


dataset <- testdata

# OLS for starting values
mod1 <- lm(y ~ z + x, dataset)

startvalue = c(mod1$coefficients, var(mod1$residuals), rep(0, 3), 
               c(LaplacesDemon::rdirichlet(n = 1, rep(1, length(dataset$z)))), 
               c(LaplacesDemon::rdirichlet(n = 1, rep(1, length(dataset$x)))),
               rep(1000, 3))

iters <- 102000

chain001 <- metropolis_Gibbs_MCMC1(startvalue, iterations= iters, 
                                   data = dataset)

# burnin: disregard the first 2000 iterations
burnin <- 2000
chain01 <- chain001[(burnin + 1):iters, ]


# thinning: keep only every 100th observation
thin <- 100
chain1 <- chain01[seq(1, nrow(chain01), thin), ]


### Posterior diagnostics

# Intercept
summary(chain1[, 1])
hist(chain1[, 1])
plot.ts(chain1[, 1])
acf(chain1[, 1])


# Endogenous regressor
summary(chain1[, 2])
hist(chain1[, 2])
plot.ts(chain1[, 2])
acf(chain1[, 2])


# Exogenous regressor
summary(chain1[, 3])
hist(chain1[, 3])
plot.ts(chain1[, 3])
acf(chain1[, 3])


# Error variance
summary(chain1[, 4])
hist(chain1[, 4])
plot.ts(chain1[, 4])
acf(chain1[, 4])


# Correlation between endogenous and exogenous regressor
summary(chain1[, 5])
hist(chain1[, 5])
plot.ts(chain1[, 5])
acf(chain1[, 5])


# Correlation between endogenous regressor and error
summary(chain1[, 6])
hist(chain1[, 6])
plot.ts(chain1[, 6])
acf(chain1[, 6])


# Correlation between exogenous regressor and error (forced to be zero)
summary(chain1[, 7])


# Probability masses formalising PDF of endogenous regressor's distribution
summary(chain1[, 8]) # first observation of z
hist(chain1[, 8])
plot.ts(chain1[, 8])
acf(chain1[, 8])

summary(chain1[, 9]) # second observation of z
hist(chain1[, 9])
plot.ts(chain1[, 9])
acf(chain1[, 9])

# ... and so on


# Take all of them and generate the Bayesian CDF (posterior mean and 95%
# posterior credible interval)
ord <- order(dataset$z)
z_sorted <- dataset$z[ord]

lam_cumsums <- apply(chain1[, 8:(N+7)], 1, function(lam) {
  cumsum(lam[ord])
})
lam_cumsums <- t(lam_cumsums)

mean_vals <- apply(lam_cumsums, 2, mean)
low_vals  <- apply(lam_cumsums, 2, quantile, 0.025)
high_vals <- apply(lam_cumsums, 2, quantile, 0.975)

plot(z_sorted, mean_vals, type = "l", lwd = 2,
     ylab = "Posterior CDF", xlab = "z")
lines(z_sorted, low_vals, lty = 2, col = "grey40")
lines(z_sorted, high_vals, lty = 2, col = "grey40")


# Same for exogenous regressor
ord <- order(dataset$x)
z_sorted <- dataset$x[ord]

lam_cumsums <- apply(chain1[, (N + 8):(2*N + 7)], 1, function(lam) {
  cumsum(lam[ord])
})
lam_cumsums <- t(lam_cumsums)

mean_vals <- apply(lam_cumsums, 2, mean)
low_vals  <- apply(lam_cumsums, 2, quantile, 0.025)
high_vals <- apply(lam_cumsums, 2, quantile, 0.975)

plot(z_sorted, mean_vals, type = "l", lwd = 2,
     ylab = "Posterior CDF", xlab = "x")
lines(z_sorted, low_vals, lty = 2, col = "grey40")
lines(z_sorted, high_vals, lty = 2, col = "grey40")


# Variance of normal prior for intercept
summary(chain1[, (7 + 2*N) + 1])
hist(chain1[, (7 + 2*N) + 1])
plot.ts(chain1[, (7 + 2*N) + 1])
acf(chain1[, (7 + 2*N) + 1])


# Variance of normal prior for delta
summary(chain1[, (7 + 2*N) + 2])
hist(chain1[, (7 + 2*N) + 2])
plot.ts(chain1[, (7 + 2*N) + 2])
acf(chain1[, (7 + 2*N) + 2])


# Variance of normal prior for beta
summary(chain1[, (7 + 2*N) + 3])
hist(chain1[, (7 + 2*N) + 3])
plot.ts(chain1[, (7 + 2*N) + 3])
acf(chain1[, (7 + 2*N) + 3])




