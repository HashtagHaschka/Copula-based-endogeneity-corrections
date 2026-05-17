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

eps <- mvtnorm::rmvnorm(n = N, mean = c(0, 0, 0), sigma = sigma_mat1, method = "chol")

testdata$const <- rep(1, N)
testdata$z <- (qlnorm(p = pnorm(eps[, 2])) - exp(.5)) / sqrt((exp(1) - 1)*exp(1))
testdata$x <- qnorm(pnorm(eps[, 3]))
testdata$e <- qnorm(pnorm(eps[, 1]), sd = sqrt(sig2))


# dependent variable
testdata$y <- cbind(testdata$const, testdata$z, testdata$x) %*% coeffs + testdata$e


### Auxiliary functions
build_margin <- function(v) {
  uv  <- sort(unique(v))
  grp <- match(v, uv)
  m   <- length(uv)
  cnt <- tabulate(grp, nbins = m)
  list(uv = uv, grp = grp, m = m, cnt = cnt)
}
xi_from_lambda <- function(lambda, mg) {
  Fmid <- (cumsum(lambda) - lambda / 2) * (mg$m / (mg$m + .01))
  qnorm(Fmid[mg$grp])
}
draw_lambda <- function(u_channel, mg) {
  g_obs   <- qgamma(u_channel, shape = 1, rate = 1)
  G_data  <- as.numeric(tapply(g_obs,
                               factor(mg$grp, levels = seq_len(mg$m)),
                               sum))
  G_data[is.na(G_data)] <- 0
  G_prior <- rgamma(mg$m, shape = 1, rate = 1)
  G <- G_prior + G_data
  G / sum(G)
}
post2 <- function(param, lamZ, lamX, mgz, mgx, data) {
  
  a  <- param[1]
  b1 <- param[2]
  b2 <- param[3]
  s  <- param[4]
  p1 <- param[5]
  p2 <- param[6]
  p3 <- param[7]
  
  sa  <- param[(7 + 2*N) + 1]
  sb1 <- param[(7 + 2*N) + 2]
  sb2 <- param[(7 + 2*N) + 3]
  
  xi_z <- xi_from_lambda(lamZ, mgz)
  xi_x <- xi_from_lambda(lamX, mgx)
  
  e <- data$y - a - b1*data$z - b2*data$x
  
  X_hat <- matrix(0, nrow = N, ncol = 3)
  X_hat[, 1] <- pnorm(xi_z)
  X_hat[, 2] <- pnorm(xi_x)
  X_hat[, 3] <- pnorm(q = e, mean = 0, sd = sqrt(s))
  X_hat[] <- pmin(pmax(X_hat, 1e-12), 1 - 1e-12)
  
  Phi1 <- matrix(nrow = 3, ncol = 3, byrow = TRUE,
                 data = c(1, p1, p2,
                          p1, 1, p3,
                          p2, p3, 1))
  
  t1 <- sum(dCopula(copula = normalCopula(param = P2p(Phi1), dim = 3,
                                          dispstr = "un"), u = X_hat,
                    log = TRUE)) +
    sum(dnorm(mean = 0, sd = sqrt(s), x = e, log = TRUE)) +
    invgamma::dinvgamma(x = s, shape = .001, rate = .001, log = TRUE) +
    dnorm(x = a, mean = 0, sd = sqrt(sa), log = TRUE) +
    dnorm(x = b1, mean = 0, sd = sqrt(sb1), log = TRUE) +
    dnorm(x = b2, mean = 0, sd = sqrt(sb2), log = TRUE)
  return(t1)
  
}
compute_logsig2_derivs <- function(a, b1, b2, sig2, p1, p2, p3,
                                   xi_z, xi_x, dat1) {
  e    <- dat1$y - a - b1 * dat1$z - b2 * dat1$x
  xi_e <- e / sqrt(sig2)
  
  Phi <- matrix(c(1, p1, p2,
                  p1, 1, p3,
                  p2, p3, 1), nrow = 3, byrow = TRUE)
  A <- solve(Phi) - diag(3)
  
  A31 <- A[3, 1]
  A32 <- A[3, 2]
  A33 <- A[3, 3]
  Xinv33 <- A33 + 1   # (Xi^{-1})_{33}
  
  # Score (first derivative)
  Score <- sum(0.5 * (A31 * xi_z * xi_e +
                        A32 * xi_x * xi_e +
                        Xinv33 * xi_e^2) - 0.5) -
    0.001 + 0.001 / sig2
  
  # Second derivative (negative near the mode; P = -1/f2 below)
  f2 <- sum(-0.25 * (A31 * xi_z * xi_e + A32 * xi_x * xi_e) -
              0.5 * Xinv33 * xi_e^2) -
    0.001 / sig2
  
  list(Score = Score, f2 = f2)
}

### Sampler (Metropolis-Hastings and Gibbs)
metropolis_Gibbs_MCMC1 <- function(startvalue, iterations, data) {
  
  dat1 <- data                                                                  # FIXED observation order for the whole run
  
  chain <- matrix(data = NA, nrow = iterations + 1, ncol = (7 + 2*N) + 3)
  chain[1, ] <- startvalue
  
  W01 <- diag(3)
  n   <- N
  K   <- 1   # dim of xi_z
  L   <- 1   # dim of xi_x
  Omega_cur <- diag(K)
  
  # Per-unique-value margin structures (fixed for the run)
  mgz <- build_margin(dat1$z)
  mgx <- build_margin(dat1$x)
  
  # current per-cell mass vectors. chain keeps its original column layout;
  # only the first mgz$m / mgx$m mass columns are used (rest stay 0).
  lamZ <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, mgz$m)))
  lamX <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, mgx$m)))
  chain[1, 8:(N + 7)]         <- 0
  chain[1, 8:(7 + mgz$m)]     <- lamZ
  chain[1, (N + 8):(2*N + 7)] <- 0
  chain[1, (N + 8):(N + 7 + mgx$m)] <- lamX
  
  for (i in 1:iterations) {
    
    if (i %% 100 == 0) {
      cat("Iteration:", i,
          "| beta =", round(chain[i,1:3], 3),
          "| sigma2 =", round(chain[i,4], 3), "\n")
    }
    
    # current latent margins from current per-cell masses (step PIT)
    xi_z <- xi_from_lambda(lamZ, mgz)
    xi_x <- xi_from_lambda(lamX, mgx)
    
    ############################################################################
    # MH step for regression coefficients - IWLS / Newton proposal with
    # Hastings correction. Sigma_prop is beta-independent in this sub-step
    # (s and the rho's are held fixed), so reusing it in the reverse density
    # is exact, not an approximation.
    ############################################################################
    
    cur_a  <- chain[i, 1]
    cur_b1 <- chain[i, 2]
    cur_b2 <- chain[i, 3]
    cur_s2 <- chain[i, 4]
    cur_p1 <- chain[i, 5]
    cur_p2 <- chain[i, 6]
    cur_p3 <- chain[i, 7]
    
    sigma <- sqrt(cur_s2)
    
    eta_cur <- cur_a + cur_b1 * dat1$z + cur_b2 * dat1$x
    e_cur   <- dat1$y - eta_cur
    
    xi_e <- e_cur / sigma
    
    Phi <- matrix(c(1, cur_p1, cur_p2,
                    cur_p1, 1, cur_p3,
                    cur_p2, cur_p3, 1), nrow = 3, byrow = TRUE)
    invPhi <- solve(Phi)
    A <- invPhi - diag(3)               # Xi^{-1} - I
    invPhi33 <- invPhi[3, 3]            # (Xi^{-1})_{33}
    
    # Score vector nu (derivative of log-lik w.r.t. eta_i)
    t3 <- A[3,1]*xi_z + A[3,2]*xi_x + A[3,3]*xi_e
    nu <- (1/sigma) * t3 + e_cur / cur_s2   # length N
    
    X <- cbind(dat1$const, dat1$z, dat1$x)
    
    # Working weight: exact observed information per obs = (Xi^{-1})_{33}/sigma^2
    w_vec <- invPhi33 / cur_s2
    
    XtWX <- w_vec * crossprod(X)
    Sigma_prop <- solve(XtWX)
    Sigma_prop <- (Sigma_prop + t(Sigma_prop)) / 2   # ensure symmetry
    
    mu_cur <- c(cur_a, cur_b1, cur_b2) + Sigma_prop %*% (t(X) %*% nu)
    mu_cur <- as.vector(mu_cur)
    
    prop_coef <- mvtnorm::rmvnorm(1, mean = mu_cur, sigma = Sigma_prop)
    prop_coef <- as.vector(prop_coef)
    
    proposal <- chain[i, ]
    proposal[1:3] <- prop_coef
    
    log_post_cur  <- post2(param = chain[i, ], lamZ, lamX, mgz, mgx, dat1)
    log_post_prop <- post2(param = proposal,   lamZ, lamX, mgz, mgx, dat1)
    
    log_q_fwd <- mvtnorm::dmvnorm(prop_coef, mean = mu_cur,
                                  sigma = Sigma_prop, log = TRUE)
    
    eta_prop  <- prop_coef[1] + prop_coef[2] * dat1$z + prop_coef[3] * dat1$x
    e_prop    <- dat1$y - eta_prop
    xi_e_prop <- e_prop / sigma
    t3_prop   <- A[3,1]*xi_z + A[3,2]*xi_x + A[3,3]*xi_e_prop
    nu_prop   <- (1/sigma) * t3_prop + e_prop / cur_s2
    mu_new    <- prop_coef + Sigma_prop %*% (t(X) %*% nu_prop)
    mu_new    <- as.vector(mu_new)
    log_q_rev <- mvtnorm::dmvnorm(c(cur_a, cur_b1, cur_b2), mean = mu_new,
                                  sigma = Sigma_prop, log = TRUE)
    
    q_ratio   <- log_q_rev - log_q_fwd
    log_alpha <- log_post_prop - log_post_cur + q_ratio
    if (is.nan(log_alpha)) log_alpha <- -Inf
    
    if (log(runif(1)) < log_alpha) {
      chain[i+1, ] <- proposal
    } else {
      chain[i+1, ] <- chain[i, ]
    }
    
    ############################################################################
    # MH step for log sigma^2 using a 2nd-order Taylor (Laplace) proposal
    ############################################################################
    
    cur_a  <- chain[i+1, 1]
    cur_b1 <- chain[i+1, 2]
    cur_b2 <- chain[i+1, 3]
    cur_s2 <- chain[i+1, 4]
    cur_p1 <- chain[i+1, 5]
    cur_p2 <- chain[i+1, 6]
    cur_p3 <- chain[i+1, 7]
    
    tau_cur <- log(cur_s2)
    
    deriv_cur <- compute_logsig2_derivs(cur_a, cur_b1, cur_b2, cur_s2,
                                        cur_p1, cur_p2, cur_p3,
                                        xi_z, xi_x, dat1)
    
    if (deriv_cur$f2 < 0) {
      P_cur <- -1 / deriv_cur$f2
    } else {
      P_cur <- 1.0   # fallback (should not happen near the mode)
    }
    mu_cur <- P_cur * deriv_cur$Score + tau_cur
    
    tau_new <- rnorm(1, mean = mu_cur, sd = sqrt(P_cur))
    s2_new  <- exp(tau_new)
    
    prop_param <- chain[i+1, ]
    prop_param[4] <- s2_new
    
    deriv_new <- compute_logsig2_derivs(cur_a, cur_b1, cur_b2, s2_new,
                                        cur_p1, cur_p2, cur_p3,
                                        xi_z, xi_x, dat1)
    
    if (deriv_new$f2 < 0) {
      P_new <- -1 / deriv_new$f2
    } else {
      P_new <- 1.0
    }
    mu_new <- P_new * deriv_new$Score + tau_new
    
    q_ratio <- dnorm(tau_cur, mean = mu_new, sd = sqrt(P_new), log = TRUE) -
      dnorm(tau_new, mean = mu_cur, sd = sqrt(P_cur), log = TRUE)
    
    log_post_cur <- post2(param = chain[i+1, ], lamZ, lamX, mgz, mgx, dat1)
    log_post_new <- post2(param = prop_param,  lamZ, lamX, mgz, mgx, dat1)
    
    log_alpha <- log_post_new - log_post_cur + q_ratio
    if (is.nan(log_alpha)) log_alpha <- -Inf
    
    if (log(runif(1)) < log_alpha) {
      chain[i+1, ] <- prop_param   # accept
    }
    # else chain[i+1, ] stays at current sigma^2
    
    ############################ Gibbs step Wishart ############################
    # Block-recursive sampler of W. xi_z / xi_x are the
    # step-PIT margins from the CURRENT masses; xi_eps from the updated
    # beta, sigma^2 (Algorithm 1 order).
    
    xi_eps <- qnorm(pnorm(dat1$y - chain[i+1, 1] -
                            chain[i+1, 2]*dat1$z -
                            chain[i+1, 3]*dat1$x,
                          sd = sqrt(chain[i+1, 4])))
    
    xi_z_m <- matrix(xi_z, n, 1)
    xi_x_m <- matrix(xi_x, n, 1)
    
    # Sufficient statistics
    S_xx   <- t(xi_x_m) %*% xi_x_m
    S_ee   <- sum(xi_eps^2)
    X_tilde <- cbind(xi_x_m, xi_eps)
    S_tilde    <- t(X_tilde) %*% X_tilde
    S_z_xtilde <- t(xi_z_m)  %*% X_tilde
    
    nu1   <- L + 2
    Psi1  <- diag(L)
    a_e   <- .001
    b_e   <- .001
    M0    <- matrix(0, nrow = K, ncol = L + 1)
    V0    <- diag(L + 1)
    V0_inv <- solve(V0)
    nu2   <- K + 2
    Psi2  <- diag(K)
    
    Sigma_xx_cur <- LaplacesDemon::rinvwishart(nu = nu1 + n, S = Psi1 + S_xx)
    sigma_e2_cur <- invgamma::rinvgamma(n = 1, shape = a_e + n / 2, rate = b_e + S_ee / 2)
    
    V_bar <- solve(V0_inv + S_tilde)
    M_bar <- (M0 %*% V0_inv + S_z_xtilde) %*% V_bar
    E_mat    <- matrix(rnorm(K * (L + 1)), nrow = K, ncol = L + 1)
    chol_Om  <- t(chol(Omega_cur))
    chol_Vb  <- chol(V_bar)
    C_cur    <- M_bar + chol_Om %*% E_mat %*% chol_Vb
    
    Resid <- xi_z_m - X_tilde %*% t(C_cur)
    S_eps <- t(Resid) %*% Resid
    Omega_cur <- LaplacesDemon::rinvwishart(nu = nu2 + n, S = Psi2 + S_eps)
    
    B_x    <- C_cur[, 1:L,   drop = FALSE]
    beta_e <- C_cur[, L + 1, drop = FALSE]
    
    W_11 <- Omega_cur + B_x %*% Sigma_xx_cur %*% t(B_x) + sigma_e2_cur * beta_e %*% t(beta_e)
    W_12 <- B_x %*% Sigma_xx_cur
    W_13 <- sigma_e2_cur * beta_e
    W_22 <- Sigma_xx_cur
    W_23 <- matrix(0, nrow = L, ncol = 1)
    W_33 <- matrix(sigma_e2_cur, 1, 1)
    
    W1 <- rbind(
      cbind(W_11, W_12, W_13),
      cbind(t(W_12), W_22, W_23),
      cbind(t(W_13), t(W_23), W_33)
    )
    
    sds    <- sqrt(diag(W1))
    D_inv  <- diag(1 / sds)
    W01 <- D_inv %*% W1 %*% D_inv
    
    chain[i+1, 5:7] <- P2p(W01)
    
    ########################### Gibbs step Dirichlet ###########################
    # one Xi-correlated triple PER OBSERVATION (row of eps_n);
    # the per-observation Gamma(1,1) variates are aggregated WITHIN value-cells
    # (plus a Gamma(1,1) prior pseudo-count per cell) to give the per-cell
    # Dirichlet lambda_w ~ Dir(1 + n_j).
    
    eps_n <- mvtnorm::rmvnorm(n = N, mean = rep(0, 3), sigma = W01, method = "eigen")
    
    lamZ <- draw_lambda(pnorm(eps_n[, 1]), mgz)
    lamX <- draw_lambda(pnorm(eps_n[, 2]), mgx)
    
    # keep chain mass columns populated for inspection (first m entries are
    # the per-cell masses; remaining entries left as 0)
    chain[i+1, 8:(N + 7)]               <- 0
    chain[i+1, 8:(7 + mgz$m)]           <- lamZ
    chain[i+1, (N + 8):(2*N + 7)]       <- 0
    chain[i+1, (N + 8):(N + 7 + mgx$m)] <- lamX
    
    ### Hyperpriors
    chain[i+1, (7+2*N)+1] <- invgamma::rinvgamma(1, shape = 0.501,
                                                 rate = (chain[i+1,1]^2 + 0.002)/2)
    chain[i+1, (7+2*N)+2] <- invgamma::rinvgamma(1, shape = 0.501,
                                                 rate = (chain[i+1,2]^2 + 0.002)/2)
    chain[i+1, (7+2*N)+3] <- invgamma::rinvgamma(1, shape = 0.501,
                                                 rate = (chain[i+1,3]^2 + 0.002)/2)
    
  }
  
  attr(chain, "mgz") <- mgz
  attr(chain, "mgx") <- mgx
  return(chain)
}


dataset <- testdata

# OLS for starting values
mod1 <- lm(y ~ z + x, dataset)

startvalue = c(mod1$coefficients,
               var(mod1$residuals), rep(0, 3),
               rep(0, N),
               rep(0, N),
               rep(1000, 3))

iters <- 12000

chain001 <- metropolis_Gibbs_MCMC1(startvalue, iterations= iters,
                                   data = dataset)
mgz <- attr(chain001, "mgz")
mgx <- attr(chain001, "mgx")

# burnin: disregard the first 200 iterations
burnin <- 2000
chain01 <- chain001[(burnin + 1):iters, ]


# thinning: keep only every 10th observation
thin <- 10
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
# (per-unique-value: first mgz$m columns of the z-mass block)
summary(chain1[, 8]) # mass for z's first unique value
hist(chain1[, 8])
plot.ts(chain1[, 8])
acf(chain1[, 8])

summary(chain1[, 9]) # mass for z's second unique value
hist(chain1[, 9])
plot.ts(chain1[, 9])
acf(chain1[, 9])

# ... and so on


# Bayesian CDF of z (posterior mean and 95% credible interval), at the
# UNIQUE values: cumsum of the per-cell masses.
lamZ_draws <- chain1[, 8:(7 + mgz$m), drop = FALSE]
cdf_z <- t(apply(lamZ_draws, 1, cumsum))
mean_vals <- apply(cdf_z, 2, mean)
low_vals  <- apply(cdf_z, 2, quantile, 0.025)
high_vals <- apply(cdf_z, 2, quantile, 0.975)

plot(mgz$uv, mean_vals, type = "s", lwd = 2,
     ylab = "Posterior CDF", xlab = "z")
lines(mgz$uv, low_vals, lty = 2, type = "s", col = "grey40")
lines(mgz$uv, high_vals, lty = 2, type = "s", col = "grey40")


# Same for exogenous regressor x
lamX_draws <- chain1[, (N + 8):(N + 7 + mgx$m), drop = FALSE]
cdf_x <- t(apply(lamX_draws, 1, cumsum))
mean_vals <- apply(cdf_x, 2, mean)
low_vals  <- apply(cdf_x, 2, quantile, 0.025)
high_vals <- apply(cdf_x, 2, quantile, 0.975)

plot(mgx$uv, mean_vals, type = "s", lwd = 2,
     ylab = "Posterior CDF", xlab = "x")
lines(mgx$uv, low_vals, lty = 2, type = "s", col = "grey40")
lines(mgx$uv, high_vals, lty = 2, type = "s", col = "grey40")


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
