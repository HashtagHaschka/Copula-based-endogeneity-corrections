# Load required packages
pacman::p_load(
  MASS,
  copula,
  mvtnorm,
  Metrics,
  matrixcalc,
  Formula,
  dplyr,
  extRC,
  nlme,
  pbapply,
  tsoutliers,
  data.table,
  nortest,
  ks,
  resample,
  parallel
)


############################## Proposed estimator ##############################


# auxiliary functions
FOD <- function(x) {
  
  if (length(x) == 1) {
    
    foo1 <- NA
    
  } else {
    
    D <- extRC::dfm(length(x))
    A <- chol(solve(D%*%t(D)))%*%D
    foo1 <- c(A%*%as.numeric(x), NA)
    
  }
  
  return(foo1)
  
}
FOD1 <- Vectorize(FOD)

likelihood1 <- function(values, M1, P1, Y1, M_star, M_star1) {
  
  betas1 <- values[1:dim(M1)[2]]
  rhos1 <- tanh(values[(dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])])
  s1 <- exp(values[length(values)])
  
  e1 <- Y1 - t(betas1%*%t(M1))
  d1 <- sum(dnorm(x = e1[, 1], mean = 0, sd = sqrt(s1), log = TRUE))
  
  e1_star1 <- pnorm(q = e1[, 1], mean = 0, sd = sqrt(s1))
  if (max(e1_star1) == 1) { e1_star1 <- pobs(e1[, 1]) }
  
  rhos2 <- copula::P2p(cor(M_star))
  rhos <- c(rhos1, rep(0, dim(M1)[2] - dim(P1)[2]), rhos2)
  
  U1 <- as.matrix(cbind(e1_star1, M_star1))
  d2 <- sum(copula::dCopula(copula = copula::normalCopula(param = rhos, dim = dim(U1)[2], 
                                                          dispstr = "un"), u = U1, 
                            log = TRUE))
  
  d3 <- (d1 + d2)*(-1)
  
  if(is.infinite(d3)) {d3 <- 1e12}
  
  return(d3)
  
}
likelihood2 <- function(values, M1, P1, Y1, M_star, M_star1) {
  
  betas1 <- values[1:(dim(M1)[2] + 1)]
  rhos1 <- tanh(values[(dim(M1)[2] + 2):(dim(M1)[2] + dim(P1)[2] + 1)])
  s1 <- exp(values[length(values)])
  
  e1 <- Y1 - t(betas1%*%t(cbind(rep(1, dim(M1)[1]), M1)))
  d1 <- sum(dnorm(x = e1[, 1], mean = 0, sd = sqrt(s1), log = TRUE))
  
  e1_star1 <- pnorm(q = e1[, 1], mean = 0, sd = sqrt(s1))
  if (max(e1_star1) == 1) { e1_star1 <- pobs(e1[, 1]) }
  
  rhos2 <- copula::P2p(cor(M_star))
  rhos <- c(rhos1, rep(0, dim(M1)[2] - dim(P1)[2]), rhos2)
  
  U1 <- as.matrix(cbind(e1_star1, M_star1))
  d2 <- sum(copula::dCopula(copula = copula::normalCopula(param = rhos, dim = dim(U1)[2], 
                                                          dispstr = "un"), u = U1, 
                            log = TRUE))
  
  d3 <- (d1 + d2)*(-1)
  
  if(is.infinite(d3)) {d3 <- 1e12}
  
  return(d3)
  
}
likelihood_b1 <- function(values, M1, P1, Y1, M_star, M_star1) {
  
  betas1 <- values[1:dim(M1)[2]]
  rhos1 <- values[(dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])]
  s1 <- values[length(values)]
  
  if ((sum(rhos1 > -1) == length(rhos1)) && ((sum(rhos1 < 1) == length(rhos1))) && (sum(s1 > 0) == length(s1))) {
    
    e1 <- Y1 - t(betas1%*%t(M1))
    d1 <- sum(dnorm(x = e1[, 1], mean = 0, sd = sqrt(s1), log = TRUE))
    
    e1_star1 <- pnorm(q = e1[, 1], mean = 0, sd = sqrt(s1))
    if (max(e1_star1) == 1) { e1_star1 <- pobs(e1[, 1]) }
    
    rhos2 <- copula::P2p(cor(M_star))
    rhos <- c(rhos1, rep(0, dim(M1)[2] - dim(P1)[2]), rhos2)
    
    U1 <- as.matrix(cbind(e1_star1, M_star1))
    d2 <- sum(copula::dCopula(copula = copula::normalCopula(param = rhos, dim = dim(U1)[2], 
                                                            dispstr = "un"), u = U1, 
                              log = TRUE))

    d3 <- (d1 + d2)*(-1)
    
  } else { d3 <- 1e12 }
  
  return(d3)
  
}
likelihood_b2 <- function(values, M1, P1, Y1, M_star, M_star1) {
  
  betas1 <- values[1:(dim(M1)[2] + 1)]
  rhos1 <- values[(dim(M1)[2] + 2):(dim(M1)[2] + dim(P1)[2] + 1)]
  s1 <- values[length(values)]
  
  if ((sum(rhos1 > -1) == length(rhos1)) && ((sum(rhos1 < 1) == length(rhos1))) && (sum(s1 > 0) == length(s1))) {
    
    e1 <- Y1 - t(betas1%*%t(cbind(rep(1, dim(M1)[1]), M1)))
    d1 <- sum(dnorm(x = e1[, 1], mean = 0, sd = sqrt(s1), log = TRUE))
    
    e1_star1 <- pnorm(q = e1[, 1], mean = 0, sd = sqrt(s1))
    if (max(e1_star1) == 1) { e1_star1 <- pobs(e1[, 1]) }
    
    rhos2 <- copula::P2p(cor(M_star))
    rhos <- c(rhos1, rep(0, dim(M1)[2] - dim(P1)[2]), rhos2)
    
    U1 <- as.matrix(cbind(e1_star1, M_star1))
    d2 <- sum(copula::dCopula(copula = copula::normalCopula(param = rhos, dim = dim(U1)[2], 
                                                            dispstr = "un"), u = U1, 
                              log = TRUE))
    
    d3 <- (d1 + d2)*(-1)
    
  } else { d3 <- 1e12 }
  
  return(d3)
  
}
bootst1 <- function(r1, tpp1) {
  
  data2 <- samples1[[r1]]

  if (length(f1) == 1) {
    
    # dependent variable
    Y1 <- with(data2, model.frame(formula, data = data2)[1])
    
    # endogenous regressor(s)
    P1 <- with(data2, stats::get_all_vars(f1Z, data = data2))
    P_star <- P1
    P_star1 <- P1
    if (ecdf == TRUE) { P_star1 <- copula::pobs(P1)} else {
      
      for (a in 1:dim(P1)[2]) {
        
        Fhat <- ks::kcde(P1[, a])
        P_star1[, a] <- predict(Fhat, x = P1[, a])
        
      }
      
    }
    P_star <- as.data.frame(qnorm(as.matrix(P_star1)))
    
    M1 <- P1
    M_star <- P_star
    M_star1 <- P_star1
    
    
    # MLE without intercept
    tpp3 <- tryCatch(optim(fn = likelihood_b1, par = Estimate, method = method,
                           control = list(maxit = 1000000),
                           M1 = M1, P1 = P1, Y1 = Y1, M_star = M_star, M_star1 = M_star1), 
                     error = function(e) NA)
    
    sds1 <- c(tpp3$par[1:dim(M1)[2]], 
              tpp3$par[(dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])],
              tpp3$par[length(tpp1$par)])
    
  } else {
    
    # dependent variable
    Y1 <- data2[, 2]
    
    # endogenous regressor(s)
    P1 <- with(data2, stats::get_all_vars(f1Z, data = data2))
    P_star <- P1
    P_star1 <- P1
    if (ecdf == TRUE) { P_star1 <- copula::pobs(P1)} else {
      
      for (a in 1:dim(P1)[2]) {
        
        Fhat <- ks::kcde(P1[, a])
        P_star1[, a] <- predict(Fhat, x = P1[, a])
        
      }
      
    }
    
    P_star <- as.data.frame(qnorm(as.matrix(P_star1)))
    
    # exogenous regressor(s)
    X1 <- data2[, -c(1:(dim(P1)[2] + 2))]
    X1 <- as.matrix(X1)
    X_star <- X1
    X_star1 <- X1
    
    M1 <- cbind(P1, X1)
    
    
    if (matrixcalc::is.singular.matrix(t(as.matrix(M1))%*%as.matrix(M1)) == FALSE) {
      
      if (ecdf == TRUE) { X_star1 <- copula::pobs(X1)} else {
        
        for (a in 1:dim(X1)[2]) {
          
          Fhat <- ks::kcde(X1[, a])
          X_star1[, a] <- predict(Fhat, x = X1[, a])
          
        }
        
      }
      
      X_star <- as.data.frame(qnorm(as.matrix(X_star1)))
      
      M_star <- cbind(P_star, X_star)
      M_star1 <- cbind(P_star1, X_star1)
      
      
      if (var(f1X_m[, 1]) != 0) {
        
        # MLE without intercept
        tpp3 <- tryCatch(optim(fn = likelihood_b1, par = Estimate, method = method,
                               control = list(maxit = 1000000),
                               M1 = M1, P1 = P1, Y1 = Y1, M_star = M_star, M_star1 = M_star1), 
                         error = function(e) NA)
        
        sds1 <- c(tpp3$par[1:dim(M1)[2]], 
                  tpp3$par[(dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])],
                  tpp3$par[length(tpp1$par)])
        
      } else {
        
        # MLE with intercept
        tpp3 <- tryCatch(optim(fn = likelihood_b2, par = Estimate, method = method,
                               control = list(maxit = 1000000),
                               M1 = M1, P1 = P1, Y1 = Y1, M_star = M_star, M_star1 = M_star1), 
                         error = function(e) NA)
        
        sds1 <- c(tpp3$par[1:(dim(M1)[2] + 1)], 
                 tpp3$par[(dim(M1)[2] + 2):(dim(M1)[2] + dim(P1)[2] + 1)],
                 tpp3$par[length(tpp1$par)])
        
      }
      
    }
    
  }
  
  return(sds1)
  
}

CopRegML_par <- function(formula, index, data, ecdf = TRUE, nboots = 199,
                         starting_values = NULL, method = "Nelder-Mead") {

  # selecting data
  panelvar <- index[1]
  timevar <- index[2]
  data <- as.data.frame(data[with(data, order(seq(from = 1, to = dim(data)[1]))), ])
  
  # seperate endogenous and exogenous regressor(s)
  f1 <- nlme::splitFormula(formula, sep = "|")
  
  # check if arguments are correctly specified
  if (is.data.frame(data) == FALSE) {
    
    stop("Data must be data.frame", call. = FALSE)
    
  }
  
  
  if (length(f1) == 1) {
    
    f1Z <- f1[[1]]
    
    # Check if all variables exist in the data
    variables <- all.vars(f1Z)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {
      
      stop(paste("The following variables are missing in the data:", paste(missing_vars, collapse=", ")))
      
    }
    
    suppressMessages(attach(data))
    data <- data %>% dplyr::select(c(all_of(panelvar), colnames(model.frame(formula)[1]), colnames(get_all_vars(f1Z))))
    data <- data[complete.cases(data), ]
    data <- as.data.frame(data)
    suppressMessages(detach(data))
    
    aux_data <- do.call(rbind.data.frame, base::split(data, data[, which(colnames(data) == panelvar)]))
    data <- do.call(rbind.data.frame, lapply(base::split(data, data[, which(colnames(data) == panelvar)]), FOD1))
    data[, which(colnames(data) == panelvar)] <- aux_data[, which(colnames(aux_data) == panelvar)]
    data <- data[complete.cases(data), ]
    
    suppressMessages(attach(data))
    
    
    # dependent variable
    Y1 <- model.frame(formula)[1]
    
    # endogenous regressor(s)
    P1 <- stats::get_all_vars(f1Z)
    P_star <- P1
    P_star1 <- P1
    if (ecdf == TRUE) { P_star1 <- copula::pobs(P1)} else {
      
      for (a in 1:dim(P1)[2]) {
        
        Fhat <- ks::kcde(P1[, a])
        P_star1[, a] <- predict(Fhat, x = P1[, a])
        
      }
      
    }
    P_star <- as.data.frame(qnorm(as.matrix(P_star1)))
    
    M1 <- P1
    M_star <- P_star
    M_star1 <- P_star1
    
    
    # MLE without intercept
    lm0 <- lm(Formula::as.Formula(formula))
    lm1 <- lm(Formula::as.Formula(stats::update.formula(stats::formula(lm0), ~ . - 1)))
    starts <- c(lm1$coefficients, rep(0, dim(P1)[2]), log(var(lm1$residuals)))
    
    if (is.null(starting_values)) {starts1 <- starts} else {
      
      if (length(starting_values) != length(starts)) {
        
        stop("Vector of starting values is of wrong length", call. = FALSE)
        
      }
      
      if (is.numeric(starting_values) == FALSE) {
        
        stop("Starting values must be a numeric vector", call. = FALSE)
        
      }
      
      starts1 <- starting_values
      names(starts1) <- names(starts)
      
    }
    
    tpp1 <- tryCatch(optim(fn = likelihood1, par = starts1, method = method,
                           control = list(maxit = 1000000),
                           M1 = M1, P1 = P1, Y1 = Y1, M_star = M_star, M_star1 = M_star1), 
                     error = function(e) NA)
    
    if (sum(is.na(tpp1)) == 1) {
      
      stop("Model cannot be evaluated at initial values.", call. = FALSE)
      
    }
    
    if (tpp1$value == 1e+12) {
      
      stop("Model cannot be evaluated at initial values.", call. = FALSE)
      
    }
    
    
    Estimate <- c(tpp1$par[1:dim(M1)[2]], 
                  tanh(tpp1$par[(dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])]),
                  exp(tpp1$par[length(tpp1$par)]))
    resids <- unlist(Y1 - t(Estimate[1:dim(M1)[2]]%*%t(M1)))
    
    normtest_res <- tsoutliers::JarqueBera.test(resids)[[2]]$p.value
    z_vars <- colnames(stats::get_all_vars(f1Z))
    
    normtest_z <- rep(NA, dim(P1)[2])
    for (b1 in 1:length(normtest_z)) {
      
      normtest_z[b1] <- suppressWarnings(nortest::ad.test(scale(P1[, b1]))$p.value)
      
    }
    
    kstest_z <- rep(NA, dim(P1)[2])
    for (b2 in 1:length(kstest_z)) {
      
      kstest_z[b1] <- suppressWarnings(stats::ks.test(scale(resids), scale(P1[, b2]))$p.value)
      
    }
    
    names(Estimate)[length(Estimate)] <- "sigma2"
    for (b3 in (dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])){
      
      names(Estimate)[b3] <- paste("rho_", paste(z_vars[b3 - dim(M1)[2]], collapse = ""), sep = "")
      
    }
    
    
  } else {
    
    f1Z <- f1[[1]]
    f1X <- f1[[2]]
    
    variables <- all.vars(f1Z)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {
      
      stop(paste("The following variables are missing in the data:", paste(missing_vars, collapse=", ")))
      
    }
    
    variables <- all.vars(f1X)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {
      
      stop(paste("The following variables are missing in the data:", paste(missing_vars, collapse=", ")))
      
    }
    
    data1 <- data[data[, paste(timevar)] != max(data[, paste(timevar)], na.rm = TRUE), ]
    data1 <- data1 %>% dplyr::select(c(colnames(stats::get_all_vars(f1X, data = data1)), colnames(stats::get_all_vars(f1Z, data = data1))))
    data1 <- data1[stats::complete.cases(data1), ]
    f1X_m <- stats::model.matrix(f1X, data = data1)
    suppressMessages(attach(data))
    
    suppressWarnings(data <- data %>% dplyr::select(c(dplyr::all_of(panelvar), colnames(stats::model.frame(formula)[1]), colnames(stats::get_all_vars(f1Z)), colnames(stats::get_all_vars(f1X)))))
    data <- data[complete.cases(data), ]
    data <- as.data.frame(data)
    suppressMessages(detach(data))
    
    suppressMessages(attach(data))
    data <- cbind(data[, 1:(length(colnames(stats::get_all_vars(f1Z))) + 2)], stats::model.matrix(f1X)[, colnames(f1X_m)])
    if (var(f1X_m[, 1]) != 0) { colnames(data)[(length(colnames(stats::get_all_vars(f1Z))) + 3):dim(data)[2]] <- colnames(stats::model.matrix(f1X)) }
    data <- data %>% dplyr::select_if(resample::colVars(.) != 0)
    suppressMessages(detach(data))
    
    
    aux_data <- do.call(rbind.data.frame, base::split(data, data[, which(colnames(data) == panelvar)]))
    data <- do.call(rbind.data.frame, lapply(base::split(data, data[, which(colnames(data) == panelvar)]), FOD1))
    data[, which(colnames(data) == panelvar)] <- aux_data[, which(colnames(aux_data) == panelvar)]
    data <- data[stats::complete.cases(data), ]
    suppressMessages(attach(data))
    
    
    # check if design matrix has full column rank
    if (matrixcalc::is.singular.matrix(t(as.matrix(data[, -1]))%*%as.matrix(data[, -1]))) {
      
      stop("Design matrix is rank deficient. Either some regressors are perfectly collinear, or increase linearly over time, or are time invariant.", call. = FALSE)
      
    }
    
    
    # dependent variable
    Y1 <- data[, 2]
    
    # endogenous regressor(s)
    P1 <- stats::get_all_vars(f1Z)
    P_star <- P1
    P_star1 <- P1
    if (ecdf == TRUE) { P_star1 <- copula::pobs(P1)} else {
      
      for (a in 1:dim(P1)[2]) {
        
        Fhat <- ks::kcde(P1[, a])
        P_star1[, a] <- predict(Fhat, x = P1[, a])
        
      }
      
    }
    P_star <- as.data.frame(qnorm(as.matrix(P_star1)))
    
    # exogenous regressor(s)
    X1 <- data[, -c(1:(dim(P1)[2] + 2))]
    X1 <- as.matrix(X1)
    X_star <- X1
    X_star1 <- X1
    if (ecdf == TRUE) { X_star1 <- copula::pobs(X1)} else {
      
      for (a in 1:dim(X1)[2]) {
        
        Fhat <- ks::kcde(X1[, a])
        X_star1[, a] <- predict(Fhat, x = X1[, a])
        
      }
      
    }
    X_star <- as.data.frame(qnorm(as.matrix(X_star1)))
    
    M1 <- cbind(P1, X1)
    M_star <- cbind(P_star, X_star)
    M_star1 <- cbind(P_star1, X_star1)
    
    
    # check for enough observations
    if (dim(data)[2] + dim(P_star)[2] - 1 > dim(data)[1]) {
      
      stop("Not enough observations", call. = FALSE)
      
    }
    
    
    
    if (var(f1X_m[, 1]) != 0) {
      
      # MLE without intercept
      lm0 <- lm(Formula::as.Formula(formula))
      lm1 <- lm(Formula::as.Formula(stats::update.formula(stats::formula(lm0), ~ . - 1)))
      starts <- c(lm1$coefficients, rep(0, dim(P1)[2]), log(var(lm1$residuals)))
      
      if (is.null(starting_values)) {starts1 <- starts} else {
        
        if (length(starting_values) != length(starts)) {
          
          stop("Vector of starting values is of wrong length", call. = FALSE)
          
        }
        
        if (is.numeric(starting_values) == FALSE) {
          
          stop("Starting values must be a numeric vector", call. = FALSE)
          
        }
        
        starts1 <- starting_values
        names(starts1) <- names(starts)
        
      }
      
      tpp1 <- tryCatch(optim(fn = likelihood1, par = starts1, method = method,
                             control = list(maxit = 1000000),
                             M1 = M1, P1 = P1, Y1 = Y1, M_star = M_star, M_star1 = M_star1), 
                       error = function(e) NA)
      
      if (sum(is.na(tpp1)) == 1) {
        
        stop("Model cannot be evaluated at initial values.", call. = FALSE)
        
      }
      
      if (tpp1$value == 1e+12) {
        
        stop("Model cannot be evaluated at initial values.", call. = FALSE)
        
      }
      
      
      Estimate <- c(tpp1$par[1:dim(M1)[2]], 
                    tanh(tpp1$par[(dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])]),
                    exp(tpp1$par[length(tpp1$par)]))
      resids <- unlist(Y1 - t(Estimate[1:dim(M1)[2]]%*%t(M1)))
      
      normtest_res <- tsoutliers::JarqueBera.test(resids)[[2]]$p.value
      z_vars <- colnames(stats::get_all_vars(f1Z))
      
      normtest_z <- rep(NA, dim(P1)[2])
      for (a1 in 1:length(normtest_z)) {
        
        normtest_z[a1] <- suppressWarnings(nortest::ad.test(scale(P1[, a1]))$p.value)
        
      }
      
      kstest_z <- rep(NA, dim(P1)[2])
      for (a2 in 1:length(kstest_z)) {
        
        kstest_z[a2] <- suppressWarnings(stats::ks.test(scale(resids), scale(P1[, a2]))$p.value)
        
      }
      
      names(Estimate)[length(Estimate)] <- "sigma2"
      for (a3 in (dim(M1)[2] + 1):(dim(M1)[2] + dim(P1)[2])){
        
        names(Estimate)[a3] <- paste("rho_", paste(z_vars[a3 - dim(M1)[2]], collapse = ""), sep = "")
        
      }
      
    } else {
      
      # MLE with intercept
      lm1 <- lm(data[, -1])
      starts <- c(lm1$coefficients, rep(0, dim(P1)[2]), log(var(lm1$residuals)))
      
      if (is.null(starting_values)) {starts1 <- starts} else {
        
        if (length(starting_values) != length(starts)) {
          
          stop("Vector of starting values is of wrong length", call. = FALSE)
          
        }
        
        if (is.numeric(starting_values) == FALSE) {
          
          stop("Starting values must be a numeric vector", call. = FALSE)
          
        }
        
        starts1 <- starting_values
        names(starts1) <- names(starts)
        
      }
      
      tpp1 <- tryCatch(optim(fn = likelihood2, par = starts1, method = method,
                             control = list(maxit = 1000000),
                             M1 = M1, P1 = P1, Y1 = Y1, M_star = M_star, M_star1 = M_star1), 
                       error = function(e) NA)
      
      if (sum(is.na(tpp1)) == 1) {
        
        stop("Model cannot be evaluated at initial values.", call. = FALSE)
        
      }
      
      if (tpp1$value == 1e+12) {
        
        stop("Model cannot be evaluated at initial values.", call. = FALSE)
        
      }
      
      
      Estimate <- c(tpp1$par[1:(dim(M1)[2] + 1)], 
                    tanh(tpp1$par[(dim(M1)[2] + 2):(dim(M1)[2] + dim(P1)[2] + 1)]),
                    exp(tpp1$par[length(tpp1$par)]))
      resids <- unlist(Y1 - t(Estimate[1:(dim(M1)[2] + 1)]%*%t(cbind(rep(1, dim(M1)[1]), M1))))
      
      normtest_res <- tsoutliers::JarqueBera.test(resids)[[2]]$p.value
      z_vars <- colnames(stats::get_all_vars(f1Z))
      
      normtest_z <- rep(NA, dim(P1)[2])
      for (a4 in 1:length(normtest_z)) {
        
        normtest_z[a4] <- suppressWarnings(nortest::ad.test(scale(P1[, a4]))$p.value)
        
      }
      
      kstest_z <- rep(NA, dim(P1)[2])
      for (a5 in 1:length(kstest_z)) {
        
        kstest_z[a5] <- suppressWarnings(stats::ks.test(scale(resids), scale(P1[, a5]))$p.value)
        
      }
      
      names(Estimate)[length(Estimate)] <- "sigma2"
      for (a6 in (dim(M1)[2] + 2):(dim(M1)[2] + dim(P1)[2] + 1)){
        
        names(Estimate)[a6] <- paste("rho_", paste(z_vars[a6 - (dim(M1)[2] + 1)], collapse = ""), sep = "")
        
      }
      
    }
    
  }
  
  suppressMessages(detach(data))
  
  
  ################################ Bootstrapping ###############################
  
  if (nboots > 2) {
    
    if (dim(data)[1] > length(base::unique(data[, 1]))) {
      
      # bootstrap standard errors
      m <- 1
      samples1 <- list()
      
      while (m <= nboots) {
        
        data1 <- data %>% dplyr::group_by(data[, which(colnames(data) == panelvar)]) %>% dplyr::slice_sample(prop = 1, replace = TRUE)
        data1 <- as.data.frame(data1[, -dim(data1)[2]])
        suppressMessages(attach(data1))
        
        
        if (length(f1) == 1) {
          
          # dependent variable
          Y1 <- model.frame(formula)[1]
          
          # endogenous regressor(s)
          P1 <- stats::get_all_vars(f1Z)
          P_star <- P1
          P_star1 <- P1
          if (ecdf == TRUE) { P_star1 <- copula::pobs(P1)} else {
            
            for (a in 1:dim(P1)[2]) {
              
              Fhat <- ks::kcde(P1[, a])
              P_star1[, a] <- predict(Fhat, x = P1[, a])
              
            }
            
          }
          P_star <- as.data.frame(qnorm(as.matrix(P_star1)))
          
          M1 <- P1
          M_star <- P_star
          M_star1 <- P_star1
          
          
          # MLE without intercept
          tpp2 <- tryCatch(likelihood_b1(values = Estimate, 
                                         M1 = M1, P1 = P1, Y1 = Y1, 
                                         M_star = M_star, M_star1 = M_star1), 
                           error = function(e) NA)
          
          suppressMessages(detach(data1))
          
          suppressWarnings(if (sum(is.na(tpp2)) != 1) {
            
            if (tpp2 != 1e12) {
              
              samples1[[m]] <- data1
              
              m <- m + 1
              
            } else { m <- m }
            
          } else { m <- m })
          
          
        } else {
          
          # dependent variable
          Y1 <- data1[, 2]
          
          # endogenous regressor(s)
          P1 <- stats::get_all_vars(f1Z)
          P_star <- P1
          P_star1 <- P1
          if (ecdf == TRUE) { P_star1 <- copula::pobs(P1)} else {
            
            for (a in 1:dim(P1)[2]) {
              
              Fhat <- ks::kcde(P1[, a])
              P_star1[, a] <- predict(Fhat, x = P1[, a])
              
            }
            
          }
          P_star <- as.data.frame(qnorm(as.matrix(P_star1)))
          
          # exogenous regressor(s)
          X1 <- data1[, -c(1:(dim(P1)[2] + 2))]
          X1 <- as.matrix(X1)
          X_star <- X1
          X_star1 <- X1
          
          M1 <- cbind(P1, X1)
          
          
          if (matrixcalc::is.singular.matrix(t(as.matrix(M1))%*%as.matrix(M1)) == FALSE) {
            
            if (ecdf == TRUE) { X_star1 <- copula::pobs(X1)} else {
              
              for (a in 1:dim(X1)[2]) {
                
                Fhat <- ks::kcde(X1[, a])
                X_star1[, a] <- predict(Fhat, x = X1[, a])
                
              }
              
            }
            X_star <- as.data.frame(qnorm(as.matrix(X_star1)))
            
            M_star <- cbind(P_star, X_star)
            M_star1 <- cbind(P_star1, X_star1)
            
            
            if (var(f1X_m[, 1]) != 0) {
              
              # MLE without intercept
              tpp2 <- tryCatch(likelihood_b1(values = Estimate,
                                             M1 = M1, P1 = P1, Y1 = Y1, 
                                             M_star = M_star, M_star1 = M_star1), 
                               error = function(e) NA)
              
              suppressMessages(detach(data1))
              
              suppressWarnings(if (sum(is.na(tpp2)) != 1) {
                
                if (tpp2 != 1e12) {
                  
                  samples1[[m]] <- data1
                  
                  m <- m + 1
                  
                } else { m <- m }
                
              } else { m <- m })
              
              
            } else {
              
              # MLE with intercept
              tpp2 <- tryCatch(likelihood_b2(values = Estimate,
                                             M1 = M1, P1 = P1, Y1 = Y1, 
                                             M_star = M_star, M_star1 = M_star1), 
                               error = function(e) NA)
              
              suppressMessages(detach(data1))
              
              suppressWarnings(if (sum(is.na(tpp2)) != 1) {
                
                if (tpp2 != 1e12) {
                  
                  samples1[[m]] <- data1
                  
                  m <- m + 1
                  
                } else { m <- m }
                
              } else { m <- m })
              
            }
            
          } else { m <- m }
          
        }
        
        try(suppressMessages(detach(data1)), silent = TRUE)
        
      }
      
      # Parallelisation
      print("calculating bootstrap standard errors")
      core1 <- detectCores()
      c1 <- makeCluster(core1 - 1)
      
      if (length(f1) == 1) {
        
        clusterExport(cl = c1, varlist = c("samples1", "f1", "bootst1", "tpp1",
                                           "likelihood_b1", "likelihood_b2",
                                           "formula", "f1Z", "ecdf", "Estimate",
                                           "pobs", "method"),
                      envir = environment())
        
      } else {
        
        clusterExport(cl = c1, varlist = c("samples1", "f1", "bootst1", 
                                           "likelihood_b1", "likelihood_b2", 
                                           "formula", "f1Z", "f1X", "f1X_m", 
                                           "ecdf", "Estimate", "tpp1", "pobs",
                                           "method"),
                      envir = environment())
      }
      
      # Calculate bootstrap SE
      sds_mat <- pbsapply(X = 1:nboots, FUN = bootst1, cl = c1, tpp1 = tpp1)
      stopCluster(c1)
      
      
      # Identification checks
      if (normtest_res < .1) {warning("Residuals may not be symmetrically distributed: Jarque-Bera skewness p = ", 
                                      paste(round(normtest_res, digits = 3), 
                                            collapse = ""), call. = FALSE)}
      
      for (a7 in 1:length(normtest_z)) {
        
        if (normtest_z[a7] > .1) {warning("Endogenous regressor ", paste(z_vars[a7], collapse = ""), " may not be sufficiently different from normality: Anderson-Darling p = ", 
                                          paste(round(normtest_z[a7], digits = 3), 
                                                collapse = ""), call. = FALSE)}
        
      }
      
      for (a8 in 1:length(kstest_z)) {
        
        if (kstest_z[a8] > .1) {warning("Difference between endogenous regressor ", paste(z_vars[a8], collapse = ""), " and error distribution may not be sufficient for identification: Kolmogorov-Smirnov p = ", 
                                        paste(round(kstest_z[a8], digits = 3), 
                                              collapse = ""), call. = FALSE)}
        
      }
      
      
      # close(progress_bar)
      Std.Error <- sapply(as.data.frame(t(sds_mat)), sd)
      return(cbind(Estimate, Std.Error))
      
      
    } else {
      
      warning("Not enough observations within panels. Cannot calculate standard errors",
              call. = FALSE)
      return(Estimate)
      
    }
    
  } else { return(Estimate) }
  
}


################################################################################
#################################### Case 1 ####################################
################################################################################


# Dimensions
N <- 100
# N <- 200
Ti <- 5
# Ti <- 10


# Panel structure
testdata <- as.data.frame(matrix(NA, nrow = N*Ti, ncol = 0))
testdata$pan <- rep(x = seq(from = 1, to = N, by = 1), times = Ti)
testdata$year <- NA
testdata$year <- rep(seq(1:Ti), each = N)


# Individual-specific effects
testdata$re <- NA
testdata$re <- rep(times = Ti, runif(n = N, min = 0, max = 1))


# Explanatory variable and error
testdata$x <- NA
testdata$e <- NA

sigma_mat1 <- matrix(c(1, .5, .5, 1), ncol = 2, nrow = 2)

eps <- rmvnorm(n = Ti*N, mean = c(0, 0),
               sigma = sigma_mat1, method = "chol")

testdata$e <- qnorm(pnorm(eps[, 1]))
testdata$x <- qlnorm(p = pnorm(eps[, 2]), meanlog = 0, sdlog = 1)
# testdata$x <- qchisq(p = pnorm(eps[, 2]), df = 2)
# testdata$x <- qt(p = pnorm(eps[, 2]), df = 4)
testdata$x <- scale(testdata$x)


# Generate Y
testdata$Y <- NA
testdata$Y <- testdata$re + testdata$x + testdata$e


# Estimation
CopRegML_par(formula = Y ~ x, index = c("pan", "year"), data = testdata)
CopRegML_par(formula = Y ~ x, index = c("pan", "year"), data = testdata,
             starting_values = c(0, 0, 1))


rm(list = ls())



################################################################################
#################################### Case 2 ####################################
################################################################################


# Dimensions
N <- 100
# N <- 200
Ti <- 5
# Ti <- 10


# Degree of endogeneity
pcorr <- 0
# pcorr <- .8


# Panel structure
testdata <- as.data.frame(matrix(NA, nrow = N*Ti, ncol = 0))
testdata$pan <- rep(x = seq(from = 1, to = N, by = 1), times = Ti)
testdata$year <- NA
testdata$year <- rep(seq(1:Ti), each = N)


# Individual-specific effects
testdata$re <- NA
testdata$re <- rep(times = Ti, runif(n = N, min = 0, max = 1))


# Explanatory variables
testdata$x <- NA
testdata$z <- NA
testdata$z <- rnorm(n = N*Ti, mean = testdata$re, sd = 1)
testdata$z <- scale(testdata$z)

testdata$e <- NA
testdata$inst1 <- NA
testdata$inst1 <- NA

sigma_mat1 <- matrix(c(1, pcorr, 0, 0, pcorr, 1, .5, .3, 0, .5, 1, 0, 0, .3, 0, 
                       1), ncol = 4, nrow = 4)

eps <- rmvnorm(n = Ti*N, mean = c(0, 0, 0, 0),
               sigma = sigma_mat1, method = "chol")

testdata$e <- qnorm(pnorm(eps[, 1]))
testdata$x <- qlnorm(p = pnorm(eps[, 2]), meanlog = 0, sdlog = 1)
testdata$x <- scale(testdata$x)
testdata$inst1 <- qnorm(pnorm(eps[, 3]))
testdata$inst2 <- qnorm(pnorm(eps[, 4]))


# Generate Y
testdata$Y <- NA
testdata$Y <- testdata$re + testdata$x + testdata$z + testdata$e


# Estimation
CopRegML_par(formula = Y ~ x | z + -1, index = c("pan", "year"), data = testdata)

CopRegML_par(formula = Y ~ x | z + -1, index = c("pan", "year"), data = testdata,
             starting_values = c(0, 0, 0, 1))


rm(list = ls())



################################################################################
#################################### Case 3 ####################################
################################################################################


# Dimensions
N <- 100
# N <- 200
Ti <- 5
# Ti <- 10


# Correlations
pcorr <- .4
qcorr <- .3
# qcorr <- .6


# Panel structure
testdata <- as.data.frame(matrix(NA, nrow = N*Ti, ncol = 0))
testdata$pan <- rep(x = seq(from = 1, to = N, by = 1), times = Ti)
testdata$year <- NA
testdata$year <- rep(seq(1:Ti), each = N)


# Individual-specific effects
testdata$re <- NA
testdata$re <- rep(times = Ti, runif(n = N, min = 0, max = 1))


# Explanatory variables
testdata$x1 <- NA
testdata$x2 <- NA
testdata$z1 <- NA
testdata$z2 <- NA

testdata$e <- NA

sigma_mat1 <- matrix(c(1, pcorr, 0, pcorr, 0, pcorr, 1, qcorr, 0, 0, 0, qcorr, 
                       1, 0, 0, pcorr, 0, 0, 1, qcorr, 0, 0, 0, qcorr, 1),
                     ncol = 5, nrow = 5)

eps <- rmvnorm(n = N*Ti, mean = c(0, 0, 0, 0, 0),
               sigma = sigma_mat1, method = "chol")

testdata$e <- qnorm(pnorm(eps[, 1]))
testdata$x1 <- qchisq(p = pnorm(eps[, 2]), df = 2)
testdata$x1 <- scale(testdata$x1)
testdata$z1 <- qnorm(pnorm(eps[, 3]))
testdata$z1 <- scale(testdata$z1)
testdata$x2 <- qchisq(p = pnorm(eps[, 4]), df = 2)
testdata$x2 <- scale(testdata$x2)
testdata$z2 <- ifelse(eps[, 5] < 0, 0, 1)
testdata$z2 <- scale(testdata$z2)


# Generate Y
testdata$Y <- NA
testdata$Y <- testdata$re + testdata$x1 + testdata$z1 + testdata$x2 + testdata$z2 + testdata$e


# Estimation
CopRegML_par(formula = Y ~ x1 + x2 | z1 + z2 -1, index = c("pan", "year"), 
             data = testdata)
CopRegML_par(formula = Y ~ x1 + x2 | z1 + z2 -1, index = c("pan", "year"), 
             data = testdata, starting_values = c(0, 0, 0, 0, 0, 0, 1))


rm(list = ls())



################################################################################
#################################### Case 4 ####################################
################################################################################


# Dimensions
N <- 100
# N <- 200
Ti <- 5
# Ti <- 10


# Serial correlation
pcorr <- 0
# pcorr <- -.4


# Panel structure
testdata <- as.data.frame(matrix(NA, nrow = N*(Ti + 1), ncol = 0))
testdata$pan <- rep(x = seq(from = 1, to = N, by = 1), times = Ti + 1)
testdata$year <- NA
testdata$year <- rep(seq(1:(Ti + 1)), each = N)


# Individual-specific effects
testdata$re <- NA
testdata$re <- rep(times = Ti + 1, runif(n = N, min = 2, max = 3))

testdata$e <- NA
testdata$x <- NA

# Initial period
testdata$e[testdata$year == 1] <- rnorm(N)

testdata$Y <- NA
testdata$Y[testdata$year == 1] <- testdata$re[testdata$year == 1] + testdata$e[testdata$year == 1]


# Dynamic process
for (h1 in 2:(Ti + 1)) {
  
  testdata$e[testdata$year == h1] <- rnorm(N) + pcorr*testdata$e[testdata$year == (h1 - 1)]
  testdata$x[testdata$year == h1] <- testdata$Y[testdata$year == (h1 - 1)]
  testdata$Y[testdata$year == h1] <- testdata$re[testdata$year == h1] + .5*testdata$x[testdata$year == h1] + testdata$e[testdata$year == h1]
  
}

testdata <- testdata[complete.cases(testdata),]


# Estimation
CopRegML_par(formula = Y ~ x, index = c("pan", "year"), data = testdata)
CopRegML_par(formula = Y ~ x, index = c("pan", "year"), data = testdata, 
             starting_values = c(0, 0, 1))


rm(list = ls())


