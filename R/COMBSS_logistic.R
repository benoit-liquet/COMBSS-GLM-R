###################################################################
# COMBSS for Binary Logistic Regression
# Matching Algorithm 1 exactly as written in the paper  :
#   - for loop i = 1, ..., N (no convergence check)
#   - delta = min(delta_min * r^i, delta_max)
#   - r = (delta_max / delta_min)^(1/N)
###################################################################

library(glmnet)

# ================================================================
# 1. COLUMN NORMALISATION (Section 3.1)  
# ================================================================

scale_design_matrix <- function(X) {
  v <- sqrt(colSums(X^2))
  v[v == 0] <- 1
  X_scaled <- sweep(X, 2, v, "/")
  list(X_scaled = X_scaled, v = v)
}

# ================================================================
# 2. POWER ITERATION AND PENALTY SCHEDULE CALIBRATION  
# ================================================================

compute_nu_max <- function(X, n_iter = 200, tol = 1e-8) {
  set.seed(42)
  p <- ncol(X)
  v <- rnorm(p)
  v <- v / sqrt(sum(v^2))
  nu_old <- 0
  for (i in 1:n_iter) {
    Xv   <- X %*% v
    XtXv <- crossprod(X, Xv)
    nu   <- sum(Xv^2)
    v    <- as.numeric(XtXv)
    v    <- v / sqrt(sum(v^2))
    if (abs(nu - nu_old) < tol * max(nu, 1)) break
    nu_old <- nu
  }
  return(nu)
}

calibrate_penalty_schedule <- function(X_u, n, family = "logistic") {
  nu_max     <- compute_nu_max(X_u)
  delta_conc <- if (family == "multinomial") nu_max / (4 * n) else nu_max / (8 * n)
  delta_max  <- delta_conc
  delta_min  <- 1e-3 * delta_conc
  list(nu_max = nu_max, delta_conc = delta_conc,
       delta_max = delta_max, delta_min = delta_min)
}

# ================================================================
# 3. INNER RIDGE SOLVER via glmnet [2][4] 
# ================================================================

solve_inner_glm <- function(t, X_no_int, y, delta, lambda,
                            mand_idx = integer(0),
                            sel_idx = 1:ncol(X_no_int)) {
  n <- nrow(X_no_int)
  p <- ncol(X_no_int)
  
  omega <- rep(0, p)
  if (length(mand_idx) > 0) {
    omega[mand_idx] <- lambda
  }
  omega[sel_idx] <- (lambda + delta) / (t^2) - delta
  omega <- pmax(omega, 1e-12)
  
  omega_sum <- sum(omega)
  if (omega_sum < 1e-12) omega_sum <- 1e-12
  penalty_factors <- p * omega / omega_sum
  lambda_glmnet <- 2 * omega_sum / p
  
  fit <- glmnet(
    x = X_no_int, y = y, family = "binomial", alpha = 0,
    lambda = lambda_glmnet, penalty.factor = penalty_factors,
    intercept = TRUE, standardize = FALSE,
    thresh = 1e-12, maxit = 1e6
  )
  
  xi <- as.numeric(coef(fit))
  return(xi)
}

# ================================================================
# 4. DANSKIN GRADIENT (Eq. 12) [2] 
# ================================================================

grad_danskin <- function(t, X_no_int, y, delta, lambda,
                         mand_idx = integer(0),
                         sel_idx = 1:ncol(X_no_int)) {
  xi <- solve_inner_glm(t, X_no_int, y, delta, lambda, mand_idx, sel_idx)
  xi_sel <- xi[1 + sel_idx]
  grad_t <- -2 * (lambda + delta) * (xi_sel^2) / (t^3)
  return(as.numeric(grad_t))
}

# ================================================================
# 5. REFIT AND EVALUATE  
# ================================================================

refit_and_evaluate <- function(X_train, y_train, X_test, y_test,
                               selected_idx) {
  if (length(selected_idx) == 0) {
    majority <- as.numeric(names(which.max(table(y_train))))
    pred <- rep(majority, length(y_test))
    return(list(accuracy = mean(pred == y_test), pred = pred))
  }
  
  X_sub_train <- as.matrix(X_train[, selected_idx, drop = FALSE])
  X_sub_test  <- as.matrix(X_test[, selected_idx, drop = FALSE])
  
  if (length(selected_idx) == 1) {
    df_train <- data.frame(y = y_train, x = X_sub_train)
    fit <- glm(y ~ ., data = df_train, family = binomial())
    df_test <- data.frame(x = X_sub_test)
    pred_prob <- predict(fit, newdata = df_test, type = "response")
    pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  } else {
    refit <- glmnet(
      x = X_sub_train, y = y_train, family = "binomial",
      alpha = 0, lambda = 1e-8,
      intercept = TRUE, standardize = FALSE
    )
    pred_prob <- predict(refit, newx = X_sub_test,
                         type = "response", s = 1e-8)
    pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  }
  
  accuracy <- mean(pred_class == y_test)
  return(list(accuracy = accuracy, pred = as.numeric(pred_class)))
}

# ================================================================
# 6. MAIN ALGORITHM — Matching Algorithm 1 in the paper  
#
# Input:  lambda >= 0, k, delta_min, delta_max, N
# Init:   t <- (k/(p-m)) * 1, r = (delta_max/delta_min)^(1/N)
# For i = 1, ..., N:
#   delta <- min(delta_min * r^i, delta_max)
#   g <- envelope gradient of f_{delta,lambda}(t)
#   s <- k smallest components of g
#   t <- (1 - alpha) * t + alpha * s
# Fit on last vertex s
# ================================================================

COMBSS_logistic <- function(X, y, Kmax,
                            delta_min = NULL, delta_max = NULL,alpha=0.01,
                            Niter = 50,
                            lambda = 0,
                            mandatory = NULL,
                            Xtest = NULL, ytest = NULL,
                            normalize = TRUE,
                            verbose = TRUE) {
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (!is.null(Xtest)) Xtest <- as.matrix(Xtest)
  if (!is.null(ytest)) ytest <- as.numeric(ytest)
  n <- nrow(X)
  p <- ncol(X)
  
  # --- Mandatory / selectable indices ---
  if (is.null(mandatory)) {
    mand_idx <- integer(0)
    sel_idx  <- 1:p
  } else {
    mand_idx <- sort(mandatory)
    sel_idx  <- setdiff(1:p, mand_idx)
  }
  p_sel <- length(sel_idx)
  
  # --- Column normalisation (Section 3.1)   ---
  if (normalize) {
    sc <- scale_design_matrix(X)
    X_scaled <- sc$X_scaled
  } else {
    X_scaled <- X
  }
  
  # --- Penalty schedule calibration   ---
  cal <- calibrate_penalty_schedule(X_scaled[, sel_idx, drop = FALSE],
                                    n, family = "logistic")
  if (is.null(delta_min)) delta_min <- cal$delta_min
  if (is.null(delta_max)) delta_max <- cal$delta_max
  
  # r = (delta_max / delta_min)^(1/N)  
  r <- (delta_max / delta_min)^(1 / Niter)
  
  if (verbose) {
    cat(sprintf("Calibration: nu_max=%.4f | delta_conc=%.6f | delta_min=%.8f | delta_max=%.6f | r=%.6f\n",
                cal$nu_max, cal$delta_conc, delta_min, delta_max, r))
  }
  
  # --- Storage ---
  selected_models   <- list()
  best_model_matrix <- matrix(0, nrow = Kmax, ncol = p_sel)
  convergence_iters <- numeric(Kmax)
  
  # --- Main loop over subset sizes k = 1, ..., Kmax ---
  for (k in 1:Kmax) {
    if (verbose) cat(sprintf("\n--- k = %d ---\n", k))
    
    # Init: t <- (k/(p-m)) * 1   
    t <- rep(k / p_sel, p_sel)
    
    # For i = 1, ..., N   
    s <- rep(0, p_sel)
    for (i in 1:(2*Niter)) {
      
      # delta <- min(delta_min * r^i, delta_max)  
      delta <- min(delta_min * r^i, delta_max)
      
      # Compute envelope gradient g  
      g <- grad_danskin(t, X_scaled, y, delta, lambda,
                        mand_idx, sel_idx)
      
      # s <- k smallest components of g 
      s <- rep(0, p_sel)
      s[order(g)[1:k]] <- 1
      
      
      #alpha <- i / Niter
      
      # t <- (1 - alpha) * t + alpha * s  
      t <- (1 - alpha) * t + alpha * s
      
      # Clamp to interior for numerical stability
      t <- pmin(pmax(t, 1e-4), 1 - 1e-4)
      
      if (verbose && i %% 10 == 0) {
        cat(sprintf("  iter %3d | delta = %10.6f",
                    i, delta))
      }
    }
    
    convergence_iters[k] <- Niter
    
    # Map selected indices back to original X columns
    selected_sel  <- which(s == 1)
    selected_full <- sel_idx[selected_sel]
    if (length(mand_idx) > 0) {
      selected_full <- sort(c(mand_idx, selected_full))
    }
    
    selected_models[[k]]    <- selected_full
    best_model_matrix[k, ]  <- s
    
    if (verbose) {
      cat(sprintf("  Selected: %s\n",
                  paste(selected_full, collapse = ", ")))
    }
  }
  
  # --- Fit on last vertex s; evaluate on test set   ---
  test_accuracy <- numeric(Kmax)
  if (!is.null(Xtest) && !is.null(ytest)) {
    for (k in 1:Kmax) {
      idx <- selected_models[[k]]
      res <- refit_and_evaluate(X, y, Xtest, ytest, idx)
      test_accuracy[k] <- res$accuracy
    }
    best_k <- which.max(test_accuracy)
  } else {
    best_k <- Kmax
  }
  
  list(
    selected_models   = selected_models,
    best_model_matrix = best_model_matrix,
    best_k            = best_k,
    best_subset       = selected_models[[best_k]],
    test_accuracy     = test_accuracy,
    convergence       = convergence_iters
  )
}