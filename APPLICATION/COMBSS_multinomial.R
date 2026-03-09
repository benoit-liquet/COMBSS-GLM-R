library(glmnet)


# ================================================================
# Run COMBSS for a SINGLE subset size k
# Matches the delta schedule and adaptive alpha from
# COMBSS_multinomial() above.
#
# delta_i = min(delta_min * r^i, delta_max), i = 1..N
#   where r = (delta_max / delta_min)^(1 / N)
# Algorithm always runs exactly N iterations.
# ================================================================

COMBSS_multinomial_single_k <- function(X, y, k,
                                        delta_min = NULL,
                                        delta_max = NULL,
                                        Niter = 50,
                                        lambda = 0,
                                        alpha=0.01,
                                        mandatory = NULL,
                                        Xtest = NULL, ytest = NULL,
                                        normalize = TRUE,
                                        lambda_refit = 0,
                                        refit_method = c("glmnet","nnet"),
                                        verbose = TRUE) {
  
  refit_method <- match.arg(refit_method)
  X <- as.matrix(X)
  y <- as.factor(y)
  n <- nrow(X)
  p <- ncol(X)
  
  # Column normalisation [2]
  if (normalize) {
    sc <- scale_design_matrix(X)
    X_norm <- sc$X_scaled
  } else {
    X_norm <- X
  }
  
  # Mandatory variables setup
  if (!is.null(mandatory)) {
    mand_idx <- mandatory
    sel_idx  <- setdiff(1:p, mand_idx)
  } else {
    mand_idx <- integer(0)
    sel_idx  <- 1:p
  }
  p_sel <- length(sel_idx)
  
  # Penalty schedule calibration [2]
  cal <- calibrate_penalty_schedule(X_norm[, sel_idx, drop = FALSE],
                                    n, family = "multinomial")
  if (is.null(delta_min)) delta_min <- cal$delta_min
  if (is.null(delta_max)) delta_max <- cal$delta_max
  r <- (delta_max / delta_min)^(1 / Niter)
  
  if (verbose) {
    cat(sprintf("  Calibration: nu_max=%.4f | delta_conc=%.6f | delta_min=%.6f | delta_max=%.6f | r=%.6f\n",
                cal$nu_max, cal$delta_conc, delta_min, delta_max, r))
    cat(sprintf("\n--- k = %d | lambda = %.6f ---\n", k, lambda))
  }
  
  # Initialise t at centroid of T_k [2]
  t <- rep(k / p_sel, p_sel)
  
  # Homotopy loop: i = 1, ..., N (Algorithm 1) [2]
  s <- rep(0, p_sel)
  for (i in 1:(2*Niter)) {
    # Geometric delta schedule clamped at delta_max [2]
    delta <- min(delta_min * r^i, delta_max)
    
    # Adaptive step size: alpha = i/N [2]
    #alpha <- i / Niter
    
    # Compute Danskin gradient [2]
    g <- grad_danskin_multinomial(t, X_norm, y, delta, lambda,
                                  mand_idx, sel_idx)
    
    # Select k smallest gradient components -> FW vertex [2]
    k_smallest <- order(g)[1:k]
    s <- rep(0, p_sel)
    s[k_smallest] <- 1
    
    # Frank-Wolfe convex combination update [2]
    t <- (1 - alpha) * t + alpha * s
    
    # Truncate t to interior for numerical stability [2]
    t <- pmin(pmax(t, 0.0001), 0.9999)
    
    if (verbose && i %% 10 == 0) {
      cat(sprintf("  iter %3d | delta = %.6f | alpha = %.4f\n",
                  i, delta, alpha))
    }
  }
  
  # Selected features (map back to original indices)
  selected_sel  <- which(s == 1)
  selected_orig <- sel_idx[selected_sel]
  if (length(mand_idx) > 0) {
    selected_orig <- sort(c(mand_idx, selected_orig))
  }
  
  if (verbose) {
    cat(sprintf("  Selected: %s\n", paste(selected_orig, collapse = ", ")))
  }
  
  # Refit on original (unnormalised) data and evaluate [2][5]
  test_acc <- NA
  if (!is.null(Xtest) && !is.null(ytest)) {
    X_refit     <- X[, selected_orig, drop = FALSE]
    Xtest_refit <- as.matrix(Xtest)[, selected_orig, drop = FALSE]
    
    if (refit_method == "nnet") {
      df_train <- data.frame(y = y, X_refit)
      df_test  <- data.frame(Xtest_refit)
      colnames(df_test) <- colnames(df_train)[-1]
      
      # When lambda_refit = 0, use a small decay (1e-4) for
      # numerical stability in near-separable settings,
      # and always use maxit = 5000 to ensure convergence [5]
      nnet_decay <- ifelse(lambda_refit > 0, lambda_refit, 1e-4)
      
      refit  <- nnet::multinom(y ~ ., data = df_train, trace = FALSE,
                               MaxNWts = 10000, maxit = 5000,
                               decay = nnet_decay)
      y_pred <- as.character(predict(refit, newdata = df_test,
                                     type = "class"))
    } else {
      # glmnet: requires >= 2 columns; add zero column for k = 1
      if (ncol(X_refit) == 1) {
        X_refit     <- cbind(X_refit, 0)
        Xtest_refit <- cbind(Xtest_refit, 0)
      }
      if (lambda_refit > 0) {
        refit <- glmnet(x = X_refit, y = y, family = "multinomial",
                        type.multinomial = "grouped", alpha = 0,
                        standardize = TRUE, intercept = TRUE,
                        lambda = lambda_refit,
                        thresh = 1e-5, maxit = 1e5)
        lam_refit <- lambda_refit
      } else {
        refit <- glmnet(x = X_refit, y = y, family = "multinomial",
                        type.multinomial = "grouped", alpha = 0,
                        standardize = TRUE, intercept = TRUE,
                        nlambda = 50, thresh = 1e-5, maxit = 1e5)
        lam_refit <- min(refit$lambda)
      }
      y_pred <- as.character(predict(refit, newx = Xtest_refit,
                                     s = lam_refit,
                                     type = "class")[, 1])
    }
    
    test_acc <- mean(y_pred == as.character(ytest))
    if (verbose) cat(sprintf("  Test accuracy: %.3f\n", test_acc))
  }
  
  list(
    k             = k,
    lambda        = lambda,
    selected      = selected_orig,
    test_accuracy = test_acc,
    convergence   = Niter,
    delta_min     = delta_min,
    delta_max     = delta_max,
    r             = r
  )
}




###################################################################
# combss_multinomial.R
#
# Reference: Algorithm 1 with Frank-Wolfe homotopy
#
# Key difference from binary logistic:
#   - Coefficient vector beta becomes matrix B (p x C-1)
#   - Danskin gradient sums squared coefficients across classes:
#     df/dt_j = -2(lambda + delta) * sum_c(xi_{j,c}^2) / t_j^3
#   - glmnet with family="multinomial", type.multinomial="grouped"
#     enforces common sparsity pattern across classes
###################################################################



# ================================================================
# 1. COLUMN NORMALISATION [2]
# ================================================================

scale_design_matrix <- function(X) {
  v <- sqrt(colSums(X^2))
  # Avoid division by zero
  v[v == 0] <- 1
  X_scaled <- sweep(X, 2, v, "/")
  list(X_scaled = X_scaled, v = v)
}

# ================================================================
# 2. POWER ITERATION AND PENALTY SCHEDULE CALIBRATION [2]
#
# From Section 3 (Penalty Schedule Calibration) of the paper:
#   nu_max    = largest eigenvalue of X_u^T X_u (power iteration)
#   delta_conc = nu_max / (4n)  [multinomial]
#              = nu_max / (8n)  [logistic]
#   delta_max = delta_conc
#   delta_min = 1e-3 * delta_conc
#
# The growth rate r is NOT a user input; it is computed inside the
# main algorithm as:
#   r = (delta_max / delta_min)^(1 / N)
# so that delta grows from delta_min * r to delta_max in exactly N steps.
# ================================================================

compute_nu_max <- function(X, n_iter = 200, tol = 1e-8) {
  # Power iteration to estimate nu_max = lambda_max(X^T X)
  set.seed(42)
  p <- ncol(X)
  v <- rnorm(p)
  v <- v / sqrt(sum(v^2))
  
  nu_old <- 0
  for (i in 1:n_iter) {
    Xv   <- X %*% v           # n-vector
    XtXv <- crossprod(X, Xv)  # p-vector = X^T X v
    nu   <- sum(Xv^2)         # = ||X v||^2 = v^T X^T X v
    v    <- as.numeric(XtXv)
    v    <- v / sqrt(sum(v^2))
    if (abs(nu - nu_old) < tol * max(nu, 1)) break
    nu_old <- nu
  }
  return(nu)
}

calibrate_penalty_schedule <- function(X_u, n, family = "multinomial") {
  # X_u : n x (p-m) normalised matrix of relaxed (uncertain) variables
  # Returns delta_max = delta_conc and delta_min = 1e-3 * delta_conc.
  # r is NOT returned; it is derived from delta_max, delta_min, N inside
  # the main algorithm as r = (delta_max / delta_min)^(1 / N).
  nu_max     <- compute_nu_max(X_u)
  delta_conc <- if (family == "multinomial") nu_max / (4 * n) else nu_max / (8 * n)
  delta_max  <- delta_conc
  delta_min  <- 1e-3 * delta_conc
  
  list(nu_max = nu_max, delta_conc = delta_conc,
       delta_max = delta_max, delta_min = delta_min)
}

# ================================================================
# 3. INNER SOLVER: Ridge-penalised multinomial GLM via glmnet [2]
#
# For interior t in (0,1)^(p-m), the inner problem is:
#   min_{xi_0, Xi} -ell/n + sum_j omega_j(t) * ||Xi[j,]||^2
#
# where omega_j(t) = (lambda + delta)/t_j^2 - delta  for relaxed variables
#       omega_j(t) = lambda                           for mandatory variables
#
# This is the correct Xi-parameterisation penalty from the paper:
#   sum_j [(lambda+delta)/t_j^2 - delta] * ||Xi[j,]||^2
# which in B-parameterisation (Xi = T_t B) gives delta*||Gamma_t B||^2.
#
# glmnet (alpha=0) minimises:
#   -ell/n + (lambda_glmnet/2) * sum_j pf_j * ||beta[j,]||_F^2
#
# Matching: lambda_glmnet/2 * pf_j = omega_j
#   => pf_j = p * omega_j / omega_sum
#   => lambda_glmnet = 2 * omega_sum / p
# ================================================================

solve_inner_multinomial <- function(t, X, y, delta, lambda,
                                    mand_idx = integer(0),
                                    sel_idx = 1:ncol(X)) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Compute omega_j for each variable [2]
  omega <- rep(lambda, p)
  
  # Mandatory variables: omega = lambda
  if (length(mand_idx) > 0) {
    omega[mand_idx] <- lambda
  }
  
  # Relaxed variables: omega_j = (lambda + delta) / t_j^2 - delta [2]
  for (j in seq_along(sel_idx)) {
    omega[sel_idx[j]] <- (lambda + delta) / (t[j]^2) - delta
  }
  
  # Map to glmnet parameterisation [2]
  omega_sum <- sum(omega)
  if (omega_sum < 1e-12) omega_sum <- 1e-12
  penalty_factors <- p * omega / omega_sum
  lambda_glmnet <- 2 * omega_sum / p
  
  # Solve inner ridge multinomial GLM [2]
  # Use a small lambda sequence to avoid single-lambda convergence failure
  lambda_seq <- c(lambda_glmnet * 1.1, lambda_glmnet)
  
  fit <- tryCatch({
    glmnet(
      x = X,
      y = y,
      family = "multinomial",
      type.multinomial = "grouped",
      alpha = 0,
      lambda = lambda_seq,
      penalty.factor = penalty_factors,
      standardize = FALSE,
      intercept = TRUE,
      thresh = 1e-5,
      maxit = 1e6
    )
  }, error = function(e) {
    # Fallback: let glmnet choose its own lambda path
    glmnet(
      x = X,
      y = y,
      family = "multinomial",
      type.multinomial = "grouped",
      alpha = 0,
      penalty.factor = penalty_factors,
      standardize = FALSE,
      intercept = TRUE,
      thresh = 1e-5,
      maxit = 1e6
    )
  })
  
  # Extract coefficients at the target lambda
  coef_list <- coef(fit, s = lambda_glmnet)
  C <- length(coef_list)
  beta0 <- sapply(coef_list, function(b) as.numeric(b[1]))
  Xi <- matrix(0, p, C)
  for (c in 1:C) {
    Xi[, c] <- as.numeric(coef_list[[c]][-1])
  }
  
  list(beta0 = beta0, Xi = Xi)
}

# ================================================================
# 4. DANSKIN GRADIENT for multinomial [2]
#
# df/dt_j = -2(lambda + delta) * sum_c(xi_{m+j,c}^2) / t_j^3
#
# This is the key multinomial extension: we sum the squared
# coefficients across all C classes for each variable j [2]
# ================================================================

grad_danskin_multinomial <- function(t, X, y, delta, lambda,
                                     mand_idx = integer(0),
                                     sel_idx = 1:ncol(X)) {
  
  result <- solve_inner_multinomial(t, X, y, delta, lambda,
                                    mand_idx, sel_idx)
  Xi <- result$Xi
  
  # Sum squared coefficients across classes for each selected variable [2]
  row_norms_sq <- rowSums(Xi[sel_idx, , drop = FALSE]^2)
  
  # Danskin gradient [2]
  grad_t <- -2 * (lambda + delta) * row_norms_sq / (t^3)
  
  return(as.numeric(grad_t))
}

# ================================================================
# 5. MAIN ALGORITHM: COMBSS for Multinomial Logistic Regression
#    Algorithm 1 with Frank-Wolfe homotopy [2]
#
# delta schedule: delta_i = min(delta_min * r^i, delta_max), i = 1..N
#   where r = (delta_max / delta_min)^(1 / N)
#   so delta grows from delta_min*r to delta_max over N steps.
#
# Adaptive step size: alpha_i = i/N  [paper Algorithm 1]
# Algorithm always runs exactly N iterations.
# ================================================================

COMBSS_multinomial <- function(X, y, Kmax,
                               delta_min = NULL, delta_max = NULL,
                               Niter = 50,
                               lambda = 0, alpha=0.01,
                               mandatory = NULL,
                               Xtest = NULL, ytest = NULL,
                               normalize = TRUE,
                               lambda_refit = 0,
                               refit_method = c("glmnet","nnet"),
                               verbose = TRUE) {
  refit_method <- match.arg(refit_method)
  
  X <- as.matrix(X)
  y <- as.factor(y)
  n <- nrow(X)
  p <- ncol(X)
  
  # Column normalisation [2]
  if (normalize) {
    sc <- scale_design_matrix(X)
    X_norm <- sc$X_scaled
    col_norms <- sc$v
  } else {
    X_norm <- X
    col_norms <- rep(1, p)
  }
  
  # Mandatory variables setup
  if (!is.null(mandatory)) {
    mand_idx <- mandatory
    sel_idx <- setdiff(1:p, mand_idx)
  } else {
    mand_idx <- integer(0)
    sel_idx <- 1:p
  }
  p_sel <- length(sel_idx)
  
  # Penalty schedule calibration [2]
  # Compute delta_max = delta_conc and delta_min via power iteration.
  # r is then derived as (delta_max/delta_min)^(1/(N-1)).
  cal <- calibrate_penalty_schedule(X_norm[, sel_idx, drop = FALSE],
                                    n, family = "multinomial")
  if (is.null(delta_min)) delta_min <- cal$delta_min
  if (is.null(delta_max)) delta_max <- cal$delta_max
  r <- (delta_max / delta_min)^(1 / Niter)
  if (verbose) {
    cat(sprintf("  Calibration: nu_max=%.4f | delta_conc=%.6f | delta_min=%.6f | delta_max=%.6f | r=%.6f\n",
                cal$nu_max, cal$delta_conc, delta_min, delta_max, r))
  }
  
  # Storage
  selected_models <- list()
  test_accuracy <- numeric(Kmax)
  convergence <- numeric(Kmax)
  best_model_matrix <- matrix(0, Kmax, p_sel)
  
  # ---- Loop over subset sizes k = 1, ..., Kmax ----
  for (k in 1:Kmax) {
    if (verbose) cat(sprintf("\n--- k = %d ---\n", k))
    
    # Initialise t at interior point: centroid of T_k
    t <- rep(k / p_sel, p_sel)
    
    # Homotopy loop: i = 1, ..., N  (Algorithm 1)
    s <- rep(0, p_sel)
    for (i in 1:(2*Niter)) {
      # Geometric delta schedule clamped at delta_max
      delta <- min(delta_min * r^i, delta_max)
      
      # Adaptive step size: alpha = i/N
      #alpha <- i / Niter
      
      # Compute Danskin gradient [2]
      g <- grad_danskin_multinomial(t, X_norm, y, delta, lambda,
                                    mand_idx, sel_idx)
      
      # Select k smallest gradient components -> FW vertex [2]
      k_smallest <- order(g)[1:k]
      s <- rep(0, p_sel)
      s[k_smallest] <- 1
      
      # Frank-Wolfe convex combination update [2]
      t <- (1 - alpha) * t + alpha * s
      
      # Truncate t to interior for numerical stability [2]
      t <- pmin(pmax(t, 0.0001), 0.9999)
      
      if (verbose && i %% 50 == 0) {
        cat(sprintf("  iter %3d | delta = %.6f | alpha = %.4f\n",
                    i, delta, alpha))
      }
    }
    
    convergence[k] <- Niter
    
    # Selected features (map back to original indices)
    selected_sel <- which(s == 1)
    selected_orig <- sel_idx[selected_sel]
    if (length(mand_idx) > 0) {
      selected_orig <- sort(c(mand_idx, selected_orig))
    }
    selected_models[[k]] <- selected_orig
    best_model_matrix[k, selected_sel] <- 1
    
    if (verbose) {
      cat(sprintf("  Selected: %s\n", paste(selected_orig, collapse = ", ")))
    }
    
    # Refit on original (unnormalised) data and evaluate [2]
    if (!is.null(Xtest) && !is.null(ytest)) {
      X_refit     <- X[, selected_orig, drop = FALSE]
      Xtest_refit <- as.matrix(Xtest)[, selected_orig, drop = FALSE]
      
      if (refit_method == "nnet") {
        df_train <- data.frame(y = y, X_refit)
        df_test  <- data.frame(Xtest_refit)
        colnames(df_test) <- colnames(df_train)[-1]
        
        # When lambda_refit = 0, use a small decay (1e-4) for
        # numerical stability in near-separable settings,
        # and always use maxit = 5000 to ensure convergence 
        nnet_decay <- ifelse(lambda_refit > 0, lambda_refit, 1e-4)
        
        refit  <- nnet::multinom(y ~ ., data = df_train, trace = FALSE,
                                 MaxNWts = 10000, maxit = 5000,
                                 decay = nnet_decay)
        y_pred <- as.character(predict(refit, newdata = df_test,
                                       type = "class"))
      } else {
        # glmnet: requires >= 2 columns; add zero column for k = 1 [2]
        if (ncol(X_refit) == 1) {
          X_refit     <- cbind(X_refit, 0)
          Xtest_refit <- cbind(Xtest_refit, 0)
        }
        if (lambda_refit > 0) {
          refit <- glmnet(x = X_refit, y = y, family = "multinomial",
                          type.multinomial = "grouped", alpha = 0,
                          standardize = TRUE, intercept = TRUE,
                          lambda = lambda_refit, thresh = 1e-5, maxit = 1e5)
          lam_refit <- lambda_refit
        } else {
          refit <- glmnet(x = X_refit, y = y, family = "multinomial",
                          type.multinomial = "grouped", alpha = 0,
                          standardize = TRUE, intercept = TRUE,
                          nlambda = 50, thresh = 1e-5, maxit = 1e5)
          lam_refit <- min(refit$lambda)
        }
        y_pred <- as.character(predict(refit, newx = Xtest_refit,
                                       s = lam_refit, type = "class")[, 1])
      }
      
      test_accuracy[k] <- mean(y_pred == as.character(ytest))
      
      if (verbose) {
        cat(sprintf("  Test accuracy: %.3f\n", test_accuracy[k]))
      }
    }
  }
  
  # Find best k
  best_k <- which.max(test_accuracy)
  
  list(
    selected_models = selected_models,
    best_model_matrix = best_model_matrix,
    best_k = best_k,
    best_subset = selected_models[[best_k]],
    test_accuracy = test_accuracy,
    convergence = convergence
  )
}
