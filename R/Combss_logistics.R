###################################################################
# COMBSS for Logistic Regression — Final Version
# Algorithm 1: Homotopy with Frank-Wolfe Steps [1][4]
# Gradient via Danskin's Envelope Theorem (Eq. 10)
###################################################################

library(glmnet)

# ================================================================
# INNER RIDGE SOLVER
#
# Solves: min_{xi_0, xi} g_{delta,lambda}(t, xi_0, xi)
#
# From Eq. (10) [1][4]:
# g = -(1/n)*ell(xi_0, xi; D)
#     + lambda * sum_{j=1}^{m} xi_j^2
#     + sum_{j=1}^{p-m} ((lambda+delta)/t_j^2 - delta) * xi_{m+j}^2
#
# For m=0 (no mandatory):
# omega_j = (lambda + delta) / t_j^2 - delta
#
# For mandatory variables (j = 1,...,m):
# omega_j = lambda
#
# glmnet (binomial, alpha=0) minimises:
#   -(1/n)*ell(xi) + (lambda_glmnet/2) * sum(pf_j * xi_j^2)
#
# Mapping: (lambda_glmnet/2) * pf_j = omega_j
# With glmnet normalisation sum(pf) = p:
#   pf_j = p * omega_j / sum(omega)
#   lambda_glmnet = 2 * sum(omega) / p
# ================================================================

solve_inner_glm <- function(t, X_no_int, y, delta, lambda,
                            mand_idx = integer(0), sel_idx = 1:ncol(X_no_int)) {
  n <- nrow(X_no_int)
  p <- ncol(X_no_int)
  
  # Build penalty weight vector omega (length p)
  omega <- rep(0, p)
  if (length(mand_idx) > 0) {
    omega[mand_idx] <- lambda
  }
  omega[sel_idx] <- (lambda + delta) / (t^2) - delta
  omega <- pmax(omega, 1e-12)
  
  # Map to glmnet parameterisation
  omega_sum <- sum(omega)
  if (omega_sum < 1e-12) omega_sum <- 1e-12
  penalty_factors <- p * omega / omega_sum
  lambda_glmnet <- 2 * omega_sum / p
  
  fit <- glmnet(
    x = X_no_int, y = y, family = "binomial", alpha = 0,
    lambda = lambda_glmnet, penalty.factor = penalty_factors,
    intercept = TRUE, standardize = FALSE, thresh = 1e-12,
    maxit = 1e6
  )
  
  xi <- as.numeric(coef(fit))  # (intercept, xi_1, ..., xi_p)
  return(xi)
}

# ================================================================
# DANSKIN GRADIENT (Eq. 10 in [1][4])
#
# d f / d t_j = -2(lambda + delta) * xi_{m+j}^2 / t_j^3
#
# Only computed for selectable variables (not mandatory)
# ================================================================

grad_danskin <- function(t, X_no_int, y, delta, lambda,
                         mand_idx = integer(0), sel_idx = 1:ncol(X_no_int)) {
  
  xi <- solve_inner_glm(t, X_no_int, y, delta, lambda, mand_idx, sel_idx)
  xi_sel <- xi[1 + sel_idx]  # +1 because coef() puts intercept first
  
  grad_t <- -2 * (lambda + delta) * (xi_sel^2) / (t^3)
  return(as.numeric(grad_t))
}

# ================================================================
# REFIT: Unpenalised (or lightly penalised) logistic regression
# on the selected subset, then evaluate on test set
# ================================================================

refit_and_evaluate <- function(X_train, y_train, X_test, y_test, selected_idx) {
  
  if (length(selected_idx) == 0) {
    majority <- as.numeric(names(which.max(table(y_train))))
    pred <- rep(majority, length(y_test))
    return(list(accuracy = mean(pred == y_test), pred = pred))
  }
  
  # Force matrix conversion
  X_sub_train <- as.matrix(X_train[, selected_idx, drop = FALSE])
  X_sub_test <- as.matrix(X_test[, selected_idx, drop = FALSE])
  
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
    pred_prob <- predict(refit, newx = X_sub_test, type = "response", s = 1e-8)
    pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  }
  
  accuracy <- mean(pred_class == y_test)
  return(list(accuracy = accuracy, pred = as.numeric(pred_class)))
}
# ================================================================
# MAIN ALGORITHM: COMBSS for Logistic Regression
# Algorithm 1 from [1][4]
#
# @param X        Design matrix (n x p), WITHOUT intercept
# @param y        Binary response (0/1)
# @param Kmax     Maximum subset size to search
# @param delta_min  Initial penalty (small, e.g. 0.1)
# @param r        Geometric growth factor for delta (r > 1)
# @param Niter    Maximum iterations per k
# @param alpha    Frank-Wolfe step size in (0,1)
# @param lambda   Ridge penalty (>= 0, default 0)
# @param epsilon  Corner tolerance for convergence
# @param mandatory  Integer vector of mandatory variable indices
#                   (1-based columns of X). NULL if none.
# @param Xtest    Test design matrix (WITHOUT intercept)
# @param ytest    Test response vector
# @param normalize  Whether to normalise columns (Section 3.1)
# ================================================================

COMBSS_logistic <- function(X, y, Kmax, delta_min = 0.1, r = 1.5,
                            Niter = 500, alpha = 0.05, lambda = 0,
                            epsilon = 0.01, mandatory = NULL,
                            Xtest = NULL, ytest = NULL,
                            normalize = TRUE) {
  
  # --- Force matrix types ---
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (!is.null(Xtest)) Xtest <- as.matrix(Xtest)
  if (!is.null(ytest)) ytest <- as.numeric(ytest)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # --- Mandatory / selectable indices ---
  if (is.null(mandatory)) {
    m <- 0
    mand_idx <- integer(0)
    sel_idx <- 1:p
  } else {
    mand_idx <- sort(mandatory)
    m <- length(mand_idx)
    sel_idx <- setdiff(1:p, mand_idx)
  }
  p_sel <- length(sel_idx)
  
  # --- Column normalisation (Section 3.1 in [1][4]) ---
  if (normalize) {
    col_norms <- apply(X, 2, function(col) sqrt(sum(col^2)))
    col_norms[col_norms == 0] <- 1
    X_scaled <- sweep(X, 2, col_norms, "/")
    if (!is.null(Xtest)) {
      Xtest_scaled <- sweep(Xtest, 2, col_norms, "/")
    }
  } else {
    X_scaled <- X
    if (!is.null(Xtest)) {
      Xtest_scaled <- Xtest
    }
  }
  
  # --- Storage ---
  selected_models <- list()
  best_model_matrix <- matrix(0, nrow = Kmax, ncol = p_sel)
  convergence_iters <- numeric(Kmax)
  
  # --- Main loop over subset sizes k = 1, ..., Kmax ---
  for (k in 1:Kmax) {
    
    # Initialisation: t = k/(p-m) * 1  (Algorithm 1, [1][4])
    t <- rep(k / p_sel, p_sel)
    
    converged <- FALSE
    iter <- 0
    
    while (!converged && iter < Niter) {
      iter <- iter + 1
      
      # Step 2: Geometric delta schedule
      delta_current <- delta_min * r^(iter - 1)
      
      # Step 3: Gradient via Danskin (Eq. 10, [1][4])
      g <- grad_danskin(t, X_scaled, y, delta_current, lambda,
                        mand_idx, sel_idx)
      
      # Step 4: FW vertex — select k smallest gradient components
      s <- rep(0, p_sel)
      s[order(g)[1:k]] <- 1
      
      # Step 5: Frank-Wolfe convex combination update
      t <- (1 - alpha) * t + alpha * s
      
      # Clamp to interior (numerical safety)
      t <- pmax(t, 1e-4)
      t <- pmin(t, 1 - 1e-4)
      
      # Step 7: Check convergence (||t - s||_inf <= epsilon)
      if (max(abs(t - s)) <= epsilon) {
        converged <- TRUE
      }
    }
    
    # Map selected indices back to original X columns
    selected_sel <- which(s == 1)
    selected_full <- sel_idx[selected_sel]
    if (m > 0) {
      selected_full <- sort(c(mand_idx, selected_full))
    }
    
    selected_models[[k]] <- selected_full
    best_model_matrix[k, ] <- s
    convergence_iters[k] <- iter
  }
  
  # --- Step 8: Refit on last vertex and evaluate on test set ---
  # Use ORIGINAL (unnormalised) X for refitting (Section 3.1, [1][4])
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
  
  # --- Return results ---
  result <- list(
    selected_models = selected_models,       # list of selected indices for each k
    best_model_matrix = best_model_matrix,   # p_sel-length binary vectors
    best_k = best_k,                         # best k by test accuracy
    best_subset = selected_models[[best_k]], # best selected variables
    test_accuracy = test_accuracy,           # test accuracy for each k
    convergence = convergence_iters          # iterations to converge for each k
  )
  return(result)
}

