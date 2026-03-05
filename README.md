# COMBSS-GLM: Best Subset Selection for Generalised Linear Models

R implementation of **COMBSS** (Continuous Optimisation for Best Subset Selection)
for generalised linear models, as described in:

> **Efficient Best Subset Selection in Sparse Generalised Linear Models for Biomedical
> Applications**
>
> Authors: Benoit Liquet, Anant Mathur, Sarat Moka, Samuel Muller
>
> Authors are listed in alphabetical order
>
> Submitted to *Biometrics*, 2026.

## Overview

Best subset selection, finding the optimal subset of *k* predictors that maximises
the likelihood, is fundamental for interpretable and parsimonious statistical
modelling, but is NP-hard in general. COMBSS overcomes this computational barrier
by reformulating the discrete combinatorial problem as a continuous optimisation
over the hypercube $[0,1]^p$ via a Boolean relaxation, and solving it using a
Frank–Wolfe homotopy algorithm.

The key innovation of this implementation is the use of **Danskin's envelope
theorem** to compute the gradient of the relaxed objective exactly, requiring
only a single ridge-penalised GLM solve (via `glmnet`) per iteration — no
Hessian assembly or linear system solve is needed. This makes the method
scalable to high-dimensional settings with p >> n.

For full methodological details, including the Boolean relaxation, the homotopy
schedule, convergence properties, and the multinomial extension, please refer to
the paper.

## Features

- **Binary logistic regression** (two-class classification)
- **Multinomial logistic regression** (multi-class, C > 2 classes)
- **Support for mandatory variables** (always included in the model)
- **Column normalisation** for numerical stability
- **Efficient inner solves** via the `glmnet` package
- **Gradient computation** via Danskin's envelope theorem (no Hessian required)
- **Frank–Wolfe homotopy algorithm** with geometric delta schedule

## Installation

No installation is required. Clone the repository and source the main function
file in R:

```{r}
source("R/COMBSS_logistic.R")
```

## Dependencies

The only dependency is glmnet:

```{r}
install.packages("glmnet")
```

## Quick Start

### Example 1: Low-Dimensional Simulation (n=200, p=30)
This example demonstrates COMBSS in a low-dimensional logistic regression
setting with n=200 observations and p=30 predictors, where the first 10
predictors are truly active.


```{r}
library(glmnet)
source("R/COMBSS_logistic.R")

# Load example data (tab-separated, no header)
# Column 1 = response (0/1), Column 2 = intercept, Columns 3+ = features
train <- read.csv("DATA/Example_small_data_train.csv", sep = "\t", header = FALSE)
test  <- read.csv("DATA/Example_small_data_test.csv",  sep = "\t", header = FALSE)

# Run COMBSS for subset sizes k = 1 to 15
result <- COMBSS_logistic(
  X         = as.matrix(train[, 3:ncol(train)]),
  y         = as.numeric(train[, 1]),
  Kmax      = 15,
  delta_min = 0.1,
  r         = 1.5,
  Niter     = 500,
  alpha     = 0.01,
  lambda    = 0,
  epsilon   = 0.01,
  Xtest     = as.matrix(test[, 3:ncol(test)]),
  ytest     = as.numeric(test[, 1])
)

# View results
cat("Best k:", result$best_k, "\n")
cat("Selected variables:", result$best_subset, "\n")
cat("Test accuracy:", result$test_accuracy[result$best_k], "\n")

```


```
Best k: 10
Best subset: 1 2 3 4 5 6 7 8 9 10

Test accuracies:
  k =  1 | accuracy = 0.704 | variables: 4
  k =  2 | accuracy = 0.756 | variables: 4,6
  k =  3 | accuracy = 0.776 | variables: 3,4,6
  k =  4 | accuracy = 0.829 | variables: 3,4,6,8
  k =  5 | accuracy = 0.837 | variables: 3,4,6,8,9
  k =  6 | accuracy = 0.843 | variables: 3,4,6,7,8,9
  k =  7 | accuracy = 0.850 | variables: 3,4,5,6,7,8,9
  k =  8 | accuracy = 0.878 | variables: 2,3,4,5,6,7,8,9
  k =  9 | accuracy = 0.889 | variables: 1,2,3,4,5,6,7,8,9
  k = 10 | accuracy = 0.898 | variables: 1,2,3,4,5,6,7,8,9,10
  k = 11 | accuracy = 0.888 | variables: 1,2,3,4,5,6,7,8,9,10,11
  k = 12 | accuracy = 0.888 | variables: 1,2,3,4,5,6,7,8,9,10,11,12
  k = 13 | accuracy = 0.888 | variables: 1,2,3,4,5,6,7,8,9,10,11,12,29
  k = 14 | accuracy = 0.888 | variables: 1,2,3,4,5,6,7,8,9,10,11,12,28,29
  k = 15 | accuracy = 0.884 | variables: 1,2,3,4,5,6,7,8,9,10,11,12,28,29,30
```

COMBSS correctly identifies the best model at **k=10**, recovering all 10
truly active variables (variables 1–10) with a test accuracy of **89.8%**.
The inclusion path shows that the strongest signals (variables 4 and 6) are
selected first, with the remaining active variables progressively added as k
increases. Beyond k=10, accuracy slightly decreases as noise variables enter
the model, confirming that the method correctly identifies the true model
size.
