# COMBSS-GLM: Best Subset Selection for Generalised Linear Models

R implementation of **COMBSS** (Continuous Optimisation for Best Subset Selection)
for generalised linear models, as described in:

> **Efficient Best Subset Selection in Sparse Generalised Linear Models for Biomedical
> Applications**
>
> Authors: [Your names here]
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


