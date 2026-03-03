# COMBSS for Multinomial Regression: Khan SRBCT Application



## The Khan SRBCT Dataset

The small round blue cell tumour (SRBCT) dataset of Khan et al. (2001)
contains expression levels of **p = 2,308 genes** measured on tissue samples
belonging to four distinct childhood tumour types: Ewing sarcoma (EWS),
Burkitt lymphoma (BL), neuroblastoma (NB), and rhabdomyosarcoma (RMS).

The dataset is split into a training set of **n_train = 63 samples** and an
independent test set of **n_test = 20 samples**, as defined in the original
study. This partition is preserved in the R package `ISLR2`. The four-class
structure (C = 4) and extreme dimensionality ratio (p/n ≈ 37) make this
dataset a challenging and widely studied benchmark for high-dimensional
classification and variable selection.

## Repository Structure

```
COMBSS-Multinomial/
├── README.md
├── R/
│   └── Benoit_function_singlek.R     # COMBSS multinomial functions
├── examples/
│   └── example_khan.R                # Khan SRBCT demo script
└── results/
    └── figures/                      # Generated figures
```



### Dependencies

```r
install.packages(c("glmnet", "nnet", "ISLR2"))
```

## Quick Start

```r
library(ISLR2)
library(glmnet)
library(nnet)

# Source the COMBSS multinomial code
source("R/Benoit_function_singlek.R")

# Load Khan SRBCT data
data(Khan)
X_train <- Khan$xtrain
y_train <- Khan$ytrain
X_test  <- Khan$xtest
y_test  <- Khan$ytest

cat(sprintf("Training: %d samples, %d genes, %d classes\n",
            nrow(X_train), ncol(X_train), length(unique(y_train))))
cat(sprintf("Test:     %d samples\n\n", nrow(X_test)))
```

**Output:**

```
Training: 63 samples, 2308 genes, 4 classes
Test:     20 samples
```

## Running COMBSS for a Single Subset Size

To run COMBSS for a specific subset size (e.g., k = 12), use the
`COMBSS_multinomial_single_k()` function. The penalty schedule
(δ_min, δ_max, r) is automatically calibrated from the data using
power iteration to estimate the largest eigenvalue of $X_u^T X_u$.

```r
result_k12 <- COMBSS_multinomial_single_k(
  X            = X_train,
  y            = y_train,
  k            = 12,
  Niter        = 50,
  lambda       = 0,
  Xtest        = X_test,
  ytest        = y_test,
  normalize    = TRUE,
  lambda_refit = 0,
  refit_method = "glmnet",
  verbose      = TRUE
)
```

**Output:**

```
Calibration: nu_max=1228.3661 | delta_conc=4.874469 | delta_min=0.004874 | delta_max=4.874469 | r=1.148154

--- k = 12 | lambda = 0.000000 ---
  iter  10 | delta = 0.019406 | alpha = 0.2000
  iter  20 | delta = 0.077255 | alpha = 0.4000
  iter  30 | delta = 0.307558 | alpha = 0.6000
  iter  40 | delta = 1.224411 | alpha = 0.8000
  iter  50 | delta = 4.874469 | alpha = 1.0000
  Selected: 187, 246, 509, 545, 910, 1074, 1319, 1389, 1645, 1954, 1955, 2050
  Test accuracy: 1.000
```

COMBSS selects **12 genes** out of 2,308 and achieves **perfect classification
(100%)** on the independent test set of 20 samples.

## Running COMBSS for All Subset Sizes k = 1 to 20

To compute the full inclusion path, use the `COMBSS_multinomial()` function
which loops over all subset sizes from k = 1 to Kmax:

```r
result <- COMBSS_multinomial(
  X            = X_train,
  y            = y_train,
  Kmax         = 20,
  Niter        = 50,
  lambda       = 0,
  Xtest        = X_test,
  ytest        = y_test,
  refit_method = "glmnet",
  lambda_refit = 0,
  verbose      = TRUE
)
```

**Output:**

```
k =  1 | acc = 0.45 | genes: 1954
k =  2 | acc = 0.70 | genes: 1954, 1955
k =  3 | acc = 0.75 | genes: 246, 1954, 1955
k =  4 | acc = 0.70 | genes: 246, 1389, 1954, 1955
k =  5 | acc = 0.85 | genes: 246, 1389, 1645, 1954, 1955
k =  6 | acc = 0.85 | genes: 187, 246, 1389, 1645, 1954, 1955
k =  7 | acc = 0.85 | genes: 187, 246, 509, 1389, 1645, 1954, 1955
k =  8 | acc = 0.95 | genes: 187, 246, 509, 1389, 1645, 1954, 1955, 2050
k =  9 | acc = 0.95 | genes: 187, 246, 509, 545, 1389, 1645, 1954, 1955, 2050
k = 10 | acc = 0.95 | genes: 187, 246, 509, 545, 1319, 1389, 1645, 1954, 1955, 2050
k = 11 | acc = 0.95 | genes: 187, 246, 509, 545, 910, 1319, 1389, 1645, 1954, 1955, 2050
k = 12 | acc = 1.00 | genes: 187, 246, 509, 545, 910, 1074, 1319, 1389, 1645, 1954, 1955, 2050
k = 13 | acc = 1.00 | genes: 187, 246, 509, 545, 910, 1074, 1319, 1389, 1645, 1954, 1955, 1980, 2050
k = 14 | acc = 0.95 | genes: 187, 246, 509, 545, 910, 1074, 1319, 1389, 1645, 1772, 1954, 1955, 1980, 2050
k = 15 | acc = 0.95 | genes: 107, 187, 246, 509, 545, 910, 1074, 1319, 1389, 1645, 1772, 1954, 1955, 1980, 2050
k = 16 | acc = 0.95 | genes: 107, 187, 246, 509, 545, 910, 1074, 1319, 1389, 1645, 1708, 1772, 1954, 1955, 1980, 2050
k = 17 | acc = 0.95 | genes: 107, 187, 246, 373, 509, 545, 910, 1074, 1319, 1389, 1645, 1708, 1772, 1954, 1955, 1980, 2050
k = 18 | acc = 0.95 | genes: 107, 187, 246, 373, 509, 545, 566, 910, 1074, 1319, 1389, 1645, 1708, 1772, 1954, 1955, 1980, 2050
k = 19 | acc = 0.95 | genes: 107, 187, 246, 373, 509, 545, 566, 910, 1074, 1207, 1319, 1389, 1645, 1708, 1772, 1954, 1955, 1980, 2050
k = 20 | acc = 0.95 | genes: 107, 187, 246, 373, 509, 545, 566, 910, 1074, 1093, 1207, 1319, 1389, 1645, 1708, 1772, 1954, 1955, 1980, 2050
```

## Results Summary

### Classification Performance on the Khan SRBCT Test Set

| Method | Genes selected | Test accuracy |
|---|:---:|:---:|
| COMBSS (k = 8) | 8 | 0.95 |
| COMBSS (k = 12) | 12 | **1.00** |
| COMBSS (k = 13) | 13 | **1.00** |
| Group Lasso (λ.1se) | 28 | 0.95 |
| Group Lasso (λ.min) | 30 | 0.95 |
| Group Lasso (best λ) | 35 | 1.00 |

COMBSS achieves perfect classification accuracy (20/20) on the independent
test set using only **12 genes** out of 2,308. By contrast, the multinomial
group Lasso requires **35 genes** to reach the same performance, and its
standard cross-validated models (λ.min and λ.1se) select 30 and 28 genes
respectively while achieving only 95% test accuracy. 

### Best-Subset Inclusion Path

The inclusion path reveals that genes 1954 and 1955 are the first to enter
the model (at k = 1 and k = 2), followed by gene 246 at k = 3. A core set of
genes is progressively built up: by k = 8 the model already achieves 95%
accuracy, and perfect classification is reached at k = 12 with genes 187,
246, 509, 545, 910, 1074, 1319, 1389, 1645, 1954, 1955, and 2050. Beyond
k = 13, accuracy drops back to 95% and remains stable.

## Algorithm Details

### Penalty Schedule Calibration

The penalty schedule is automatically calibrated from the data. The largest
eigenvalue ν_max of X_u^T X_u (where X_u is the column-normalised design
matrix of relaxed variables) is estimated via power iteration. The
concentration penalty is then set as δ_conc = ν_max / (4n) for multinomial
models, and the schedule parameters are:

```
δ_min = 10⁻³ × δ_conc
δ_max = δ_conc
r = (δ_max / δ_min)^(1/N)
```

This ensures the penalty traverses the full range [δ_min, δ_conc] uniformly
on a log-scale within N steps.

### Inner Solver

The inner ridge-penalised multinomial GLM is solved via `glmnet` with
`family = "multinomial"` and `type.multinomial = "grouped"`, which enforces
a common sparsity pattern across the C−1 = 3 class-specific coefficient
vectors. The per-variable penalty weights ω_j(t) = (λ + δ)/t_j² − δ are
mapped to `glmnet`'s `penalty.factor` parameterisation as:

```
pf_j = p × ω_j / Σ_j ω_j
λ_glmnet = 2 × Σ_j ω_j / p
```

### Danskin Gradient (Multinomial Extension)

The gradient of the relaxed objective with respect to t_j is computed via
Danskin's envelope theorem. For the multinomial case, this sums the squared
coefficients across all C−1 classes for each variable j:

```
∂f/∂t_j = −2(λ + δ) ‖Ξ_{m+j,:}‖² / t_j³
```

where Ξ is the minimiser of the inner ridge problem. This avoids any Hessian
computation and requires only one call to `glmnet` per iteration.

### Refitting

After the algorithm selects a support set for each k, the final model is
refit on the original (unnormalised) training data. The refitting can be
done via `glmnet` (ridge with a small or zero penalty) or `nnet::multinom`
(with sufficient iterations and a small decay for numerical stability in
near-separable settings).

## Function Reference

### `COMBSS_multinomial_single_k()`

Runs the Frank–Wolfe homotopy algorithm for a **single** subset size k.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `X` | matrix | — | Design matrix (n × p), without intercept |
| `y` | factor/numeric | — | Response vector (class labels 1, ..., C) |
| `k` | integer | — | Target subset size |
| `Niter` | integer | 50 | Number of Frank–Wolfe iterations |
| `lambda` | numeric | 0 | Ridge parameter λ ≥ 0 |
| `delta_min` | numeric | NULL | Initial δ (auto-calibrated if NULL) |
| `delta_max` | numeric | NULL | Maximum δ (auto-calibrated if NULL) |
| `mandatory` | integer | NULL | Indices of mandatory variables |
| `Xtest` | matrix | NULL | Test design matrix |
| `ytest` | factor/numeric | NULL | Test response vector |
| `normalize` | logical | TRUE | Column-normalise X to unit ℓ₂ norm |
| `lambda_refit` | numeric | 0 | Ridge penalty for refitting (0 = near MLE) |
| `refit_method` | character | "glmnet" | Refitting method: `"nnet"` or `"glmnet"` |

### `COMBSS_multinomial()`

Runs the algorithm for **all** subset sizes k = 1, ..., Kmax sequentially.
Same parameters as above, with `Kmax` replacing `k`.


