# COMBSS for Generalised Linear Models

R implementation of **COMBSS** (Continuous Optimisation for Best Subset Selection)
for logistic and multinomial regression, as described in:

> *Efficient Best Subset Selection in Sparse Generalized Linear Models
> for Biomedical Applications* (2026)

COMBSS reformulates the NP-hard best subset selection problem as a continuous
optimisation over the hypercube [0,1]^p via a Boolean relaxation, and solves it
using a Frank–Wolfe homotopy algorithm. The gradient of the relaxed objective
is computed exactly using Danskin's envelope theorem, requiring only a single
ridge-penalised GLM solve per iteration.

## Features

- Binary logistic regression (two-class)
- Multinomial logistic regression (multi-class, via Appendix A of the paper)
- Support for mandatory variables (always included in the model)
- Column normalisation for numerical stability
- Efficient inner solves via `glmnet`

## Installation

No installation required. Clone the repository and source the main function file:

```r
source("R/COMBSS_logistic.R")
```
