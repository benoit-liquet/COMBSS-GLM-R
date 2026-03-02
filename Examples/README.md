# Example of COMBSS-glm for logistic model
## High dimensional data set with p=10000 and n=200

The first 10 features are related to the outcome (see the paper regarding the true generated model)

## Source the model 
```{r}
library(glmnet)
source("../R/COMBSS_logistics.R")
```

## Load example data (tab-separated, no header)
### Column 1 = response (0/1), Column 2 = intercept, Columns 3+ = features
```{r}
train <- read.csv("../DATA/n-200-p1000Replica1.csv", sep = "\t", header = FALSE)
test <- read.table("../DATA/n-200-p1000Test-Replica1-500.txt", sep = "\t", header = FALSE)
dim(test)
```{r}

The test data here is only 1000 data points for available github storage requirement.


## Run COMBSS for subset sizes k = 1 to 15
```{r}
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
```

## View results

```{r}
cat("Best k:", result$best_k, "\n")
cat("Best subset:", result$best_subset, "\n")
cat("Test accuracies:\n")
for (k in 1:15) {
  cat(sprintf("  k = %2d | accuracy = %.3f | variables: %s\n",
              k, result$test_accuracy[k],
              paste(result$selected_models[[k]], collapse = ",")))
}
```



```

Best k: 12
Best subset: 1 2 3 4 5 6 7 8 9 10 92 760

Test accuracies:
  k =  1 | accuracy = 0.707 | variables: 7
  k =  2 | accuracy = 0.786 | variables: 4,7
  k =  3 | accuracy = 0.791 | variables: 4,6,7
  k =  4 | accuracy = 0.789 | variables: 4,6,7,8
  k =  5 | accuracy = 0.834 | variables: 2,4,6,7,8
  k =  6 | accuracy = 0.837 | variables: 2,4,5,6,7,8
  k =  7 | accuracy = 0.842 | variables: 2,3,4,5,6,7,8
  k =  8 | accuracy = 0.858 | variables: 1,2,3,4,5,6,7,8
  k =  9 | accuracy = 0.872 | variables: 1,2,3,4,5,6,7,8,9
  k = 10 | accuracy = 0.879 | variables: 1,2,3,4,5,6,7,8,9,10
  k = 11 | accuracy = 0.880 | variables: 1,2,3,4,5,6,7,8,9,10,92
  k = 12 | accuracy = 0.882 | variables: 1,2,3,4,5,6,7,8,9,10,92,760
  k = 13 | accuracy = 0.879 | variables: 1,2,3,4,5,6,7,8,9,10,92,486,760
  k = 14 | accuracy = 0.880 | variables: 1,2,3,4,5,6,7,8,9,10,92,486,760,825
  k = 15 | accuracy = 0.880 | variables: 1,2,3,4,5,6,7,8,9,10,92,486,760,825,978
```

COMBSS successfully recovers all 10 truly active variables (variables 1–10)
by **k=10** with a test accuracy of **87.9%**, despite the challenging $p >> n$
setting (p=1000, n=200)
