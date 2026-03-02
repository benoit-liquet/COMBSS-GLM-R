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
train <- read.csv("../data/n-200-p1000Replica1.csv", sep = "\t", header = FALSE)
test <- read.table("../data/n-200-p1000Test-Replica1-500.txt", sep = "\t", header = FALSE)
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

cat("Best k:", result$best_k, "\n")
cat("Best subset:", result$best_subset, "\n")
cat("Test accuracies:\n")
for (k in 1:15) {
  cat(sprintf("  k = %2d | accuracy = %.3f | variables: %s\n",
              k, result$test_accuracy[k],
              paste(result$selected_models[[k]], collapse = ",")))
}
