# lsReg

`lsReg` is an R package for large-scale regression. It efficiently tests a large
number of candidate covariates against a fixed base model without re-fitting the
base model each time. The canonical use case is genome-wide association studies
(GWAS), where thousands or millions of genetic variants are screened one at a
time against a common set of baseline covariates.

## Supported models and tests

| Family | `glm` family argument |
|---|---|
| Linear | `gaussian` |
| Logistic | `binomial` |
| Poisson | `poisson` |

| Test | `testtype` |
|---|---|
| Likelihood ratio | `"lrt"` |
| Score (Rao) | `"score"` |
| Wald | `"wald"` |
| Robust score (Huber-White) | `"robustscore"` |
| Robust Wald (Huber-White) | `"robustwald"` |

## How it works

The workflow has two steps:

1. **Allocate** — fit the base model with `glm()`, then call `lsReg()` to
   pre-compute and cache expensive quantities (QR decomposition, residuals,
   information matrices) from the base model.
2. **Test** — call `addcovar()` once per candidate covariate. The cached
   quantities are reused for every test, avoiding redundant computation.

## Installation

```r
# Install from GitHub
remotes::install_github("USCbiostats/lsReg")
```

## Example: large-scale linear regression

The example below screens 10 candidate covariates (`z1`–`z10`) against a
Gaussian base model using the Wald test. The outcome `ylin` was simulated as a
linear function of `x1`, `x2`, `z5`, and `z9`.

```r
library(lsReg)

# Load example data
datafile <- system.file("extdata", "simulated_data.rds", package = "lsReg")
dat <- readRDS(datafile)

# Step 1: fit the base model
basemdl <- glm(ylin ~ x1 + x2, data = dat, family = gaussian)

# Step 2: allocate memory (done once)
mem <- lsReg(basemdl, colstoadd = 1, testtype = "wald")

# Step 3: test each candidate covariate
zvars <- paste0("z", 1:10)

results <- data.frame(
  variable  = zvars,
  estimate  = NA_real_,
  statistic = NA_real_
)

for (i in seq_along(zvars)) {
  xr <- as.matrix(dat[, zvars[i], drop = FALSE])
  addcovar(mem, xr)
  results$estimate[i]  <- mem$fitdata$betab[1]
  results$statistic[i] <- mem$testvalue[1, 1]
}

results$pvalue <- 2 * pnorm(-abs(results$statistic))
print(results, digits = 4)
#>    variable  estimate statistic    pvalue
#> 1        z1 -0.469022 -1.405562 0.1598543
#> 2        z2 -0.327854 -1.875926 0.0606655
#> 3        z3 -0.096312 -0.532620 0.5942964
#> 4        z4 -0.096532 -0.409402 0.6822450
#> 5        z5  0.744905  3.565794 0.0003628
#> 6        z6 -0.024623 -0.082187 0.9344980
#> 7        z7  0.682225  2.241809 0.0249737
#> 8        z8  0.001659  0.005704 0.9954488
#> 9        z9 -0.748608 -2.471511 0.0134543
#> 10      z10 -0.533338 -1.774089 0.0760484
```

`z5` and `z9` show the strongest signals, consistent with the data-generating model.

## Vignettes

- [Large-Scale Linear Regression](vignettes/linear-regression.Rmd)
- [Large-Scale Logistic Regression](vignettes/logistic-regression.Rmd)
- [Large-Scale Poisson Regression](vignettes/poisson-regression.Rmd)
- [Hypothesis Tests in lsReg](vignettes/test-statistics.Rmd)
