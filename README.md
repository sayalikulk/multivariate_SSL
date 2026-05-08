## mSSL: The Multivariate Spike-and-Slab LASSO

[Sameer K. Deshpande](https://people.csail.mit.edu/sameerd/)

An R package for simultaneous variable and covariance selection in multivariate linear regression with spike-and-slab LASSO (SSL) regularization.
For further details about the method, please see [Deshpande, Rockova, and George (2018)](https://arxiv.org/abs/1708.08911).


### Details

#### Installation

The package source files are contained in the sub-directory mSSL/.
To install, you can either download that directory and then build and install the package from the command line (e.g. `R CMD BUILD mSSL` followed by `R CMD INSTALL mSSL_1.0.tar.gz`).
You can also install the package using `devtools::install_github` as follows.

```r
library(devtools)
devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")
```

#### Examples

The following code chunk shows the basic usage of the basic functions.
```r
library(SSLASSO)
library(MASS)
library(mSSL)

n <- 100
p <- 50
q <- 25
rho <- 0.9

set.seed(129)
SigmaX <- 0.7^(abs(outer(1:p, 1:p, "-")))
X <- mvrnorm(n, rep(0, times = p), SigmaX)
B_orig <- matrix(sample(c(runif(floor(p*q*0.2), -2,2),
                          rep(0, times = p*q - floor(p*q*0.2)))), nrow = p, ncol = q)
XB <- X %*% B_orig
Sigma_orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega_orig <- solve(Sigma_orig)
Omega_orig[abs(Omega_orig) < 1e-3] <- 0 # need to truncate the essentially negligible terms introduced by solve
set.seed(130)
E <- mvrnorm(n, rep(0, times = q), Sigma_orig)
Y <- XB + E

fit_mSSL_dpe <- mSSL_dpe(X,Y)
fit_mSSL_dcpe <- mSSL_dcpe(X,Y)
fit_sslg <- sep_sslg(X,Y)


# Compare variable selection performance
round(error_B(fit_mSSL_dpe$B, B_orig), digits = 4)
round(error_B(fit_mSSL_dcpe$B, B_orig), digits = 4)
round(error_B(fit_sslg$B, B_orig), digits = 4)

# Compare covariance selection performance
round(error_Omega(fit_mSSL_dpe$Omega, Omega_orig), digits = 4)
round(error_Omega(fit_mSSL_dcpe$Omega, Omega_orig), digits = 4)
round(error_Omega(fit_sslg$Omega, Omega_orig), digits = 4)
```


The sub-directory examples/ contains R scripts to replicate some of the examples in our paper.

#### Benchmarking the QUIC / B_theta improvements

The `tests/` directory contains a script that compares the **pre-improvements** build of mSSL (commit `0b0f157` on `master`) against the current `quic-improvements` branch on identical synthetic data. It checks two things at once: that the improved solver still produces the same answer, and how much faster (if any) it is.

To run it:

```bash
bash tests/run_benchmark.sh
```

What the script does:

1. Creates a git worktree at `.bench-master/` checked out to `master` (the pre-improvements baseline).
2. Builds the **base** package from that worktree into `tests/lib_base/`  (`R CMD INSTALL --library=tests/lib_base mSSL`).
3. Builds the **improved** package from the current branch into `tests/lib_improved/`.
4. Runs `tests/benchmark.R`, which launches a fresh R process per build (using `.libPaths(...)` so the two namespaces never collide), fits `mSSL_dcpe` on the same simulated `(X, Y)` in each, and writes a report.

The report (`tests/results/comparison.txt`) contains:

- **Timing** for each build and the speedup ratio.
- **Output agreement**: relative Frobenius differences `‖B_new − B_base‖_F / ‖B_base‖_F` and the same for `Ω`, plus Jaccard agreement on the sparsity pattern (i.e. whether both versions zero out the same entries).
- **Recovery error vs. ground truth** (`B_true`, `Ω_true`) for each build, so you can confirm that any small numerical drift hasn't degraded statistical accuracy.

Bump `n`, `p`, or `q` in the `make_data()` call near the top of `tests/benchmark.R` for a heavier test.

Requirements: a working R install (`R` and `Rscript` on `PATH`) with `Rcpp`, `RcppArmadillo`, `SSLASSO`, and `MASS`. The build trees in `tests/lib_base/` and `tests/lib_improved/` plus the worktree at `.bench-master/` and `tests/results/` are gitignored — they're regenerated on each run.

