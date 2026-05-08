# tests/benchmark.R
# Compare the master ("base") and quic-improvements ("improved") builds of mSSL
# on identical synthetic data. Times the two solvers and reports differences
# in the recovered B and Omega.
#
# Invoked by tests/run_benchmark.sh — but you can also run it manually:
#   Rscript tests/benchmark.R tests/lib_base tests/lib_improved tests/results

suppressMessages({
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript benchmark.R <lib_base> <lib_improved> <results_dir>")
}
LIB_BASE    <- normalizePath(args[1])
LIB_NEW     <- normalizePath(args[2])
RESULTS_DIR <- normalizePath(args[3], mustWork = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Build a reproducible synthetic dataset -------------------------------
# Match the scale of examples/example.R but a touch smaller so a single run
# takes seconds, not minutes. Bump (n, p, q) up if you want a heavier test.
make_data <- function(seed = 129, n = 200, p = 200, q = 15, rho = 0.7,
                      sparsity = 0.2) {
  set.seed(seed)
  SigmaX <- 0.7 ^ (abs(outer(1:p, 1:p, "-")))
  X <- mvrnorm(n, rep(0, p), SigmaX)

  nz <- floor(p * q * sparsity)
  B  <- matrix(sample(c(runif(nz, -2, 2), rep(0, p * q - nz))), p, q)

  Sigma <- rho ^ (abs(outer(1:q, 1:q, "-")))
  Omega <- solve(Sigma); Omega[abs(Omega) < 1e-3] <- 0

  set.seed(seed + 1L)
  E <- mvrnorm(n, rep(0, q), Sigma)
  list(X = X, Y = X %*% B + E, B_true = B, Omega_true = Omega)
}

dat <- make_data()
cat(sprintf("Data: n=%d, p=%d, q=%d\n",
            nrow(dat$X), ncol(dat$X), ncol(dat$Y)))

# ---- 2. Run one fitter from a specified library tree -------------------------
# We re-launch a fresh R process for each version so the two builds of mSSL
# never share a namespace. That avoids any DLL-collision oddness.
run_fit <- function(libpath, label, dat, fitter = "mSSL_dcpe") {
  rds_in  <- tempfile(fileext = ".rds")
  rds_out <- tempfile(fileext = ".rds")
  saveRDS(dat, rds_in)

  script <- sprintf('
    .libPaths(c("%s", .libPaths()))
    suppressMessages(library(mSSL))
    dat <- readRDS("%s")
    gc(reset = TRUE, full = TRUE)
    t0  <- proc.time()
    fit <- %s(dat$X, dat$Y)
    elapsed <- (proc.time() - t0)[["elapsed"]]
    saveRDS(list(fit = fit, elapsed = elapsed,
                 sessionInfo = utils::sessionInfo()), "%s")
  ', libpath, rds_in, fitter, rds_out)

  script_file <- tempfile(fileext = ".R")
  writeLines(script, script_file)

  status <- system2("Rscript", c("--vanilla", script_file), stdout = TRUE,
                    stderr = TRUE)
  if (!file.exists(rds_out)) {
    cat(paste(status, collapse = "\n"), "\n")
    stop("Fit failed for ", label)
  }
  readRDS(rds_out)
}

# ---- 3. Run both versions ----------------------------------------------------
cat("\n==> Running BASE (master)\n")
res_base <- run_fit(LIB_BASE, "base", dat)
cat(sprintf("    elapsed: %.2fs\n", res_base$elapsed))

cat("\n==> Running IMPROVED (quic-improvements)\n")
res_new  <- run_fit(LIB_NEW, "improved", dat)
cat(sprintf("    elapsed: %.2fs\n", res_new$elapsed))

# ---- 4. Compare outputs ------------------------------------------------------
# mSSL_dcpe returns (at least) B and Omega — pull those and diff.
B_base  <- res_base$fit$B;  O_base <- res_base$fit$Omega
B_new   <- res_new$fit$B;   O_new  <- res_new$fit$Omega

frob <- function(A) sqrt(sum(A * A))
rel_diff <- function(A, B) frob(A - B) / max(frob(A), 1e-12)

# Sparsity-pattern agreement
support_match <- function(A, B, tol = 1e-8) {
  sA <- abs(A) > tol; sB <- abs(B) > tol
  list(jaccard = sum(sA & sB) / max(sum(sA | sB), 1L),
       n_diff  = sum(sA != sB))
}

sup_B <- support_match(B_base, B_new)
sup_O <- support_match(O_base, O_new)

# ---- 5. Write report ---------------------------------------------------------
report <- c(
  "mSSL benchmark — base vs. improved",
  paste(rep("=", 60), collapse = ""),
  sprintf("Data: n=%d, p=%d, q=%d", nrow(dat$X), ncol(dat$X), ncol(dat$Y)),
  "",
  "Timing (seconds, mSSL_dcpe):",
  sprintf("  base     : %.3f", res_base$elapsed),
  sprintf("  improved : %.3f", res_new$elapsed),
  sprintf("  speedup  : %.2fx",
          res_base$elapsed / max(res_new$elapsed, 1e-9)),
  "",
  "Output agreement (improved vs. base):",
  sprintf("  ||B_new - B_base||_F / ||B_base||_F     = %.3e",
          rel_diff(B_base, B_new)),
  sprintf("  ||O_new - O_base||_F / ||O_base||_F     = %.3e",
          rel_diff(O_base, O_new)),
  sprintf("  B support Jaccard = %.4f   (mismatched entries: %d)",
          sup_B$jaccard, sup_B$n_diff),
  sprintf("  O support Jaccard = %.4f   (mismatched entries: %d)",
          sup_O$jaccard, sup_O$n_diff),
  "",
  "Recovery vs. ground truth (lower is better):",
  sprintf("  base     B err = %.3e   O err = %.3e",
          rel_diff(dat$B_true,     B_base),
          rel_diff(dat$Omega_true, O_base)),
  sprintf("  improved B err = %.3e   O err = %.3e",
          rel_diff(dat$B_true,     B_new),
          rel_diff(dat$Omega_true, O_new))
)

out_path <- file.path(RESULTS_DIR, "comparison.txt")
writeLines(report, out_path)
saveRDS(list(base = res_base, improved = res_new, data = dat),
        file.path(RESULTS_DIR, "raw_results.rds"))

cat("\n", paste(report, collapse = "\n"), "\n", sep = "")
cat(sprintf("\nWrote %s\n", out_path))
