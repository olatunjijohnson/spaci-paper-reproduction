## =====================================================================
## Redesigned simulation study (task #8) -- headline pieces:
##  (A) COVERAGE: iid vs spatial-HAC vs block-bootstrap intervals for recoverU+
##      -- closes out T4. Reported for the consistent regime (delta_u=0, coverage
##      of the true ATT is meaningful) and a confounded regime (coverage of the
##      estimator's probability limit, isolating SE calibration from bias).
##  (B) WEIGHTS: mean iDAPS weights (pi1,pi2,pi3) over (u, gamma) -- should track
##      the true confounding/interference mix (pi3 up in gamma, pi2 up in u).
## =====================================================================
suppressMessages(library(spaci))
OUT <- "results"
dir.create(OUT, showWarnings = FALSE)
library(parallel); RNGkind("L'Ecuyer-CMRG")

cover_rep <- function(n, seed, delta_u, theta_spatial = 1.5, gamma = 1.5, tau = 0.1) {
  sim <- simulate_spatial_causal(n = n, delta_u = delta_u, theta_spatial = theta_spatial,
             gamma_interference = gamma, tau_exp = tau, seed = seed)
  fit <- tryCatch(recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, tau = tau),
                  error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$att)) return(NULL)
  hac <- tryCatch(vcov_hac(fit), error = function(e) NULL)
  bb  <- tryCatch(boot_spatial(sim$Y, sim$Z, sim$X, sim$coords,
                    method = "recoverUplus", B = 400, seed = seed), error = function(e) NULL)
  if (is.null(hac) || is.null(bb)) return(NULL)
  c(att = fit$att,
    iid_lo = fit$ci[["lower"]],  iid_hi = fit$ci[["upper"]],
    hac_lo = hac$ci[["lower"]],  hac_hi = hac$ci[["upper"]],
    boot_lo = bb$ci[["lower"]],  boot_hi = bb$ci[["upper"]])
}

weight_rep <- function(n, seed, u, gamma, tau = 0.1) {
  sim <- simulate_spatial_causal(n = n, delta_u = u, theta_spatial = 1.5,
             gamma_interference = gamma, tau_exp = tau, seed = seed)
  fit <- tryCatch(idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = tau, seed = seed),
                  error = function(e) NULL)
  if (is.null(fit) || is.null(fit$weights)) return(NULL)
  w <- unname(fit$weights)
  c(u = u, gamma = gamma, pi1 = w[1], pi2 = w[2], pi3 = w[3])
}

args <- commandArgs(trailingOnly = TRUE); MODE <- if (length(args) >= 1) args[1] else "coverage"

if (MODE == "coverage") {
  cells <- list(list(du = 0, n = 150), list(du = 0, n = 300), list(du = 2, n = 150))
  R <- 400
  res <- lapply(seq_along(cells), function(k) {
    cc <- cells[[k]]
    M <- do.call(rbind, mclapply(1:R, function(r)
      tryCatch(cover_rep(cc$n, seed = 40000 * k + r, delta_u = cc$du),
               error = function(e) NULL), mc.cores = 11))
    list(delta_u = cc$du, n = cc$n, M = M)
  })
  saveRDS(res, file.path(OUT, "coverage.rds")); cat("COVERAGE_DONE\n")
}

if (MODE == "coverage2") {
  ## fill the confounding axis at n=150 (the regime where recoverU+ is used)
  cells <- list(list(du = 0.5, n = 150), list(du = 1.0, n = 150))
  R <- 300
  res <- lapply(seq_along(cells), function(k) {
    cc <- cells[[k]]
    M <- do.call(rbind, mclapply(1:R, function(r)
      tryCatch(cover_rep(cc$n, seed = 50000 * k + r, delta_u = cc$du),
               error = function(e) NULL), mc.cores = 11))
    list(delta_u = cc$du, n = cc$n, M = M)
  })
  saveRDS(res, file.path(OUT, "coverage2.rds")); cat("COVERAGE2_DONE\n")
}

if (MODE == "weights") {
  cells <- expand.grid(u = c(0, 1, 2), gamma = c(0, 1.5, 3))
  R <- 160; n <- 250
  res <- lapply(seq_len(nrow(cells)), function(k) {
    u <- cells$u[k]; g <- cells$gamma[k]
    M <- do.call(rbind, mclapply(1:R, function(r)
      tryCatch(weight_rep(n, seed = 60000 * k + r, u = u, gamma = g),
               error = function(e) NULL), mc.cores = 11))
    list(u = u, gamma = g, M = M)
  })
  saveRDS(res, file.path(OUT, "weights.rds")); cat("WEIGHTS_DONE\n")
}
