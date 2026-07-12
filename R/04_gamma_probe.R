## =====================================================================
## Gamma orthogonality probe for recoverU+ (paper-development, T5 §4 / T6 §2)
##
## Q: is the DR ATT first-order INSENSITIVE to the recovered-confounder
##    ESTIMATION step? If yes (Gamma = 0), Theorem 1 uses the simple influence
##    function (no generated-regressor correction).
##
## Primary test: compare the real estimator (Matern params phi estimated per
##   dataset) to one with phi held FIXED at its population value phi0 -- same
##   kriging construction & target, differing only by first-stage estimation
##   error. Delta_phi = tau(phi_hat) - tau(phi0). If sd(Delta_phi) shrinks like
##   ~n^-1 (faster than tau's ~n^-1/2), phi-estimation is first-order negligible
##   => Gamma_phi ~ 0. If ~n^-1/2 (constant ratio), Gamma_phi != 0.
##
## Secondary: true-U and U_R-proxy oracles (cost-of-recovery / bias view).
## =====================================================================

suppressMessages(library(spaci))
OUT <- "results"
dir.create(OUT, showWarnings = FALSE)

## ---- stable exponential (nu = 0.5 fixed) Matern MLE on residual field ----
## The full 4-param MLE degenerates (nu -> boundary) on weak residual fields;
## fixing nu = 0.5 matches the exponential DGP and stabilises identification.
fit_matern_exp <- function(resid, coords) {
  d <- as.matrix(stats::dist(coords)); r <- as.numeric(resid); n <- length(r)
  v <- stats::var(r); if (!is.finite(v) || v <= 0) v <- 1
  nll <- function(p) {
    s2 <- exp(p[1]); th <- exp(p[2]); ng <- exp(p[3])
    C <- matern_cov_matrix(d, s2, th, 0.5) + ng * diag(n)
    ch <- tryCatch(chol(C), error = function(e) NULL); if (is.null(ch)) return(1e10)
    one <- rep(1, n)
    Ci_one <- backsolve(ch, forwardsolve(t(ch), one))
    Ci_r   <- backsolve(ch, forwardsolve(t(ch), r))
    mu <- sum(one * Ci_r) / sum(one * Ci_one); res <- r - mu
    Ci_res <- backsolve(ch, forwardsolve(t(ch), res))
    val <- 0.5 * (2 * sum(log(diag(ch))) + sum(res * Ci_res) + n * log(2 * pi))
    if (is.finite(val)) val else 1e10
  }
  st <- log(c(0.7 * v, 0.2, 0.3 * v))
  o <- tryCatch(stats::optim(st, nll, method = "Nelder-Mead",
                             control = list(maxit = 400, reltol = 1e-8)),
                error = function(e) NULL)
  if (is.null(o)) return(list(sigma2 = 0.7 * v, theta = 0.2, nu = 0.5, sigma2_eps = 0.3 * v))
  p <- exp(o$par); list(sigma2 = p[1], theta = p[2], nu = 0.5, sigma2_eps = p[3])
}

## ---- recovery step (mirrors recoverU_core lines 32-63); phi optional ----
recover_UR <- function(Y, Z, X, coords, E, phi = NULL) {
  n <- length(Y); Xdf <- as.data.frame(X); xnames <- colnames(X)
  dat <- data.frame(Y = Y, Z = Z, Xdf, G = E)
  init_form <- stats::as.formula(paste("Y ~ Z +", paste(c(xnames, "G"), collapse = " + ")))
  fit_initial <- stats::lm(init_form, data = dat)
  resid_initial <- stats::residuals(fit_initial)
  d_space_raw <- as.matrix(stats::dist(coords))

  par_hat <- if (is.null(phi)) fit_matern_exp(resid_initial, coords) else phi
  Sigma_U <- matern_cov_matrix(d_space_raw, par_hat$sigma2, par_hat$theta, par_hat$nu)
  V <- Sigma_U + max(par_hat$sigma2_eps, 1e-8) * diag(n) + 1e-6 * diag(n)
  W <- stats::model.matrix(init_form, data = dat)
  theta_gls <- solve(t(W) %*% solve(V, W), t(W) %*% solve(V, Y))
  resid_gls <- as.numeric(Y - W %*% theta_gls)
  Uhat <- as.vector(Sigma_U %*% solve(V, resid_gls))
  list(Uhat = Uhat, phi = par_hat)
}

## ---- faithful recoverU+ downstream with an injectable confounder ----
dr_downstream <- function(Y, Z, X, E, Uconf) {
  n_obs <- length(Y); n1 <- sum(Z == 1)
  Xdf <- as.data.frame(X); xnames <- colnames(X)
  dat <- data.frame(Y = Y, Z = Z, Xdf, G = E, Uhat = safe_scale(Uconf))
  rhs <- c(xnames, "G", "Uhat")
  ps_fit <- suppressWarnings(stats::glm(
    stats::as.formula(paste("Z ~", paste(rhs, collapse = " + "))),
    family = stats::binomial(), data = dat,
    control = stats::glm.control(maxit = 100)))
  ehat <- clip_ps(stats::fitted(ps_fit))
  m0_fit <- stats::lm(stats::as.formula(paste("Y ~", paste(rhs, collapse = " + "))),
                      data = dat[dat$Z == 0, , drop = FALSE])
  m0hat <- as.numeric(stats::predict(m0_fit, newdata = dat))
  w <- ehat / (1 - ehat)
  psi <- (Z - (1 - Z) * w) * (Y - m0hat)
  att <- sum(psi) / n1
  se <- stats::sd(psi / (n1 / n_obs) - att, na.rm = TRUE) / sqrt(n_obs)
  list(att = att, se = se)
}

## ---- one replicate ----
## theta_spatial = 1.5 gives a residual spatial signal that dominates nugget,
## so the recovery is well-identified (else the probe is uninformative).
one_rep <- function(n, seed, phi0 = NULL, tau = 0.1, delta_u = 2, gamma = 1.5,
                    theta_spatial = 1.5) {
  sim <- simulate_spatial_causal(n = n, delta_u = delta_u,
                                 gamma_interference = gamma, tau_exp = tau,
                                 theta_spatial = theta_spatial, seed = seed)
  Y <- sim$Y; Z <- sim$Z; X <- sim$X; S <- sim$coords; U <- sim$U
  E <- neighbourhood_exposure(S, Z, tau = tau, normalize = TRUE)$E

  rec <- tryCatch(recover_UR(Y, Z, X, S, E, phi = NULL), error = function(e) NULL)
  if (is.null(rec)) return(NULL)
  d_rec <- dr_downstream(Y, Z, X, E, rec$Uhat)

  out <- c(att_rec = d_rec$att, se_rec = d_rec$se,
           phi_sigma2 = rec$phi$sigma2, phi_theta = rec$phi$theta,
           phi_nu = rec$phi$nu, phi_eps = rec$phi$sigma2_eps)

  if (!is.null(phi0)) {
    fix <- tryCatch(recover_UR(Y, Z, X, S, E, phi = phi0), error = function(e) NULL)
    if (!is.null(fix)) out["att_fixphi"] <- dr_downstream(Y, Z, X, E, fix$Uhat)$att
    ## secondary oracles
    out["att_oroU"]  <- dr_downstream(Y, Z, X, E, U)$att
    UR2 <- stats::residuals(stats::lm(U ~ Z))
    out["att_oroUR"] <- dr_downstream(Y, Z, X, E, UR2)$att
  }
  out
}

## ---------------- MODE dispatch ----------------
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "validate"

if (MODE == "validate") {
  cat("=== VALIDATION ===\n")
  cat("(a) dr_downstream fed the PACKAGE's recovered Uhat == package att:\n")
  for (s in 1:3) {
    sim <- simulate_spatial_causal(n = 200, delta_u = 2, gamma_interference = 1.5,
                                   theta_spatial = 1.5, tau_exp = 0.1, seed = 100 + s)
    pkg <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, matern_method = "mle")
    E <- neighbourhood_exposure(sim$coords, sim$Z, tau = 0.1)$E
    our <- dr_downstream(sim$Y, sim$Z, sim$X, E, pkg$extras$Uhat)$att
    cat(sprintf("   seed %d: pkg=%.8f our=%.8f |diff|=%.2e\n",
                100 + s, pkg$att, our, abs(pkg$att - our)))
  }
  cat("(b) exponential-nu recovery: does Uhat correlate with true U? phi sane?\n")
  for (s in 1:3) {
    sim <- simulate_spatial_causal(n = 250, delta_u = 2, gamma_interference = 1.5,
                                   theta_spatial = 1.5, tau_exp = 0.1, seed = 200 + s)
    E <- neighbourhood_exposure(sim$coords, sim$Z, tau = 0.1)$E
    rec <- recover_UR(sim$Y, sim$Z, sim$X, sim$coords, E, phi = NULL)
    cat(sprintf("   seed %d: cor(Uhat,U)=%.3f  phi=(s2=%.3f,theta=%.3f,nug=%.3f)\n",
                200 + s, cor(rec$Uhat, sim$U), rec$phi$sigma2, rec$phi$theta, rec$phi$sigma2_eps))
  }
}

if (MODE == "calib") {
  ## estimate population phi0 = mean phi_hat over reps at moderate n
  library(parallel)
  R <- 40; n <- 300
  M <- do.call(rbind, mclapply(1:R, function(r) {
    x <- tryCatch(one_rep(n, seed = 5000 + r), error = function(e) NULL)
    if (is.null(x)) NULL else x[c("phi_sigma2","phi_theta","phi_nu","phi_eps")]
  }, mc.cores = 10))
  phi0 <- list(sigma2 = median(M[,1]), theta = median(M[,2]),
               nu = 0.5, sigma2_eps = median(M[,4]))
  saveRDS(phi0, file.path(OUT, "phi0.rds"))
  cat("phi0 (median phi_hat over", nrow(M), "reps at n=300):\n"); str(phi0)
}

if (MODE == "run") {
  library(parallel)
  RNGkind("L'Ecuyer-CMRG")
  phi0 <- readRDS(file.path(OUT, "phi0.rds"))
  grid <- list(list(n = 100, R = 160), list(n = 150, R = 160),
               list(n = 225, R = 160), list(n = 350, R = 130),
               list(n = 500, R = 90))
  for (g in grid) {
    f <- file.path(OUT, sprintf("res_n%04d.rds", g$n))
    if (file.exists(f)) { cat("skip n=", g$n, "\n"); next }
    t0 <- Sys.time()
    M <- do.call(rbind, mclapply(1:g$R, function(r)
      tryCatch(one_rep(g$n, seed = 1e6 * g$n + r, phi0 = phi0), error = function(e) NULL),
      mc.cores = 11))
    saveRDS(list(n = g$n, M = M), f)
    cat(sprintf("n=%4d done: %d reps, %.1f min\n", g$n, nrow(M),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  cat("ALL_DONE\n")
}

## ---- Regime B: correct specification (delta_u = 0 => U fully recoverable) ----
if (MODE == "calibB") {
  library(parallel)
  R <- 40; n <- 300
  M <- do.call(rbind, mclapply(1:R, function(r) {
    x <- tryCatch(one_rep(n, seed = 6000 + r, delta_u = 0), error = function(e) NULL)
    if (is.null(x)) NULL else x[c("phi_sigma2","phi_theta","phi_nu","phi_eps")]
  }, mc.cores = 10))
  phi0 <- list(sigma2 = median(M[,1]), theta = median(M[,2]),
               nu = 0.5, sigma2_eps = median(M[,4]))
  saveRDS(phi0, file.path(OUT, "phi0B.rds"))
  cat("phi0B (delta_u=0):\n"); str(phi0)
}

if (MODE == "runB") {
  library(parallel)
  RNGkind("L'Ecuyer-CMRG")
  phi0 <- readRDS(file.path(OUT, "phi0B.rds"))
  grid <- list(list(n = 100, R = 160), list(n = 150, R = 160),
               list(n = 225, R = 160), list(n = 350, R = 130),
               list(n = 500, R = 90))
  for (g in grid) {
    f <- file.path(OUT, sprintf("resB_n%04d.rds", g$n))
    if (file.exists(f)) { cat("skip n=", g$n, "\n"); next }
    t0 <- Sys.time()
    M <- do.call(rbind, mclapply(1:g$R, function(r)
      tryCatch(one_rep(g$n, seed = 2e6 * g$n + r, phi0 = phi0, delta_u = 0),
               error = function(e) NULL), mc.cores = 11))
    saveRDS(list(n = g$n, M = M), f)
    cat(sprintf("B n=%4d done: %d reps, %.1f min\n", g$n, nrow(M),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  cat("ALL_DONE_B\n")
}
