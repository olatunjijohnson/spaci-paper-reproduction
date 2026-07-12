## =====================================================================
## T2: variogram bias bound for iDAPS (Proposition 2) -- validation.
##
## iDAPS matched-pair difference (linear working model):
##   Y_i - Y_j = th1 + th2'(X_i-X_j) + th3(E_i-E_j) + thU(U_i-U_j) + (eps_i-eps_j)
## => EXACT decomposition of the iDAPS bias:
##   bias = th2'E[dX] + th3 E[dE] + thU E[dU]                       (identity)
## => BOUND (Cauchy-Schwarz on the U-term, gamma_U = semivariogram):
##   |bias| <= ||th2|| E||dX|| + |th3| E|dE| + |thU| E sqrt(2 gamma_U(d_match))
## Validate: (a) identity holds; (b) bound >= |bias|; (c) confounding term
## shrinks as matched distances shrink (caliper sweep => consistency corollary).
## =====================================================================
suppressMessages(library(spaci))
OUT <- "results"
dir.create(OUT, showWarnings = FALSE)

## stable exponential (nu=0.5) Matern MLE on residuals (from T5 probe; see T8)
fit_matern_exp <- function(resid, coords) {
  d <- as.matrix(stats::dist(coords)); r <- as.numeric(resid); n <- length(r)
  v <- stats::var(r); if (!is.finite(v) || v <= 0) v <- 1
  nll <- function(p) {
    s2 <- exp(p[1]); th <- exp(p[2]); ng <- exp(p[3])
    C <- matern_cov_matrix(d, s2, th, 0.5) + ng * diag(n)
    ch <- tryCatch(chol(C), error = function(e) NULL); if (is.null(ch)) return(1e10)
    one <- rep(1, n); Ci1 <- backsolve(ch, forwardsolve(t(ch), one))
    Cir <- backsolve(ch, forwardsolve(t(ch), r)); mu <- sum(one*Cir)/sum(one*Ci1)
    res <- r - mu; Cires <- backsolve(ch, forwardsolve(t(ch), res))
    val <- 0.5*(2*sum(log(diag(ch))) + sum(res*Cires) + n*log(2*pi))
    if (is.finite(val)) val else 1e10
  }
  o <- tryCatch(stats::optim(log(c(0.7*v,0.2,0.3*v)), nll, method="Nelder-Mead",
                control=list(maxit=400, reltol=1e-8)), error=function(e) NULL)
  if (is.null(o)) return(list(sigma2=0.7*v, theta=0.2, sigma2_eps=0.3*v))
  p <- exp(o$par); list(sigma2=p[1], theta=p[2], sigma2_eps=p[3])
}

## ---- PLUG-IN bias bound (the real diagnostic for applications) ----
## Coefficients from an outcome working model; confounding term from the fitted
## SPATIAL variogram of the residual field (sill already carries thU^2 Var(U)).
bias_bound <- function(Y, Z, X, coords, pairs, tau = 0.1, E = NULL) {
  if (is.null(E)) E <- neighbourhood_exposure(coords, Z, tau = tau)$E
  Xdf <- as.data.frame(X); xn <- colnames(X)
  dat <- data.frame(Y = Y, Z = Z, Xdf, G = E)
  fit <- stats::lm(stats::as.formula(paste("Y ~ Z +", paste(c(xn,"G"), collapse="+"))), data=dat)
  co <- stats::coef(fit); th2 <- co[xn]; th3 <- unname(co["G"])
  phi <- fit_matern_exp(stats::residuals(fit), coords)
  ti <- pairs[,1]; ci <- pairs[,2]
  dX <- X[ti,,drop=FALSE] - X[ci,,drop=FALSE]; dE <- E[ti] - E[ci]
  d_pair <- sqrt(rowSums((coords[ti,,drop=FALSE] - coords[ci,,drop=FALSE])^2))
  rho <- exp(-sqrt(2) * d_pair / phi$theta)          # Matern nu=0.5 correlation
  gU  <- phi$sigma2 * (1 - rho)                       # spatial semivariance at matched d
  bd_X <- sqrt(sum(th2^2)) * mean(sqrt(rowSums(dX^2)))
  bd_E <- abs(th3) * mean(abs(dE))
  bd_U <- mean(sqrt(2 * pmax(gU, 0)))
  list(bd_X = bd_X, bd_E = bd_E, bd_U = bd_U, total = bd_X + bd_E + bd_U,
       mean_d = mean(d_pair), sigma2 = phi$sigma2, theta = phi$theta)
}

## ---- one replicate: identity + true-bound + plugin-bound + realized bias ----
one_rep <- function(n, seed, theta_spatial, gamma, caliper, tau = 0.1,
                    delta_u = 2, u_phi = 0.2) {
  sim <- simulate_spatial_causal(n = n, delta_u = delta_u, theta_spatial = theta_spatial,
             gamma_interference = gamma, tau_exp = tau, u_phi = u_phi, seed = seed)
  Y <- sim$Y; Z <- sim$Z; X <- sim$X; S <- sim$coords; U <- sim$U
  fit <- tryCatch(idaps(Y, Z, X, S, tau = tau, caliper = caliper, seed = seed),
                  error = function(e) NULL)
  if (is.null(fit) || is.null(fit$extras$pairs) || nrow(fit$extras$pairs) < 2 ||
      !is.finite(fit$att)) return(NULL)
  E <- fit$extras$E; pr <- fit$extras$pairs; ti <- pr[,1]; ci <- pr[,2]
  dX <- X[ti,,drop=FALSE] - X[ci,,drop=FALSE]; dE <- E[ti]-E[ci]; dU <- U[ti]-U[ci]
  d_pair <- sqrt(rowSums((S[ti,,drop=FALSE]-S[ci,,drop=FALSE])^2))
  th2 <- c(1, 0.5); th3 <- gamma; thU <- theta_spatial     # TRUE DGP coefficients

  bias <- fit$att - 2
  ## exact signed decomposition (true coeffs + true U)
  tX <- sum(th2 * colMeans(dX)); tE <- th3 * mean(dE); tU <- thU * mean(dU)
  ## true bound terms
  gU_true <- 1 - exp(-d_pair / u_phi)                      # gamma_U for exp field, Var=1
  bdX <- sqrt(sum(th2^2)) * mean(sqrt(rowSums(dX^2)))
  bdE <- abs(th3) * mean(abs(dE))
  bdU <- abs(thU) * mean(sqrt(2 * gU_true))
  ## plug-in bound
  pb <- bias_bound(Y, Z, X, S, pr, tau = tau, E = E)

  c(att = fit$att, bias = bias, n_match = nrow(pr), mean_d = mean(d_pair),
    dec_X = tX, dec_E = tE, dec_U = tU, dec_sum = tX + tE + tU,
    bd_X = bdX, bd_E = bdE, bd_U = bdU, bound = bdX + bdE + bdU,
    pbd_X = pb$bd_X, pbd_E = pb$bd_E, pbd_U = pb$bd_U, pbound = pb$total)
}

## ---------------- MODE ----------------
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "grid"
library(parallel); RNGkind("L'Ecuyer-CMRG")

if (MODE == "test") {
  cat("=== single-rep sanity (theta_spatial=2, gamma=1.5, caliper=0.25) ===\n")
  for (s in c(42, 43, 44)) {
    r <- one_rep(n = 250, seed = s, theta_spatial = 2.0, gamma = 1.5, caliper = 0.25)
    cat(sprintf("seed %d: att=%.3f bias=%.3f dec_sum=%.3f noise=%.3f | bound=%.3f  bound>=|bias|:%s\n",
        s, r["att"], r["bias"], r["dec_sum"], r["bias"] - r["dec_sum"], r["bound"],
        r["bound"] >= abs(r["bias"])))
    cat(sprintf("   terms: dec(X=%.3f,E=%.3f,U=%.3f)  bd(X=%.3f,E=%.3f,U=%.3f)  plugin_bd_U=%.3f\n",
        r["dec_X"], r["dec_E"], r["dec_U"], r["bd_X"], r["bd_E"], r["bd_U"], r["pbd_U"]))
  }
}

if (MODE == "grid") {
  cells <- expand.grid(theta_spatial = c(0.4, 1.0, 2.0), gamma = c(0, 1.5, 3.0))
  R <- 200; n <- 250
  res <- lapply(seq_len(nrow(cells)), function(k) {
    th <- cells$theta_spatial[k]; gm <- cells$gamma[k]
    M <- do.call(rbind, mclapply(1:R, function(r)
      tryCatch(one_rep(n, seed = 7000*k + r, theta_spatial = th, gamma = gm, caliper = 0.25),
               error = function(e) NULL), mc.cores = 11))
    list(theta_spatial = th, gamma = gm, M = M)
  })
  saveRDS(res, file.path(OUT, "grid.rds")); cat("GRID_DONE\n")
}

if (MODE == "caliper") {
  cals <- c(0.05, 0.1, 0.15, 0.25, 0.5, Inf); R <- 200; n <- 250
  res <- lapply(seq_along(cals), function(k) {
    M <- do.call(rbind, mclapply(1:R, function(r)
      tryCatch(one_rep(n, seed = 9000*k + r, theta_spatial = 2.0, gamma = 1.5, caliper = cals[k]),
               error = function(e) NULL), mc.cores = 11))
    list(caliper = cals[k], M = M)
  })
  saveRDS(res, file.path(OUT, "caliper.rds")); cat("CALIPER_DONE\n")
}
