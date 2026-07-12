## =====================================================================
## Missing simulation components (protocol §4 / F5):
##  (A) NUISANCE MISSPECIFICATION -- the actual double-robustness test.
##      Truth: outcome has +0.8*X1^2; PS has +0.5*X1*X2 (inside logit).
##      Fits: CC (both augmented), MC (PS omits interaction), CM (outcome
##      omits X1^2), MM (both omit). At u=0 (no unrecoverable confounding)
##      partial DR predicts: CC, MC, CM unbiased; MM biased.
##  (B) EXPOSURE-MAPPING MISSPECIFICATION. DGP exposure: exp(tau=0.1)
##      [correct], exp(tau=0.05) [bandwidth wrong], kNN k=5 [form wrong];
##      estimators always use exp(tau=0.1). u=0, gamma=1.5.
## Standalone (does not touch the reverted package). n=250, R=150.
## =====================================================================
suppressMessages(library(spaci))
OUT <- "results"
dir.create(OUT, showWarnings = FALSE)
library(parallel); RNGkind("L'Ecuyer-CMRG")

fit_matern_exp <- function(resid, coords) {
  d <- as.matrix(stats::dist(coords)); r <- as.numeric(resid); n <- length(r)
  v <- stats::var(r); if (!is.finite(v) || v <= 0) v <- 1
  nll <- function(p) {
    C <- matern_cov_matrix(d, exp(p[1]), exp(p[2]), 0.5) + exp(p[3]) * diag(n)
    ch <- tryCatch(chol(C), error = function(e) NULL); if (is.null(ch)) return(1e10)
    one <- rep(1, n); Ci1 <- backsolve(ch, forwardsolve(t(ch), one))
    Cir <- backsolve(ch, forwardsolve(t(ch), r)); mu <- sum(one*Cir)/sum(one*Ci1)
    rr <- r - mu; Crr <- backsolve(ch, forwardsolve(t(ch), rr))
    val <- 0.5*(2*sum(log(diag(ch))) + sum(rr*Crr) + n*log(2*pi))
    if (is.finite(val)) val else 1e10
  }
  o <- tryCatch(stats::optim(log(c(0.7*v, 0.2, 0.3*v)), nll, method = "Nelder-Mead",
                control = list(maxit = 400, reltol = 1e-8)), error = function(e) NULL)
  if (is.null(o)) return(list(sigma2 = 0.7*v, theta = 0.2, sigma2_eps = 0.3*v))
  p <- exp(o$par); list(sigma2 = p[1], theta = p[2], sigma2_eps = p[3])
}

## DGP with optional nonlinear truth and exposure-mapping variants
gen <- function(n, seed, u, gamma = 1.5, theta_s = 1.5, nl = FALSE,
                expo = c("exp01", "exp005", "knn5")) {
  expo <- match.arg(expo); set.seed(seed)
  S <- cbind(runif(n), runif(n)); d <- as.matrix(dist(S))
  X1 <- rnorm(n); X2 <- rbinom(n, 1, 0.5)
  L <- chol(exp(-d / 0.2) + 1e-8 * diag(n))
  U <- as.vector(scale(crossprod(L, rnorm(n)))); U[!is.finite(U)] <- 0
  lin <- 0.1 + 0.1*X1 + 0.2*X2 + u*U + if (nl) 0.45*(X1^2 - 1) else 0
  A <- rbinom(n, 1, pmin(pmax(plogis(lin), 1e-6), 1-1e-6))
  if (length(unique(A)) < 2) A <- rbinom(n, 1, 0.5)
  E_true <- switch(expo,
    exp01 = { K <- exp(-d/0.10); diag(K) <- 0; rs <- rowSums(K); rs[rs==0] <- 1; as.vector((K/rs) %*% A) },
    exp005 = { K <- exp(-d/0.05); diag(K) <- 0; rs <- rowSums(K); rs[rs==0] <- 1; as.vector((K/rs) %*% A) },
    knn5 = { nn <- apply(d, 1, function(x) order(x)[2:6]); sapply(1:n, function(i) mean(A[nn[,i]])) })
  Y <- 2.5 + 2*A + 1.0*X1 + 0.5*X2 + (if (nl) 0.8*X1^2 else 0) +
       gamma*E_true + theta_s*U + rnorm(n)
  list(Y = Y, A = A, X1 = X1, X2 = X2, S = S, U = U)
}

## recoverU+ with configurable nuisance formulas (working exposure: exp tau=0.1)
rup <- function(dd, ps_aug = FALSE, m0_aug = FALSE) {
  Y <- dd$Y; A <- dd$A; S <- dd$S
  E <- neighbourhood_exposure(S, A, tau = 0.1)$E   # working exposure, always
  dat <- data.frame(Y = Y, Z = A, X1 = dd$X1, X2 = dd$X2, G = E)
  ## the analyst's outcome specification governs BOTH the recovery-stage
  ## initial model and m0: "outcome model correct" means both include X1^2
  init_f <- if (m0_aug) Y ~ Z + X1 + X2 + I(X1^2) + G else Y ~ Z + X1 + X2 + G
  init <- lm(init_f, data = dat)
  ph <- fit_matern_exp(residuals(init), S)
  dm <- as.matrix(dist(S)); n <- length(Y)
  Sig <- matern_cov_matrix(dm, ph$sigma2, ph$theta, 0.5)
  V <- Sig + max(ph$sigma2_eps, 1e-8)*diag(n) + 1e-6*diag(n)
  W <- model.matrix(init); tg <- solve(t(W) %*% solve(V, W), t(W) %*% solve(V, Y))
  dat$Uhat <- spaci:::safe_scale(as.vector(Sig %*% solve(V, as.numeric(Y - W %*% tg))))
  ps_rhs <- paste("X1 + X2", if (ps_aug) "+ I(X1^2)" else "", "+ G + Uhat")
  m0_rhs <- paste("X1 + X2", if (m0_aug) "+ I(X1^2)" else "", "+ G + Uhat")
  ps <- suppressWarnings(glm(as.formula(paste("Z ~", ps_rhs)), family = binomial(),
        data = dat, control = glm.control(maxit = 100)))
  ehat <- spaci:::clip_ps(fitted(ps))
  m0 <- lm(as.formula(paste("Y ~", m0_rhs)), data = dat[dat$Z == 0, , drop = FALSE])
  m0hat <- as.numeric(predict(m0, newdata = dat))
  w <- ehat / (1 - ehat)
  sum((A - (1 - A) * w) * (Y - m0hat)) / sum(A)
}

args <- commandArgs(trailingOnly = TRUE); MODE <- if (length(args)) args[1] else "all"
R <- 150; n <- 250

if (MODE %in% c("all", "dr")) {
  cells <- list(CC = c(TRUE, TRUE), MC = c(FALSE, TRUE),
                CM = c(TRUE, FALSE), MM = c(FALSE, FALSE))
  res <- lapply(names(cells), function(nm) {
    a <- cells[[nm]]
    v <- unlist(mclapply(1:R, function(r) tryCatch({
      dd <- gen(n, seed = 70000 + r, u = 0, nl = TRUE)   # u = 0: clean DR test
      rup(dd, ps_aug = a[1], m0_aug = a[2])
    }, error = function(e) NA_real_), mc.cores = 11))
    v <- v[is.finite(v)]
    data.frame(cell = nm, ps_correct = a[1], m0_correct = a[2],
               reps = length(v), bias = mean(v) - 2, sd = sd(v))
  })
  saveRDS(do.call(rbind, res), file.path(OUT, "dr_misspec.rds"))
  cat("DR_MISSPEC_DONE\n"); print(do.call(rbind, res), row.names = FALSE)
}

if (MODE %in% c("all", "expo")) {
  res <- lapply(c("exp01", "exp005", "knn5"), function(ex) {
    M <- do.call(rbind, mclapply(1:R, function(r) tryCatch({
      dd <- gen(n, seed = 80000 + r, u = 0, expo = ex)
      ru <- rup(dd, ps_aug = FALSE, m0_aug = FALSE)   # linear truth here (nl=FALSE)
      fi <- idaps(dd$Y, dd$A, cbind(X1 = dd$X1, X2 = dd$X2), dd$S,
                  tau = 0.1, caliper = 0.25, seed = r)$att
      c(ru = ru, id = fi)
    }, error = function(e) c(ru = NA_real_, id = NA_real_)), mc.cores = 11))
    data.frame(dgp_exposure = ex,
               rup_bias = mean(M[,"ru"], na.rm = TRUE) - 2, rup_sd = sd(M[,"ru"], na.rm = TRUE),
               idaps_bias = mean(M[,"id"], na.rm = TRUE) - 2, idaps_sd = sd(M[,"id"], na.rm = TRUE))
  })
  saveRDS(do.call(rbind, res), file.path(OUT, "expo_misspec.rds"))
  cat("EXPO_MISSPEC_DONE\n"); print(do.call(rbind, res), row.names = FALSE)
}
cat("ALL_MISSPEC_DONE\n")
