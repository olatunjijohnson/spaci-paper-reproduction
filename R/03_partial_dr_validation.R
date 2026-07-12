## =====================================================================
## T3: partial double robustness of recoverU+ (Proposition 3) -- validation.
##
## Prop 3: with the PS-model odds weight w = p(V)/(1-p(V)), V=(X,E,U_R),
##   bias(recoverU+) = theta_U * ( E[U_A | A=1] - E^w[U_A | A=0] ) =: b_UA
## where U_A is the confounder component the adjustment set cannot see.
## Operational (definition-robust) b_UA: r_U = residual of U on the actual
##   adjustment set (X, E, Uhat); imbalance of r_U between treated and
##   odds-weighted controls, times theta_U.
## Checks: (a) realized bias tracks b_UA (45 deg); (b) bias -> 0 as delta_u -> 0
##   (U_A -> 0, classic DR recovered); (c) bias ∝ theta_U; oracle-full-U ~ 0.
## =====================================================================
suppressMessages(library(spaci))
OUT <- "results"
dir.create(OUT, showWarnings = FALSE)

fit_matern_exp <- function(resid, coords) {
  d <- as.matrix(stats::dist(coords)); r <- as.numeric(resid); n <- length(r)
  v <- stats::var(r); if (!is.finite(v) || v <= 0) v <- 1
  nll <- function(p) {
    C <- matern_cov_matrix(d, exp(p[1]), exp(p[2]), 0.5) + exp(p[3]) * diag(n)
    ch <- tryCatch(chol(C), error = function(e) NULL); if (is.null(ch)) return(1e10)
    one <- rep(1, n); Ci1 <- backsolve(ch, forwardsolve(t(ch), one))
    Cir <- backsolve(ch, forwardsolve(t(ch), r)); mu <- sum(one*Cir)/sum(one*Ci1)
    res <- r - mu; Cires <- backsolve(ch, forwardsolve(t(ch), res))
    val <- 0.5*(2*sum(log(diag(ch))) + sum(res*Cires) + n*log(2*pi)); if (is.finite(val)) val else 1e10
  }
  o <- tryCatch(stats::optim(log(c(0.7*v,0.2,0.3*v)), nll, method="Nelder-Mead",
                control=list(maxit=400, reltol=1e-8)), error=function(e) NULL)
  if (is.null(o)) return(list(sigma2=0.7*v, theta=0.2, sigma2_eps=0.3*v))
  p <- exp(o$par); list(sigma2=p[1], theta=p[2], sigma2_eps=p[3])
}

recover_UR <- function(Y, Z, X, coords, E) {
  n <- length(Y); Xdf <- as.data.frame(X); xn <- colnames(X)
  dat <- data.frame(Y=Y, Z=Z, Xdf, G=E)
  init_form <- stats::as.formula(paste("Y ~ Z +", paste(c(xn,"G"), collapse="+")))
  fit <- stats::lm(init_form, data=dat); resid <- stats::residuals(fit)
  dr <- as.matrix(stats::dist(coords))
  ph <- fit_matern_exp(resid, coords)
  Sig <- matern_cov_matrix(dr, ph$sigma2, ph$theta, 0.5)
  V <- Sig + max(ph$sigma2_eps,1e-8)*diag(n) + 1e-6*diag(n)
  W <- stats::model.matrix(init_form, data=dat)
  tg <- solve(t(W)%*%solve(V,W), t(W)%*%solve(V,Y))
  as.vector(Sig %*% solve(V, as.numeric(Y - W%*%tg)))
}

## recoverU+ with an injectable confounder; returns att, ehat, m0hat, Uconf used
dr_full <- function(Y, Z, X, E, Uconf) {
  n <- length(Y); n1 <- sum(Z==1); Xdf <- as.data.frame(X); xn <- colnames(X)
  dat <- data.frame(Y=Y, Z=Z, Xdf, G=E, Uhat=safe_scale(Uconf))
  rhs <- c(xn,"G","Uhat")
  ps <- suppressWarnings(stats::glm(stats::as.formula(paste("Z ~", paste(rhs,collapse="+"))),
        family=stats::binomial(), data=dat, control=stats::glm.control(maxit=100)))
  ehat <- clip_ps(stats::fitted(ps))
  m0 <- stats::lm(stats::as.formula(paste("Y ~", paste(rhs,collapse="+"))), data=dat[dat$Z==0,,drop=FALSE])
  m0hat <- as.numeric(stats::predict(m0, newdata=dat))
  w <- ehat/(1-ehat); psi <- (Z-(1-Z)*w)*(Y-m0hat)
  list(att=sum(psi)/n1, ehat=ehat, uhat=dat$Uhat)
}

one_rep <- function(n, seed, theta_spatial, delta_u, gamma=1.5, tau=0.1) {
  sim <- simulate_spatial_causal(n=n, delta_u=delta_u, theta_spatial=theta_spatial,
             gamma_interference=gamma, tau_exp=tau, seed=seed)
  Y<-sim$Y; Z<-sim$Z; X<-sim$X; S<-sim$coords; U<-sim$U
  E <- neighbourhood_exposure(S, Z, tau=tau)$E
  Uhat <- tryCatch(recover_UR(Y,Z,X,S,E), error=function(e) NULL); if (is.null(Uhat)) return(NULL)

  d <- dr_full(Y, Z, X, E, Uhat)                       # real recoverU+
  bias <- d$att - 2
  ## plug-in b_UA: r_U = part of U not seen by adjustment set (X, E, Uhat)
  r_U <- stats::residuals(stats::lm(U ~ X[,1] + X[,2] + E + d$uhat))
  w <- d$ehat/(1-d$ehat)
  eE_treat <- mean(r_U[Z==1])
  eE_ctrl_w <- sum((1-Z)*w*r_U)/sum((1-Z)*w)
  b_UA <- theta_spatial * (eE_treat - eE_ctrl_w)
  ## oracle: plug true full U (no unrecoverable part) -> bias ~ 0
  bias_oracle <- dr_full(Y, Z, X, E, U)$att - 2

  c(theta_spatial=theta_spatial, delta_u=delta_u, bias=bias, b_UA=b_UA,
    bias_oracle=bias_oracle, cor_UhatU=cor(Uhat, U))
}

args <- commandArgs(trailingOnly=TRUE); MODE <- if (length(args)>=1) args[1] else "run"
library(parallel); RNGkind("L'Ecuyer-CMRG")

if (MODE == "test") {
  for (du in c(0, 1, 2)) for (th in c(1, 2)) {
    r <- one_rep(250, 11, theta_spatial=th, delta_u=du)
    cat(sprintf("theta=%.1f delta_u=%.1f: bias=%.3f b_UA=%.3f oracle=%.3f cor=%.2f\n",
        th, du, r["bias"], r["b_UA"], r["bias_oracle"], r["cor_UhatU"]))
  }
}

if (MODE == "run") {
  cells <- expand.grid(theta_spatial=c(0.4,1.0,2.0), delta_u=c(0,0.5,1.0,2.0,3.0))
  R <- 200; n <- 250
  res <- lapply(seq_len(nrow(cells)), function(k) {
    th<-cells$theta_spatial[k]; du<-cells$delta_u[k]
    M <- do.call(rbind, mclapply(1:R, function(r)
      tryCatch(one_rep(n, seed=13000*k+r, theta_spatial=th, delta_u=du), error=function(e) NULL),
      mc.cores=11))
    list(theta_spatial=th, delta_u=du, M=M)
  })
  saveRDS(res, file.path(OUT,"grid.rds")); cat("T3_DONE\n")
}
