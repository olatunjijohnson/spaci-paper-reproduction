OUT <- "results"
TRUE_ATT <- 2; z <- qnorm(0.975)

## ---------- (A) COVERAGE / SE CALIBRATION ----------
cf <- c(file.path(OUT, "coverage.rds"), file.path(OUT, "coverage2.rds"))
cells <- unlist(lapply(cf[file.exists(cf)], readRDS), recursive = FALSE)
if (length(cells) > 0) {
  S <- do.call(rbind, lapply(cells, function(cell) {
    M <- cell$M; M <- M[is.finite(M[,"att"]), , drop = FALSE]
    plim <- mean(M[,"att"]); mcsd <- sd(M[,"att"])
    se <- function(lo, hi) mean((M[,hi] - M[,lo]) / (2 * z))
    cov_of <- function(lo, hi, t) mean(M[,lo] <= t & t <= M[,hi])
    data.frame(delta_u = cell$delta_u, n = cell$n, reps = nrow(M),
      bias = plim - TRUE_ATT, mcsd = mcsd,
      r_iid = se("iid_lo","iid_hi")/mcsd, r_hac = se("hac_lo","hac_hi")/mcsd,
      r_boot = se("boot_lo","boot_hi")/mcsd,
      cov_iid = cov_of("iid_lo","iid_hi", plim),
      cov_hac = cov_of("hac_lo","hac_hi", plim),
      cov_boot = cov_of("boot_lo","boot_hi", plim))
  }))
  S <- S[order(S$n, S$delta_u), ]
  cat("=== (A) recoverU+ SE calibration (SE / true Monte-Carlo SD; 1.0 = calibrated) ===\n")
  cat("    and coverage of the estimator's probability limit (nominal 0.95)\n\n")
  print(round(S[, c("delta_u","n","reps","bias","mcsd","r_iid","r_hac","r_boot",
                    "cov_iid","cov_hac","cov_boot")], 3), row.names = FALSE)

  ## figure: SE/MCsd vs delta_u at n=150
  s150 <- S[S$n == 150, ]; s150 <- s150[order(s150$delta_u), ]
  png(file.path(OUT, "coverage.png"), width = 820, height = 520, res = 118)
  op <- par(mar = c(4.5, 4.5, 3, 1))
  yl <- range(0.4, s150$r_iid, s150$r_hac, s150$r_boot, 1.4)
  plot(s150$delta_u, s150$r_iid, type = "b", pch = 19, col = "#C44E52", ylim = yl,
       xlab = "spatial confounding strength (delta_u)",
       ylab = "estimated SE / true SD", main = "recoverU+ standard-error calibration (n=150)")
  lines(s150$delta_u, s150$r_hac, type = "b", pch = 17, col = "#DD8452")
  lines(s150$delta_u, s150$r_boot, type = "b", pch = 15, col = "#4C72B0")
  abline(h = 1, lty = 2, col = "grey40")
  legend("topleft", bty = "n", pch = c(19,17,15), col = c("#C44E52","#DD8452","#4C72B0"),
         legend = c("i.i.d.","spatial HAC","block bootstrap"), cex = 0.9)
  text(max(s150$delta_u), 1.02, "calibrated (=1)", pos = 2, cex = 0.75, col = "grey40")
  par(op); invisible(dev.off())
  cat("\nFigure:", file.path(OUT, "coverage.png"), "\n")
  saveRDS(S, file.path(OUT, "coverage_summary.rds"))
}

## ---------- (B) WEIGHTS ----------
if (file.exists(file.path(OUT, "weights.rds"))) {
  wt <- readRDS(file.path(OUT, "weights.rds"))
  W <- do.call(rbind, lapply(wt, function(c) {
    M <- c$M; M <- M[is.finite(M[,"pi1"]), , drop = FALSE]
    data.frame(u = c$u, gamma = c$gamma,
               pi1 = mean(M[,"pi1"]), pi2 = mean(M[,"pi2"]), pi3 = mean(M[,"pi3"]))
  }))
  W <- W[order(W$u, W$gamma), ]
  cat("\n=== (B) mean iDAPS weights over (u, gamma):  pi1=PS, pi2=spatial, pi3=interference ===\n")
  print(round(W, 3), row.names = FALSE)
  ## the weights adapt to CONFOUNDING strength u (a design quantity), not to the
  ## interference outcome-coefficient gamma
  cat(sprintf("\npi1 (PS weight)        vs u:  cor = %+.2f  (expect < 0: PS becomes inadequate)\n", cor(W$u, W$pi1)))
  cat(sprintf("pi3 (interference wt)  vs u:  cor = %+.2f  (expect > 0)\n", cor(W$u, W$pi3)))
  cat(sprintf("pi3 (interference wt)  vs gamma: cor = %+.2f  (design weight ~ independent of outcome coef)\n", cor(W$gamma, W$pi3)))

  ## average over gamma; grouped bars of (pi1,pi2,pi3) by confounding strength u
  agg <- aggregate(cbind(pi1, pi2, pi3) ~ u, data = W, FUN = mean)
  png(file.path(OUT, "weights.png"), width = 760, height = 500, res = 118)
  op <- par(mar = c(4.5, 4.5, 3, 1))
  bp <- barplot(t(as.matrix(agg[, c("pi1","pi2","pi3")])), beside = TRUE,
                names.arg = paste0("u=", agg$u), ylim = c(0, 0.6),
                col = c("#4C72B0","#DD8452","#C44E52"),
                ylab = "mean iDAPS weight",
                main = "iDAPS weights adapt to spatial confounding strength")
  legend("topright", bty = "n", fill = c("#4C72B0","#DD8452","#C44E52"),
         legend = c(expression(pi[1]~"(propensity score)"),
                    expression(pi[2]~"(spatial proximity)"),
                    expression(pi[3]~"(neighbourhood exposure)")), cex = 0.85)
  par(op); invisible(dev.off())
  cat("Figure:", file.path(OUT, "weights.png"), "\n")
  saveRDS(W, file.path(OUT, "weights_summary.rds"))
}
