## Analyze the T2 bias-bound validation: identity, bound validity, C_overlap,
## plug-in agreement, and the caliper-sweep consistency corollary.
OUT <- "results"

## ---------- GRID ----------
grid <- readRDS(file.path(OUT, "grid.rds"))
G <- do.call(rbind, lapply(grid, function(cell) {
  M <- cell$M; M <- M[is.finite(M[,"bias"]), , drop = FALSE]
  data.frame(
    theta_spatial = cell$theta_spatial, gamma = cell$gamma, reps = nrow(M),
    mean_bias = mean(M[,"bias"]), mean_decsum = mean(M[,"dec_sum"]),
    mean_bound = mean(M[,"bound"]),
    frac_bound_ok = mean(M[,"bound"] >= abs(M[,"bias"])),
    mean_bdU = mean(M[,"bd_U"]), mean_absdecU = mean(abs(M[,"dec_U"])),
    frac_bdU_ok = mean(M[,"bd_U"] >= abs(M[,"dec_U"])),
    C_overlap = mean(abs(M[,"dec_U"])) / mean(M[,"bd_U"]),
    C_overlap_q95 = quantile(abs(M[,"dec_U"]) / pmax(M[,"bd_U"], 1e-8), 0.95),
    plugin_bdU = mean(M[,"pbd_U"]),                 # plug-in vs true bd_U
    mean_d = mean(M[,"mean_d"]), n_match = mean(M[,"n_match"])
  )
}))
cat("=== GRID: identity + bound validity ===\n")
print(round(G[, c("theta_spatial","gamma","mean_bias","mean_decsum",
                  "mean_bound","frac_bound_ok")], 3), row.names = FALSE)
cat("\n=== U-term: raw bound validity + selection factor C_overlap ===\n")
print(round(G[, c("theta_spatial","gamma","mean_bdU","mean_absdecU",
                  "frac_bdU_ok","C_overlap","C_overlap_q95","plugin_bdU")], 3), row.names = FALSE)
cat(sprintf("\nIdentity check: max |mean_bias - mean_decsum| over cells = %.4f (should be ~0)\n",
            max(abs(G$mean_bias - G$mean_decsum))))
cat(sprintf("Total-bound validity: min frac(bound>=|bias|) = %.3f\n", min(G$frac_bound_ok)))
cat(sprintf("Raw U-term bound: min frac(bd_U>=|dec_U|) = %.3f  =>  V2 needed if <1\n",
            min(G$frac_bdU_ok)))
cat(sprintf("Selection factor C_overlap: range %.2f - %.2f (mean %.2f); q95 up to %.2f\n",
            min(G$C_overlap), max(G$C_overlap), mean(G$C_overlap), max(G$C_overlap_q95)))
cat(sprintf("Plug-in vs true bd_U: mean ratio %.3f\n", mean(G$plugin_bdU / G$mean_bdU)))

## ---------- CALIPER SWEEP ----------
if (file.exists(file.path(OUT, "caliper.rds"))) {
  cal <- readRDS(file.path(OUT, "caliper.rds"))
  Cw <- do.call(rbind, lapply(cal, function(cc) {
    M <- cc$M; M <- M[is.finite(M[,"bias"]), , drop = FALSE]
    data.frame(caliper = cc$caliper, reps = nrow(M), n_match = mean(M[,"n_match"]),
      mean_d = mean(M[,"mean_d"]), mean_absbias = mean(abs(M[,"bias"])),
      mean_bdU = mean(M[,"bd_U"]), mean_absdecU = mean(abs(M[,"dec_U"])),
      mean_bound = mean(M[,"bound"]))
  }))
  Cw <- Cw[order(Cw$mean_d), ]
  cat("\n=== CALIPER SWEEP (consistency corollary): finer matching -> smaller confounding ===\n")
  print(round(Cw, 3), row.names = FALSE)

  png(file.path(OUT, "t2_caliper.png"), width = 1000, height = 460, res = 118)
  op <- par(mfrow = c(1,2), mar = c(4.3,4.3,3,1))
  ## panel 1: confounding bound term & realized |bias| shrink with matched distance
  plot(Cw$mean_d, Cw$mean_bdU, type="b", pch=19, col="#C44E52",
       xlab="mean matched distance", ylab="magnitude",
       main="Confounding bias shrinks with matched distance", cex.main=0.9,
       ylim=range(0, Cw$mean_bdU, Cw$mean_absbias, Cw$mean_absdecU))
  lines(Cw$mean_d, Cw$mean_absdecU, type="b", pch=17, col="#DD8452")
  lines(Cw$mean_d, Cw$mean_absbias, type="b", pch=15, col="#4C72B0")
  legend("topleft", bty="n", cex=0.8,
    legend=c("bound term |thU|E sqrt(2 gamma_U(d))","actual |thU E[dU]|","realized |iDAPS bias|"),
    pch=c(19,17,15), col=c("#C44E52","#DD8452","#4C72B0"))
  ## panel 2: total bound dominates |bias| across calipers
  plot(Cw$mean_d, Cw$mean_bound, type="b", pch=19, col="#1A1A1A",
       xlab="mean matched distance", ylab="magnitude",
       main="Total bound dominates realized bias", cex.main=0.9,
       ylim=range(0, Cw$mean_bound, Cw$mean_absbias))
  lines(Cw$mean_d, Cw$mean_absbias, type="b", pch=15, col="#4C72B0")
  legend("topleft", bty="n", cex=0.8, legend=c("estimated bound","realized |bias|"),
    pch=c(19,15), col=c("#1A1A1A","#4C72B0"))
  par(op); invisible(dev.off())
  cat("\nFigure:", file.path(OUT, "t2_caliper.png"), "\n")
}
saveRDS(list(grid = G, caliper = if (exists("Cw")) Cw else NULL),
        file.path(OUT, "summary.rds"))
