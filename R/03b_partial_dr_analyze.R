## Analyze T3 partial-double-robustness validation.
OUT <- "results"
grid <- readRDS(file.path(OUT, "grid.rds"))

pool <- do.call(rbind, lapply(grid, function(c) {
  M <- c$M; M[is.finite(M[,"bias"]) & is.finite(M[,"b_UA"]), , drop=FALSE]
}))
G <- do.call(rbind, lapply(grid, function(c) {
  M <- c$M; M <- M[is.finite(M[,"bias"]), , drop=FALSE]
  data.frame(theta_spatial=c$theta_spatial, delta_u=c$delta_u, reps=nrow(M),
    mean_bias=mean(M[,"bias"]), mean_bUA=mean(M[,"b_UA"]),
    mean_oracle=mean(M[,"bias_oracle"]), cor_UhatU=mean(M[,"cor_UhatU"]))
}))
G <- G[order(G$theta_spatial, G$delta_u), ]

cat("=== Cell means: realized bias vs analytic b_UA vs oracle(true U) ===\n")
print(round(G, 3), row.names=FALSE)

## (a) 45-degree tracking (pooled per-rep + design-point means)
r_rep  <- cor(pool[,"bias"], pool[,"b_UA"])
fit_cell <- lm(mean_bias ~ mean_bUA, data=G)
cat(sprintf("\n(a) FORMULA TRACKS: per-rep cor(bias, b_UA) = %.3f;  design-point-mean slope = %.3f, intercept = %.3f (ideal 1, 0)\n",
    r_rep, coef(fit_cell)[2], coef(fit_cell)[1]))

## (b) delta_u -> 0 : classic DR recovered
b0 <- G[G$delta_u==0, ]
cat(sprintf("(b) delta_u=0 (U not confounding): mean|bias| = %.3f (should ~0; vs delta_u=2: %.3f)\n",
    mean(abs(b0$mean_bias)), mean(abs(G$mean_bias[G$delta_u==2]))))

## (c) oracle (full U) ~ unbiased
cat(sprintf("(c) oracle plug true U: mean|bias_oracle| over design points = %.3f (should ~0)\n",
    mean(abs(G$mean_oracle))))

## theta scaling at delta_u=2
d2 <- G[G$delta_u==2, ]
cat(sprintf("    bias ∝ theta_U at delta_u=2: theta 0.4/1/2 -> bias %.2f/%.2f/%.2f\n",
    d2$mean_bias[1], d2$mean_bias[2], d2$mean_bias[3]))

## ---- figures ----
png(file.path(OUT, "t3_dr.png"), width=1000, height=460, res=118)
op <- par(mfrow=c(1,2), mar=c(4.3,4.3,3,1))
## panel 1: 45-degree scatter (pooled reps, colored light) + design-point means (bold)
plot(pool[,"b_UA"], pool[,"bias"], pch=16, col="#B0B0B033", cex=0.5,
     xlab="analytic residual bias  b_UA", ylab="realized recoverU+ bias",
     main="Prop 3: bias tracks b_UA", cex.main=0.95,
     xlim=range(pool[,"b_UA"]), ylim=range(pool[,"bias"]))
abline(0,1,col="#C44E52",lwd=2,lty=2)
points(G$mean_bUA, G$mean_bias, pch=19, col="#1A1A1A", cex=1.2)
legend("topleft", bty="n", cex=0.8, legend=c("per-rep","design-point mean","y = x"),
       pch=c(16,19,NA), lty=c(NA,NA,2), col=c("#B0B0B0","#1A1A1A","#C44E52"))
## panel 2: bias vs delta_u for each theta_spatial (-> 0 as delta_u->0)
cols <- c("#4C72B0","#DD8452","#C44E52"); ths <- sort(unique(G$theta_spatial))
plot(NA, xlim=range(G$delta_u), ylim=range(0,G$mean_bias), xlab="delta_u (confounding of treatment)",
     ylab="mean recoverU+ bias", main="Bias -> 0 as U stops confounding", cex.main=0.95)
abline(h=0, col="grey70")
for (i in seq_along(ths)) {
  g <- G[G$theta_spatial==ths[i], ]; g <- g[order(g$delta_u),]
  lines(g$delta_u, g$mean_bias, type="b", pch=19, col=cols[i])
}
legend("topleft", bty="n", cex=0.8, legend=paste0("theta_U=",ths), pch=19, col=cols)
par(op); invisible(dev.off())
cat("\nFigure:", file.path(OUT,"t3_dr.png"), "\n")
saveRDS(G, file.path(OUT,"summary.rds"))
