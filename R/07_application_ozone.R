## Table 6 + Figure 7 (forest) + bias-bound diagnostic + sensitivity.
## The ozone data ship with the spaci package (ppm; multiplied to ppb here).
suppressMessages(library(spaci)); suppressMessages(library(readxl))
SEED <- 115   # declared matching seed used in the paper (see README)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
p <- system.file("extdata", "analysis_dat.xlsx", package = "spaci")
dat <- as.data.frame(read_excel(p))
out <- "mean4maxOzone"; trt <- "SnCR"; no2 <- "totNOxemissions"
cn <- c("Fac.Longitude", "Fac.Latitude")
cov <- setdiff(names(dat), c(out, trt, no2, cn))
dat <- dat[complete.cases(dat[, c(out, trt, cn, cov)]), ]
Y <- dat[[out]] * 1000; Z <- dat[[trt]]
S <- as.matrix(dat[, cn]); X <- as.matrix(dat[, cov])

## --- Table 6: matching rows (declared seed) + DR rows (geoR engine) ---------
f_nv <- naive_ps(Y, Z, X, caliper = 0.25, seed = SEED)
f_dp <- daps(Y, Z, X, S, caliper = 0.25, seed = SEED)
f_id <- idaps(Y, Z, X, S, tau = 0.2, caliper = 0.25, seed = SEED)
f_ru <- recoverU(Y, Z, X, S, tau = 0.2, matern_method = "geoR")
f_rp <- recoverUplus(Y, Z, X, S, tau = 0.2, matern_method = "geoR")
h_ru <- vcov_hac(f_ru); h_rp <- vcov_hac(f_rp)   # spatial HAC intervals
tab <- data.frame(
  method = c("Naive PS", "DAPS", "iDAPS", "recoverU", "recoverU+"),
  est = c(f_nv$att, f_dp$att, f_id$att, f_ru$att, f_rp$att),
  lo  = c(f_nv$ci[1], f_dp$ci[1], f_id$ci[1], h_ru$ci[1], h_rp$ci[1]),
  hi  = c(f_nv$ci[2], f_dp$ci[2], f_id$ci[2], h_ru$ci[2], h_rp$ci[2]))
print(tab, row.names = FALSE); saveRDS(tab, "results/table6_ozone.rds")

## --- 40-seed sensitivity ranges for the matching estimators -----------------
rng <- sapply(1:40, function(s) c(
  naive_ps(Y, Z, X, caliper = 0.25, seed = s)$att,
  daps(Y, Z, X, S, caliper = 0.25, seed = s)$att,
  idaps(Y, Z, X, S, tau = 0.2, caliper = 0.25, seed = s)$att))
cat("40-seed ranges (naive, DAPS, iDAPS):\n"); print(t(apply(rng, 1, range)))

## --- Prop. 2 bias-bound diagnostic for the iDAPS match -----------------------
print(bias_bound(f_id, Y, Z, X, S, tau = 0.2))

## --- sensitivity: tau (recoverU+, iDAPS) and caliper (iDAPS) -----------------
for (tt in c(0.1, 0.2, 0.4)) {
  cat(sprintf("recoverU+ tau=%.1f: %.3f\n", tt,
      recoverUplus(Y, Z, X, S, tau = tt, matern_method = "geoR")$att))
  cat(sprintf("iDAPS     tau=%.1f: %.3f\n", tt,
      idaps(Y, Z, X, S, tau = tt, caliper = 0.25, seed = SEED)$att))
}
for (cc in c(0.1, 0.25, 0.5)) cat(sprintf("iDAPS caliper=%.2f: %.3f\n", cc,
      idaps(Y, Z, X, S, tau = 0.2, caliper = cc, seed = SEED)$att))

## --- Figure 7: forest plot ---------------------------------------------------
d <- tab[rev(seq_len(nrow(tab))), ]
png("figures/ozone-forest.png", width = 880, height = 520, res = 118)
op <- par(mar = c(4.4, 7.2, 2.5, 1.5)); yy <- seq_len(nrow(d))
plot(d$est, yy, xlim = range(-1.4, d$lo, d$hi), ylim = c(0.6, nrow(d) + 0.4),
     pch = NA, yaxt = "n", ylab = "", bty = "n",
     xlab = "Effect of SCR/SNCR on ambient ozone (ppb)",
     main = "Estimated ATT-type direct effect with 95% intervals")
axis(2, at = yy, labels = d$method, las = 1, tick = FALSE)
abline(v = 0, lty = 2, col = "grey45")
segments(d$lo, yy, d$hi, yy, col = "#4C72B0", lwd = 2)
segments(d$lo, yy - .1, d$lo, yy + .1, col = "#4C72B0", lwd = 2)
segments(d$hi, yy - .1, d$hi, yy + .1, col = "#4C72B0", lwd = 2)
points(d$est, yy, pch = 19, col = "#1A1A1A", cex = 1.15)
text(d$est, yy + 0.27, sprintf("%.2f", d$est), cex = 0.82)
par(op); invisible(dev.off())
cat("figures/ozone-forest.png written\n")
