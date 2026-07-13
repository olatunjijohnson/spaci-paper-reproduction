## Figures 8-9: matched-pairs map and propensity-score overlap, restyled.
## Same data, seed and models as Table 6 / script 07; only the plotting
## changes (US basemap, muted colour-blind-safe palette, no in-figure
## titles, thinner connecting lines, panel annotations).
suppressMessages({
  library(spaci); library(readxl); library(ggplot2); library(maps)
})
SEED <- 115
dir.create("figures", showWarnings = FALSE)

p <- system.file("extdata", "analysis_dat.xlsx", package = "spaci")
dat <- as.data.frame(read_excel(p))
out <- "mean4maxOzone"; trt <- "SnCR"; no2 <- "totNOxemissions"
cn <- c("Fac.Longitude", "Fac.Latitude")
cov <- setdiff(names(dat), c(out, trt, no2, cn))
dat <- dat[complete.cases(dat[, c(out, trt, cn, cov)]), ]
Y <- dat[[out]] * 1000; Z <- dat[[trt]]
S <- as.matrix(dat[, cn]); X <- as.matrix(dat[, cov])

f_nv <- naive_ps(Y, Z, X, caliper = 0.25, seed = SEED)
f_dp <- daps(Y, Z, X, S, caliper = 0.25, seed = SEED)
f_id <- idaps(Y, Z, X, S, tau = 0.2, caliper = 0.25, seed = SEED)
f_ru <- recoverU(Y, Z, X, S, tau = 0.2, matern_method = "geoR")
f_rp <- recoverUplus(Y, Z, X, S, tau = 0.2, matern_method = "geoR")

pal <- c(Control = "#D98449", Treated = "#2C5F8A")

## ---------------------------------------------------------------- Figure 8
prep_match <- function(fit, coords, Z, method) {
  pr <- fit$extras$pairs
  seg <- data.frame(Method = method, x = coords[pr[, 1], 1], y = coords[pr[, 1], 2],
                     xend = coords[pr[, 2], 1], yend = coords[pr[, 2], 2])
  pts <- data.frame(Method = method, x = coords[, 1], y = coords[, 2],
                     Treatment = factor(Z, levels = c(0, 1), labels = c("Control", "Treated")))
  d_pair <- sqrt(rowSums((coords[pr[, 1], , drop = FALSE] - coords[pr[, 2], , drop = FALSE])^2))
  lbl <- sprintf("%s  (%d pairs, mean %.1f° apart)", method, nrow(pr), mean(d_pair))
  list(seg = seg, pts = pts, lbl = lbl)
}
mm_naive <- prep_match(f_nv, S, Z, "Naive PS")
mm_daps  <- prep_match(f_dp, S, Z, "DAPS")
mm_idaps <- prep_match(f_id, S, Z, "iDAPS")

lbls <- c(mm_naive$lbl, mm_daps$lbl, mm_idaps$lbl)
names(lbls) <- c("Naive PS", "DAPS", "iDAPS")
seg_df <- rbind(mm_naive$seg, mm_daps$seg, mm_idaps$seg)
pt_df  <- rbind(mm_naive$pts, mm_daps$pts, mm_idaps$pts)
seg_df$Method <- factor(seg_df$Method, levels = c("Naive PS", "DAPS", "iDAPS"))
pt_df$Method  <- factor(pt_df$Method,  levels = c("Naive PS", "DAPS", "iDAPS"))

us <- map_data("state")

p8 <- ggplot() +
  geom_polygon(data = us, aes(long, lat, group = group),
               fill = "grey97", colour = "grey80", linewidth = 0.25) +
  geom_segment(data = seg_df, aes(x, y, xend = xend, yend = yend),
               colour = "grey55", linewidth = 0.22, alpha = 0.4) +
  geom_point(data = pt_df, aes(x, y, colour = Treatment), size = 1.6, alpha = 0.85) +
  facet_wrap(~Method, ncol = 1, labeller = as_labeller(lbls)) +
  coord_fixed(ratio = 1.3, xlim = range(S[, 1]) + c(-2, 2), ylim = range(S[, 2]) + c(-2, 2)) +
  scale_colour_manual(values = pal) +
  labs(x = "Longitude", y = "Latitude", colour = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 9.5, hjust = 0),
        legend.position = "bottom",
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3),
        panel.spacing = unit(0.6, "lines"))

ggsave("figures/matched-pairs.png", p8, width = 5.4, height = 8.4, dpi = 300, bg = "white")
cat("figures/matched-pairs.png written\n")

## ---------------------------------------------------------------- Figure 9
mk_overlap <- function(ps, Z, method) {
  data.frame(PropensityScore = ps, Method = method,
             Treatment = factor(Z, levels = c(0, 1), labels = c("Control", "Treated")))
}
overlap_dat <- rbind(
  mk_overlap(f_nv$extras$ps, Z, "Naive PS"),
  mk_overlap(f_ru$extras$ps, Z, "recoverU"),
  mk_overlap(f_rp$extras$ps, Z, "recoverU+")
)
overlap_dat$Method <- factor(overlap_dat$Method, levels = c("Naive PS", "recoverU", "recoverU+"))

p9 <- ggplot(overlap_dat, aes(PropensityScore, fill = Treatment, colour = Treatment)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  facet_wrap(~Method, ncol = 1) +
  scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
  labs(x = "Estimated propensity score", y = "Density", fill = "", colour = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 9.5, hjust = 0),
        legend.position = "bottom",
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3))

ggsave("figures/ps-overlap.png", p9, width = 5.6, height = 5.6, dpi = 300, bg = "white")
cat("figures/ps-overlap.png written\n")
