# Reproduction materials: "Unified estimation of direct causal effects with spatial confounding and interference"

This repository contains R scripts that reproduce **every table and figure**
in the paper (Ogunsola & Johnson). The methods themselves live in the R
package [`spaci`](https://github.com/Ogunsolaia/spaci); this repository holds
only the experiment and analysis code.

## Setup

```r
install.packages("remotes")
remotes::install_github("Ogunsolaia/spaci")   # methods + bundled ozone data
install.packages(c("readxl", "geoR"))          # geoR: application's Matern engine
```

R >= 4.1. Scripts write intermediate results to `results/` and figures to
`figures/`. Where a script is Monte-Carlo heavy, the replicate count can be
reduced via the environment variable noted in its header.

## Script -> output map

| Script | Reproduces | Approx. runtime* |
|---|---|---|
| `R/01_accuracy_table.R` | Table 3 (accuracy across confounding strength) | hours (NSIM=1000); set `NSIM=100` for a check |
| `R/02_bias_bound_validation.R` + `02b` | Fig. 2 and the Prop. 2 validation numbers (bound >= bias in 100% of reps; caliper sweep) | ~1 h |
| `R/03_partial_dr_validation.R` + `03b` | Fig. 4 and the Prop. 3 validation (bias tracks b_UA, cor 0.94) | ~30 min |
| `R/04_gamma_probe.R` + `04b` | Fig. 6 (recovery orthogonality probe) | ~30 min |
| `R/05_coverage_weights.R` + `05b` | Table 2 + Fig. 5 (SE calibration/coverage) and Table 1 + Fig. 3 (weight adaptivity) | ~1 h |
| `R/06_misspecification.R` | Tables 4-5 (double-robustness stress test; exposure-mapping misspecification) | ~30 min |
| `R/07_application_ozone.R` | Table 6, Fig. 7 (forest), bias-bound diagnostic, tau/caliper/seed sensitivity | ~5 min |

*on ~10 cores; scripts use `parallel::mclapply`.

## Notes

- The ozone data (Papadogeorgou 2016, Harvard Dataverse doi:10.7910/DVN/DKXXSN)
  ship inside the `spaci` package (`inst/extdata/`); the raw outcome is in ppm
  and all analyses convert to ppb.
- Matching estimators use greedy 1:1 matching whose result depends on the
  declared seed; the paper uses seed 115 throughout and reports the 40-seed
  range. `spaci` also offers deterministic optimal matching
  (`match_method = "optimal"`).
- Scripts 02-06 are the *validation* experiments; each corresponds to a
  numbered proposition or stress test in the paper.
