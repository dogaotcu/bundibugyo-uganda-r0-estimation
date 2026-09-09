# Estimating the reproduction number of Bundibugyo virus from imported cases in Uganda

This repository contains the data and R code used to estimate the basic
reproduction number (R0) of the 2026 Bundibugyo virus disease epidemic in
the Democratic Republic of the Congo (DRC), using laboratory-confirmed
imported cases reported in Uganda.

Two complementary models are fitted to the same daily imported-case series:

- **Model A** — a simple exponential-growth traveler model fitted to the
  pre-border-closure period.
- **Model B** — a whole-period model that additionally accounts for the
  border closure on 27 May 2026 and the incubation-period delay between
  border crossing and symptom onset.

R0 is obtained from the estimated exponential growth rate `r` via the
Euler–Lotka equation, assuming a gamma-distributed generation interval.
95% confidence intervals and simultaneous confidence bands are obtained by
parametric bootstrap (1,000 replications).

## Repository structure

```
.
├── data/
│   └── imported_cases_uganda.csv   # daily/cumulative imported cases, 15 May-29 Jun 2026
├── R/
│   ├── 00_utils.R                  # shared setup, theme, R0 conversion, bootstrap/CI helpers
│   ├── 01_model_A.R                # Model A: fit + bootstrap (end days 22-24)
│   ├── 02_model_B.R                # Model B: fit + bootstrap (closure days 12-14)
│   ├── 03_figure1_epicurve_map.R   # Figure 1: epidemic curve + DRC/Uganda map
│   ├── 04_figure2_model_fit.R      # Figure 2: observed vs fitted cases (Models A & B)
│   ├── 05_figure3_r0_sensitivity.R # Figure 3: R0 vs assumed end date / closure date
│   ├── 06_figure4_GI_sensitivity.R # Figure 4: R0 sensitivity to generation-time assumptions
│   └── run_all.R                   # master script: reproduces everything, in order
└── output/                         # all figures and .rds results
```

`R/00_utils.R` holds every function shared by both models (the publication
theme, the `R0` conversion, and the parametric-bootstrap / simultaneous
confidence-band routines), so the estimation logic that is common to both
models is defined once rather than duplicated.

## Data

`data/imported_cases_uganda.csv` is the single source of truth for the
observed case series and is read by every downstream script. Columns:

| column             | description                                   |
|---------------------|-----------------------------------------------|
| `day`               | analysis day (1 = 15 May 2026)                |
| `date`              | calendar date                                 |
| `daily_cases`       | newly confirmed imported cases on that day    |
| `cumulative_cases`  | cumulative confirmed imported cases           |

Source: WHO AFRO and Africa CDC situation reports (15 May–29 June 2026).

## Reproducing the analysis

Requirements: R >= 4.3, with the packages `dplyr`, `tibble`, `ggplot2`
(>= 3.5.0), `future`, `future.apply`, `parallelly`, `scales`, `patchwork`,
and, for Figure 1 only, `sf`, `rnaturalearth`, and `geodata`.

```r
install.packages(c(
  "dplyr", "tibble", "ggplot2", "future", "future.apply",
  "parallelly", "scales", "patchwork", "sf", "rnaturalearth", "geodata"
))
```

From the repository root:

```bash
Rscript R/run_all.R
```

This fits both models, saves `output/model_A_results.rds` and
`output/model_B_results.rds`, and writes Figures 1–4 to `output/` as PNG
and TIFF. Each script can also be run individually — Figures 2–4 read the
saved `.rds` files, so Model A and Model B only need to be re-fitted when
their assumptions change.

## Model summary

**Model A** (pre-closure exponential growth):

```
E[c_t] = N/(m*q) * i0 * exp(r*t),   c_t ~ Poisson(E[c_t])
```

**Model B** (whole-period, border closure + incubation delay):

```
E[c_t] = N/(m*q) * i0 * exp(r*t)                                  for t < t1
E[c_t ]= ∫_0^(t_1)〖N/mq i_0  e^r(t_1-s)  f((t-t_1 )+s)ds 〗+ (1-ζ) ∫_(t_1)^t〖N/mq i_0  e^ru f(t-u)du〗                                      for t >= t1
```
![Uploading image.png…]()

`log(i0)` and `r` are estimated jointly by maximum likelihood (no closed
form is available once the incubation convolution is introduced).

**R0 conversion** (Euler–Lotka equation with a gamma generation interval):

```
R0 = (1 + r * sigma^2 / mu) ^ (mu^2 / sigma^2)
```

where `mu` and `sigma` are the mean and SD of the generation-interval
distribution (13.7 and 4.5 days, respectively, in the baseline analysis).

## Citation

If you use this code or data, please cite the accompanying manuscript.
(Add the full citation here once available, e.g. journal, DOI.)

## License

Released under the MIT License (see `LICENSE`).
