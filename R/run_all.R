###############################################################################
## run_all.R
##
## Reproduces all analyses and figures from scratch. Run from the repository
## root, e.g.:
##
##   Rscript R/run_all.R
##
## or, from within R with the working directory set to the repo root:
##
##   source("R/run_all.R")
###############################################################################

message(">> Fitting Model A (pre-closure exponential-growth model)...")
source("R/01_model_A.R")

message(">> Fitting Model B (hybrid traveler / border-closure model)...")
source("R/02_model_B.R")

message(">> Building Figure 1 (epidemic curve and map)...")
source("R/03_figure1_epicurve_map.R")

message(">> Building Figure 2 (observed vs fitted cases)...")
source("R/04_figure2_model_fit.R")

message(">> Building Figure 3 (R0 sensitivity to assumed dates)...")
source("R/05_figure3_r0_sensitivity.R")

message(">> Building Figure 4 (R0 sensitivity to generation time)...")
source("R/06_figure4_GI_sensitivity.R")

message(">> Done. See the output/ directory for all figures and .rds results.")
