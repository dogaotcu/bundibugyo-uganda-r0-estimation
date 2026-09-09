###############################################################################
## 06_figure4_GI_sensitivity.R
##
## Figure 4. Sensitivity of R0 to the assumed generation-time distribution.
##   (A) Model A (headline end-day fit).
##   (B) Model B (headline closure-day fit).
## Each tile shows R0 recomputed at a fixed growth-rate estimate r-hat
## across combinations of the generation-time mean and SD.
## Requires output/model_A_results.rds and output/model_B_results.rds.
###############################################################################

source("R/00_utils.R")

model_A_results <- readRDS("output/model_A_results.rds")
model_B_results <- readRDS("output/model_B_results.rds")

r_hat_A <- model_A_results$summary |>
  dplyr::filter(end_day == model_A_results$main_end_day) |> dplyr::pull(r_hat)
r_hat_B <- model_B_results$summary |>
  dplyr::filter(closure_day == model_B_results$main_closure_day) |> dplyr::pull(r_hat)

if (length(r_hat_A) != 1L) stop("Reference r for Model A GI sensitivity is unavailable.")
if (length(r_hat_B) != 1L) stop("Reference r for Model B GI sensitivity is unavailable.")

gi_heatmap <- function(r_hat, tag) {
  results <- expand.grid(GI_mean = GI_MEAN_RANGE, GI_sd = GI_SD_RANGE) |>
    dplyr::mutate(R0 = calculate_R0(r = r_hat, GI_mean = GI_mean, GI_sd = GI_sd))

  ggplot(results, aes(factor(sprintf("%.1f", GI_mean)), factor(sprintf("%.1f", GI_sd)), fill = R0)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f", R0)), size = 3.5, color = "black") +
    scale_fill_gradient(low = "#D6EAF8", high = "#C0392B", name = expression(R[0]),
                         breaks = scales::breaks_pretty(n = 5), labels = \(x) sprintf("%.1f", x)) +
    labs(tag = tag, x = "Mean generation time (days)", y = "SD of the generation time (days)") +
    theme_publication +
    theme(
      plot.tag = element_text(size = 18, face = "bold"), plot.tag.position = c(0, 1),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, colour = "black"),
      legend.position = "right", legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.height = unit(1.4, "cm"), legend.key.width = unit(0.5, "cm")
    )
}

panel_4A <- gi_heatmap(r_hat_A, "A")
panel_4B <- gi_heatmap(r_hat_B, "B")

Figure_4 <- panel_4A + panel_4B + plot_layout(widths = c(1, 1))

print(Figure_4)

dir.create("output", showWarnings = FALSE)
ggsave("output/Figure4_GI_sensitivity.tiff", Figure_4, device = ragg::agg_tiff,
       width = 16, height = 6, units = "in", dpi = 600)
cat("[Figure 4] saved to output/\n")
