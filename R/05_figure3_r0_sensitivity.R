###############################################################################
## 05_figure3_r0_sensitivity.R
##
## Figure 3. R0 estimates under alternative assumptions.
##   (A) Model A: assumed latest date of exponential growth (5, 6, 7 June).
##   (B) Model B: assumed effective date of border closure (26, 27, 28 May).
## Grey bands show the systematically reviewed R0 range for Zaire/Sudan
## ebolavirus across all studies and among studies conducted in Uganda.
## Requires output/model_A_results.rds and output/model_B_results.rds.
###############################################################################

source("R/00_utils.R")

model_A_results <- readRDS("output/model_A_results.rds")
model_B_results <- readRDS("output/model_B_results.rds")

y_max <- 15
y_breaks <- seq(0, y_max, 1)

reference_ranges <- data.frame(
  range = c("All countries", "Uganda"),
  lower = c(PUBLISHED_R0_ALL_LOWER, PUBLISHED_R0_UGANDA_LOWER),
  upper = c(PUBLISHED_R0_ALL_UPPER, PUBLISHED_R0_UGANDA_UPPER),
  fill  = c("grey80", "grey60")
)

reference_bands <- list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = PUBLISHED_R0_ALL_LOWER, ymax = PUBLISHED_R0_ALL_UPPER,
           fill = "grey80", alpha = 0.8),
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = PUBLISHED_R0_UGANDA_LOWER, ymax = PUBLISHED_R0_UGANDA_UPPER,
           fill = "grey60", alpha = 0.8)
)

theme_figure3 <- theme_publication +
  theme(
    axis.text.x = element_text(size = 16, hjust = 0.5, margin = margin(t = 6)),
    axis.title.x = element_text(size = 20, margin = margin(t = 5))
  )

# --------------------------------------------------------------------------
# Panel A: Model A, R0 by assumed end date
# --------------------------------------------------------------------------
fig3A_data <- model_A_results$summary |>
  dplyr::arrange(end_day) |>
  dplyr::mutate(
    date_label = factor(format_calendar_date(START_DATE + end_day - 1L),
                         levels = format_calendar_date(START_DATE + sort(unique(end_day)) - 1L)),
    R0_lower_plot = pmax(0, R0_lower)
  )

panel_3A <- ggplot(fig3A_data, aes(date_label, R0_hat)) +
  reference_bands +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  geom_errorbar(aes(ymin = R0_lower_plot, ymax = R0_upper), width = 0.14, linewidth = 1.25, colour = "black") +
  geom_point(size = 5, colour = "black") +
  scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "The latest date of exponential growth", y = "Basic reproduction number") +
  theme_figure3

# --------------------------------------------------------------------------
# Panel B: Model B, R0 by assumed border-closure date
# --------------------------------------------------------------------------
fig3B_data <- model_B_results$summary |>
  dplyr::arrange(closure_day) |>
  dplyr::mutate(
    date_label = factor(format_calendar_date(START_DATE + closure_day - 1L),
                         levels = format_calendar_date(START_DATE + sort(unique(closure_day)) - 1L)),
    R0_lower_plot = pmax(0, R0_lower)
  )

panel_3B <- ggplot(fig3B_data, aes(date_label, R0_hat)) +
  reference_bands +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  geom_errorbar(aes(ymin = R0_lower_plot, ymax = R0_upper), width = 0.14, linewidth = 1.25, colour = "black") +
  geom_point(size = 5, colour = "black") +
  scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "The effective date of border closure", y = "Basic reproduction number") +
  theme_figure3

Figure_3 <- (
  (panel_3A | panel_3B) + plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "A")
) & theme(plot.tag = element_text(size = 28, face = "bold", colour = "black"), plot.tag.position = c(0, 1))

print(Figure_3)

dir.create("output", showWarnings = FALSE)
ggsave("output/Figure3_R0_sensitivity.png", Figure_3, width = 14, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure3_R0_sensitivity.tiff", Figure_3, width = 14, height = 6, dpi = 600, bg = "white", compression = "lzw")
cat("[Figure 3] saved to output/\n")
