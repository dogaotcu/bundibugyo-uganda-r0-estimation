###############################################################################
## 04_figure2_model_fit.R
##
## Figure 2. Observed vs model-predicted daily imported cases.
##   (A) Model A, fitted through the day-23 (6 June 2026) end date.
##   (B) Model B, assuming the reported border closure on day 13 (27 May).
## Requires output/model_A_results.rds and output/model_B_results.rds
## (run R/01_model_A.R and R/02_model_B.R first, or use R/run_all.R).
###############################################################################

source("R/00_utils.R")

model_A_results <- readRDS("output/model_A_results.rds")
model_B_results <- readRDS("output/model_B_results.rds")

reported_closure_date <- as.Date("2026-05-27")

# --------------------------------------------------------------------------
# Assemble a common observed series and shared axis limits
# --------------------------------------------------------------------------
fig2A_data <- model_A_results$fitted |>
  dplyr::filter(end_day == model_A_results$main_end_day) |>
  dplyr::mutate(calendar_date = START_DATE + day - 1L)

fig2B_data <- model_B_results$fitted |>
  dplyr::filter(closure_day == model_B_results$main_closure_day) |>
  dplyr::mutate(calendar_date = as.Date(date))

full_observed <- fig2B_data |>
  dplyr::select(calendar_date, observed_daily) |>
  dplyr::arrange(calendar_date)

date_breaks <- make_weekly_breaks(min(full_observed$calendar_date), max(full_observed$calendar_date))
y_max <- max(1, ceiling(max(full_observed$observed_daily, fig2A_data$fitted_daily_upper,
                             fig2B_data$fitted_upper, na.rm = TRUE)))

shared_x_scale <- scale_x_date(breaks = date_breaks, labels = format_calendar_date,
                                limits = c(min(full_observed$calendar_date), max(full_observed$calendar_date)),
                                expand = expansion(mult = c(0.01, 0.01)))
shared_y_scale <- scale_y_continuous(breaks = seq(0, y_max, 1), limits = c(-0.2, y_max),
                                      expand = expansion(mult = c(0, 0.05)))

theme_figure2 <- theme_publication +
  theme(
    axis.text.x = element_text(size = 16, hjust = 0.5, margin = margin(t = 6)),
    axis.title.x = element_text(size = 22, margin = margin(t = 3)),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    legend.key = element_rect(fill = "white", color = NA)
  )

# --------------------------------------------------------------------------
# Panel A: Model A fit (day 23)
# --------------------------------------------------------------------------
panel_2A <- ggplot() +
  geom_ribbon(data = fig2A_data, aes(calendar_date, ymin = fitted_daily_lower, ymax = fitted_daily_upper),
              fill = "grey75", alpha = 0.55, show.legend = FALSE) +
  geom_line(data = fig2A_data, aes(calendar_date, fitted_daily, colour = "Estimated cases"), linewidth = 1.25) +
  geom_point(data = full_observed, aes(calendar_date, observed_daily, colour = "Observed cases"), size = 3.1) +
  scale_colour_manual(name = NULL, breaks = c("Observed cases", "Estimated cases"),
                       values = c("Observed cases" = "black", "Estimated cases" = "red")) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, NA), linetype = c(0, 1), linewidth = c(0, 1.4)))) +
  shared_x_scale + shared_y_scale +
  labs(x = "Calendar time", y = "Daily number of imported cases") +
  theme_figure2 +
  theme(legend.position = "inside", legend.position.inside = c(0.07, 0.93), legend.justification = c(0, 1))

# --------------------------------------------------------------------------
# Panel B: Model B fit (closure day 13, 27 May)
# --------------------------------------------------------------------------
panel_2B <- ggplot(fig2B_data, aes(calendar_date)) +
  geom_ribbon(aes(ymin = fitted_daily_lower, ymax = fitted_daily_upper), fill = "grey75", alpha = 0.55) +
  geom_line(aes(y = fitted_daily), linewidth = 1.25, colour = "red") +
  geom_point(aes(y = observed_daily), size = 3.1, colour = "black") +
  geom_vline(xintercept = reported_closure_date, linetype = "dashed", linewidth = 1, colour = "black") +
  shared_x_scale + shared_y_scale +
  labs(x = "Calendar time", y = "Daily number of imported cases") +
  theme_figure2 +
  theme(legend.position = "none")

Figure_2 <- (
  (panel_2A | panel_2B) +
    plot_layout(widths = c(1, 1), axis_titles = "collect_x") +
    plot_annotation(tag_levels = "A")
) & theme(plot.tag = element_text(size = 28, face = "bold", colour = "black"), plot.tag.position = c(0, 1))

print(Figure_2)

dir.create("output", showWarnings = FALSE)
ggsave("output/Figure2_model_fit.png", Figure_2, width = 14, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure2_model_fit.tiff", Figure_2, width = 14, height = 6, dpi = 600, bg = "white", compression = "lzw")
cat("[Figure 2] saved to output/\n")
