###############################################################################
## 03_figure1_epicurve_map.R
##
## Figure 1. (A) Daily imported cases in Uganda with vertical markers for
## the first imported case, the PHEIC declaration, and border closure.
## (B) Map of the DRC and Uganda, with Ituri province and Kampala marked.
###############################################################################

source("R/00_utils.R")
library(sf)
library(rnaturalearth)
library(geodata)

observed <- load_observed_cases()

events <- data.frame(
  date  = as.Date(c("2026-05-15", "2026-05-17", "2026-05-27")),
  label = c("First case reported in Uganda", "Declared PHEIC", "Border closure")
)

# --------------------------------------------------------------------------
# Panel A: epidemic curve
# --------------------------------------------------------------------------
panel_A <- ggplot(observed, aes(x = date, y = daily_cases)) +
  geom_col(fill = "black", width = 0.9) +
  geom_segment(data = events, aes(x = date, xend = date, y = 0, yend = 2.5),
               color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_text(data = events, aes(x = date, y = 2.5, label = label),
            angle = 90, hjust = 0, vjust = 0.5, size = 3.2, fontface = "bold") +
  scale_x_date(
    breaks = c(seq(min(observed$date), max(observed$date) - 3, by = "2 days"), max(observed$date)),
    date_labels = "%d %b"
  ) +
  scale_y_continuous(limits = c(0, 5), breaks = 0:4, expand = c(0, 0)) +
  labs(tag = "A", x = "Calendar date", y = "Daily number of reported cases") +
  theme_classic(base_size = 14) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

# --------------------------------------------------------------------------
# Panel B: DRC / Uganda map with Ituri province and Kampala
# --------------------------------------------------------------------------
world  <- ne_countries(scale = "medium", returnclass = "sf")
africa <- world[world$continent == "Africa", ]
uganda <- africa[africa$name == "Uganda", ]
drc    <- africa[africa$name == "Dem. Rep. Congo", ]

drc_provinces <- st_as_sf(gadm(country = "COD", level = 1, path = tempdir()))
ituri <- drc_provinces[drc_provinces$NAME_1 == "Ituri", ]

kampala <- data.frame(city = "Kampala", lon = 32.5825, lat = 0.3136)

panel_B_main <- ggplot() +
  geom_sf(data = africa, fill = "grey90", color = "white", linewidth = 0.2) +
  geom_sf(data = drc, fill = "#FCAE91", color = "white", linewidth = 0.2) +
  geom_sf(data = drc_provinces, fill = NA, color = "grey40", linewidth = 0.25) +
  geom_sf(data = ituri, fill = "#CB181D", color = "white", linewidth = 0.45) +
  geom_sf(data = uganda, fill = "#FEE5D9", color = "white", linewidth = 0.2) +
  geom_sf_text(data = ituri, aes(label = NAME_1), size = 3.8, fontface = "bold",
               color = "white", nudge_y = -0.25) +
  geom_point(data = kampala, aes(x = lon, y = lat), color = "black", size = 2) +
  geom_text(data = kampala, aes(x = lon, y = lat, label = city),
            hjust = -0.2, vjust = 0, size = 3.5, fontface = "bold", color = "black") +
  coord_sf(xlim = c(22, 36), ylim = c(-5, 5.85), expand = FALSE) +
  labs(tag = "B") +
  theme_void(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))

panel_B_inset <- ggplot() +
  geom_sf(data = africa, fill = "grey90", color = "white", linewidth = 0.1) +
  geom_sf(data = drc, fill = "#FCAE91", color = NA) +
  geom_sf(data = uganda, fill = "#FEE5D9", color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

panel_B <- panel_B_main +
  inset_element(panel_B_inset, left = 0.65, bottom = 0.02, right = 0.98, top = 0.35)

legend_panel <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.06, ymin = 0.72, ymax = 0.92, fill = "#FCAE91", color = "black", linewidth = 0.3) +
  annotate("text", x = 0.08, y = 0.82, label = "DRC", hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold") +
  annotate("rect", xmin = 0.02, xmax = 0.06, ymin = 0.44, ymax = 0.64, fill = "#CB181D", color = "black", linewidth = 0.3) +
  annotate("text", x = 0.08, y = 0.54, label = "Ituri province where the transmission is most intense",
           hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold") +
  annotate("rect", xmin = 0.02, xmax = 0.06, ymin = 0.16, ymax = 0.36, fill = "#FEE5D9", color = "black", linewidth = 0.3) +
  annotate("text", x = 0.08, y = 0.26, label = "Uganda", hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, b = 0))

right_column <- (panel_B / legend_panel) + plot_layout(ncol = 1, heights = c(1, 0.22))

Figure_1 <- panel_A + right_column + plot_layout(ncol = 2, widths = c(1.2, 1))

print(Figure_1)

dir.create("output", showWarnings = FALSE)
ggsave("output/Figure1_epicurve_map.png", Figure_1, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave("output/Figure1_epicurve_map.tiff", Figure_1, width = 12, height = 5.5, dpi = 300, bg = "white", compression = "lzw")
cat("[Figure 1] saved to output/\n")
