##############################################################################
# plot_results.R — Visualize simulation results across scenarios
##############################################################################
rm(list = ls())

library(data.table)
library(plyr)
library(ggplot2)
library(ggpubr)

# ── 1. Load all result files ────────────────────────────────────────────────

files   <- list.files("outMod/", pattern = "Results_", full.names = TRUE)
rawData <- lapply(files, function(f) {
  data.frame(Year = 1:50, readRDS(f), stringsAsFactors = FALSE)
})

rawData <- rbindlist(rawData)

names(rawData) <- c("Year", "rep", "Scenario",
                    "ParMeanT1", "ParMeanT2",
                    "ParVarT1", "ParVarT2",
                    "varyMeanT1", "varyMeanT2",
                    "Inbreeding", "IbdDeltaF",
                    "h2MeanT1", "h2MeanT2")

# ── 2. Summarize across reps ───────────────────────────────────────────────

t1_summary <- ddply(rawData, c("Year", "Scenario"), summarize,
                    mean = mean(ParMeanT1, na.rm = TRUE),
                    se   = sd(ParMeanT1, na.rm = TRUE) / sqrt(length(ParMeanT1)))

t2_summary <- ddply(rawData, c("Year", "Scenario"), summarize,
                    mean = mean(ParMeanT2, na.rm = TRUE),
                    se   = sd(ParMeanT2, na.rm = TRUE) / sqrt(length(ParMeanT2)))

# ── 3. Colors ──────────────────────────────────────────────────────────────

scenario_colors <- c("Conventional"  = "cornflowerblue",
                     "Strategy1"     = "tomato2",
                     "Strategy2"     = "seagreen",
                     "Strategy3"     = "gold2")

# ── 4. Yield plot ───────────────────────────────────────────────────────────

p_yield <- ggplot(t1_summary, aes(x = Year, y = mean, color = Scenario)) +
  theme_minimal() +
  geom_point(size = 3.5) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous("Year", limits = c(19, 50)) +
  scale_y_continuous("Parent mean", expand = c(0, 0),
                     limits = c(40, 70), breaks = seq(40, 70, by = 5)) +
  ggtitle("Yield") +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_blank(),
    plot.title           = element_text(size = 26, hjust = 0.5, face = "bold"),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 20, face = "bold"),
    legend.key           = element_blank(),
    legend.key.width     = unit(2, "cm"),
    axis.text.x          = element_text(size = 20, colour = "black", face = "bold",
                                        angle = 45, hjust = 1),
    axis.text.y          = element_text(size = 20, colour = "black", face = "bold"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_text(size = 24, face = "bold"),
    axis.ticks           = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length    = unit(0.2, "cm"),
    axis.line            = element_line(color = "black", linewidth = 0.7)
  )

# ── 5. Specific gravity plot ───────────────────────────────────────────────

p_sg <- ggplot(t2_summary, aes(x = Year, y = mean, color = Scenario)) +
  theme_minimal() +
  geom_point(size = 3.5) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous("Year", limits = c(19, 50)) +
  scale_y_continuous("Parent mean", expand = c(0, 0),
                     limits = c(1.05, 1.16), breaks = seq(1.05, 1.16, by = 0.02)) +
  ggtitle("Specific Gravity") +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_blank(),
    plot.title           = element_text(size = 26, hjust = 0.5, face = "bold"),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 20, face = "bold"),
    legend.key           = element_blank(),
    legend.key.width     = unit(2, "cm"),
    axis.text.x          = element_text(size = 20, colour = "black", face = "bold",
                                        angle = 45, hjust = 1),
    axis.text.y          = element_text(size = 20, colour = "black", face = "bold"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_text(size = 24, face = "bold"),
    axis.ticks           = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length    = unit(0.2, "cm"),
    axis.line            = element_line(color = "black", linewidth = 0.7)
  )

# ── 6. Combined figure ─────────────────────────────────────────────────────

combined_plot <- ggarrange(
  p_yield, p_sg,
  labels        = c("(A)", "(B)"),
  font.label    = list(size = 20, face = "bold"),
  ncol          = 2,
  nrow          = 1,
  common.legend = TRUE,
  legend        = "bottom",
  align         = "hv"
)

tiff("Yield_SG_simulation_results.tiff",
     width = 16, height = 8, units = "in", res = 450)
print(combined_plot)
dev.off()

message("Plot saved: Yield_SG_simulation_results.tiff")

