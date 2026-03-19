##############################################################################
# plot_results.R — Visualize simulation results across scenarios
#
# Reads all Results_*.rds files from outMod/, aggregates across reps,
# and plots parent genetic mean over years for Yield and SG.
#
# Usage: Rscript plot_results.R
# (run from the simulation directory containing outMod/)
##############################################################################

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

# Yield (Trait 1)
t1_summary <- ddply(rawData, c("Year", "Scenario"), summarize,
                    mean = mean(ParMeanT1, na.rm = TRUE),
                    se   = sd(ParMeanT1, na.rm = TRUE) / sqrt(length(ParMeanT1)))

# Specific Gravity (Trait 2)
t2_summary <- ddply(rawData, c("Year", "Scenario"), summarize,
                    mean = mean(ParMeanT2, na.rm = TRUE),
                    se   = sd(ParMeanT2, na.rm = TRUE) / sqrt(length(ParMeanT2)))

# ── 3. Plotting theme ──────────────────────────────────────────────────────

scenario_colors <- c("Conventional" = "cornflowerblue",
                     "Strategy1"    = "tomato2",
                     "Strategy2"    = "seagreen",
                     "Strategy3"    = "gold2")

plot_theme <- theme(
  panel.grid.minor     = element_blank(),
  panel.grid.major     = element_blank(),
  plot.title           = element_text(size = 20, hjust = 0.5, face = "bold"),
  legend.title         = element_blank(),
  legend.text          = element_text(size = 14),
  legend.key           = element_blank(),
  legend.key.width     = unit(2, "cm"),
  axis.text.x          = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
  axis.text.y          = element_text(size = 14, colour = "black"),
  axis.title.x         = element_blank(),
  axis.title.y         = element_text(size = 18, face = "bold")
)

# ── 4. Yield plot ───────────────────────────────────────────────────────────

p_yield <- ggplot(t1_summary, aes(x = Year, y = mean, color = Scenario)) +
  geom_point(size = 2) +
  plot_theme +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous("Year", limits = c(19, 50)) +
  scale_y_continuous("Parent mean", expand = c(0, 0)) +
  ggtitle("Yield")

# ── 5. Specific gravity plot ───────────────────────────────────────────────

p_sg <- ggplot(t2_summary, aes(x = Year, y = mean, color = Scenario)) +
  geom_point(size = 2) +
  plot_theme +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous("Year", limits = c(19, 50)) +
  scale_y_continuous("Parent mean", expand = c(0, 0)) +
  ggtitle("Specific Gravity")

# ── 6. Combined figure ─────────────────────────────────────────────────────

combined_plot <- ggarrange(
  p_yield, p_sg,
  labels        = c("(A)", "(B)"),
  ncol          = 2,
  nrow          = 1,
  common.legend = TRUE,
  legend        = "bottom",
  align         = "hv"
)

# Save
tiff("Yield_SG_simulation_results.tiff",
     width = 16, height = 10, units = "in", res = 450)
print(combined_plot)
dev.off()

message("Plot saved: Yield_SG_simulation_results.tiff")
