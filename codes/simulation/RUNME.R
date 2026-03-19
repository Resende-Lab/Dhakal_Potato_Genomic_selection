##############################################################################
# RUNME.R — Potato breeding simulation with AlphaSimR
#
# Simulates a 6-stage potato breeding pipeline over 50 years:
#   - 20-year burn-in (phenotypic selection)
#   - 30-year future period testing 4 scenarios:
#       Conventional : Standard phenotypic pipeline
#       Strategy1    : Genotype 650 seedlings
#       Strategy2    : Genotype 1,600 seedlings
#       Strategy3    : Genotype all F1 individuals
#
# Usage:
#   Rscript RUNME.R <rep> <ddMeanSG>
#   e.g., Rscript RUNME.R 1 0.3
#
# Dependencies: AlphaSimR, AGHmatrix
##############################################################################

require("AlphaSimR")
require("AGHmatrix")

options(echo = TRUE)
args     <- commandArgs(trailingOnly = TRUE)
rep      <- as.numeric(args[1])
ddMeanSG <- as.numeric(args[2])

# ── 1. Global parameters & base population ──────────────────────────────────

source("GlobalParameters.R")

nYears  <- burninYears + futureYears
nTraits <- 2

parentMean  <- matrix(NA, nYears, nTraits)
parentVar   <- matrix(NA, nYears, nTraits)
varietyMean  <- matrix(NA, nYears, nTraits)
Inbreeding  <- matrix(NA, nYears)
IbdDeltaF   <- matrix(NA, nYears)
heritab     <- matrix(NA, nYears, nTraits)
weight      <- c(0.5, 0.5)

# ── 2. Fill initial pipeline ────────────────────────────────────────────────

source("FillPipeline.R")

P <- runif(burninYears + futureYears)  # p-values for G×Y effects

# ── 3. Burn-in period (20 years, phenotypic selection) ──────────────────────

for (year in 1:burninYears) {
  cat("Burn-in year:", year, "\n")
  source("UpdateParents_Burnin.R")
  source("AdvanceBurnin.R")
  source("UpdateResults.R")
}

save.image(paste0("outMod/BURNIN_", rep, "_", ddMean[1], ".RData"))

# ══════════════════════════════════════════════════════════════════════════════
# FUTURE SCENARIOS — Each reloads burn-in state to start from the same point
# ══════════════════════════════════════════════════════════════════════════════

# ── 4. Conventional — Standard phenotypic selection ─────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMean[1], ".RData"))

cat("Working on: Conventional\n")
for (year in (burninYears + 1):(burninYears + futureYears)) {
  cat("  Year:", year, "\n")
  p <- P[year]
  source("UpdateParents_Scenario.R")
  source("Advance_Conventional.R")
  source("UpdateResults.R")
}

output_conv <- data.frame(
  rep          = rep(rep, nYears),
  scenario     = "Conventional",
  ParMeanT1    = parentMean[, 1],
  ParMeanT2    = parentMean[, 2],
  ParVarT1     = parentVar[, 1],
  ParVarT2     = parentVar[, 2],
  varyMeanT1   = varietyMean[, 1],
  varyMeanT2   = varietyMean[, 2],
  Inbreeding   = Inbreeding,
  IbdDeltaF    = IbdDeltaF,
  h2MeanT1     = heritab[, 1],
  h2MeanT2     = heritab[, 2],
  stringsAsFactors = FALSE
)
saveRDS(output_conv, paste0("outMod/Results_Conventional_", rep, "_", ddMean[1], ".rds"))

# ── 5. Strategy 1 — Genotype 650 seedlings ─────────────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMean[1], ".RData"))

cat("Working on: Strategy1\n")
for (year in (burninYears + 1):(burninYears + futureYears)) {
  cat("  Year:", year, "\n")
  p <- P[year]
  source("UpdateParents_Scenario.R")
  source("Advance_Strategy1.R")
  source("UpdateResults.R")
}

output_s1 <- data.frame(
  rep          = rep(rep, nYears),
  scenario     = "Strategy1",
  ParMeanT1    = parentMean[, 1],
  ParMeanT2    = parentMean[, 2],
  ParVarT1     = parentVar[, 1],
  ParVarT2     = parentVar[, 2],
  varyMeanT1   = varietyMean[, 1],
  varyMeanT2   = varietyMean[, 2],
  Inbreeding   = Inbreeding,
  IbdDeltaF    = IbdDeltaF,
  h2MeanT1     = heritab[, 1],
  h2MeanT2     = heritab[, 2],
  stringsAsFactors = FALSE
)
saveRDS(output_s1, paste0("outMod/Results_Strategy1_", rep, "_", ddMean[1], ".rds"))

# ── 6. Strategy 2 — Genotype 1,600 seedlings ───────────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMean[1], ".RData"))

cat("Working on: Strategy2\n")
for (year in (burninYears + 1):(burninYears + futureYears)) {
  cat("  Year:", year, "\n")
  p <- P[year]
  source("UpdateParents_Scenario.R")
  source("Advance_Strategy2.R")
  source("UpdateResults.R")
}

output_s2 <- data.frame(
  rep          = rep(rep, nYears),
  scenario     = "Strategy2",
  ParMeanT1    = parentMean[, 1],
  ParMeanT2    = parentMean[, 2],
  ParVarT1     = parentVar[, 1],
  ParVarT2     = parentVar[, 2],
  varyMeanT1   = varietyMean[, 1],
  varyMeanT2   = varietyMean[, 2],
  Inbreeding   = Inbreeding,
  IbdDeltaF    = IbdDeltaF,
  h2MeanT1     = heritab[, 1],
  h2MeanT2     = heritab[, 2],
  stringsAsFactors = FALSE
)
saveRDS(output_s2, paste0("outMod/Results_Strategy2_", rep, "_", ddMean[1], ".rds"))

# ── 7. Strategy 3 — Genotype all F1 individuals ────────────────────────────
# NOTE: ddMean is reset to 0 here before loading burn-in.
# This loads a different burn-in file (BURNIN_<rep>_0.RData).

ddMean <- 0
load(paste0("outMod/BURNIN_", rep, "_", ddMean[1], ".RData"))

cat("Working on: Strategy3\n")
for (year in (burninYears + 1):(burninYears + futureYears)) {
  cat("  Year:", year, "\n")
  p <- P[year]
  source("UpdateParents_Scenario.R")
  source("Advance_Strategy3.R")
  source("UpdateResults.R")
}

output_s3 <- data.frame(
  rep          = rep(rep, nYears),
  scenario     = "Strategy3",
  ParMeanT1    = parentMean[, 1],
  ParMeanT2    = parentMean[, 2],
  ParVarT1     = parentVar[, 1],
  ParVarT2     = parentVar[, 2],
  varyMeanT1   = varietyMean[, 1],
  varyMeanT2   = varietyMean[, 2],
  Inbreeding   = Inbreeding,
  IbdDeltaF    = IbdDeltaF,
  h2MeanT1     = heritab[, 1],
  h2MeanT2     = heritab[, 2],
  stringsAsFactors = FALSE
)
saveRDS(output_s3, paste0("outMod/Results_Strategy3_", rep, "_", ddMean[1], ".rds"))

cat("All scenarios complete for rep", rep, "\n")
