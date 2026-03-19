##############################################################################
# GlobalParameters.R — Trait architecture, base population, simulation settings
##############################################################################

# ── 1. Trait architecture ───────────────────────────────────────────────────
# Trait 1: Yield       Trait 2: Specific Gravity (SG)

mean  <- c(32, 1.05)
var   <- c(5, 4.362632e-05)
h2    <- c(0.1, 0.002)
varGE <- c(26, 4.362632e-05)
varE  <- c(30, 4.578333e-05)

# Dominance parameters
ddVar  <- c(0, 0.2)
ddMean <- c(0, ddMeanSG)        # ddMeanSG passed from command line via RUNME.R

# Marker/QTL counts
nSnp <- 400
nQtl <- 300

# ── 2. Base population ──────────────────────────────────────────────────────

set.seed(9862)
founderPop <- runMacs(nInd     = 200,
                      nChr     = 12,
                      segSites = nSnp + nQtl,
                      inbred   = FALSE,
                      species  = "GENERIC",
                      ploidy   = 4L)

# ── 3. Simulation parameters (AlphaSimR) ───────────────────────────────────

SP <- SimParam$new(founderPop)
SP$restrSegSites(nQtl, nSnp)

if (nSnp > 0) {
  SP$addSnpChip(nSnp)
}

SP$addTraitADG(nQtl,
               mean   = mean,
               var    = var,
               meanDD = ddMean,
               varDD  = ddVar,
               varGxE = varGE,
               varEnv = 0)$
  setVarE(varE = varE)

# ── 4. Create breeding population & initial parents ─────────────────────────

new_pop <- newPop(founderPop)
Parents <- selectInd(new_pop, 100)
Parents <- setPheno(Parents, varE = varE)

# ── 5. Year settings ────────────────────────────────────────────────────────

burninYears   <- 20
futureYears   <- 30
startTrainPop <- 11
