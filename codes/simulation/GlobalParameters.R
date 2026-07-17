##############################################################################
# GlobalParameters.R — Trait architecture, base population, simulation settings
##############################################################################

# ── 1. Trait architecture ───────────────────────────────────────────────────
# Trait 1: Yield       Trait 2: Specific Gravity (SG)

mean = c(32,1.05)
var = c(5, 4.362632e-05)
h2=c(0.1, 0.002)
varGE = c(1.25,4.362632e-05)
varE = c(30, 4.578333e-05)

# Dominance parameters
ddVar  = c(0.2, 0)         
ddMean = c(ddMeanYLD, 0)        

# Marker/QTL counts
nSnp <- 1000
nQtl <- 300

# ── 2. Base population ──────────────────────────────────────────────────────

set.seed(9862)
founderPop <- AlphaSimR::runMacs(nInd = 2*100, nChr = 12, segSites = nSnp+nQtl, inbred = FALSE,
                                 species = "GENERIC", split = NULL, ploidy = 4L,
                                 manualCommand = NULL, manualGenLen = NULL)


# Traits
SP = AlphaSimR::SimParam$new(founderPop)

SP$restrSegSites(nQtl,nSnp)

if(nSnp>0){
  SP$addSnpChip(nSnp)
}

# ── 3. Simulation parameters (AlphaSimR) ───────────────────────────────────



# Trait 1 = YIELD: additive + directional dominance (meanDD, varDD both > 0)
SP$addTraitAG(nQtl,
               mean   = mean[1],
               var    = var[1],
               varGxE = varGE[1],
               varEnv = 0)


# Trait 2 = SPECIFIC GRAVITY: purely additive (meanDD = 0, varDD = 0)
SP$addTraitAG(nQtl,
               mean   = mean[2],
               var    = var[2],
               varGxE = varGE[2],
               varEnv = 0)


# ── 4. Create breeding population & initial parents ─────────────────────────

new_pop <- newPop(founderPop)
Parents = AlphaSimR::selectInd(new_pop, 100)
Parents = AlphaSimR::setPheno(Parents, varE = varE)

# ── 5. Year settings ────────────────────────────────────────────────────────

burninYears   <- 20
futureYears   <- 30
startTrainPop <- 15

