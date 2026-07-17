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

source("auxFunctions.R")

options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
rep       <- as.numeric(args[1])
ddMeanYLD <- as.numeric(args[2])

# ── 1. Global parameters & base population ──────────────────────────────────

source("GlobalParameters.R")

# Matrices to record results (one row per year, one column per trait)
nYears  = burninYears + futureYears
nTraits = 2
parentMean  = matrix(NA, nYears, nTraits)
parentVar   = matrix(NA, nYears, nTraits)
varietyMean = matrix(NA, nYears, nTraits)
Inbreeding  = matrix(NA, nYears)
IbdDeltaF   = matrix(NA, nYears)
heritab     = matrix(NA, nYears, nTraits)
h2_seed = matrix(NA, nYears, nTraits)
h2_1C   = matrix(NA, nYears, nTraits)
h2_2C   = matrix(NA, nYears, nTraits)
h2_OT1  = matrix(NA, nYears, nTraits)
h2_OT2  = matrix(NA, nYears, nTraits)
Accuracy    = matrix(NA, nYears, nTraits)   # accuracy at the F1 / year-1 selection step
AccPar      = matrix(NA, nYears, nTraits)   # accuracy at the parent-update step
weight      = matrix(c(0.5,0.5), ncol = 1)

# ── 2. Fill initial pipeline ────────────────────────────────────────────────

source("FillPipeline.R")

P = runif(nYears)  # p-values for G×Y effects

# ── 3. Burn-in period (20 years, phenotypic selection) ──────────────────────

for(year in 1:burninYears){
  cat("Burn-in year:", year, "\n")
  source("UpdateParents.R")    # pick parents from last year's data
  source("AdvanceSC1.R")       # advance the pipeline one year
  source("UpdateResults.R")    # record metrics
  source("WriteRecords_GS.R")  # build the GS training population
}

# Save burn-in state so every scenario starts from the same point
save.image(paste0("outMod/BURNIN_", rep, "_", ddMeanYLD, ".RData"))

# ══════════════════════════════════════════════════════════════════════════════
# FUTURE SCENARIOS — Each reloads burn-in state to start from the same point
# ══════════════════════════════════════════════════════════════════════════════

# ── 4. Conventional — Standard phenotypic selection ─────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMeanYLD, ".RData"))

cat("Working on Conventional (SC2)\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  cat("Working on year:", year, "\n")
  p = P[year]
  source("UpdateParents_SC2.R")
  source("AdvanceSC2.R")
  source("UpdateResults.R")
}

output2 = data.frame(rep=rep(rep, nYears),
                     scenario=rep("Conv", nYears),
                     ParMeanT1 = parentMean[,1],
                     ParMeanT2 = parentMean[,2],
                     ParVarT1  = parentVar[,1],
                     ParVarT2  = parentVar[,2],
                     varyMeanT1= varietyMean[,1],
                     varyMeanT2= varietyMean[,2],
                     Inbreeding,
                     IbdDeltaF,
                     h2MeanT1  = heritab[,1],
                     h2MeanT2  = heritab[,2],
                     AccuracyT1= Accuracy[,1],
                     AccuracyT2= Accuracy[,2],
                     AccParT1  = AccPar[,1],
                     AccParT2  = AccPar[,2],
                     h2_seed_T1 = h2_seed[,1], h2_seed_T2 = h2_seed[,2],
                     h2_1C_T1   = h2_1C[,1],   h2_1C_T2   = h2_1C[,2],
                     h2_2C_T1   = h2_2C[,1],   h2_2C_T2   = h2_2C[,2],
                     h2_OT1_T1  = h2_OT1[,1],  h2_OT1_T2  = h2_OT1[,2],
                     h2_OT2_T1  = h2_OT2[,1],  h2_OT2_T2  = h2_OT2[,2],
                     stringsAsFactors=FALSE)

saveRDS(output2, paste0("outMod/Results_Conventional_", rep, "_", ddMeanYLD, ".rds"))

# ── 5. Strategy 1 — Genotype 650 seedlings ─────────────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMeanYLD, ".RData"))

source("auxFunctions.R")

cat("Working on Strategy 1 (SC3)\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  if (year == burninYears+1){
    
    ####<---------------- YEAR 1: Cross
    Clones <- randCross(Parents, nCrosses= 80, nProgeny = 200, ignoreSexes = TRUE)
    
    ###<--------------- YEAR 2/3: seedlings + GS
    Field_seedling <- setPheno(Clones, h2 = c(0.05, 0.000001))
    Field_seedling_Sel <- selectInd(Field_seedling, 650, use = "pheno")
    
    Field_seedling_Sel@ebv = setEBV_mtm(Target = Field_seedling_Sel,
                                        Training = trainPop,
                                        Dominance = FALSE,     # yield = additive + dominance
                                        ploidy = Field_seedling_Sel@ploidy,
                                        tarTraits = NULL)
    
    First_Clonal_Sel <- selectInd(Field_seedling_Sel, 200, use = "ebv",
                                  trait = selIndex, b = weight, scale = TRUE)
  }
  
  cat("Working on year:", year, "\n")
  p = P[year]
  source("UpdateParents_GS.R")   # picks parents on EBV; also records AccPar
  source("AdvanceSC3.R")
  source("UpdateResults_GS.R")   # F1 accuracy on Field_seedling_Sel (correct for S1)
  source("WriteRecords_GS.R")
}

output3 = data.frame(rep=rep(rep, nYears),
                     scenario=rep("Strategy1", nYears),
                     ParMeanT1 = parentMean[,1],
                     ParMeanT2 = parentMean[,2],
                     ParVarT1  = parentVar[,1],
                     ParVarT2  = parentVar[,2],
                     varyMeanT1= varietyMean[,1],
                     varyMeanT2= varietyMean[,2],
                     Inbreeding,
                     IbdDeltaF,
                     h2MeanT1  = heritab[,1],
                     h2MeanT2  = heritab[,2],
                     AccuracyT1= Accuracy[,1],
                     AccuracyT2= Accuracy[,2],
                     AccParT1  = AccPar[,1],
                     AccParT2  = AccPar[,2],
                     h2_seed_T1 = h2_seed[,1], h2_seed_T2 = h2_seed[,2],
                     h2_1C_T1   = h2_1C[,1],   h2_1C_T2   = h2_1C[,2],
                     h2_2C_T1   = h2_2C[,1],   h2_2C_T2   = h2_2C[,2],
                     h2_OT1_T1  = h2_OT1[,1],  h2_OT1_T2  = h2_OT1[,2],
                     h2_OT2_T1  = h2_OT2[,1],  h2_OT2_T2  = h2_OT2[,2],
                     stringsAsFactors=FALSE)

saveRDS(output3, paste0("outMod/Results_Strategy1_", rep, "_", ddMeanYLD, ".rds"))

# ── 6. Strategy 2 — Genotype 1,600 seedlings ───────────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMeanYLD, ".RData"))

source("auxFunctions.R")

cat("Working on Strategy 2 (SC4)\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  if (year == burninYears+1){
    
    ####<---------------- YEAR 1: Cross
    Clones <- randCross(Parents, nCrosses= 80, nProgeny = 200, ignoreSexes = TRUE)
    
    ###<--------------- YEAR 2/3: seedlings + GS
    Field_seedling <- setPheno(Clones, h2 = c(0.05, 0.000001))
    Field_seedling_Sel <- selectInd(Field_seedling, 1600, use = "pheno")
    
    Field_seedling_Sel@ebv = setEBV_mtm(Target = Field_seedling_Sel,
                                        Training = trainPop,
                                        Dominance = FALSE,     # yield = additive + dominance
                                        ploidy = Field_seedling_Sel@ploidy,
                                        tarTraits = NULL)
    
    First_Clonal_Sel <- selectInd(Field_seedling_Sel, 200, use = "ebv",
                                  trait = selIndex, b = weight, scale = TRUE)
  }
  cat("Working on year:", year, "\n")
  p = P[year]
  source("UpdateParents_GS.R")
  source("AdvanceSC4.R")
  source("UpdateResults_GS.R")   # FIXED: S2 genotypes Field_seedling_Sel (was UpdateResults_GS_SC5.R)
  source("WriteRecords_GS.R")
}

output3 = data.frame(rep=rep(rep, nYears),
                     scenario=rep("Strategy2", nYears),
                     ParMeanT1 = parentMean[,1],
                     ParMeanT2 = parentMean[,2],
                     ParVarT1  = parentVar[,1],
                     ParVarT2  = parentVar[,2],
                     varyMeanT1= varietyMean[,1],
                     varyMeanT2= varietyMean[,2],
                     Inbreeding,
                     IbdDeltaF,
                     h2MeanT1  = heritab[,1],
                     h2MeanT2  = heritab[,2],
                     AccuracyT1= Accuracy[,1],
                     AccuracyT2= Accuracy[,2],
                     AccParT1  = AccPar[,1],
                     AccParT2  = AccPar[,2],
                     h2_seed_T1 = h2_seed[,1], h2_seed_T2 = h2_seed[,2],
                     h2_1C_T1   = h2_1C[,1],   h2_1C_T2   = h2_1C[,2],
                     h2_2C_T1   = h2_2C[,1],   h2_2C_T2   = h2_2C[,2],
                     h2_OT1_T1  = h2_OT1[,1],  h2_OT1_T2  = h2_OT1[,2],
                     h2_OT2_T1  = h2_OT2[,1],  h2_OT2_T2  = h2_OT2[,2],
                     stringsAsFactors=FALSE)

saveRDS(output3, paste0("outMod/Results_Strategy2_", rep, "_", ddMeanYLD, ".rds"))

# ── 7. Strategy 3 — Genotype all F1 individuals ────────────────────────────

load(paste0("outMod/BURNIN_", rep, "_", ddMeanYLD, ".RData"))

source("auxFunctions.R")

cat("Working on Strategy 3 (SC5)\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  if (year == burninYears+1){
    
    ####<---------------- YEAR 1: Cross
    Clones <- randCross(Parents, nCrosses= 80, nProgeny = 200, ignoreSexes = TRUE)
    
    Clones@ebv = setEBV_mtm(Target = Clones,
                            Training = trainPop,
                            Dominance = FALSE,     # yield = additive + dominance
                            ploidy = Clones@ploidy,
                            tarTraits = NULL)
    
    # GS: select directly on EBV among all genotyped F1s
    First_Clonal_Sel <- selectInd(Clones, 200, use = "ebv",
                                  trait = selIndex, b = weight, scale = TRUE)
    
    # Third clonal generation
    Second_Clonal <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
    Second_Clonal_Sel <- selectInd(Second_Clonal, 60, use = "pheno",
                                   trait = selIndex, b = weight, scale = TRUE)
  }
  cat("Working on year:", year, "\n")
  p = P[year]
  source("UpdateParents_GS.R")
  source("AdvanceSC5.R")
  source("UpdateResults_GS_SC5.R")   # FIXED: S3 genotypes Clones (was UpdateResults_GS.R)
  source("WriteRecords_GS.R")
}

output3 = data.frame(rep=rep(rep, nYears),
                     scenario=rep("Strategy3", nYears),
                     ParMeanT1 = parentMean[,1],
                     ParMeanT2 = parentMean[,2],
                     ParVarT1  = parentVar[,1],
                     ParVarT2  = parentVar[,2],
                     varyMeanT1= varietyMean[,1],
                     varyMeanT2= varietyMean[,2],
                     Inbreeding,
                     IbdDeltaF,
                     h2MeanT1  = heritab[,1],
                     h2MeanT2  = heritab[,2],
                     AccuracyT1= Accuracy[,1],
                     AccuracyT2= Accuracy[,2],
                     AccParT1  = AccPar[,1],
                     AccParT2  = AccPar[,2],
                     h2_seed_T1 = h2_seed[,1], h2_seed_T2 = h2_seed[,2],
                     h2_1C_T1   = h2_1C[,1],   h2_1C_T2   = h2_1C[,2],
                     h2_2C_T1   = h2_2C[,1],   h2_2C_T2   = h2_2C[,2],
                     h2_OT1_T1  = h2_OT1[,1],  h2_OT1_T2  = h2_OT1[,2],
                     h2_OT2_T1  = h2_OT2[,1],  h2_OT2_T2  = h2_OT2[,2],
                     stringsAsFactors=FALSE)

saveRDS(output3, paste0("outMod/Results_Strategy3_", rep, "_", ddMeanYLD, ".rds"))

cat("All scenarios complete for rep", rep, "\n")
