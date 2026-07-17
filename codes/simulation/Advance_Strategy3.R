##############################################################################
# Advance_Strategy3.R — Genotype all F1 individuals
#
# Most aggressive GS strategy
# skipping the seedling field stage entirely. Selection goes directly from
# F1 evaluation to clonal trials.
##############################################################################

####<---------------YEAR 6: Official trials 2
Fouth_Clonal <- setPheno(Third_Clonal_Sel, reps = 30, varE = varE)
h2_OT2[year,] <- diag(varG(Fouth_Clonal)) / diag(varP(Fouth_Clonal))
Fouth_Clonal_Sel <- selectInd(Fouth_Clonal, 3, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)
Third_Clonal_Pheno = Fouth_Clonal

####<---------------YEAR 5: Official trials 1
Third_Clonal <- setPheno(Second_Clonal_Sel, reps = 15, varE = varE)
h2_OT1[year,] <- diag(varG(Third_Clonal)) / diag(varP(Third_Clonal))
Third_Clonal_Sel <- selectInd(Third_Clonal, 20, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)
Second_Clonal_Pheno = Third_Clonal

####<---------------YEAR 2/3/4: Genotype ALL F1s + GS
Clones@ebv = setEBV_mtm(Target = Clones, Training = trainPop, Dominance = FALSE,
                        ploidy = Clones@ploidy, tarTraits = 1)
Clones@ebv = cbind(Clones@ebv,
                   setEBV_mtm(Target = Clones, Training = trainPop, Dominance = FALSE,
                              ploidy = Clones@ploidy, tarTraits = 2))

# Snapshot the genotyped candidates (GV + EBV) before Clones is re-crossed below,
# so UpdateResults_GS_SC5.R can measure F1 accuracy on the correct population.
Clones_GS <- Clones

First_Clonal_Sel <- selectInd(Clones,  200, use = "ebv",
                              trait = selIndex, b = weight, scale = TRUE)

# Third clonal generation
Second_Clonal <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
h2_2C[year,] <- diag(varG(Second_Clonal)) / diag(varP(Second_Clonal))
Second_Clonal_Sel <- selectInd(Second_Clonal, 60, use = "pheno",
                               trait = selIndex, b = weight, scale = TRUE)
First_Clonal_Pheno = Second_Clonal

####<----------------YEAR 1: Cross
Clones <- randCross(Parents, nCrosses= 80, nProgeny = 200)
