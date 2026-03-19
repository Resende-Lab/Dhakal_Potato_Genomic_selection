##############################################################################
# Advance_Strategy3.R — Genotype all F1 individuals
#
# Most aggressive GS strategy: all F1 clones are phenotyped at higher H2,
# skipping the seedling field stage entirely. Selection goes directly from
# F1 evaluation to clonal trials.
##############################################################################

# ── Year 6: Official trials 2 (30 reps → select 3) ──
Fouth_Clonal     <- setPheno(Third_Clonal_Sel, reps = 30, varE = varE)
Fouth_Clonal_Sel <- selectInd(Fouth_Clonal, 3, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 5: Official trials 1 (15 reps → select 20) ──
Third_Clonal     <- setPheno(Second_Clonal_Sel, reps = 15, varE = varE)
Third_Clonal_Sel <- selectInd(Third_Clonal, 20, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 3: Evaluate all F1 clones at higher H2 (select 200) ──
Clones1          <- setPheno(Clones, H2 = c(0.3, 0.5))
First_Clonal_Sel <- selectInd(Clones1, 200, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 4: Third clonal generation (6 reps → select 60) ──
Second_Clonal     <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
Second_Clonal_Sel <- selectInd(Second_Clonal, 60, use = "pheno",
                               trait = selIndex, b = weight, scale = TRUE)

# ── Year 1: Crossing ──
Clones <- randCross(Parents, nCrosses = 80, nProgeny = 200)
