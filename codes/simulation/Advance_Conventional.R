##############################################################################
# Advance_Conventional.R — Conventional phenotypic selection
#
# Standard potato breeding pipeline: phenotypic selection at all stages
# with selection index (yield + SG).
##############################################################################

# ── Year 6: Official trials 2 (30 reps → select 3) ──
Fouth_Clonal     <- setPheno(Third_Clonal_Sel, reps = 30, varE = varE)
Fouth_Clonal_Sel <- selectInd(Fouth_Clonal, 3, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 5: Official trials 1 (15 reps → select 20) ──
Third_Clonal     <- setPheno(Second_Clonal_Sel, reps = 15, varE = varE)
Third_Clonal_Sel <- selectInd(Third_Clonal, 20, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 4: Third clonal generation (6 reps → select 60) ──
Second_Clonal     <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
Second_Clonal_Sel <- selectInd(Second_Clonal, 60, use = "pheno",
                               trait = selIndex, b = weight, scale = TRUE)

# ── Year 3: Second clonal generation (2 reps → select 200) ──
First_Clonal     <- setPheno(Field_seedling_Sel, H2 = c(0.1, 0.2), reps = 2)
First_Clonal_Sel <- selectInd(First_Clonal, 200, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 2: Seedling field (select 650) ──
Field_seedling     <- setPheno(Clones, h2 = c(0.05, 0))
Field_seedling_Sel <- selectInd(Field_seedling, 650, use = "pheno")

# ── Year 1: Crossing ──
Clones <- randCross(Parents, nCrosses = 80, nProgeny = 200)
