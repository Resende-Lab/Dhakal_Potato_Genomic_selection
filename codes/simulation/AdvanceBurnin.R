##############################################################################
# AdvanceBurnin.R — Advance pipeline one year during burn-in
#
# Phenotypic selection at all stages. Works backward through the pipeline
# to avoid overwriting data needed by later stages.
##############################################################################

# ── Year 6: Official trials 2 ──
Fouth_Clonal     <- setPheno(Third_Clonal_Sel, reps = 30, varE = varE)
Fouth_Clonal_Sel <- selectInd(Fouth_Clonal, 3, use = "pheno")

# ── Year 5: Official trials 1 ──
Third_Clonal     <- setPheno(Second_Clonal_Sel, reps = 15, varE = varE)
Third_Clonal_Sel <- selectInd(Third_Clonal, 20, use = "pheno")

# ── Year 4: Third clonal generation ──
Second_Clonal     <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
Second_Clonal_Sel <- selectInd(Second_Clonal, 60)

# ── Year 3: Second clonal generation ──
First_Clonal     <- setPheno(Field_seedling_Sel, reps = 2, varE = varE)
First_Clonal_Sel <- selectInd(First_Clonal, 200, use = "pheno")

# ── Year 2: Seedling field ──
Field_seedling     <- setPheno(Clones, h2 = c(0.05, 0.5))
Field_seedling_Sel <- selectInd(Field_seedling, 650, use = "pheno")

# ── Year 1: Crossing ──
Clones <- randCross(Parents, nCrosses = 80, nProgeny = 200)
