##############################################################################
# AdvanceBurnin.R — Advance pipeline one year during burn-in
#
# Phenotypic selection at all stages. Works backward through the pipeline
# to avoid overwriting data needed by later stages.
##############################################################################

# ── Year 6: Official trials 2 ──
Fouth_Clonal <- setPheno(Third_Clonal_Sel, reps = 30, varE = varE)
h2_OT2[year,] <- diag(varG(Fouth_Clonal)) / diag(varP(Fouth_Clonal))
Fouth_Clonal_Sel <- selectInd(Fouth_Clonal, 3, use = "pheno")

# ── Year 5: Official trials 1 ──
Third_Clonal <- setPheno(Second_Clonal_Sel, reps = 15, varE = varE)
h2_OT1[year,] <- diag(varG(Third_Clonal)) / diag(varP(Third_Clonal))
Third_Clonal_Sel <- selectInd(Third_Clonal, 20, use = "pheno")
Third_Clonal_Pheno = Third_Clonal_Sel

# ── Year 4: Third clonal generation ──
Second_Clonal <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
h2_2C[year,] <- diag(varG(Second_Clonal)) / diag(varP(Second_Clonal))
Second_Clonal_Sel <- selectInd(Second_Clonal, 60)
Second_Clonal_Pheno = Second_Clonal_Sel

# ── Year 3: Second clonal generation ──
First_Clonal <- setPheno(Field_seedling_Sel, reps = 2, varE = varE)
h2_1C[year,] <- diag(varG(First_Clonal)) / diag(varP(First_Clonal))
First_Clonal_Sel <- selectInd(First_Clonal, 200, use = "pheno")
First_Clonal_Pheno = First_Clonal_Sel

# ── Year 2: Seedling field ──
Field_seedling <- setPheno(Clones, h2 = c(0.02, 0.00001))
h2_seed[year,] <- diag(varG(Field_seedling)) / diag(varP(Field_seedling))
Field_seedling_Sel <- selectInd(Field_seedling, 650, use = "pheno")

# ── Year 1: Crossing ──
Clones <- randCross(Parents, nCrosses = 80, nProgeny = 200)
