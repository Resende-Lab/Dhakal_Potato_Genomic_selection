##############################################################################
# Advance_Strategy1.R — Genotype 650 seedlings
#
# Modified pipeline: seedling field selection (650) followed by phenotyping
# at higher H2 in the second clonal stage. Skips the conventional Year 3
# order — seedling field is evaluated first, then clonal trials proceed.
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

# ── Year 3: Seedling field (select 650 phenotypically) ──
Field_seedling     <- setPheno(Clones, h2 = c(0.05, 0))
Field_seedling_Sel <- selectInd(Field_seedling, 650, use = "pheno")

# ── Year 2: Second clonal with higher H2 (2 reps → select 200) ──
First_Clonal     <- setPheno(Field_seedling_Sel, H2 = c(0.3, 0.5), reps = 2)
First_Clonal_Sel <- selectInd(First_Clonal, 200, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

# ── Year 1: Crossing ──
Clones <- randCross(Parents, nCrosses = 80, nProgeny = 200)
