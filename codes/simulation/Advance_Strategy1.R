##############################################################################
# Advance_Strategy1.R — Genotype 650 seedlings
##############################################################################

####<---------------YEAR 6: Official trials 2
Third_Clonal_Pheno <- setPheno(Third_Clonal_Sel, reps = 30, varE = varE)
h2_OT2[year,] <- diag(varG(Third_Clonal_Pheno)) / diag(varP(Third_Clonal_Pheno))
Fouth_Clonal_Sel <- selectInd(Third_Clonal_Pheno, 3, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

####<---------------YEAR 5: Official trials 1
Second_Clonal_Pheno <- setPheno(Second_Clonal_Sel, reps = 15, varE = varE)
h2_OT1[year,] <- diag(varG(Second_Clonal_Pheno)) / diag(varP(Second_Clonal_Pheno))
Third_Clonal_Sel <- selectInd(Second_Clonal_Pheno, 20, use = "pheno",
                              trait = selIndex, b = weight, scale = TRUE)

####<---------------YEAR 4: Third clonal generation
First_Clonal_Pheno <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
h2_2C[year,] <- diag(varG(First_Clonal_Pheno)) / diag(varP(First_Clonal_Pheno))
Second_Clonal_Sel <- selectInd(First_Clonal_Pheno, 60, use = "pheno",
                               trait = selIndex, b = weight, scale = TRUE)

###<---------------YEAR 2/3: Second clonal generation (seedling + GS)
Field_seedling <- setPheno(Clones, h2 = c(0.02, 0.00001))
h2_seed[year,] <- diag(varG(Field_seedling)) / diag(varP(Field_seedling))
Field_seedling_Sel <- selectInd(Field_seedling, 650, use = "pheno")

# Both traits additive -> single additive multi-trait EBV (2-column: yield, SG)
Field_seedling_Sel@ebv = setEBV_mtm(Target = Field_seedling_Sel,
                                    Training = trainPop,
                                    Dominance = FALSE,
                                    ploidy = Field_seedling_Sel@ploidy,
                                    tarTraits = 1)
Field_seedling_Sel@ebv = cbind(Field_seedling_Sel@ebv,
                               setEBV_mtm(Target = Field_seedling_Sel,
                                          Training = trainPop,
                                          Dominance = FALSE,
                                          ploidy = Field_seedling_Sel@ploidy,
                                          tarTraits = 2))

First_Clonal_Sel <- selectInd(Field_seedling_Sel, 200,  use = "ebv",
                              trait = selIndex, b = weight, scale = TRUE) 

####<----------------YEAR 1: Cross
Clones <- randCross(Parents, nCrosses= 80, nProgeny = 200)


