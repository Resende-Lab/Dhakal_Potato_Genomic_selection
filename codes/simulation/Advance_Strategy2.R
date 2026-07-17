##############################################################################
# Advance_Strategy2.R — Genotype 1,600 seedlings
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

####<---------------YEAR 4: Third clonal generation
Second_Clonal <- setPheno(First_Clonal_Sel, reps = 6, varE = varE)
h2_2C[year,] <- diag(varG(Second_Clonal)) / diag(varP(Second_Clonal))
Second_Clonal_Sel <- selectInd(Second_Clonal, 60, use = "pheno",
                               trait = selIndex, b = weight, scale = TRUE)

First_Clonal_Pheno = Second_Clonal

####<---------------YEAR 2/3: First Clonal (seedling + GS)
Field_seedling <- setPheno(Clones, h2 = c(0.05, 0.00001))
h2_seed[year,] <- diag(varG(Field_seedling)) / diag(varP(Field_seedling))
Field_seedling_Sel <- selectInd(Field_seedling, 1600 , use = "pheno")

# Both traits additive -> single additive multi-trait EBV (2-column: yield, SG)
# Field_seedling_Sel@ebv = setEBV_mtm(Target = Field_seedling_Sel,
#                                     Training = trainPop,
#                                     Dominance = FALSE,
#                                     ploidy = Field_seedling_Sel@ploidy,
#                                     tarTraits = NULL)

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

First_Clonal_Sel <- selectInd(Field_seedling_Sel, 200, use = "ebv",
                              trait = selIndex, b = weight, scale = TRUE) #*** Parents from here plus higher h2

####<----------------YEAR 1: Cross
Clones <- randCross(Parents, nCrosses= 80, nProgeny = 200)
