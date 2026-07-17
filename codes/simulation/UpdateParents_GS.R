
### ---------------
First_Clonal_Sel@ebv = setEBV_mtm(Target = First_Clonal_Sel,
                                    Training = trainPop,
                                  Dominance = FALSE,
                                  ploidy = First_Clonal_Sel@ploidy,
                                  tarTraits = NULL)


AccPar[year,] <- diag(cor(First_Clonal_Sel@gv, First_Clonal_Sel@ebv))   

# Update of the parents at each cycle
Parents = c(selectInd(First_Clonal_Sel, 50, use = "ebv",
                      trait = selIndex, b = weight, scale = TRUE))
)
