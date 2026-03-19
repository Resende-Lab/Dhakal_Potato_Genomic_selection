##############################################################################
# UpdateParents_Burnin.R — Select new parents during burn-in
#
# Selects top 50 from Third Clonal stage based on selection index.
##############################################################################

Parents <- selectInd(Second_Clonal, 50,
                     trait = selIndex, b = weight, scale = TRUE)
