##############################################################################
# UpdateParents_Scenario.R — Select new parents during future scenarios
#
# Multi-stage parent selection: combines top individuals from Second Clonal,
# Third Clonal, and Fourth Clonal stages.
##############################################################################

Parents <- c(
  selectInd(Second_Clonal, 12,
            trait = selIndex, b = weight, scale = TRUE),
  selectInd(Third_Clonal, 6,
            trait = selIndex, b = weight, scale = TRUE),
  selectInd(Fouth_Clonal, 3,
            trait = selIndex, b = weight, scale = TRUE)
)
