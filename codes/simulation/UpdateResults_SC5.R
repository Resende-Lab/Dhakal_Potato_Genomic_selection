# Measure the results from each cycle

#>>---- Parents
parentMean[year,] = meanG(Parents)
parentVar[year,] = diag(varG(Parents))

#>>---- Released material
varietyMean[year,] = meanG(Fouth_Clonal_Sel)

#>>---- Heritability
heritab[year,] <- diag(varG(Parents))/diag(varP(Parents))

# F1 / year-1 selection accuracy, measured on the genotyped candidates (all F1s).
# Uses Clones_GS (snapshot taken in AdvanceSC5.R) because Clones itself is
# re-crossed at the end of the advance and no longer carries EBVs.
Accuracy[year,] <- diag(cor(Clones_GS@gv, Clones_GS@ebv))

#>>---- Inbreeding


#' Calculates the inbreeding rates of a target population
#'
#' @description
#' The function implements the estimation for inbreeding rates based on the
#' proposition of Falconer and Mackey (1996). The amount of heterozygous is
#' measured as a proxy for the inbreeding rate for a target population.
#'
#' @param W Marker matrix with the SNPs coded as 0,1,2.
#'
#' @return Inbreeding rate for a target population
#'
#' @export

inbr_rate <- function(W){
  het=1-abs(W-1)
  fi=rowSums(het)/(ncol(W))
  inbreeding=1-fi
  return(inbreeding)
}


MarkersA = pullSnpGeno(Parents)
Inbreeding[year] = mean(inbr_rate(MarkersA))



#####>>>>----- 2. Based on Falconer equation

# Inbreeding falconer
inbCoef = function(Pop){
  
  Markers = pullSnpGeno(Pop)
  p = colMeans(Markers)/4
  
  f = sum(p^2)/ncol(Markers)
  
  return(f)
}

IbdDeltaF[year] = inbCoef(Parents)
