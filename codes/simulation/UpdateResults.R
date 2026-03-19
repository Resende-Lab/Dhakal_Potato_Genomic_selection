##############################################################################
# UpdateResults.R — Record metrics for each breeding cycle year
##############################################################################

# ── Parent genetic mean & variance ──
parentMean[year, ] <- meanG(Parents)
parentVar[year, ]  <- diag(varG(Parents))

# ── Released variety genetic mean ──
varietyMean[year, ] <- meanG(Fouth_Clonal_Sel)

# ── Heritability (narrow-sense proxy) ──
heritab[year, ] <- diag(varG(Parents)) / diag(varP(Parents))

# ── Inbreeding rate (heterozygosity-based, Falconer & Mackay 1996) ──
inbr_rate <- function(W) {
  het <- 1 - abs(W - 1)
  fi  <- rowSums(het) / ncol(W)
  return(1 - fi)
}

MarkersA <- pullSnpGeno(Parents)
Inbreeding[year] <- mean(inbr_rate(MarkersA))

# ── Inbreeding coefficient (allele-frequency-based) ──
inbCoef <- function(Pop) {
  Markers <- pullSnpGeno(Pop)
  p <- colMeans(Markers) / 2
  return(sum(p^2) / ncol(Markers))
}

IbdDeltaF[year] <- inbCoef(Parents)
