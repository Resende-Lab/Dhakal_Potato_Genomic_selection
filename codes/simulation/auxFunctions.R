#' Estimate Breeding Values via Multi-Trait Regression (bWGR)
#'
#' Computes estimated breeding values (EBVs) for a target population using
#' an external multivariate solver from the \pkg{bWGR} package. This scales
#' better to large datasets than AlphaSimR's internal genomic prediction
#' model, at the cost of slightly less accurate variance component
#' estimation.
#'
#' @param Target An AlphaSimR population object for which EBVs will be
#'   predicted (must have genotypes available via \code{pullSnpGeno}).
#' @param Training An AlphaSimR population object used to train the model
#'   (must have phenotypes available via \code{pheno} and genotypes via
#'   \code{pullSnpGeno}).
#' @param Dominance Logical. If \code{TRUE}, fits an additive + dominance
#'   model using \code{bWGR::mrr2X}. If \code{FALSE} (default), fits an
#'   additive-only model using \code{bWGR::mrr}.
#' @param ploidy Integer. Ploidy level of the species, used to construct
#'   the dominance incidence matrix when \code{Dominance = TRUE}. Default
#'   is \code{4} (autotetraploid).
#'
#' @details
#' For \code{Dominance = FALSE}, the additive marker effects are estimated
#' via \code{bWGR::mrr(Y, X)} and EBVs are computed as \eqn{M \%*\% b}.
#'
#' For \code{Dominance = TRUE}, a dominance incidence matrix is constructed
#' as \code{(0 < M < ploidy) * 1} (heterozygote indicator), and additive
#' plus dominance effects are jointly estimated via
#' \code{bWGR::mrr2X(Y, X1, X2)}. EBVs are then computed as
#' \eqn{M \%*\% b_1 + M_{dom} \%*\% b_2}.
#'
#' Note: marker matrices are not centered prior to fitting. This shifts
#' all EBVs by a constant but has no adverse effect on selection decisions,
#' since ranking is preserved.
#'
#' @return A matrix of estimated breeding values for individuals in
#'   \code{Target}, with one column per trait.
#'
#' @importFrom bWGR mrr mrr2X
#'
#' @examples
#' \dontrun{
#' ebv_add    <- setEBV_mtm(Target = TargetPop, Training = TrainPop)
#' ebv_addom  <- setEBV_mtm(Target = TargetPop, Training = TrainPop,
#'                           Dominance = TRUE, ploidy = 4)
#' }
#'
#' @export
setEBV_mtm <- function(Target = NULL,
                       Training = NULL,
                       Dominance = FALSE,
                       ploidy = 4,
                       tarTraits = NULL) {
  
  ## ---- Input validation -------------------------------------------------
  if (is.null(Target) || is.null(Training)) {
    stop("Both 'Target' and 'Training' populations must be provided.")
  }
  if (!is.logical(Dominance) || length(Dominance) != 1) {
    stop("'Dominance' must be a single logical value (TRUE/FALSE).")
  }
  if (!is.numeric(ploidy) || ploidy <= 0) {
    stop("'ploidy' must be a positive numeric value.")
  }
  
  ## ---- Extract training data ---------------------------------------------
  Y <- pheno(Training)
  M <- pullSnpGeno(Training)
  
  if(is.null(tarTraits)){
    if (isFALSE(Dominance)) {
      
      ## ---- Additive-only model (bWGR::mrr) ----------------------------------
      fit <- bWGR::mrr(Y = Y, X = M)
      
      ## Assign EBVs to target population.
      ## Note: markers are not centered here -- this only shifts all EBVs by a
      ## constant and does not affect selection (ranking is preserved).
      M_target  <- pullSnpGeno(Target)
      targetEbv <- M_target %*% fit$b
      
    } else {
      
      ## ---- Additive + dominance model (bWGR::mrr2X) -------------------------
      Mdom <- (M > 0 & M < ploidy) * 1
      fit  <- bWGR::mrr2X(Y = as.matrix(Y[,2]), X1 = M, X2 = Mdom)
      
      M_target    <- pullSnpGeno(Target)
      Mdom_target <- (M_target > 0 & M_target < ploidy) * 1
      targetEbv <- M_target %*% fit$b1 + Mdom_target %*% fit$b2
    }
  
     } else { 
    if (isFALSE(Dominance)) {
    
      ## ---- Additive-only model (bWGR::mrr) ----------------------------------
      yMat = as.matrix(Y[,tarTraits])
      fit <- bWGR::mrr(Y = yMat, X = M)
    
      ## Assign EBVs to target population.
      ## Note: markers are not centered here -- this only shifts all EBVs by a
      ## constant and does not affect selection (ranking is preserved).
      M_target  <- pullSnpGeno(Target)
      targetEbv <- M_target %*% fit$b
    
    } else {
    
      ## ---- Additive + dominance model (bWGR::mrr2X) -------------------------
      Mdom <- (M > 0 & M < ploidy) * 1
      yMat = as.matrix(Y[,tarTraits])
      fit  <- bWGR::mrr2X(Y = yMat, X1 = M, X2 = Mdom)
      
      M_target    <- pullSnpGeno(Target)
      Mdom_target <- (M_target > 0 & M_target < ploidy) * 1
      targetEbv <- M_target %*% fit$b1 + Mdom_target %*% fit$b2
    }
  }
  return(targetEbv)
}





