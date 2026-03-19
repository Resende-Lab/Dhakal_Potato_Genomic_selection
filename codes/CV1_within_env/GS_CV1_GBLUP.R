##############################################################################
# Genomic Selection — CV1: Predicting unobserved genotypes
# 
# Description:
#   5-fold cross-validation with 10 repetitions to evaluate genomic prediction
#   accuracy (GBLUP via RKHS) for potato yield and quality traits across
#   multiple environments.
#
# Inputs:
#   - GS_BLUES_all_FINAL_2.xlsx : BLUEs for Total_yield, Mkt_yield, SG by Env
#   - Dart_ADM_filtered_4K.csv  : DArT SNP dosage matrix (genotypes as columns)
#
# Output:
#   - predictions_CV1_GBLUP.csv : per-fold prediction accuracies (correlations)
#
# Dependencies:
#   AGHmatrix, BGLR, dplyr, readxl
##############################################################################

library(AGHmatrix)
library(BGLR)
library(dplyr)
library(readxl)

# ── 1. Load & prepare phenotype data ────────────────────────────────────────

GS_BLUEs <- read_excel("GS_BLUES_all_FINAL_2.xlsx",
                        sheet = "Sheet1", col_names = TRUE)

# Scale traits within each environment (zero mean, unit variance)
GS_BLUEs_scaled <- GS_BLUEs %>%
  group_by(Env) %>%
  mutate(Total_yield = scale(as.numeric(Total_yield)),
         Mkt_yield   = scale(as.numeric(Mkt_yield)),
         SG          = scale(as.numeric(SG))) %>%
  ungroup()

# ── 2. Define environments and traits ───────────────────────────────────────

selected_envs <- c("BG_2023", "FL_M1_2024", "TRS_2024", "BG_2024")
environments  <- intersect(unique(GS_BLUEs_scaled$Env), selected_envs)
traits        <- c("Total_yield", "Mkt_yield", "SG")

# ── 3. Load & format SNP data ──────────────────────────────────────────────

Genotypes <- as.matrix(read.csv("Dart_ADM_filtered_4K.csv"))
colnames(Genotypes) <- gsub("\\_", "-", colnames(Genotypes))

# ── 4. Cross-validation parameters ─────────────────────────────────────────

nReps  <- 10
nFolds <- 5
set.seed(928761)

# ── 5. Run CV1 across environments × traits ────────────────────────────────

all_predictions <- data.frame()

for (env in environments) {
  for (trait in traits) {

    # --- Subset phenotype data for this env × trait ---
    phenodata <- GS_BLUEs_scaled %>%
      filter(Env == env) %>%
      select(Genotype, all_of(trait)) %>%
      na.omit()

    # --- Align genotypes between pheno and geno ---
    common_genotypes <- intersect(phenodata$Genotype, colnames(Genotypes))
    phenodata        <- phenodata %>% filter(Genotype %in% common_genotypes)
    Genotypes_filt   <- Genotypes[, common_genotypes]

    # Sort alphabetically to ensure consistent ordering
    common_sorted    <- sort(common_genotypes)
    phenodata        <- phenodata %>%
      arrange(factor(Genotype, levels = common_sorted))
    Genotypes_filt   <- Genotypes_filt[, common_sorted]

    # --- Build genomic relationship matrix (VanRaden, tetraploid) ---
    snpdata <- t(Genotypes_filt)
    snpdata <- apply(snpdata, 2, as.numeric)

    G_vr <- Gmatrix(SNPmatrix    = snpdata,
                    maf           = 0.05,
                    thresh.missing = 0.2,
                    method        = "VanRaden",
                    ploidy        = 4,
                    missingValue  = "NA")

    Y <- as.numeric(phenodata[[trait]])

    # --- Repeated k-fold CV ---
    for (Rep in 1:nReps) {
      folds <- sample(1:nFolds, size = length(Y), replace = TRUE)

      for (fold in 1:nFolds) {
        tst     <- which(folds == fold)
        yNA     <- Y
        yNA[tst] <- NA

        # Fit GBLUP (RKHS kernel)
        fm <- BGLR(y       = yNA,
                   ETA     = list(list(K = G_vr, model = "RKHS")),
                   nIter   = 30000,
                   burnIn  = 3000,
                   thin    = 10,
                   verbose = FALSE)

        # Record prediction accuracy
        acc <- cor(Y[tst], fm$yHat[tst], use = "complete.obs")

        all_predictions <- rbind(all_predictions,
                                 data.frame(Environment = env,
                                            Trait       = trait,
                                            Repetition  = Rep,
                                            Fold        = fold,
                                            Cor         = acc))
      }
    }

    mean_cor <- mean(all_predictions$Cor[all_predictions$Environment == env &
                                          all_predictions$Trait == trait],
                     na.rm = TRUE)
    message(sprintf("Env: %-15s Trait: %-12s Mean r = %.3f", env, trait, mean_cor))
  }
}

# ── 6. Save results ────────────────────────────────────────────────────────

write.csv(all_predictions, file = "predictions_CV1_GBLUP.csv", row.names = FALSE)
message("Done. Results saved to predictions_CV1_GBLUP.csv")
