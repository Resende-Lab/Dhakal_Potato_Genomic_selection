##############################################################################
# GS_CV1_DArT_4K.R — CV1 with DArT 4K markers (platform comparison)
#
# Same CV1 framework as 01_CV1_within_env but using only the DArT 4K
# marker set. Part of the platform comparison (Section 3.4).
#
# Inputs:
#   - GS_BLUES_all_FINAL_2.xlsx  : BLUEs across environments
#   - Dart_ADM_filtered_4K.csv   : DArT 4K SNP dosage matrix
#
# Output:
#   - predictions_CV1_DArT_4K.csv : Per-fold prediction accuracies
##############################################################################

library(dplyr)
library(readxl)
library(AGHmatrix)
library(BGLR)

# ── 1. Load & scale phenotype data ─────────────────────────────────────────

GS_BLUEs <- read_excel("GS_BLUES_all_FINAL_2.xlsx",
                        sheet = "Sheet1", col_names = TRUE)

GS_BLUEs_scaled <- GS_BLUEs %>%
  group_by(Env) %>%
  mutate(Total_yield = scale(as.numeric(Total_yield)),
         Mkt_yield   = scale(as.numeric(Mkt_yield)),
         SG          = scale(as.numeric(SG))) %>%
  ungroup()

# ── 2. Define environments and traits ──────────────────────────────────────

selected_envs <- c("TRS_2023", "BG_2023", "FL_M1_2024", "TRS_2024", "BG_2024")
environments  <- intersect(selected_envs, unique(GS_BLUEs_scaled$Env))
traits        <- c("Total_yield", "Mkt_yield", "SG")

# ── 3. Load SNP data (DArT 4K) ────────────────────────────────────────────

Genotypes <- as.matrix(read.csv("Dart_ADM_filtered_4K.csv"))
colnames(Genotypes) <- gsub("_", "-", colnames(Genotypes))

# ── 4. Cross-validation parameters ────────────────────────────────────────

nReps  <- 10
nFolds <- 5
set.seed(928761)

# ── 5. CV1 loop ───────────────────────────────────────────────────────────

all_preds <- data.frame()

for (env in environments) {
  for (trait in traits) {

    # Subset phenotype
    phenodata <- GS_BLUEs_scaled %>%
      filter(Env == env) %>%
      select(Genotype, all_of(trait)) %>%
      na.omit()

    # Align genotypes
    common_genotypes <- intersect(phenodata$Genotype, colnames(Genotypes))
    phenodata        <- phenodata %>% filter(Genotype %in% common_genotypes)
    Genotypes_filt   <- Genotypes[, common_genotypes]

    common_sorted  <- sort(common_genotypes)
    phenodata      <- phenodata %>% arrange(factor(Genotype, levels = common_sorted))
    Genotypes_filt <- Genotypes_filt[, common_sorted]

    # G matrix
    snpdata <- apply(t(Genotypes_filt), 2, as.numeric)
    G_vr <- Gmatrix(SNPmatrix    = snpdata,
                    maf           = 0.05,
                    thresh.missing = 0.2,
                    method        = "VanRaden",
                    ploidy        = 4,
                    missingValue  = "NA")

    Y <- as.numeric(phenodata[[trait]])

    # Repeated k-fold CV
    for (Rep in 1:nReps) {
      folds <- sample(1:nFolds, size = length(Y), replace = TRUE)

      for (fold in 1:nFolds) {
        tst      <- which(folds == fold)
        yNA      <- Y
        yNA[tst] <- NA

        fm <- BGLR(y       = yNA,
                   ETA     = list(list(K = G_vr, model = "RKHS")),
                   nIter   = 30000,
                   burnIn  = 3000,
                   thin    = 10,
                   verbose = FALSE)

        all_preds <- rbind(all_preds,
                           data.frame(Environment = env,
                                      Trait       = trait,
                                      Repetition  = Rep,
                                      Fold        = fold,
                                      Cor         = cor(Y[tst], fm$yHat[tst],
                                                        use = "complete.obs")))
      }
    }

    message(sprintf("Env: %-15s Trait: %-12s Mean r = %.3f",
                    env, trait,
                    mean(all_preds$Cor[all_preds$Environment == env &
                                       all_preds$Trait == trait], na.rm = TRUE)))
  }
}

# ── 6. Save ────────────────────────────────────────────────────────────────

write.csv(all_preds, "predictions_CV1_DArT_4K.csv", row.names = FALSE)
message("Done. Results saved to predictions_CV1_DArT_4K.csv")
