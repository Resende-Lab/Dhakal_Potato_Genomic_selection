##############################################################################
# Genomic Selection — CV0: Predicting across environments (LOEO)
#
# Description:
#   Leave-one-environment-out (LOEO) cross-validation comparing four
#   multi-environment GBLUP models of increasing complexity:
#     M1: G                       (genomic only)
#     M2: E + G                   (environment + genomic)
#     M3: E + G + G×E             (+ genomic × environment interaction)
#     M4: E + G + G×W             (+ genomic × enviromic interaction)
#
# Inputs:
#   - GS_BLUES_all_FINAL_2.xlsx          : BLUEs across environments
#   - M_LGC_DART_4k_ADM_merged.csv       : Merged SNP dosage matrix
#   - ECData(VarEnv_raw_allDates_100).RData : Enviromic covariance (W matrix)
#
# Output:
#   - prediction_results_CV0_all_models.csv : LOEO accuracies for all models
#
# Dependencies:
#   AGHmatrix, BGLR, dplyr, tibble, ggplot2, readxl
##############################################################################

library(AGHmatrix)
library(BGLR)
library(dplyr)
library(tibble)
library(ggplot2)
library(readxl)

# ── 1. Load & prepare phenotype data ────────────────────────────────────────

GS_BLUEs <- read_excel("GS_BLUES_all_FINAL_2.xlsx",
                        sheet = "Sheet1", col_names = TRUE)

BLUEsPheno <- GS_BLUEs %>%
  select(Genotype, Env, Total_yield, Mkt_yield, SG) %>%
  mutate(Genotype = gsub("[-.]", "_", Genotype))

# Scale traits within each environment
BLUEsPheno_scaled <- BLUEsPheno %>%
  group_by(Env) %>%
  mutate(Total_yield = scale(as.numeric(Total_yield)),
         Mkt_yield   = scale(as.numeric(Mkt_yield)),
         SG          = scale(as.numeric(SG))) %>%
  ungroup()

# ── 2. Load & align SNP data ───────────────────────────────────────────────

Genotypes <- read.csv("M_LGC_DART_4k_ADM_merged.csv")

common_genotypes <- intersect(BLUEsPheno_scaled$Genotype, colnames(Genotypes))
common_sorted    <- sort(common_genotypes)

BLUEsPheno_scaled <- BLUEsPheno_scaled %>%
  filter(Genotype %in% common_sorted) %>%
  arrange(factor(Genotype, levels = common_sorted))

Genotypes <- Genotypes[, common_sorted]
snpdata   <- t(Genotypes)

# ── 3. Build genomic relationship matrix ────────────────────────────────────

G_mat <- Gmatrix(SNPmatrix    = snpdata,
                 maf           = 0.05,
                 thresh.missing = 0.2,
                 method        = "VanRaden",
                 ploidy        = 4,
                 missingValue  = "NA")

# Subset and order G matrix to match phenotype genotypes
geno_ids <- sort(unique(BLUEsPheno_scaled$Genotype))
G_mat    <- G_mat[geno_ids, geno_ids]

# ── 4. Build design matrices & kernels ──────────────────────────────────────

# Genotype incidence matrix → genomic kernel
Z_G <- model.matrix(~ 0 + Genotype, data = BLUEsPheno_scaled)
KG  <- Z_G %*% G_mat %*% t(Z_G)

# Environment incidence matrix → G×E kernel
Z_E <- model.matrix(~ 0 + Env, data = BLUEsPheno_scaled)
KGE <- KG * (Z_E %*% t(Z_E))

# Enviromic covariance → G×W kernel
load("ECData(VarEnv_raw_allDates_100).RData")
KW  <- Z_E %*% WWt[[2]] %*% t(Z_E)
KGW <- KG * KW

# ── 5. Define models ───────────────────────────────────────────────────────

models <- list(
  M1 = list(list(K = KG, model = "RKHS")),

  M2 = list(list(model = "BRR", X = Z_E),
            list(K = KG, model = "RKHS")),

  M3 = list(list(model = "BRR", X = Z_E),
            list(K = KG, model = "RKHS"),
            list(K = KGE, model = "RKHS")),

  M4 = list(list(model = "BRR", X = Z_E),
            list(K = KG, model = "RKHS"),
            list(K = KGW, model = "RKHS"))
)

model_descriptions <- c(
  M1 = "G",
  M2 = "E + G",
  M3 = "E + G + GxE",
  M4 = "E + G + GxW"
)

# ── 6. LOEO cross-validation ───────────────────────────────────────────────

envs   <- c("FL_2023", "FL_M2_2024", "TRS_2023", "BG_2023",
            "FL_M1_2024", "TRS_2024", "BG_2024")
traits <- c("Total_yield", "Mkt_yield", "SG")

all_results <- data.frame()

for (mod_name in names(models)) {
  for (test_env in envs) {
    for (trait in traits) {

      # Mask the test environment
      yNA <- BLUEsPheno_scaled[[trait]]
      yNA[BLUEsPheno_scaled$Env == test_env] <- NA

      # Fit model
      set.seed(123)
      fm <- BGLR(y       = yNA,
                 ETA     = models[[mod_name]],
                 nIter   = 30000,
                 burnIn  = 3000,
                 thin    = 10,
                 verbose = FALSE)

      # Evaluate on held-out environment
      idx_test <- BLUEsPheno_scaled$Env == test_env
      acc <- cor(fm$yHat[idx_test],
                 BLUEsPheno_scaled[[trait]][idx_test],
                 use = "complete.obs")

      all_results <- rbind(all_results,
                           data.frame(Model    = mod_name,
                                      Model_description = model_descriptions[mod_name],
                                      Env_Test = test_env,
                                      Trait    = trait,
                                      Cor      = acc))

      message(sprintf("%-3s | Env: %-12s | Trait: %-12s | r = %.3f",
                      mod_name, test_env, trait, acc))
    }
  }
}

# ── 7. Save results ────────────────────────────────────────────────────────

write.csv(all_results, "prediction_results_CV0_all_models.csv", row.names = FALSE)
message("Done. Results saved to prediction_results_CV0_all_models.csv")
