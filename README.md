# Potato_GS_pipeline

Genomic selection (GS) pipeline for potato breeding: within- and across-environment prediction, genotyping platform comparison, and breeding program simulation using AlphaSimR.

---

## Repository Structure

```
Potato_GS_pipeline/
│
├── data/                            # Placeholder (data may be added after PI approval)
│
├── codes/
│   ├── CV1_within_env/              # Section 3.2 — Within-environment prediction
│   │   └── GS_CV1_GBLUP.R
│   │
│   ├── CV0_across_env/              # Section 3.3 — Across-environment prediction
│   │   └── GS_CV0_LOEO.R
│   │
│   ├── platform_comparison/         # Section 3.4 — Genotyping platform comparison
│   │   └── GS_CV1_DArT_4K.R
│   │
│   ├── simulation/                  # Section 3.5 — AlphaSimR breeding simulation
│   │   ├── RUNME.R                  #   Main orchestrator
│   │   ├── GlobalParameters.R       #   Trait architecture & base population
│   │   ├── FillPipeline.R           #   Initialize breeding pipeline
│   │   ├── AdvanceBurnin.R          #   Burn-in advancement (phenotypic)
│   │   ├── Advance_Conventional.R   #   Scenario: Conventional selection
│   │   ├── Advance_Strategy1.R      #   Scenario: Genotype 650 seedlings
│   │   ├── Advance_Strategy2.R      #   Scenario: Genotype 1,600 seedlings
│   │   ├── Advance_Strategy3.R      #   Scenario: Genotype all F1s
│   │   ├── UpdateParents_Burnin.R   #   Parent selection (burn-in)
│   │   ├── UpdateParents_Scenario.R #   Parent selection (scenarios)
│   │   ├── UpdateResults.R          #   Record metrics per year
│   │   ├── plot_results.R           #   Visualization
│   │   └── rep.txt                  #   SLURM array input
│   │
│   └── SLURM/                       # HPC job submission scripts
│       ├── run_CV1.sbatch
│       ├── run_CV0.sbatch
│       └── run_simulation.sbatch
│
├── .gitignore
└── README.md
```

---

## Analyses

### 1. CV1 — Within-Environment Prediction (Section 3.2)

Predicts **unobserved genotypes** within environments where other genotypes have been phenotyped.

- **Script**: `codes/CV1_within_env/GS_CV1_GBLUP.R`
- **Model**: GBLUP via RKHS kernel (`BGLR`)
- **G matrix**: VanRaden for autotetraploids (`AGHmatrix`, ploidy = 4)
- **Markers**: ~4K DArT SNPs (dosage)
- **CV**: 5-fold, 10 repetitions per environment × trait
- **Traits**: Total yield, marketable yield, specific gravity
- **Environments**: BG_2023, FL_M1_2024, TRS_2024, BG_2024
- **Phenotypes**: BLUEs scaled within environment (zero mean, unit variance)

### 2. CV0 — Across-Environment Prediction (Section 3.3)

Predicts genotype performance in **entirely unobserved environments** using leave-one-environment-out (LOEO).

- **Script**: `codes/CV0_across_env/GS_CV0_LOEO.R`
- **Models compared**:

| Model | Components | Description |
|-------|-----------|-------------|
| M1 | G | Genomic kernel only |
| M2 | E + G | Environment main effect + genomic |
| M3 | E + G + G×E | + Genomic × environment interaction |
| M4 | E + G + G×W | + Genomic × enviromic interaction |

- **Environments**: FL_2023, FL_M2_2024, TRS_2023, BG_2023, FL_M1_2024, TRS_2024, BG_2024
- **Enviromic data**: Weather covariance matrix from `ECData(VarEnv_raw_allDates_100).RData`

### 3. Platform Comparison (Section 3.4)

Compares prediction accuracy across genotyping platforms and marker densities.

- **Script**: `codes/platform_comparison/GS_CV1_DArT_4K.R` — CV1 with DArT 4K markers
- Additional scripts for Flex-seq 4K vs 105K and imputation to be added.

### 4. Breeding Program Simulation (Section 3.5)

Stochastic simulation of a potato breeding pipeline comparing conventional phenotypic selection against genomic selection strategies with varying genotyping intensity.

- **Scripts**: `codes/simulation/` (14 files, orchestrated by `RUNME.R`)
- **Entry point**: `Rscript RUNME.R <rep> <ddMeanSG>`

**Pipeline (6 stages per year):**

| Stage | Description | Reps | Selected |
|-------|------------|------|----------|
| Year 1 | Crossing (80 crosses × 200 progeny) | — | 16,000 |
| Year 2 | Seedling field | — | varies by scenario |
| Year 3 | 2nd clonal generation | 2 | 200 |
| Year 4 | 3rd clonal generation | 6 | 60 |
| Year 5 | Official trial 1 | 15 | 20 |
| Year 6 | Official trial 2 | 30 | 3 |

**Scenarios:**

| Scenario | Description |
|----------|-------------|
| Conventional | Standard phenotypic selection at all stages |
| Strategy 1 | Genotype 650 seedlings, higher H2 at 2nd clonal |
| Strategy 2 | Genotype 1,600 seedlings, higher H2 at 2nd clonal |
| Strategy 3 | Genotype all F1s, skip seedling stage |

**Design**: 200 tetraploid founders, 12 chromosomes, 300 QTL, 400 SNPs, 20-year burn-in + 30-year future, 30 replicates via SLURM array.

---

## Requirements

- R (≥ 4.x)
- R packages: `BGLR`, `AGHmatrix`, `AlphaSimR`, `dplyr`, `readxl`, `ggplot2`, `ggpubr`, `data.table`, `plyr`
- HPC with SLURM (scripts written for UF HiPerGator; edit SLURM headers before use)

## Usage

```bash
# Clone the repo
git clone https://github.com/<your-username>/Potato_GS_pipeline.git
cd Potato_GS_pipeline

# --- CV1 ---
cd codes/CV1_within_env
# Place input data files here, then:
sbatch ../../codes/SLURM/run_CV1.sbatch

# --- CV0 ---
cd codes/CV0_across_env
sbatch ../../codes/SLURM/run_CV0.sbatch

# --- Simulation ---
cd codes/simulation
mkdir -p outMod logs
sbatch ../../codes/SLURM/run_simulation.sbatch
# After all jobs complete:
Rscript plot_results.R
```

## Input Data (not included)

| File | Used by | Description |
|------|---------|-------------|
| `GS_BLUES_all_FINAL_2.xlsx` | CV1, CV0, Platform | BLUEs across environments |
| `Dart_ADM_filtered_4K.csv` | CV1, Platform | DArT 4K SNP dosage matrix |
| `M_LGC_DART_4k_ADM_merged.csv` | CV0 | Merged DArT + LGC SNP dosage matrix |
| `ECData(VarEnv_raw_allDates_100).RData` | CV0 | Enviromic covariance matrices |

## Citation

*TODO: Add citation / DOI once published.*

## License

*TODO: Choose license (e.g., MIT, GPL-3).*
