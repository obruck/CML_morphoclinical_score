rm(list = ls())

# Load necessary libraries
source("~/mounts/research/src/Rfunctions/library.R")
library(cowplot)

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")

dir.create(paste0(results, "_response"))

# Read data
selected_features = readxl::read_xlsx(paste0(results, "_response/univariate_cox_results_top_features.xlsx"))
clinical = readxl::read_xlsx(paste0(export, "_response/full_data.xlsx")) %>%
  dplyr::select(gender, age_at_dg,
                sokal_class, hasford_class, eutos_class, elts_class, spleen_size,
                b_leuk, b_trom, l_baso, l_blast, l_eos, l_lymf, bm_blast, selected_features$names
  )

# Export
writexl::write_xlsx(clinical, paste0(results, "_response/univariate_cox_results_top_features_clinical_balloonplot.xlsx"))
