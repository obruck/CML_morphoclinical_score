rm(list = ls())

# Load necessary libraries
source("~/library.R")
library(cowplot)

# Read data
selected_features = readxl::read_xlsx("./univariate_cox_results_top_features.xlsx")
clinical = readxl::read_xlsx("./full_data.xlsx") %>%
  dplyr::select(gender, age_at_dg,
                sokal_class, hasford_class, eutos_class, elts_class, spleen_size,
                b_leuk, b_trom, l_baso, l_blast, l_eos, l_lymf, bm_blast, selected_features$names
  )

# Export
writexl::write_xlsx(clinical, "/univariate_cox_results_top_features_clinical_balloonplot.xlsx")
