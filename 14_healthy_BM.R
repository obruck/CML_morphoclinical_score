rm(list = ls())

# Load necessary libraries
source("~/mounts/research/src/Rfunctions/library.R")

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")


# Healthy data
healthy_sample_statistics = readxl::read_xlsx("mounts/research/husdatalake/disease/processed_data/CML/CML_HUS_healthy_hematoscope_100x_sample_statistics_25-08-2025.xlsx") %>%
  dplyr::mutate(album_id = as.integer(album_id)) %>%
  dplyr::select(database, album_id, contains("monocyte", ignore.case = TRUE))


healthy_differential_counts = readxl::read_xlsx("mounts/research/husdatalake/disease/processed_data/CML/CML_HUS_healthy_hematoscope_100x_differential_counts_25-08-2025.xlsx") %>%
  dplyr::mutate(album_id = as.integer(album_id)) %>%
  dplyr::select(database, album_id, contains("monocyte", ignore.case = TRUE))


# Join data
healthy = full_join(healthy_sample_statistics, healthy_differential_counts)

rm(healthy_sample_statistics)
rm(healthy_differential_counts)

# Read data
df = read_xlsx(paste0(results, "_response/univariate_cox_results_top_features_clinical_balloonplot.xlsx"))
top_features = read_xlsx(paste0(results, "_response/univariate_cox_results_top_features.xlsx"))

# Read model coefficients
coefficients_df = read_xlsx(paste0(results, "_response/coefficients.xlsx"))
coefficients = coefficients_df$coef
names(coefficients) = coefficients_df$covariates

# Cutpoint
cutpoint <- median(df$monocytes_nuclei_perimeter_median, na.rm = TRUE)

# Divide patients by monocyte nuclei size
df$group = ifelse(df$monocytes_nuclei_perimeter_median > cutpoint, "High", "Low")
df$group = factor(df$group, levels = c("Low", "High"))

df = df %>%
  dplyr::rename("Gender, female (yes/no)" = "gender",
                "Age" = "age_at_dg",
                "ELTS, high (yes/no)" = "elts_class",
                "Spleen size (cm)" = "spleen_size",
                "PB PLT (E9/L)" = "b_trom",
                "PB Baso (%)" = "l_baso",
                "PB Blast (%)" = "l_blast",
                "PB Eos (%)" = "l_eos",
                "PB Lymph (%)" = "l_lymf",
                "BM Blast (%)" = "bm_blast")
colnames(df) = gsub("_", " ",
                    gsub("_percentage", " (%)",
                         gsub("b_leuk", "PB WBC (E9/L)",
                              gsub("_mean", " (mean)",
                                   gsub("_percentile_25", " (25%)",
                                        gsub("_percentile_75", " (75%)",
                                             gsub("_median", " (median)",
                                                  gsub("_total", "",
                                                       gsub("_proportion", " (%/Area)",
                                                            gsub("sec_gen_tki", "2GTKI (yes/no)",
                                                                 gsub("imatinib", "Imatinib (yes/no)",
                                                                      gsub("mega", "Mega",
                                                                           gsub("mono", "Mono",
                                                                                gsub("proe", "Proe",
                                                                                     gsub("lipid", "Lipid",
                                                                                          gsub("dasatinib", "Dasatinib (yes/no)", colnames(df) ))))))))))))))))




# Apply the same column name transformation to healthy data
colnames(healthy) = gsub("_", " ",
                         gsub("_percentage", " (%)",
                              gsub("b_leuk", "PB WBC (E9/L)",
                                   gsub("_mean", " (mean)",
                                        gsub("_percentile_25", " (25%)",
                                             gsub("_percentile_75", " (75%)",
                                                  gsub("_median", " (median)",
                                                       gsub("_total", "",
                                                            gsub("_proportion", " (%/Area)",
                                                                 gsub("sec_gen_tki", "2nd gen TKI (yes/no)",
                                                                      gsub("imatinib", "Imatinib (yes/no)",
                                                                           gsub("mega", "Mega",
                                                                                gsub("mono", "Mono",
                                                                                     gsub("proe", "Proe",
                                                                                          gsub("lipid", "Lipid",
                                                                                               gsub("dasatinib", "Dasatinib (yes/no)", colnames(healthy) ))))))))))))))))

# Convert "nuc" to "nuclei" in healthy data to match df
colnames(healthy) = gsub("nuc", "nuclei", colnames(healthy))

mono_cols = grep("Monocytes nuclei perimeter", names(healthy), value = TRUE)
for (cn in mono_cols) {
  if (is.numeric(healthy[[cn]])) {
    healthy[[cn]] = healthy[[cn]] / 10
  }
}

# Morphological feature names
clin <- c("group")
top_features <- c("Monocytes nuclei perimeter (median)")


# Combine data 
healthy$group <- "Healthy"
df_subset = df %>% dplyr::select("Monocytes nuclei perimeter (median)", group)
healthy_subset = healthy %>% dplyr::select("Monocytes nuclei perimeter (median)", group)

combined_data = rbind(df_subset, healthy_subset)


# Remove rows with NA values
df_clean <- na.omit(combined_data)

# Run pairwise Wilcoxon
pw <- pairwise.wilcox.test(df_clean$`Monocytes nuclei perimeter (median)`, df_clean$group,
                           p.adjust.method = "BH")

# Extract results
pvals <- pw$p.value
pvals_df <- as.data.frame(as.table(pvals))
pvals_df <- na.omit(pvals_df)
colnames(pvals_df) <- c("Group1", "Group2", "p.value")


# Set group order for plotting: Healthy, Low, High
df_clean$group <- factor(df_clean$group, levels = c("Healthy", "Low", "High"))

# Helper to fetch p-value regardless of order in pvals_df
get_p <- function(a, b, tbl) {
  hit <- dplyr::filter(tbl, (Group1 == a & Group2 == b) | (Group1 == b & Group2 == a))
  if (nrow(hit) == 0) return(NA_real_) else return(hit$p.value[1])
}

pv_hl <- get_p("Healthy", "Low", pvals_df)
pv_hh <- get_p("Healthy", "High", pvals_df)

# Build annotation rows explicitly to guarantee two brackets
y_top <- suppressWarnings(max(df_clean$`Monocytes nuclei perimeter (median)`, na.rm = TRUE))
ann <- data.frame(
  group1 = c("Healthy", "Healthy"),
  group2 = c("Low", "High"),
  y.position = y_top + c(0.03, 0.08) * (ifelse(y_top == 0, 1, y_top)),
  label = c(paste0("Wilcoxon, p = ", format.pval(pv_hl, digits = 3)),
            paste0("Wilcoxon, p = ", format.pval(pv_hh, digits = 3)))
)


g <- ggplot(df_clean, aes(x = group, y = `Monocytes nuclei perimeter (median)`)) +
  geom_jitter(size=3, width = 0.2, aes(fill=group), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x = "Group", y = "Monocytes nuclei perimeter (median)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  ggpubr::stat_pvalue_manual(ann, label = "label", y.position = "y.position",
                             tip.length = 0, bracket.size = 0.7, size = 5.5,
                             inherit.aes = FALSE); g


# Save the plot
ggsave(plot = g, filename = paste0(results, "_response/Scatterplots/new/", "monocyte_nuclei_perimeter", "_with_healthy_updated.png"), width = 5, height = 5, units = "in", dpi = 300)

