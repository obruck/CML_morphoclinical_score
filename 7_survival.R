# Survival plots

# Library
source("~/library.R")
library(ggpubr)
library(survminer)
library(cowplot)


# Read data
df1 = read_xlsx("./univariate_cox_results_top_features_clinical_balloonplot.xlsx")
df = read_xlsx("./full_data.xlsx") %>%
  dplyr::select(mmr, mmr_time, names(df1))
top_features = read_xlsx("./univariate_cox_results_top_features.xlsx")


# Categorize these
df = df %>%
  dplyr::mutate(lipid_total_periphery_proportion = lipid_total_periphery_proportion*100) %>%
  dplyr::mutate_at(vars(matches("(?i)(Age|Lipid|^bm_|^b_|Mono|Mega|Proe|^l_|spleen)")),
                   .funs = function(x) ntile(x, 3)) %>%
  dplyr::mutate_at(vars(matches("(?i)(Age|Lipid|^bm_|^b_|Mono|Mega|Proe|^l_|spleen)")),
                   .funs = factor) %>%
  dplyr::mutate_at(vars(matches("(?i)(Age|Lipid|^bm_|^b_|Mono|Mega|Proe|^l_|spleen)")),
                   .funs = function(x) x = ifelse(is.na(x), NA, ifelse(x==1, "Low", ifelse(x==2, "Int", "High")))) %>%
  dplyr::mutate_at(.vars = c("sokal_class", "hasford_class", "elts_class"),
                   .funs = function(x) ifelse(is.na(x), NA,
                                              ifelse(x < 2, 0, 1))) %>%
  dplyr::mutate_at(vars(matches("(?i)(gender|inib|sec|class)")),
                   .funs = function(x) x = ifelse(is.na(x), NA, ifelse(x==0, "No", "Yes")))

df = df %>%
  dplyr::rename("Gender, female (yes/no)" = "gender",
                "Age" = "age_at_dg",
                "Sokal, high (yes/no)" = "sokal_class",
                "Hasford, high (yes/no)" = "hasford_class",
                "EUTOS, high (yes/no)" = "eutos_class",
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
                                                            gsub("sec_gen_tki", "2nd gen TKI (yes/no)",
                                                                 gsub("imatinib", "Imatinib (yes/no)",
                                                                      gsub("mega", "Mega",
                                                                           gsub("mono", "Mono",
                                                                                gsub("proe", "Proe",
                                                                                     gsub("lipid", "Lipid",
                                                                                          gsub("dasatinib", "Dasatinib (yes/no)", colnames(df) ))))))))))))))))
colnames(df) = gsub("mmr time", "mmr_time", colnames(df))


# Define variable sets
top_features1 = df
colnames(top_features1) = gsub("cell ", "cell\n", colnames(top_features1))
colnames(top_features1) = gsub("nuclei ", "nuclei\n", colnames(top_features1))
colnames(top_features1) = gsub("periphery ", "periphery\n", colnames(top_features1))
colnames(top_features1) = gsub(" \\(yes\\/no\\)", "", colnames(top_features1))


# Plot
for (feature1 in names(top_features1)[3:(length(names(top_features1))-2)]) {
  
  # Survival plots
  top_features1$tmp = top_features1[[feature1]]
  top_features1$tmp = factor(top_features1$tmp, levels = unique(top_features1[order(top_features1$tmp, decreasing = TRUE),]$tmp))
  fit1 = survfit(Surv(mmr_time, mmr) ~ tmp, data = top_features1)
  
  # Plot
  g <- ggsurvplot(fit1,
                  fun = "event",
                  data=top_features1,
                  palette = "Set1",
                  size = 2.5,   #line thickness
                  ggtheme = theme_minimal(), #theme
                  font.main = c(15, "black"), #title font
                  font.x = c(15, "black"), #x-axis font
                  font.y = c(15, "black"), #y-axis font
                  font.tickslab = c(12, "black"), #axis numbering font
                  conf.int = FALSE, #confidence interval
                  pval = TRUE, #p-value
                  pval.size = 5, #p-value size
                  break.x.by = 12,
                  risk.table.fontsize  = 5,
                  pval.coord = c(48, 0.1),
                  risk.table.pos = "out",
                  tables.y.text = FALSE,
                  tables.theme = theme_cleantable(),
                  risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                  risk.table.col = "strata", #risk table color
                  risk.table.height = 0.25,
                  ylab = "MMR (%)",
                  xlab = "Time (years)",
                  surv.scale = "percent",
                  ylim = c(0,1),
                  xlim = c(0,60),
                  censor = TRUE,
                  censor.shape = 108,
                  censor.size = 4,
                  font.legend = c(14, "bold", "black"),     #font voi olla esim. "bold" tai "plain"
                  legend.title = feature1,
                  legend.labs = gsub("tmp = ", "", levels(top_features1[!is.na(top_features1$tmp),]$tmp)))
  
  g = plot_grid(g$plot + guides(colour = guide_legend(nrow = 1)), g$table +  guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
  
  # Save
  ggsave(plot = g, filename = paste0("./Survival/KM_", janitor::make_clean_names(feature1), ".png"),
         width = 5, height = 5, units = "in", dpi = 300, bg = "white") 
  
}
