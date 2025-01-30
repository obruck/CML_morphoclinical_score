# Balloonplot

# Library
source("~/mounts/research/src/Rfunctions/library.R")
library(ggpubr)

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")

dir.create(paste0(results, "_response/Scatterplots"))


# Read data
df = read_xlsx(paste0(results, "_response/univariate_cox_results_top_features_clinical_balloonplot.xlsx"))
top_features = read_xlsx(paste0(results, "_response/univariate_cox_results_top_features.xlsx"))


# Categorize these
df = df %>%
  dplyr::mutate_at(.vars = c("age_at_dg", "spleen_size", "b_trom", "l_baso", "l_blast", "l_eos", "l_lymf", "bm_blast"),
                   .funs = function(x) ifelse(x > median(x, na.rm=TRUE), 1, 0)) %>%
  dplyr::mutate_at(.vars = c("sokal_class", "hasford_class", "elts_class"),
                   .funs = function(x) ifelse(is.na(x), NA,
                                              ifelse(x < 2, 0, 1))) %>%
  dplyr::mutate(lipid_total_periphery_proportion = lipid_total_periphery_proportion*100)


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


# Define variable sets
top_features1 = df %>%
  dplyr::select(contains("Lipid"), contains("WBC"), contains("Mono"), contains("Mega"), contains("Proe")) %>% #, contains("inib"), contains("2nd")) %>%
  names()
clin = names(df)[!names(df) %in% top_features1]

pvalue_df = data.frame()

for (variable1 in unique(clin)) {
  
  print(variable1)
  
  df$tmp = df[[variable1]]
  
  # Wilcoxon
  ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
  multiple_t_tests_p_value <- lapply(df[top_features1], function(x) wilcox.test(x ~ df$tmp, na.rm=TRUE, exact=FALSE))
  ### P-values can be extracted from the result object
  pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
  ### Create a matrix and dataframe of the p-values
  pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
    #### Add the p values to a new dataframe
    p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"),
    #### Add also the t values, 95%CI to the same dataframe
    median1 = unlist(lapply(df[top_features1], function(x) median(x[df$tmp==1], na.rm=TRUE))),  # 1 = High
    median2 = unlist(lapply(df[top_features1], function(x) median(x[df$tmp==0], na.rm=TRUE)))   # 0 = Low
  ) %>%
    ### Rownames to column
    rownames_to_column() %>%
    mutate(clin_var = variable1) %>%
    rename(pvalue = p.value)
  # arrange(pvalue)
  rownames(pvalue_df1) = NULL
  
  pvalue_df <- rbind(pvalue_df, pvalue_df1)
}


# Rename features
# pvalue_df1 = pvalue_df %>%
#   dplyr::filter(!str_detect(rowname, "percentile|perimeter")) %>%
#   dplyr::mutate(clin_var = ifelse(clin_var == "gender", "Gender, female (yes/no)", clin_var),
#                 clin_var = ifelse(clin_var == "age_at_dg", "Age", clin_var),
#                 clin_var = ifelse(clin_var == "sokal_class", "Sokal, high (yes/no)", clin_var),
#                 clin_var = ifelse(clin_var == "hasford_class", "Hasford, high (yes/no)", clin_var),
#                 clin_var = ifelse(clin_var == "eutos_class", "EUTOS, high (yes/no)", clin_var),
#                 clin_var = ifelse(clin_var == "elts_class", "ELTS, high (yes/no)", clin_var),
#                 clin_var = ifelse(clin_var == "spleen_size", "Spleen size (cm)", clin_var),
#                 clin_var = ifelse(clin_var == "b_trom", "PB PLT (E9/L)", clin_var),
#                 clin_var = ifelse(clin_var == "l_baso", "PB Baso (%)", clin_var),
#                 clin_var = ifelse(clin_var == "l_blast", "PB Blast (%)", clin_var),
#                 clin_var = ifelse(clin_var == "l_eos", "PB Eos (%)", clin_var),
#                 clin_var = ifelse(clin_var == "l_lymf", "PB Lymph (%)", clin_var),
#                 clin_var = ifelse(clin_var == "bm_blast", "BM Blast (%)", clin_var),
#                 rowname = gsub("_", " ",
#                                gsub("_percentage", " (%)",
#                                     gsub("b_leuk", "PB WBC (E9/L)",
#                                          gsub("_mean", " (mean)",
#                                               gsub("_percentile_25", " (25%)",
#                                                    gsub("_percentile_75", " (75%)",
#                                                         gsub("_median", " (median)",
#                                                              gsub("_total", "",
#                                                                   gsub("_proportion", " (%/Area)",
#                                                                        gsub("sec_gen_tki", "2nd gen TKI (yes/no)",
#                                                                             gsub("imatinib", "Imatinib (yes/no)",
#                                                                                  gsub("dasatinib", "Dasatinib (yes/no)", rowname
#                                                                                  )))))))))))),
#                 rowname = str_to_sentence(rowname),
#                 rowname = gsub("Pb wbc \\(e9\\/l\\)", "PB WBC (E9/L)", 
#                                gsub("2nd gen tki \\(yes\\/no\\)", "2nd gen TKI (yes/no)",
#                                     gsub("cell ", "cell\n", 
#                                          gsub("nuclei ", "nuclei\n", rowname))))
#   )

pvalue_df1 = pvalue_df %>%
  dplyr::filter(!str_detect(rowname, "5%")) %>%
  dplyr::mutate(rowname = gsub("cell ", "cell\n",
                               gsub("nuclei ", "nuclei\n", rowname))
  )
unique(pvalue_df1$rowname)
unique(pvalue_df1$clin_var)

## Calculate FC and -log10 P value
pvalue_df1 <- pvalue_df1 %>%
  mutate(FC = ifelse(median1 == 0 & median2 == 0, 1, ifelse(median2 == 0, median1+1, ifelse(median1 == 0, 1/median2, median1 / median2)))) %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    neg_log10_p_adj = -log10(p_adjusted),
    p_adjusted_cat = ifelse(p_adjusted<0.001, 0.001, 
                            ifelse(p_adjusted<0.01, 0.01,
                                   ifelse(p_adjusted<0.05, 0.05,
                                          ifelse(p_adjusted<0.1, 0.1, "ns")))),
    pvalue_cat = ifelse(pvalue<0.001, 0.001, 
                        ifelse(pvalue<0.01, 0.01,
                               ifelse(pvalue<0.05, 0.05,
                                      ifelse(pvalue<0.1, 0.1, "ns")))),
    pvalue_cat_size = ifelse(pvalue<0.001, 0.001, 
                             ifelse(pvalue<0.01, 0.01,
                                    ifelse(pvalue<0.05, 0.05,
                                           ifelse(pvalue<0.1, 0.1, 0.5)))),
    pvalue_cat_size = -log10(pvalue_cat_size),
    FC_log = log(FC) )


# Balloonplot
p <- ggballoonplot(pvalue_df1, x = "rowname", y = "clin_var",
                   fill = "FC_log",
                   size = "neg_log10_p",
                   # size.range = c(1, 10),
                   ggtheme = theme_bw()) +
  scale_size(breaks = c(0, 1, 2, 3), range = c(1, 10), limits = c(1, max(pvalue_df1$neg_log10_p))) +
  # scale_fill_viridis_c(option = "C") +
  # gradient_fill(c("blue", "white", "red")) +
  # ylab("Clinical variables") +
  # xlab("MMR biomarkers") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       limits = c(min(pvalue_df1$FC_log), max(pvalue_df1$FC_log)), na.value = "gray") +
  guides(size = guide_legend(title="LOG10 P", nrow = 1, title.vjust = 0.5),
         fill = guide_colorbar(title="LOG10 FC", title.vjust = 0.75)) +
  font("xy.text", size = 10, color = "black", face="plain") +
  theme_bw() +
  theme(#axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
    axis.text.x = element_text(colour="black", angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour="black"),
    legend.position = "bottom", legend.margin=margin()); p
ggsave(plot = p, filename = paste0(results, "_response/Balloonplot_top_features.png"),
       width = 6, height = 5, dpi = 300, units = "in", bg = "white")

# Balloonplot adjusted p
p <- ggballoonplot(pvalue_df1, x = "rowname", y = "clin_var",
                   fill = "FC_log",
                   size = "neg_log10_p_adj",
                   # size.range = c(1, 10),
                   ggtheme = theme_bw()) +
  scale_size(breaks = c(0, 1, 2, 3), range = c(1, 10), limits = c(1, max(pvalue_df1$neg_log10_p_adj))) +
  # scale_fill_viridis_c(option = "C") +
  # gradient_fill(c("blue", "white", "red")) +
  # ylab("Clinical variables") +
  # xlab("MMR biomarkers") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       limits = c(min(pvalue_df1$FC_log), max(pvalue_df1$FC_log)), na.value = "gray") +
  guides(size = guide_legend(title="LOG10 AdjP", nrow = 1, title.vjust = 0.5),
         fill = guide_colorbar(title="LOG10 FC", title.vjust = 0.75)) +
  font("xy.text", size = 10, color = "black", face="plain") +
  theme_bw() +
  theme(#axis.title.y = element_text(size=12, colour="black", face="bold", angle = 90),
    axis.text.x = element_text(colour="black", angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour="black"),
    legend.position = "bottom", legend.margin=margin()); p
ggsave(plot = p, filename = paste0(results, "_response/Balloonplot_top_features_adjp.png"),
       width = 6, height = 5, dpi = 300, units = "in", bg = "white")


top_features1 = top_features1[grep(x = top_features1, pattern = "yes", invert = TRUE)]
# Scatter plots
for (i in top_features1) {
  for (j in clin[3:6]) {
    
    # Reset
    df1 = df
    df1$y1 = df1[[i]]
    df1$x1 = factor(df1[[j]])
    
    # Remove NA
    df1 = df1 %>%
      dplyr::filter(!(is.na(x1) | is.na(y1)))
    
    # Edit variables
    df1$x1 = ifelse(j == "Gender, female (yes/no)" & df1$x1 == 1, "Female",
                    ifelse(j == "Gender, female (yes/no)" & df1$x1 == 0, "Male",
                           ifelse(str_detect(j, "ELTS|Sokal|Has") & df1$x1 == 0, "Low-Int",
                                  ifelse(str_detect(j, "EUTOS") & df1$x1 == 0, "Low",
                                         ifelse(str_detect(j, "ELTS|EUTOS|Sokal|Has") & df1$x1 == 1, "High",
                                                ifelse(str_detect(j, "yes") & df1$x1 == 0, "FALSE",
                                                       ifelse(str_detect(j, "yes") & df1$x1 == 1, "TRUE",
                                                              ifelse(df1$x1 == 1, "High", "Low"))))))))
    df1$x1 = factor(df1$x1, levels = unique(df1$x1)[order(unique(df1$x1), decreasing = TRUE)])
    i = gsub(",[[:print:]]*", "", i)
    j = gsub(" \\([[:print:]]*", "", j)
    j = gsub(" female| high", "", j)
    j = gsub(",", "", j)
    
    
    g <- ggplot(data = df1, aes(y = y1, x = x1)) +
      geom_jitter(size=5, width = 0.2, aes(fill=x1), shape = 21, color = "black") +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      labs(x=j, y=i) +
      # ylim(0, 105) +
      theme_bw() +
      theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
            axis.text.y = element_text(size=13, colour = "black", face="bold"),
            axis.title=element_text(size=14, face="bold", colour = "black"),
            legend.position = "none") +
      scale_fill_brewer(palette = c("Set1"), direction = -1) +
      # scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
      stat_compare_means(method = "wilcox.test",
                         # label = "p.signif",
                         # label.y = 12,
                         label.x = 1.2,
                         size = 6); g
    ggsave(plot = g, filename = paste0(results, "_response/Scatterplots/", janitor::make_clean_names(j), "_", janitor::make_clean_names(i), ".png"),
           width = 5, height = 5, units = "in", dpi = 300)
    
  }
}

