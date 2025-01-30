# Balloonplot

# Library
source("~/mounts/research/src/Rfunctions/library.R")
library(ggpubr)
library(rstatix)

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")

dir.create(paste0(results, "_response/Scatterplots_WBCdiff"))


# Data
df0 = readRDS(paste0(export, "_response/data_for_modelling1.rds")) %>%
  dplyr::filter(!(is.na(MMR_time) | MMR_time == 0)) %>%
  mutate(
    l_mono = ifelse(l_mono==70, NA, l_mono),
    imatinib = ifelse(first_TKI=="imatinib", 1, 0),
    nilotinib = ifelse(first_TKI=="nilotinib", 1, 0),
    dasatinib = ifelse(first_TKI=="dasatinib", 1, 0),
    sec_gen_tki = ifelse(imatinib == 1, 0, 1),
    Sokal_class = ifelse(Sokal_class == "Low", 0,
                         ifelse(Sokal_class == "Intermediate", 1,
                                ifelse(Sokal_class == "High", 2, NA))),
    Hasford_class = ifelse(Hasford_class == "Low", 0,
                            ifelse(Hasford_class == "Intermediate", 1,
                                   ifelse(Hasford_class == "High", 2, NA))),
    ELTS_class = ifelse(ELTS_class == "Low", 0,
                        ifelse(ELTS_class == "Intermediate", 1,
                               ifelse(ELTS_class == "High", 2, NA))),
    EUTOS_class = ifelse(EUTOS_class == "Low", 0,
                         ifelse(EUTOS_class == "High", 1, NA)),
    MMR_time = as.numeric(MMR_time)/30.44,
    MR4.0_time = as.numeric(MR4.0_time)/30.44,
    MR4.5_time = as.numeric(MR4.5_time)/30.44,
    response_at_24 = ifelse(!is.na(MMR) & MMR_time<24 & MMR == 0, NA,
                            ifelse(MMR_time >= 24, FALSE, TRUE)),
    HUS_or_AUS_pt = ifelse(str_detect(henkilotunnus, '^(02139|AUS)'), 1, 0)) %>%
  mutate(across(c(dasatinib, imatinib, nilotinib), as.integer))

table(df0$response_at_24)

# Rename
df0 = df0 %>%
  janitor::clean_names()

# Define endpoints
df0 = df0 %>%
  dplyr::mutate(
    mmr_time = ifelse(is.na(mmr_time), 0, mmr_time),
    mmr_time1=mmr_time,
    mmr = response_at_24
    # mmr_time=ifelse(!is.na(response_at_24), 1/mmr_time, NA)
  ) %>%
  dplyr::filter(!is.na(response_at_24)) %>%
  dplyr::filter(!(response_at_24==FALSE & mmr_time1<24))


# Divide some variables with 10 to make the HR more interpretable
df0 = df0 %>%
  dplyr::mutate(b_leuk = b_leuk / 10) %>%
  dplyr::mutate_at(vars(matches('perimeter')), function(x) x = x/10)


## Replace NaN and Inf with NA
df0[df0 == "NaN"] <- NA
df0[df0 == "Inf"] <- NA
df0[df0 == "-Inf"] <- NA

df = df0

cor.test(df$b_leuk, df$l_neut, method="spearman")



# Scatter plots
## Categorize diff
for (i in unique(names(df)[grep(pattern = "^l_", x=names(df))])) {
  
  # Reset
  df1 = df
  df1$x1 = ifelse(df1[[i]] > median(df1[[i]], na.rm = TRUE), "High", "Low")
  df1$x1 = factor(df1$x1, levels = c("Low", "High"))
  
  
  # Remove NA
  df1 = df1 %>%
    dplyr::filter(!(is.na(b_leuk) | is.na(x1)))
  
  # Plot
  g <- ggplot(data = df1, aes(y = b_leuk, x = x1)) +
    geom_jitter(size=5, width = 0.2, aes(fill=x1), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    labs(x=gsub("Erblast", "Erythroblasts",
                gsub("Baso", "Basophils",
                     gsub("Eos", "Eosinophils",
                          gsub("Lymf", "Lymphocytes",
                               gsub("Mono", "Monocytes",
                                    gsub("Neut", "Neutrophils", paste0(paste0("PB ", (stringr::str_to_sentence(gsub("l_", "", i)))), " (%)"))))))),
         y="PB WBC (%)") +
    # ylim(0, 105) +
    theme_bw() +
    theme(axis.text.x = element_text(size=13, colour = "black"),
          axis.text.y = element_text(size=13, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1"), direction = -1) +
    # scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
    stat_compare_means(method = "wilcox.test",
                       # label = "p.signif",
                       # label.y = 12,
                       label.x = 1.2,
                       size = 6); g
  ggsave(plot = g, filename = paste0(results, "_response/Scatterplots_WBCdiff/WBC_", janitor::make_clean_names(i), ".png"),
         width = 5, height = 5, units = "in", dpi = 300)
  
}

## Categorize b_leuk
for (i in unique(names(df)[grep(pattern = "^l_", x=names(df))])) {
  
  # Reset
  df1 = df
  df1$y1 = df1[[i]]
  df1$b_leuk = ifelse(df1$b_leuk > median(df1$b_leuk, na.rm = TRUE), "High", "Low")
  df1$b_leuk = factor(df1$b_leuk, levels = c("Low", "High"))
  
  
  # Remove NA
  df1 = df1 %>%
    dplyr::filter(!(is.na(b_leuk) | is.na(y1)))
  
  # Plot
  g <- ggplot(data = df1, aes(x = b_leuk, y = y1)) +
    geom_jitter(size=5, width = 0.2, aes(fill=b_leuk), shape = 21, color = "black") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    labs(y=gsub("Erblast", "Erythroblasts",
                gsub("Baso", "Basophils",
                     gsub("Eos", "Eosinophils",
                          gsub("Lymf", "Lymphocytes",
                               gsub("Mono", "Monocytes",
                                    gsub("Neut", "Neutrophils", paste0(paste0("PB ", (stringr::str_to_sentence(gsub("l_", "", i)))), " (%)"))))))),
         x="PB WBC (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(size=13, colour = "black"),
          axis.text.y = element_text(size=13, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          legend.position = "none") +
    scale_fill_brewer(palette = c("Set1"), direction = -1) +
    # scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
    stat_compare_means(method = "wilcox.test",
                       # label = "p.signif",
                       # label.y = 12,
                       label.x = 1.2,
                       size = 6); g
  ggsave(plot = g, filename = paste0(results, "_response/Scatterplots_WBCdiff/WBC_", janitor::make_clean_names(i), "_1.png"),
         width = 5, height = 5, units = "in", dpi = 300)
  
}


# Scatter plot for correlation
## Categorize b_leuk
for (i in unique(names(df)[grep(pattern = "^l_", x=names(df))])) {
  
  # Reset
  df1 = df
  df1$x1 = df1[[i]]
  
  if (str_detect(i, "eos")) {
    df1$y1 = df1$eosinophils_percentage
    x2 = "PB Eosinophils (%)"
    y2 = "BM Eosinophils (%)"
    tick_breaks = seq(0, 25, by = 5)
  } else if (str_detect(i, "bas")) {
    df1$y1 = df1$basophils_percentage
    x2 = "PB Basophils (%)"
    y2 = "BM Basophils (%)"
    tick_breaks = seq(0, 30, by = 10)
  } else if (str_detect(i, "mono")) {
    df1$y1 = df1$monocytes_percentage
    x2 = "PB Monocytes (%)"
    y2 = "BM Monocytes (%)"
    tick_breaks = seq(0, 20, by = 5)
  } else if (str_detect(i, "neut")) {
    df1$y1 = df1$neutrophils_percentage
    x2 = "PB Neutrophils (%)"
    y2 = "BM Neutrophils (%)"
    tick_breaks = seq(0, 100, by = 25)
  } else if (str_detect(i, "lym")) {
    df1$y1 = df1$lymphocytes_percentage
    x2 = "PB Lymphocytes (%)"
    y2 = "BM Lymphocytes (%)"
    tick_breaks = seq(0, 55, by = 10)
  } else if (str_detect(i, "l_blas")) {
    df1$y1 = df1$blasts_percentage
    x2 = "PB Blasts (%)"
    y2 = "BM Blasts (%)"
    tick_breaks = seq(0, 20, by = 5)
  } else if (str_detect(i, "erblas")) {
    df1$y1 = df1$erythroblasts_percentage
    x2 = "PB Erythroblasts (%)"
    y2 = "BM Erythroblasts (%)"
    tick_breaks = seq(0, 25, by = 5)
  }
  
  
  # Remove NA
  df1 = df1 %>%
    dplyr::filter(!(is.na(y1) | is.na(x1)))
  
  # Plot
  g <- ggplot(data = df1, aes(x = x1, y = y1)) +
    geom_point(size=3, color = "black") +
    # geom_smooth(size=2, method='lm', formula= y~x) +
    labs(x = x2,
         y = y2) +
    # ylim(0, 105) +
    theme_bw() +
    coord_equal(ratio = 1) +
    expand_limits(x = 0, y = 0) +
    xlim(min(tick_breaks), max(tick_breaks)) +
    ylim(min(tick_breaks), max(tick_breaks)) +
    # scale_x_continuous(breaks = tick_breaks) +
    # scale_y_continuous(breaks = tick_breaks) +
    theme(axis.text.x = element_text(size=13, colour = "black"),
          axis.text.y = element_text(size=13, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          legend.position = "none") +
    # scale_fill_brewer(palette = c("Set1"), direction = -1) +
    stat_cor(method = "spearman",
             label.x = 1.2,
             size = 6); g
  ggsave(plot = g, filename = paste0(results, "_response/Scatterplots_WBCdiff/", janitor::make_clean_names(x2), "_", janitor::make_clean_names(y2), ".png"),
         width = 5, height = 5, units = "in", dpi = 300)
  
}

# Balloonplot
tmp = df %>%
  dplyr::mutate(
    b_leuk = ifelse(b_leuk > median(b_leuk, na.rm = TRUE), "High", "Low"),
    b_leuk = factor(b_leuk, levels = c("Low", "High"))
  ) %>%
  dplyr::select(b_leuk, starts_with("l_")) %>%
  melt(id.vars = "b_leuk") %>%
  ungroup()
tmp1 = tmp %>%
  dplyr::group_by(variable) %>%
  rstatix::wilcox_test(value~b_leuk, p.adjust.method = "BH") %>%
  ungroup()
tmp2 = tmp %>%
  ungroup() %>%
  dplyr::filter(b_leuk=="Low") %>%
  dplyr::group_by(variable) %>%
  summarise(median2 = mean(value, na.rm=TRUE)) %>%
  ungroup()
tmp3 = tmp %>%
  dplyr::filter(b_leuk=="High") %>%
  dplyr::group_by(variable) %>%
  summarise(median1 = mean(value, na.rm=TRUE)) %>%
  ungroup()
tmp1 = tmp1 %>%
  dplyr::left_join(tmp2) %>%
  dplyr::left_join(tmp3)
## Calculate FC and -log10 P value
pvalue_df1 <- tmp1 %>%
  mutate(FC = ifelse(median1 == 0 & median2 == 0, 1, ifelse(median2 == 0, median1+1, ifelse(median1 == 0, 1/median2, median1 / median2)))) %>%
  mutate(
    neg_log10_p = -log10(p),
    neg_log10_p_adj = -log10(p.adj),
    p_adjusted_cat = ifelse(p.adj<0.001, 0.001, 
                            ifelse(p.adj<0.01, 0.01,
                                   ifelse(p.adj<0.05, 0.05,
                                          ifelse(p.adj<0.1, 0.1, "ns")))),
    pvalue_cat = ifelse(p<0.001, 0.001, 
                        ifelse(p<0.01, 0.01,
                               ifelse(p<0.05, 0.05,
                                      ifelse(p<0.1, 0.1, "ns")))),
    pvalue_cat_size = ifelse(p<0.001, 0.001, 
                             ifelse(p<0.01, 0.01,
                                    ifelse(p<0.05, 0.05,
                                           ifelse(p<0.1, 0.1, 0.5)))),
    pvalue_cat_size = -log10(pvalue_cat_size),
    FC_log = log(FC),
    variable = paste0(paste0("PB ", (stringr::str_to_sentence(gsub("l_", "", variable)))), " (%)"),
    variable = gsub("Erblast", "Erythroblasts",
                    gsub("Baso", "Basophils",
                                    gsub("Eos", "Eosinophils",
                                         gsub("Lymf", "Lymphocytes",
                                              gsub("Mono", "Monocytes",
                                                   gsub("Neut", "Neutrophils", variable)))))),
    group1 = "PB WBC (E9/L)\nHigh vs. low")


# Balloonplot
p <- ggballoonplot(pvalue_df1, x = "variable", y = "group1",
                   fill = "FC",
                   size = "neg_log10_p_adj",
                   # size.range = c(1, 10),
                   ggtheme = theme_bw()) +
  scale_size(breaks = c(0, 1, 2, 3), range = c(1, 10), limits = c(1, max(pvalue_df1$neg_log10_p_adj))) +
  # scale_fill_viridis_c(option = "C") +
  # gradient_fill(c("blue", "white", "red")) +
  # ylab("Clinical variables") +
  # xlab("MMR biomarkers") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 1,
                       # breaks = c(1, 2, 3),
                       limits = c(min(pvalue_df1$FC), max(pvalue_df1$FC)), na.value = "gray") +
  guides(size = guide_legend(title="LOG10 AdjP", nrow = 1, title.vjust = 0.5),
         fill = guide_colorbar(title="LOG10 FC", title.vjust = 0.75)) +
  font("xy.text", size = 10, color = "black", face="plain") +
  theme_bw() +
  theme(#axis.title.y = element_text(size=12, colour="black", angle = 90),
    axis.text.x = element_text(colour="black", angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour="black", hjust = 0.5),
    legend.position = "bottom", legend.margin=margin()); p
ggsave(plot = p, filename = paste0(results, "_response/Balloonplot_WBC_diff.png"),
       width = 5.5, height = 2.5, dpi = 300, units = "in", bg = "white")

