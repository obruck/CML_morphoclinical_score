# Validation with the AUS set #2


# Library
source("~/library.R")
library(ggpubr)

# Read data
df = read_xlsx("./Helsinki Full TFR cohort data _de-identified.xlsx", sheet = "CML Full TFR cohort data")
df1 = read_xlsx("./Helsinki Full TFR cohort data _de-identified.xlsx", sheet = "TFR cohort BM data")
df = df %>%
  dplyr::left_join(df1) %>%
  janitor::clean_names() %>%
  dplyr::filter(!id_number == "12") %>%
  dplyr::mutate(mmr_time = (as.Date(date_of_0_1_percent) - as.Date(date_of_tki_commencement))/365.24,
                ccyr_time = (as.Date(date_of_1_percent) - as.Date(date_of_tki_commencement))/365.24)

# Process data
df = df %>%
  janitor::clean_names() %>%
  dplyr::mutate(time_to_mr4 = ifelse(is.na(time_to_mr4_yrs), NA,
                                     ifelse(time_to_mr4_yrs < 3, "<3 years", "≥3 years")),
                wbc = ifelse(pb_white_cell_count_e9_l > 50, ">50 (E9/L)", "≤50 (E9/L)"),
                wbc = factor(wbc, levels = c("≤50 (E9/L)", ">50 (E9/L)")),
                erytroid = as.numeric(erythroid_precursors_erythroblasts_percent),
                erytroid1 = ifelse(erytroid >= median(erytroid, na.rm=TRUE), "High", "Low"),
                erytroid1 = factor(erytroid1, levels = c("Low", "High")))

# WBC and CCYR
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(wbc)), aes(y = ccyr_time, x = wbc)) +
  geom_jitter(size=5, width = 0.2, aes(fill=wbc), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="PB WBC (E9/L)", y="Time to CCyR (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 1.5,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./CCYR_time_WBC.png", width = 5, height = 5, units = "in", dpi = 300)


# Erytroid and CCYR
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(erytroid1)), aes(y = ccyr_time, x = erytroid1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=erytroid1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="BM Erythroid precursors (%)", y="Time to CCyR (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 1.5,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./CCYR_time_erytroid.png", width = 5, height = 5, units = "in", dpi = 300)


# WBC and MMR
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(wbc)), aes(y = mmr_time, x = wbc)) +
  geom_jitter(size=5, width = 0.2, aes(fill=wbc), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="PB WBC (E9/L)", y="Time to MMR (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 7.5,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MMR_time_WBC.png", width = 5, height = 5, units = "in", dpi = 300)


# Erytroid and MMR
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(erytroid1)), aes(y = mmr_time, x = erytroid1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=erytroid1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="BM Erythroid precursors (%)", y="Time to MMR (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 5.5,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MMR_time_erytroid.png", width = 5, height = 5, units = "in", dpi = 300)



# WBC and MR4.0
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(wbc)), aes(y = time_to_mr4_yrs, x = wbc)) +
  geom_jitter(size=5, width = 0.2, aes(fill=wbc), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="PB WBC (E9/L)", y="Time to MR4.0 (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 15,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MR4_time_WBC.png", width = 5, height = 5, units = "in", dpi = 300)


# Erytroid and MR4.0
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(erytroid1)), aes(y = time_to_mr4_yrs, x = erytroid1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=erytroid1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="BM Erythroid precursors (%)", y="Time to MR4.0 (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 13,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MR4_time_erytroid.png", width = 5, height = 5, units = "in", dpi = 300)
