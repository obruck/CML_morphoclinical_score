# Validation with the FIN set


# Library
source("~/library.R")
library(ggpubr)

# Read data
df1 = read_xlsx("./univariate_cox_results_top_features_clinical_balloonplot.xlsx")
df = read_xlsx("./full_data.xlsx") %>%
  dplyr::select(henkilotunnus, mmr, mmr_time, mr4_0, mr4_0_time, names(df1))
response = fread("./CML_response.csv") %>%
  dplyr::filter(henkilotunnus %in% df$henkilotunnus) %>%
  dplyr::filter(CCyR==TRUE) %>%
  dplyr::select(henkilotunnus, CCyR_time)
df = df %>%
  dplyr::left_join(response) %>%
  dplyr::mutate(mmr_time = mmr_time / 12,
                time_to_mr4_yrs = mr4_0_time / 12,
                ccyr_time = CCyR_time / 365.24)

# Process data
df = df %>%
  janitor::clean_names() %>%
  dplyr::mutate(time_to_mr4 = ifelse(is.na(time_to_mr4_yrs), NA,
                                     ifelse(time_to_mr4_yrs < 3, "<3 years", "≥3 years")),
                wbc = ifelse(b_leuk > 50, ">50 (E9/L)", "≤50 (E9/L)"),
                wbc = factor(wbc, levels = c("≤50 (E9/L)", ">50 (E9/L)")),
                erytroid = as.numeric(proerythroblasts_percentage),
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
                     label.y = 3.5,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./CCYR_time_WBC_FIN.png", width = 5, height = 5, units = "in", dpi = 300)


# Erytroid and CCYR
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(erytroid1)), aes(y = ccyr_time, x = erytroid1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=erytroid1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="BM Proerythroblasts (%)", y="Time to CCyR (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 3.5,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./CCYR_time_erytroid_FIN.png", width = 5, height = 5, units = "in", dpi = 300)


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
                     label.y = 10,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MMR_time_WBC_FIN.png", width = 5, height = 5, units = "in", dpi = 300)


# Erytroid and MMR
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(erytroid1)), aes(y = mmr_time, x = erytroid1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=erytroid1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="BM Proerythroblasts (%)", y="Time to MMR (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 10,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MMR_time_erytroid_FIN.png", width = 5, height = 5, units = "in", dpi = 300)



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
                     label.y = 14,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MR4_time_WBC_FIN.png", width = 5, height = 5, units = "in", dpi = 300)


# Erytroid and MR4.0
g <- ggplot(data = df %>%
              dplyr::filter(!is.na(erytroid1)), aes(y = time_to_mr4_yrs, x = erytroid1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=erytroid1), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(x="BM Proerythroblasts (%)", y="Time to MR4.0 (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  stat_compare_means(method = "wilcox.test",
                     label.y = 14,
                     label.x = 1.2,
                     size = 6); g
ggsave(plot = g, filename = "./MR4_time_erytroid_FIN.png", width = 5, height = 5, units = "in", dpi = 300)
