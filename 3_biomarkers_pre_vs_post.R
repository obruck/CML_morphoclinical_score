rm(list = ls())

# Load necessary libraries
source("~/mounts/research/src/Rfunctions/library.R")
library(survminer)
library(MASS)
library(cowplot)

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")


###### Process data ########


# Diagnostic data
df0 = readxl::read_xlsx(paste0(export, "_response/full_data.xlsx"))

# All data
df = readRDS(paste0(export, "_response/data_for_modelling_all_cml.rds")) %>%
  dplyr::filter(henkilotunnus %in% df0$henkilotunnus) %>%
  dplyr::left_join(df0 %>% dplyr::select(henkilotunnus, mmr_dg = mmr, b_leuk_dg = b_leuk)) %>%
  dplyr::mutate(b_leuk_dg1 = ifelse(b_leuk_dg > median(b_leuk_dg, na.rm=TRUE), 1, 0)) %>%
  dplyr::filter(!is.na(mgg_time_diff)) %>%
  dplyr::filter(!paste0(database, "_", album_id) %in% paste0(df0$database, "_", df0$album_id)) %>%
  mutate(
    imatinib = ifelse(first_TKI=="imatinib", 1, 0),
    nilotinib = ifelse(first_TKI=="nilotinib", 1, 0),
    dasatinib = ifelse(first_TKI=="dasatinib", 1, 0),
    sec_gen_tki = ifelse(imatinib == 1, 0, 1),
    Sokal_class = ifelse(Sokal_class == "Low", 0,
                         ifelse(Sokal_class == "Intermediate", 1,
                                ifelse(Sokal_class == "High", 2, NA))),
    Hashford_class = ifelse(Hashford_class == "Low", 0,
                            ifelse(Hashford_class == "Intermediate", 1,
                                   ifelse(Hashford_class == "High", 2, NA))),
    ELTS_class = ifelse(ELTS_class == "Low", 0,
                        ifelse(ELTS_class == "Intermediate", 1,
                               ifelse(ELTS_class == "High", 2, NA))),
    EUTOS_class = ifelse(EUTOS_class == "Low", 0,
                         ifelse(EUTOS_class == "High", 1, NA)),
    MMR_time = as.numeric(MMR_time)/30.44,
    MR4.0_time = as.numeric(MR4.0_time)/30.44,
    MR4.5_time = as.numeric(MR4.5_time)/30.44,
    response_at_48 = ifelse(!is.na(MMR) & MMR_time<24 & MMR == 0, NA,
                            ifelse(MMR_time >= 24, FALSE, TRUE)),
    HUS_or_AUS_pt = ifelse(str_detect(henkilotunnus, '^(02139|AUS)'), 1, 0)) %>%
  mutate(across(c(dasatinib, imatinib, nilotinib), as.integer)) %>%
  dplyr::select(-c(first_TKI, MR4.0, MR4.5, MR4.0_time, MR4.5_time)) %>%
  janitor::clean_names() %>%
  # Define endpoints
  dplyr::mutate(
    mmr_time = ifelse(is.na(mmr_time), 0, mmr_time),
    mmr_time1=mmr_time,
    mmr = response_at_48
  ) %>%
  dplyr::filter(!is.na(response_at_48)) %>%
  dplyr::filter(!(response_at_48==FALSE & mmr_time1<24))


## Replace NaN and Inf with NA
df[df == "NaN"] <- NA
df[df == "Inf"] <- NA
df[df == "-Inf"] <- NA

# Increase probability of MMR
# Response
summary(lm(proerythroblasts_percentage ~ mgg_time_diff, data = df[df$imatinib==TRUE & df$mmr == TRUE & df$mgg_time_diff>90 & df$mgg_time_diff<365,]))
summary(lm(proerythroblasts_percentage ~ mgg_time_diff, data = df[df$imatinib==TRUE & df$mmr == FALSE & df$mgg_time_diff>90 & df$mgg_time_diff<365,]))
summary(lm(monocytes_nuclei_perimeter_median ~ mgg_time_diff, data = df[df$imatinib==TRUE & df$mmr == TRUE & df$mgg_time_diff>90 & df$mgg_time_diff<365,]))
summary(lm(monocytes_nuclei_perimeter_median ~ mgg_time_diff, data = df[df$imatinib==TRUE & df$mmr == FALSE & df$mgg_time_diff>90 & df$mgg_time_diff<365,]))

ggplot(df[df$imatinib==TRUE & df$mgg_time_diff>50  & abs(df$mgg_time_diff)<365,]) +
  geom_point(aes(x = mgg_time_diff, y = monocytes_nuclei_perimeter_median, color = mmr))

median(df[df$imatinib==TRUE & df$mgg_time_diff>90  & abs(df$mgg_time_diff)<365 & df$mmr == FALSE,]$proerythroblasts_percentage, na.rm = TRUE)
median(df[df$imatinib==TRUE & df$mgg_time_diff>90  & abs(df$mgg_time_diff)<365 & df$mmr == TRUE,]$proerythroblasts_percentage, na.rm = TRUE)
wilcox.test(df[df$imatinib==TRUE & df$mgg_time_diff>150  & abs(df$mgg_time_diff)<365 & df$mmr == FALSE,]$proerythroblasts_percentage,
            df[df$imatinib==TRUE & df$mgg_time_diff>150  & abs(df$mgg_time_diff)<365 & df$mmr == TRUE,]$proerythroblasts_percentage)
wilcox.test(df[df$mgg_time_diff>90  & abs(df$mgg_time_diff)<270 & df$mmr == FALSE,]$monocytes_nuclei_perimeter_median,
            df[df$mgg_time_diff>90  & abs(df$mgg_time_diff)<270 & df$mmr == TRUE,]$monocytes_nuclei_perimeter_median)

tmp = ifelse(df[df$mgg_time_diff>30  & abs(df$mgg_time_diff)<190,]$monocytes_nuclei_perimeter_median>1, 1, 0)
chisq.test(tmp, df[df$mgg_time_diff>30  & abs(df$mgg_time_diff)<190,]$mmr)
table(tmp, df[df$mgg_time_diff>30  & abs(df$mgg_time_diff)<190,]$mmr)
nrow(df[df$mgg_time_diff>30  & abs(df$mgg_time_diff)<190 & df$mmr == FALSE,]) #df$imatinib==TRUE & 
nrow(df[df$mgg_time_diff>30  & abs(df$mgg_time_diff)<190 & df$mmr == TRUE,]) #df$imatinib==TRUE & 


df1 = df %>%
  dplyr::select(henkilotunnus, mmr, mgg_time_diff, monocytes_nuclei_perimeter_median, proerythroblasts_percentage) %>%
  dplyr::filter(mgg_time_diff > 50 & mgg_time_diff < 200) %>%
  group_by(henkilotunnus, mmr) %>%
  summarise(monocytes_nuclei_perimeter_median = median(monocytes_nuclei_perimeter_median, na.rm=TRUE),
            proerythroblasts_percentage = median(proerythroblasts_percentage, na.rm=TRUE))
wilcox.test(df1$monocytes_nuclei_perimeter_median ~ df1$mmr, na.rm = TRUE)
wilcox.test(df1$proerythroblasts_percentage ~ df1$mmr, na.rm = TRUE)

ggplot(data = df) +
  geom_point(aes(x = mgg_time_diff, y = proerythroblasts_percentage, color=mmr_dg))

ggplot(data = df) +
  geom_point(aes(x = mmr_dg, y = proerythroblasts_percentage, color=mmr_dg))

g <- ggplot(data = df %>%
              dplyr::mutate(b_leuk1 = ifelse(b_leuk > median(b_leuk, na.rm=TRUE), 1, 0),
                            b_leuk1 = factor(b_leuk1),
                            b_leuk_dg1 = factor(b_leuk_dg1),
                            mmr_time2 = mmr_time1-mgg_time_diff) %>%
              dplyr::filter(mgg_time_diff > 180 & mgg_time_diff < 365 & !is.na(b_leuk_dg1)),
            aes(y = proerythroblasts_percentage, x = b_leuk_dg1)) +
  geom_jitter(size=5, width = 0.2, aes(fill=mmr), shape = 21, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point() +
  # labs(x="BM Erythroid precursors (%)", y="Time to MR4.0 (years)") +
  # ylim(0, 10) +
  theme_bw() +
  theme(axis.text.x = element_text(size=13, colour = "black", face="bold"),
        axis.text.y = element_text(size=13, colour = "black", face="bold"),
        axis.title=element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  scale_fill_brewer(palette = c("Set1"), direction = -1) +
  # scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
  stat_compare_means(method = "wilcox.test",
                     # label = "p.signif",
                     # label.y = 13,
                     label.x = 1.2,
                     size = 6); g
