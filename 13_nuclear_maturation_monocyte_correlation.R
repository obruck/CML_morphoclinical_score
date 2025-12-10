# Script: 13_nuclear_maturation_monocyte_correlation.R

# Load necessary libraries
rm(list = ls())

# Load necessary libraries
source("~/mounts/research/src/Rfunctions/library.R")

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")


###### Process data ########
mono_cols <- c(
  "Monocytes_Percentage",
  "Monocytes_cell-perimeter_median",
  "Monocytes_nuclei-area_mean",
  "Monocytes_nuclei-area_median",
  "Monocytes_nuclei-perimeter_median"
)

# Data
df = readRDS(paste0(export, "_response/data_for_modelling1.rds")) %>%
  dplyr::filter(!(is.na(MMR_time) | MMR_time == 0)) %>%
  mutate(
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

df <- df %>% select(henkilotunnus, HUS_or_AUS_pt, l_mono, mono_cols)

# Process data
lab = readRDS("mounts/research/husdatalake/disease/processed_data/CML/lab_dg.rds") %>%
  dplyr::filter(tutkimus_lyhenne %in% c("B -Monos", "L -Mono", "B -Leuk"))

lab <- lab[, c("henkilotunnus", "tutkimus_lyhenne", "tulos_max")]
lab <- lab %>%
  pivot_wider(
    names_from = tutkimus_lyhenne,
    values_from = tulos_max
  )


df <- df %>%
  left_join(lab, by = "henkilotunnus") %>%
  mutate(l_mono = ifelse(is.na(l_mono), `L -Mono`, l_mono)) %>%
  select(-`L -Mono`)

# Print the column names
cat('Columns to look at:\n')
print(mono_cols)

# Check if PB monocyte values are NA 
sum(!is.na(df$l_mono))

colnames(df)[colnames(df) == "Monocytes_nuclei-perimeter_median"] <- "Monocytes_nuclei_perimeter_median"

ggscatter(df, x = "Monocytes_nuclei_perimeter_median", y = "l_mono",
          add = "reg.line",        # add linear regression line
          conf.int = TRUE,         # add confidence interval
          cor.coef = TRUE,         # add correlation coefficient
          cor.method = "spearman", # or "pearson"
          xlab = "Monocytes nuclei perimeter (median)", ylab = "l_mono") +
  theme_minimal()

# Optional: scatterplot with loess smoother (for non-linear trends)
ggplot(df, aes(x = Monocytes_nuclei_perimeter_median, y = l_mono)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(x = "Monocytes nuclei-perimeter (median)",
       y = "l_mono",
       title = "Scatterplot with LOESS smoother")

# Histograms
hist(df$Monocytes_nuclei_perimeter_median, main = "Histogram of Monocytes_nuclei_perimeter_median", xlab = "Monocytes nuclei_perimeter (median)")
hist(df$l_mono, main = "Histogram of l_mono", xlab = "l_mono")

# Q-Q plots
qqnorm(df$Monocytes_nuclei_perimeter_median); qqline(df$Monocytes_nuclei_perimeter_median)
qqnorm(df$l_mono); qqline(df$l_mono)

# Shapiro-Wilk tests
shapiro.test(df$Monocytes_nuclei_perimeter_median)
shapiro.test(df$l_mono)

# Spearman correlation (safe for non-normal data)
spearman_cor <- cor.test(df$Monocytes_nuclei_perimeter_median, df$l_mono, method = "spearman", exact = FALSE)
spearman_cor

# Pearson correlation (optional, if relationship looks roughly linear)
pearson_cor <- cor.test(df$Monocytes_nuclei_perimeter_median, df$l_mono, method = "pearson")
pearson_cor

# Function to compute Spearman correlation and p-value
cor_fun <- function(x, y) {
  res <- cor.test(x, y, method = 'spearman', use = 'complete.obs')
  return(c(cor = res$estimate, p = res$p.value))
}

# Correlation with l_mono
df = as.data.frame(df)
df[3:8] = sapply(df[3:8], as.numeric)
cor_df = data.frame()
for (i in 4:8){
  res = cor.test(df[[i]], df$l_mono, method = 'spearman', use = 'complete.obs')
  cor_df = rbind(cor_df, data.frame(feature=mono_cols[i-3],
                                    cor.rho=res$estimate,
                                    p=res$p.value))
}

print(cor_df)

# Make sure columns are numeric/factor
cor_df <- cor_df %>%
  mutate(
    cor.rho = as.numeric(cor.rho),
    p = as.numeric(p),
    feature = as.factor(feature),
    significant = ifelse(p < 0.05, "*", "")
  )

cor_results <- bind_rows(cor_df)

# Balloon plot of all correlations
cor_results$cor.rho <- as.numeric(cor_results$cor.rho)
cor_results$p <- as.numeric(cor_results$p)

# Prepare for balloon plot
balloon_data <- cor_df
balloon_data$plot_color <- ifelse(
  balloon_data$p > 0.05, "p > 0.05",
  ifelse(balloon_data$cor > 0, "Positive", "Negative")
)
balloon_data$target = "PB Mono (%)"
balloon_data <- balloon_data %>%
  mutate(feature = factor(
    feature,
    levels = c(
      "Monocytes_Percentage",
      "Monocytes_nuclei-perimeter_median",
      "Monocytes_nuclei-area_median",
      "Monocytes_nuclei-area_mean",
      "Monocytes_cell-perimeter_median"
    ),
    labels = c(
      "Monocytes (%)",
      "Monocyte nuclear perimeter (median)",
      "Monocyte nuclear area (median)",
      "Monocyte nuclear area (mean)",
      "Monocyte cell perimeter (median)"
    )
  ))

# Balloonplot
balloon_data$corrho = abs(balloon_data$cor.rho)
balloon_data$plot_color1 = factor(balloon_data$plot_color, levels = c("Positive", "Negative", "p > 0.05"))
g = ggballoonplot(balloon_data, x = "target", y = "feature",
              fill = "plot_color1",
              size = "corrho",
              ggtheme = theme_bw()) +
  scale_size(breaks = c(0.1, 0.2, 0.3)) +
  scale_fill_manual(values = c("red", "grey80", "blue")) +
  guides(size = guide_legend(title="LOG10\nAdjP", nrow = 3, title.vjust = 0.5),
         fill = guide_legend(title="LOG10 FC", title.vjust = 0.75, override.aes = list(size=3))) +
  font("xy.text", size = 10, color = "black", face="plain") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour="black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour="black"),
    legend.position = "right", legend.direction = "vertical", legend.margin=margin()); g

ggsave(plot = g, filename = paste0(results, "_response/Balloonplot_monocyte.png"),
       width = 4, height = 5, dpi = 300, units = "in", bg = "white")

