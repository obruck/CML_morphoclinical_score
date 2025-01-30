# Barplots

# Library
source("~/library.R")
library(ggpubr)

# Read data
df1 = read_xlsx("./coefficients_confint.xlsx")
df = read_xlsx("./coefficients.xlsx") %>%
  bind_cols(df1) %>%
  janitor::clean_names() %>%
  dplyr::mutate(covariates = ifelse(covariates == "elts_class", "ELTS",
                                    ifelse(covariates == "monocytes_nuclei_perimeter_median", "Monocytes nuclei\nperimeter (median)",
                                           ifelse(covariates == "imatinib", "Imatinib (yes/no)",
                                                  ifelse(covariates == "b_leuk", "PB WBC (E9/L)",
                                                         ifelse(covariates == "proerythroblasts_percentage", "Proerythroblasts (%)", covariates))))),
                HR = ifelse(exp_coef_2<1, -exp_coef_2, exp_coef_2),
                lower_95_1 = ifelse(lower_95<1, -lower_95, lower_95),
                upper_95_1 = ifelse(upper_95<1, -upper_95, upper_95),
                HRqual = ifelse(exp_coef<1, "L", "H"),
                p_value = ifelse(round(pr_z, 2)!=0, round(pr_z, 2),
                                 ifelse(round(pr_z, 3)!=0, round(pr_z, 3), round(pr_z, 4))))

# Barplot
g = ggplot(data = df, aes(x = covariates, y = HR, fill = HRqual)) +
  geom_errorbar(aes(ymin=lower_95_1, ymax=upper_95_1), size = 0.8, width = 0.3) +
  geom_col(color="black") +
  geom_text(aes(x = covariates, hjust = ifelse(exp_coef_2<1, 1.1, -0.1), label = paste0("p=", p_value)),
            y = 0, color="white") +
  theme_bw() +
  xlab("") +
  ylab("Hazard ratio of MMR in 2 years") +
  scale_fill_brewer(palette = c("Set1")) +
  theme(axis.text = element_text(size=11, colour = "black"),
        axis.title=element_text(size=12, colour = "black"),
        legend.position = "None") +
  coord_flip(); g
ggsave(plot = g, filename = "./Barplot_coefficients.png",
       width = 6, height = 3, units = "in", dpi = 300)

