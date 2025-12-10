# Barplots

# Library
source("~/library.R")
library(ggpubr)
library(cowplot)

c_index_model <- 0.67366533

# Read C-index data
c_index_df = read_xlsx(paste0(results, "_response/c_indexes.xlsx")) %>%
  janitor::clean_names() %>%
  dplyr::mutate(covariates = ifelse(covariates == "elts_class", "ELTS",
                                    ifelse(covariates == "monocytes_nuclei_perimeter_median", "Monocytes nuclei\nperimeter (median)",
                                           ifelse(covariates == "imatinib", "Imatinib (yes/no)",
                                                  ifelse(covariates == "b_leuk", "PB WBC (E9/L)",
                                                         ifelse(covariates == "proerythroblasts_percentage", "Proerythroblasts (%)", covariates))))))

# Read data
df1 = read_xlsx(paste0(results, "_response/coefficients_confint.xlsx"))
df = read_xlsx(paste0(results, "_response/coefficients.xlsx")) %>%
  bind_cols(df1) %>%
  janitor::clean_names() %>%
  dplyr::mutate(covariates = ifelse(covariates == "elts_class", "ELTS",
                                    ifelse(covariates == "monocytes_nuclei_perimeter_median", "Monocytes nuclei\nperimeter (median)",
                                           ifelse(covariates == "imatinib", "Imatinib (yes/no)",
                                                  ifelse(covariates == "b_leuk", "PB WBC (E9/L)",
                                                         ifelse(covariates == "proerythroblasts_percentage", "Proerythroblasts (%)", covariates))))),
                # HR = ifelse(exp_coef_2<1, -exp_coef_2, exp_coef_2),
                HR = exp_coef_2,
                lower_95_1 = ifelse(lower_95<1, -lower_95, lower_95),
                upper_95_1 = ifelse(upper_95<1, -upper_95, upper_95),
                HRqual = ifelse(exp_coef<1, "L", "H"),
                p_value = ifelse(round(pr_z, 2)!=0, round(pr_z, 2),
                                 ifelse(round(pr_z, 3)!=0, round(pr_z, 3), round(pr_z, 4))))

desired_order <- c("Proerythroblasts (%)",
                   "Imatinib (yes/no)",
                   "PB WBC (E9/L)",
                   "Monocytes nuclei\nperimeter (median)",
                   "ELTS")

df <- df %>% dplyr::mutate(covariates = factor(covariates, levels = desired_order))
c_index_df <- c_index_df %>% dplyr::mutate(
  covariates = factor(covariates, levels = levels(df$covariates))
)

c_index_df$change <- c_index_df$c_indexes - c_index_model

# Barplot
g1 = ggplot(data = df, aes(x = covariates, y = HR)) +
  geom_linerange(aes(ymin=lower_95, ymax=upper_95, color = HRqual), size = 1.5) +
  geom_point(size = 5, shape = 18, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_label(aes(x = covariates, vjust = -0.36, #hjust = 0.5,
                 label = paste0("HR=", substr(round(HR, 3), start=1, stop=4), ", p=", p_value)),
             fill = "white", color="black") +
  theme_bw() +
  xlab("") +
  ylab("Hazard ratio of MMR in 2 years") +
  scale_color_brewer(palette = "Set1") +
  ylim(0,2) +
  theme(axis.text = element_text(size=11, colour = "black"),
        axis.title=element_text(size=12, colour = "black"),
        legend.position = "None",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 6, 5.5, 5),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_flip()

# C-index Barplot (Right)
g2 = ggplot(data = c_index_df, aes(x = covariates, y = change)) +
  geom_linerange(aes(ymin = change, ymax = 0), color = "steelblue", size = 1.5) +
  geom_point(aes(y = change), size = 5, shape = 18, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label(aes(x = covariates, y = change, label = round(change, 3)),
             fill = "white", color = "black", vjust = -0.6) +
  theme_bw() +
  xlab("") +
  ylab("ΔC-index") +
  scale_color_brewer(palette = "Set1") +
  ylim(-0.04,0) +
  theme(axis.text = element_text(size=11, colour = "black"),
        axis.title=element_text(size=12, colour = "black"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 6, 5.5, 5),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_flip()

# Combine plots tightly and draw a crisp separator on the canvas
combined_plot = plot_grid(g1, g2, ncol = 2, rel_widths = c(1.5, 1), align = "h")

print(combined_plot)

ggsave(plot = combined_plot, filename = paste0(results, "_response/Barplot_coefficients_with_c_indexes.png"),
       width = 8, height = 3.5, units = "in", dpi = 300)
