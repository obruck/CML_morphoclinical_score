# Validation with the AUS set #2


# Library
source("~/library.R")
library(ggpubr)

# Read data
res = readxl::read_xlsx("./univariate_cox_results.xlsx")

# Prepare data for volcanoplot
res1 <- res %>%
  dplyr::filter(!str_detect(names, "(?i)(plasma|lympho|promonocytes|green|red|blue|solidity|compactness|macro)")) %>%
  dplyr::arrange(p.value) %>%
  mutate(
    fold.change.log10 = log(HR, 10),
    sig = ifelse(p.value>=0.05, "ns",
                 ifelse(HR >= 1, "HR≥1", "HR<1")),
    p.value1 = ifelse(HR > 1, p.value, -p.value),
    p.value2 = ifelse(HR > 1, -log(1/p.value, 10), log(1/p.value, 10)),
    order1 = 1:nrow(.),
    names = gsub("_", " ",
                 gsub("^blasts", "Blasts",
                      gsub("^eryt", "Eryt",
                           gsub("^spleen", "Spleen",
                                gsub("elts_class", "ELTS",
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
                                                                                                           gsub("dasatinib", "Dasatinib (yes/no)", names )))))))))))))))))))))

# Edit the point size
res1 = res1 %>%
  dplyr::mutate(p.value3 = ifelse(p.value < 0.05, 0.00001, p.value))

# Plot
g = ggplot(res1, aes(order1, p.value2)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
  annotate(label="p=0.05", geom = "text", x = 550, y = -log10(0.035), color = "black", size = 4) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
  annotate(label="p=0.05", geom = "text", x = 550, y = log10(0.035), color = "black", size = 4) +
  geom_point(aes(fill=sig, size=1/p.value3), shape = 21, color="black") +
  labs(x="Rank", y="P-value") +
  guides(fill=guide_legend("P-value",
                           nrow = 1,
                           # title.position = "top",
                           title.hjust = 0.5,
                           override.aes = list(size = 5)),
         size = FALSE) +
  theme_bw() +
  theme(legend.direction = 'horizontal', 
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=12, face="bold", colour = "black"),
        legend.title=element_text(size=12, face="bold", colour = "black"),
        legend.text=element_text(size=12, colour = "black"),
        legend.position = "bottom") +
  scale_fill_manual(values=c("blue", "red", "black")) +
  scale_y_continuous(breaks = c(-2:2), labels=c(0.01,0.1,0,0.1,0.01)) +
  geom_label_repel(data=res1 %>%
                     dplyr::slice(1:10) %>%
                     dplyr::filter(!names == "Monocytes nuclei area (median)"), size=5, aes(label=names),
                   box.padding = 0.7, max.overlaps=Inf, force=20, force_pull=20);g
# Remove manually collinear variables (Monocytes nuclei area (median))
ggsave(plot = g, filename = "./volcanoplot.png", width = 8, height = 5.5, units = "in", dpi = 300)
