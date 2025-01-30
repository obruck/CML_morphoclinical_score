rm(list = ls())

# Load necessary libraries
source("~/library.R")
library(uwot)



###### Process data ########


# Data
df0 = readRDS("./data_for_modelling_all_cml1.rds") %>%
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
    HUS_or_AUS_pt = ifelse(str_detect(henkilotunnus, '^(02139|AUS)'), 1, 0)) %>%
  mutate(across(c(dasatinib, imatinib, nilotinib), as.integer))

# Rename
df0 = df0 %>%
  janitor::clean_names()


## Replace NaN and Inf with NA
df0[df0 == "NaN"] <- NA
df0[df0 == "Inf"] <- NA
df0[df0 == "-Inf"] <- NA

df = df0


###### Split data ########


# Set seed
seed1=10
set.seed(seed1)

# Select features for UMAP
umap_data = df %>%
  dplyr::select(which(names(df)=="rbc_sample_proportion"):which(names(df)=="macrophages_nuclei_compactness_percentile_75"))
umap_metadata = df
## Remove rows where all image values are NA
remove_rows = rowSums(is.na(umap_data)) != ncol(umap_data)
umap_data = umap_data[remove_rows,]
umap_metadata = umap_metadata[remove_rows,]
## Fill remaining NA with median
umap_data = umap_data %>% 
  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))

# Scale
image_10x = umap_data %>%
  dplyr::select(which(names(umap_data)=="rbc_sample_proportion"):which(names(umap_data)=="lipid_total_periphery_proportion"))
image_100x_cells = umap_data %>%
  dplyr::select(!which(names(umap_data)=="rbc_sample_proportion"):which(names(umap_data)=="lipid_total_periphery_proportion")) %>%
  dplyr::select(contains("percentage")) %>%
  dplyr::select(-c(contains("_class"), contains("_class"), contains("quality"), contains("granularity"), contains("vacuolated"), contains("multinuclear")))
image_100x_cellstatistics = umap_data %>%
  dplyr::select(!which(names(umap_data)=="rbc_sample_proportion"):which(names(umap_data)=="lipid_total_periphery_proportion")) %>%
  dplyr::select(-names(image_100x_cells))

image_10x_scaled = sapply(image_10x, function(x) 2*(x - mean(x)) / sd(x))
image_100x_cells_scaled = sapply(image_100x_cells, function(x) 10*(x - mean(x)) / sd(x))
image_100x_cellstatistics_scaled = sapply(image_100x_cellstatistics, function(x) (x - mean(x)) / sd(x))
## Join
image_100x_cells_scaled = as.data.frame(image_100x_cells_scaled)
image_100x_cells_scaled$blasts_percentage = image_100x_cells_scaled$blasts_percentage*10; umap_data1 = image_100x_cells_scaled
umap_data1 = cbind(cbind(image_10x_scaled, image_100x_cells_scaled), image_100x_cellstatistics_scaled)


# UMAP
umap1 = umap(umap_data1, scale = FALSE, n_neighbors = 15, seed = 10, metric = "cosine") %>% # metric = "correlation", "cosine"
  as.data.frame()

# Add ids
umap1 = umap1 %>%
  dplyr::rename(UMAP_1 = V1, UMAP_2 = V2) %>%
  cbind(umap_metadata)

# Plot
## BM Blast %
umap1$b_leuk1 = ifelse(abs(umap1$bm_blast) > 20, 20, umap1$bm_blast)
umap1$order = ifelse(is.na(umap1$b_leuk1), 0, umap1$b_leuk1)
g = ggplot(data = umap1 %>%
             dplyr::arrange(order)) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = b_leuk1), size = 2, alpha=1) +
  scale_color_distiller(palette = "RdBu", name = "BM Blast (%)") + 
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    legend.position = "none",
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.background=element_blank()); g
ggsave(plot = g, filename = "./UMAP_bmblast.png", bg = "white", width = 6, height = 6, units = "in", dpi = 300)

## WBC
umap1$b_leuk1 = ifelse(abs(umap1$b_leuk) > 100, 100, umap1$b_leuk) 
umap1$order = ifelse(is.na(umap1$b_leuk1), 0, umap1$b_leuk1)
g = ggplot(data = umap1 %>%
             dplyr::arrange(order)) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = b_leuk1), size = 2, alpha=1) +#, size = 1) +
  scale_color_distiller(palette = "RdBu", name = "PB WBC (E9/L)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        axis.title.y=element_blank(),
        plot.background=element_blank()); g
ggsave(plot = g, filename = "./UMAP_wbc.png", bg = "white", width = 6, height = 6, units = "in", dpi = 300)

## WBC
umap1 = umap1 %>%
  dplyr::mutate(rown = 1:nrow(.))
umap2 = umap1 %>%
  dplyr::arrange(mgg_time_diff) %>%
  dplyr::filter(abs(mgg_time_diff) <= 60) %>%
  dplyr::group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup()
umap1$b_leuk1 = ifelse(umap1$rown %in% umap2$rown, "DG", "FU") 
umap1$b_leuk1 = ifelse(is.na(umap1$b_leuk1), "FU", umap1$b_leuk1)
umap1$order = ifelse(is.na(umap1$b_leuk1), 0, umap1$b_leuk1)
g = ggplot(data = umap1 %>%
             dplyr::arrange(order)) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = b_leuk1), size = 2, alpha=1) +#, size = 1) +
  scale_color_brewer(palette = "Set1", name = "Time from\ndiagnosis\n(days)") + 
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "none",
    plot.background=element_blank()); g
ggsave(plot = g, filename = "./UMAP_baseline_followup.png", bg = "white", width = 6, height = 6, units = "in", dpi = 300)

