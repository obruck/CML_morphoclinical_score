rm(list = ls())

# Load necessary libraries
source("~/mounts/research/src/Rfunctions/library.R")

# General parameters
source("~/mounts/research/husdatalake/disease/scripts/CML/R/parameters")

dir.create(paste0(results, "_response"))



# Clinical data

# Abroad
abroad = readxl::read_xlsx("mounts/research/husdatalake/disease/processed_data/CML/CML_clinical_data_v1.15.xlsx") %>%
  dplyr::rename(henkilotunnus = `Subject ID`,
                album_id = `Hemavision ID`,
                first_TKI = `First-line tyrosine kinase inhibitor generic name before discontinuation`) %>%
  dplyr::mutate(album_id = as.integer(album_id),
                first_TKI = tolower(first_TKI))

## Remove unnecessary data
columns_to_drop <- grep("(?i)(discontinuation|cytogenetic|BCR-ABL|pauses|last |12 months|relapse|original|Major|comments|trial|transcript level|mutations|year of|Disease stage|Organ|sample type)", names(abroad), value = TRUE)
abroad = abroad %>%
  dplyr::select(-all_of(columns_to_drop))#%>%

## Rename variables
abroad = abroad %>%
  dplyr::rename(age_at_dg = "Patient age at diagnosis (completed years)",
                gender = "Patient gender",
                first_TKI_class = "First-line tyrosine kinase inhibitor is second-generation TKI",
                IFNa = "Usage of IFNa (dose, duration)",
                MMR_time = "Time to MMR (months)",
                MR4.0_time = "Time to MR4.0 (months)",
                MR4.5_time = "Time to MR4.5 (months)",
                spleen_size = "Spleen size at diagnosis (max. distance from costal margin)",
                b_leuk = "PB white blood cells (E9/L)",
                l_baso = "PB basophils at diagnosis (%)",
                b_trom = "PB platelet count at diagnosis (E9/L)",
                l_eos = "PB eosinophils at diagnosis (%)",
                l_lymf = "PB lymphocytes at diagnosis (%)",
                l_blast = "PB blasts at diagnosis (%)",
                bm_blast = "BM blasts at diagnosis (%)",
                Sokal_class = "Sokal score",
                Hasford_class = "Hasford score",
                EUTOS_class = "EUTOS score",
                ELTS_class = "ELTS score") %>%
  dplyr::mutate(age_at_dg = as.numeric(age_at_dg),
                gender = ifelse(gender == "female", 1, 0),
                IFNa = ifelse(!is.na(IFNa), FALSE, TRUE),
                spleen_size = as.numeric(spleen_size),
                first_TKI_class = ifelse(first_TKI_class == "yes", 0, 1),
                l_baso = round(as.numeric(l_baso), 1),
                l_eos = round(as.numeric(l_eos), 1),
                l_lymf = round(as.numeric(l_lymf), 1),
                l_blast = round(as.numeric(l_blast), 1),
                bm_blast = round(as.numeric(bm_blast), 1),
                spleen_size = as.numeric(spleen_size),
                MMR = ifelse(is.na(MMR_time), 0, 1),
                MR4.0 = ifelse(is.na(MR4.0_time), 0, 1),
                MR4.5 = ifelse(is.na(MR4.5_time), 0, 1),
                MMR_time = MMR_time * 30.4,
                MR4.0_time = MR4.0_time * 30.4,
                MR4.5_time = MR4.5_time * 30.4,
                Sokal_class = ifelse(is.na(Sokal_class), NA,
                                     ifelse(Sokal_class < 0.8, "Low",
                                            ifelse(Sokal_class <= 1.2, "Intermediate", "High"))),
                Hasford_class = ifelse(is.na(Hasford_class), NA,
                                        ifelse(Hasford_class <= 780, "Low",
                                               ifelse(Hasford_class < 1480, "Intermediate", "High"))),
                EUTOS_class = ifelse(is.na(EUTOS_class), NA,
                                     ifelse(EUTOS_class <= 87, "Low", "High")),
                ELTS_class = ifelse(is.na(ELTS_class), NA,
                                    ifelse(ELTS_class <= 1.5680, "Low",
                                           ifelse(ELTS_class <= 2.2185, "Intermediate", "High"))))

# HUS
## MGG samples
hus = fread(paste0(export, "/mgg_samples_cml_dg.csv"))

## Demo
demo = readRDS("mounts/research/husdatalake/disease/processed_data/CML/demographics.rds") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::select(henkilotunnus, dg_date_combined, sukupuoli_selite, age_at_dg, first_TKI_manual)

## Risk classes
cml_pat = fread("mounts/research/husdatalake/disease/processed_data/CML/CML_pathology.csv") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::select(henkilotunnus, ends_with("_class"))

## Responses
response = fread("mounts/research/husdatalake/disease/processed_data/CML/CML_response.csv") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::select(henkilotunnus, starts_with("MMR"), starts_with("MR"), ends_with("_time"), -starts_with("CCyR"))

## Spleen size
spleen1 = fread("mounts/research/husdatalake/data/radiology/csv/spleen_df.csv") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::left_join(demo %>% dplyr::select(henkilotunnus, dg_date_combined)) %>%
  dplyr::mutate(time_diff = abs(as.Date(dg_date_combined) - as.Date(tutkimus_pvm))) %>%
  group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::filter(time_diff < 90) %>%
  dplyr::select(henkilotunnus, pernan_mitta)
spleen2 = fread("mounts/research/husdatalake/disease/processed_data/CML/CML_spleen.csv") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::mutate(time_diff = abs(as.Date(dg_date_combined) - as.Date(pvm))) %>%
  group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::filter(time_diff < 90) %>%
  dplyr::select(henkilotunnus, pernan_mitta)
spleen = bind_rows(spleen1, spleen2) %>%
  group_by(henkilotunnus) %>%
  summarise(pernan_mitta = max(pernan_mitta, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(pernan_mitta = ifelse(pernan_mitta < 11, 0, pernan_mitta - 11)) %>%
  dplyr::rename(spleen_size = pernan_mitta)
rm(spleen1, spleen2)


## Lab
lab = readRDS("mounts/research/husdatalake/disease/processed_data/CML/lab_dg.rds") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::filter(tutkimus_lyhenne %in% c("B -Erblast", "B -Trom", "B -Eos", "B -Baso", "B -Lymf", "B -Neut", "B -Mono", "B -Leuk",
                                        "L -Neut", "L -Mono", "L -Eos", "B -Leuk", "L -Blast", "L -Baso", "L -Lymf")) %>%
  group_by(henkilotunnus) %>%
  mutate(leuk = ifelse(tutkimus_lyhenne == "B -Leuk", tulos_max, NA)) %>%
  fill(leuk, .direction = "downup") %>%
  ungroup() %>%
  dplyr::mutate(tulos_max = ifelse(tutkimus_lyhenne %in% c("B -Erblast", "B -Eos", "B -Baso", "B -Lymf", "B -Neut", "B -Mono"), round(100*tulos_max/leuk), tulos_max),
                tutkimus_lyhenne = ifelse(tutkimus_lyhenne %in% c("B -Erblast", "B -Eos", "B -Baso", "B -Lymf", "B -Neut", "B -Mono"), gsub("B ", "L ", tutkimus_lyhenne), tutkimus_lyhenne)) %>%
  group_by(henkilotunnus, tutkimus_lyhenne) %>%
  summarise(tulos_max = max(tulos_max, na.rm = TRUE)) %>%
  ungroup() %>%
  # dplyr::mutate(tulos = ifelse(tutkimus_lyhenne %in% c("b_leuk", "b_trom", "l_baso", "l_eos", "l_blast"), tulos_max, tulos_mean)) %>%
  # dplyr::select(-c(tulos_min, tulos_mean)) %>%
  reshape2::dcast(henkilotunnus~tutkimus_lyhenne, value.var = "tulos_max") %>%
  janitor::clean_names() %>%
  dplyr::mutate(l_blast = ifelse(is.na(l_blast), 0, l_blast))

## Treatment
treatment = fread("mounts/research/husdatalake/disease/processed_data/CML/CML_treatment.csv") %>%
  dplyr::filter(henkilotunnus %in% hus$henkilotunnus) %>%
  dplyr::select(henkilotunnus, IFNa, first_TKI)

## Blasts
## MGG
pat1 = readRDS(paste0(import, "/pathology_myplus.rds")) %>%
  dplyr::filter(henkilotunnus %in% demo$henkilotunnus)
pat2 = readRDS(paste0(import, "/pathology_qpati.rds")) %>%
  dplyr::filter(henkilotunnus %in% demo$henkilotunnus)
pathology = rbind(pat1 %>% distinct(henkilotunnus, samplenumber, sampletaken, examinations, statement, diagnose), pat2 %>% distinct(henkilotunnus, samplenumber, sampletaken, examinations, statement, diagnose)); rm(pat1, pat2)
mgg = readRDS(paste0(import, "/mgg.rds")) %>%
  dplyr::left_join(pathology) %>%
  dplyr::inner_join(demo)
mgg1 = mgg %>%
  dplyr::filter(mgg_type_string == "blast") %>%
  dplyr::mutate(time_diff = abs(as.Date(sampletaken) - as.Date(dg_date_combined)),
                mgg_section1 = ifelse(mgg_section == "bone_marrow", 1, ifelse(mgg_section == "summary", 2, 3))) %>%
  dplyr::arrange(henkilotunnus, time_diff, mgg_section1) %>%
  group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(henkilotunnus, bm_blast_mgg = mgg_value_percentage)
## BM-PAD
pad = pathology %>%
  dplyr::filter(str_detect(examinations, "(?i)(Bm-PAD|kri|cri|Lu-PAD)")) %>%
  dplyr::inner_join(demo) %>%
  dplyr::mutate(time_diff = abs(as.Date(sampletaken) - as.Date(dg_date_combined))) %>%
  dplyr::arrange(henkilotunnus, time_diff) %>%
  group_by(henkilotunnus, sampletaken) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(bm_blast_pad = gsub(" ", "",
                             gsub(",", ".",
                                  gsub("%", "",
                                       str_extract(diagnose, "[[:digit:]][[:print:]]{0,8}%")))),
         bm_blast_pad0 = str_extract(statement, "(?i)(blast|kantaso)[^.]*%"),
         bm_blast_pad0 = gsub(",", ".",
                              gsub(" ", "",
                                   gsub("%", "", str_extract(bm_blast_pad0, "[[:digit:]][[:print:]]{0,8}%")))),
         bm_blast_pad = ifelse(!is.na(bm_blast_pad), bm_blast_pad, bm_blast_pad0), 
         bm_blast_pad1 = ifelse(str_detect(bm_blast_pad, "-"), gsub("-[[:print:]]*", "", bm_blast_pad), bm_blast_pad),
         bm_blast_pad2 = ifelse(str_detect(bm_blast_pad, "-"), gsub("[[:print:]]*-", "", bm_blast_pad), bm_blast_pad),
         bm_blast_pad = ifelse(!str_detect(bm_blast_pad, "-"), bm_blast_pad, (as.numeric(bm_blast_pad1)+as.numeric(bm_blast_pad2))/2),
         bm_blast_pad = ifelse(!is.na(bm_blast_pad), bm_blast_pad,
                               ifelse(str_detect(statement, "(?i)(ei |eikä|ilman)[^.]{0,30}blast") | str_detect(statement, "(?i)blast[^.]{0,30}alle"), 3, NA)),
         bm_blast_pad = as.numeric(bm_blast_pad)
  ) %>%
  dplyr::select(-c(bm_blast_pad1, bm_blast_pad2)) %>%
  # dplyr::select(henkilotunnus, statement, diagnose, bm_blast_pad)
  dplyr::group_by(henkilotunnus) %>%
  dplyr::summarise(bm_blast_pad = first(bm_blast_pad, na.rm = TRUE),
                   bm_blast_pad_date = first(sampletaken, na.rm = TRUE),
                   time_diff = first(time_diff, na.rm=TRUE)) %>%
  dplyr::filter(time_diff < 100) %>%
  mutate(bm_blast_pad = ifelse(is.na(bm_blast_pad), 3, bm_blast_pad)) %>%
  dplyr::select(henkilotunnus, bm_blast_pad)
bmblast = full_join(mgg1, pad) %>%
  dplyr::mutate(bm_blast = ifelse(!is.na(bm_blast_mgg), bm_blast_mgg, bm_blast_pad)) %>%
  dplyr::select(henkilotunnus, bm_blast)
rm(mgg, mgg1, pad, pathology)

# Join
hus1 = hus %>%
  dplyr::inner_join(demo %>%
                     dplyr::select(-c(dg_date_combined)) %>%
                     dplyr::mutate(age_at_dg = as.numeric(age_at_dg))) %>%
  dplyr::rename(gender = sukupuoli_selite) %>%
  dplyr::mutate(gender = ifelse(gender == "Nainen", 1, 0)) %>%
  dplyr::left_join(treatment) %>%
  dplyr::mutate(first_TKI = ifelse(!is.na(first_TKI_manual), first_TKI_manual, first_TKI),
                first_TKI = ifelse(is.na(first_TKI) | first_TKI == "", "imatinib", gsub("ibi$", "ib", first_TKI)),
                first_TKI_class = ifelse(first_TKI == "imatinib", 0, 1)) %>%
  dplyr::select(-first_TKI_manual) %>%
  dplyr::left_join(cml_pat) %>%
  dplyr::left_join(spleen) %>%
  dplyr::left_join(response) %>%
  dplyr::left_join(lab) %>%
  dplyr::left_join(bmblast) %>%
  dplyr::select(-c(samplenumber, sampletaken, anamnesis, diagnose, statement, dg_date_combined, time_diff))



# Image data
## 100x cell statistics
cml100x_statistics_hus = readxl::read_xlsx("mounts/research/import/image/CML/CML_HUS_hematoscope_100x_sample_statistics_09-10-2024.xlsx") %>%
  dplyr::rename(album_id = Sample) %>%
  dplyr::mutate(album_id = as.integer(album_id))
cml100x_statistics_abroad = readxl::read_xlsx("mounts/research/import/image/CML/CML_multi-center_hematoscope_100x_sample_statistics_09-10-2024.xlsx") %>%
  dplyr::rename(album_id = `Hemavision ID`) %>%
  dplyr::mutate(album_id = as.integer(album_id))
## 100x cell counts
cml100x_cells_hus = readxl::read_xlsx("mounts/research/import/image/CML/CML_HUS_hematoscope_100x_differential_counts_09-10-2024.xlsx") %>%
  dplyr::rename(album_id = Sample) %>%
  dplyr::mutate(album_id = as.integer(album_id))
cml100x_cells_abroad = readxl::read_xlsx("mounts/research/import/image/CML/CML_multi-center_hematoscope_100x_cell_differential_counts_09-10-2024.xlsx") %>%
  dplyr::rename(album_id = `Hemavision ID`) %>%
  dplyr::mutate(album_id = as.integer(album_id))
## 10x statistics
cml10x_hus = readxl::read_xlsx("mounts/research/import/image/CML/CML_HUS_hematoscope_10x_wsi_results_18-09-2024.xlsx") %>%
  dplyr::rename(album_id = ID,
                database = Dataset) %>%
  dplyr::mutate(album_id = as.integer(album_id),
                database = gsub("_csc|_06\\.07\\.21", "", database),
                database = gsub("vimages2", "vimages", database))
cml10x_abroad = readxl::read_xlsx("mounts/research/import/image/CML/CML_Abroad_hematoscope_10x_wsi_results_31-05-2024.xlsx") %>%
  dplyr::rename(album_id = `Hemavision ID`) %>%
  dplyr::mutate(album_id = as.integer(album_id))

unique(cml10x_hus$database)[!unique(cml10x_hus$database) %in% unique(cml100x_cells_hus$database)]
unique(cml10x_hus$database)[!unique(cml10x_hus$database) %in% unique(cml100x_statistics_hus$database)]
unique(cml100x_statistics_hus$database)[!unique(cml100x_statistics_hus$database) %in% unique(cml10x_hus$database)]
unique(cml100x_cells_hus$database)[!unique(cml100x_cells_hus$database) %in% unique(cml10x_hus$database)]

# Join
abroad = abroad %>%
  dplyr::left_join(cml10x_abroad) %>%
  dplyr::left_join(cml100x_cells_abroad) %>%
  dplyr::left_join(cml100x_statistics_abroad)
hus1 = hus1 %>%
  dplyr::left_join(cml10x_hus) %>%
  dplyr::left_join(cml100x_cells_hus) %>%
  dplyr::left_join(cml100x_statistics_hus) %>%
  dplyr::mutate(Center="HUS")
df = full_join(hus1, abroad)

# Remove ids with no image data in both 10x and 100x
rem_ids =  which(is.na(df$Megakaryocyte_Center_Area))[which(is.na(df$Megakaryocyte_Center_Area)) %in% which(is.na(df$Neutrophils))]
df = df[-rem_ids]

# Verify that columns are identical between HUS and abroad
names(cml10x_abroad)[!names(cml10x_abroad) %in% names(cml10x_hus)]
names(cml100x_cells_abroad)[!names(cml100x_cells_abroad) %in% names(cml100x_cells_hus)]
names(cml100x_statistics_abroad)[!names(cml100x_statistics_abroad) %in% names(cml100x_statistics_hus)]

names(cml10x_hus)[!names(cml10x_hus) %in% names(cml10x_abroad)]
names(cml100x_cells_hus)[!names(cml100x_cells_hus) %in% names(cml100x_cells_abroad)]
names(cml100x_statistics_hus)[!names(cml100x_statistics_hus) %in% names(cml100x_statistics_abroad)]

# Remove unnecessary variables at 10x
columns_to_drop <- grep("(?i)(^periphery|Periphery_Area$|focus|empty|small_|large_|medium_|sample_area$|center_area$)", names(df), value = TRUE)
df <- df %>% dplyr::select(-all_of(columns_to_drop), -Sample_Proportion, -Center_Proportion)
## Discard non-percentage variables (100x)
df = df %>%
  dplyr::select(-all_of(gsub("_Percentage", "", names(df)[grep(pattern = "Percentage", names(df), ignore.case = TRUE)])))
## Other redundant 100x variables
columns_to_drop <- grep("(?i)(broken|unclear|artefac|all cells|unknown|^living cells_percentage)", names(df), value = TRUE)
df = df %>%
  dplyr::select(-all_of(columns_to_drop))

# Save
saveRDS(df, paste0(export, "_response/data_for_modelling1.rds"))
fwrite(df, paste0(export, "_response/data_for_modelling1.csv"))
