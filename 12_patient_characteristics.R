rm(list = ls())

# Load necessary libraries
source("~/library.R")


##### Descriptive ##### 


# Read data
df0 = readRDS("./data_for_modelling_all_cml1.rds") %>%
  dplyr::arrange(mgg_time_diff) %>%
  dplyr::filter(abs(mgg_time_diff) <= 60) %>%
  dplyr::group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(
    imatinib = ifelse(first_TKI=="imatinib", TRUE, FALSE),
    nilotinib = ifelse(first_TKI=="nilotinib", TRUE, FALSE),
    dasatinib = ifelse(first_TKI=="dasatinib", TRUE, FALSE),
    sec_gen_tki = ifelse(imatinib == TRUE, FALSE, TRUE),
    first_TKI_class = ifelse(is.na(first_TKI), NA,
                             ifelse(first_TKI == "imatinib", TRUE, FALSE)),
    MR4.0 = ifelse(MR4.0 == 1, TRUE, FALSE),
    MMR = ifelse(MMR == 1, TRUE, FALSE),
    MMR_time = as.numeric(MMR_time)/365.24,
    MR4.0_time = as.numeric(MR4.0_time)/365.24,
    MR4.5_time = as.numeric(MR4.5_time)/365.24,
    response_at_24 = ifelse(!is.na(MMR) & MMR_time<2 & MMR == 0, NA,
                            ifelse(MMR_time >= 2, FALSE, TRUE)),
    HUS_or_AUS_pt = ifelse(str_detect(henkilotunnus, '^(02139|AUS)'), 1, 0),
    Cohort = "Descriptive") %>%
  mutate(across(c(dasatinib, imatinib, nilotinib), as.integer))

# Rename
df0 = df0 %>%
  janitor::clean_names()

## Replace NaN and Inf with NA
df0[df0 == "NaN"] <- NA
df0[df0 == "Inf"] <- NA
df0[df0 == "-Inf"] <- NA

df_descriptive = df0


##### Modelling ##### 


# Read data
df0 = readRDS("./data_for_modelling.rds") %>%
  dplyr::filter(!(is.na(MMR_time) | MMR_time == 0)) %>%
  mutate(
    imatinib = ifelse(first_TKI=="imatinib", TRUE, FALSE),
    nilotinib = ifelse(first_TKI=="nilotinib", TRUE, FALSE),
    dasatinib = ifelse(first_TKI=="dasatinib", TRUE, FALSE),
    sec_gen_tki = ifelse(imatinib == TRUE, FALSE, TRUE),
    first_TKI_class = ifelse(is.na(first_TKI), NA,
                             ifelse(first_TKI == "imatinib", TRUE, FALSE)),
    MR4.0 = ifelse(MR4.0 == 1, TRUE, FALSE),
    MMR = ifelse(MMR == 1, TRUE, FALSE),
    MMR_time = as.numeric(MMR_time)/365.24,
    MR4.0_time = as.numeric(MR4.0_time)/365.24,
    MR4.5_time = as.numeric(MR4.5_time)/365.24,
    response_at_24 = ifelse(!is.na(MMR) & MMR_time<2 & MMR == 0, NA,
                            ifelse(MMR_time >= 2, FALSE, TRUE)),
    HUS_or_AUS_pt = ifelse(str_detect(henkilotunnus, '^(02139|AUS)'), 1, 0),
    Cohort = "Modelling") %>%
  mutate(across(c(dasatinib, imatinib, nilotinib), as.integer))


# Rename
df0 = df0 %>%
  janitor::clean_names()


## Replace NaN and Inf with NA
df0[df0 == "NaN"] <- NA
df0[df0 == "Inf"] <- NA
df0[df0 == "-Inf"] <- NA

df_modelling = df0


##### Test #3 ##### 


# Read data
df = read_xlsx("./Helsinki Full TFR cohort data _de-identified.xlsx", sheet = "CML Full TFR cohort data")
df1 = read_xlsx("./Helsinki Full TFR cohort data _de-identified.xlsx", sheet = "TFR cohort BM data")
df_test = df %>%
  dplyr::left_join(df1) %>%
  janitor::clean_names() %>%
  dplyr::filter(!id_number == "12") %>%
  dplyr::mutate(mmr_time = (as.Date(date_of_0_1_percent) - as.Date(date_of_tki_commencement))/365.24,
                ccyr_time = (as.Date(date_of_1_percent) - as.Date(date_of_tki_commencement))/365.24)


##### FINAL PROCESSING #####


# Age, gender, tki, sokal, Hasford, elts, eutos
df_descriptive = df_descriptive %>%
  dplyr::select(cohort, gender, age_at_dg, first_tki, ends_with("class"), mmr, mmr_time, response_at_24, mr4_0, mr4_0_time, b_leuk) %>%
  dplyr::mutate(gender = ifelse(gender==1, "Female", "Male"))
df_modelling = df_modelling %>%
  dplyr::select(cohort, gender, age_at_dg, first_tki, ends_with("class"), mmr, mmr_time, response_at_24, mr4_0, mr4_0_time, b_leuk) %>%
  dplyr::mutate(gender = ifelse(gender==1, "Female", "Male"))
df_test = df_test %>%
  dplyr::select(gender, age_at_dg = age_at_diagnosis_yrs, first_tki, mmr_time,
                mr4_0_time=time_to_mr4_yrs, sokal_class=sokal_score, hasford_class=hasford, eutos_class=eutos,
                elts_class=elts, mr4_0 = mr4_achievement_0_no_1_yes_2_never, b_leuk = pb_white_cell_count_e9_l,
  ) %>%
  dplyr::mutate(cohort = "Test #3",
                mmr = TRUE,
                mmr_time = as.numeric(mmr_time),
                first_tki = gsub(" [[:print:]]*", "", first_tki),
                first_tki = ifelse(is.na(first_tki), NA,
                                   ifelse(str_detect(first_tki, "(?i)IM"), "imatinib",
                                          ifelse(str_detect(first_tki, "(?i)NIL"), "nilotinib",
                                                 ifelse(str_detect(first_tki, "(?i)DAS"), "dasatinib", first_tki)))),
                first_tki_class = ifelse(is.na(first_tki), NA,
                                         ifelse(first_tki=="imatinib", TRUE, FALSE)),
                elts_class = ifelse(is.na(elts_class), NA,
                                    ifelse(elts_class<=1.5680, "Low",
                                           ifelse(elts_class>2.2185, "High", "Intermediate"))),
                sokal_class = ifelse(is.na(sokal_class), NA,
                                     ifelse(sokal_class<0.8, "Low",
                                            ifelse(sokal_class>1.2, "High", "Intermediate"))),
                hasford_class = ifelse(is.na(hasford_class), NA,
                                       ifelse(hasford_class<=780, "Low",
                                              ifelse(hasford_class>1480, "High", "Intermediate"))),
                eutos_class = ifelse(is.na(eutos_class), NA,
                                     ifelse(eutos_class<=87, "Low", "High")),
                response_at_24 = ifelse(!is.na(mmr) & mmr_time<2 & mmr == 0, NA,
                                        ifelse(mmr_time >= 2, FALSE, TRUE)),
  )


##### Patient characteristics #####


# Bind rows
df = bind_rows(df_modelling, as.data.frame(df_test)) %>%
  dplyr::select(-first_tki_class)
df = df %>%
  dplyr::mutate(first_tki = ifelse(str_detect(first_tki, "(?i)(bosu|pona)"), "Bosutinib/Ponatinib", str_to_sentence(first_tki)))

# Table
table1 <-
  tbl_summary(
    df,
    statistic = list(
      all_continuous() ~ "{median} [{p25}-{p75}]",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    by = cohort, # split table by group
    digits = all_continuous() ~ 2,
    label = c(gender ~ "Gender",
              age_at_dg ~ "Age at diagnosis",
              first_tki ~ "First TKI",
              sokal_class ~ "Sokal",
              hasford_class ~ "Hasford",
              eutos_class ~ "EUTOS",
              elts_class ~ "ELTS",
              mmr ~ "MMR",
              mmr_time ~ "MMR time",
              response_at_24 ~ "Response at 24 months",
              mr4_0 ~ "MR4.0",
              mr4_0_time ~ "MR4.0 time",
              b_leuk ~ "PB WBC at diagnosis (E9/L)"),
    missing_text = "(Missing)"
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(
    list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test")
  ) %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()
table1
table1 %>% as_gt() %>%
  gt::gtsave(filename = "./Patient_table.png") # use extensions .png, .html, .docx, .rtf, .tex, .ltx
table1 %>% as_gt() %>%
  gt::gtsave(filename = "./Patient_table.docx")
table1 %>% as_gt() %>%
  gt::gtsave(filename = "./Patient_table.html")
