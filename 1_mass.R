rm(list = ls())


# Load necessary libraries
source("~/library.R")
library(survminer)
library(caret)
library(pROC)
library(MLmetrics)
library(PRROC)
library(cutpointr)
library(MASS)
library(cowplot)
library(patchwork)


###### Process data ########


# Data
df0 = readRDS("./data_for_modelling.rds") %>%
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


# Rename
df0 = df0 %>%
  janitor::clean_names()

# Define endpoints
df0 = df0 %>%
  dplyr::mutate(
    mmr_time = ifelse(is.na(mmr_time), 0, mmr_time),
    mmr_time1=mmr_time,
    mmr = response_at_24
  ) %>%
  dplyr::filter(!is.na(response_at_24)) %>%
  dplyr::filter(!(response_at_24==FALSE & mmr_time1<24))


# Divide some variables with 10 to make the HR more interpretable
df0 = df0 %>%
  dplyr::mutate(b_leuk = b_leuk / 10) %>%
  dplyr::mutate_at(vars(matches('perimeter')), function(x) x = x/10)


## Replace NaN and Inf with NA
df0[df0 == "NaN"] <- NA
df0[df0 == "Inf"] <- NA
df0[df0 == "-Inf"] <- NA


# Store
df = df0

# Export
writexl::write_xlsx(df, "./full_data.xlsx")

# Remove unnecessary variables
df = df %>%
  dplyr::select(-c(first_tki, mr4_0, mr4_5, mr4_0_time, mr4_5_time))


###### Split data ########


# Set seed
seed1=10
set.seed(seed1)


# Split data
df_train1 <- df %>%
  dplyr::filter(!center %in% c("Australia"))
train_index <- createDataPartition(y = df_train1$mmr, p = 0.8, list = FALSE) # df$hus_or_aus_pt
df_train <- df_train1[train_index, ]
df_val <- df_train1[-train_index, ]
df_test <- df %>%
  dplyr::filter(center %in% c("Australia"))


####### Select features for Cox model ########


# Univariate Cox regression to select features
covariates <- names(df_train)[!(names(df_train) %in% c("henkilotunnus", "database", "album_id", "mmr_time", "mmr", "hus_or_aus_pt", "center", "mmr_time1", "response_at_24"))]


# Function to apply
## Function
cox_fun <- sapply(covariates,
                  function(x) as.formula(paste('Surv(df_train$mmr_time, df_train$response_at_24)~', x, sep="")))
## Apply
model_list <- lapply(cox_fun, function(x){coxph(x, data = df_train)})

# Extract data
univ_results <- lapply(model_list,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$coefficients[5], digits=5)
                         beta<-signif(x$coef[1], digits=5);#coeficient beta
                         HR <-signif(x$coef[2], digits=5);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 5)
                         HR.confint.upper <- signif(x$conf.int[1,"upper .95"], 5)
                         HR1 <- paste0(HR, " (",
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR1, p.value)
                         names(res)<-c("beta", "HR", "HR (95% CI for HR)", "p.value")
                         return(res)
                       })

# Results
res <- t(as.data.frame(univ_results)) %>%
  as.data.frame %>%
  tibble::rownames_to_column("names") %>%
  mutate(
    beta = as.numeric(as.character(beta)),
    HR = as.numeric(as.character(HR)),
    p.value = as.numeric(as.character(p.value)),
    p.value_adj = p.adjust(p.value, method = "BH", n = nrow(.)))

# Save
if (!exists("res1")) {
  res1 = res
} else {
  res1 = rbind(res1, res)
}
writexl::write_xlsx(res1, "./univariate_cox_results.xlsx")


# Univariate
selected_features <- res1 %>%
  dplyr::arrange(p.value) %>%
  dplyr::filter(p.value<0.05) %>%
  distinct(names)
## Covariates to remove
covariates_to_remove = covariates[!covariates %in% selected_features$names]


selected_features2 = data.frame()
for (j in 1:10) {
  
  print(j)
  
  # Round 2
  df_train1 = df_train %>%
    dplyr::slice_sample(prop = 0.5)
  
  # Function to apply
  ## Function
  cox_fun <- sapply(selected_features$names,
                    function(x) as.formula(paste('Surv(df_train1$mmr_time, df_train1$response_at_24)~', x, sep="")))
  ## Apply
  model_list <- lapply(cox_fun, function(x){coxph(x, data = df_train1)})
  
  # Extract data
  univ_results <- lapply(model_list,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$coefficients[5], digits=5)
                           beta<-signif(x$coef[1], digits=5);#coeficient beta
                           HR <-signif(x$coef[2], digits=5);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 5)
                           HR.confint.upper <- signif(x$conf.int[1,"upper .95"], 5)
                           HR1 <- paste0(HR, " (",
                                         HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, HR1, p.value)
                           names(res)<-c("beta", "HR", "HR (95% CI for HR)", "p.value")
                           return(res)
                         })
  
  # Results
  res <- t(as.data.frame(univ_results)) %>%
    as.data.frame %>%
    tibble::rownames_to_column("names") %>%
    mutate(
      beta = as.numeric(as.character(beta)),
      HR = as.numeric(as.character(HR)),
      p.value = as.numeric(as.character(p.value)),
      p.value_adj = p.adjust(p.value, method = "BH", n = nrow(.)))
  
  # Univariate
  selected_features1 <- res %>%
    dplyr::arrange(p.value) %>%
    dplyr::filter(p.value<0.05) %>%
    dplyr::mutate(loop = j)
  ## Covariates to keep
  selected_features2 = rbind(selected_features2, selected_features1)
  
}


# Proportion of cells
df %>%
  dplyr::select(dplyr::any_of(ends_with("percentage")) & !dplyr::any_of(contains(c("vacu", "multi", "mega", "class", "gran", "apl")))) %>%
  summarise_all(.funs = median, na.rm=TRUE)
df %>%
  dplyr::group_by(center) %>%
  summarise(median = median(myelocytes_cell_solidity_mean, na.rm=TRUE),
            mean = mean(myelocytes_cell_solidity_mean, na.rm=TRUE))


# Filter most common variables
selected_features3 = selected_features2 %>%
  dplyr::group_by(names) %>%
  summarise(n = n()) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(n>3); selected_features3



# Calculate HR for training set

# Function to apply
## Function
cox_fun <- sapply(selected_features3$names,
                  function(x) as.formula(paste('Surv(df_train$mmr_time, df_train$response_at_24)~', x, sep="")))
## Apply
model_list <- lapply(cox_fun, function(x){coxph(x, data = df_train)})

# Extract data
univ_results <- lapply(model_list,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$coefficients[5], digits=5)
                         beta<-signif(x$coef[1], digits=5);#coeficient beta
                         HR <-signif(x$coef[2], digits=5);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 5)
                         HR.confint.upper <- signif(x$conf.int[1,"upper .95"], 5)
                         HR1 <- paste0(HR, " (",
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR1, p.value)
                         names(res)<-c("beta", "HR", "HR (95% CI for HR)", "p.value")
                         return(res)
                       })

# Results
res <- t(as.data.frame(univ_results)) %>%
  as.data.frame %>%
  tibble::rownames_to_column("names") %>%
  mutate(
    beta = as.numeric(as.character(beta)),
    HR = as.numeric(as.character(HR)),
    p.value = as.numeric(as.character(p.value)))

# Process
selected_features4 <- res %>%
  dplyr::arrange(p.value) %>%
  dplyr::left_join(selected_features3)

# Export
writexl::write_xlsx(selected_features4, "./univariate_cox_results_top_features.xlsx")



# Forward feature selection
df_train2 = df_train
df_train2$mmr = ifelse(df_train2$mmr == FALSE, 0, 1)

# Define a function to replace NA values with the median of each column
replace_na_with_median <- function(x) {
  # Check if the column has any NA values
  if (any(is.na(x))) {
    # Replace NAs with the median, ignoring NAs when calculating the median
    x[is.na(x)] <- median(x, na.rm = TRUE)
  }
  return(x)
}
df_train2 <- as.data.frame(lapply(df_train2, replace_na_with_median))


# Define the full model with all potential predictors
top_features = selected_features3[order(selected_features3$n, decreasing = TRUE),]$names
if (!"imatinib" %in% top_features) {
  top_features = c(top_features, "imatinib")
}
if (!"nilotinib" %in% top_features) {
  top_features = c(top_features, "nilotinib")
}
if (!"elts_class" %in% top_features) {
  top_features = c(top_features, "elts_class")
}
top_features = unique(top_features)
formula_string <- paste("Surv(mmr_time, mmr) ~", paste(top_features, collapse = "+"))
cox_formula <- as.formula(formula_string)

# Fit the Cox proportional hazards model
full_model <- coxph(cox_formula, data = df_train2[,c("mmr_time", "mmr", top_features)])

# Perform stepwise forward selection using AIC
model.null = coxph(Surv(mmr_time, mmr) ~ elts_class, data = df_train2[,c("mmr_time", "mmr", top_features)])
cox_model = MASS::stepAIC(model.null, direction = "forward",
                          scope = list(lower = model.null,
                                       upper = full_model))
# steps = 10)

# Summary of the Cox model
cox_summary <- summary(cox_model)


# Extract the coefficients from the model
coefficients <- cox_summary$coefficients[, "coef"]

# View the coefficients
print(coefficients)
tmp = data.frame(cox_summary$coefficients)
tmp$covariates = rownames(cox_summary$coefficients)
writexl::write_xlsx(tmp, "./coefficients.xlsx")
writexl::write_xlsx(data.frame(cox_summary$conf.int), "./coefficients_confint.xlsx")
writexl::write_xlsx(data.frame(cox_summary$sctest), "./coefficients_logranktest.xlsx")


# Compute the risk score for each individual
# Multiply each individual's predictor values by the corresponding coefficients
df_train3 = df_train2
df_train3$score <- as.matrix(df_train3[, gsub("TRUE", "", rownames(cox_summary$coefficients))]) %*% coefficients



# VAL
df_val1 <- as.data.frame(lapply(df_val, replace_na_with_median))
df_val1 = df_val1 %>%
  dplyr::mutate(mmr = ifelse(mmr == FALSE, 0, 1))
df_val1$score <- as.matrix(df_val1[, gsub("TRUE", "", rownames(cox_summary$coefficients))]) %*% coefficients


# CUTOFF
cp <- cutpointr(df_val1, score, mmr, pos_class = 1, neg_class = 0, 
                method = maximize_metric, metric = youden)  #prod_sens_spec, prod_ppv_npv, risk_ratio, youden, F1_score, odds_ratio
predicted_labels_val <- ifelse(df_val1$score >= cp$optimal_cutpoint, 1, 0) #cp$optimal_cutpoint


# METRICS
f1_score <- F1_Score(y_pred = predicted_labels_val, y_true = df_val1$mmr, positive = "0")
print(paste("F1-Score:", f1_score))

# Calculate the AUC
df1 = data.frame(V1 = as.factor(df_val1$mmr), predicted_labels = predicted_labels_val)
auc = pROC::roc(df1$V1, df1$predicted_labels, plot=TRUE)
prauc = PRAUC(y_true = df1$V1, y_pred = df1$predicted_labels); prauc


# Generate the confusion matrix
df_val1$pred = as.factor(predicted_labels_val)
conf_matrix <- confusionMatrix(as.factor(predicted_labels_val), as.factor(df_val1$mmr), positive = "1")
conf_matrix1 = as.data.frame(conf_matrix$table) %>%
  group_by(Reference) %>%
  dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
  ungroup() %>%
  dplyr::mutate(Prediction = ifelse(is.na(Prediction), NA,
                                    ifelse(Prediction == 0, "No MMR", "MMR")),
                Prediction = factor(Prediction, levels = c("No MMR", "MMR")),
                Reference = ifelse(is.na(Reference), NA,
                                   ifelse(Reference == 0, "No MMR", "MMR")),
                Reference = factor(Reference, levels = c("No MMR", "MMR")))

# Plot
# Create a nice plot using ggplot2
g = ggplot(data = conf_matrix1, aes(x = reorder(Reference, as.integer(Reference)), y = reorder(Prediction, -as.integer(Prediction)), fill = Prop)) +
  geom_tile(color="black", linewidth = 0.6) +
  geom_text(aes(label=paste0(round(Prop, 1), "%\nn=", round(Freq, 0))), size=4.5) +
  scale_fill_gradient(low="white", high="red") +
  labs(y = "Predicted response", x = "True response", fill = "Proportion (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_text(colour = "black", size = 11),
        legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
ggsave(plot = g, filename = paste0("./confusion_matrix_", seed1, "_val.png"),
       bg = "white", width = 3.5, height = 3.5, units = "in", dpi = 300)


# Combine metrics
TP = conf_matrix$table[2, 2]
FN = conf_matrix$table[1, 2]
TN = conf_matrix$table[1, 1]
FP = conf_matrix$table[2, 1]

# Print the results
metrics_val = data.frame(
  analysis = "val",
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  auc = round(auc$auc, 3),
  precision = round(TP / (TP + FP), 3),
  recall = round(TP / (TP + FN), 3)
) %>%
  dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))

# Save
if (!exists("metrics1")) {
  metrics1 = metrics_val
} else {
  metrics1 = rbind(metrics1, metrics_val)
}


# TEST
df_test1 = df_test %>%
  dplyr::mutate(mmr = ifelse(mmr == FALSE, 0, 1))

df_test1 <- as.data.frame(lapply(df_test1, replace_na_with_median))
df_test1$score <- as.matrix(df_test1[, gsub("TRUE", "", rownames(cox_summary$coefficients))]) %*% coefficients

predicted_labels_test <- ifelse(df_test1$score >= cp$optimal_cutpoint, 1, 0) #cp$optimal_cutpoint

# Calculate the F1-scores
f1_score <- F1_Score(y_pred = predicted_labels_test, y_true = df_test1$mmr, positive = "1")
print(paste("F1-Score:", f1_score))

# Calculate the AUC
df1 = data.frame(V1 = as.factor(df_test1$mmr), predicted_labels = predicted_labels_test)
auc = pROC::roc(df1$V1, df1$predicted_labels, plot=TRUE); auc
prauc = PRAUC(y_true = df1$V1, y_pred = df1$predicted_labels); prauc

# Generate the confusion matrix
df_test1$pred = as.factor(predicted_labels_test)
conf_matrix <- confusionMatrix(as.factor(predicted_labels_test), as.factor(df_test1$mmr), positive = "1")
conf_matrix1 = as.data.frame(conf_matrix$table) %>%
  group_by(Reference) %>%
  dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
  ungroup() %>%
  dplyr::mutate(Prediction = ifelse(is.na(Prediction), NA,
                                    ifelse(Prediction == 0, "No MMR", "MMR")),
                Prediction = factor(Prediction, levels = c("No MMR", "MMR")),
                Reference = ifelse(is.na(Reference), NA,
                                   ifelse(Reference == 0, "No MMR", "MMR")),
                Reference = factor(Reference, levels = c("No MMR", "MMR")))

# Plot
# Create a nice plot using ggplot2
g = ggplot(data = conf_matrix1, aes(x = reorder(Reference, as.integer(Reference)), y = reorder(Prediction, -as.integer(Prediction)), fill = Prop)) +
  geom_tile(color="black", linewidth = 0.6) +
  geom_text(aes(label=paste0(round(Prop, 1), "%\nn=", round(Freq, 0))), size=4.5) +
  scale_fill_gradient(low="white", high="red") +
  labs(y = "Predicted response", x = "True response", fill = "Proportion (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_text(colour = "black", size = 11),
        legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
ggsave(plot = g, filename = paste0("./confusion_matrix_", seed1, "_test.png"),
       bg = "white", width = 3.5, height = 3.5, units = "in", dpi = 300)


# Combine metrics
TP = conf_matrix$table[2, 2]
FN = conf_matrix$table[1, 2]
TN = conf_matrix$table[1, 1]
FP = conf_matrix$table[2, 1]

# Print the results
metrics_test = data.frame(
  analysis = "test",
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  auc = round(auc$auc, 3),
  precision = round(TP / (TP + FP), 3),
  recall = round(TP / (TP + FN), 3)
) %>%
  dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))

# Save
if (!exists("metrics1")) {
  metrics1 = metrics_test
} else {
  metrics1 = rbind(metrics1, metrics_test)
}


# TRAIN
df_train1 <- as.data.frame(lapply(df_train, replace_na_with_median))
df_train1 = df_train1 %>%
  dplyr::mutate(mmr = ifelse(mmr == FALSE, 0, 1))
df_train1$score <- as.matrix(df_train1[, gsub("TRUE", "", rownames(cox_summary$coefficients))]) %*% coefficients

predicted_labels_train <- ifelse(df_train1$score >= cp$optimal_cutpoint, 1, 0)


# METRICS
f1_score <- F1_Score(y_pred = predicted_labels_train, y_true = df_train1$mmr, positive = "0")
print(paste("F1-Score:", f1_score))

# Calculate the AUC
df1 = data.frame(V1 = as.factor(df_train1$mmr), predicted_labels = predicted_labels_train)
auc = pROC::roc(df1$V1, df1$predicted_labels, plot=TRUE)
prauc = PRAUC(y_true = df1$V1, y_pred = df1$predicted_labels); prauc


# Generate the confusion matrix
df_train1$pred = as.factor(predicted_labels_train)
conf_matrix <- confusionMatrix(as.factor(predicted_labels_train), as.factor(df_train1$mmr), positive = "1")
conf_matrix1 = as.data.frame(conf_matrix$table) %>%
  group_by(Reference) %>%
  dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
  ungroup() %>%
  dplyr::mutate(Prediction = ifelse(is.na(Prediction), NA,
                                    ifelse(Prediction == 0, "No MMR", "MMR")),
                Prediction = factor(Prediction, levels = c("No MMR", "MMR")),
                Reference = ifelse(is.na(Reference), NA,
                                   ifelse(Reference == 0, "No MMR", "MMR")),
                Reference = factor(Reference, levels = c("No MMR", "MMR")))

# Plot
# Create a nice plot using ggplot2
g = ggplot(data = conf_matrix1, aes(x = reorder(Reference, as.integer(Reference)), y = reorder(Prediction, -as.integer(Prediction)), fill = Prop)) +
  geom_tile(color="black", linewidth = 0.6) +
  geom_text(aes(label=paste0(round(Prop, 1), "%\nn=", round(Freq, 0))), size=4.5) +
  scale_fill_gradient(low="white", high="red") +
  labs(y = "Predicted response", x = "True response", fill = "Proportion (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_text(colour = "black", size = 11),
        legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
ggsave(plot = g, filename = paste0("./confusion_matrix_", seed1, "_train.png"),
       bg = "white", width = 3.5, height = 3.5, units = "in", dpi = 300)



# Combine metrics
TP = conf_matrix$table[2, 2]
FN = conf_matrix$table[1, 2]
TN = conf_matrix$table[1, 1]
FP = conf_matrix$table[2, 1]

# Print the results
metrics_train = data.frame(
  analysis = "train",
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  auc = round(auc$auc, 3),
  precision = round(TP / (TP + FP), 3),
  recall = round(TP / (TP + FN), 3)
  # seed = seed1
) %>%
  dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))

# Save
if (!exists("metrics1")) {
  metrics1 = metrics_train
} else {
  metrics1 = rbind(metrics1, metrics_train)
}


# Loop by site in the training set
# Generate the confusion matrix
for (center1 in unique(df_train1$center)) {
  
  print(center1)
  
  # Filter
  df_train_center = df_train1 %>%
    dplyr::filter(center == center1) %>%
    dplyr::select(pred, mmr)
  
  # Add dummy row and remove it before plotting as some centers do not have MMR == 0
  df_train_center = rbind(df_train_center, data.frame("pred" = 0, "mmr" = 0))
  
  conf_matrix <- confusionMatrix(as.factor(df_train_center$pred), as.factor(df_train_center$mmr), positive = "1")
  conf_matrix1 = as.data.frame(conf_matrix$table) %>%
    dplyr::mutate(Freq = ifelse(Prediction == 0 & Reference == 0, Freq-1, Freq)) %>%
    group_by(Reference) %>%
    dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
    ungroup() %>%
    dplyr::mutate(Prop = ifelse(is.na(Prop), 0, Prop),
                  Prediction = ifelse(is.na(Prediction), NA,
                                      ifelse(Prediction == 0, "No MMR", "MMR")),
                  Prediction = factor(Prediction, levels = c("No MMR", "MMR")),
                  Reference = ifelse(is.na(Reference), NA,
                                     ifelse(Reference == 0, "No MMR", "MMR")),
                  Reference = factor(Reference, levels = c("No MMR", "MMR")))
  
  # Plot
  # Create a nice plot using ggplot2
  g = ggplot(data = conf_matrix1, aes(x = reorder(Reference, as.integer(Reference)), y = reorder(Prediction, -as.integer(Prediction)), fill = Prop)) +
    geom_tile(color="black", linewidth = 0.6) +
    geom_text(aes(label=paste0(round(Prop, 1), "%\nn=", round(Freq, 0))), size=4.5) +
    scale_fill_gradient(low="white", high="red") +
    labs(y = "Predicted response", x = "True response", fill = "Proportion (%)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 11),
          axis.text.x = element_text(colour = "black", size = 11),
          legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); print(g)
  ggsave(plot = g, filename = paste0("./site/confusion_matrix_", seed1, "_train_", center1, "_center.png"),
         bg = "white", width = 3.5, height = 3.5, units = "in", dpi = 300)
  
}


##### BENCHMARK #####


# Subset data
df_benchmark = df

# ELTS high treated
df_benchmark$elts_class_high = ifelse(is.na(df_benchmark$elts_class), NA, ifelse(df_benchmark$elts_class==2, 0, 1))
benchmark_summary = summary(coxph(Surv(mmr_time, mmr) ~ elts_class_high, data = df_benchmark))

df_benchmark = as.data.frame(lapply(df_benchmark[,c("mmr", "mmr_time", "elts_class_high")], replace_na_with_median))
df_benchmark$mmr = ifelse(df_benchmark$mmr == TRUE, 1, 0)

# METRICS
f1_score <- F1_Score(y_pred = df_benchmark$elts_class_high, y_true = df_benchmark$mmr, positive = "0")
print(paste("F1-Score:", f1_score))

# Calculate the AUC
df1 = data.frame(V1 = as.factor(df_benchmark$mmr), predicted_labels = df_benchmark$elts_class_high)
auc = pROC::roc(df1$V1, df1$predicted_labels, plot=TRUE)
prauc = PRAUC(y_true = df1$V1, y_pred = df1$predicted_labels); prauc

conf_matrix <- confusionMatrix(as.factor(df_benchmark$elts_class_high), as.factor(df_benchmark$mmr), positive = "1")
conf_matrix1 = as.data.frame(conf_matrix$table) %>%
  group_by(Reference) %>%
  dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
  ungroup() %>%
  dplyr::mutate(Prediction = ifelse(is.na(Prediction), NA,
                                    ifelse(Prediction == 0, "No MMR", "MMR")),
                Prediction = factor(Prediction, levels = c("No MMR", "MMR")),
                Reference = ifelse(is.na(Reference), NA,
                                   ifelse(Reference == 0, "No MMR", "MMR")),
                Reference = factor(Reference, levels = c("No MMR", "MMR")))

# Plot
# Create a nice plot using ggplot2
g = ggplot(data = conf_matrix1, aes(x = reorder(Reference, as.integer(Reference)), y = reorder(Prediction, -as.integer(Prediction)), fill = Prop)) +
  geom_tile(color="black", linewidth = 0.6) +
  geom_text(aes(label=paste0(round(Prop, 1), "%\nn=", round(Freq, 0))), size=4.5) +
  scale_fill_gradient(low="white", high="red") +
  labs(y = "Predicted response", x = "True response", fill = "Proportion (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_text(colour = "black", size = 11),
        legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
ggsave(plot = g, filename = paste0("./confusion_matrix_", seed1, "_benchmark.png"),
       bg = "white", width = 3.5, height = 3.5, units = "in", dpi = 300)

# Combine metrics
TP = conf_matrix$table[2, 2]
FN = conf_matrix$table[1, 2]
TN = conf_matrix$table[1, 1]
FP = conf_matrix$table[2, 1]

# Print the results
metrics_benchmark = data.frame(
  analysis = "benchmark",
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  auc = round(auc$auc, 3),
  precision = round(TP / (TP + FP), 3),
  recall = round(TP / (TP + FN), 3)
) %>%
  dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))

# Save
if (!exists("metrics1")) {
  metrics1 = metrics_benchmark
} else {
  metrics1 = rbind(metrics1, metrics_benchmark)
}


##### HALVING TIME #####


# Load data
halving_times_hus = readRDS("./kml_halvingtimes.rds") %>%
  dplyr::select(henkilotunnus, ht=ht1)
halving_times_abroad = readxl::read_xlsx("./halving_time_abroad.xlsx") %>%
  dplyr::select(album_id=`Hemavision ID`, ht=`3-month halving time`)

# Process data
ht_index = 76 # Halving time according to 10.1182/blood-2014-03-566323
df_hus = df %>%
  dplyr::filter(center == "HUS") %>%
  dplyr::select(henkilotunnus, album_id) %>%
  dplyr::left_join(halving_times_hus) %>%
  dplyr::mutate(ht1 = ifelse(ht < ht_index, 1, 0))
df_aus = df %>%
  dplyr::filter(center == "Australia") %>%
  dplyr::select(henkilotunnus, album_id) %>%
  dplyr::left_join(halving_times_abroad) %>%
  dplyr::mutate(ht1 = ifelse(ht < ht_index, 1, 0))
## Join
df_halving_times = rbind(df_hus, df_aus) %>%
  dplyr::left_join(df %>% dplyr::select(henkilotunnus, mmr, mmr_time)) %>%
  distinct()

# Join df and halving time
df = df %>%
  dplyr::left_join(df_halving_times) %>%
  distinct()
chisq.test(df$ht1, df$mmr)
table(df$ht1, df$mmr)

# Calculate the F1-scores
df$mmr = ifelse(df$mmr == TRUE, 1, 0)
f1_score <- F1_Score(y_pred = df$ht1, y_true = df$mmr, positive = "1")
print(paste("F1-Score:", f1_score))

# Calculate the AUC
df1 = df[!is.na(df$ht1)]
df1 = data.frame(V1 = as.factor(df1$mmr), predicted_labels = df1$ht1, ht = df1$ht)
auc = pROC::roc(df1$V1, df1$predicted_labels, plot=TRUE); auc
prauc = PRAUC(y_true = df1$V1, y_pred = df1$predicted_labels); prauc

# Generate the confusion matrix
conf_matrix <- confusionMatrix(as.factor(df1$predicted_labels), as.factor(df1$V1), positive = "1")
conf_matrix1 = as.data.frame(conf_matrix$table) %>%
  group_by(Reference) %>%
  dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
  ungroup() %>%
  dplyr::mutate(Prediction = ifelse(is.na(Prediction), NA,
                                    ifelse(Prediction == 0, "No MMR", "MMR")),
                Prediction = factor(Prediction, levels = c("No MMR", "MMR")),
                Reference = ifelse(is.na(Reference), NA,
                                   ifelse(Reference == 0, "No MMR", "MMR")),
                Reference = factor(Reference, levels = c("No MMR", "MMR")))

# Plot
# Create a nice plot using ggplot2
g = ggplot(data = conf_matrix1, aes(x = reorder(Reference, as.integer(Reference)), y = reorder(Prediction, -as.integer(Prediction)), fill = Prop)) +
  geom_tile(color="black", linewidth = 0.6) +
  geom_text(aes(label=paste0(round(Prop, 1), "%\nn=", round(Freq, 0))), size=4.5) +
  scale_fill_gradient(low="white", high="red") +
  labs(y = "Predicted response", x = "True response", fill = "Proportion (%)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_text(colour = "black", size = 11),
        legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
ggsave(plot = g, filename = paste0("./confusion_matrix_", seed1, "_halving_time.png"),
       bg = "white", width = 3.5, height = 3.5, units = "in", dpi = 300)


# Combine metrics
TP = conf_matrix$table[2, 2]
FN = conf_matrix$table[1, 2]
TN = conf_matrix$table[1, 1]
FP = conf_matrix$table[2, 1]

# Print the results
metrics_ht = data.frame(
  analysis = "halving_time",
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  auc = round(auc$auc, 3),
  precision = round(TP / (TP + FP), 3),
  recall = round(TP / (TP + FN), 3)
) %>%
  dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))

# Save
if (!exists("metrics1")) {
  metrics1 = metrics_ht
} else {
  metrics1 = rbind(metrics1, metrics_ht)
}
writexl::write_xlsx(metrics1, "./model_metrics.xlsx")


# Combine score and halving time
## Val
df_val2 = df_val1 %>% 
  dplyr::select(henkilotunnus, album_id, score, mmr) %>%
  dplyr::mutate(predicted_labels = predicted_labels_val) %>%
  dplyr::inner_join(df_hus)
table(df_val2[df_val2$predicted_labels==0,]$ht1, df_val2[df_val2$predicted_labels==0,]$mmr)

## Test
df_test2 = df_test1 %>% 
  dplyr::select(henkilotunnus, album_id, score, mmr) %>%
  dplyr::mutate(predicted_labels = predicted_labels_test) %>%
  dplyr::left_join(df_aus)
table(df_test2[df_test2$predicted_labels==0,]$ht1, df_test2[df_test2$predicted_labels==0,]$mmr)


##### KAPLAN MEIER ANALYSES #####


# BENCHMARK
# Fit
fit <- survfit(Surv(mmr_time, mmr) ~ elts_class_high, data = df_benchmark)

## Plot
g = ggsurvplot(fit,
               fun = "event",
               data = df_benchmark,
               palette = "Set1",
               size = 2.5,   #line thickness
               ggtheme = theme_minimal(), #theme
               font.main = c(15, "black"), #title font
               font.x = c(15, "bold", "black"), #x-axis font
               font.y = c(15, "bold", "black"), #y-axis font
               font.tickslab = c(15, "bold", "black"), #axis numbering font
               conf.int = FALSE, #confidence interval
               pval = TRUE, #p-value
               pval.size = 5, #p-value size
               pval.coord = c(48, 0.1),
               risk.table.pos = "out",
               tables.y.text = FALSE,
               tables.theme = theme_cleantable(),
               break.x.by = 12,
               risk.table.fontsize = 5,
               risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
               risk.table.col = "strata", #risk table color
               risk.table.height = 0.25,
               ylab = "MMR (%)",
               xlab = "Time (month)",
               xlim = c(0, 60),
               surv.scale = "percent",
               ylim = c(0,1),
               censor = TRUE,
               censor.shape = 108,
               censor.size = 4,
               font.legend = c(15, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
               legend.title = "ELTS high",
               legend.labs = c("TRUE", "FALSE")); g

g1 = plot_grid(g$plot +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size=12, colour = "black"),
                       axis.text.y = element_text(size=12, colour = "black"),
                       axis.title=element_text(size=14, colour = "black"),
                       legend.text = element_text(size=12, colour = "black"),
                       legend.title = element_text(size=12, face="bold", colour = "black"),
                       legend.position = "top"),
               g$table +
                 guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
ggsave(plot = g1, filename = "./kaplan_meier_mmr_benchmark.png", width = 5, height = 5, bg = "white", units = 'in', dpi = 300)





# HALVING TIME
# Fit
df_halving_times1 = df_halving_times %>%
  dplyr::mutate(ht1 = ifelse(is.na(ht1), NA,
                             ifelse(ht1 == 1, "<76 days", "≥76 days")),
                ht1 = factor(ht1, levels = c("<76 days", "≥76 days")))
fit <- survfit(Surv(mmr_time, mmr) ~ ht1, data = df_halving_times1)

## Plot
rhg_cols <- c("#377eb8", "#e41a1c")
g = ggsurvplot(fit,
               fun = "event",
               data = df_halving_times,
               palette = rhg_cols,
               size = 2.5,   #line thickness
               ggtheme = theme_minimal(), #theme
               font.main = c(15, "black"), #title font
               font.x = c(15, "bold", "black"), #x-axis font
               font.y = c(15, "bold", "black"), #y-axis font
               font.tickslab = c(15, "bold", "black"), #axis numbering font
               conf.int = FALSE, #confidence interval
               pval = TRUE, #p-value
               pval.size = 5, #p-value size
               pval.coord = c(48, 0.1),
               risk.table.pos = "out",
               tables.y.text = FALSE,
               tables.theme = theme_cleantable(),
               break.x.by = 12,
               risk.table.fontsize  = 5,
               risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
               risk.table.col = "strata", #risk table color
               risk.table.height = 0.25,
               ylab = "MMR (%)",
               xlab = "Time (month)",
               surv.scale = "percent",
               ylim = c(0,1),
               xlim = c(0, 60),
               censor = TRUE,
               censor.shape = 108,
               censor.size = 4,
               font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
               legend.title = "BCR-ABL1 halving time",
               legend.labs = c("<76 days", "≥76 days")); g

g1 = plot_grid(g$plot +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size=12, colour = "black"),
                       axis.text.y = element_text(size=12, colour = "black"),
                       axis.title=element_text(size=14, colour = "black"),
                       legend.text = element_text(size=12, colour = "black"),
                       legend.title = element_text(size=12, face="bold", colour = "black"),
                       legend.position = "top"),
               g$table +
                 guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
ggsave(plot = g1, filename = "./kaplan_meier_mmr_halving_time.png",
       width = 5, height = 5, bg = "white", units = 'in', dpi = 300)



# SCORE BY TKI

# Fit
df5_ima = rbind(df_val1, df_test1) %>%
  dplyr::filter(imatinib == 1) %>%
  dplyr::mutate(score_cat = ifelse(score > median(score, na.rm=TRUE), "Ima high score", "Ima low score"))
df5_sectki = rbind(df_val1, df_test1) %>%
  dplyr::filter(imatinib == 0) %>%
  dplyr::mutate(score_cat = ifelse(score > median(score, na.rm=TRUE), "2Gen-TKI high score", "2Gen-TKI low score"))
df5 = rbind(df5_ima, df5_sectki)
df5$score_cat = factor(df5$score_cat)
fit <- survfit(Surv(mmr_time, mmr) ~ score_cat, data = df5)

## Plot
g = ggsurvplot(fit,
               fun = "event",
               data = df5,
               palette = "Set1",
               size = 2.5,   #line thickness
               ggtheme = theme_minimal(), #theme
               font.main = c(15, "black"), #title font
               font.x = c(15, "bold", "black"), #x-axis font
               font.y = c(15, "bold", "black"), #y-axis font
               font.tickslab = c(15, "bold", "black"), #axis numbering font
               conf.int = FALSE, #confidence interval
               pval = TRUE, #p-value
               pval.size = 5, #p-value size
               pval.coord = c(48, 0.1),
               risk.table.pos = "out",
               tables.y.text = FALSE,
               tables.theme = theme_cleantable(),
               break.x.by = 12,
               risk.table.fontsize  = 5,
               risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
               risk.table.col = "strata", #risk table color
               risk.table.height = 0.25,
               ylab = "MMR (%)",
               xlab = "Time (month)",
               surv.scale = "percent",
               ylim = c(0,1),
               xlim = c(0, 60),
               censor = TRUE,
               censor.shape = 108,
               censor.size = 4,
               font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
               legend.title = "Morphoclinical score",
               legend.labs = c("2Gen-TKI high", "2Gen-TKI low", "Ima high", "Ima low")); g

g1 = plot_grid(g$plot +
                 guides(colour = guide_legend(nrow = 2)) +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size=12, colour = "black"),
                       axis.text.y = element_text(size=12, colour = "black"),
                       axis.title=element_text(size=14, colour = "black"),
                       legend.text = element_text(size=12, colour = "black"),
                       legend.title = element_text(size=12, face="bold", colour = "black"),
                       legend.position = "top"),
               g$table +
                 guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
ggsave(plot = g1, filename = "./kaplan_meier_mmr_imaging_by_tki_val_test.png",
       width = 5, height = 5, bg = "white", units = 'in', dpi = 300)



# Score by ELTS group

# print(unique(df_test1$elts_class))
# Fit
df6_ELTS_low = rbind(df_val1, df_test1) %>%
  dplyr::filter(elts_class == 0) %>%
  dplyr::mutate(score_cat = ifelse(score > median(score, na.rm=TRUE), "ELTS low-high score", "ELTS low-low score"))
df6_ELTS_intermediate = rbind(df_val1, df_test1) %>%
  dplyr::filter(elts_class == 1) %>%
  dplyr::mutate(score_cat = ifelse(score > median(score, na.rm=TRUE), "ELTS intermed-high score", "ELTS intermed-low score"))
df6_ELTS_high = rbind(df_val1, df_test1) %>%
  dplyr::filter(elts_class == 2) %>%
  dplyr::mutate(score_cat = ifelse(score > median(score, na.rm=TRUE), "ELTS high-high score", "ELTS high-low score"))
df6 = rbind(df6_ELTS_low, df6_ELTS_intermediate, df6_ELTS_high)
df6$score_cat = factor(df6$score_cat)
fit <- survfit(Surv(mmr_time, mmr) ~ score_cat, data = df6)
fit_low <- survfit(Surv(mmr_time, mmr) ~ score_cat, data = df6_ELTS_low)
fit_intermediate <- survfit(Surv(mmr_time, mmr) ~ score_cat, data = df6_ELTS_intermed)
fit_high <- survfit(Surv(mmr_time, mmr) ~ score_cat, data = df6_ELTS_high)

my_colors <- brewer.pal(8, "Set1")   # 9 colors in Set1
my_colors[6] <- "#E69F00" 

## Plots
q_low = ggsurvplot(fit_low,
                   fun = "event",
                   title = "ELTS Low-risk", 
                   data = df6_ELTS_low,
                   palette = my_colors,
                   size = 2,   #line thickness
                   ggtheme = theme_minimal(), #theme
                   font.main = c(11, "black"), #title font
                   font.x = c(11, "black"), #x-axis font
                   font.y = c(11, "black"), #y-axis font
                   font.tickslab = c(11, "black"), #axis numbering font
                   conf.int = FALSE, #confidence interval
                   pval = FALSE, #p-value
                   risk.table.pos = "out",
                   tables.y.text = FALSE,
                   tables.theme = theme_cleantable(),
                   break.x.by = 12,
                   risk.table.fontsize = 4,
                   risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                   risk.table.col = "strata", #risk table color
                   risk.table.height = 0.25,
                   ylab = "MMR (%)",
                   xlab = "Time (month)",
                   surv.scale = "percent",
                   ylim = c(0,1),
                   xlim = c(0, 48),
                   censor = TRUE,
                   censor.shape = 108,
                   censor.size = 4,
                   font.legend = c(11, "black"),    #font voi olla esim. "bold" tai "plain"
                   legend.title = "Morphoclinical score",
                   legend.labs = c("High", "Low"))


q_intermediate = ggsurvplot(fit_intermediate,
                            fun = "event",
                            title = "ELTS Intermediate-risk",
                            data = df6_ELTS_intermediate,
                            palette = my_colors,
                            size = 2,   #line thickness
                            ggtheme = theme_minimal(), #theme
                            font.main = c(11, "black"), #title font
                            font.x = c(11, "black"), #x-axis font
                            font.y = c(11, "black"), #y-axis font
                            font.tickslab = c(11, "black"), #axis numbering font
                            conf.int = FALSE, #confidence interval
                            pval = FALSE, #p-value
                            risk.table.pos = "out",
                            tables.y.text = FALSE,
                            tables.theme = theme_cleantable(),
                            break.x.by = 12,
                            risk.table.fontsize  = 4,
                            risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                            risk.table.col = "strata", #risk table color
                            risk.table.height = 0.25,
                            ylab = "MMR (%)",
                            xlab = "Time (month)",
                            surv.scale = "percent",
                            ylim = c(0,1),
                            xlim = c(0, 48),
                            censor = TRUE,
                            censor.shape = 108,
                            censor.size = 4,
                            font.legend = c(11, "black"),    #font voi olla esim. "bold" tai "plain"
                            legend.title = "Morphoclinical score",
                            legend.labs = c("High", "Low"))
q_high = ggsurvplot(fit_high,
                    fun = "event",
                    title = "ELTS High-risk",
                    data = df6_ELTS_high,
                    palette = my_colors,
                    size = 2,   #line thickness
                    ggtheme = theme_minimal(), #theme
                    font.main = c(11, "black"), #title font
                    font.x = c(11, "black"), #x-axis font
                    font.y = c(11, "black"), #y-axis font
                    font.tickslab = c(11, "black"), #axis numbering font
                    conf.int = FALSE, #confidence interval
                    pval = FALSE, #p-value
                    risk.table.pos = "out",
                    tables.y.text = FALSE,
                    tables.theme = theme_cleantable(),
                    break.x.by = 12,
                    risk.table.fontsize  = 4,
                    risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                    risk.table.col = "strata", #risk table color
                    risk.table.height = 0.25,
                    ylab = "MMR (%)",
                    xlab = "Time (month)",
                    surv.scale = "percent",
                    ylim = c(0,1),
                    xlim = c(0, 48),
                    censor = TRUE,
                    censor.shape = 108,
                    censor.size = 4,
                    font.legend = c(11, "black"),    #font voi olla esim. "bold" tai "plain"
                    legend.title = "Morphoclinical score",
                    legend.labs = c("High", "Low"))



final_plot = (q_low$plot/q_low$table + guides(colour = "none"))|(q_intermediate$plot/q_intermediate$table + guides(colour = "none"))|(q_high$plot/q_high$table + guides(colour = "none"))
final_plot
ggsave(plot = final_plot, filename = paste0(results, "_response/KM_plots_ELTS_risk_groups.png"),
       width = 11, height = 5, bg = "white", units = 'in', dpi = 300)


## Edits ELTS risk groups
# Combine validation and test sets
df6 <- rbind(df_val1, df_test1)

# Create ELTS group factor and model risk groups within each ELTS group
df6 <- df6 %>%
  mutate(ELTS_group = factor(elts_class, labels = c("ELTS Low", "ELTS Intermediate", "ELTS High"))) %>%
  group_by(ELTS_group) %>%
  mutate(model_risk = ifelse(score > median(score, na.rm=TRUE), "High model risk", "Low model risk")) %>%
  ungroup()

# Fit KM curves stratified by model risk and faceted by ELTS
fit <- survfit(Surv(mmr_time, mmr) ~ model_risk, data = df6)

# Plot with facets by ELTS
g <- ggsurvplot_facet(
  fit,
  data = df6,
  facet.by = "ELTS_group",
  palette = c("firebrick", "steelblue"),
  fun = "event",
  size = 2.5,
  ggtheme = theme_minimal(),
  conf.int = FALSE,
  pval = TRUE,
  pval.size = 5,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.col = "strata",
  tables.y.text = FALSE,
  break.x.by = 12,
  ylim = c(0,1),
  xlim = c(0, 60),
  ylab = "MMR (%)",
  xlab = "Time (months)",
  surv.scale = "percent",
  censor = TRUE,
  censor.shape = 108,
  censor.size = 4,
  font.legend = c(12, "bold", "black"),
  legend.title = "Model-predicted risk",
  legend.labs = c("Low", "High")
)

# Clean theme
g$plot <- g$plot +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size=12, colour="black"),
    axis.text.y = element_text(size=12, colour="black"),
    axis.title=element_text(size=14, colour="black"),
    legend.text = element_text(size=12, colour="black"),
    legend.title = element_text(size=12, face="bold", colour="black"),
    legend.position = "top"
  )

# Combine plot + risk table
g1 <- plot_grid(
  g$plot + guides(colour = guide_legend(nrow=1)), 
  g$table + guides(colour = "none"),
  ncol = 1, align = "v", rel_heights = c(2,1)
)
g1
# Save
ggsave(plot = g1, filename = "./kaplan_meier_mmr_by_ELTS_2.png",
       width = 7, height = 5, bg = "white", units = 'in', dpi = 300)


# SCORE TRAIN
# Fit
fit <- survfit(Surv(mmr_time, mmr) ~ pred, data = df_train1)

## Plot
g = ggsurvplot(fit,
               fun = "event",
               data = df_train1,
               palette = "Set1",
               size = 2.5,   #line thickness
               ggtheme = theme_minimal(), #theme
               font.main = c(15, "black"), #title font
               font.x = c(15, "bold", "black"), #x-axis font
               font.y = c(15, "bold", "black"), #y-axis font
               font.tickslab = c(15, "bold", "black"), #axis numbering font
               conf.int = FALSE, #confidence interval
               pval = TRUE, #p-value
               pval.size = 5, #p-value size
               pval.coord = c(48, 0.1),
               risk.table.pos = "out",
               tables.y.text = FALSE,
               tables.theme = theme_cleantable(),
               break.x.by = 12,
               risk.table.fontsize  = 5,
               risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
               risk.table.col = "strata", #risk table color
               risk.table.height = 0.25,
               ylab = "MMR (%)",
               xlab = "Time (month)",
               surv.scale = "percent",
               ylim = c(0,1),
               xlim = c(0, 60),
               censor = TRUE,
               censor.shape = 108,
               censor.size = 4,
               font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
               legend.title = "Morphoclinical score",
               legend.labs = c("Low", "High")); g

g1 = plot_grid(g$plot +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size=12, colour = "black"),
                       axis.text.y = element_text(size=12, colour = "black"),
                       axis.title=element_text(size=14, colour = "black"),
                       legend.text = element_text(size=12, colour = "black"),
                       legend.title = element_text(size=12, face="bold", colour = "black"),
                       legend.position = "top"),
               g$table +
                 guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
ggsave(plot = g1, filename = "./kaplan_meier_mmr_imaging_train.png",
       width = 5, height = 5, bg = "white", units = 'in', dpi = 300)


# SCORE VAL
# Fit
fit <- survfit(Surv(mmr_time, mmr) ~ pred, data = df_val1)

## Plot
g = ggsurvplot(fit,
               fun = "event",
               data = df_val1,
               palette = "Set1",
               size = 2.5,   #line thickness
               ggtheme = theme_minimal(), #theme
               font.main = c(15, "black"), #title font
               font.x = c(15, "bold", "black"), #x-axis font
               font.y = c(15, "bold", "black"), #y-axis font
               font.tickslab = c(15, "bold", "black"), #axis numbering font
               conf.int = FALSE, #confidence interval
               pval = TRUE, #p-value
               pval.size = 5, #p-value size
               pval.coord = c(48, 0.1),
               risk.table.pos = "out",
               tables.y.text = FALSE,
               tables.theme = theme_cleantable(),
               break.x.by = 12,
               risk.table.fontsize  = 5,
               risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
               risk.table.col = "strata", #risk table color
               risk.table.height = 0.25,
               ylab = "MMR (%)",
               xlab = "Time (month)",
               surv.scale = "percent",
               ylim = c(0,1),
               xlim = c(0, 60),
               censor = TRUE,
               censor.shape = 108,
               censor.size = 4,
               font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
               legend.title = "Morphoclinical score",
               legend.labs = c("Low", "High")); g

g1 = plot_grid(g$plot +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size=12, colour = "black"),
                       axis.text.y = element_text(size=12, colour = "black"),
                       axis.title=element_text(size=14, colour = "black"),
                       legend.text = element_text(size=12, colour = "black"),
                       legend.title = element_text(size=12, face="bold", colour = "black"),
                       legend.position = "top"),
               g$table +
                 guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
ggsave(plot = g1, filename = "./kaplan_meier_mmr_imaging_val.png",
       width = 5, height = 5, bg = "white", units = 'in', dpi = 300)


# SCORE TEST
# Fit
fit <- survfit(Surv(mmr_time, mmr) ~ pred, data = df_test1)

## Plot
g = ggsurvplot(fit,
               fun = "event",
               data = df_test1,
               palette = "Set1",
               size = 2.5,   #line thickness
               ggtheme = theme_minimal(), #theme
               font.main = c(15, "black"), #title font
               font.x = c(15, "bold", "black"), #x-axis font
               font.y = c(15, "bold", "black"), #y-axis font
               font.tickslab = c(15, "bold", "black"), #axis numbering font
               conf.int = FALSE, #confidence interval
               pval = TRUE, #p-value
               pval.size = 5, #p-value size
               pval.coord = c(48, 0.1),
               risk.table.pos = "out",
               tables.y.text = FALSE,
               tables.theme = theme_cleantable(),
               break.x.by = 12,
               risk.table.fontsize  = 5,
               risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
               risk.table.col = "strata", #risk table color
               risk.table.height = 0.25,
               ylab = "MMR (%)",
               xlab = "Time (month)",
               surv.scale = "percent",
               ylim = c(0,1),
               xlim = c(0, 60),
               censor = TRUE,
               censor.shape = 108,
               censor.size = 4,
               font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
               legend.title = "Morphoclinical score",
               legend.labs = c("Low", "High")); g

g1 = plot_grid(g$plot +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size=12, colour = "black"),
                       axis.text.y = element_text(size=12, colour = "black"),
                       axis.title=element_text(size=14, colour = "black"),
                       legend.text = element_text(size=12, colour = "black"),
                       legend.title = element_text(size=12, face="bold", colour = "black"),
                       legend.position = "top"),
               g$table +
                 guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
ggsave(plot = g1, filename = "./kaplan_meier_mmr_imaging_test.png",
       width = 5, height = 5, bg = "white", units = 'in', dpi = 300)



##### AUROC and PRROC #####


# First use halving time and later binary halving time

# Continuous halving time

# Calculate and plot the ROC Curve
df_test3 = df_test1 %>%
  dplyr::mutate(response_at_24 = as.numeric(ifelse(response_at_24==TRUE, 1, 0)),
                score = as.numeric(score))
df_benchmark1 = df_benchmark %>%
  dplyr::mutate(response_at_24 = mmr,
                score = as.numeric(elts_class_high))
df1_continuous = df1 %>%
  dplyr::mutate(response_at_24 = V1,
                score = as.numeric(ht))
roc_curve1 <- pROC::roc(df_test3$response_at_24, df_test3$score)
roc_curve2 <- pROC::roc(df_benchmark1$response_at_24, df_benchmark1$score)
roc_curve3 <- pROC::roc(df1_continuous$response_at_24, df1_continuous$score)

# Convert the ROC curve to a data frame for ggplot2
roc_df1 <- data.frame(
  FPR = rev(1 - roc_curve1$specificities),  # False Positive Rate
  TPR = rev(roc_curve1$sensitivities),      # True Positive Rate
  Thresholds = rev(roc_curve1$thresholds),   # Thresholds
  Model = paste0("Morphoclinical (AUROC=", round(roc_curve1$auc, 2), ")")
)
roc_df2 <- data.frame(
  FPR = rev(1 - roc_curve2$specificities),  # False Positive Rate
  TPR = rev(roc_curve2$sensitivities),      # True Positive Rate
  Thresholds = rev(roc_curve2$thresholds),   # Thresholds
  Model = paste0("ELTS high-risk (AUROC=", round(roc_curve2$auc, 2), ")")
)
roc_df3 <- data.frame(
  FPR = rev(1 - roc_curve3$specificities),  # False Positive Rate
  TPR = rev(roc_curve3$sensitivities),      # True Positive Rate
  Thresholds = rev(roc_curve3$thresholds),   # Thresholds
  Model = paste0("BCR-ABL1 halving time (AUROC=", round(roc_curve3$auc, 2), ")")
)

# Join
roc_df = rbind(rbind(roc_df1, roc_df2), roc_df3)

# Order legend
roc_df = roc_df %>%
  dplyr::mutate(Model = factor(Model, levels = c(
    paste0("Morphoclinical (AUROC=", round(roc_curve1$auc, 2), ")"),
    paste0("ELTS high-risk (AUROC=", round(roc_curve2$auc, 2), ")"),
    paste0("BCR-ABL1 halving time (AUROC=", round(roc_curve3$auc, 2), ")")
  )))

# Delong
pROC::roc.test(roc_curve1, roc_curve2, method="delong")
pROC::roc.test(roc_curve1, roc_curve3, method="delong")

# Plot ROC Curve using ggplot2
g = ggplot(roc_df, aes(x = FPR, y = TPR, color=Model)) +
  geom_line(size = 3) +            # Line for ROC curve
  geom_abline(linetype = "dashed") +                  # Diagonal line for random classifier
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Reds", direction = -1) +
  theme(legend.position = c(0.65, 0.15),
        legend.title = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)); g
ggsave(plot = g, filename = "./AUROC_continuous_halving_time.png", bg = "white", width = 5, height = 5, units = "in", dpi = 300)


# Calculate and plot the PR Curve
pr_curve1 <- pr.curve(scores.class0 = df_test3[df_test3$response_at_24 == 1,]$score,
                      scores.class1 = df_test3[df_test3$response_at_24 == 0,]$score,
                      curve = TRUE)
pr_curve2 <- pr.curve(scores.class0 = df_benchmark1[df_benchmark1$response_at_24 == 1,]$score,
                      scores.class1 = df_benchmark1[df_benchmark1$response_at_24 == 0,]$score,
                      curve = TRUE)
pr_curve3 <- pr.curve(scores.class0 = df1_continuous[df1_continuous$response_at_24 == 1,]$score,
                      scores.class1 = df1_continuous[df1_continuous$response_at_24 == 0,]$score,
                      curve = TRUE)
# Extract Precision-Recall data into a data frame
pr_df1 <- data.frame(
  Recall = pr_curve1$curve[, 1],
  Precision = pr_curve1$curve[, 2],
  Model = paste0("Morphoclinical (PRAUC=", round(pr_curve1$auc.integral, 2), ")")
)
pr_df2 <- data.frame(
  Recall = pr_curve2$curve[, 1],
  Precision = pr_curve2$curve[, 2],
  Model = paste0("ELTS high-risk (PRAUC=", round(pr_curve2$auc.integral, 2), ")")
)
pr_df3 <- data.frame(
  Recall = pr_curve3$curve[, 1],
  Precision = pr_curve3$curve[, 2],
  Model = paste0("BCR-ABL1 halving time (PRAUC=", round(pr_curve3$auc.integral, 2), ")")
)

# Join
pr_df = rbind(rbind(pr_df1, pr_df2), pr_df3)

# Order legend
pr_df = pr_df %>%
  dplyr::mutate(Model = factor(Model, levels = c(
    paste0("Morphoclinical (PRAUC=", round(pr_curve1$auc.integral, 2), ")"),
    paste0("ELTS high-risk (PRAUC=", round(pr_curve2$auc.integral, 2), ")"),
    paste0("BCR-ABL1 halving time (PRAUC=", round(pr_curve3$auc.integral, 2), ")")
  )))

# Plot PR Curve using ggplot2
g = ggplot(pr_df, aes(x = Recall, y = Precision, color=Model)) +
  geom_line(size = 3) +
  labs(
    x = "Recall",
    y = "Precision"
  ) +
  scale_color_brewer(palette = "Reds", direction = -1) +
  theme_bw() +
  theme(legend.position = c(0.65, 0.15),
        legend.title = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0.4, 1), expand = c(0.01, 0.01)); g
ggsave(plot = g, filename = "./PRAUC_continuous_halving_time.png", bg = "white", width = 5, height = 5, units = "in", dpi = 300)



# Binary halving time

# Calculate and plot the ROC Curve
df_test3 = df_test1 %>%
  dplyr::mutate(response_at_24 = as.numeric(ifelse(response_at_24==TRUE, 1, 0)),
                score = as.numeric(score))
df_benchmark1 = df_benchmark %>%
  dplyr::mutate(response_at_24 = mmr,
                score = as.numeric(elts_class_high))
df1_binary = df1 %>%
  dplyr::mutate(response_at_24 = V1,
                score = as.numeric(predicted_labels))
roc_curve1 <- pROC::roc(df_test3$response_at_24, df_test3$score)
roc_curve2 <- pROC::roc(df_benchmark1$response_at_24, df_benchmark1$score)
roc_curve3 <- pROC::roc(df1_binary$response_at_24, df1_binary$score)

# Convert the ROC curve to a data frame for ggplot2
roc_df1 <- data.frame(
  FPR = rev(1 - roc_curve1$specificities),  # False Positive Rate
  TPR = rev(roc_curve1$sensitivities),      # True Positive Rate
  Thresholds = rev(roc_curve1$thresholds),   # Thresholds
  Model = paste0("Morphoclinical (AUROC=", round(roc_curve1$auc, 2), ")")
)
roc_df2 <- data.frame(
  FPR = rev(1 - roc_curve2$specificities),  # False Positive Rate
  TPR = rev(roc_curve2$sensitivities),      # True Positive Rate
  Thresholds = rev(roc_curve2$thresholds),   # Thresholds
  Model = paste0("ELTS high-risk (AUROC=", round(roc_curve2$auc, 2), ")")
)
roc_df3 <- data.frame(
  FPR = rev(1 - roc_curve3$specificities),  # False Positive Rate
  TPR = rev(roc_curve3$sensitivities),      # True Positive Rate
  Thresholds = rev(roc_curve3$thresholds),   # Thresholds
  Model = paste0("BCR-ABL1 halving time (AUROC=", round(roc_curve3$auc, 2), ")")
)

# Join
roc_df = rbind(rbind(roc_df1, roc_df2), roc_df3)

# Order legend
roc_df = roc_df %>%
  dplyr::mutate(Model = factor(Model, levels = c(
    paste0("Morphoclinical (AUROC=", round(roc_curve1$auc, 2), ")"),
    paste0("ELTS high-risk (AUROC=", round(roc_curve2$auc, 2), ")"),
    paste0("BCR-ABL1 halving time (AUROC=", round(roc_curve3$auc, 2), ")")
  )))

# Delong
## roc_curve1 = score at test set (n=62)
## roc_curve2 = ELTS at entire set (n=238)
## roc_curve3 = BCR-ABL1 halving time (n=166)
pROC::roc.test(roc_curve1, roc_curve2, method="delong")
pROC::roc.test(roc_curve1, roc_curve3, method="delong")

# Plot ROC Curve using ggplot2
g = ggplot(roc_df, aes(x = FPR, y = TPR, color=Model)) +
  geom_line(size = 3) +            # Line for ROC curve
  geom_abline(linetype = "dashed") +                  # Diagonal line for random classifier
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Reds", direction = -1) +
  theme(legend.position = c(0.65, 0.15),
        legend.title = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)); g
ggsave(plot = g, filename = "./AUROC_binary_halving_time.png", bg = "white", width = 5, height = 5, units = "in", dpi = 300)

# Calculate and plot the PR Curve
pr_curve1 <- pr.curve(scores.class0 = df_test3[df_test3$response_at_24 == 1,]$score,
                      scores.class1 = df_test3[df_test3$response_at_24 == 0,]$score,
                      curve = TRUE)
pr_curve2 <- pr.curve(scores.class0 = df_benchmark1[df_benchmark1$response_at_24 == 1,]$score,
                      scores.class1 = df_benchmark1[df_benchmark1$response_at_24 == 0,]$score,
                      curve = TRUE)
pr_curve3 <- pr.curve(scores.class0 = df1_binary[df1_binary$response_at_24 == 1,]$score,
                      scores.class1 = df1_binary[df1_binary$response_at_24 == 0,]$score,
                      curve = TRUE)
# Extract Precision-Recall data into a data frame
pr_df1 <- data.frame(
  Recall = pr_curve1$curve[, 1],
  Precision = pr_curve1$curve[, 2],
  Model = paste0("Morphoclinical (PRAUC=", round(pr_curve1$auc.integral, 2), ")")
)
pr_df2 <- data.frame(
  Recall = pr_curve2$curve[, 1],
  Precision = pr_curve2$curve[, 2],
  Model = paste0("ELTS high-risk (PRAUC=", round(pr_curve2$auc.integral, 2), ")")
)
pr_df3 <- data.frame(
  Recall = pr_curve3$curve[, 1],
  Precision = pr_curve3$curve[, 2],
  Model = paste0("BCR-ABL1 halving time (PRAUC=", round(pr_curve3$auc.integral, 2), ")")
)

# Join
pr_df = rbind(rbind(pr_df1, pr_df2), pr_df3)

# Order legend
pr_df = pr_df %>%
  dplyr::mutate(Model = factor(Model, levels = c(
    paste0("Morphoclinical (PRAUC=", round(pr_curve1$auc.integral, 2), ")"),
    paste0("ELTS high-risk (PRAUC=", round(pr_curve2$auc.integral, 2), ")"),
    paste0("BCR-ABL1 halving time (PRAUC=", round(pr_curve3$auc.integral, 2), ")")
  )))

# Plot PR Curve using ggplot2
g = ggplot(pr_df, aes(x = Recall, y = Precision, color=Model)) +
  geom_line(size = 3) +            # Line for PR curve
  labs(
    x = "Recall",
    y = "Precision"
  ) +
  scale_color_brewer(palette = "Reds", direction = -1) +
  theme_bw() +
  theme(legend.position = c(0.65, 0.15),
        legend.title = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0.01, 0.01)); g
ggsave(plot = g, filename = "./PRAUC_binary_halving_time.png", bg = "white", width = 5, height = 5, units = "in", dpi = 300)

