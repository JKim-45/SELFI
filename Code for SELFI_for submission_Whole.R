

# This software is licensed under the MIT License.
# See the license.txt file for more information.

# Package loading
library(readxl)
library(pROC)
library(openxlsx)
library(dplyr)
library(caret)    # kind of dummy package

# set random seeds
set.seed(123)

# Raw data reading
tdata <- read_excel("Whole.xlsx", sheet = "Sheet1", col_names = TRUE)

# Data separation
roc_obj_C <- tdata$'SELFI'   # C means CLFA (= SELFI)
roc_obj_E <- tdata$'ELISA'  # E means ELISA
roc_obj_A <- tdata$'AuNPs'  # A means single AuNPs
roc_obj_stat <- tdata$'stats' # stat means presence of PDAC (0 : healthy control, 1 : PDAC patient)

# number of repeatation
n_of_exp <- 100

# size of each sample group
n_of_samples_H <- 50    # H means Healthy control
n_of_samples_E <- 40    # E means Early stage patient
n_of_samples_L <- 60    # L means Late stage patient

# ratio of exclude
r_of_exc <- 0.1

# variable for average calculation
auc_CLFA <- c()
auc_ELISA <- c()
auc_AuNPs <- c()

# variable for ROC curve stacking
roc_data_C <- list()
roc_data_E <- list()
roc_data_A <- list()

# variable for average ROC of LFIA
max_sen <- c()
min_sen <- c()

# get the current folder address
current_dir <- getwd()

# get the YYMMDD_hhmm (Y - year, M - month, D - day, h - hour, m - minute)
folder_name <- format(Sys.time(), "%y%m%d_%H%M_Whole samples")

# generate new folder with YYMMDD_hhmm
new_dir <- file.path(current_dir, folder_name)
if (!dir.exists(new_dir)) {
  dir.create(new_dir) }

# change the path to new dir
setwd(new_dir)


# Starting the repeatation
for(i in 1:n_of_exp)
{
  
  # Random generation of exclude number of each sample
  exclude_indices1 <- sample(1:n_of_samples_H,n_of_samples_H * r_of_exc)
  exclude_indices2 <- sample(1:n_of_samples_E,n_of_samples_E * r_of_exc)
  exclude_indices3 <- sample(1:n_of_samples_L,n_of_samples_L * r_of_exc)
  exclude_indices <- c(exclude_indices1, exclude_indices2 + n_of_samples_H, exclude_indices3 + n_of_samples_H + n_of_samples_E)
  
  # Generation of new vector variable, after excluding from original data
  roc_CLFA <- roc_obj_C[-exclude_indices]
  roc_ELISA <- roc_obj_E[-exclude_indices]
  roc_AuNPs <- roc_obj_A[-exclude_indices]
  roc_stat <- roc_obj_stat[-exclude_indices]

  # ROC curve set generation
  roc_curve_C <- roc(roc_stat, roc_CLFA)
  roc_curve_E <- roc(roc_stat, roc_ELISA)
  roc_curve_A <- roc(roc_stat, roc_AuNPs)

  # Making xlsx file with ROC curve sets
  roc_C_temp <- data.frame(
    threshold = roc_curve_C$thresholds,
    sensitivity = roc_curve_C$sensitivities,
    specificity = roc_curve_C$specificities
  )
  file_name <- paste("CLFA_ROC_",as.character(i))
  file_name <- paste(file_name,".xlsx")
  write.xlsx(roc_C_temp, file = file_name, sheetName = "temp", rownames = FALSE)
  
  roc_E_temp <- data.frame(
    threshold = roc_curve_E$thresholds,
    sensitivity = roc_curve_E$sensitivities,
    specificity = roc_curve_E$specificities
  )
  file_name <- paste("ELISA_ROC_",as.character(i))
  file_name <- paste(file_name,".xlsx")
  write.xlsx(roc_E_temp, file = file_name, sheetName = "temp", rownames = FALSE)
  
  roc_A_temp <- data.frame(
    threshold = roc_curve_A$thresholds,
    sensitivity = roc_curve_A$sensitivities,
    specificity = roc_curve_A$specificities
  )
  file_name <- paste("AuNPs_ROC_",as.character(i))
  file_name <- paste(file_name,".xlsx")
  write.xlsx(roc_A_temp, file = file_name, sheetName = "temp", rownames = FALSE)
  
  # stacking the ROC curve results
  roc_data_C <- append(roc_data_C, list(roc_curve_C))
  roc_data_E <- append(roc_data_E, list(roc_curve_E))
  roc_data_A <- append(roc_data_A, list(roc_curve_A))
  
  # stacking the AUC value
  auc_CLFA <- c(auc_CLFA, auc(roc_curve_C))
  auc_ELISA <- c(auc_ELISA, auc(roc_curve_E))
  auc_AuNPs <- c(auc_AuNPs, auc(roc_curve_A))
  
  # finding the maximum sensitivity at 1-specificity = 0 in LFIA
  max_sen <- c(max_sen, roc_A_temp %>%
    filter(specificity == 1) %>%
    summarise(max_sensitivity = max(sensitivity)) %>%
    pull(max_sensitivity))
  
  # finding the minimum sensitivity at 1-specificity = 1 in LFIA
  min_sen <- c(min_sen, roc_A_temp %>%
    filter(specificity == 0) %>%
    summarise(min_sensitivity = min(sensitivity)) %>%
    pull(min_sensitivity))
  
}

# print summary of AUC (for confirmation)
#summary (auc_CLFA)
#summary (auc_ELISA)

# Drawing the average ROC curve with linear interpolation
all_fpr_C <- sort(unique(unlist(lapply(roc_data_C, function(roc) roc$specificities))))
mean_tpr_C <- sapply(all_fpr_C, function(fpr) {
  tprs <- sapply(roc_data_C, function(roc) {
    approx(roc$specificities, roc$sensitivities, xout = fpr, rule = 2)$y
  })
  mean(tprs)
})

all_fpr_E <- sort(unique(unlist(lapply(roc_data_E, function(roc) roc$specificities))))
mean_tpr_E <- sapply(all_fpr_E, function(fpr) {
  tprs <- sapply(roc_data_E, function(roc) {
    approx(roc$specificities, roc$sensitivities, xout = fpr, rule = 2)$y
  })
  mean(tprs)
})

## Because ROC curves of LFIA were not suitable for linear interpolation, we didn't use linear interpolation for LFIA
##all_fpr_A <- sort(unique(unlist(lapply(roc_data_A, function(roc) roc$specificities))))
##mean_tpr_A <- sapply(all_fpr_A, function(fpr) {
##  tprs <- sapply(roc_data_A, function(roc) {
##    approx(roc$specificities, roc$sensitivities, xout = fpr, rule = 2)$y
##  })
##  mean(tprs)
##})

mean_fpr_C <- all_fpr_C
mean_fpr_E <- all_fpr_E
##mean_fpr_A <- all_fpr_A

# Dataframe for making xlsx files
roc_avg_C <- data.frame(
  criterion = mean_fpr_C,
  FPR = 1 - mean_fpr_C,   
  TPR = mean_tpr_C
)

roc_avg_E <- data.frame(
  criterion = mean_fpr_E,
  FPR = 1 - mean_fpr_E,   
  TPR = mean_tpr_E
)

roc_avg_A <- data.frame(
  FPR = c(1, 1, 0, 0),
  TPR = c(1, mean(min_sen), mean(max_sen), 0)
)

# AUC dataframe for xlsx files
C_data_new <- data.frame(group = "CLFA", value = auc_CLFA)
E_data_new <- data.frame(group = "ELISA", value = auc_ELISA)
A_data_new <- data.frame(group = "AuNPs", value = auc_AuNPs)
df_auc <- cbind(ELISA = auc_ELISA, SELFI = auc_CLFA, AuNPs = auc_AuNPs)

# Average AUC and 95% CI
avg_AUC <- colMeans(df_auc)
lower_ci <- c()
upper_ci <- c()
for (i in 1:3) {
  ci_value <- quantile(df_auc[,i], probs = c(0.025, 0.975))
  lower_ci[i] <- ci_value[1]
  upper_ci[i] <- ci_value[2]
}
roc_ci_df <- data.frame ("Average AUC value" = avg_AUC, "Lower CI" = lower_ci, "Upper CI" = upper_ci)
rownames(roc_ci_df) <- c("ELISA", "SELFI", "AuNPs")

# xlsx file generation
write.xlsx(roc_avg_C, file = "Average_ROC_data_SELFI.xlsx", sheetName = "Average_ROC_SELFI", rowNames = FALSE)
write.xlsx(roc_avg_E, file = "Average_ROC_data_ELISA.xlsx", sheetName = "Average_ROC_ELISA", rowNames = FALSE)
write.xlsx(roc_avg_A, file = "Average_ROC_data_AuNPs.xlsx", sheetName = "Average_ROC_AuNPs", rowNames = FALSE)
write.xlsx(df_auc, file = "AUC values.xlsx", sheetName = "AUC_values", rowNames = FALSE)
write.xlsx(roc_ci_df, file = "auc_average.xlsx", sheetName = "AUC_CI", rowNames = TRUE)

# change the path to original folder
setwd(current_dir)