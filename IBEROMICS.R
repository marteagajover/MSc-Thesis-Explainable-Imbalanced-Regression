#################################################################################
# IBEROMICS PRE-PROCESS #
#################################################################################

# packages
packages <- c("dplyr", "openxlsx", "lubridate")
for (package in packages) {
  if (!require(package, character.only = TRUE)) 
  {install.packages(package)}}

# import data
DB_path <- "C:/Users/Maria/Desktop/TFM/BDs/2024_02_14_Base_integrada_OMICS.csv"
DB_sep = ";"
DB_dec = ","
DB_header <- TRUE

# create new variables
AL <- c()
BMI_father <- c(45, 46)
BMI_mother <- c(47, 48)

# obmetrics
obmetrics_input_directory <- "C:/Users/Maria/Desktop/Bioinformática/4 semestre/Datasets/OBMETRICS/ZIBEROMICS_obmetrics.xlsx"
obmetrics_output_BMI_directory <- "C:/Users/Maria/Desktop/Bioinformática/4 semestre/Datasets/OBMETRICS/ZIBEROMICS_obmetrics_output_BMI.xlsx"
obmetrics_output_WC_directory <- "C:/Users/Maria/Desktop/Bioinformática/4 semestre/Datasets/OBMETRICS/ZIBEROMICS_obmetrics_output_WC.xlsx"

# filter subjects by 1 variable (e.g., acc)
# filter <- 563
filter <- 0

# filer variables
vars <- c(1, 5, 7, 8, 9, # Metadata
         10:12, 14, 18, 19, 21:26, 30, 31, 37, 38, 673:674, 49:51, 57, 58, 67, # Background_information
         104, 108, 113, 115, 117, # Physical_examination
         119, 121, 122, # Pre-examination_questions
         # 123:127, 130:132, 134:137, 141:147, 150, 152, 156, 164, 167, 173, 176:181, # Biochemical_analysis
         123:127, 130:132, 134:137, 141:147, 150, 152, 156, 164, 167, 176:181, # Biochemical_analysis
         431:433, # Dual_energy_X_ray_absorptiometry
         # 565:585, # Accelerometer
         602, 607, 628, 629, 631, 633, 634) # Outcomes

# correct OZ136 BMI (wrong age)
real_zbmi <- 5.06462

# new variables names
var_names <- c(
  "Code", "Origin", "Analysis_month", "Age", "Sex", "Father_education",
  "Mother_education", "Father_working", "Mother_working", "Parents_marital_status",
  "Who_spends_more_time_with_child", "Obesity_parents", "Diabetes_parents",
  "Myocardial_infarction_parents", "Vascular_cerebral_problems_parents", "Hypercholesterolemia_parents",
  "Elevated_TG_parents", "Obesity_maternal_grandparents", "Diabetes_maternal_grandparents",
  "Obesity_paternal_grandparents", "Diabetes_paternal_grandparents", "BMI_father", "BMI_mother",
  "Pregnancy_weight_gain", "Gestational_diabetes", "DM_treatment", "Gestational_age",
  "Birth_weight", "Exclusive_Breastfeeding_Duration", "Tanner_stage", "Acantosis_nigricans",
  "Tan_FM_percent", "Tan_LM_percent", "Tan_ACT_percent", "Illnes_AINEs_last_week",
  "Fasting_hours", "Child_intense_exercise_last_24_hours",
  "Glucose", "Uric_acid", "Urea", "Creatinine", "Protein", "LDLc", "APO_A", "APO_B",
  "AST", "ALT", "GGT", "ALP", "LH", "FSH", "TSH", "T4", "Cortisol", "Testosterone",
  "Estradiol", "Vitamin_D", "Haemoglobin", "Leukocytes", "Iron",
  "CRP", "TNF", "IL8", "IL6", "Adiponectin", "Leptin", "ALR",
  "Lean_mass", "Fat_mass", "Fat_mass_percentage", "zBMI_Orbegozo",
  "zWC_Ferrandez", "zTG_Stavnsbo", "zHDL_Stavnsbo", "zHOMA_Stavnsbo",
  "zDBP_NHBPEP", "zSBP_NHBPEP"
)

# factorize label-encoded-variables
factor_var <- c(
  # metadata
  "Code", "Origin", "Analysis_month", "Sex", 
  
  # Background_information
  "Father_education", "Mother_education", "Father_working", "Mother_working", 
  "Parents_marital_status", "Who_spends_more_time_with_child", "Obesity_parents", 
  "Diabetes_parents", "Myocardial_infarction_parents", "Vascular_cerebral_problems_parents", 
  "Hypercholesterolemia_parents", "Elevated_TG_parents", 
  "Obesity_maternal_grandparents", "Diabetes_maternal_grandparents", 
  "Obesity_paternal_grandparents", "Diabetes_paternal_grandparents", 
  "Gestational_diabetes", "DM_treatment", 
  
  # Physical_examination
  "Tanner_stage", "Acantosis_nigricans",
  
  # Pre-examination_questions
  "Illnes_AINEs_last_week", 
  "Child_intense_exercise_last_24_hours"
)

# import PGS scores
scores_path <- "C:/Users/Maria/Desktop/TFM/BDs/PGS/IBEROMICS_PGS/scores/scores.txt"
scores_var <- c(
  "sample",    # ID
  "PGS000306", # BMI-adjusted fasting blood glucose measurement
  "PGS000308", #	BMI-adjusted fasting blood insulin measurement
  "PGS000299", #	BMI-adjusted waist-hip ratio
  "PGS000843", #	BMI-adjusted waist-hip ratio
  "PGS000027", #	body mass index
  "PGS002313", #	body mass index
  "PGS000716", #	Early life body size, body mass index, comparative body size at age 10, self-reported
  "PGS000305", #	fasting blood glucose measurement
  "PGS000307", #	fasting blood insulin measurement
  "PGS001350", #	fasting blood glucose measurement
  "PGS000839", #	glucose tolerance test
  "PGS003469", #	HOMA-B
  "PGS003470", #	HOMA-IR
  "PGS000877", #	insulin resistance
  "PGS000834", #	insulin response measurement
  "PGS000871", #	insulin secretion measurement
  "PGS000837", #	insulin sensitivity measurement
  "PGS003400", #	Obesity
  "PGS000840", #	proinsulin measurement
  "PGS003124", #	Thinness, body mass index
  "PGS000021", #	type 1 diabetes mellitus
  "PGS000014", #	type 2 diabetes mellitus
  "PGS000020", #	type 2 diabetes mellitus
  "PGS000330", #	type 2 diabetes mellitus
  "PGS000729", #	type 2 diabetes mellitus
  "PGS002243", #	type 2 diabetes mellitus
  "PGS000844", #	visceral adipose tissue measurement
  "PGS000828", #	Waist circumference (female)
  "PGS000827", #	Waist circumference (male)
  "PGS000842" #	waist-hip ratio
) # PGS

# objective selection
targets <- c("zHOMA_Stavnsbo", "zWC_Ferrandez") # target variables
# del_vars: variables to delete for a certain variable (e.g., delete glucose for zHOMA-IR)
del_var <- list(zHOMA_Stavnsbo=c("Glucose"))

# select variables < X% missing data
pmissing <- 15

# impute with MissForest
imputation <- TRUE

# save results
results_directory <- "C:/Users/Maria/Desktop/Bioinformática/4 semestre/Datasets/BDs/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# IMPORT DATA

# import data in CSV
if (DB_path != ""){
  DB_total <- read.table(file = DB_path, sep = DB_sep, dec = DB_dec, 
                         header = DB_header, stringsAsFactors = TRUE)
} else {
  warning("DB path not provided.")
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# CREATE NEW VARIABLES

# calculate the adipoectin-leptin ratio
if (length(AL) == 2){
  DB_total["ALR"] <- DB_total[AL[1]] / DB_total[AL[2]]
}

# calculate the parents BMI
if (length(BMI_father) == 2){
  DB_total["BMI_father"] <- DB_total[BMI_father[1]] / (DB_total[BMI_father[2]]/100)^2}

if (length(BMI_mother) == 2){
  DB_total["BMI_mother"] <- DB_total[BMI_mother[1]] / (DB_total[BMI_mother[2]]/100)^2
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# transform age to REAL decimal age in Origen=2 (not real decimal age for Zaragoza)
DB_total[DB_total$Origen == 2,]$Age <- round(time_length(interval(dmy(DB_total[DB_total$Origen == 2,]$Fecha_de_nacimiento),
                                                                  dmy(DB_total[DB_total$Origen == 2,]$Fecha_analitica)), "years"), 4)

## ----------------------------------------------------------------------------------------------------------------------------------------
## OBMETRICS
# these variables are created in the total DB so they can be selected later, in the correct order

# save a CSV with the data necessary to calculate the z-scores via Obmetrics
DB_Obmetrics <- DB_total[DB_total$Origen == 2, c("Code_omics", "Age", "Sex", "Height",
                           "Peso_Kg", "Perimetro_de_cintura", "TD", "TS",
                           "TAG__mg_dl_", "HDLc__mg_dl_","Glucose__mg_dl_", 
                           "Insulin__mU_ml_", "Estadio_tanner")]
DB_Obmetrics$Height <- DB_Obmetrics$Height / 100
DB_Obmetrics$Sex[DB_Obmetrics$Sex == 1] <- 0
DB_Obmetrics$Sex[DB_Obmetrics$Sex == 2] <- 1
DB_Obmetrics$Code_omics <- seq(1, nrow(DB_Obmetrics))

# change names to the proper ones
names(DB_Obmetrics) <- c("id", "decimal_age", "sex", "height_m", "weight_kg", 
                         "wc_cm", "dbp_mmHg", "sbp_mmHg", "tg_mg_dl", "hdl_mg_dl",
                         "glucose_mg_dl", "insulin_microU_ml", "tanner_index")

# save results
# write.xlsx(DB_Obmetrics, file = obmetrics_input_directory,
#            sep = ";", dec = ",", rowNames = FALSE)

# USE OBMETRICS

# load the calculated z-scores
obmetrics <- read.xlsx(xlsxFile = obmetrics_output_BMI_directory)
DB_total[DB_total$Origen ==2, ]["BMI_zscore_Orbegozo_longi"] <- obmetrics["Obesity"]
DB_total[DB_total$Origen ==2, ]["Tryglicerides_zscore_Stavnsbo"] <- obmetrics["Tryglicerides"]
DB_total[DB_total$Origen ==2, ]["HDL_zscore_Stavnsbo"] <- obmetrics["HDL"]
DB_total[DB_total$Origen ==2, ]["HOMA_zscore_Stavnsbo"] <- obmetrics["Insulin_resistance"]
DB_total[DB_total$Origen ==2, ]["DBP_zscore_Stavnsbo"] <- obmetrics["DBP"]
DB_total[DB_total$Origen ==2, ]["SBP_zscore_Stavnsbo"] <- obmetrics["SBP"]
obmetrics_WC <- read.xlsx(xlsxFile = obmetrics_output_WC_directory)
DB_total[DB_total$Origen ==2, ]["WC_zscore_Ferrandez"] <- obmetrics_WC["Obesity"]

## ----------------------------------------------------------------------------------------------------------------------------------------
# FILTER SUBJECTS BY 1 NON-MISSING VARIABLE (e.g., ACC VARIABLES)

# DB filtered by 1 variable (e.g., those subjects with data in GWAS, EWAS, acc, etc.)
if (filter > 0){
  DB <- DB_total[!is.na(DB_total[filter]), ]
} else {
  DB <- DB_total
}


## ----------------------------------------------------------------------------------------------------------------------------------------
# FILTER VARIABLES

# DB filtered (variables to maintain)
if (length(vars) > 0){
  DB <- DB[vars]
}


## ----------------------------------------------------------------------------------------------------------------------------------------
# RENAME VARIABLES

# new variables names
names(DB) <- var_names

## ----------------------------------------------------------------------------------------------------------------------------------------
# TRANSFORM VARIABLES TO A MORE SUITABLE FORMAT

# transform a date to month
DB$Analysis_month <- as.factor(format(as.Date(as.character(DB$Analysis_month)), "%m"))

# recode marital status (simplification)
DB$Parents_marital_status <- factor(
  ifelse(DB$Parents_marital_status %in% c(0, 4), 0,
         ifelse(DB$Parents_marital_status %in% c(1, 2), 1,
                ifelse(DB$Parents_marital_status == 3, 2, NA))),
  levels = 0:2,
  labels = c("0", "1", "2")
)

# recode sex (from 1: boy, 1: girl to 0: boy, 1: girl)
DB$Sex <- ifelse(DB$Sex == 1, 0, 1)

# transform urea from g/L to mg/dL
DB$Urea <- DB$Urea * 100

# transform urea from thousands/mm3 to units/mm3
DB$Leukocytes <- DB$Leukocytes * 1000

## ----------------------------------------------------------------------------------------------------------------------------------------
# FACTORIZE LABEL-ENCODED-VARIABLES 

if (length(factor_var) > 0){
  DB[factor_var] <- lapply(DB[factor_var], as.factor)
}


## ----------------------------------------------------------------------------------------------------------------------------------------
# IMPORT GPS SCORES

# subjects with genomics data
scores <- read.table(scores_path, sep=",", header = TRUE)

if (length(scores_var) > 0){
  scores <- scores[scores_var]
}

# merge previous data and PGS
DB <- merge(DB, scores, by.x = "Code", by.y = "sample", all = FALSE)
DB$Code <- NULL   # delete previous data ID (neccesary to merge)
DB$sample <- NULL # delete PGS data ID (neccesary to merge)

# missing_percentages <- sapply(DB, function(x) sum(is.na(x)) / length(x) * 100)
# variables_with_missing <- names(missing_percentages[missing_percentages > 15])
# variables_with_missing


## ----------------------------------------------------------------------------------------------------------------------------------------
# OBJECTIVE SELECTION

# function to create DBs as different list based on the target variable
objective_data <- function(data, vars, del_vars){
  data_list <- list()
  
  objective <- function(df, obj_var, data){  
    data_obj_var <- df[obj_var]
    df <- df[, !names(df) %in% obj_var]
    df <- cbind(df, data_obj_var)
    df <- df[!is.na(df[obj_var]),]
    return(df)
  }
  
  for (var in vars){
    df <- objective(df = data, obj_var = var, data = data)  
    for (del_var in del_vars[[var]]){
      df[, del_var] <- NULL
    }
    data_list[[var]] <- df
  }
  
  return(data_list)
}

# list of DBs (1 db for target variable)
DB_list <- objective_data(data = DB, vars = targets, del_vars = del_var)

# for (df in names(DB_list[[1]])){
#   missing_percentage <- sapply(DB_list[[df]], function(col) sum(is.na(col)) 
#                                / length(col) * 100)
#   missing_percentage[missing_percentage >= pmissing]
#   selected_columns <- names(missing_percentage[missing_percentage <= pmissing])
#   DB_vars <- DB_list[[df]][selected_columns]
# }
# 
# for (df in names(DB_list[[2]])){
#   missing_percentage <- sapply(DB_list[[df]], function(col) sum(is.na(col))
#                                / length(col) * 100)
#   missing_percentage[missing_percentage >= pmissing]
#   selected_columns <- names(missing_percentage[missing_percentage <= pmissing])
#   DB_vars <- DB_list[[df]][selected_columns]
# }

## ----------------------------------------------------------------------------------------------------------------------------------------
# SELECT VARIABLES < X% MISSING DATA
var_sel <- function(data_list, missing){
  for (df in names(data_list)){
    missing_percentage <- sapply(data_list[[df]], function(col) sum(is.na(col)) 
                                 / length(col) * 100)
    selected_columns <- names(missing_percentage[missing_percentage <= missing])
    data_list[[df]] <- data_list[[df]][selected_columns]
  }
  return(data_list)
}

# DB without missing
if (pmissing > 0){
  DB_list <- var_sel(data_list = DB_list, missing = pmissing)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# IMPUTE WITH MISSFOREST

# function to impute data
imp <- function(data_list){
  if (!require("missForest")){
    install.packages("missForest")}
  library(missForest)
  for (df in names(data_list)){
    data_list[[df]] <- missForest(data_list[[df]])$ximp
  }
  return(data_list)
}

# impute data
if (imputation == TRUE){
  DB_list <- imp(DB_list)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# SAVE RESULTS

csv <- function(data_list, directory, end=""){
  for (df in names(data_list)){
    write.csv(x = data_list[[df]], 
              paste0(file = directory, "/", df, end, ".csv"),
              row.names = FALSE)
  }
}

# save results as csv
#csv(data_list = DB_list, directory = results_directory)
