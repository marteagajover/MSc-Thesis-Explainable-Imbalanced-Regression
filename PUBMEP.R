#################################################################################
# PUBMEP PRE-PROCESS #
#################################################################################

# packages
packages <- c("dplyr", "openxlsx")
for (package in packages) {
  if (!require(package, character.only = TRUE)) 
    {install.packages(package)}}

# import data
DB_path <- "C:/Users/Maria/Desktop/TFM/BDs/BASE_PUBMEP_LONGITUDINAL_LONGformat_NINAS_Y_NINOS_213inds_31_01_2022_EPIC.csv"
DB_sep = ";"
DB_dec = ","
DB_header <- TRUE


# create new variables
AL <- c(3354, 3361)
BMI_father <- c(3048, 3049)
BMI_mother <- c(3050, 3051)

# obmetrics
obmetrics_input_directory <- "C:/Users/Maria/Desktop/Bioinformática/4 semestre/Datasets/OBMETRICS/PUBMEP_obmetrics.xlsx"
obmetrics_output_directory <- "C:/Users/Maria/Desktop/Bioinformática/4 semestre/Datasets/OBMETRICS/PUBMEP_obmetrics_output.xlsx"

# DEXA T1
GENOBOX_path <- "C:/Users/Maria/Desktop/TFM/BDs/Genobox_integrada_2023_02_01.csv"

# DEXA T2
DEXA_T2_path <- "C:/Users/Maria/Desktop/TFM/BDs/DEXA/originales/2024_03_21_tanita_dexa_pubmep (2).csv"
DEXA_T2_sep = ","
DEXA_T2_dec = "."
DEXA_T2_header <- TRUE

# Vit D T2
VitD_T2_path <- "C:/Users/Maria/Desktop/TFM/BDs/Vit D/Codigos_Base_PUBMEP_Santiago_VitD 03 02 20.xlsx"


# filter subjects by 1 variable (e.g., acc)
filter <- 0
# acc_vars <- 3363:3383
# max_var <- acc_vars[which.max(sapply(DB_total[acc_vars], function(var) sum(!is.na(var))))]
# filter <- max_var

# filer variables
vars <- c(
  # Metadata
  1, 3, 3130, 3010, 2563,
  
  # Background_information
  3013:3015, 3017, 3021:3022, 3024:3029, 3033:3034, 3040:3041, 3533:3534,
  3052:3053, 3055, 3060:3061, 3070,
  
  # Physical_examination
  # 3012, 3105, 3117, 3119, 3121,
  3545, 3105, 3117, 3119, 3121,
  
  # Pre-examination_questions
  3126, 3128:3129,
  
  # Biochemical_analysis
  3307:3311, 3315:3321, 3324:3328, 3526, 3329, 3544, 3331, 3333,
  3322, 3334, 3360, 3358, 3357, 3354, 3361, 3532,
  
  # Dual_energy_X_ray_absorptiometry
  3541:3543, 
  
  # Accelerometer
  # 3363:3383,
  
  # Outcomes
  3349,
  
  # z-scores
  3535:3540,
  
  # Cole
  3350)

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
  "zDBP_NHBPEP", "zSBP_NHBPEP", "Cole"
)

# factorize label-encoded-variables
factor_var <- c(
  # metadata
  "Origin", "Analysis_month", "Sex", 
  
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

# correct variables 
NA_value <- 98 # some NAs are encoded as 98

# import PGS scores
scores_path <- "C:/Users/Maria/Desktop/TFM/BDs/PGS/PUBMEP_PGS/scores/scores.txt"
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
  DB_total <- read.table(file = DB_path, sep = DB_sep, dec = DB_dec, header = DB_header, stringsAsFactors = TRUE)
} else {
  warning("DB path not provided.")
}


## ----------------------------------------------------------------------------------------------------------------------------------------
# CREATE NEW VARIABLES
# these variables are created in the total DB so they can be selected later, in the correct order

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
## OBMETRICS
# these variables are created in the total DB so they can be selected later, in the correct order

# save a CSV with the data necessary to calculate the z-scores via Obmetrics
DB_Obmetrics <- DB_total[c("Code", "Age_Time", "Sex.1_T1", "Talla_m_Time",
                           "Peso_Kg_Time", "WC_Time", "DBP_Time", "SBP_Time",
                           "TAG__mg_dl__Time", "HDLc__mg_dl__Time",
                           "Glucose__mg_dl__Time", "Insulin__mU_l__Time")]
DB_Obmetrics$Tanner_Time <- c(DB_total$Tanner2_T1[1:213], DB_total$Pubarquia_T2[1:213])
DB_Obmetrics$Code <- seq(1, 426)


# change names to the proper ones
names(DB_Obmetrics) <- c("id", "decimal_age", "sex", "height_m", "weight_kg", 
                         "wc_cm", "dbp_mmHg", "sbp_mmHg", "tg_mg_dl", "hdl_mg_dl",
                         "glucose_mg_dl", "insulin_microU_ml", "tanner_index")

# save results
# write.xlsx(DB_Obmetrics, file = obmetrics_input_directory,
#            sep = ";", dec = ",", rowNames = FALSE)

# USE OBMETRICS

# load the calculated z-scores
obmetrics <- read.xlsx(xlsxFile = obmetrics_output_directory)
DB_total["zWC_Ferrandez"] <- obmetrics["Obesity"]
DB_total["zTAG_Stavnsbo"] <- obmetrics["Tryglicerides"]
DB_total["zHDL_Stavnsbo"] <- obmetrics["HDL"]
DB_total["zHOMA_IR_Stavnsbo"] <- obmetrics["Insulin_resistance"]
DB_total["zDPB_NHBPEP"] <- obmetrics["DBP"]
DB_total["zSBP_NHBPEP"] <- obmetrics["SBP"]

## ----------------------------------------------------------------------------------------------------------------------------------------
# DEXA T1
# DB with T1 values (GENOBOX)
if (GENOBOX_path != ""){
  GENOBOX <- read.table(file = GENOBOX_path, sep = DB_sep, dec = DB_dec, 
                        header = DB_header, stringsAsFactors = TRUE)
} else {
  warning("DEXA T1 DB path not provided.")
}

# add "Mmagratg" and "Mgrasatg"
DB_total <- DB_total %>% mutate(row_id = row_number()) # row ID
DB_merged <- merge(DB_total, GENOBOX[, c("Code", "Mmagratg", "Mgrasatg")], 
                   by = "Code", all.x = TRUE, sort = FALSE) # merge
DB_merged <- unique(DB_merged) # delete duplicated
DB_merged <- DB_merged[order(DB_merged$row_id), ] # order as DB_total
DB_total$row_id <- NULL # delete row id
DB_merged$Mmagratg[214:426] <- NA # delete T1 values wrongly added to T2
DB_merged$Mgrasatg[214:426] <- NA # delete T1 values wrongly added to T2
DB_total[c("Mmagratg", "Mgrasatg")] <- DB_merged[c("Mmagratg", "Mgrasatg")] # add columns
rm(DB_merged) # delete DB_merged


## ----------------------------------------------------------------------------------------------------------------------------------------
# DEXA T2
if (DEXA_T2_path != ""){
  DEXA_T2 <- read.table(file = DEXA_T2_path, sep = DEXA_T2_sep, dec = DEXA_T2_dec,
                        header = DEXA_T2_header, stringsAsFactors = TRUE)
} else {
  warning("DEXA T2 DB path not provided.")
}

names(DEXA_T2) <- lapply(names(DEXA_T2), function(x) {
  paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep = "")
})

# outliers corrections (correct mistakes)

# DEXA

# Mmagratg → DXA_Base_Datos_MMagra(kg)
DEXA_T2[DEXA_T2$Code_new_t2 == "MS252",]$Mmagratg <- 47148
DEXA_T2[DEXA_T2$Code_new_t2 == "MS292",]$Mmagratg <- 10954

# Mgrasatg → DXA_Base_Datos_MGrasa(kg)
DEXA_T2[DEXA_T2$Code_new_t2 == "MS252",]$Mgrasatg <- 14574
DEXA_T2[DEXA_T2$Code_new_t2 == "MS304",]$Mgrasatg <- 16432

# porcgrasat → Masa_Grasa_%_DXA

# tanita

# Tan_LM_percent → Masa_magra_%_tan

# Tan_FM_percent → Masa_grasa_%_tan
DEXA_T2[DEXA_T2$Code_new_t2 == "MS228",]$MG_percent_Tan <- 71080
DEXA_T2[DEXA_T2$Code_new_t2 == "MS205",]$MG_percent_Tan <- 55500
DEXA_T2[DEXA_T2$Code_new_t2 == "MS241",]$MG_percent_Tan <- 73970
DEXA_T2[DEXA_T2$Code_new_t2 == "MS244",]$MG_percent_Tan <- 82740

# ---

# save T1 values
Mmagratg_T1 <- DB_total$Mmagratg[1:213]
Mgrasatg_T1 <- DB_total$Mgrasatg[1:213]
Porcgrasat_T1 <- DB_total$percent_Grasa_corporal_T1[1:213]

# delete variables
DB_total$Mmagratg<- NULL
DB_total$Mgrasatg <- NULL

# add T2 values
DB_total <- DB_total %>% mutate(row_id = row_number()) # row ID
DB_merged <- merge(DB_total, DEXA_T2[, c("Code", "Mmagratg", "Mgrasatg", "Porcgrasat")], 
                   by = "Code", all.x = TRUE, sort = FALSE) # merge
DB_merged <- unique(DB_merged) # delete duplicated
DB_merged <- DB_merged[order(DB_merged$row_id), ] # order as DB_total
DB_total$row_id <- NULL # delete row id
DB_merged$Mmagratg[1:213] <- NA # delete T2 values wrongly added to T1
DB_merged$Mgrasatg[1:213] <- NA # delete T2 values wrongly added to T1
DB_total[c("Mmagratg", "Mgrasatg", "Porcgrasat")] <- DB_merged[c("Mmagratg", 
                                                                 "Mgrasatg", "Porcgrasat")] # add columns
rm(DB_merged) # delete DB_merged

# restore T1 values
DB_total$Mmagratg[1:213] <- Mmagratg_T1 
DB_total$Mgrasatg[1:213] <- Mgrasatg_T1 
DB_total$Porcgrasat[1:213] <- Porcgrasat_T1 

## ----------------------------------------------------------------------------------------------------------------------------------------
# Vitamin D T2
vitD_T2 <- read.xlsx(VitD_T2_path)

# save T1 value
vitD_T1 <- DB_total$Vit.D_T1[1:213]
names(vitD_T2) <- c("Code_new_T2", "Code", "VitD")

# add T2 values
DB_total <- DB_total %>% mutate(row_id = row_number()) # row ID
DB_merged <- merge(DB_total, vitD_T2[, c("Code_new_T2", "VitD")], 
                   by = "Code_new_T2", all.x = TRUE, sort = FALSE) # merge
DB_merged <- unique(DB_merged) # delete duplicated
DB_merged <- DB_merged[order(DB_merged$row_id), ] # order as DB_total
DB_total$row_id <- NULL # delete row id
DB_merged$VitD[1:213] <- NA # delete T2 values wrongly added to T1
DB_total[c("VitD")] <- DB_merged[c("VitD")] # add column
rm(DB_merged) # delete DB_merged

# restore T1 values
DB_total$VitD[1:213] <- vitD_T1 

## ----------------------------------------------------------------------------------------------------------------------------------------
# estadio tanner
DB_total$Tanner_stage <- c(DB_total$Tanner2_T1[1:213], DB_total$Pubarquia_T2[1:213])

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

# DB$Analysis_month <- as.factor(format(as.Date(as.character(DB$Analysis_month)), "%m"))
DB$Analysis_month <- as.factor(format(as.Date(as.character(DB$Analysis_month), 
                                              format = "%m/%d/%y"), "%m"))

# recode marital status (simplification)

DB$Parents_marital_status <- factor(
  ifelse(DB$Parents_marital_status %in% c(0, 4), 0,
         ifelse(DB$Parents_marital_status %in% c(1, 2), 1,
                ifelse(DB$Parents_marital_status == 3, 2, NA))),
  levels = 0:2,
  labels = c("0", "1", "2")
)

## ----------------------------------------------------------------------------------------------------------------------------------------
# FACTORIZE LABEL-ENCODED-VARIABLES
if (length(factor_var) > 0){
  DB[factor_var] <- lapply(DB[factor_var], as.factor)
}

DB$Code <- as.factor(DB$Code)


## ----------------------------------------------------------------------------------------------------------------------------------------
# DIVIDE DB IN T1 AND T2
# T1
DB_T1 <- DB[1:(nrow(DB)/2),]
DB_T1 <- DB_T1[DB_T1$Cole %in% c("Obeso_T1", "Sobrepeso_T1"), ]
DB_T1$Cole <- NULL # delete Cole (neccesary to mselect obese/overweight)

# T2
DB_T2 <- DB[(nrow(DB)/2+1):nrow(DB),]
DB_T2 <- DB_T2[DB_T2$Cole %in% c("Obeso_T2", "Sobrepeso_T2"), ]
DB_T2$Cole <- NULL # delete Cole (neccesary to mselect obese/overweight)

## ----------------------------------------------------------------------------------------------------------------------------------------
# IMPORT GPS SCORES

# subjects with genomics data
scores <- read.table(scores_path, sep=",", header = TRUE)

if (length(scores_var) > 0){
  scores <- scores[scores_var]
}

# merge previous data and PGS
pgs <- function(database){
  database <- database %>% mutate(row_id = row_number()) # row ID
  database <- merge(database, scores, by.x = "Code", by.y = "sample", all = FALSE, sort = FALSE)
  database <- database[order(database$row_id), ] # order as DB_total
  database$Code <- NULL   # delete previous data ID (neccesary to merge)
  database$sample <- NULL # delete PGS data ID (neccesary to merge)
  database$row_id <- NULL # delete PGS data ID (neccesary to merge)
  return(database)
}

DB_T1 <- pgs(DB_T1)
DB_T2 <- pgs(DB_T2)


## ----------------------------------------------------------------------------------------------------------------------------------------
# exclude Córdoba
DB_T1 <- DB_T1[DB_T1$Origin != 0,]
DB_T2 <- DB_T2[DB_T2$Origin != 0,] 


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
DB_list_T1 <- objective_data(data = DB_T1, vars = targets, del_vars = del_var)
DB_list_T2 <- objective_data(data = DB_T2, vars = targets, del_vars = del_var)

for (df in names(DB_list_T1[[1]])){
  missing_percentage <- sapply(DB_list_T1[[df]], function(col) sum(is.na(col))
                               / length(col) * 100)
  DB_list_T1_1_missing <- missing_percentage[missing_percentage > pmissing]
  DB_list_T1_1_non_selected_columns <- names(missing_percentage[missing_percentage > pmissing])
  DB_list_T1_1_selected_columns <- names(missing_percentage[missing_percentage <= pmissing])
}
aa <- intersect(names(DB_list[[1]]), DB_list_T1_1_non_selected_columns)

for (df in names(DB_list_T1[[2]])){
  missing_percentage <- sapply(DB_list_T1[[df]], function(col) sum(is.na(col))
                               / length(col) * 100)
  DB_list_T1_2_missing <- missing_percentage[missing_percentage > pmissing]
  DB_list_T1_2_non_selected_columns <- names(missing_percentage[missing_percentage > pmissing])
  DB_list_T1_2_selected_columns <- names(missing_percentage[missing_percentage <= pmissing])
}
intersect(names(DB_list[[2]]), DB_list_T1_2_non_selected_columns)

for (df in names(DB_list_T2[[1]])){
  missing_percentage <- sapply(DB_list_T2[[df]], function(col) sum(is.na(col))
                               / length(col) * 100)
  DB_list_T2_1_missing <- missing_percentage[missing_percentage > pmissing]
  DB_list_T2_1_non_selected_columns <- names(missing_percentage[missing_percentage > pmissing])
  DB_list_T2_1_selected_columns <- names(missing_percentage[missing_percentage <= pmissing])
}
intersect(names(DB_list[[1]]), DB_list_T2_1_non_selected_columns)


for (df in names(DB_list_T2[[2]])){
  missing_percentage <- sapply(DB_list_T2[[df]], function(col) sum(is.na(col))
                               / length(col) * 100)
  DB_list_T2_2_missing <- missing_percentage[missing_percentage > pmissing]
  DB_list_T2_2_non_selected_columns <- names(missing_percentage[missing_percentage > pmissing])
  DB_list_T2_2_selected_columns <- names(missing_percentage[missing_percentage <= pmissing])
}
intersect(names(DB_list[[2]]), DB_list_T2_2_non_selected_columns)


DB_list_T1_nm <- DB_list_T1[[1]][complete.cases(DB_list_T1[[1]][, aa]), ]
DB_all_missing <- DB_list_T1[rowSums(is.na(DB_list_T1[, aa])) == length(aa), ]


for (df in names(DB_list_T1_nm[[2]])){
  missing_percentage <- sapply(DB_list_T2[[df]], function(col) sum(is.na(col))
                               / length(col) * 100)
  m <- missing_percentage[missing_percentage > pmissing]
  nm <- names(missing_percentage[missing_percentage > pmissing])
  nnm <- names(missing_percentage[missing_percentage <= pmissing])
}

DB_list_T1_nm <- DB_list[[1]][complete.cases(DB_list[[1]][, names(DB_list[[1]])]), ]


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
  DB_list_T1 <- var_sel(data_list = DB_list_T1, missing = pmissing)
  DB_list_T2 <- var_sel(data_list = DB_list_T2, missing = pmissing)
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
  DB_list_T1 <- imp(DB_list_T1)
  DB_list_T2 <- imp(DB_list_T2)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# SAVE RESULTS

# function to save results as different csv
csv <- function(data_list, directory, end = ""){
  for (df in names(data_list)){
    write.csv(x = data_list[[df]], 
              paste0(file = directory, "/", df, end, ".csv"),
              row.names = FALSE)
  }
}

# save results as csv
# csv(data_list = DB_list_T1, directory = results_directory, end = "PM_T1.csv")
# csv(data_list = DB_list_T1, directory = results_directory, end = "PM_T2.csv")
