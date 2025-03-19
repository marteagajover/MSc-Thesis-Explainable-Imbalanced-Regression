# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 17:07:07 2024

@author: Maria

Instructions

The working directory must include two folders:
- Datasets: It MUST contain one folder for each dataset (with all the partitions)
- Results: It will contain the results as *.txt
- Graphs: It will contain the graphs as *.txt

Algorithm: Random Forest Regressor

Description: A random forest is a meta estimator that fits a number of 
classifying decision trees on various sub-samples of the dataset and uses 
averaging to improve the predictive accuracy and control over-fitting. 

Available: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html
"""

# Previous:
    # Create 3 folders: Datasets (with the datasets), "Results" and "Graphs"
    # script "calculo_f1.R" in this same directory
    # rpy2 (in the cmd) pip install rpy2==3.5.1
    # R_HOME (in Windows, in "Editing System Environment Variables", define R_HOME as C:/Program Files/R/R-4.2.2/) 
        # alternative (not working) (change path) os.environ["R_HOME"] = r"C:/Program Files/R/R-4.2.2/"
    # RStudio (in RStudio)
        # install.packages("remotes")
        # remotes::install_github("cran/DMwR")
        # remotes::install_github("rpribeiro/uba")
        # (download uba, github: https://github.com/paobranco/Pre-processingApproachesImbalanceRegression?tab=readme-ov-file)
        # install.packages("uba_0.7.7.tar.gz",repos=NULL, dependencies=T)
        # The datasets folders and archives names bust be the same
        # (e.g., folder "zHOMA", partition "zHOMA-1tra.dat")

# import
from sklearn.metrics import r2_score
import os
import numpy as np
from time import time
from sklearn.ensemble import RandomForestRegressor
import rpy2.robjects as robjects
import pandas as pd
import shap
import matplotlib.pyplot as plt


# previous
plt.close('all')
cwd = "C:/Users/Maria/Desktop/TFM/Exp"
f1file = "/calculo_f1.R"
kfold = 5
times = 6

# directory
algorithm = "RandomForest"
os.chdir(cwd)
directory = cwd + "/Datasets"
directoryR = cwd + f1file
directorygraphs = cwd + "/Graphs"
dataset = os.listdir(directory)
seeds = [12345678]

# Loading datasets

# zHOMA_acc
# dataset = ["zHOMA_acc"] 
# features_names = [["Origin", "Analysis_month", "Age", "Sex", "Gestational_age", 
#                     "Birth_weight", "Tanner_stage", "%FM_tan", "%LM_tan", "%BW_tan", 
#                     "Uric_acid", "Urea", "Creatinine", "Protein", "LDLc", "Bilirubin", 
#                     "AST", "ALT", "GGT", "ALP", "Calcium", "Sodium", "Potassium", 
#                     "LH", "FSH", "TSH", "Cortisol", "Testosterone", "Estradiol", 
#                     "PTH", "Vitamin_D", "Erythrocytes", "Haemoglobin", "MCV", "MCH", 
#                     "MCHC", "Leukocytes", "%Neutrophils", "%Lymphocytes", "%Monocytes", 
#                     "%Eosinophils", "%Basophils", "Platelets", "MPV", "Iron", 
#                     "Transferrin", "Ferritin", "Adiponectin", "Leptin", "ALR", 
#                     "Number_days", "SedentaryWd", "Sedentary_weekend", "Sedentary", 
#                     "LightWEd", "LightWEd", "Light", "ModerateWd", "ModerateWEd", 
#                     "Moderate", "VigorousWd", "VigorousWEd", "Vigorous", "cpmWd", 
#                     "cpmWEd", "cpm", "StepsWd", "StepsWEd", "Steps", "mvpaWd", 
#                     "mvpaWd", "mvpa", "zBMI", "zWC", "zTG", 
#                     "zHDL", "zDBP", "zSBP", "zHOMA"]]

# zHOMA_PGS
# dataset = ["zHOMA_pgs"]
# features_names = [["Origin", "Analysis_month", "Age", "Sex", "Civil_status_par", 
#                     "Obesity_par", "Diabetes_par", "MI_par", "Cerebrovascular_par", 
#                     "Hyperchol_par", "Obesity_mat_gp", "Diabetes_mat_gp", 
#                     "Gestational_age", "Birth_weight", "Tanner_stage", "%FM_tan", 
#                     "%LM_tan", "%BW_tan", "Uric_acid", "Urea", "Creatinine", "Protein", 
#                     "LDLc", "AST", "ALT", "GGT", "ALP", "LH", "FSH", "TSH", "Testosterone", 
#                     "Estradiol", "Vitamin_D", "Haemoglobin", "Leukocytes", "Iron", 
#                     "%FM_DXA", "%FM_DXA", "zBMI", "zWC", "zTG", "zHDL", "zDBP", "zSBP", 
#                     "PGS000306", "PGS000308", "PGS000299", "PGS000843", "PGS000027", 
#                     "PGS002313", "PGS000716", "PGS000305", "PGS000307", "PGS001350", 
#                     "PGS000839", "PGS003469", "PGS003470", "PGS000877", "PGS000834", 
#                     "PGS000871", "PGS000837", "PGS003400", "PGS000840", "PGS003124", 
#                     "PGS000021", "PGS000014", "PGS000020", "PGS000330", "PGS000729", 
#                     "PGS002243", "PGS000844", "PGS000828", "PGS000827", "PGS000842", 
#                     "zHOMA"]]

# zWC_acc
# dataset = ["zWC_acc"]
# features_names = [["Origin", "Analysis_month", "Age", "Sex", "Gestational_age", 
#                     "Birth_weight", "Tanner_stage", "%Tan_FM", "%Tan_LM", "%BW_tan", 
#                     "Glucose", "Uric_acid", "Urea", "Creatinine", "Protein", "LDLc", 
#                     "AST", "ALT", "GGT", "ALP", "Calcium", "Sodium", "Potassium", 
#                     "LH", "FSH", "TSH", "Cortisol", "Testosterone", "Estradiol", 
#                     "PTH", "Vitamin_D", "Erythrocytes", "Haemoglobin", "MCV", "HCM", 
#                     "CHCM", "%Neutrophils", "%Lymphocytes", "%Monocytes", "%Eosinophils", 
#                     "%Basophils", "Platelets", "VPM", "Iron", "Transferrin", "Ferritin", 
#                     "Adiponectin", "Leptin", "ALR", "number_days", "SedentaryWd", 
#                     "SedentaryWEd", "Sedentary", "LightWd", "LightWEd", "Light", 
#                     "ModerateWd", "ModerateWEd", "Moderate", "VigorousWd", "VigorousWEd", 
#                     "Vigorous", "cpmWd", "cpmWEd", "cpm", "StepsWd", "StepsWEd", 
#                     "Steps", "mvpaWd", "mvpaWEd", "mvpa", "zBMI", 
#                     "zTG", "zHDL", "zHOMA", "zDBP", 
#                     "zSBP", "zWC"]]

# zWC_pgs
dataset = ["zWC_pgs"] 
features_names = [["Origin", "Analysis_month", "Age", "Sex", "Diabetes_par", 
                    "MI_par", "Cerebrovascular_par", "Hyperchol_par", "Obesity_mat_gp", 
                    "Diabetes_mat_gp", "Gestational_age", "Birth_weight", 
                    "Tanner_stage", "%FM_tan", "%LM_tan", "%BW_tan", "Glucose", 
                    "Uric_acid", "Urea", "Creatinine", "Protein", "LDLc", "AST", 
                    "ALT", "GGT", "ALP", "LH", "FSH", "TSH", "Testosterone", "Estradiol", 
                    "Vitamin_D", "Haemoglobin", "Leukocytes", "Iron", "IL8", "%LM_tan", 
                    "%FM_DXA", "%FM_DXA", "zBMI", "zTG", "zHDL", "zHOMA", "zDBP", "zSBP", 
                    "PGS000306", "PGS000308", "PGS000299", "PGS000843", "PGS000027", 
                    "PGS002313", "PGS000716", "PGS000305", "PGS000307", "PGS001350", 
                    "PGS000839", "PGS003469", "PGS003470", "PGS000877", "PGS000834", 
                    "PGS000871", "PGS000837", "PGS003400", "PGS000840", "PGS003124", 
                    "PGS000021", "PGS000014", "PGS000020", "PGS000330", "PGS000729", 
                    "PGS002243", "PGS000844", "PGS000828", "PGS000827", "PGS000842", 
                    "zWC"]]

# results
header=("Algorithm").rjust(len(algorithm)+1) + ("Dataset").rjust(12) + ("ExeTime").rjust(12) + ("MSE-tra").rjust(12) + ("SD-tra").rjust(12) + ("MSE-tst").rjust(12) + ("SD-tst").rjust(12)+ ("F1-tst").rjust(12)+ ("MAE-tra").rjust(12)+ ("MAE-tst").rjust(12)+ ("R2-tst").rjust(12)+("MAPE-tst").rjust(12)
rfile = open(cwd + "/Results/" + algorithm + "_Results.txt", 'w+')
rfile.write(header)
rfile.write("\n")
print(header)


for k in range(len(dataset)):
    rtime = []
    resulttraX = []
    resulttestX = []
    resulttraXmae = []
    resulttestXmae = []
    resulttestXmape = []
    resultr2 = []
    resultf1 = []

    shap_array = pd.DataFrame()
    for i in range(1, kfold +1):
        for t in range(1, times +1): 
            
            # training file path
            
            # normal:
            # train_file = directory + "/" + dataset[k] + "/" + dataset[k] + str(t)+ "-" + str(i) + ".tra"
            train_file = directory + "/" + dataset[k] + "/" + dataset[k] + "-" + str(t) + "-" + str(i) + "tra.dat"

            # SMOGN:
            # train_file = directory + "/" + dataset[k] + "/" + dataset[k] + str(t) + "X-" + str(i) + ".tra"
            
            # test file path
            # test_file = directory + "/" + dataset[k] + "/" + dataset[k] + str(t) + "-" + str(i) + ".tst"
            test_file = directory + "/" + dataset[k] + "/" + dataset[k] + "-" + str(t) + "-" + str(i) + "tst.dat"

            # training
            # datasettra = pd.read_csv(train_file, delimiter="\t", skiprows=2, comment='@', header=None)
            # not  delimiter=", " (that space causes an error)
            datasettra = pd.read_csv(train_file, delimiter=",", comment='@', header=None)
            # attributes = datasettra.shape[1]-1
            x = datasettra.iloc[:, :-1]
            y = datasettra.iloc[:, -1]
                        
            # test
            # datasettest = pd.read_csv(test_file, delimiter="\t", skiprows=2, comment='@', header=None)
            # not  delimiter=", " (that space causes an error)
            datasettest = pd.read_csv(test_file, delimiter=",", comment='@', header=None)
            xtest = datasettest.iloc[:, :-1]
            ytest = datasettest.iloc[:, -1]
            
            # features names
            x.columns = features_names[k][0:-1]
            y.columns = features_names[k][-1]
            xtest.columns = features_names[k][0:-1]
            ytest.columns = features_names[k][-1] 
            
            # execution
            for seed in seeds:
                initial_time = time() 
     
    			# model
                regr = RandomForestRegressor(n_estimators = 500, random_state = seed)
    
    			# training the model using the training sets
                regr.fit(x, y)
                 

                # SHAP
                
                # explainer
                # explainer = shap.Explainer(regr, x, feature_names=x.columns, check_additivity=False)
                explainer = shap.Explainer(regr, x, feature_names=x.columns)
                
                # SHAP values
                shap_values = explainer.shap_values(x, check_additivity=False)
                
                shap_val = pd.DataFrame(shap_values, columns=features_names[k][0:-1])
                shap_abs = shap_val.abs().sum()
                shap_array[f"shap_abs_{t}_{i}_{seed}"] = shap_abs
                  
                # summary plot
                # Calculates the average importance of each feature across all samples in the data set.
                plt.figure(figsize=(6, 10))
                shap.summary_plot(shap_values, features=x, feature_names=x.columns)
                plt.title("Summary plot (" + dataset[k] + ")")
                plt.savefig(directorygraphs + '/summary_plot_' + dataset[k] + '-' + str(t) + '-' + str(i) + '.png')
                                            
                training_output = np.mean((regr.predict(x) - y) ** 2)
                test_output = np.mean((regr.predict(xtest) - ytest) ** 2)
                
                training_outputmae = np.mean(abs(regr.predict(x) - y))
                test_outputmae = np.mean(abs(regr.predict(xtest) - ytest))
                
               
                # f1
                np.savetxt('y.txt', ytest)
                np.savetxt('y_pred.txt', regr.predict(xtest))
                               
                r_source = robjects.r['source']
                r_source(directoryR)
                
                f1=np.loadtxt('f1.txt')
                final_time = time() 
                execution_time = final_time - initial_time
     
                resultf1.append(f1)                          
                rtime.append(execution_time)
                resulttestX.append(test_output)
                resulttraX.append(training_output)
                
                #mae_value = mean_absolute_error(ytest, regr.predict(xtest))# alternative to calculate mae
                #resulttestXmae.append(mae_value)
                resulttestXmae.append(test_outputmae)
                resulttraXmae.append(training_outputmae)


                # Calcular MAPE
                mape = np.mean(np.abs((np.array(ytest) - np.array(regr.predict(xtest))) / np.maximum(np.array(ytest), 1))) * 100
                resulttestXmape.append(mape)
                
                # Calcular R^2
                r_squared = r2_score(ytest, regr.predict(xtest))
                if r_squared < 0:
                    r_squared = 0
                resultr2.append(r_squared)

    exestime = np.mean(rtime)
    finalF1 = np.mean(resultf1)
    
    ecmtra = np.mean(resulttraX)
    sdtra = np.std(resulttraX)
    ecmtst = np.mean(resulttestX)
    sdtst = np.std(resulttestX)

    ecmtramae= np.mean(resulttraXmae)
    ecmtstmae = np.mean(resulttestXmae)
    r2 = np.mean(resultr2)
    ecmtstmape = np.mean(resulttestXmape)


    # results
    results = (algorithm).rjust(len(algorithm)+1) + (dataset[k]).rjust(12) + str(round(exestime, 6)).rjust(11) +" "+ str(round(ecmtra, 6)).rjust(11) +" " + str(round(sdtra, 6)).rjust(11) +" "  + str(round(ecmtst, 6)).rjust(11) +" " + str(round(sdtst, 6)).rjust(11) +" " + str(round(finalF1, 6)).rjust(12) +" "  + str(round(ecmtramae, 6)).rjust(11)+" "  + str(round(ecmtstmae, 6)).rjust(11)+" "  + str(round(r2, 6)).rjust(11)+" "  + str(round(ecmtstmape, 6)).rjust(11)+" " 
    print(results)
    rfile = open(cwd + "/Results/" + algorithm + "_Results.txt", 'w+')
    rfile.write(results)
    rfile.write("\n")
    
    # shap_array_abs = shap_array.abs().sum()
    shap_array_abs = shap_array.sum(axis=1) # axis1 = sum by rows
    shap_array_sort = shap_array_abs.sort_values(ascending=False) # order (descending)

    del resulttraX[:]
    del resulttestX[:]
    del rtime[:]
    
# SHAP values to xlsx
shap_array_sort.reset_index().to_excel('SHAP_values.xlsx')