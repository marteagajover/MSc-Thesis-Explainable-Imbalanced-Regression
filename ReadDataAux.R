# # Function to read datasets in KEEL format (*.dat).
# 
# readData = function(name, location = ""){
#   
#   r = read.csv(paste(".", location,"/",name, sep = ""), comment.char = "@", header = FALSE)
#   colnames(r)[ncol(r)] = "y"
#   r
#   # add name class:
# }
# 
# 
# # Cross Validation with stratification (generating the folds)
# cvStrat = function (DS, folds = 10) {   
#   permutation = sample(1:nrow(DS))
#   y = DS[permutation,]$y
#   names(y) = permutation
#   if (is.numeric(y)) {
#     cuts <- floor(length(y)/folds)
#     if (cuts < 2) 
#       cuts <- 2
#     if (cuts > 5) 
#       cuts <- 5
#     breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
#     y <- cut(y, breaks, include.lowest = TRUE)
#   }
#   if (folds < length(y)) {
#     y <- factor(as.character(y))
#     numInClass <- table(y)
#     foldVector <- vector(mode = "integer", length(y))
#     for (i in 1:length(numInClass)) {
#       min_reps <- numInClass[i]%/%folds
#       if (min_reps > 0) {
#         spares <- numInClass[i]%%folds
#         seqVector <- rep(1:folds, min_reps)
#         if (spares > 0) 
#           seqVector <- c(seqVector, sample(1:folds, spares))
#         foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
#       }
#       else {
#         foldVector[which(y == names(numInClass)[i])] <- sample(1:folds, 
#                                                                size = numInClass[i])
#       }
#     }
#   }
#   else
#     foldVector <- seq(along = y)
#   out = split(seq(along = y), foldVector)
#   
#   lapply(out, function(x)permutation[x])
# }
# 
# 
# # Generating training and test sets from dataset 
# 
# cvPartitions = function(DS_name, times = 2, folds = 10, strat = FALSE){
#   DS = readData(DS_name, DSlocation)
#   name_extension = strsplit(DS_name, "[.]")[[1]]
#   
#   createFile = function(it){
#     CV = cvStrat(DS, fold = 10) # Folds with stratification
#     CV_without = list()
#     if(!strat){
#       CV_without = normalCVFolds(DS$y, k = 10)
#     }
#     
#     
#     cv_iter = function(cv_it){ 
#       
#       # Directory where data will be stored
#       new_dir = file.path(getwd(), "OUTPUT")
#       dir.create(new_dir, showWarnings = FALSE)
#       new_dir = file.path(getwd(), paste("OUTPUT", name_extension[1], sep = "/"))
#       dir.create(new_dir, showWarnings = FALSE)
#       
#       # Setting the workign diretory as the current dataset directory
#       setwd(new_dir)
#       
#       
#       # Names for training and test sets (with and without SMOGN)
#       train_name = paste(name_extension[1], it, "-", cv_it, ".tra", sep = "")
#       train_name_SMOGN = paste(name_extension[1], it, "X", "-", cv_it, ".tra", sep = "")
#       test_name = paste(name_extension[1], it, "-", cv_it, ".tst", sep = "")
#       train_name_without_strat = paste(name_extension[1], it,"N", "-", cv_it, ".tra", sep = "")
#       test_name_without_strat = paste(name_extension[1], it,"N", "-", cv_it, ".tst", sep = "")
#       train_name_SMOGN_without_strat = paste(name_extension[1], it,"XN", "-", cv_it, ".tra", sep = "")
# 
#       
#       # File with training set
#       write_in_file(DS[- CV[[cv_it]],], train_name)
#       
#       # File with training with SMOGN set
#       pc <- UBL::phi.control(as.numeric(DS[- CV[[cv_it]],]$y), method="extremes") # relevance parameters
#       SMOGN_res = DIBSRegress(y~., dat = DS[- CV[[cv_it]],], method = "extremes", 
#                               npts = pc$npts, control.pts = pc$control.pts, thr.rel = 0.8, pert = 0.01)
#       write_in_file(SMOGN_res, train_name_SMOGN)
#       
#       # File with test set
#       write_in_file(DS[CV[[cv_it]],], test_name)
#       
#       ## IF ADITIONAL NORMAL CV:
#       if(!strat){
#         # SMOGN without estratification:
#         pc_sin <- UBL::phi.control(as.numeric(DS[- CV_without[[cv_it]],]$y), method="extremes") # relevance parameters
#         SMOGN_res_sin = DIBSRegress(y~., dat = DS[- CV_without[[cv_it]],], method = "extremes", 
#                                 npts = pc_sin$npts, control.pts = pc_sin$control.pts, thr.rel = 0.8, pert = 0.01)
#         write_in_file(SMOGN_res, train_name_SMOGN_without_strat)
#         
#         # File with training random CV
#         write_in_file(DS[- CV_without[[cv_it]],], train_name_without_strat)
#         # File with test random CV
#         write_in_file(DS[CV_without[[cv_it]],], test_name_without_strat)
#       }
#       
#       setwd("../..")
#     }
#     lapply(1:folds, cv_iter)
#   }
#   lapply(1:times, createFile)
# }
# 
# 
# # Writing data in file
# write_in_file = function(DS, file_name){
#   cat(paste(toString(nrow(DS)),toString(ncol(DS)), "", sep = "\n"), file=file_name)
#   write.table(DS, file_name, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
# }
# 
# 
# # Obtaining percentages of rares in datasets
# getPercentages = function(DSnames, DSlocation){
#   percentages = c()
#   for (i in DSnames){
#     d = readData(i, DSlocation)
#     pc <- UBL::phi.control(as.numeric(d$y), method="extremes")
#     
#     p = UBL::phi(as.numeric(d$y), pc)
#     
#     porcentaje = length(p[p  >=0.8])/length(d$y) *100
#     percentages = rbind(percentages, c(i, porcentaje))
#     
#   }
#   
#   percentages
# }
# 
# 
# # Cross-validation
# normalCVFolds = function(y, kfold = 10){
#   # Randomly shuffling the data
#   permutacion = sample(length(y))
#   y <- y[sample(length(y))]
#   # Creating 10 equally sized folds
#   folds <- cut(seq(1,length(y)),breaks=kfold,labels=FALSE)
#   # Performing 10-fold cross-validation
#   cv = list()
#   for(i in 1:kfold){
#     # Segmenting the data by folds using the which() function 
#     testIndexes <- which(folds==i,arr.ind=TRUE)
#     # Using the test and training data partitions as desired...
#     cv[[i]] = permutacion[testIndexes]
#     names(cv[i]) = toString(i)
#   }
#   cv
# }

# Function to read datasets in KEEL format (*.dat).

readData = function(name, location = ""){
  
  r = read.csv(paste(".", location,"/",name, sep = ""), comment.char = "@", header = FALSE)
  colnames(r)[ncol(r)] = "y"
  r
  # add name class:
}


# Cross Validation with stratification (generating the folds)
cvStrat = function (DS, folds = 10) {   
  permutation = sample(1:nrow(DS))
  y = DS[permutation,]$y
  names(y) = permutation
  if (is.numeric(y)) {
    cuts <- floor(length(y)/folds)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (folds < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%folds
      if (min_reps > 0) {
        spares <- numInClass[i]%%folds
        seqVector <- rep(1:folds, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:folds, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:folds, 
                                                               size = numInClass[i])
      }
    }
  }
  else
    foldVector <- seq(along = y)
  out = split(seq(along = y), foldVector)
  
  lapply(out, function(x)permutation[x])
}


# Generating training and test sets from dataset 

cvPartitions = function(DS_name, times = 2, folds = 10, strat = FALSE){
  DS = readData(DS_name, DSlocation)
  name_extension = strsplit(DS_name, "[.]")[[1]]
  
  createFile = function(it){
    CV = cvStrat(DS, fold = 10) # Folds with stratification
    CV_without = list()
    if(!strat){
      CV_without = normalCVFolds(DS$y, k = 10)
    }
    
    
    cv_iter = function(cv_it){ 
      
      # Directory where data will be stored
      new_dir = file.path(getwd(), "OUTPUT")
      dir.create(new_dir, showWarnings = FALSE)
      new_dir = file.path(getwd(), paste("OUTPUT", name_extension[1], sep = "/"))
      dir.create(new_dir, showWarnings = FALSE)
      
      # Setting the workign diretory as the current dataset directory
      setwd(new_dir)
      
      
      # Names for training and test sets (with and without SMOGN)
      train_name = paste(name_extension[1], it, "-", cv_it, ".tra", sep = "")
      train_name_SMOGN = paste(name_extension[1], it, "X", "-", cv_it, ".tra", sep = "")
      test_name = paste(name_extension[1], it, "-", cv_it, ".tst", sep = "")
      train_name_without_strat = paste(name_extension[1], it,"N", "-", cv_it, ".tra", sep = "")
      test_name_without_strat = paste(name_extension[1], it,"N", "-", cv_it, ".tst", sep = "")
      train_name_SMOGN_without_strat = paste(name_extension[1], it,"XN", "-", cv_it, ".tra", sep = "")
      
      
      # File with training set
      write_in_file(DS[- CV[[cv_it]],], train_name)
      
      # File with training with SMOGN set
      pc <- UBL::phi.control(as.numeric(DS[- CV[[cv_it]],]$y), method="extremes") # relevance parameters
      SMOGN_res = DIBSRegress(y~., dat = DS[- CV[[cv_it]],], method = "extremes", 
                              npts = pc$npts, control.pts = pc$control.pts, thr.rel = 0.8, pert = 0.01)
      write_in_file(SMOGN_res, train_name_SMOGN)
      
      # File with test set
      write_in_file(DS[CV[[cv_it]],], test_name)
      
      ## IF ADITIONAL NORMAL CV:
      if(!strat){
        # SMOGN without estratification:
        pc_sin <- UBL::phi.control(as.numeric(DS[- CV_without[[cv_it]],]$y), method="extremes") # relevance parameters
        SMOGN_res_sin = DIBSRegress(y~., dat = DS[- CV_without[[cv_it]],], method = "extremes", 
                                    npts = pc_sin$npts, control.pts = pc_sin$control.pts, thr.rel = 0.8, pert = 0.01)
        write_in_file(SMOGN_res, train_name_SMOGN_without_strat)
        
        # File with training random CV
        write_in_file(DS[- CV_without[[cv_it]],], train_name_without_strat)
        # File with test random CV
        write_in_file(DS[CV_without[[cv_it]],], test_name_without_strat)
      }
      
      setwd("../..")
    }
    lapply(1:folds, cv_iter)
  }
  lapply(1:times, createFile)
}


# Writing data in file
write_in_file = function(DS, file_name){
  cat(paste(toString(nrow(DS)),toString(ncol(DS)), "", sep = "\n"), file=file_name)
  write.table(DS, file_name, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
}


# Obtaining percentages of rares in datasets
getPercentages = function(DSnames, DSlocation){
  percentages = c()
  for (i in DSnames){
    d = readData(i, DSlocation)
    pc <- UBL::phi.control(as.numeric(d$y), method="extremes")
    
    p = UBL::phi(as.numeric(d$y), pc)
    
    porcentaje = length(p[p  >=0.8])/length(d$y) *100
    percentages = rbind(percentages, c(i, porcentaje))
    
  }
  
  percentages
}


# Cross-validation
normalCVFolds = function(y, kfold = 10){
  # Randomly shuffling the data
  permutacion = sample(length(y))
  y <- y[sample(length(y))]
  # Creating 10 equally sized folds
  folds <- cut(seq(1,length(y)),breaks=kfold,labels=FALSE)
  # Performing 10-fold cross-validation
  cv = list()
  for(i in 1:kfold){
    # Segmenting the data by folds using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    # Using the test and training data partitions as desired...
    cv[[i]] = permutacion[testIndexes]
    names(cv[i]) = toString(i)
  }
  cv
}

