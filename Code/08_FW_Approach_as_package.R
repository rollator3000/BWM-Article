"Functions for R. Hornung for the implementation of a R-Package for the FW-RF

 FW-Approach:
  > Fit a seperate RF on each fold of the train-set
  > Prune these foldwise fitted RFs in regard to the test-set
  > Internally evaluate pruned FW-RF with their oob-AUC 
  > Each RF is then used to predict on the observations of the test-set 
  > Create an overall prediciton by calculating an weighted average of the fold-wise
    predicitons - as weights we use the oob-AUC of the pruned FW-RFs
  > Evaluate the results with common metrics (AUC, Accuracy, F-1 Score, ...)
"
# [0] SetWD, load packages, define variables and functions                  ----
# 0-1 Set WD
setwd("/Users/frederik/Desktop/BWM-Article/")             # Mac
setwd("C:/Users/kuche/Desktop/BWM-Paper")                 # Windows
setwd("/dss/dsshome1/lxc0B/ru68kiq3/Project/BWM-Article") # Server

# 0-2 Load packages
library(checkmate)
library(randomForestSRC)
library(caret)
library(pROC)
library(doParallel)

# 0-3 Define fix variables

# 0-4 Define functions
# 0-4-1 Load functions from 'code/01_Create_BWM_Pattern"
source("./Code/01_Create_BWM_Pattern.R")

# 0-4-2 Load functions so we can fit and prune a RF 
#       (not possible with 'randomForestSRC'-library...)
source("./code/06_1_simpleRF_adaption.R")

# 0-4-3 Check whether all trees are grown correctly
all_trees_grown_correctly <- function(trees) {
  "Check, whether 'trees', were grown correctly & if not grow these trees again,
   as long, as they are grown correctly! 
      --> Growning not correctly: No 'childNodeIDs', no split variables etc...
  
   Args:
      trees (list) : list filled with object of the class 'Tree'! 
                     For each object in there check, whether it was grown 
                     correctly (has childNodeIDs) - if not grown it again!
                     
   Return: 
      list of trees, where all of these trees were grown correctly
      --> each tree has at least one 'childNodeID'
  "
  # [0] Check Inputs  ----------------------------------------------------------
  # 0-1 All Objects in 'trees' of class 'Tree'
  trees_classes <- sapply(trees, function(x) class(x))
  if (any(!grepl("Tree", trees_classes))) {
    stop("Not all Objects in 'trees' are of class 'Tree'")
  }
  
  # [1] Get the entrance of the objects, that miss child node IDs  -------------
  wrong_trees  <- unlist(lapply(1:length(trees), 
                                FUN = function(x) {
                                  if (length(trees[[x]]$child_nodeIDs) == 0) x
                                }))
  
  # [2] Regrow the trees, that were not grown correctly  -----------------------
  #     If there are any trees not grown correctly, grow them again until all
  #     of the trees were grown correctly!
  while (length(wrong_trees) > 0) {
    
    # grow the errours trees again
    trees[wrong_trees] <- lapply(trees[wrong_trees], 
                                 function(x) {
                                   x$grow(replace = TRUE)
                                   x
                                 })
    
    # check whether any of the trees is not grown correctly!
    wrong_trees  <- unlist(lapply(1:length(trees), 
                                  FUN = function(x) {
                                    if (length(trees[[x]]$child_nodeIDs) == 0) x
                                  }))
  }
  
  # [3] Return the correclty grown trees  --------------------------------------
  return(trees)
}

# 0-4-4 Calculate the oob AUC of single Forest
get_oob_AUC <- function(trees) {
  "Calculate OOB AUC of a list of trees! 
   For this we go  through all OOB Predictions and obtain aggregated predicitons
   from all trees, that have the same observation as out-of-bag!
   In case the metrics are 'NA' [not defined] they are repalced by 0 [worst value]
    
    Args: 
      - trees (list)        : list filled w/ objects of class 'Tree'
      
    Return:
      - Average oob-Acc & oob-F1 metric for 'trees'!
  "
  # [0] Check Input ------------------------------------------------------------
  # 0-1 Make sure 'trees' is a list filled with 'Trees'
  assert_list(trees, min.len = 1, any.missing = FALSE)
  if (any(sapply(trees, function(x) 'Tree' %in% class(x)))) {
    stop("not all elements in 'trees' are of class 'Tree'")
  }
  
  # [1] Get the OOB Predicitons ------------------------------------------------
  # 1-1 Get the trees, that are usable - not pruned in the first split_variable!
  usable_trees <- sapply(1:length(trees), function(x) {
    
    # Check whether the first split_var was pruned!
    if (trees[[x]]$child_nodeIDs[[1]][1] != "pruned") {
      return(x)
    } 
  })
  
  # 1-1-1 If all the trees were pruned in the first node, we can not do
  #       any OOB predictions --> OOB-Accuracy = 0
  if (length(usable_trees) < 1) {
    return(0)
  } else {
    usable_trees <- unlist(usable_trees)
  }
  
  # 1-2 For all usable trees [not pruned in first split variable] get the IDs
  #     of the OOB observations. Then collect all and only keep unique IDs!
  oob_ids_all_trees <- sapply(usable_trees, function(x) trees[[x]]$oob_sampleIDs)
  unique_oob_ids    <- unique(unlist(oob_ids_all_trees))
  
  # 1-3 Loop over all 'unique_oob_ids' and get the oob predictions from all 
  #     trees, that have the same OOB observation!
  all_oob_preds_class0 <- c()
  for (curr_oob in unique_oob_ids) {
    
    # 1-3-1 Get all trees that have 'curr_oob' as OOB observation!
    trees_same_oob <- unlist(sapply(usable_trees, function(x) {
      if (curr_oob %in% trees[[x]]$oob_sampleIDs) x
    }))
    
    # 1-3-2 Get the feas of the observation that is OOB for 'trees_same_oob' 
    curr_oob_feas <- Data$new(data = trees[[1]]$data$data[curr_oob,])
    
    # 1-3-3 Get a Prediciton for the 'curr_oob' from all trees!
    predicted_probs <- sapply(trees_same_oob, 
                              function(x) trees[[x]]$predict(curr_oob_feas))
    
    # 1-3-4 Aggregate the predictions from the different trees and
    #       Get the probability for class 0! [First Row is class 0]
    all_oob_preds_class0 <- c(all_oob_preds_class0, 
                              sum(predicted_probs[1,]) / length(trees_same_oob))
  }
  
  # [2] Compare predicted probabilites with the true classes 
  # 2-1 Convert 'all_oob_preds_class0' to 'all_oob_preds_class1'
  all_oob_preds_class1 <- ((all_oob_preds_class0 - 1) * - 1)
  
  # 2-2 Get the AUC
  AUC <- tryCatch(expr = pROC::auc(trees[[1]]$data$data[unique_oob_ids, 1], 
                                   all_oob_preds_class1, quiet = T),
                  error = function(c) 0)
  
  # 2-3 Return it
  return(AUC)
}


data          = train_test_bwm$Train$data
folds         = train_test_bwm$Train$fold_index
num_trees     = 500
mtry          = NULL
min_node_size = 1


# 0-4-5 Train a foldwise RF
train <- function(data, folds, num_trees = 500, mtry = NULL, min_node_size = 1) {
  " Train a seperate RF on each fold available in 'data', such that the final RF
    consists of as many FW-RFs as the data contains unique folds.
    
    Args:
      data (data.frame)  : Dataframe with dimensions n*p. Must contain a binary 
                           factor-column 'ytarget' (0 = neg. class / 1 = pos. class)
      folds (vec)        : Vector of length 'n' filled with integers. Indicates the 
                           for each row in 'data' to which fold it belongs.
                           (!! Obs. in the same fold need to be observed in the 
                               features !!)
      num_trees (int)    : Amount of trees that are used to grow a foldwise RF
      mtry (int)         : Amount of variables to be checked as split-variables
                           at every split. Default = ceiling(sqrt(p))
      min_node_size (int): Amount of Observations, that, at least, need 
                           to be in a terminal node!
                           
    Return: Return a list of fold-wise fitted RFs
  "
  # [0] Check inputs
  # 0-1 'data' has to be a data.frame & contain 'ytarget' as binary factor variable
  assert_data_frame(data)
  if (!('ytarget' %in% colnames(data))) {
    stop("'data' must contain 'ytarget' as column")
  } 
  assert_factor(data$ytarget, levels = c('0', '1')
  
  # 0-2 'folds' must be of the same length as 'data' & must only contain integers
  if (nrow(data) != length(folds)) stop("'folds' & 'data' differ in lenght")
  
  # 0-3 'num_trees' & 'min_node_size' must be integers >= 1
  assert_int(num_trees, lower = 1)
  assert_int(min_node_size, lower = 1)
  
  # 0-4 'mtry' must be an int <= amount of cols in data - if not 'NULL'
  if is.null(mtry) {
    mtry = ceiling(sqrt(ncol(data)))
  } else {
    assert_int(mtry, upper = ncol(data), lower = 1)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}