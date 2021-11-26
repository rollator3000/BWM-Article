"Script to evaluate the Single-Block-Case approach on data with blockwise missingness

  > Use only a single block to train a RF and predict on the test-set then 
  > 1. Remove all blocks from the training-set that are not featured in the test-set
  > 2. Train a seperate RF on each of the blocks
  > 3. Rate the performance of the RF with the AUC based on the out-of-bag predicitons
  > 4. The RF with the best out-of-bag AUC is then used to predict on the test-set
"
# [0] SetWD, load packages, define fix variables and fuctions                ----
# 0-1 Set WD (currently out-commented, as we need to load the script)
setwd("/Users/frederik/Desktop/BWM-Article/")             # Mac
setwd("C:/Users/kuche/Desktop/BWM-Paper")                 # Windows
setwd("/dss/dsshome1/lxc0B/ru68kiq3/Project/BWM-Article") # Server

# 0-2 Load packages
library(checkmate)
library(randomForestSRC)
library(parallel)
library(doParallel)
library(caret)
library(pROC)

# 0-3 Define fixed variables
# 0-3-1 Define amount of usable cores (parallel computing)
detectCores()
registerDoParallel(cores = 2)

# 0-4 Load functions from 'code/01_Create_BWM_Pattern"
source("./Code/01_Create_BWM_Pattern.R")

# 0-5 Define functions

# 0-5-1 Function to fit an RF & return its oob-AUC
fit_RF_get_oob_AUC <- function(data) {
  "Fit an RF (with its standard settings) to 'data' ('ytarget' must be in there & is used as response).
   Only return the oob-AUC of the fit RF.
   
   Args:
    > data (data.frame): Data with at least two columns & observations. 
                         Must contain the column 'ytarget'
                         
   Return:
    > Return the AUC-metric that was calculated based on the oob-observations of the RF
  "
  # [0] Check Inputs
  
  # [1] Fit RF on the data
  
  # [2] Get the AUC based on the oob-observations & return it then
}

# 0-5-2 Evaluate a RF with the complete-case approach & get its metrics
path = './Data/Raw/BLCA.Rda'
frac_train = 0.75
split_seed = 1312
block_seed_train = 1234
block_seed_test = 1312
train_pattern = 2 
train_pattern_seed = 12
test_pattern = 2
eval_sb_appr <- function(path = './Data/Raw/BLCA.Rda', frac_train = 0.75, split_seed = 1312,
                         block_seed_train = 1234, block_seed_test = 1312, train_pattern = 2, 
                         train_pattern_seed = 12, test_pattern = 2) {
  "Evaluate the SB-Approach on the data 'path' points to.
   On each block that the test- & train-set have in commom, a seperate RF is trained & evaluated 
   with the AUC based on the out-of-bag predictions - the RFs are trained with their standard settings
   (e.g. 'ntree', 'mtry' & 'min_node_size'). The block that leads to the best out-of-bag AUC is then
   used to predict on the test-set.
   Finally return a DF with the the AUC, the Brier-Score and  the standard metrics Precision, Recall,
   Sensitivity, Specificity, F-1 Score & Accuracy + all the settings for the evaluation 
   (e.g. path, seeds, train_pattern, settings for RF, block_order, ...).
   
   Args:
      > path               (str): Path to a dataset - must contain 'Data/Raw'
      > frac_train       (float): Fraction of observations for the train-set - ]0;1[
      > split_seed         (int): Seed for the split of the data to train & test
      > block_seed_train   (int): Seed for the shuffeling of the block-order in train
      > block_seed_test    (int): Seed for the shuffeling of the block-order in test
      > train_pattern      (int): Seed for the induction of the pattern for train
                                  (obs. are assigned to different folds!)
      > train_pattern_seed (int): Pattern to induce into train (1, 2, 3, 4, 5)
      > test_pattern       (int): Pattern to induce into test (1, 2, 3, 4)
   
   Return:
      > A DF with the settings of the experiment (path to the data, train pattern, ...)
        as well as the settings of the RF (ntree, mtry, ...) and the results
        of the evaluation (AUC; Brier-Score; Accuracy)
  "
  # [0] Check Inputs
  #     --> All arguments are checked in the functions 'get_train_test()' &
  #         'get_predicition()' that werde loaded from 'Code/01_Create_BWM_Pattern.R'
  
  # [1] Load & prepare the data 
  # 1-1 Load the data from 'path', split it to test- & train & induce block-wise 
  #     missingness into both of them according to 'train_pattern' & 'test_pattern'
  train_test_bwm <- get_train_test(path = path,                             # Path to the data
                                   frac_train = frac_train,                 # Fraction of data used for Training (rest for test)
                                   split_seed = split_seed,                 # Seed for the split of the data into test- & train
                                   block_seed_train = block_seed_train,     # Seed to shuffle the block-order in train
                                   block_seed_test = block_seed_test,       # Seed to shuffle the block-order in test
                                   train_pattern = train_pattern,           # Pattern to introduce to train
                                   train_pattern_seed = train_pattern_seed, # Seed for the introduction of the BWM into train
                                   test_pattern = test_pattern)             # Pattern for the test-set
  
  # 1-2 Extract the observed blocks from test & store the corresponding DFs seperatly in the list 
  #     'test_blocks' (each fully observed block is an own entry in the list w/ its block-name)
  test_blocks <- list()
  for (curr_block in train_test_bwm$Test$block_names) {
    
    # --1 Which Index has the current block
    curr_block_idx <- which(train_test_bwm$Test$block_names == curr_block)
    
    # --2 Get the corresponding columsn to 'curr_block'
    curr_block_cols <- which(train_test_bwm$Test$block_index == curr_block_idx)
    
    # --3 Check whether the whole block is observed (no NA's)
    curr_block_observed <- all(!is.na(train_test_bwm$Test$data[,curr_block_cols] <= 0))
    
    # --4 If the current block is fully observed, add it to 'test_blocks' & add the response as 
    #     additional column to the DF then
    if (curr_block_observed) {
      test_blocks[[curr_block]] <- train_test_bwm$Test$data[,curr_block_cols]
      test_blocks[[curr_block]]$ytarget <- train_test_bwm$Test$data$ytarget
    }
  }
  
  # 1-3 Extract the blocks from train-set that are fully observed in the test-set & store the 
  #     corresponding DFs in the list 'train_blocks' 
  train_blocks <- list()
  for (curr_test_block in names(test_blocks)) {
    
    # --1 Get the index of the current test-block from the train-set
    curr_train_block_idx <- which(train_test_bwm$Train$block_names == curr_test_block)
    
    # --2 Extract the corresponding columns from train for current train-block 
    curr_train_block_cols <- which(train_test_bwm$Train$block_index == curr_train_block_idx)
    
    # --3 Check whether the whole block has (at least some) observed values (no NA's), then we can use it!
    curr_train_block_partly_observed <- any(!is.na(train_test_bwm$Train$data[,curr_train_block_cols] <= 0))
    
    # --4 If the current_test_block is patzly observed in the train-set, add it to 'train_blocks'
    #     & additionaly add the corresponding response to DF then
    if (curr_train_block_partly_observed) {
      
      #--4-1 Get the rows that are fully observed
      observed_rows <- which(rowSums(is.na(train_test_bwm$Train$data[,curr_train_block_cols])) == 0)
      
      # --4-2 Add the block with only observed values to 'train_blocks'
      train_blocks[[curr_test_block]] <- train_test_bwm$Train$data[observed_rows, curr_train_block_cols]
      
      # --4-3 Add the corresponding response variable
      train_blocks[[curr_test_block]]$ytarget <- train_test_bwm$Train$data$ytarget[observed_rows]
    }
  }
  
  # 1-4 If 'train_blocks' is empty the SB approach can not be applied & return an DF w/o metrics
  if (length(names(train_blocks)) <= 0) {
    return(data.frame("path"               = curr_path, 
                      "frac_train"         = 0.75, 
                      "split_seed"         = curr_split_seed, 
                      "block_seed_train"   = curr_block_seed_test,
                      "block_seed_test"    = curr_block_seed_test, 
                      "block_order_train_for_BWM" = '---',
                      "block_order_test_for_BWM"  = '---',
                      "train_pattern"      = curr_train_pattern, 
                      "train_pattern_seed" = curr_train_pattern_seed, 
                      "test_pattern"       = curr_test_pattern, 
                      "ntree"              = '---', 
                      "mtry"               = '---', 
                      "min_node_size"      = '---', 
                      "AUC"                = '---',
                      "Accuracy"           = '---', 
                      "Sensitivity"        = '---', 
                      "Specificity"        = '---', 
                      "Precision"          = '---', 
                      "Recall"             = '---', 
                      "F1"                 = '---', 
                      "BrierScore"         = '---'))
  }
  
  # [2] Train is not empty: Fit an RF on each block of the train-set & evaluate it with the oob-AUC
  # 2-1 Loop over each block and fit an RF on it seperatly & get its oob-AUC
  for (names(train_blocks)) {
    
  }
  
  
  
  # ----------- STOPPED HERE -----------------------------------------------------------------------
  
  # , and we& evaluate a RF (with its standard-settings) 
  # 2-1 Get predictions for the test-set from a RF that is fitted with its 
  #     standard settings to the processed train-set
  preds_test_set <- get_predicition(train = train_test_bwm$Train$data, 
                                    test = train_test_bwm$Test$data)
  
  # 2-2 Calculate the metrics based on the true & predicted labels
  # 2-2-1  Confusion Matrix & all corresponding metrics (Acc, F1, Precision, ....)
  metrics_1 <- caret::confusionMatrix(preds_test_set$pred_classes, 
                                      factor(train_test_bwm$Test$data$ytarget, 
                                             levels = c(0, 1)),
                                      positive = "1")
  
  # 2-2-2 Calculate the AUC
  AUC <- pROC::auc(factor(train_test_bwm$Test$data$ytarget, levels = c(0, 1)), 
                   preds_test_set$pred_prob_pos_class, quiet = T)
  
  # 2-2-3 Calculate the Brier-Score
  brier <- mean((preds_test_set$pred_prob_pos_class - train_test_bwm$Test$data$ytarget)  ^ 2)
  
  # [3] Return the results as DF
  return(data.frame("path"               = path, 
                    "frac_train"         = frac_train, 
                    "split_seed"         = split_seed, 
                    "block_seed_train"   = block_seed_train,
                    "block_seed_test"    = block_seed_test, 
                    "block_order_train_for_BWM" = paste(train_test_bwm$Train$block_names, collapse = ' - '),
                    "block_order_test_for_BWM"  = paste(train_test_bwm$Test$block_names, collapse = ' - '),
                    "train_pattern"      = train_pattern, 
                    "train_pattern_seed" = train_pattern_seed, 
                    "test_pattern"       = test_pattern, 
                    "ntree"              = preds_test_set$RF_ntree, 
                    "mtry"               = preds_test_set$RF_mtry, 
                    "min_node_size"      = preds_test_set$RF_min_node_size, 
                    "AUC"                = AUC,
                    "Accuracy"           = metrics_1$overall['Accuracy'], 
                    "Sensitivity"        = metrics_1$byClass['Sensitivity'], 
                    "Specificity"        = metrics_1$byClass['Specificity'], 
                    "Precision"          = metrics_1$byClass['Precision'], 
                    "Recall"             = metrics_1$byClass['Recall'], 
                    "F1"                 = metrics_1$byClass['F1'], 
                    "BrierScore"         = brier))
}