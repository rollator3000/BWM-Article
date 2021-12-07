"Script to evaluate the blockwise approach on the data with blockwise missingness.

  > Get the blocks that the test- & train-set have in common
  > Train a seperate RF on each of these train-blocks (only on the complete-case observations)
  > Each RF predicts for the test-set
  > Get an overall prediciton by doing an weighted average of the block-wise predicitons
  > Use the oob-AUC of the block-wise fitted RFs to determine their weight in the weighted average! 
"
# [0] SetWD, load packages, define fix variables and fuctions                ----
# 0-1 Set WD
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
# 0-5-1 Fit an RF to data, return the RF & its oob AUC
fit_RF_get_oob_AUC <- function(data) {
  "Fit an RF (with its standard settings) to 'data' ('ytarget' must be in there & is used as response).
   Only return the oob-AUC of the fit RF.
   
   Args:
    > data (data.frame): Data with at least two columns & observations. 
                         Must contain the column 'ytarget' & no missing values.
                         
   Return:
    > Return the AUC-metric that was calculated based on the oob-observations of the RF
  "
  # [0] Check Inputs
  # 0-1 'data' has to be a DF, with at least 2 observations & w/ columns
  assert_data_frame(data, any.missing = F, min.rows = 2, min.cols = 2)
  if (!'ytarget' %in% colnames(data)) stop("'data' must contain 'ytarget' as column")
  
  # [1] Fit RF on the data
  # 1-1 Train a RF on 'train'
  # --1 Convert the response to a factor
  data[,'ytarget'] <- as.factor(data[,'ytarget'])
  
  # --2 Create a formula to pass to the RF 
  #     (define response & use remaining variables as features)
  formula_all <- as.formula(paste('ytarget', " ~ ."))
  
  # --3 Fit the actual RF (only use standard-settings)
  RF <- rfsrc(formula = formula_all, data = data, samptype = "swr", 
              seed = 12345678, var.used = 'all.trees')
  
  # [2] Get the AUC based on the oob-observations & return it then
  # 2-1 Get the predicted probabilities for the oob-observations
  pred_prob_oob <- RF$predicted.oob[,'1']
  
  # 2-2 Compare the predicited class-prob. with the true classes & calc the AUC
  AUC <- pROC::auc(data$ytarget, pred_prob_oob, quiet = T)
  
  # 2-3 Return the AUC
  return(list("RF" = RF,
              "AUC" = AUC))
}

# 0-5-2 Function to evaluate the BW-Apprach
eval_bw_approach <- function(path = './Data/Raw/BLCA.Rda', frac_train = 0.75, split_seed = 1312,
                             block_seed_train = 1234, block_seed_test = 1312, train_pattern = 2, 
                             train_pattern_seed = 12, test_pattern = 2) {
  "Evaluate the BW-Approach on the data 'path' points to. 
   Get the blocks 'Train' & 'Test' have in common. Then train a seperate RF on these blocks 
   (use standard settings 'ntree', 'mtry' & 'min_node_size' & only observations w/o NA). Each of 
   these RFs predict on the test-set then. The final prediciton equals a weighted average of the
   block-wise predicitons. The weights in the weighted average equal the oob-auc of the block-wise
   fitted RFs.
   Finally return a DF with the the AUC, the Brier-Score and the standard metrics Precision, Recall, 
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
  
  # 1-2 Get the observed blocks of test- & train-set
  # --1 Test: Get the fully observed blocks and store their names in 'test_blocks' 
  test_blocks <- list()
  for (curr_block in train_test_bwm$Test$block_names) {
    
    # --1 Which Index has the current block
    curr_block_idx <- which(train_test_bwm$Test$block_names == curr_block)
    
    # --2 Get the corresponding columns to 'curr_block'
    curr_block_cols <- which(train_test_bwm$Test$block_index == curr_block_idx)
    
    # --3 Check whether the whole block is observed (no NA's)
    curr_block_observed <- all(!is.na(train_test_bwm$Test$data[,curr_block_cols] <= 0))
    
    # --4 If the current block is fully observed, add it to 'test_blocks'
    if (curr_block_observed) {
      test_blocks[[curr_block]] <- train_test_bwm$Test$data[,curr_block_cols]
    }
  }
  
  # --2 For each block in 'test_blocks', check whether train has observations in the corresponding 
  #     block. If this is the case, add the observed part of the block to 'train_blocks'
  train_blocks <- list()
  for (curr_test_block in names(test_blocks)) {
    
    # --1 Get the index of the current test-block from the train-set
    curr_train_block_idx <- which(train_test_bwm$Train$block_names == curr_test_block)
    
    # --2 Extract the corresponding columns from train for current train-block 
    curr_train_block_cols <- which(train_test_bwm$Train$block_index == curr_train_block_idx)
    
    # --3 Check whether the whole block has (at least some) observed values (no NA's), then we can use it!
    curr_train_block_partly_observed <- any(!is.na(train_test_bwm$Train$data[,curr_train_block_cols] <= 0))
    
    # --4 If the current_test_block is (partly) observed in the train, add it to 'train_blocks' &
    #     additionally add the corresponding response to the DF then
    if (curr_train_block_partly_observed) {
      
      #--4-1 Get the rows that are fully observed
      observed_rows <- which(rowSums(is.na(train_test_bwm$Train$data[,curr_train_block_cols])) == 0)
      
      # --4-2 Add the block with only observed values to 'train_blocks'
      train_blocks[[curr_test_block]] <- train_test_bwm$Train$data[observed_rows, curr_train_block_cols]
      
      # --4-3 Add the corresponding response variable
      train_blocks[[curr_test_block]]$ytarget <- train_test_bwm$Train$data$ytarget[observed_rows]
    }
  }
  
  # [2] Train & evaluate a seperate RF (with its standard-settings) on each block in 'train_blocks'
  # 2-1 Fit a seperate RF on each of the blocks in 'train_blocks'
  fitted_RFs <- list()
  for (curr_train_block in names(train_blocks)) {
    
    # --1 Fit an RF on the curr_train_block, save the DF & its oob AUC
    fitted_RFs[[curr_train_block]] <- fit_RF_get_oob_AUC(train_blocks[[curr_train_block]])
  }
  
  # 2-2 Get predicitons (prob. for class 1) for the test-set from each RF
  preds_test_set <- list()
  for (curr_block in names(train_blocks)) {
    
    # --1 Fit an RF on the curr_block of the test-set
    curr_pred <- predict(fitted_RFs[[curr_block]]$RF, 
                         test_blocks[[curr_block]])
    
    # --2 Save the predicted probs for class 1 to 'preds_test_set'
    preds_test_set[[curr_block]] <- curr_pred$predicted[,2]
  }
  
  # 2-3 Get the weighted average of the predicitons on the test-set
  # --1 Collect oob-AUCs of the blockwise fitted RFs
  oob_weights <- sapply(names(fitted_RFs), function(x) fitted_RFs[[x]]$AUC)
  
  # --2 Get the weighted sum of the predicted probabilities 
  weighted_predicitons <- c()
  for (curr_block in names(train_blocks)) {
    
    if (length(weighted_predicitons) == 0) {
      weighted_predicitons <- (preds_test_set[[curr_block]] * oob_weights[which(names(fitted_RFs) == curr_block)])
    } else {
      weighted_predicitons <- weighted_predicitons + (preds_test_set[[curr_block]] * oob_weights[which(names(fitted_RFs) == curr_block)]) 
    }
  }
  
  # --3 Scale them down by dividing with the sum of 'oob_weights'
  weighted_predicitons <- weighted_predicitons / sum(oob_weights)
  
  # --4 Get the predicted class
  classes_predicted <- factor(as.numeric(weighted_predicitons >= 0.5), levels = c(0, 1))
  
  # [3] Calculate the metrics based on the true & predicted labels
  # 3-1  Confusion Matrix & all corresponding metrics (Acc, F1, Precision, ....)
  metrics_1 <- caret::confusionMatrix(classes_predicted, 
                                      factor(train_test_bwm$Test$data$ytarget, 
                                             levels = c(0, 1)),
                                      positive = "1")
  
  # 3-2 Calculate the AUC
  AUC <- pROC::auc(factor(train_test_bwm$Test$data$ytarget, levels = c(0, 1)), 
                   weighted_predicitons, quiet = T)
  
  # 3-3 Calculate the Brier-Score
  brier <- mean((weighted_predicitons - train_test_bwm$Test$data$ytarget)  ^ 2)
  
  # [4] Return the results as DF
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
                    "common_blocks"      = paste(names(train_blocks), collapse = ' - '),
                    "AUC"                = AUC,
                    "Accuracy"           = metrics_1$overall['Accuracy'], 
                    "Sensitivity"        = metrics_1$byClass['Sensitivity'], 
                    "Specificity"        = metrics_1$byClass['Specificity'], 
                    "Precision"          = metrics_1$byClass['Precision'], 
                    "Recall"             = metrics_1$byClass['Recall'], 
                    "F1"                 = metrics_1$byClass['F1'], 
                    "BrierScore"         = brier))
}

# [1] Run the experiments                                                    ----
# 1-1 Initalize a empty DF to store the results
BW_res <- data.frame()

# 1-2 Define a list with the paths to the availabe DFs
df_paths <- paste0("./Data/Raw/", list.files("./Data/Raw/"))

# 1-3 Loop over all the possible settings for the evaluation of the SB-Approach
#     each setting is evaluated 5-times!
for (curr_path in df_paths) {
  for (curr_train_pattern in c(1, 2, 3, 4, 5)) {
    for (curr_test_pattern in c(1, 2, 3, 4)) {
      for (curr_repetition in c(1, 2, 3, 4, 5)) {
        
        cat('-----------------------------------------------\n',
            "Current Path:          >", curr_path, '\n',
            "Current Train Pattern: >", curr_train_pattern, '\n',
            "Current Test Patter:   >", curr_test_pattern, '\n',
            "Current Repetition:    >", curr_repetition, '\n')
        
        # Set the seed for the 'split'
        curr_split_seed = 12345678 + curr_repetition
        
        # Set the seed for the shuffling of the blocks of train & test
        curr_block_seed_train = 1234567 + curr_repetition
        curr_block_seed_test  = 7654321 + curr_repetition
        
        # Set the seed for the train_pattern (shuffling of observations)
        curr_train_pattern_seed = 12345 + curr_repetition
        
        # Run the evaluation with current settings
        curr_res <- tryCatch(eval_bw_approach(path               = curr_path, 
                                              frac_train         = 0.75, 
                                              split_seed         = curr_split_seed,
                                              block_seed_train   = curr_block_seed_train, 
                                              block_seed_test    = curr_block_seed_test,
                                              train_pattern      = curr_train_pattern,
                                              train_pattern_seed = curr_train_pattern_seed, 
                                              test_pattern       = curr_test_pattern),
                             error = function(c) {
                               data.frame("path"               = curr_path, 
                                          "frac_train"         = 0.75, 
                                          "split_seed"         = curr_split_seed, 
                                          "block_seed_train"   = curr_block_seed_test,
                                          "block_seed_test"    = curr_block_seed_test, 
                                          "block_order_train_for_BWM" = '---',
                                          "block_order_test_for_BWM"  = '---',
                                          "train_pattern"      = curr_train_pattern, 
                                          "train_pattern_seed" = curr_train_pattern_seed, 
                                          "test_pattern"       = curr_test_pattern,
                                          "common_blocks"      = "---",
                                          "AUC"                = '---',
                                          "Accuracy"           = '---', 
                                          "Sensitivity"        = '---', 
                                          "Specificity"        = '---', 
                                          "Precision"          = '---', 
                                          "Recall"             = '---', 
                                          "F1"                 = '---', 
                                          "BrierScore"         = '---')
                             }
        ) 
        
        # Add the curr_repetition to 'curr_res', before adding it to 'SB_res'
        curr_res$repetition <- curr_repetition
        curr_res$approach   <- 'Blockwise'
        
        # Add the results of the setting to 'SB_res' & save it
        BW_res <- rbind(BW_res, curr_res)
        write.csv(BW_res, './Docs/Evaluation_Results/BW_Approach/BW_Eval.csv')
      }
    }
  }
}
