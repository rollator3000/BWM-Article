# Create an temporary script to track down the error 'unkown factors' in the 
# Complete-Case-Apporoach

# [0] SetWD, load packages, define variables and functions
# 0-1 Set WD
setwd("/Users/frederik/Desktop/BWM-Article/")             # Mac
setwd("C:/Users/kuche/Desktop/BWM-Paper")                 # Windows

# 0-2 Load packages
library(checkmate)
library(randomForestSRC)
library(caret)
library(pROC)

# 0-3 Define variables

# 0-4 Define functions
# 0-4-1 Load functions from 'code/01_Create_BWM_Pattern"
source("./Code/01_Create_BWM_Pattern.R")

# 0-4-2 Function to fit a RF & return its OOB-AUC
fit_RF_get_oob_AUC <- function(data) {
  "Fit an RF (with its standard settings) to 'data' ('ytarget' must be in there &
   is used as response) - only return the oob-AUC of the fit RF.
   
   Args:
    > data (data.frame): Data with at least two columns & observations. 
                         Must contain the column 'ytarget' & no missing values.
                         
   Return:
    > Return the AUC-metric that was calculated based on the oob-observations of the RF
  "
  # [0] Check Inputs
  # 0-1 'data' has to be a DF, with at least 2 observations & w/ columns
  assert_data_frame(data, any.missing = F, min.rows = 2, min.cols = 2)
  
  # 0-2 'data' must contain 'ytarget' as column
  if (!'ytarget' %in% colnames(data)) stop("'data' must contain 'ytarget' as column")

  # [1] Fit RF on the data
  # 1-1 Train a RF
  # --1 Create a formula to pass to the RF 
  #     (define response & use remaining variables as features)
  formula_all <- as.formula(paste('ytarget', " ~ ."))
  
  # --2 Fit the actual RF (only use standard-settings)
  RF <- rfsrc(formula = formula_all, data = data, samptype = "swr", 
              seed = 12345678, var.used = 'all.trees')
  
  # [2] Get the AUC based on the oob-observations & return it then
  # 2-1 Get the predicted probabilities for the oob-observations
  pred_prob_oob <- RF$predicted.oob[,'1']
  
  # 2-2 Compare the predicted class-prob. with the true classes & calc the AUC
  #  -> in case RF predict all OOB as 0/ 1 error will arise -> set AUC to 0
  AUC <- tryCatch(expr = pROC::auc(data$ytarget, pred_prob_oob, quiet = T),
                  error = function(c) 0)
  
  # 2-3 Return the AUC
  return(AUC)
}

# 0-4-3 Function to fit an RF on train & predict on test then --> get all the metrics
get_predicition <- function(train, test) {
  " Fit a RF on 'train' and create predicitons for 'test' then 

    Args:
      - train  (DF): DF that only contains variables also availabe for 'test'.
                     Must not contain any NA values! 
      - test   (DF): DF that is completly observed in all observations
                     
    Return:
      - List with: > 'pred_classes' = predicted class for each observation in 'test'
                   > 'pred_prob_pos_class' = predicted probability for each obs. 
                                             to be in class 1
                   > settings of the RF for 'mtry', 'min_node_size' & 'ntree'
  "
  # [0] Check Inputs
  # 0-1 'train' & 'test' must be dataframes w/o missing values
  assert_data_frame(train, any.missing = FALSE)
  assert_data_frame(test, min.rows = 1, any.missing = FALSE)
  
  # 0-2 'train' must not contain any col-names not available in 'test' & vic versa
  if (!all((colnames(train) %in% colnames(test)))) {
    stop("Train-Set has different features than the Test-Set!")
  }
  
  if (!all((colnames(test) %in% colnames(train)))) {
    stop("Test-Set has different features than the Train-Set!")
  }
  
  # [1] Train a RF & create predictions for the test-set
  # 1-1 Train a RF on 'train'
  # --1 Create a formula to pass to the RF 
  #     (define response & use remaining variables as features)
  formula_all <- as.formula(paste('ytarget', " ~ ."))
  
  # --2 Fit the RF (use standard-settings)
  RF <- rfsrc(formula = formula_all, data = train, samptype = "swr", 
              seed = 12345678, var.used = 'all.trees')
  
  # 1-2 Get predictions for the test-set
  predicitons <- predict(RF, test)
  
  # 1-3 Return the predicted classes & the predicted probabilities for class '1'
  #     as well as the settings of the RF
  return(list('pred_classes'        = predicitons$class,
              'pred_prob_pos_class' = predicitons$predicted[,'1'],
              'RF_ntree'            = RF$ntree,
              'RF_mtry'             = RF$mtry,
              'RF_min_node_size'    = RF$nodesize))
}

# [1] Evaluate the CC-Approach
# 1-1 Set arguments
path          = './Data/Raw/ESCA.Rda'
frac_train    = 0.75
train_pattern = 2
test_pattern  = 4
int_seed      = 7315704

# 1-2 Get seeds for the split of train- & test-set, shuffeling block-order, ...
set.seed(int_seed)    
seeds <- round(runif(4, 0, 100000))

#     1. Seed to split data into test & train
curr_split_seed <- seeds[1]
split_seed      <- curr_split_seed

#     2. Seed to shuffle the block order in 'train'
curr_block_seed_train <- seeds[2]
block_seed_train      <- curr_block_seed_train

#     3. Seed to shuffle the block order in 'test'
curr_block_seed_test <- seeds[3]
block_seed_test      <- curr_block_seed_test

#     4. Seed for the train-pattern (assignment of obs. in train to folds)
curr_train_pattern_seed <- seeds[4]
train_pattern_seed      <- curr_train_pattern_seed

##### SB APPROACH ------------------------------------------------------RAW  ---


# [0] Check Inputs
#     --> All arguments are checked in the functions 'get_train_test()' &
#         'get_predicition()' that were loaded from 'Code/01_Create_BWM_Pattern.R'

# [1] Load the data & prepare them for the SB-Approach
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

# 1-2 Extract the observed blocks from test & store the corresponding DFs separately in the list 
#     'test_blocks' (each fully observed block is an own entry in the list w/ its block-name)
test_blocks <- list()
for (curr_block in train_test_bwm$Test$block_names) {
  
  # --1 Which Index has the current block
  curr_block_idx <- which(train_test_bwm$Test$block_names == curr_block)
  
  # --2 Get the corresponding columns to 'curr_block'
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

# 1-3 For each block in 'test_blocks', check whether train has observations in the corresponding 
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

# 1-4 If 'train_blocks' is empty the SB approach can not be applied & return an DF w/o metrics
if (length(names(train_blocks)) <= 0) {
  return(data.frame("path"               = curr_path, 
                    "frac_train"         = 0.75, 
                    "split_seed"         = curr_split_seed, 
                    "block_seed_train"   = curr_block_seed_train,
                    "block_seed_test"    = curr_block_seed_test, 
                    "block_order_train_for_BWM" = paste(train_test_bwm$Train$block_names, collapse = ' - '),
                    "block_order_test_for_BWM"  = paste(train_test_bwm$Train$block_names, collapse = ' - '),
                    "train_pattern"      = curr_train_pattern, 
                    "train_pattern_seed" = curr_train_pattern_seed, 
                    "test_pattern"       = curr_test_pattern,
                    "common_blocks"      = "---",
                    "block_best_oob"     = "---",
                    "block_best_oob_auc" = "---",
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

# [2] If 'train_blocks' is not empty: Fit an RF on each block of the train-set & evaluate them 
#     with their OOB-AUC
# 2-1 Loop over each block in 'train_blocks', fit an RF on it & get the corresponding oob-AUC
res_df <- data.frame('block' = character(), 
                     'auc'   = numeric())

for (curr_train_block in names(train_blocks)) {
  
  # --1 Extract the data of the 'curr_train_block' from 'train_blocks'
  curr_data <- train_blocks[[curr_train_block]]
  
  # --2 Get the oob-AUC of a RF fit on 'curr_data'
  curr_AUC <- fit_RF_get_oob_AUC(data = curr_data)
  
  # --3 Collect 'block_name' & corresponding AUC
  res_df[nrow(res_df) + 1, ] <- c(curr_train_block, curr_AUC)
}

# 2-2 Get the name of the block that led to the highest oob-AUC
best_block <- res_df$block[which(res_df$auc == max(res_df$auc))]

# 2-3 In case 'best_block' has more than one entrance, randomly sample one
#     (this might happen due to multiple RFs with the same oob-AUC in 2-1)
set.seed(block_seed_train)
if (length(best_block) > 1) best_block <- sample(best_block, 1)

# [3] Get predictions for the test-set (based on 'best_block') & get the corresponding metrics
# 3-1 Train a RF on the 'best_block' of train-set & use it to predict on the test-set then!
preds_test_set <- get_predicition(train = train_blocks[[best_block]],
                                  test = test_blocks[[best_block]])

# 3-2 Calculate the metrics based on the true & predicted labels
# 3-2-1  Confusion Matrix & all corresponding metrics (Acc, F1, Precision, ....)
metrics_1 <- caret::confusionMatrix(preds_test_set$pred_classes, 
                                    train_test_bwm$Test$data$ytarget, 
                                    positive = "1")

# 3-2-2 Calculate the AUC
AUC <- pROC::auc(train_test_bwm$Test$data$ytarget, 
                 preds_test_set$pred_prob_pos_class, quiet = T)

# 3-2-3 Calculate the Brier-Score
brier <- mean((preds_test_set$pred_prob_pos_class - as.numeric(levels(train_test_bwm$Test$data$ytarget))[train_test_bwm$Test$data$ytarget]) ^ 2)

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
                  "common_blocks"      = paste(names(train_blocks), collapse = ' - '),
                  "block_best_oob"     = best_block,
                  "block_best_oob_auc" = res_df$auc[res_df$block == best_block],
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