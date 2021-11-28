"Script to evaluate the Imputation approach on data with blockwise missingness

  > 1. Impute thes missing values in the train-set with TOMBI
  > 2. Remove all blocks from train that are not available in test
  > 3. Fit a RF on the remaining train-set
  > 4. Use this RF to predict on the test-set & get the corresponding metrics
  
  > Functions for the imputation with TOMBI have been supplied by Dr. Roman Hornung.
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
# 0-5-1 Altered TOBMI function with greater speed than the original TOMBI
TOBMIfast <- function(x = cpg, y = exp) {
  
  ##Calculating the distances among un-/complete cases using auxiliary dataset
  dist.matrix <- as.matrix(dist( x ))
  
  ##Neighbors list for every uncomplete cases
  missing_num <- length(which(complete.cases(y) == F)) 
  donors <- list()
  for(i in 1:missing_num){
    donors[[i]] <- as.matrix(sort(dist.matrix[i,c(c(missing_num + 1):dim(x)[1])])[1 : floor(sqrt(dim(x)[1] - missing_num))])
    ## NEW: If the Mahalanobis distance was zero, the weights were NaN. --> Replace
    ## weights of zero by the smallest observed distance greater than zero:
    donors[[i]][,1][donors[[i]][,1]==0] <- min(donors[[i]][,1][donors[[i]][,1]!=0])
  }
  
  ##Neighbors will be weighted by distance 
  donors.w<-list()		
  for(i in 1:missing_num){
    donors.w[[i]]<-(1/donors[[i]][,1])/sum((1/donors[[i]][,1]))
  }
  
  neighbourindices <- lapply(donors.w, function(x1) sapply(names(x1), function(x2) which(rownames(y)==x2)))
  
  ##Imputation process
  for(j in 1:missing_num){
    sweep(as.matrix(y[neighbourindices[[j]],]), 1, donors.w[[j]], "*")->donors.calculate
    y[j,]<-apply(donors.calculate, MARGIN = 2,sum)
  }
  
  imputed.data<-y
}

# 0-5-2 Function to do imputation a la TOMBI
ImputeWithTOBMI <- function(omicsdata, blockind) {
  
  rownamessafe        <- rownames(omicsdata)
  rownames(omicsdata) <- 1:nrow(omicsdata)
  
  blockscompl <- unique(blockind[apply(omicsdata, 2, function(x) !any(is.na(x)))])
  blocksres   <- setdiff(unique(blockind), blockscompl)
  
  omicsdatacompl <- omicsdata[,blockind==blockscompl]
  
  for(i in seq(along=blocksres))
    omicsdata[,blockind==blocksres[i]] <- ImputeTwo(omicsdatacompl, omicsdata[,blockind==blocksres[i]])
  
  rownames(omicsdata) <- rownamessafe
  
  return(omicsdata)
}

# 0-5-3 HelpFunction for 0-5-1
ImputeTwo <- function(omicsdatacompl, omicsdatamiss) {
  
  reorderind <- order(complete.cases(omicsdatamiss))
  rereorderind <- order(reorderind)
  
  imputed <- TOBMIfast(x = omicsdatacompl[reorderind,], y = omicsdatamiss[reorderind,])
  imputed <- imputed[rereorderind,]
  
  return(imputed)
}

# 0-5-4 Fit a RF on train & get the predicitons on test
get_predicition <- function(train, test) {
  " Get predictions from a RF-Model for 'test', whereby the RF is trained on 'train'..
    Important: > All obs. in 'test' are fully observed!
               > 'train' only consits of features that are availabe in 'test'
                  & contains not a single NA.
                  
      --> Train a RF on 'train' & use this to generate predicitons for 'test' then.
                
    Args:
      - train  (DF): DF that only contains variables that are also availabe for 'test'.
                     Must not contain any NA values! 
      - test   (DF): DF that is completly observed for all its observations
                     
    Return:
      - List with: > 'pred_classes' = predicted class for each observation in 'test'
                   > 'pred_prob_pos_class' = predicted probability for a obs. 
                                             to be in class 1
                   > settings of the RF for 'mtry', 'min_node_size' & 'ntree'
  "
  # [0] Check Inputs
  # 0-1 'train' & 'test' must be dataframes w/o missing values
  assert_data_frame(train, any.missing = FALSE)
  assert_data_frame(test, min.rows = 1, any.missing = FALSE)
  
  # 0-2 'train' must not contain any colnames not avaible in 'test' & vic versa
  if (!all((colnames(train) %in% colnames(test)))) {
    stop("Train-Set has different features than the Test-Set!")
  }
  
  if (!all((colnames(test) %in% colnames(train)))) {
    stop("Test-Set has different features than the Train-Set!")
  }
  
  # [1] Train a RF & create predictions for the test-set
  # 1-1 Train a RF on 'train'
  # --1 Convert the response to a factor
  train[,'ytarget'] <- as.factor(train[,'ytarget'])
  
  # --2 Create a formula to pass to the RF 
  #     (define response & use remaining variables as features)
  formula_all <- as.formula(paste('ytarget', " ~ ."))
  
  # --3 Fit the actual RF (only use standard-settings)
  RF <- rfsrc(formula = formula_all, data = train, samptype = "swr", 
              seed = 12345678, var.used = 'all.trees')
  
  # 1-2 Get Prediciton on the testset from the RF
  predicitons <- predict(RF, test)
  
  # 1-3 Return the predicted classes & the predicted probabilites for class '1'
  #     aswell as the settings of the RF
  return(list('pred_classes'        = predicitons$class,
              'pred_prob_pos_class' = predicitons$predicted[,'1'],
              'RF_ntree'            = RF$ntree,
              'RF_mtry'             = RF$mtry,
              'RF_min_node_size'    = RF$nodesize))
}

# 0-5-5 Evaluate the imputation approach
evaal_imp_approach <- function(path = './Data/Raw/BLCA.Rda', frac_train = 0.75, split_seed = 1312,
                               block_seed_train = 1234, block_seed_test = 1312, train_pattern = 2, 
                               train_pattern_seed = 12, test_pattern = 2) {
  "Evaluate the Imputation-Approach on the data 'path' points to.
   The (block-wise) missing values in the train-set are imputed with TOMBI. Then all the blocks that
   are not available for the test-set are then removed from the train-set. A RF is then fit on the
   reamining train-set & used to predict on the test-set then. Based on the predicitions for the 
   test-set, various metrics are calculated.
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
      > A DF with the settings of the experiment (path to the data, train pattern, ...), 
        hte common blocks between train- & test-set, as well as the settings of the RF 
        (ntree, mtry, ...) and the results of the evaluation (AUC; Brier-Score; Accuracy)
  "
  # [0] Check Inputs
  #     --> All arguments are checked in the functions 'get_train_test()' &
  #         'get_predicition()' that werde loaded from 'Code/01_Create_BWM_Pattern.R'
  
  # [1] Load the data & prepare them for the Imputation-Approach
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
  
  # 1-2 Impute the missing values in the train-set
  sum(is.na(train_test_bwm$Train$data)) # - 9496920 NAs  -------------------------------------------------------- TO DELETE
  dim(train_test_bwm$Train$data)        # 232 x 81875
  
  # 1-2-1 Remove the response variable 'ytarget' from the train-set (temporary for imputation)
  train_ytarget                     <- train_test_bwm$Train$data$ytarget
  train_test_bwm$Train$data$ytarget <- NULL
  
  # 1-2-2 Do the imputation 
  train_test_bwm$Train$data_imputed <- ImputeWithTOBMI(omicsdata = train_test_bwm$Train$data, 
                                                       blockind  = train_test_bwm$Train$block_index) 
  
  sum(is.na(train_test_bwm$Train$data_imputed)) # - 0 NAs  -------------------------------------------------------- TO DELETE
  dim(train_test_bwm$Train$data_imputed)        # 232 x 81874 (as 'y_target' was removed )
  
  # 1-3 Remove the blocks from the imputed train-set that are not available for the test-set
  # 1-3-1 Get the names of the observed blocks in test
  observed_test_blocks <- c()
  for (curr_test_block in train_test_bwm$Test$block_names) {
    
    # --1 Which block-index has 'curr_test_block' 
    curr_test_block_idx <- which(train_test_bwm$Test$block_names == curr_test_block)
    
    # --2 Get the corresponding columns to 'curr_test_block'
    curr_test_block_cols <- which(train_test_bwm$Test$block_index == curr_test_block_idx)
    
    # --3 Check if 'curr_test_block' is fully observed in the test-set, if so, add it to 'observed_test_blocks'
    if (sum(is.na(train_test_bwm$Test$data[,curr_test_block_cols])) == 0) {
      observed_test_blocks <- c(observed_test_blocks, curr_test_block)
    }
  }
  
  # 1-3-2 Only keep the observed blocks in the test-set
  blocks_to_keep_ <- which(train_test_bwm$Test$block_names %in% observed_test_blocks)
  cols_to_keep_   <- which(train_test_bwm$Test$block_index %in% blocks_to_keep_)
  y_target_idx_   <- which(colnames(train_test_bwm$Test$data) == "ytarget")
  train_test_bwm$Test$data <- train_test_bwm$Test$data[,c(cols_to_keep_, y_target_idx_)]
  
  # 1-3-3 Remove all features from the imputed train-set that are not observed in test
  blocks_to_keep <- which(train_test_bwm$Train$block_names %in% observed_test_blocks)
  cols_to_keep   <- which(train_test_bwm$Train$block_index %in% blocks_to_keep)
  train_test_bwm$Train$data_imputed <- train_test_bwm$Train$data_imputed[,cols_to_keep]
  
  # 1-4 Add the original response to the train-set
  train_test_bwm$Train$data_imputed$ytarget <- train_ytarget

  # [2] Fit an RF on the imputed train-set, get predictions for the test-set & get metrics
  # 2-1 Train a RF on the imputed train set & use it to predict on the test-set then!
  preds_test_set <- get_predicition(train = train_test_bwm$Train$data_imputed,
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
                    "common_blocks"      = paste(observed_test_blocks, collapse = ' - '),
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

# [1] Run the experiments                                                    ----
# 1-1 Initalize a empty DF to store the results
IMP_res <- data.frame()

# 1-2 Define a list with the paths to the availabe DFs
df_paths <- paste0("./Data/Raw/", list.files("./Data/Raw/"))

# 1-3 Loop over all the possible settings for the evaluation of the IMP-Approach
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
        curr_res <- tryCatch(evaal_imp_approach(path               = curr_path, 
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
                                          "common_blocks"      = '---',
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
                                          "BrierScore"         = '---')
                             }
        ) 
        
        # Add the curr_repetition to 'curr_res', before adding it to 'IMP_res'
        curr_res$repetition <- curr_repetition
        curr_res$approach   <- 'Imputation'
        
        # Add the results of the setting to 'IMP_res' & save it
        IMP_res <- rbind(IMP_res, curr_res)
        write.csv(IMP_res, './Docs/Evaluation_Results/IMP_Approach/IMP_Eval.csv')
      }
    }
  }
}