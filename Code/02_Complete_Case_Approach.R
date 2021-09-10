"Script to evaluate the Complete-Case approach on data with blockwise missingness

  > All those blocks from the training-data that are not available in the test-data
    are removed
  > Then remove all observations from the (remaining) training data that contain 
    missing values
  > Train a RF on the resulting DF & use it then to create predicitons for the 
    test-set (structure of test-data has to be known before trainng a RF)
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

# 0-3 Define fixed variables
# 0-3-1 Define amount of usable cores (parallel computing)
detectCores()
registerDoParallel(cores = 2)

# 0-4 Load functions from 'code/01_Create_BWM_Pattern"
source("./Code/01_Create_BWM_Pattern.R")

# 0-5 Define functions
# 0-5-1 Get predicitons for the test-set from a RF trained on train-set
get_predicition <- function(train, test, ntree = NULL, mtry = NULL, min_node_size = NULL) {
  " Get predictions from a RF-Model for 'test', whereby the RF is trained on 'train'..
    Important: > All obs. in 'test' are fully observed!
               > 'train' only consits of features that are availabe in 'test'
                  & contains not a single NA.
                  
      --> Train a RF on 'train' & use this to generate predicitons for 'test' then.
                
    Args:
      - train  (DF): DF that only contains variables that are also availabe for 'test'.
                     Must not contain any NA values! 
      - test   (DF): DF that is completly observed for all its observations
      - ntree (int): Amount of trees to be fit on 'train' - 300 is the default!
      - mtry  (int): Amount of split-variables to try, when looking for a split. 
                     - If it is 'NULL' it is set to sqrt(p) [p = ncol of train]
      - mi_no (int): Minimum node size - amount of obs. a node must at least 
                     contain, so the model keeps on trying to split them 
                     - 1 is the default.
                     
    Return:
      - vector with the predicted class for each observation in 'test'!
  "
  # [0] Check Inputs
  # 0-1 'train' & 'test' must be dataframes w/o missing values
  assert_data_frame(train, any.missing = FALSE)
  assert_data_frame(test, min.rows = 1, any.missing = FALSE)
  
  # 0-2 'train' must not contain any colnames not avaible in 'test'
  if (!all((colnames(train) %in% colnames(test)))) {
    stop("Train-Set has different features than the Test-Set!")
  }
  
  # 0-3 'ntree', 'min_node_size' & 'mtry' must be integer > 10 / 5 / 1 if not NULL
  if (!is.null(ntree)) assert_int(ntree, lower = 10)
  if (is.null(ntree))  ntree = 300
  if (!is.null(mtry))  assert_int(mtry, lower = 5)
  if (is.null(mtry))   mtry = ceiling(sqrt(ncol(train)))
  if (!is.null(min_node_size)) assert_int(min_node_size, lower = 1)
  if (is.null(min_node_size))  min_node_size = 1
  
  # [1] Train a RF & create predicitons for the test-set
  # 1-1 Train a RF on 'train'
  # --1 Convert the response to a factor
  train[,'ytarget'] <- as.factor(train[,'ytarget'])
  
  # --2 Create a formula to pass to the RF 
  #     (define response & use remaining variables as features)
  formula_all <- as.formula(paste('ytarget', " ~ ."))
  
  # --3 Fit the actual RF
  RF <- rfsrc(formula = formula_all, data = train, 
              ntree = ntree, mtry = mtry, nodesize = min_node_size, 
              samptype = "swr", seed = 12345678, var.used = 'all.trees')
  
  # 1-2 Get Prediciton on the testset from the RF
  predicitons <- predict(RF, test)
  
  # 1-3 Return the predicted classes & the predicted probabilites for class '1'
  return(list('pred_classes' = predicitons$class,
              'pred_prob_pos_class' = predicitons$predicted[,'1']))
}

# [1] Evaluate the approach                                                  ----
# 1-1 Load & process the data, such that we can use it for evaluation
# --1 Load the data with the induced BWM
train_test_bwm <- get_train_test(path = './Data/Raw/BLCA.Rda',  # Path to the data
                                 frac_train = 0.75,             # Fraction of data used for Training (rest for test)
                                 split_seed = 1312,             # Seed for the split of the data into test- & train
                                 block_seed_train = 1234,       # Seed to shuffle the block-order in train
                                 block_seed_test = 1342,        # Seed to shuffle the block-order in test
                                 train_pattern = 2,             # Pattern to introduce to train
                                 train_pattern_seed = 12,       # Seed for the introduction of the BWM into train
                                 test_pattern = 2)              # Pattern for the test-set

# --2 Extract Train- & Test-Set
test_set  <- train_test_bwm$Test
train_set <- train_test_bwm$Train

# --3 Extract the observed features from the test-set
# --3-1 Get the observed features from the test-set
obs_test_fea <- names(which(colSums(is.na(test_set$data)) <= 0))

# --3-2 Drop all other features from the test-set (only contains observed features)
test_set$data <- test_set$data[,obs_test_fea]

# --4 Prepare the training data 
# --4-1 Remove all blocks from the training-data, that are not availabe for test
#       (--> only contains features that are available for test)
train_set$data <- train_set$data[,obs_test_fea]

# --4-2 Remove all observations with missing values
train_set$data <- train_set$data[complete.cases(train_set$data), ]

# 1-2 Train a RF on the processed train_set & get predicitons for the test-set
preds_test_set <- get_predicition(train = train_set$data, test = test_set$data, 
                                  ntree = NULL, mtry = NULL, min_node_size = NULL)

# 1-3 Analyse the results
table(preds_test_set$pred_classes, test_set$data$ytarget)





