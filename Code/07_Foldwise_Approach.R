"Start with the implementation of the FW-Approach"
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

# 0-3 Define variables

# 0-4 Define functions
# 0-4-1 Load functions from 'code/01_Create_BWM_Pattern"
source("./Code/01_Create_BWM_Pattern.R")

# 0-4-2 Fit an RF to data, return the RF & its oob-AUC
fit_RF_get_oob_AUC <- function(data) {
  "Fit an RF (w/ its standard settings) to 'data' - 'ytarget' is used as response.
   Return the oob-AUC of the fit RF, as well as the RF itself.
   
   Args:
    > data (data.frame): Data with at least two columns & observations. 
                         Must contain the column 'ytarget' & no have missing values.
                         
   Return:
    > List with the two entrances:
      > AUC: AUC-metric that was calculated based on the oob-observations of the RF
      > RF:  The RF itself
  "
  # [0] Check Inputs
  # 0-1 'data' has to be a DF, with at least 2 observations & columns
  assert_data_frame(data, any.missing = F, min.rows = 2, min.cols = 2)
  if (!'ytarget' %in% colnames(data)) stop("'data' must contain 'ytarget' as column")
  
  # [1] Fit RF on the data
  # 1-1 Train a RF on 'train'
  # --1 Create a formula to pass to the RF 
  #     (define response & use remaining variables as features)
  formula_all <- as.formula(paste('ytarget', " ~ ."))
  
  # --2 Fit the actual RF (only use standard-settings)
  RF <- rfsrc(formula = formula_all, data = data, samptype = "swr", 
              seed = 12345678, var.used = 'all.trees')
  
  # [2] Get the AUC based on the oob-observations 
  # 2-1 Get the predicted probabilities for the oob-observations
  pred_prob_oob <- RF$predicted.oob[,'1']
  
  # 2-2 Compare the predicted class-prob. with the true classes & calc the AUC
  #  -> in case RF predict all OOB as 0/ 1 error will arise -> set AUC to 0
  AUC <- tryCatch(expr = pROC::auc(data$ytarget, pred_prob_oob, quiet = T),
                  error = function(c) 0)
  
  # [3] Return a list with the AUC & the RF itself
  return(list("RF" = RF,
              "AUC" = AUC))
}


#### ----  USe 'code/05...' as template for implementation of the FW-Approach

# ----------- 1 Set argumets for evaluation
path               = './Data/Raw/BLCA.Rda'
frac_train         = 0.75
split_seed         = 1312
block_seed_train   = 1234
block_seed_test    = 1312
train_pattern      = 2
train_pattern_seed = 12
test_pattern       = 2

# ----------- 2 Start running the function 'eval_bw_approach()'
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
"--> Get a single Test-Train-Set with the pattern 'train_pattern' & 'test_pattern'
    > A list with the lists 'Train' & 'Test'. Both of the lists contain the entrances:
       - 'data' (w/ BWM according to 'train_pattern' / 'test_pattern')
       - 'block_index' (index of the blocks after the order of the blocks has been shuffled)
       - 'block_names' (names of the blocks after they have been shuffled)
       - 'Train' also has a additional entrance 'fold_index' with the assigned fold for each obs.

"
# shape of the data
dim(train_test_bwm$Train$data)
# Folds
train_test_bwm$Train$fold_index
length(unique(train_test_bwm$Train$fold_index))

# Select a current fold (all obs. are observed in the same features)
curr_fold <- train_test_bwm$Train$data[which(train_test_bwm$Train$fold_index == 1),]

# Remove the unknown variables
# - sum of variables that are unknown (for all obs.)
sum(sapply(curr_fold, function(x) sum(is.na(x)) == nrow(curr_fold)))
curr_fold <- curr_fold[,-which(sapply(curr_fold, function(x) sum(is.na(x)) == nrow(curr_fold)))]


