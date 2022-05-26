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
library(doParallel)

# 0-3 Define variables

# 0-4 Define functions
# 0-4-1 Load functions from 'code/01_Create_BWM_Pattern"
source("./Code/01_Create_BWM_Pattern.R")

# 0-4-2 Load functions so we can fit and prune a RF 
#       (not possible with 'randomForestSRC'-library...)
source("./code/07_1_simpleRF_adaption.R")

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

#### ----  USe 'code/05...' as template for implementation of the FW-Approach

# ----------- 1 Set arguments for evaluation 
#            (this all will be packed into a single Eval-Function)
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

# TEMPORATILY --- REMOVE EVERY SECOND COLUMN IN TRAIN TO HAVE LOWER COMP. NEEDS
cols_to_keep_ex           <- colnames(train_test_bwm$Train$data)[seq(from = 1, to = 81875, by = 4)]
cols_to_keep_ex           <- c(cols_to_keep_ex, 'ytarget')
train_test_bwm$Train$data <- train_test_bwm$Train$data[cols_to_keep_ex]
train_test_bwm$Test$data  <- train_test_bwm$Test$data[cols_to_keep_ex]

# [2] Check the data, extract a single fold & fit an RF on it
# 2-1 Shape of the train-data & get its availabe folds
dim(train_test_bwm$Train$data)
train_test_bwm$Train$fold_index
length(unique(train_test_bwm$Train$fold_index))

# 2-2 Fit the single trees of the whole RandomForest
Forest <- list()
for (j_ in 1:length(unique(train_test_bwm$Train$fold_index))) {
  
  # Extract the current fold 'j_' and remove all columns w/ NAs 
  curr_fold <- train_test_bwm$Train$data[which(train_test_bwm$Train$fold_index == j_),]
  curr_fold <- curr_fold[,-which(sapply(curr_fold, function(x) sum(is.na(x)) == nrow(curr_fold)))]
  
  # Define formula
  formula_all <- as.formula(paste("ytarget ~ ."))
  
  # Fit a RF to the corresponding folds
  fold_RF <- simpleRF(formula           = formula_all, 
                      data              = curr_fold, 
                      num_trees         = 5,     # Same settings as for 'rfSCR' ! 500 !
                      mtry              = NULL, 
                      min_node_size     = 1,
                      replace           = TRUE,  # always TRUE, as we need OOB!
                      splitrule         = NULL,  # always NULL!
                      unordered_factors = "ignore")
  
  # Grow the single trees of the RF
  fold_RF <- lapply(fold_RF, function(x) {
    x$grow(replace = TRUE)
    x
  })
  
  # Add grown fold_RF to 'Forest'
  Forest[[j_]] <- fold_RF
}

# [3] Get predictions for the test-set




# 


