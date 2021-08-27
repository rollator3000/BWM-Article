"Script to define functions to split a file in data/raw into a train- and test-set
 (3:1) and induce the various BWM-Pattern into the train- and test-set.
 
 As the files are extremly large, the BWM needs to be induced on the fly and can
 not be saved in the repository!
"
# [0] SetWD, load packages, define fix variables and fuctions                ----
# 0-1 Set WD
setwd("/Users/frederik/Desktop/BWM-Article/")             # Mac
setwd("C:/Users/kuche/Desktop/BWM-Paper")                 # Windows
setwd("/dss/dsshome1/lxc0B/ru68kiq3/Project/BWM-Article") # Server

# 0-2 Load packages
library(checkmate)

# 0-3 Define fix variables

# 0-4 Define functions
# 0-4-1 Load the multi-omics blocks from 'Data/Raw', store each of them in a list
#       (w/ corresponding name) and return the list!
load_data <- function(path) {
  "Load the files 'path' points to (has to be '.Rda' in 'Data/Raw'), store them 
   in a list (with corresponding names) and return the list then!
  "
  # [0] Check arguments
  # 0-1 'path' has to be a string & contain '.Rda' & 'Data/Raw'
  assert_string(path, pattern = 'Data/Raw')
  assert_string(path, pattern = '.Rda')
  
  # [1] Load the data from 'path' & store the omcis-blocks into a list 
  # 1-1 Create a local environment (-> don't load the variables to the global env)
  local_env <- new.env()
  load(file = path, local_env)
  
  # 1-2 Check whether 'local_env' contains the six necessary omics-blocks 
  #     ("clin", "mirna", "mutation", "cnv", "rna")
  if (any(sapply(c("clin", "mirna", "mutation", "cnv", "rna"), function(x) !x %in% names(local_env)))) {
    stop('"data_path" led to a DF that does not have the six blocks ("clin", "mirna", "mutation", "cnv", "rna")')
  }
  
  # 1-3 Store the omics-blocks into a list and return it
  return(
    list('clin'     = local_env[['clin']],
         'cnv'      = local_env[['cnv']],
         'mirna'    = local_env[['mirna']],
         'rna'      = local_env[['rna']],
         'mutation' = local_env[['mutation']])
  )
}

# 0-4-2 Process the loaded data & merge the various omics-blocks to a single DF
process_loaded_data <- function(raw_data) {
  " Process the list of omics-blocks laoded with 'load_data()' & create a single
    DF from it! The processing includes:
      - extract the response 'TP53' from the mutation block
      - remove the mutation block (not relevant any more then, as we only neeed
                                   the response 'TP53' from this block)
      - merge the remaining blocks 'clin', 'cnv', 'rna' & 'mirna' to a single DF
      - add the response 'TP53' as 'ytarget' to the DF

    Args:
      > raw_data (list): List filled with 6 data-frames (one for each omics-block).
                         Must contain 'clin', 'mirna', 'mutation', 'cnv' & 'rna'!
                         
    Return:
    A list filled with:
      > 'data': A single DF (fully obserbed) made of the four-blocks 'clin', 'cnv', 
        'mirna' & 'rna' (also in this order), as well as the response variable
        'ytarget'.
      > 'block_index: A vector with the index of which variable belongs to which 
         block (e.g. [1, 1, 2, 2, 2, 2] - > first 2-variables 1. block, rest in 2. block)
      > 'block_names': A vector with the names of the blocks in the correct order
  "
  # [0] Check Inputs
  # 0-1 'raw_data' must be a list & contain the relevant blocks
  assert_list(raw_data)
  if (any(sapply(c("clin", "mirna", "mutation", "cnv", "rna"), function(x) !x %in% names(raw_data)))) {
    stop('"raw_data" must contain of 6 blocks ("clin", "mirna", "mutation", "cnv", "rna")')
  }
  
  # 0-2 Each element in 'raw_data' must be a dataframe/ matrix
  if (! all(sapply(raw_data, function(x) class(x) == "data.frame" || class(x) == "matrix"))) {
    stop("'raw_data' must only contain data.frames/ matrices!")
  }
  
  # [1] Process the data 
  # 1-1 Extract the target-variable ('TP53' from the 'mutation'-Block)
  ytarget <- raw_data$mutation[,"TP53"]
  
  # 1-2 Merge the blocks to create a single DF & create a block index. The 
  #     response-variable is added as 'ytarget' to the data. The mutation-block 
  #     is removed as we extracted the response from it!
  dataset         <- data.frame(cbind(raw_data$clin, raw_data$cnv, raw_data$mirna, 
                                      raw_data$rna))
  dataset$ytarget <- ytarget
  blockind        <- rep(1:4, times = c(ncol(raw_data$clin), ncol(raw_data$cnv), 
                                        ncol(raw_data$mirna), ncol(raw_data$rna)))
  
  # [2] Return 'dataset', 'blockind' & 'block_names' in a list
  return(list('data'        = dataset,
              'block_index' = blockind,
              'block_names' = c('clin', 'cnv', 'mirna', 'rna')))
}

# 0-4-3 Split the processed data into a test- & train-set
split_processed_data <- function(data, fraction_train = 0.75, seed = 1312) {
  " Split the processed data (from 'process_loaded_data()') into a train- & 
    test-set. 
    
    Args: 
      > data             (list): List filled with 'data', 'block_index' & 
                                 'block_names' coming from 'process_loaded_data()'
      > fraction_train (double): Fraction of observations that shall be used
                                 for the training-set. 1 - fraction_train equals
                                 the fraction of the test-set!
      > seed              (int): Seed for reproducibility
      
    Return:
      > A list containing two further lists (once train, once test). Each of
        these lists is filled with: 
        > 'data': A single DF (fully obserbed) made of the four-blocks 'clin', 
          'cnv', 'mirna' & 'rna' (also in this order), as well as the response
          variable 'ytarget'.
        > 'block_index: A vector with the index of which variable belongs to  
           which block (e.g. [1, 1, 2, 2, 2, 2] - > first 2-variables 1. block, 
                                                    rest in 2. block)
        > 'block_names': A vector with the names of the blocks in the correct order
  "
  # [0] Check inputs
  # 0-1 'data' has to be list with the entrances 'data', 'block_index' & 'block_names'
  assert_list(data, len = 3)
  if (!all(sapply(names(data), function(x) x %in% c('data', 'block_index', 'block_names')))) {
    stop("'data' must contain 'data', 'block_index' & 'block_names' as entrances")
  }
  
  # 0-2 'data' has to be data.frame, 'block_index' & 'block_names' must be a vector
  assert_data_frame(data$data)
  assert_vector(data$block_index)
  assert_vector(data$block_names)
  
  # 0-3 'block_index'/ 'block_names' must only contain int/ char
  if (!all(sapply(data$block_index, function(x) is.integer(x)))) {
    stop("'data$block_index' must only contain integers")
  }
  if (!all(sapply(data$block_names, function(x) is.character(x)))) {
    stop("'data$block_names' must only contain strings")
  }
  
  # 0-4 'fraction_train' must be float in ]0;1[ & 'seed' an integer
  assert_number(fraction_train, lower = 0, upper = 1)
  assert_int(seed)
  
  # [1] Split the data into test- & train-set 
  #     (incl all corresponding entrances from 'block_index' & 'block_names')
  # 1-1 Get the amount of data-points for the train-set
  amount_train = round(fraction_train * nrow(data$data))
  
  # 1-2 Sample 'amount_train' data-points between 1-amount of observations
  #     & get the row indeces for the test-obs aswell 
  set.seed(seed)
  train_obs <- sample(1:nrow(data$data), amount_train)
  test_obs  <- which(! c(1:nrow(data$data)) %in% train_obs)
  
  # 1-3 Split the list and all its entrances to 'train' & 'test'
  # 1-3-1 TRAIN
  train_list <- list('data' = data$data[train_obs,],
                     'block_index' = data$block_index,
                     'block_names' = data$block_names)
  # 1-3-2 TEST
  test_list <- list('data' = data$data[test_obs,],
                    'block_index' = data$block_index,
                    'block_names' = data$block_names)
  
  
  # [2] Return the Train- & Test-Set in a list with corresponding entrances
  return(list('train_set' = train_list,
              'test_set'  = test_list))
}

# 0-4-4 Shuffle the order of the blocks randomly
shuffle_block_order <- function(data, seed) {
  "Shuffle the block order of 'data' - created in 'split_processed_data()'.
   The order of all blocks is shuffled, except for the 'clin' block, which
   will always be the first block. But instead of shuffeling the data, we only shuffle
   the 'block_index' & 'block_names' (which we need to access the corresponding variables)
  
  Args:
    > data (list): List filled with 'data', 'block_index' & 'block_names' coming
                   from 'split_processed_data()'
    > seed  (int): Seed to make the results reproducible
    
  Return:
    > The original 'data' list, but with changed block-order! For that, only
      entrances 'block_index' & 'block_names' are updated (as these are used
      to access the variables, hence no need to change the order in the data itself!)
  "
  # [0] Check Inputs
  # 0-1 'data' has to be list with the entrances 'data', 'block_index' & 'block_names'
  assert_list(data, len = 3)
  if (!all(sapply(names(data), function(x) x %in% c('data', 'block_index', 'block_names')))) {
    stop("'data' must contain 'data', 'block_index' & 'block_names' as entrances")
  }
  
  # 0-2 'data' has to be data.frame, 'block_index' & 'block_names' must be a vector
  assert_data_frame(data$data)
  assert_vector(data$block_index)
  assert_vector(data$block_names)
  
  # 0-3 'seed' has to be an integer
  assert_int(seed)
  
  # [1] Shuffle the order of the blocks
  # 1-1 Randomly shuffle the order of all blocks except for 'clin'
  #     (stays in first place all the time!)
  set.seed(seed)
  new_order <- c('clin', sample(data$block_names[data$block_names != 'clin']))
  
  # 1-2 Get the new index for the various blocks - according to 'new_order'
  new_idx <- c()
  
  for (curr_block in data$block_names) {
    
    # --1 Get the 'block_index' for 'curr_block' in the original data
    #     (which variables belong to 'curr_block')
    org_idx <- which(data$block_index == which(data$block_names == curr_block))
    
    # --2 Get the position of 'curr_block' in the original 'block_order'
    new_block_position <- which(new_order == curr_block)
    
    # --3 Repeat 'old_position' as often, as curr_block has variables
    new_idx <- c(new_idx, rep(new_block_position, times = length(org_idx)))
  } 
  
  # 1-3 Overwrite the 'block_index' & the 'block_names' in data
  # 1-3-1 Overwrite the 'block_names' in data with 'new_order'
  data$block_names <- new_order
  
  # 1-3-2 Overwrite the 'block_index' in data with 'new_block_order_idx'
  data$block_index <- new_idx
  
  # [2] Return the data-set with shuffled 'block_names' & 'block_index'
  return(data)
}

# 0-4-5 Induce block-wise missingness pattern to Train
induce_bwm_train <- function(data, pattern, seed) {
  "Induce the pattern of block-wise missingness into the train-data - 
   totally there are 5 different patterns (1. pattern, doesn't even have bwm)!
  
   Args: 
    > data   (list): List filled with 'data', 'block_index' & 'block_names' 
                     coming from 'shuffle_block_order()'
    > pattern (int): Which pattern of BWM is induced into the data
                     (must be int in [1-5])
    > seed    (int): Seed to keep results reproducible
  
   Return:
    > data (list) - 'block_index' & 'block_names' are untouched.
      The 'data' is induced with block-wise missingness pattern.
      Additionally add the 'fold_index' to the list - showing
      which observations belong to which fold.
  "
  # [0] Check Inputs
  # 0-1 'data' has to be list with the entrances 'data', 'block_index' & 'block_names'
  assert_list(data, len = 3)
  if (!all(sapply(names(data), function(x) x %in% c('data', 'block_index', 'block_names')))) {
    stop("'data' must contain 'data', 'block_index' & 'block_names' as entrances")
  }
  
  # 0-2 'data' has to be data.frame, 'block_index' & 'block_names' must be a vector
  assert_data_frame(data$data)
  assert_vector(data$block_index)
  assert_vector(data$block_names)
  
  # 0-3 'pattern' has to be a integer in [1;5] & 'seed' an integer
  assert_int(pattern, lower = 1, upper = 5)
  assert_int(seed)
  
  # [1] Induce block-wise missingness
  # 1-1 Pattern 1 - all blocks are observed for all observations
  if (pattern == 1) {
    
    # --1 Add the fold-index to 'data'
    #     (in this pattern all obs. are in fold 1)
    data$fold_index <- rep(1, times = nrow(data$data))
    
    # --2 Return data
    return(data)
  }
  
  # 1-2 Pattern 2 - data is split into 2 folds, whereby one fold is 
  #                 observed in 'clin' + 2 blocks & the other fold is 
  #                 'clin' + 1 observed block
  if (pattern == 2) {
    
    # --1 Mix the original row order (incl. seed for reproducibility)
    set.seed(seed)
    data$data <- data$data[sample(1:nrow(data$data)),]
    
    # --2 Assign the observations to one of the two fold
    # --2-1 Amount of obs. per fold
    obs_per_fold <- round(nrow(data$data) / 2) 
    
    # --2-2 Get fold-index for each observation
    folds <- rep(1:2, times=c(obs_per_fold, 
                              nrow(data$data) - obs_per_fold)) 
    
    # --3 Get the variables for the blocks
    block_1_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[1]))]
    block_2_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[2]))]
    block_3_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[3]))]
    block_4_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[4]))]
    
    # --3 Remove the values of the blocks & folds -according to pattern 2
    data$data[folds==1, block_3_var]                 <- NA
    data$data[folds==2, c(block_2_var, block_4_var)] <- NA
    
    # --4 Add the fold-index to 'data'
    data$fold_index <- folds
    
    # --5 Return the data with induced BWM
    return(data)
  }
  
  # 1-3 Pattern 3 - data is split into 2 folds, whereby each fold is 
  #                 'clin' + 2 blocks
  if (pattern == 3) {
    
    # --1 Mix the original row order (incl. seed for reproducibility)
    set.seed(seed)
    data$data <- data$data[sample(1:nrow(data$data)),]
    
    # --2 Assign the observations to one of the two fold
    # --2-1 Amount of obs. per fold
    obs_per_fold <- round(nrow(data$data) / 2) 
    
    # --2-2 Get fold-index for each observation
    folds <- rep(1:2, times=c(obs_per_fold, 
                              nrow(data$data) - obs_per_fold)) 
    
    # --3 Get the variables for the blocks
    block_1_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[1]))]
    block_2_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[2]))]
    block_3_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[3]))]
    block_4_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[4]))]
    
    # --3 Remove the values of the blocks & folds -according to pattern 3
    data$data[folds==1, block_4_var] <- NA
    data$data[folds==2, block_3_var] <- NA
    
    # --4 Add the fold-index to 'data'
    data$fold_index <- folds
    
    # --5 Return the data with induced BWM
    return(data)
  }
  
  # 1-4 Pattern 4 - data is split into 3 folds, whereby each is observed in 
  #                 'clin' and one additional block 
  if (pattern == 4) {
    
    # --1 Mix the original row order (incl. seed for reproducibility)
    set.seed(seed)
    data$data <- data$data[sample(1:nrow(data$data)),]
    
    # --2 Assign the observations to one of the two fold
    # --2-1 Amount of obs. per fold
    obs_per_fold <- round(nrow(data$data) / 3) 
    
    # --2-2 Get fold-index for each observation
    folds <- rep(1:3, times=c(obs_per_fold,
                              obs_per_fold,
                              nrow(data$data) - 2 * obs_per_fold)) 
    
    # --3 Get the variables for the blocks
    block_1_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[1]))]
    block_2_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[2]))]
    block_3_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[3]))]
    block_4_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[4]))]
    
    # --3 Remove the values of the blocks & folds -according to pattern 4
    data$data[folds==1, c(block_3_var, block_4_var)] <- NA
    data$data[folds==2, c(block_2_var, block_4_var)] <- NA
    data$data[folds==3, c(block_2_var, block_3_var)] <- NA
    
    # --4 Add the fold-index to 'data'
    data$fold_index <- folds
    
    # --5 Return the data with induced BWM
    return(data)
  }
  
  # 1-5 Pattern 5 - data is split into 8 folds, whereby each is observed
  #                 in different combination of blocks ('clin' always!)
  if (pattern == 5) {
    
    # --1 Mix the original row order (incl. seed for reproducibility)
    set.seed(seed)
    data$data <- data$data[sample(1:nrow(data$data)),]
    
    # --2 Assign the observations to one of the two fold
    # --2-1 Amount of obs. per fold
    obs_per_fold <- round(nrow(data$data) / 8) 
    
    # --2-2 Get fold-index for each observation
    folds <- rep(1:8, times=c(obs_per_fold, obs_per_fold, obs_per_fold, obs_per_fold, obs_per_fold, 
                              obs_per_fold, obs_per_fold, nrow(data$data) - 7 * obs_per_fold)) 
    
    # --3 Get the variables for the blocks
    block_1_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[1]))]
    block_2_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[2]))]
    block_3_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[3]))]
    block_4_var <- colnames(data$data)[which(data$block_index == which(data$block_names == data$block_names[4]))]
    
    # --3 Remove the values of the blocks & folds -according to pattern 5
    data$data[folds==2, block_2_var] <- NA
    data$data[folds==3, block_3_var] <- NA
    data$data[folds==4, block_4_var] <- NA
    data$data[folds==5, c(block_3_var, block_4_var)] <- NA
    data$data[folds==6, c(block_2_var, block_4_var)] <- NA    
    data$data[folds==7, c(block_2_var, block_3_var)] <- NA
    data$data[folds==8, c(block_2_var, block_3_var, block_4_var)] <- NA
    
    # --4 Add the fold-index to 'data'
    data$fold_index <- folds
    
    # --5 Return the data with induced BWM
    return(data)
  }
}

# 0-4-6 Induce block-wise missingness patter to Test
induce_bwm_test <- function(data, pattern) {
  "Induce the pattern of block-wise missingness into the test-data - 
   totally there are 4 different patterns (4. pattern, doesn't even have bwm)!
  
   Args: 
    > data   (list): List filled with 'data', 'block_index' & 'block_names' 
                     coming from 'shuffle_block_order()'
    > pattern (int): Which pattern of BWM is induced into the data
                     (must be int in [1-5])
  
   Return:
    > data (list) - 'block_index' & 'block_names' are untouched.
      The 'data' is induced with block-wise missingness pattern.
      Additionally add the 'fold_index' to the list - showing
      which observations belong to which fold.
  "
  # [0] Check Inputs
  # 0-1 'data' has to be list with the entrances 'data', 'block_index' & 'block_names'
  assert_list(data, len = 3)
  if (!all(sapply(names(data), function(x) x %in% c('data', 'block_index', 'block_names')))) {
    stop("'data' must contain 'data', 'block_index' & 'block_names' as entrances")
  }
  
  # 0-2 'data' has to be data.frame, 'block_index' & 'block_names' must be a vector
  assert_data_frame(data$data)
  assert_vector(data$block_index)
  assert_vector(data$block_names)
  
  # 0-3 'pattern' has to be a integer in [1;4]
  assert_int(pattern, lower = 1, upper = 4)
  
  # [1] Insert BWM-Pattern 'pattern'
  # 1-1 Pattern 1 - Only block 1 is observed
  if (pattern == 1) {
    
    # --1 Get the blocks & the corresponding variables that get censored
    na_blocks <- unique(data$block_index)[2:length(unique(data$block_index))]
    na_cols   <- which(data$block_index %in% na_blocks)
    
    # --2 Set the variables of 'na_cols' to NA
    data$data[,na_cols] <- NA
    
    # --3 Return data
    return(data)
  }
  
  # 1-2 Pattern 2 - Only block 1 & 2is observed
  if (pattern == 2) {
    
    # --1 Get the blocks & the corresponding variables that get censored
    na_blocks <- unique(data$block_index)[3:length(unique(data$block_index))]
    na_cols   <- which(data$block_index %in% na_blocks)
    
    # --2 Set the variables of 'na_cols' to NA
    data$data[,na_cols] <- NA
    
    # --3 Return data
    return(data)
  }
  
  # 1-3 Pattern 3 - Only block 1, 2 & 3 is observed
  if (pattern == 3) {
    
    # --1 Get the blocks & the corresponding variables that get censored
    na_blocks <- unique(data$block_index)[4:length(unique(data$block_index))]
    na_cols   <- which(data$block_index %in% na_blocks)
    
    # --2 Set the variables of 'na_cols' to NA
    data$data[,na_cols] <- NA
    
    # --3 Return data
    return(data)
  }
  
  # 1-4 Pattern 4 - all blocks are observed
  if (pattern == 4) {
    
    # --1 Return data
    return(data)
  }
}

# 0-4-7 Wrap-Function that combines all of the above functions!
get_train_test <- function(path, frac_train = 0.75, split_seed = 1312,
                           block_seed = 1234, train_pattern = 1, 
                           train_pattern_seed = 12, test_pattern = 2) {
  "Wrap up the functions from 0-4-1 to 0-4-7.
   Load the data, process it to a single DF, split it to test- & train-set, 
   shuffle the order of the blocks in test- & train-set & induce BWM into them
   according to 'train_pattern' & 'test_pattern'
  
  Args: 
    > path               (str): Path to a dataset - must contain 'Data/Raw'
    > frac_train       (float): Fraction of observations for the train-set (]0;1[)
    > split_seed         (int): Seed for the split of the data to train & test
    > block_seed         (int): Seed for the shuffeling of the block-order
    > train_pattern      (int): Seed for the induction of the pattern for train
                                (obs. are assigned to different folds!)
    > train_pattern_seed (int): Pattern to induce into train (1, 2, 3, 4, 5)
    > test_pattern       (int): Pattern to induce into test (1, 2, 3, 4)
    
  Return:
    > A list with the lists 'Train' & 'Test'. Both of the lists contain the 
      entrances 'data', 'block_index' & 'block_names'. The 'Train'-List also has
      a additional entrance 'fold_index' with the assigned fold for each obs.
  "
  # [0] Check Inputs
  # 0-1 'path' must be a string with 'Data/Raw' in it
  assert_string(path, pattern = 'Data/Raw')
  
  # 0-2 'frac_train' must be a float between 0 & 1
  assert_numeric(frac_train, lower = 0, upper = 1)
  
  # 0-3 'split_seed', 'block_seed' & 'train_pattern_seed' must be integers
  assert_integer(split_seed)
  assert_integer(block_seed)
  assert_integer(train_pattern_seed)
  
  # 0-4 'train_pattern'/ 'test_pattern' has to be a int in [1-5]/ [1-4]
  assert_integer(train_pattern, lower = 1, upper = 5)
  assert_integer(test_pattern, lower = 1, upper = 4)
}

# [1] Test the implementations                                               ----
# 1-1 Load a raw DF
data_raw <- load_data('./Data/Raw/BLCA.Rda')

# 1-2 Process it
data_processed <- process_loaded_data(data_raw)

# 1-3 Split 'data_processed' to Train- & Test-Set
train_test <- split_processed_data(data_processed, fraction_train = 0.75, seed = 1312)

# 1-4 Shuffle the block-order of test & train
train_shuffled <- shuffle_block_order(train_test$train_set, seed = 1312)
test_shuffled  <- shuffle_block_order(train_test$test_set, seed = 1312)

# 1-5 Induce BWM-Pattern to Train
train_bwm <- induce_bwm_train(data = train_shuffled, pattern = 2, seed = 1312)

# 1-6 Induce BWM-Pattern to Test
test_bwm <- induce_bwm_test(data = test_shuffled, pattern = 1)


#  TEST THE RESULTS OF THE FUNCTION -----------------------------------------------------------------------------------------------
# --> Compare the old & new order + check whether we can access the right variables
train_test$train_set$block_names
train_shuffled$block_names

colnames(train_test$train_set$data)[which(train_test$train_set$block_index == which(train_test$train_set$block_names == "clin"))]
colnames(train_shuffled$data)[which(train_shuffled$block_index == which(train_shuffled$block_names == "clin"))]

colnames(train_test$train_set$data)[which(train_test$train_set$block_index == which(train_test$train_set$block_names == "cnv"))]
colnames(train_shuffled$data)[which(train_shuffled$block_index == which(train_shuffled$block_names == "cnv"))]

colnames(train_test$train_set$data)[which(train_test$train_set$block_index == which(train_test$train_set$block_names == "mirna"))]
colnames(train_shuffled$data)[which(train_shuffled$block_index == which(train_shuffled$block_names == "mirna"))]

colnames(train_test$train_set$data)[which(train_test$train_set$block_index == which(train_test$train_set$block_names == "rna"))]
colnames(train_shuffled$data)[which(train_shuffled$block_index == which(train_shuffled$block_names == "rna"))]
# --> Everything seems to be correct!



# Compare 'train_shuffled' & 'train_bwm'
train_shuffled$data[c(1, nrow(train_shuffled$data)), 
                    seq(from = 1, to = ncol(train_shuffled$data), by  = 500)]

train_bwm$data[c(1, nrow(train_bwm$data)), 
               seq(from = 1, to = ncol(train_bwm$data), by  = 500)]
# --> Seems to be correct

# Access the observations that are observed in a certain block
fold_1      <- train_bwm$data[which(train_bwm$fold_index == 1),]
fold_1_cols <- names(which(colSums(is.na(fold_1)) == 0))

fold_2      <- train_bwm$data[which(train_bwm$fold_index == 2),]
fold_2_cols <- names(which(colSums(is.na(fold_2)) == 0))