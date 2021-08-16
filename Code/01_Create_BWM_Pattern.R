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
data <- train_test$train_set
seed <- 1234
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
      to access the variables, hence no need to change the order in the data itself!
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
  # 1-1-1 Get the blocks we want to shuffle
  blocks_to_shuffle <- data$block_names[data$block_names != 'clin']
  
  # 1-1-2 Randomly shuffle the order of 'blocks_to_shuffle'
  set.seed(seed)
  new_order <- sample(blocks_to_shuffle)
  
  # 1-3 Get the new block_index (according to 'new_order')
  #     (-> so we know which variables in 'data$data' belong to which block)
  # 1-3-1 Start of with the indeces of 'clin', as it is always the first block
  new_block_order_idx <- c(which(data$block_index == which(data$block_names == 'clin')))
  
  # 1-3-2 Add the indices of the remaining blocks according to the order in
  #      'new_order' & add them to 'new_block_order_idx'
  for (curr_block in new_order) {
    
    # --1 Get the index of the variables from 'curr_block'
    curr_block_idx <- which(data$block_index == which(data$block_names == curr_block))
    
    #--2 Add it to 'new_block_order_idx'
    new_block_order_idx <- c(new_block_order_idx, curr_block_idx)
  }
  
  # 1-4 Overwrite the 'block_index' & the 'block_names' in data
  # 1-4-1 Overwrite the 'block_names' in data with 'new_order'
  data$block_names <-  c('clin', new_order)
  
  # 1-4-2 Overwrite the 'block_index' in data with 'new_block_order_idx'
  data$block_index <- new_block_order_idx
  
  # [2] Return the data-set with shuffled 'block_names' & 'block_index'
  return(data)
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

# --> Next Step: Shuffle the block order of the test- & train-set!
#   --> Get the names of the block & corresponding indexes
train_test$train_set$block_names
train_test$train_set$block_index

 #  --> Get the corresponding entrances from the data & the response
train_test$train_set$data[1:5, which(train_test$train_set$block_index == 1)]
train_test$train_set$data$ytarget[1:5]
