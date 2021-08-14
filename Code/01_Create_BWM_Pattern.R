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

# [1] Test the implementations                                                ----
# 1-1 Load a raw DF
data_raw <- load_data('./Data/Raw/BLCA.Rda')

# 1-2 Process it
data_processed <- process_loaded_data(data_raw)


#   --> Get the names of the block & corresponding indexes
data_processed$block_names
data_processed$block_index

 #  --> Get the corresponding entrances from the data & the response
data_processed$data[1:5, which(data_processed$block_index == 1)]
data_processed$data$ytarget[1:5]