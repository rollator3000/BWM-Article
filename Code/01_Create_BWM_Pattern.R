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
# 0-4-1 Load the data and process it to a single DF
load_data <- function(data_path) {
  " Load the raw TCGA-Data and process it. 
      - extract the response 'TP53' from the mutation block
      - mutation block is then not relevant any more
      - merge the blocks 'clin', 'cnv', 'rna' & 'mirna' to a single DF
      - add the response 'TP53' as 'ytarget' to the DF
    On this data, block-wise missingness will be induced later, such that the 
    various approaches can be evaluated on it.
    
    Args:
      > data_path (str): Path to data-set in 'data/raw'. If the argument doesn't 
                         contain 'data/raw' it will throw an error.
                         
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
  # 0-1 'data_path' must be string and contain 'data/raw'
  assert_string(data_path, pattern = 'Data/Raw')
  
  # [1] Load the data and check it
  # 1-1 Load the data 
  tmp <- new.env()
  load(file = data_path, envir = tmp)
  
  # 1-2 Check whether it contains six blocks ("clin", "mirna", "mutation", "cnv", "rna")
  if (any(sapply(c("clin", "mirna", "mutation", "cnv", "rna"), 
                 function(x) !x %in% names(tmp)))) {
    stop('"data_path" led to a DF that does not have the six blocks ("clin", "mirna", "mutation", "cnv", "rna")')
  }
  
  # [2] Process the data 
  # 2-1 Extract the target-variable ('TP53' from the 'mutation'-Block)
  ytarget <- tmp$mutation[,"TP53"]
  
  # 2-2 Assign the various multi-omics blocks to 'block'-variables
  #     (except for 'Mutation', as we use a feature from it as response)
  block1 <- tmp$clin
  block2 <- tmp$cnv
  block3 <- tmp$mirna
  block4 <- tmp$rna
  
  # 2-3 Merge the blocks to create a single DF & create a block index
  #     (response-variable is added as 'ytarget') to the data.
  #     (except for 'mutation', as we have the response from it)
  dataset         <- data.frame(cbind(block1, block2, block3, block4))
  dataset$ytarget <- ytarget
  blockind        <- rep(1:4, times = c(ncol(block1), ncol(block2), 
                                        ncol(block3), ncol(block4)))
  
  print(dim(dataset))
  
  # [3] Return 'dataset', 'blockind' & 'block_names' in a list
  return(list('data'        = dataset,
              'block_index' = blockind,
              'block_names' = c('clin', 'cnv', 'mirna', 'rna')))
}

# [1] Do shit                                                                ----
# 1-1 Set Arguments
data_path <- './Data/Raw/BLCA.Rda'
train_pattern <- '2'
test_pattern  <- '4'
seed          <- 1234

# 1-2 Load the data
loool <- load(data_path)
