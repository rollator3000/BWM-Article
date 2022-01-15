"Script to do a quality-check on the evaluation results! 

  > Ensure the settings of the evaluations were always the same 
    for the various approaches (same seed, same block-order, ...)
  > Ensure the approaches were evaluated correctly & do not miss results
    in certain settings due to an unexpected behavior
"
# [0] SetWD, load packages, define fix variables and functions                ----
# 0-1 Set WD
setwd("/Users/frederik/Desktop/BWM-Article/")             # Mac
setwd("C:/Users/kuche/Desktop/BWM-Paper")                 # Windows
setwd("/dss/dsshome1/lxc0B/ru68kiq3/Project/BWM-Article") # Server

# 0-2 Load packages

# 0-3 Define variables
# 0-3-1 Define paths to the evaluation-results
path_cc_res  <- "./Docs/Evaluation_Results/CC_Approach/CC_Eval.csv"
path_sb_res  <- "./Docs/Evaluation_Results/SB_Approach/SB_Eval.csv"
path_bw_res  <- "./Docs/Evaluation_Results/BW_Approach/BW_Eval.csv"
path_imp_res <- "./Docs/Evaluation_Results/IMP_Approach/IMP_Eval.csv"

# 0-3-2 Columns that must be available in all eval-results
#       (- settings of the evaluation [seeds, block-order in train & test, ...]
#        - resulting metrics of evaluation [F1, Acc, ROC, BrierScore, ...])
nec_cols <- c('path', 'frac_train', 'split_seed', 'block_seed_train', 
              'block_seed_test', 'block_order_train_for_BWM', 
              'block_order_test_for_BWM', 'train_pattern', 'train_pattern_seed', 
              'test_pattern', 'AUC', 'Accuracy', 'Sensitivity', 'Specificity', 
              'Precision', 'Recall', 'F1', 'BrierScore', 'int_seed', 
              'repetition', 'approach')

# 0-4 Define functions

# [1] Start investigation of the results                                      ----
# 1-1 Load the results of the evaluation
cc_res  <- read.csv(path_cc_res)
sb_res  <- read.csv(path_sb_res)
bw_res  <- read.csv(path_bw_res)
imp_res <- read.csv(path_imp_res)

# 1-2 Ensure every result-DF has necessary columns
# 1-2-1 Ensure all laoded DFs have the necessary colums
col_check <- all(sapply(list(cc_res, sb_res, bw_res, imp_res), 
                        function(x) all(nec_cols %in% colnames(x))))

if (!col_check) {
  stop("At least one of the DFs from 'path_cc_res', 'path_sb_res', 'path_bw_res',
        'imp_res' misses a column from 'nec_cols'")
}

# 1-3 Ensure all DFs have the same settings for their evaluation
# 1-3-1 Ensure the values for a given evaluation-settings are identical in all DFs
# ----1 Get all columns w/ settings for the evaluation
eval_set_cols <- nec_cols[c(1:5, 8:10, 19:20)]

# ----2 Get the settings of each approach as a single vector
settings_cc <- sapply(1:nrow(cc_res), 
                      function(x) paste(cc_res[x,eval_set_cols], collapse = " - "))
settings_bw <- sapply(1:nrow(bw_res), 
                      function(x) paste(bw_res[x,eval_set_cols], collapse = " - "))
settings_sb <- sapply(1:nrow(sb_res), 
                      function(x) paste(sb_res[x,eval_set_cols], collapse = " - "))
settings_imp <- sapply(1:nrow(imp_res), 
                       function(x) paste(imp_res[x,eval_set_cols], collapse = " - "))

# ----3 Ensure the settings are the same
which(settings_cc != settings_bw)
# --> Wrong seeds entered with CC-approach...