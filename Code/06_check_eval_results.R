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

# 0-4 Define functions