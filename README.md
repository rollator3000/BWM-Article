# BMW-Paper
This is the repository for the article 'Prediction approaches for partly missing multi-omics covariate data: An empirical comparison study'.  
This article is written in colobration with Anne-Laure Boulesteix, Roman Hornung *(both IBE @LMU)* & Jonas Hagenberg.  
My personal contribution to the article is the implementation of various RF-approaches and the evaluation of these.  


## Project description
This project compares different approaches capable to deal with block-wise missingness in Multi-Omics data.  
The approaches are either random forest based, or based on penalised regression. 

## Block-wise missingness:
- Block-wise missingness is a special type of missingness that appears frequently in the context of Multi-Omics data. It can, for example, arise when concatenating multiple clinical studies with the same target variable. Even though the datasets from the different studies have the same target variable, the observed features can still differ! The concatenation of such datasets results then in a DF with block-wise missingness!  

- Regular model fitting on data with block-wise missingness is for most approaches not directly possible, such that either the method needs to be adjusted or the data processed! 

- Besides the training, the test data can also consist of block-wise missingness. Therefore the approaches must be able to deal with block-wise missing data in the test data as well as in the train data. <br>

### Example for data with blockwise missingness:
Data with blockwise missingness always consists of different **folds** and **blocks**.  
  - A **block** describes a set of covariates containing all features collected based on a characteristic  
    -> All covariates that are related in content (e.g. *physical properties*: Height & Weight | *educational properties*: Income & Education').  
  - A **fold** represents a set of observations with the same observed blocks.  
    -> All observations with the same observed features. Each fold is unique and every observation belongs to exactly one of them.
  
| ID  | Weight  | Height  | Income  | Education   | g1      | ...   | g100    | Y   |
|---- |-------- |-------- |-------- |-----------  |-------  |-----  |-------  |---  |
| 1   | 65.4    | 187     | 2.536   | Upper       |         |       |         | 1   |
| 2   | 83.9    | 192     | 1.342   | Lower       |         |       |         | 0   |
| 3   | 67.4    | 167     | 5.332   | Upper       |         |       |         | 1   |
| 4   |         |         | 743     | Lower       | -0.42   | ...   | 1.43    | 1   |
| 5   |         |         | 2.125   | Lower       | 0.52    | ...   | -1.37   | 0   |
| 6   | 105.2   | 175     |         |             | -1.53   | ...   | 2.01    | 0   |
| 7   | 71.5    | 173     |         |             | 0.93    | ...   | 0.53    | 0   |
| 8   | 73.0    | 169     |         |             | 0.31    | ...   | -0.07   | 1   |

  - The data consits of three feature-blocks:
     - **Physical properties:**     Weight, Height
     - **Educational properties:**  Income, Education
     - **Biological properties:**   g1, ..., g100
  -  The data consits of three folds:
     - **Fold1:** All observations with observed Physical & Educational properties
     - **Fold2:** All observations with observed Educational & Biological properties
     - **Fold3:** All observations with observed Physical & Biological properties
   

## Data   
* The data comes from the TCGA *(The Cancer Genome Atlas)* and each dataset consits of multiple omics-blocks
* The data was provided by Roman Hornung, who has worked with these multi-omics data-sets already  
* The provided data doesn't contain any missing values, such that the blockwise-missingness needs to be induced  
* Each data-set uses the 'TP53'-Mutation as response and consits of four further (used) blocks 'clinical', 'copy number variation', 'miRNA' & 'RNA'

## Code  
This section contains short descriptions to the scripts in 'Code/' - there is an logical order in these scripts!  
To the '.py' scripts there is always a note to the enviroment to run the script from.  
Details to the enviroments at the end of the READ-ME. 

#### [1] 00_Inspect_raw_data.R
    - Get a overview to the files in 'Data/Raw'
      - do they contain all necessary blocks (clin, mirna, mutation, cnv, rna)
      - is the response 'TP53' part of the 'mutation' block
      - the amount of features & observations
      - amount of variables with missing values
    - Collect the amount of features per block for each DF, get the amount of
      observations, the fraction of reponse-classes that are positive and collect
       it all in a DF (resulting DF saved to Docs/DataInfo)
    - Get the average amount of features, observations and positive classes in the
      target-var over all DFs

#### [2] 01_Create_BWM_Pattern.R
    - All files in data/raw are fully observed & do not contain missing values
    - Define functions to induce the different BWM-Pattern in the raw data 
    - Based on this d

## Folder-Structure  
```
├── README.md <- Top-level README for devs working with this repository
│ 
├── Data <- All the data for this repository
│   │   
│   ├─── raw <- All 13 Multi-Omics DFs (as they were provided & w/o further processing)
│   │         
│   └─── Example_Data <- Examplary data needed for the implementations
│  
├── Docs <- Sources, Results and everything else documenting the repository  
│   │  
│   ├─── Article_Versions <- Different Versions of the article (shall be published in the end)
│   ├─── DataInfo         <- Overview to amount of rows & features per block for each DF in Data/Raw
│   └─── PhenomeImpute    <- Package to use the R-Function 'PhenomeImpute'
│  
├── envs <- the various enviroments of the project
│
└── code <- Code of the repository
    │
    ├── Example <- Templates & Examples
    └── 00_Inspect_raw_data.R
```
## Enviroments
To run the R-Scripts, you need R-Version 4.0.2/ 4.0.3  

Run the PY-Scripts from the correct Conda-Environments!  
Create the enviroment manually or install it via:  
```
conda env create -f .\envs\XYZ.yml
```  