# BMW-Paper
This is the repository for the article 'Prediction approaches for partly missing multi-omics covariate data: An empirical comparison study'.  
This article is written in colobration with Dr. Hornung & Jonas Hagenberg.  
My personal contribution to the article is the implementation of various RF-approaches and the evaluation of these.


## Project description
This project compares different approaches capable to deal with block-wise missingness in Multi-Omics data.  
The approaches are either random forest based, or based on penalised regression. 

## Block-wise missingness:
- Block-wise missingness is a special type of missingness that appears frequently in the context of Multi-Omics data. 
It can, for example, arise when concatenating multiple clinical studies with the same target variable. Even though the datasets from the different studies have the same target variable, the observed features can still differ! The concatenation of such datasets results then in a DF with block-wise missingness!  

- Regular model fitting on data with block-wise missingness is for most statistical approaches not directly possible, such that either the method needs to be adjusted or the data processed! Besides the training, the test data can also consist of block-wise missingness. Therefore the approaches must be able to deal with block-wise missing data in the test data as well as in the train data. <br>

#### Example for data with blockwise missingness:
Data with blockwise missingness always consists of different **folds** and **blocks**.
  - A **block** describes a set of covariates containing all features collected based on a characteristic.  
    All covariates that are related in content (e.g. *physical properties*: Height & Weight | *educational properties*: Income & Education').  
  - A **fold** represents a set of observations with the same observed blocks.  
    All observations with the same observed features. Each fold is unique and every observation belongs to exactly one of them.
  
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
   
