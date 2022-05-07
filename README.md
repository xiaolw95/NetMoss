# NetMoss   
[NetMoss2](https://github.com/xiaolw95/NetMoss2) now is available!       
Additional functions have been added in this new version:
1. Single file is supported as input.   
2. Network can be constructed automatically from abundance matrix.   
3. A P value is provided for each NetMoss score.   
4. NetMoss results are visualized as paired networks.  


NetMoss is a tool developed for integrating large-scale data and identifying disease associated biomarkers based on network algorithm.    
Here we provide a R package to acheive this goal.     

For more information, please see paper "Large-scale microbiome data integration enables robust biomarker identification" published in Nature Computational Science.     


## Contents  
- [Installation](#installation)     
- [Basic Usage](#basic-usage)     
- [Input](#input)     
- [Output](#output)     
- [Classification](#classification)           

## Installation    
Installation with `devtools`     
```
library(devtools)
install_github("xiaolw95/NetMoss")
library(NetMoss)
```

## Basic Usage     
The NetMoss function is used to calculate NetMoss score of significant bacteria between case and control groups. Users are demanded to provide four directories as follows:      
```
NetMoss(case_dir = case_dir,    
        control_dir = control_dir,    
        net_case_dir = net_case_dir,   
        net_control_dir = net_control_dir)   
```
`case_dir:`  the directory of case datasets.     
`control_dir:`  the directory of control datasets.      
`net_case_dir:`  the directory of case network datasets.      
`net_control_dir:`  the directory of control network datasets.      

We have provided a small dataset to test the function.     
1. Download from the data directory (https://github.com/xiaolw95/NetMoss/tree/main/data) directly. 
Or get the dataset using `git clone` commond in `Linux`:      
```
git clone https://github.com/xiaolw95/NetMoss.git     
cd NetMoss/data
```

2. After getting the dataset, the NetMoss score can be easily calculated using the `NetMoss` function:       
```
##setwd('path-to-data-directory')
case_dir = paste0(getwd(),"/case_dir")
control_dir = paste0(getwd(),"/control_dir")
net_case_dir = paste0(getwd(),"/net_case_dir")
net_control_dir = paste0(getwd(),"/net_control_dir")
result = NetMoss(case_dir = case_dir,    
        control_dir = control_dir,    
        net_case_dir = net_case_dir,   
        net_control_dir = net_control_dir) 
```   

## Input     
Abundance or network matrix should be included in the directory of the input.    

##### Abundance Table
`case_dir` or `control_dir` includes abundance matrix which refers to the relative abundance of case or contol samples, with the row as bacteria and the column as samples. Abundance file can be processed from raw sequence using [QIIME2](https://qiime2.org/), [MetaPhlAn3](https://github.com/biobakery/MetaPhlAn) or other tools.       
| taxon_names   | sample1 | sample2 | sample3 |    
|  ---  |  ---  |  ---  |  ---  |       
|   taxon1    |    60   |    20   |   10    |       
|   taxon2    |    30   |    77   |   89    |    
|   taxon3    |    0    |    23   |   15    |      
|   ... ...   |         |         |         |          

##### Network Matrix
`net_case_dir` or `net_control_dir` includes network matrix which refers to the adjacency matrix of correltaion between the bacteria. Microbial correlation can be deduced from any tools for which [SparCC](https://github.com/bio-developer/sparcc) or [SPIEC-EASI](https://github.com/zdk123/SpiecEasi) are especially recommended.     

|          | taxon1 | taxon2 | taxon3 |      
|  ------  | -----  | -----  | -----  |      
|  taxon1  |    1   |  -0.3  |  0.5   |      
|  taxon2  |  -0.3  |    1   |  0.67  |      
|  taxon3  |   0.5  |  0.67  |    1   |      
|  ... ... |        |        |        |        

## Output
The output of the NetMoss is a table of NetMoss score for each taxon:     
| taxon_names | control_mod | case_mod | NetMoss_score |      
|  ------  | -----  | -----  | -----  |      
|    taxon1   |    1     |      1     |      0.98     |      
|    taxon2   |    2     |      1     |      0.7      |      
|    taxon3   |    3     |      2     |      0.32     |      
|    ... ...  |        |        |        |       

`taxon_names:` the name of the bacteria.      
`control_mod:`  the control module of the bacteria belongs to.      
`case_mod:`  the case module of the bacteria belongs to.     
`NetMoss_score:`  the NetMoss of the bacteria gets.      

## Classification       
In this section, we provide a pipeline to classify case and control groups based on the NetMoss markers. Iterative training and 10-fold cross validation stpes are implemented to guarantee the markers contain network and abundance informations. For this reason, it will take a long time to process the real datasets which contain large samples. Please be patient.
```
netROC(case_dir = case_dir,
      control_dir = control_dir,
      marker = marker,
      metadata = metadata,
      plot.roc = T,
      train.num = 20)
```
`case_dir:` the directory of case datasets.     
`control_dir:` the directory of control datasets.    
`marker:` a table of combined markers identified by NetMoss.     
`metadata:`  a table of clinical informations for all studies.     
`plot.roc:`  a logical parameter. If TRUE then the combined ROC of the result of classification will be plotted.     
`train.num:`  a numerical parameter which refers to trainning times of the model. By default, it is set to 20.        

First of all, efficient markers should be selected manually from the NetMoss result by users. Generally, we recommend a less strict threshold for the sparse network.
Also, a metadata file contains disease or health information for each sample needs to be inculded. The format should be like this:     
|  sample_id |   type  | study |     
|  ------  | -----  | -----  |     
|  SRRXXXXX  | disease | study1 |      
|  SRRXXXXX  | disease | study2 |       
|  SRRXXXXX  | healthy | study1 |        
|  ... ... |        |        |  

After preparing the two files, classification can be realized using the function `netROC`:     
```
marker = data.frame(result[which(result$NetMoss_Score > 0.3),])       
rownames(marker) = marker$taxon_names        
metadata = read.table("metadata.txt",header = T,sep = '\t',row.names = 1)     
myROC = netROC(case_dir,control_dir,marker,metadata)     
```

The result of the classfication is a table includes true positive rate and false positive rate:     
| threhold |  TPR  |  FPR  |      
|  ------  | ----- | ----- |      
|     0    |   1   |   1   |       
|    0.01  |  0.97 | 0.99  |       
|    0.03  |  0.9  | 0.87  |        
|  ... ... |       |       |  

A combined ROC will be ploted if the parameter `plot.roc` is set to be true.     

<img src="https://github.com/xiaolw95/NetMoss/blob/main/NetMoss_ROC.png" width = "500px">     
