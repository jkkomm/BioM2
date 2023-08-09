# BioMLF: Biologically explainable Machine Learning Framework for phenotype prediction using omics data

<!-- Add buttons here -->

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/navendu-pottekkat/awesome-readme?include_prereleases)


<!-- Remove this note if you plan to copy this README -->

> **Note**: At present, the level of refinement in the details of this R package is suboptimal. Nonetheless, a series of novel functionalities are slated for inclusion to augment its capabilities.


## Motivation
Identifying reproducible and interpretable biological patterns from high-dimensional omics data is a critical factor in understanding the risk mechanism of complex disease. As such, explainable machine learning can offer biological insight in addition to personalized risk scoring.

## Deliverables
We have implemented a biologically informed multi-stage machine learning framework termed __BioMLF__ specifically for phenotype prediction using omics-scale data based on prior biological information including gene ontological (GO) and/or KEGG pathways.   

**Features of BioMLF in a nutshell**:   

1. Applicability for multiple omics data modalities (e.g. methylome, transcriptome). 
2. Various biological stratification strategies.    
3. Prioritizing outcome-associated functional patterns.   
4. Personalized scoring based on biological stratified patterns.   
5. Possibility for an extension to learning models of interest.   

## :writing_hand: Authors

Shunjie Zhang  ----  (E-mail: 1838076459@qq.com)

# Tutorial


## Data requirements

BioMLF requires:
- Genome-wide DNA methylation data / Bulk RNA-seq data
- Feature annotation data
- Pathway annotation data  (Gene ontological and KEGG pathways are used in this tutorial)
  
```
                             <<< Data requirements >>>

-----------------------------------------------------------------------------
## [Genome-wide DNA methylation data]  (data.frame)
# First column name must be 'label', and the rest are the features (e.g., CpGs).

label cg21870274 cg09499020 cg16535257 cg00168193
     0     0.0057     0.0002    -0.0313     0.0002
     0    -0.0317    -0.0444    -0.0578    -0.0160
     1    -0.0341    -0.0541    -0.0056    -0.0230
     1     0.0811    -0.0029     0.0049     0.0274
     1    -0.0187     0.0475     0.1168     0.0169
     0    -0.0158     0.0032    -0.0173     0.0133

--------------------------------------------------------------------------------
## [Feature annotation data]  (data.frame)
# The data frame must contain the two column names 'ID' and 'entrezID' .

         ID entrezID symbol
 cg00000029     5934   RBL2
 cg00000109    64778 FNDC3B
 cg00000155   221927  BRAT1
 cg00000221   162282 ANKFN1
 cg00000236     7419  VDAC3
 cg00000289       87  ACTN1

-----------------------------------------------------------------------------
## [Pathway annotation data] (list)
# The name of each subset of the list is the ID of the pathway, and each subset contains a vector of gene entrezIDs.

List of 15719
 $ GO:0000002: chr [1:31] "142" "291" "1763" "1890" ...
 $ GO:0000003: chr [1:1513] "2" "18" "49" "51" ...
 $ GO:0000012: chr [1:12] "142" "1161" "2074" "3981" ...
 $ GO:0000017: chr [1:2] "6523" "6524"
 $ GO:0000018: chr [1:131] "60" "86" "142" "604" ...
 $ GO:0000019: chr [1:7] "2068" "4292" "4361" "7014" ...
 $ GO:0000022: chr [1:9] "4926" "6795" "9055" "9212" ...
 $ GO:0000023: chr [1:3] "2548" "2595" "8972"

-----------------------------------------------------------------------------

```


## Prediction
if you want to predict the phenotype.
```
result=BioMLF ( TrainData = data , TestData = NULL ,       ## If you only have one dataset
         pathlistDB = pathlistDB ,                         ## ==>> [Pathway annotation data]
         FeatureAnno = FeatureAnno ,                       ## ==>> [Feature annotation data]
         classifier = 'liblinear' , nfolds = 5 ,           ## Choose your learner( use "lrns()" ) , currently only cross-validation is supported
         Inner_CV = 'None' , inner_folds=10 ,               ## Whether to use nested resampling 
         cores = 5                                         ## Parallel support
         )
                                                           ...(More Detail)

[1] "-----------------------------------------------------------"
[1] "------------========<<<<  Completed!  >>>>======-----------"
[1] "-----------------------------------------------------------"
[1] "{|>>>=====Learner: liblinear---Performance Metric---==>> AUC:0.953 ACC:0.876 PCCs:0.785 ======<<<|}"
  resampling_id learner_name       AUC       ACC      PCCs
1             1    liblinear 0.9642857 0.9285714 0.8309358
2             2    liblinear 0.9387755 0.7857143 0.7114370
3             3    liblinear 0.9132653 0.8214286 0.7304734
4             4    liblinear 0.9940828 0.9615385 0.9101721
5             5    liblinear 0.9526627 0.8846154 0.7415218
Time difference of 10.99379 mins
[1] "######——————  Well Done！！！——————######"

> str(result)

```
##  Biological interpretability
If you want to explore which biological pathways have a potential impact on the disease/phenotype,
  please set the parameter ( target = 'pathways') .Show the association between each biological pathway used for prediction and the phenotype.
  
```
result=BioMLF ( TrainData = data , TestData = NULL ,              
         pathlistDB = pathlistDB ,                         
         FeatureAnno = FeatureAnno ,                       
         classifier = 'liblinear' , nfolds = 5 ,          
         target='pathways',                           ##==>>  [ target = 'pathways']
         cores = 5                                        
         )

[1] "-----------------------------------------------------------"
[1] "------------========<<<<  Completed!  >>>>======-----------"
[1] "-----------------------------------------------------------"
             id       cor      p.value adjust_p.value
254    hsa05226 0.6821436 6.059108e-20   2.200668e-16
165    hsa04814 0.6819709 6.242139e-20   2.267145e-16
1141 GO:0015698 0.6633595 1.374259e-18   4.991309e-15
692  GO:0006821 0.6517819 8.435684e-18   3.063840e-14
557  GO:0003341 0.6499000 1.124604e-17   4.084562e-14
                                      term
254  Gastric cancer - Homo sapiens (human)
165  Motor proteins - Homo sapiens (human)
1141             inorganic anion transport
692                     chloride transport
557                        cilium movement
Time difference of 8.649476 mins
[1] "######——————  Well Done！！！——————######"

> str(result)

```

## Pathways Module

A pathway matrix can be obtained by using BioMLF(, target = 'pathways'). The WGCNA-based method aggregates pathways with similar expression patterns into a module, and uses biological semantic information to assist in screening modules with high biological interpretability, and compares these biological pathway modules association with phenotype.

### FindPara_Module（）：Using Biological Semantic Information to Assist in Selecting Optimal Parameters
```
result=BioMLF ( TrainData = data , TestData = NULL ,              
         pathlistDB = pathlistDB ,                         
         FeatureAnno = FeatureAnno ,                       
         classifier = 'liblinear' , nfolds = 5 ,          
         target='pathways',                           ##==>>  [ target = 'pathways']
         cores = 5                                        
         )
result$
Para=FindPara_Module(pathways_matrix=,control_label=0,minModuleSize = c(10,20),mergeCutHeight=c(0.2,0.25))

```


## Installation


# Contribute
[(Back to top)](#table-of-contents)

You can use this section to highlight how people can contribute to your project.

You can add information on how they can open issues or how they can sponsor the project.

# License
[(Back to top)](#table-of-contents)

You can also mention what license the project uses. I usually add it like this:

[MIT license](./LICENSE)
