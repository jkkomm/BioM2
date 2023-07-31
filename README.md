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


## Usage

```
BioMLF ( TrainData = data , TestData = NULL ,              ## If you only have one dataset
         pathlistDB = pathlistDB ,                         ## ==>> [Pathway annotation data]
         FeatureAnno = FeatureAnno ,                       ## ==>> [Feature annotation data]
         classifier = 'liblinear' , nfolds = 5 ,           ## Choose your learner( use "lrns()" ) , currently only cross-validation is supported
         Inner_CV = 'Yes' , inner_folds=10 ,               ## Whether to use nested resampling 
         cores = 5                                         ## Parallel support
         )
                                                           ...(More Detail)


[1] "===================BioMLF=================="
[1] "<<<<<-----Start-----Resampling(CV,folds=5)-No.1----->>>>>"
[1] "Step1: ReadData"
[1] "     |>Total number of pathways==>>3970"
[1] "Step2: FeartureSelection-features"
[1] "     |> Total number of selected features==>>38060"
[1] "Step3: MergeData"
[1] "     |> Total number of selected pathways==>>3969"
[1] "     |> Min features number of pathways==>>15.......Max features number of pathways==>>502"
[1] "Step4: Reconstruction"
[1] "     <<< Reconstruction Done! >>>     "
[1] "Step5: FeartureSelection-pathways"
[1] "     |> Final number of pathways >>>3969......Min correlation of pathways>>>0.256"
[1] "Step6: Predict and Metric"
[1] "######Resampling NO.1~~~~liblinear==>AUC:0.296 ACC:0.231 PCCs:-0.535"
Time difference of 1.552451 mins
[1] "---------------------####################------------------"
[1] "<<<<<-----Start-----Resampling(CV,folds=5)-No.2----->>>>>"
[1] "Step1: ReadData"
[1] "     |>Total number of pathways==>>3970"
.....
.....


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
