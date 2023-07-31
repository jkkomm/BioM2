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
##[Genome-wide DNA methylation data]  (data.frame)
# First column name must be 'label', and the rest are the features (e.g., CpGs).

label cg21870274 cg09499020 cg16535257 cg00168193
     0     0.0057     0.0002    -0.0313     0.0002
     0    -0.0317    -0.0444    -0.0578    -0.0160
     1    -0.0341    -0.0541    -0.0056    -0.0230
     1     0.0811    -0.0029     0.0049     0.0274
     1    -0.0187     0.0475     0.1168     0.0169
     0    -0.0158     0.0032    -0.0173     0.0133

--------------------------------------------------------------------------------
##[Feature annotation data]  (data.frame)
#The data frame must contain the two column names 'ID', 'entrezID' .

         ID entrezID symbol
 cg00000029     5934   RBL2
 cg00000109    64778 FNDC3B
 cg00000155   221927  BRAT1
 cg00000221   162282 ANKFN1
 cg00000236     7419  VDAC3
 cg00000289       87  ACTN1

-----------------------------------------------------------------------------


```


## Usage


## Installation


# Contribute
[(Back to top)](#table-of-contents)

You can use this section to highlight how people can contribute to your project.

You can add information on how they can open issues or how they can sponsor the project.

# License
[(Back to top)](#table-of-contents)

You can also mention what license the project uses. I usually add it like this:

[MIT license](./LICENSE)
