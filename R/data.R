#' @docType data
#' @name GO_Ancestor
#'
#' @title Pathways in the GO database and their Ancestor
#' @description  Inclusion relationships between pathways
#' @details In the GO database, each pathway will have its own ancestor pathway.
#' Map pathways in GO database to about 20 common ancestor pathways.
#'
#' @format A data frame :
#' \describe{
#'   ...
#' }
#' @source From GO.db
NULL




#' @docType data
#' @name GO_Ancestor_exact
#'
#' @title Pathways in the GO database and their Ancestor
#' @description  Inclusion relationships between pathways
#' @details In the GO database, each pathway will have its own ancestor pathway.
#' Map pathways in GO database to about 400 common ancestor pathways.
#'
#' @format A data frame :
#' \describe{
#'   ...
#' }
#' @source From GO.db
NULL

#' @docType data
#' @name MethylAnno
#'
#' @title An example about FeatureAnno for methylation data
#' @description  An example about FeatureAnno for methylation data
#' @details The annotation data stored in a data.frame for probe
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.
#'
#' @format A data frame :
#' \describe{
#'   ...
#' }
NULL

#' @docType data
#' @name TransAnno
#'
#' @title An example about FeatureAnno for gene expression
#' @description  An example about FeatureAnno for gene expression
#' @details The annotation data stored in a data.frame for probe
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.
#'
#' @format A data frame :
#' \describe{
#'   ...
#' }
NULL

#' @docType data
#' @name MethylData_Test
#'
#' @title An example about TrainData/TestData for methylation data
#' @description  An example about TrainData/TestData for methylation data
#' MethylData_Test.
#' @details The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @format A data frame :
#' \describe{
#'   ...
#' }
NULL

#' @docType data
#' @name TransData_Test
#'
#' @title An example about TrainData/TestData for gene expression
#' @description  An example about TrainData/TestData for gene expression
#' MethylData_Test.
#' @details The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @format A data frame :
#' \describe{
#'   ...
#' }
NULL

#' @docType data
#' @name GO2ALLEGS_BP
#'
#' @title An example about pathlistDB
#' @description  An example about pathlistDB
#' @details A list of pathways with pathway IDs and their
#' corresponding genes ('entrezID' is used).
#' @format A list :
#' \describe{
#'   ...
#' }
NULL
