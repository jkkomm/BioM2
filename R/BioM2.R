
#' Prediction by Machine Learning
#'
#' @description
#' Prediction by Machine Learning with different learners ( From 'mlr3' )
#' @param trainData The input training dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param testData The input test dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param predMode The prediction mode. Available options are
#' c('probability', 'classification').
#' @param classifier  Learners in mlr3
#' @param paramlist   Learner parameters
#' @param inner_folds k-fold cross validation ( Only supported when testData = NULL )
#'
#' @return The predicted output for the test data.
#' @import mlr3verse
#' @importFrom mlr3 as_task_classif lrn rsmp msr resample as_task_regr
#' @export
#' @author Shunjie Zhang
#' @examples
#'library(mlr3verse)
#'library(caret)
#'library(BioM2)
#'data=MethylData_Test
#'set.seed(1)
#'part=unlist(createDataPartition(data$label,p=0.8))#Split data
#'predict=baseModel(trainData=data[part,1:10],
#'                  testData=data[-part,1:10],
#'                  classifier = 'svm')#Use 10 features to make predictions,Learner uses svm
#'
#'
baseModel=function ( trainData, testData, predMode = "probability",
                     classifier,paramlist=NULL, inner_folds=10){
  predMode=match.arg(predMode)
  if(!is.null(testData)){
    if (colnames(trainData)[1] != "label") {
      stop("The first column of the 'trainData' must be the 'label'!")
    }
    if (colnames(testData)[1] != "label") {
      stop("The first column of the 'testData' must be the 'label'!")
    }

    if( predMode == 'probability'){
      classifier=paste0('classif.',classifier,'')
      trainData[,1]=as.factor(trainData[,1])
      testData[,1]=as.factor(testData[,1])
      trainData=as_task_classif(trainData,target='label')
      testData=as_task_classif(testData,target='label')
      model=lrn(classifier,predict_type = "prob")

      if(!is.null(paramlist)){
        at = auto_tuner(
          tuner = tnr("grid_search", resolution = 10, batch_size = 5),
          learner =  model,
          search_space = paramlist,
          resampling =rsmp("cv", folds =5),
          measure = msr("classif.acc")
        )
        at$train(trainData)
        model$param_set$values = at$tuning_result$learner_param_vals[[1]]
        model$train(trainData)
        predict=model$predict(testData)$prob[,2]
        return(predict)
      }else{
        sink(nullfile())
        model$train(trainData)
        sink()
        predict=model$predict(testData)$prob[,2]
        return(predict)
      }
    }else if( predMode == 'regression'){
      classifier=paste0('regr.',classifier,'')
      trainData[,1]=as.numeric(trainData[,1])
      testData[,1]=as.numeric(testData[,1])
      trainData=as_task_regr(trainData,target='label')
      testData=as_task_regr(testData,target='label')
      model=lrn(classifier)
      if(!is.null(paramlist)){
        at = auto_tuner(
          tuner = tnr("grid_search", resolution = 5, batch_size = 5),
          learner =  model,
          search_space = paramlist,
          resampling = rsmp("cv", folds =5),
          measure = msr("regr.mae")
        )
        at$train(trainData)
        model$param_set$values = at$tuning_result$learner_param_vals[[1]]
        model$train(trainData)
        predict=model$predict(testData)$response
        return(predict)
      }else{
        model$train(trainData)
        predict=model$predict(testData)$response
        return(predict)
      }
    }
  }
  else if(is.null(testData)){
    if (colnames(trainData)[1] != "label") {
      stop("The first column of the 'trainData' must be the 'label'!")
    }
    if( predMode == 'probability'){
      classifier=paste0('classif.',classifier,'')
      trainData[,1]=as.factor(trainData[,1])
      trainData=as_task_classif(trainData,target='label')
      model=lrn(classifier,predict_type = "prob")
      #set.seed(seed)
      sink(nullfile())
      rr=resample(trainData, model, rsmp("cv", folds = inner_folds))$prediction()
      sink()
      re=as.data.frame(as.data.table(rr))[,c(1,5)]
      re=re[order(re$row_ids),][,2]
      return(re)
    }else if(predMode == 'regression'){
      classifier=paste0('regr.',classifier,'')
      trainData=as_task_regr(trainData,target='label')
      model=lrn(classifier)
      #set.seed(seed)
      sink(nullfile())
      rr=resample(trainData, model, rsmp("cv", folds = inner_folds))$prediction()
      sink()
      re=as.data.frame(as.data.table(rr))[,c(1,3)]
      re=re[order(re$row_ids),][,2]
      return(re)

    }
  }
}

#' Stage 1 Fearture Selection
#'
#' @param Stage1_FeartureSelection_Method Feature selection methods. Available options are
#' c(NULL, 'cor', 'wilcox.test', 'cor_rank', 'wilcox.test_rank').
#' @param data The input training dataset. The first column is the label.
#' @param cutoff The cutoff used for feature selection threshold. It can be any value
#' between 0 and 1. Commonly used cutoffs are c(0.5, 0.1, 0.05, 0.01, etc.).
#' @param featureAnno The annotation data stored in a data.frame for probe
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.
#' (For details, please refer to data( data("MethylAnno") )
#' @param pathlistDB_sub A list of pathways with pathway IDs and their
#' corresponding genes ('entrezID' is used).
#' For details, please refer to ( data("GO2ALLEGS_BP") )
#' @param cores The number of cores used for computation.
#' @param verbose Whether to print running process information to the console
#'
#' @return  A list of matrices with pathway IDs as the associated list member
#' names.
#' @import parallel
#' @importFrom stats wilcox.test
#' @export
#' @author Shunjie Zhang
#' @examples
#'
#' library(parallel)
#' data=MethylData_Test
#' feature_pathways=Stage1_FeartureSelection(Stage1_FeartureSelection_Method='cor',
#'                  data=data,cutoff=0,
#'                  featureAnno=MethylAnno,pathlistDB_sub=GO2ALLEGS_BP,cores=1)
#'
Stage1_FeartureSelection=function(Stage1_FeartureSelection_Method='cor',data=NULL,cutoff=NULL,
                                 featureAnno=NULL,pathlistDB_sub=NULL,cores=1,verbose=TRUE){
  if(Sys.info()[1]=="Windows"){
    cores=1
  }
  if(Stage1_FeartureSelection_Method=='cor'){

    if(verbose)print(paste0('      Using <<  correlation  >>',' ,and you choose cutoff:',cutoff))
    Cor=stats::cor(data$label,data)
    Cor=ifelse(Cor>0,Cor,-Cor)
    names(Cor)=colnames(data)
    Cor_names=names(Cor)
    Cor_cutoff=Cor[which(Cor>cutoff)]
    Cor_cutoff_names=names(Cor_cutoff)
    feature_pathways=mclapply(1:length(pathlistDB_sub),function(x){
      id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
      if(length(id)>10){
        id2=id[which(id %in% Cor_cutoff_names)]
        if(length(id2)<11){
          a=Cor[id]
          id2=names(a)[order(a,decreasing = T)[1:11]]
          return(id2)
        }else{
          return(id2)
        }
      }else{
        return(id)
      }
    } ,mc.cores=cores)
    return(feature_pathways)

  }else if(Stage1_FeartureSelection_Method=='wilcox.test'){

    if(verbose)print(paste0('      Using <<  wilcox.test  >>',' ,and you choose cutoff:',cutoff))
    data_0=data[which(data$label==unique(data$label)[1]),]
    data_1=data[which(data$label==unique(data$label)[2]),]
    Cor=unlist(mclapply(1:ncol(data),function(x) wilcox.test(data_0[,x],data_1[,x])$p.value,mc.cores=cores))
    names(Cor)=colnames(data)
    Cor_names=names(Cor)
    Cor_cutoff=Cor[which(Cor<cutoff)]
    Cor_cutoff_names=names(Cor_cutoff)
    feature_pathways=mclapply(1:length(pathlistDB_sub),function(x){
      id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
      if(length(id)>10){
        id2=id[which(id %in% Cor_cutoff_names)]
        if(length(id2)<11){
          a=Cor[id]
          id2=names(a)[order(a,decreasing = T)[1:11]]
          return(id2)
        }else{
          return(id2)
        }
      }else{
        return(id)
      }
    } ,mc.cores=cores)
    return(feature_pathways)

  }else if(Stage1_FeartureSelection_Method=='cor_rank'){

    if(verbose)print(paste0('      Using <<  correlation_rank  >>',' ,and you choose cutoff:',cutoff))
    Cor=stats::cor(data$label,data)
    Cor=ifelse(Cor>0,Cor,-Cor)
    len=length(Cor)*cutoff
    names(Cor)=colnames(data)
    Cor_names=names(Cor)
    Cor_cutoff=order(Cor,decreasing=T)[1:len]
    Cor_cutoff_names=names(Cor_cutoff)
    feature_pathways=mclapply(1:length(pathlistDB_sub),function(x){
      id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
      if(length(id)>10){
        id2=id[which(id %in% Cor_cutoff_names)]
        if(length(id2)<11){
          a=Cor[id]
          id2=names(a)[order(a,decreasing = T)[1:11]]
          return(id2)
        }else{
          return(id2)
        }
      }else{
        return(id)
      }
    } ,mc.cores=cores)
    return(feature_pathways)
  }else if(Stage1_FeartureSelection_Method=='wilcox.test_rank'){

    if(verbose)print(paste0('      Using <<  wilcox.test_rank  >>',' ,and you choose cutoff:',cutoff))
    data_0=data[which(data$label==unique(data$label)[1]),]
    data_1=data[which(data$label==unique(data$label)[2]),]
    Cor=unlist(mclapply(1:ncol(data),function(x) wilcox.test(data_0[,x],data_1[,x])$p.value,mc.cores=cores))
    len=length(Cor)*cutoff
    names(Cor)=colnames(data)
    Cor_names=names(Cor)
    Cor_cutoff=order(Cor)[1:len]
    Cor_cutoff_names=names(Cor_cutoff)
    feature_pathways=mclapply(1:length(pathlistDB_sub),function(x){
      id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
      if(length(id)>10){
        id2=id[which(id %in% Cor_cutoff_names)]
        if(length(id2)<11){
          a=Cor[id]
          id2=names(a)[order(a,decreasing = T)[1:11]]
          return(id2)
        }else{
          return(id2)
        }
      }else{
        return(id)
      }
    } ,mc.cores=cores)
    return(feature_pathways)

  }else{

    feature_pathways=mclapply(1:length(pathlistDB_sub),function(x){
      id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
      return(id)
    } ,mc.cores=cores)
    return(feature_pathways)

  }

}



#' Stage 2 Fearture Selection
#'
#' @param Stage2_FeartureSelection_Method Feature selection methods. Available options are
#' c(NULL, 'cor', 'wilcox.test', 'RemoveHighcor', 'RemoveLinear').
#' @param data The input training dataset. The first column is the label.
#' @param label The label of dataset
#' @param cutoff The cutoff used for feature selection threshold. It can be any value
#' between 0 and 1.
#' @param preMode The prediction mode. Available options are
#' c('probability', 'classification').
#' @param classifier  Learners in mlr3
#' @param cores The number of cores used for computation.
#' @param verbose Whether to print running process information to the console
#'
#' @return Column index of feature
#' @export
#' @import parallel
#' @importFrom stats  wilcox.test
#' @import caret
#' @author Shunjie Zhang
#'
#'
Stage2_FeartureSelection=function(Stage2_FeartureSelection_Method='RemoveHighcor',data=NULL,label=NULL,cutoff=NULL,preMode=NULL,classifier=NULL,verbose=TRUE,cores=1){

  if(Sys.info()[1]=="Windows"){
    cores=1
  }

  if(preMode=='probability' | preMode=='classification'){

    if(Stage2_FeartureSelection_Method=='cor'){
      up=ifelse(classifier=='lda',0.999,100)
      corr=sapply(1:length(data),function(x) stats::cor(data[[x]],label,method='pearson'))

      if(verbose)print(paste0('     |> Final number of pathways >>>',length(which(corr>cutoff &  corr < up)),'......Min correlation of pathways>>>',round(min(corr[which(corr > cutoff & corr < up)]),digits = 3)))
      index=which(corr>cutoff & corr < up )
      return(index)
    }else if(Stage2_FeartureSelection_Method=='wilcox.test'){

      data=do.call(cbind,data)
      data=cbind(label=label,data)
      data=as.data.frame(data)

      data_0=data[which(data$label==unique(data$label)[1]),]
      data_1=data[which(data$label==unique(data$label)[2]),]
      pvalue=unlist(mclapply(2:ncol(data),function(x) wilcox.test(data_0[,x],data_1[,x])$p.value,mc.cores=60))
      if(cutoff < length(pvalue)){
        index=order(pvalue)[1:cutoff]
      }else{
        index=order(pvalue)
      }

      #index=which(pvalue < cutoff)
      if(verbose)print(paste0('     |> Final number of pathways >>>',length(index),'......Max p-value of pathways>>>',round(max(pvalue[index]),digits = 3)))
      return(index)

    }else if(Stage2_FeartureSelection_Method=='RemoveHighcor'){
      if(!is.null(label)){
        corr=sapply(1:length(data),function(x) stats::cor(data[[x]],label,method='pearson'))
        data=do.call(cbind,data)
        corm=stats::cor(data)
        unindex=caret::findCorrelation(corm,cutoff =cutoff)
        index=which(corr > 0)
        index=setdiff(index,unindex)
        if(verbose)print(paste0('     |> Final number of pathways >>> ',length(index),'......Min correlation of pathways>>>',round(min(corr[index]),digits = 3)))
        return(index)
      }else{
        data=do.call(rbind,data)
        label=data[,1]
        corr=sapply(2:ncol(data),function(x) stats::cor(data[,x],label,method='pearson'))
        index=which(corr > 0)
        data=data[,-1]
        corm=stats::cor(data)
        unindex=caret::findCorrelation(corm,cutoff =cutoff)
        index=setdiff(index,unindex)
        index=index+1
        return(index)
      }
    }else if(Stage2_FeartureSelection_Method=='RemoveLinear'){
      if(!is.null(label)){
        corr=sapply(1:length(data),function(x) stats::cor(data[[x]],label,method='pearson'))
        data=do.call(cbind,data)
        unindex=caret::findLinearCombos(data)$remove
        index=which(corr > 0)
        index=setdiff(index,unindex)
        if(verbose)print(paste0('     |> Final number of pathways >>>',length(index),'......Min correlation of pathways>>>',round(min(corr[index]),digits = 3)))
        return(index)
      }else{
        data=do.call(rbind,data)
        label=data[,1]
        corr=sapply(2:ncol(data),function(x) stats::cor(data[,x],label,method='pearson'))
        index=which(corr > 0)
        data=data[,-1]
        unindex=caret::findLinearCombos(data)$remove
        index=setdiff(index,unindex)
        index=index+1
        return(index)
      }

    }else{
      if(!is.null(label)){
        up=ifelse(classifier=='lda',0.99,100)
        corr=sapply(1:length(data),function(x) stats::cor(data[[x]],label,method='pearson'))

        if(verbose)print(paste0('     |> Final number of pathways >>>',length(order(corr,decreasing=T)[which(corr > 0 & corr < up)]),'......Min correlation of pathways>>>',round(min(corr[which(corr > 0 & corr < up)]),digits = 3)))
        index=which(corr > 0 & corr < up )
        return(index)
      }else{
        data=do.call(rbind,data)
        label=data[,1]
        corr=sapply(2:ncol(data),function(x) stats::cor(data[,x],label,method='pearson'))
        index=which(corr > 0)
        index=index+1
        return(index)
      }

    }
  }
  if(preMode=='regression'){

  }
}



#' Add unmapped probe
#'
#' @param train The input training dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param test The input test dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param Unmapped_num The number of unmapped probes.
#' @param Add_FeartureSelection_Method Feature selection methods. Available options are
#' c('cor', 'wilcox.test').
#' @param anno The annotation data stored in a data.frame for probe
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.
#' (For details, please refer to data( data("MethylAnno") )
#' @param len The number of unmapped probes
#' @param cores The number of cores used for computation.
#' @param verbose Whether to print running process information to the console
#'
#' @return Matrix of unmapped probes
#' @export
#' @import parallel
#' @importFrom stats  wilcox.test
#'
AddUnmapped=function(train=NULL,test=NULL,Unmapped_num=NULL,Add_FeartureSelection_Method='wilcox.test',anno=NULL,len=NULL,verbose=TRUE,cores=1){
  requireNamespace("parallel")
  if(Sys.info()[1]=="Windows"){
    cores=1
  }

  if(is.null(Add_FeartureSelection_Method) | Add_FeartureSelection_Method=='wilcox.test'){
    Unmapped_Train=train[,setdiff(colnames(train),anno$ID)]
    Unmapped_Test=test[,setdiff(colnames(train),anno$ID)]
    Unmapped_0=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[1]),]
    Unmapped_1=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[2]),]
    Unmapped_pvalue=unlist(mclapply(2:ncol(Unmapped_Train),function(x) wilcox.test(Unmapped_0[,x],Unmapped_1[,x])$p.value,mc.cores=cores))
    if(is.null(Unmapped_num)){
      Unmapped_num=len-1
    }
    if(length(Unmapped_pvalue)<Unmapped_num){
      Unmapped_num=length(Unmapped_pvalue)
    }
    Unmapped_id=order(Unmapped_pvalue)[1:Unmapped_num]
    Unmapped_id=Unmapped_id+1
    Unmapped_Train=Unmapped_Train[,Unmapped_id]
    Unmapped_Test=Unmapped_Test[,Unmapped_id]
    Unmapped_Data=list('train'= Unmapped_Train,'test'= Unmapped_Test)
    if(verbose)print(paste0('     |> Add Unmapped features==>>',length(Unmapped_id)))
    return(Unmapped_Data)
  }else if(Add_FeartureSelection_Method=='cor'){
    Unmapped_Train=train[,setdiff(colnames(train),anno$ID)]
    Unmapped_Test=test[,setdiff(colnames(train),anno$ID)]
    Unmapped_pvalue=unlist(mclapply(2:ncol(Unmapped_Train),function(x) stats::cor(Unmapped_Train$label,Unmapped_Train[,x]),mc.cores=cores))
    Unmapped_pvalue=ifelse(Unmapped_pvalue>0,Unmapped_pvalue,-Unmapped_pvalue)
    if(is.null(Unmapped_num)){
      Unmapped_num=len-1
    }
    if(length(Unmapped_pvalue)<Unmapped_num){
      Unmapped_num=length(Unmapped_pvalue)
    }
    Unmapped_id=order(Unmapped_pvalue)[1:Unmapped_num]
    Unmapped_id=Unmapped_id+1
    Unmapped_Train=Unmapped_Train[,Unmapped_id]
    Unmapped_Test=Unmapped_Test[,Unmapped_id]
    Unmapped_Data=list('train'= Unmapped_Train,'test'= Unmapped_Test)
    if(verbose)print(paste0('     |> Add Unmapped features==>>',length(Unmapped_id)))
    return(Unmapped_Data)
  }

}



#' Biologically Explainable Machine Learning Framework
#'
#' @param TrainData The input training dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param TestData The input test dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param pathlistDB A list of pathways with pathway IDs and their
#' corresponding genes ('entrezID' is used).
#' For details, please refer to ( data("GO2ALLEGS_BP") )
#' @param FeatureAnno The annotation data stored in a data.frame for probe
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.
#' (For details, please refer to data( data("MethylAnno") )
#' @param resampling Resampling in mlr3verse.
#' @param nfolds k-fold cross validation ( Only supported when TestData = NULL )
#' @param classifier Learners in mlr3
#' @param predMode The prediction mode. Available options are
#' c('probability', 'classification').
#' @param PathwaySizeUp The upper-bound of the number of genes in each
#' biological pathways.
#' @param PathwaySizeDown The lower-bound of the number of genes in each
#' biological pathways.
#' @param MinfeatureNum_pathways The minimal defined pathway size after mapping your
#' own data to pathlistDB(KEGG database/GO database).
#' @param Add_UnMapped Whether to add unmapped probes for prediction
#' @param Unmapped_num The number of unmapped probes
#' @param Add_FeartureSelection_Method Feature selection methods.
#' @param Inner_CV Whether to perform a k-fold verification on the training set.
#' @param inner_folds k-fold verification on the training set.
#' @param Stage1_FeartureSelection_Method Feature selection methods.
#' @param cutoff The cutoff used for feature selection threshold. It can be any value
#' between 0 and 1.
#' @param Stage2_FeartureSelection_Method Feature selection methods.
#' @param cutoff2 The cutoff used for feature selection threshold. It can be any value
#' between 0 and 1.
#' @param classifier2 Learner for stage 2 prediction(if classifier2==NULL,then it is the same as the learner in stage 1.)
#' @param target Is it used to predict or explore potential biological mechanisms?
#' Available options are c('predict', 'pathways').
#' @param p.adjust.method p-value adjustment method.(holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' @param save_pathways_matrix Whether to output the path matrix file
#' @param cores The number of cores used for computation.
#' @param verbose Whether to print running process information to the console
#'
#' @return A list containing prediction results and prediction result evaluation
#' @export
#' @import ROCR
#' @import caret
#' @importFrom utils head data
#' @importFrom stats  wilcox.test p.adjust cor.test
#' @examples
#'
#'
#'
#' library(mlr3verse)
#' library(caret)
#' library(parallel)
#' library(BioM2)
#' data=MethylData_Test
#' set.seed(1)
#' part=unlist(createDataPartition(data$label,p=0.8))
#' Train=data[part,]
#' Test=data[-part,]
#' pathlistDB=GO2ALLEGS_BP
#' FeatureAnno=MethylAnno
#'
#'
#' pred=BioM2(TrainData = Train,TestData = Test,
#'            pathlistDB=pathlistDB,FeatureAnno=FeatureAnno,
#'            classifier='svm',nfolds=5,
#'            PathwaySizeUp=25,PathwaySizeDown=20,MinfeatureNum_pathways=10,
#'            Add_UnMapped='Yes',Unmapped_num=300,
#'            Inner_CV='None',inner_folds=5,
#'            Stage1_FeartureSelection_Method='cor',cutoff=0.3,
#'            Stage2_FeartureSelection_Method='None',
#'            target='predict',cores=1
#' )#(To explore biological mechanisms, set target=‘pathways’)
#'
#'
#'
BioM2=function(TrainData=NULL,TestData=NULL,pathlistDB=NULL,FeatureAnno=NULL,resampling=NULL,nfolds=5,classifier='liblinear', predMode = "probability",
                PathwaySizeUp=200,PathwaySizeDown=20,MinfeatureNum_pathways=10,
                Add_UnMapped=TRUE,Unmapped_num=300,Add_FeartureSelection_Method='wilcox.test',
                Inner_CV=TRUE,inner_folds=10,
                Stage1_FeartureSelection_Method='cor',cutoff=0.3,
                Stage2_FeartureSelection_Method='RemoveHighcor',cutoff2=0.85,classifier2=NULL,
                target='predict',p.adjust.method='fdr',save_pathways_matrix=FALSE,cores=1,verbose=TRUE){
  if(verbose)print('===================BioM2==================')

  if(Sys.info()[1]=="Windows"){
    cores=1
  }

  prediction=list()
  FeatureAnno$ID=gsub('[\\.\\_\\-]','',FeatureAnno$ID)
  colnames(TrainData)=gsub('[\\.\\_\\-]','',colnames(TrainData))

  if(is.null(TestData)){
    list_pathways=list()
    Record=data.frame(resampling_id=1:nfolds,learner_name=1:nfolds,AUC=1:nfolds,ACC=1:nfolds,PCCs=1:nfolds)
    Resampling=caret::createFolds(TrainData$label,k=nfolds)
    if(!is.null(resampling)){
      nfolds=resampling$param_set$values$folds
    }
    T1=Sys.time()
    for(xxx in 1:nfolds){
      if(verbose)print(paste0('<<<<<-----Start-----','Resampling(CV',',folds=',nfolds,')-No.',xxx,'----->>>>>'))

      if(verbose)print('Step1: ReadData')
      t1=Sys.time()
      if(is.null(resampling)){
        trainData=TrainData[unlist(Resampling[-xxx]),]
        testData=TrainData[unlist(Resampling[xxx]),]
      }else{
        trainData=TrainData[resampling$train_set(xxx),]
        testData=TrainData[resampling$test_set(xxx),]
      }
      geneNum_pathways=sapply(1:length(pathlistDB),function(i) length(pathlistDB[[i]]))
      pathlistDB_sub=pathlistDB[which(geneNum_pathways > PathwaySizeDown & geneNum_pathways < PathwaySizeUp )]
      featureAnno=FeatureAnno[FeatureAnno$ID %in% colnames(trainData),]
      if(verbose)print(paste0('     |>Total number of pathways==>>',length(pathlistDB_sub)))



      if(verbose)print('Step2: FeartureSelection-features')
      feature_pathways=Stage1_FeartureSelection(Stage1_FeartureSelection_Method=Stage1_FeartureSelection_Method,data=trainData,cutoff=cutoff,
                                               featureAnno=featureAnno,pathlistDB_sub=pathlistDB_sub,cores=cores,verbose=verbose)

      lens=sapply(1:length(feature_pathways),function(x) length(feature_pathways[[x]]))

      if(verbose)print('Step3: MergeData')
      trainDataList=mclapply(1:length(feature_pathways),function(x) trainData[,feature_pathways[[x]]] ,mc.cores=cores)
      testDataList=mclapply(1:length(feature_pathways),function(x) testData[,feature_pathways[[x]]] ,mc.cores=cores)
      names(trainDataList)=names(pathlistDB_sub)
      names(testDataList)=names(pathlistDB_sub)
      trainDataList=trainDataList[which(lens>MinfeatureNum_pathways)]
      testDataList=testDataList[which(lens>MinfeatureNum_pathways)]
      featureNum_pathways=sapply(1:length(trainDataList),function(i2) length(trainDataList[[i2]]))
      if(verbose)print(paste0('     |> Total number of selected pathways==>>',length(trainDataList)))
      if(verbose)print(paste0('     |> Min features number of pathways==>>',min(featureNum_pathways)-1,'.......','Max features number of pathways==>>',max(featureNum_pathways)-1))



      #(PredictPathways)
      if(target=='pathways'){
        if(verbose)print('Step4: PredictPathways')
        test=mclapply(1:length(testDataList),function(i6) baseModel(trainData =trainDataList[[i6]],testData =testDataList[[i6]],predMode = predMode,classifier = classifier),mc.cores=cores)
        corr=sapply(1:length(testDataList),function(x) stats::cor(test[[x]],testDataList[[x]]$label,method='pearson'))
        newtest=do.call(cbind, test)
        colnames(newtest)=names(testDataList)
        newtest=cbind(label=testDataList[[1]]$label,newtest)
        list_pathways[[xxx]]=newtest
        if(verbose)print(paste0('     |>min correlation of pathways=====>>>',round(min(corr),digits = 3),'......','max correlation of pathways===>>>',round(max(corr),digits = 3)))
        if(verbose)print('     <<< PredictPathways Done! >>>     ')
        t2=Sys.time()
        if(verbose)print(t2-t1)
        if(verbose)print('---------------------####################------------------')
      }else{
        #(Reconstruction )
        if(verbose)print('Step4: Reconstruction')
        if(Inner_CV==TRUE){
          if(verbose)print('     |> Using Inner CV ~ ~ ~')
          train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =NULL,predMode =predMode,classifier = classifier,inner_folds=inner_folds),mc.cores=cores)
          test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode =predMode,classifier = classifier),mc.cores=cores)
        }else{

          train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =trainDataList[[i4]],predMode = predMode ,classifier = classifier),mc.cores=cores)
          test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode = predMode ,classifier = classifier),mc.cores=cores)


        }

        if(verbose)print('     <<< Reconstruction Done! >>>     ')

        #(FeartureSelection-pathways)
        if(verbose)print('Step5: FeartureSelection-pathways')
        index=Stage2_FeartureSelection(Stage2_FeartureSelection_Method=Stage2_FeartureSelection_Method,data=train,
                                      label=trainDataList[[1]]$label,cutoff=cutoff2,preMode='probability',classifier =classifier,verbose=verbose,cores=cores)


        newtrain=do.call(cbind,train[index])
        colnames(newtrain)=names(trainDataList)[index]
        newtrain=cbind(label=trainDataList[[1]]$label,newtrain)

        newtest=do.call(cbind,test[index])
        colnames(newtest)=names(trainDataList)[index]
        newtest=cbind(label=testDataList[[1]]$label,newtest)
        colnames(newtest)=gsub(':','',colnames(newtest))
        colnames(newtrain)=gsub(':','',colnames(newtrain))
        if(Add_UnMapped==TRUE){
          Unmapped_Data=AddUnmapped(train=trainData,test=testData,Add_FeartureSelection_Method=Add_FeartureSelection_Method,
                                    Unmapped_num=Unmapped_num,len=ncol(newtrain),anno=featureAnno,verbose=verbose,cores=cores)
          newtrain=cbind(newtrain,Unmapped_Data$train)
          newtest=cbind(newtest,Unmapped_Data$test)
          if(verbose)print(paste0('     |> Merge PathwayFeature and AddFeature ==>>',ncol(newtrain)))
        }

        #(Predict and Metric)
        if(verbose)print('Step6: Predict and Metric')

        if(is.null(classifier2)){
          classifier2=classifier
        }
        result=baseModel(trainData=newtrain,testData=newtest,predMode ='probability',classifier = classifier2)
        prediction_part=data.frame(sample=rownames(testData),prediction=result)
        prediction[[xxx]]=prediction_part
        names(prediction)[xxx]=paste('Resample No.',xxx)
        testDataY=testDataList[[1]]$label
        pre=ifelse(result>0.5,1,0)
        Record[xxx,1]=xxx
        Record[xxx,2]=classifier2
        Record[xxx,5]=stats::cor(testDataY,result,method='pearson')
        Record[xxx,3]=ROCR::performance(ROCR::prediction(result,testDataY),'auc')@y.values[[1]]
        testDataY=as.factor(testDataY)
        pre=as.factor(pre)
        Record[xxx,4]=confusionMatrix(pre, testDataY)$overall['Accuracy'][[1]]
        if(verbose)print(paste0('######Resampling NO.',xxx,'~~~~',classifier2,'==>','AUC:',round(Record[xxx,3],digits = 3),' ','ACC:',round(Record[xxx,4],digits = 3),' ','PCCs:',round(Record[xxx,5],digits = 3)))
        t2=Sys.time()
        if(verbose)print(t2-t1)
        if(verbose)print('---------------------####################------------------')

      }

      #(Comprehensive Assessment)
      if(xxx==nfolds){
        if(verbose)print('-----------------------------------------------------------')
        if(verbose)print('------------========<<<<  Completed!  >>>>======-----------')
        if(verbose)print('-----------------------------------------------------------')
        if(target=='pathways'){
          ColN=Reduce(intersect,lapply(1:nfolds,function(x) colnames(list_pathways[[x]])))
          list_pathways=lapply(1:nfolds,function(x) list_pathways[[x]][,ColN])

          index=Stage2_FeartureSelection(Stage2_FeartureSelection_Method=Stage2_FeartureSelection_Method,data=list_pathways,
                                        label=NULL,cutoff=cutoff2,preMode='probability',classifier =classifier,cores=cores)
          matrix_pathways=do.call(rbind,list_pathways)
          matrix_pathways=matrix_pathways[,c(1,index)]
          rownames(matrix_pathways)=rownames(TrainData)[unlist(Resampling)]
          if(save_pathways_matrix==TRUE){
            saveRDS(matrix_pathways,'pathways_matrix.rds')
            if(verbose)print('     |||==>>> Save the Pathways-Matrix ')
          }
          pathways_result=data.frame(
            id=colnames(matrix_pathways)[2:ncol(matrix_pathways)],
            cor=sapply(2:ncol(matrix_pathways),function(i) cor.test(matrix_pathways[,1],matrix_pathways[,i])$estimate),
            p.value=sapply(2:ncol(matrix_pathways),function(i) cor.test(matrix_pathways[,1],matrix_pathways[,i])$p.value)
          )
          pathways_result$adjust_p.value=p.adjust(pathways_result$p.value,method=p.adjust.method)
          GO_Ancestor=NA
          data("GO_Ancestor",envir = environment())
          GO_anno=GO_Ancestor[,1:2]
          GO_anno=GO_anno[-which(duplicated(GO_anno)),]
          colnames(GO_anno)=c('id','term')
          pathways_result=merge(pathways_result,GO_anno,by='id')
          pathways_result=pathways_result[order(pathways_result$cor,decreasing = TRUE),]
          if(verbose)print(head(pathways_result))
          final=list('PathwaysMatrix'= matrix_pathways,'PathwaysResult'= pathways_result)
          T2=Sys.time()
          if(verbose)print(T2-T1)
          if(verbose)print('######-------  Well Done!!!-------######')
          return(final)
        }else{
          T2=Sys.time()
          if(verbose)print(paste0('{|>>>=====','Learner: ',classifier,'---Performance Metric---==>>','AUC:',round(mean(Record$AUC),digits = 3),' ','ACC:',round(mean(Record$ACC),digits = 3),' ','PCCs:',round(mean(Record$PCCs),digits = 3),'======<<<|}'))
          if(verbose)print(Record)
          final=list('Prediction'=prediction,'Metric'=Record,'TotalMetric'=c('AUC'=round(mean(Record$AUC),digits = 3),'ACC'=round(mean(Record$ACC),digits = 3),'PCCs'=round(mean(Record$PCCs),digits = 3)))
          if(verbose)print(T2-T1)
          if(verbose)print('######-------  Well Done!!!-------######')
          return(final)
        }
      }
    }
  }else{
    colnames(TestData)=colnames(TrainData)
    Record=data.frame(learner_name=1,AUC=1,ACC=1,PCCs=1)
    if(verbose)print('Step1: ReadData')
    T1=Sys.time()
    trainData=TrainData
    testData=TestData
    geneNum_pathways=sapply(1:length(pathlistDB),function(i) length(pathlistDB[[i]]))
    pathlistDB_sub=pathlistDB[which(geneNum_pathways > PathwaySizeDown & geneNum_pathways < PathwaySizeUp )]
    featureAnno=FeatureAnno[FeatureAnno$ID %in% colnames(trainData),]
    if(verbose)print(paste0('     |>Total number of pathways==>>',length(pathlistDB_sub)))

    if(verbose)print('Step2: FeartureSelection-features')
    feature_pathways=Stage1_FeartureSelection(Stage1_FeartureSelection_Method=Stage1_FeartureSelection_Method,data=trainData,cutoff=cutoff,
                                             featureAnno=featureAnno,pathlistDB_sub=pathlistDB_sub,cores=cores,verbose=verbose)

    lens=sapply(1:length(feature_pathways),function(x) length(feature_pathways[[x]]))

    if(verbose)print('Step3: MergeData')
    trainDataList=mclapply(1:length(feature_pathways),function(x) trainData[,feature_pathways[[x]]] ,mc.cores=cores)
    testDataList=mclapply(1:length(feature_pathways),function(x) testData[,feature_pathways[[x]]] ,mc.cores=cores)
    names(trainDataList)=names(pathlistDB_sub)
    names(testDataList)=names(pathlistDB_sub)
    trainDataList=trainDataList[which(lens>MinfeatureNum_pathways)]
    testDataList=testDataList[which(lens>MinfeatureNum_pathways)]
    featureNum_pathways=sapply(1:length(trainDataList),function(i2) length(trainDataList[[i2]]))
    if(verbose)print(paste0('     |> Total number of selected pathways==>>',length(trainDataList)))
    if(verbose)print(paste0('     |> Min features number of pathways==>>',min(featureNum_pathways)-1,'.......','Max features number of pathways==>>',max(featureNum_pathways)-1))



    #(PredictPathways)
    if(target=='pathways'){
      if(verbose)print('Step4: PredictPathways')
      if(Inner_CV ==TRUE){
        if(verbose)print('     |> Using Inner CV ~ ~ ~')
        train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =NULL,predMode =predMode,classifier = classifier,inner_folds=inner_folds),mc.cores=cores)
        test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode =predMode,classifier = classifier),mc.cores=cores)
      }else{
        train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =trainDataList[[i4]],predMode = predMode ,classifier = classifier),mc.cores=cores)
        test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode = predMode ,classifier = classifier),mc.cores=cores)
      }
      if(verbose)print('Step5: FeartureSelection-pathways')
      index=Stage2_FeartureSelection(Stage2_FeartureSelection_Method=Stage2_FeartureSelection_Method,data=train,
                                    label=trainDataList[[1]]$label,cutoff=cutoff2,preMode='probability',classifier =classifier,verbose=verbose,cores=cores)
      corr=sapply(1:length(testDataList),function(x) stats::cor(test[[x]],testDataList[[x]]$label,method='pearson'))
      newtest=do.call(cbind, test)
      colnames(newtest)=names(testDataList)
      newtest=cbind(label=testDataList[[1]]$label,newtest)
      matrix_pathways=newtest
      rownames(matrix_pathways)=rownames(TestData)
      if(verbose)print(paste0('     |>min correlation of pathways=====>>>',round(min(corr),digits = 3),'......','max correlation of pathways===>>>',round(max(corr),digits = 3)))
      if(verbose)print('     <<< PredictPathways Done! >>>     ')
      if(save_pathways_matrix==TRUE){
        saveRDS(matrix_pathways,'pathways_matrix.rds')
        if(verbose)print('     |||==>>> Save the Pathways-Matrix ')
      }
      pathways_result=data.frame(
        id=colnames(matrix_pathways)[2:ncol(matrix_pathways)],
        cor=sapply(2:ncol(matrix_pathways),function(i) cor.test(matrix_pathways[,1],matrix_pathways[,i])$estimate),
        p.value=sapply(2:ncol(matrix_pathways),function(i) cor.test(matrix_pathways[,1],matrix_pathways[,i])$p.value)
      )
      pathways_result$adjust_p.value=p.adjust(pathways_result$p.value,method=p.adjust.method)
      GO_Ancestor=NA
      data("GO_Ancestor",envir = environment())
      GO_anno=GO_Ancestor[,1:2]
      GO_anno=GO_anno[-which(duplicated(GO_anno)),]
      colnames(GO_anno)=c('id','term')
      pathways_result=merge(pathways_result,GO_anno,by='id')
      pathways_result=pathways_result[order(pathways_result$cor,decreasing = TRUE),]
      if(verbose)print(head(pathways_result))
      final=list('PathwaysMatrix'= matrix_pathways,'PathwaysResult'= pathways_result)
      T2=Sys.time()
      if(verbose)print(T2-T1)
      if(verbose)print('######-------  Well Done!!!-------######')
      return(final)
    }else{
      #(Reconstruction )
      if(verbose)print('Step4: Reconstruction')
      if(Inner_CV ==TRUE){
        if(verbose)print('     |> Using Inner CV ~ ~ ~')
        train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =NULL,predMode =predMode,classifier = classifier,inner_folds=inner_folds),mc.cores=cores)
        test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode =predMode,classifier = classifier),mc.cores=cores)
      }else{
        train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =trainDataList[[i4]],predMode = predMode ,classifier = classifier),mc.cores=cores)
        test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode = predMode ,classifier = classifier),mc.cores=cores)
      }
      if(verbose)print('     <<< Reconstruction Done! >>>     ')

      #(FeartureSelection-pathways)
      if(verbose)print('Step5: FeartureSelection-pathways')
      index=Stage2_FeartureSelection(Stage2_FeartureSelection_Method=Stage2_FeartureSelection_Method,data=train,
                                    label=trainDataList[[1]]$label,cutoff=cutoff2,preMode='probability',classifier =classifier,verbose=verbose,cores=cores)


      newtrain=do.call(cbind,train[index])
      colnames(newtrain)=names(trainDataList)[index]
      newtrain=cbind(label=trainDataList[[1]]$label,newtrain)

      newtest=do.call(cbind,test[index])
      colnames(newtest)=names(trainDataList)[index]
      newtest=cbind(label=testDataList[[1]]$label,newtest)
      colnames(newtest)=gsub(':','',colnames(newtest))
      colnames(newtrain)=gsub(':','',colnames(newtrain))
      if(Add_UnMapped==TRUE){
        Unmapped_Data=AddUnmapped(train=trainData,test=testData,,Add_FeartureSelection_Method=Add_FeartureSelection_Method,
                                  Unmapped_num=Unmapped_num,anno=featureAnno,verbose=verbose,cores=cores)
        newtrain=cbind(newtrain,Unmapped_Data$train)
        newtest=cbind(newtest,Unmapped_Data$test)
        if(verbose)print(paste0('     |> Merge PathwayFeature and AddFeature ==>>',ncol(newtrain)))
      }

      #(Predict and Metric)
      if(verbose)print('Step6: Predict and Metric')
      if(is.null(classifier2)){
        classifier2=classifier
      }
      result=baseModel(trainData=newtrain,testData=newtest,predMode ='probability',classifier = classifier2)
      predict=data.frame(sample=rownames(testData),prediction=result)
      testDataY=testDataList[[1]]$label
      pre=ifelse(result>0.5,1,0)
      Record[1,1]=classifier
      Record[1,4]=stats::cor(testDataY,result,method='pearson')
      Record[1,2]=ROCR::performance(ROCR::prediction(result,testDataY),'auc')@y.values[[1]]
      testDataY=as.factor(testDataY)
      pre=as.factor(pre)
      Record[1,3]=confusionMatrix(pre, testDataY)$overall['Accuracy'][[1]]
      if(verbose)print(paste0('######~~~~',classifier2,'==>','AUC:',round(Record[1,2],digits = 3),' ','ACC:',round(Record[1,3],digits = 3),' ','PCCs:',round(Record[1,4],digits = 3)))
      final=list('Prediction'=predict,'Metric'=Record)
      T2=Sys.time()
      if(verbose)print(T2-T1)
      if(verbose)print('######-------  Well Done!!!-------######')
      return(final)
    }
  }
}



#' Find suitable parameters for partitioning pathways modules
#'
#' @param pathways_matrix A pathway matrix generated by the BioM2( target='pathways') function.
#' @param control_label The label of the control group ( A single number, factor, or character )
#' @param minModuleSize minimum module size for module detection. Detail for WGCNA::blockwiseModules()
#' @param mergeCutHeight dendrogram cut height for module merging. Detail for WGCNA::blockwiseModules()
#' @param minModuleNum Minimum total number of modules detected
#' @param power 	soft-thresholding power for network construction. Detail for WGCNA::blockwiseModules()
#' @param exact   Whether to divide GO pathways more accurately
#'
#' @return A list containing recommended parameters
#' @export
#' @importFrom utils data
#' @importFrom stats sd aggregate
#' @importFrom WGCNA pickSoftThreshold blockwiseModules
#'
#'
#'
#'
#'
FindParaModule=function(pathways_matrix=NULL,control_label=NULL,minModuleSize = seq(10,20,5),mergeCutHeight=seq(0,0.3,0.1),minModuleNum=20,power=NULL,exact=TRUE){
  if('package:WGCNA' %in% search()){
    final=list()
    if(exact==FALSE){
      GO_Ancestor=NA
      data("GO_Ancestor",envir = environment())
      anno=GO_Ancestor
    }else{
      GO_Ancestor_exact=NA
      data("GO_Ancestor_exact",envir = environment())
      anno=GO_Ancestor_exact
    }
    data=as.data.frame(pathways_matrix)
    label=data.frame(ID=rownames(data),label=data$label)
    data=data[data$label==control_label,]
    data=data[,which(colnames(data) %in% anno$GO)]

    if(is.null(power)){
      powers = c(c(1:10), seq(from = 12, to=20, by=1))
      sink(nullfile())
      sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
      sink()
      if(is.na(sft$powerEstimate)){
        stop('Could not find a proper powers , Please give a power by youself .')
      }else{
        power=sft$powerEstimate
        message('Find the proper power!')
      }
    }
    if(length(minModuleSize)==1 & length(mergeCutHeight)==1){
      stop('Length of both minModuleSize and mergeCutHeight must greater than 1')

    }else{
      Num_module=minModuleSize
      result_list=list()
      for(xxx in 1:length(Num_module)){
        cutoff=mergeCutHeight
        n=length(cutoff)
        result=data.frame(mergeCutHeight=1:n,Number_clusters=1:n,Mean_number_pathways=1:n,Mean_Fraction=1:n,Sd_Fraction=1:n,minModuleSize=1:n)
        for(ii in 1:length(cutoff)){
          sink(nullfile())
          net = blockwiseModules(data, power = power,
                                 TOMType = "unsigned", minModuleSize = Num_module[xxx],
                                 reassignThreshold = 0, mergeCutHeight = cutoff[ii],
                                 numericLabels = TRUE, pamRespectsDendro = FALSE,
                                 saveTOMs = F,
                                 saveTOMFileBase = "femaleMouseTOM",
                                 verbose = 3)
          sink()
          cluster=data.frame(ID=names(net$colors),cluster=net$colors)
          cluster$cluster=cluster$cluster+1
          cluster_list=list()
          faction=vector()
          numder_pathways=vector()
          for(i in 1:max(cluster$cluster)){
            new_anno=anno[which(anno$GO %in% cluster[which(cluster$cluster==i),]$ID),]
            sum=length(which(cluster$cluster==i))
            numder_pathways[i]=sum
            term=unique(new_anno$Ancestor)
            prop=sapply(1:length(term),function(x) length(which(new_anno$Ancestor==term[x]))*100/sum)
            a=data.frame(Term=term[which.max(prop)],Fraction=max(prop))
            faction[i]=a$Fraction
            a$Fraction=paste0(round(a$Fraction,2),'%')
            if(length(which(cluster$cluster==i))>5 & length(which(cluster$cluster==i))<300 ){
              cluster_list[[i]]=a
            }else{
              cluster_list[[i]]=NA
            }
            names(cluster_list)[i]=paste0('cluster_NO.',i,'==>> ','There are ',sum,' pathways')
          }
          result[ii,1]=cutoff[ii]
          result[ii,4]=mean(faction[which(!is.na(cluster_list))])
          result[ii,5]=sd(faction[which(!is.na(cluster_list))])
          result[ii,3]=mean(numder_pathways[which(numder_pathways<300)])
          result[ii,2]=length(which(!is.na(cluster_list)))
          result[ii,6]= Num_module[xxx]
        }
        result_list[[xxx]]=result
      }
      result=do.call(rbind,result_list)
      result=result[which(result$Number_clusters > minModuleNum),]
      d=aggregate(result$Mean_Fraction,by=list(minModuleSize=result$minModuleSize),max)
      colnames(d)[2]='Mean_Fraction'
    }

    best_minModuleSize=d[which.max(d$Mean_Fraction),]$minModuleSize
    message(paste0('The best minModuleSize is:',best_minModuleSize))

    size=result[result$minModuleSize==best_minModuleSize,]
    best_mergeCutHeight=size[which.max(size$Mean_Fraction),]$mergeCutHeight
    message(paste0('The best mergeCutHeight is:', best_mergeCutHeight))
    Para=c(power,best_minModuleSize,best_mergeCutHeight)
    names(Para)=c('power','ModuleSize','mergeCutHeight')
    final[[1]]=result
    final[[2]]=Para
    names(final)=c('TotalResult','BestParameter')
    return(final)
    message('Completed!')
  }else{
    message('If you want to use this function, please install and load the WGCNA package')
  }
}





#' Delineate differential pathway modules with high biological interpretability
#'
#' @param pathways_matrix A pathway matrix generated by the BioM2( target='pathways') function.
#' @param control_label The label of the control group ( A single number, factor, or character )
#' @param power soft-thresholding power for network construction. Detail for WGCNA::blockwiseModules()
#' @param minModuleSize minimum module size for module detection. Detail for WGCNA::blockwiseModules()
#' @param mergeCutHeight dendrogram cut height for module merging. Detail for WGCNA::blockwiseModules()
#' @param cutoff Thresholds for Biological Interpretability Difference Modules
#' @param MinNumPathways Minimum number of pathways included in the biologically interpretable difference module
#' @param p.adjust.method p-value adjustment method.(holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' @param exact Whether to divide GO pathways more accurately
#'
#' @return A list containing differential module results that are highly biologically interpretable
#' @export
#' @importFrom WGCNA pickSoftThreshold blockwiseModules moduleEigengenes
#' @importFrom stats  wilcox.test p.adjust
#' @importFrom utils data

PathwaysModule=function(pathways_matrix=NULL,control_label=NULL,power=NULL,minModuleSize=NULL,mergeCutHeight=NULL,
                        cutoff=70,MinNumPathways=5,p.adjust.method='fdr',exact=TRUE){
  if('package:WGCNA' %in% search()){
    if(exact==FALSE){
      GO_Ancestor=NA
      data("GO_Ancestor",envir = environment())
      anno=GO_Ancestor
    }else{
      GO_Ancestor_exact=NA
      data("GO_Ancestor_exact",envir = environment())
      anno=GO_Ancestor_exact
    }
    final=list()
    data=as.data.frame(pathways_matrix)
    label=data.frame(ID=rownames(data),label=data$label)
    data=data[data$label==control_label,]
    data=data[,which(colnames(data) %in% anno$GO)]

    if(is.null(power)){
      powers = c(1:30)
      sink(nullfile())
      sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
      sink()
      if(is.na(sft$powerEstimate)){
        stop('Could not find a proper powers , Please give a power by youself .')
      }else{
        power=sft$powerEstimate
        message('Find the proper power!')
      }
    }
    sink(nullfile())
    net = blockwiseModules(data, power = power,
                           TOMType = "unsigned", minModuleSize = minModuleSize,
                           reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = F,
                           verbose = 3)
    sink()
    cluster=data.frame(ID=names(net$colors),cluster=net$colors)
    cluster_list=list()
    faction=vector()
    for(i in 0:max(cluster$cluster)){
      new_anno=anno[which(anno$GO %in% cluster[which(cluster$cluster==i),]$ID),]
      sum=length(which(cluster$cluster==i))
      term=unique(new_anno$Ancestor)
      prop=sapply(1:length(term),function(x) length(which(new_anno$Ancestor==term[x]))*100/sum)
      a=data.frame(Term=term[which.max(prop)],Fraction=max(prop))
      faction[i+1]=a$Fraction
      a$Fraction=paste0(round(a$Fraction,2),'%')
      if(length(which(cluster$cluster==i))>5 & length(which(cluster$cluster==i))<1000 ){
        cluster_list[[i+1]]=a
      }else{
        cluster_list[[i+1]]=NA
      }
      names(cluster_list)[i+1]=paste0('Module_NO.',i,'==>> ','There are ',sum,' pathways')
    }
    data=as.data.frame(pathways_matrix)
    label=data.frame(ID=rownames(data),label=data$label)
    data=data[,c(1,which(colnames(data) %in% anno$GO))]
    ALL_eigengene=moduleEigengenes(data[,-1],net$colors)$eigengenes
    ALL_eigengene$label=data$label
    ALL_eigengene=ALL_eigengene[,c(ncol(ALL_eigengene),(1:(ncol(ALL_eigengene)-1)))]
    Cor=data.frame(module=rownames(as.data.frame(stats::cor(ALL_eigengene)))[-1],cor=as.data.frame(stats::cor(ALL_eigengene))[-1,1])
    ALL_eigengene_0=ALL_eigengene[ALL_eigengene$label==0,]
    ALL_eigengene_1=ALL_eigengene[ALL_eigengene$label==1,]
    pvalue=unlist(lapply(1:ncol(ALL_eigengene),function(x) wilcox.test(ALL_eigengene_0[,x],ALL_eigengene_1[,x])$p.value))
    n=ncol(ALL_eigengene)-1

    adjust=data.frame(module=colnames(ALL_eigengene),pvalue=pvalue,adjust_pvalue=p.adjust(pvalue,method=p.adjust.method))
    adjust=adjust[-which(adjust$module=='label'),]


    Num_pathways=as.vector(table(net$colors))
    names(Num_pathways)=NULL

    module=data.frame(module=paste0('ME',names(table(net$colors))),
                      Num_pathways=Num_pathways,
                      Fraction=faction)

    result=merge(module,adjust,by='module')
    result=merge(result,Cor,by='module')
    message(paste0('Filter Module those Fraction less than ',cutoff))

    id=which(result$adjust_pvalue < 0.05 & result$Fraction >= cutoff & result$Num_pathways >= MinNumPathways )
    Result=result[id,]
    Result=Result[order(Result$adjust_pvalue),]
    final[[1]]=cluster
    final[[2]]=result
    final[[3]]=Result
    final[[4]]=pathways_matrix
    names(final)=c('ModuleResult','RAW_PathwaysModule','DE_PathwaysModule','Matrix')
    message('Completed!')
    return(final)
  }else{
    message('If you want to use this function, please install and load the WGCNA package')
  }

}



#' Display biological information within each pathway module
#'
#' @param obj Results produced by PathwaysModule()
#' @param ID_Module ID of the diff module
#' @param exact Whether to divide GO pathways more accurately
#'
#' @return List containing biologically specific information within the module
#' @export
#' @importFrom utils data
#'

ShowModule=function(obj=NULL,ID_Module=NULL,exact=TRUE){
  i=ID_Module
  if(exact==FALSE){
    GO_Ancestor=NA
    data("GO_Ancestor",envir = environment())
    anno=GO_Ancestor
  }else{
    GO_Ancestor_exact=NA
    data("GO_Ancestor_exact",envir = environment())
    anno=GO_Ancestor_exact
  }
  cluster=obj$ModuleResult
  final=list()
  for(x in 1:length(i)){
    new_anno=anno[which(anno$GO %in% cluster[which(cluster$cluster==i[x]),]$ID),]
    sum=length(which(cluster$cluster==i[x]))
    term=unique(new_anno$Ancestor)
    prop=sapply(1:length(term),function(x) length(which(new_anno$Ancestor==term[x]))*100/sum)
    a=new_anno[new_anno$Ancestor==term[which.max(prop)],]
    b=setdiff(cluster[which(cluster$cluster==i[x]),]$ID,a$GO)
    if(length(b)>0){
      b2=anno[anno$GO %in% b,]
      if(length(which(duplicated(b2$GO)))==0){
        a=rbind(a,b2)
      }else{
        b2=b2[-which(duplicated(b2$GO)),]
        a=rbind(a,b2)
      }
    }
    final[[x]]=a
    names(final)[x]=paste0('ME',i[x])
  }
  return(final)
}




#' Visualisation of the results of the analysis of the pathway modules
#'
#' @param BioM2_pathways_obj Results produced by BioM2(,target='pathways')
#' @param FindParaModule_obj Results produced by FindParaModule()
#' @param ShowModule_obj Results produced by ShowModule()
#' @param PathwaysModule_obj Results produced by PathwaysModule()
#' @param exact Whether to divide GO pathways more accurately
#' @param type_text_table Whether to display it in a table
#' @param text_table_theme The topic of this table.Detail for ggtexttable()
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
#' Larger values result in more global views of the manifold, while smaller values result in more local data being preserved.
#' In general values should be in the range 2 to 100.
#' @param spread The effective scale of embedded points. In combination with min_dist, this determines how clustered/clumped the embedded points are.
#' @param min_dist The effective minimum distance between embedded points.
#' Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together,
#' while larger values will result on a more even dispersal of points.
#' The value should be set relative to the spread value,
#' which determines the scale at which embedded points will be spread out.
#' @param size Scatter plot point size
#' @param target_weight Weighting factor between data topology and target topology.
#' A value of 0.0 weights entirely on data, a value of 1.0 weights entirely on target.
#' The default of 0.5 balances the weighting equally between data and target.
#' Only applies if y is non-NULL.
#' @param alpha Alpha for ellipse specifying the transparency level of fill color. Use alpha = 0 for no fill color.
#' @param ellipse logical value. If TRUE, draws ellipses around points.
#' @param ellipse.alpha Alpha for ellipse specifying the transparency level of fill color. Use alpha = 0 for no fill color.
#' @param theme  Default:theme_base(base_family = "serif")
#' @param width  image width
#' @param height image height
#' @param save_pdf Whether to save images in PDF format
#' @param volin Can only be used when PathwaysModule_obj exists. ( Violin diagram )
#' @param control_label Can only be used when PathwaysModule_obj exists. ( Control group label )
#' @param module Can only be used when PathwaysModule_obj exists.( PathwaysModule ID )
#' @param cols palette (vector of colour names)
#'
#' @return a ggplot2 object
#' @export
#' @import ggplot2
#' @import htmlwidgets
#' @import jiebaR
#' @import ggsci
#' @import CMplot
#' @import uwot
#' @import webshot
#' @import wordcloud2
#' @import ggpubr
#' @import ggthemes
#' @importFrom utils data
#' @importFrom  stats aggregate quantile
#' @importFrom  ggstatsplot ggbetweenstats
#'
VisMultiModule=function(BioM2_pathways_obj=NULL,FindParaModule_obj=NULL,ShowModule_obj=NULL,PathwaysModule_obj=NULL,exact=TRUE,
                  type_text_table=FALSE,text_table_theme=ttheme('mOrange'),
                  volin=FALSE,control_label=0,module=NULL,cols=NULL,
                  n_neighbors = 8,spread=1,min_dist =2,target_weight = 0.5,
                  size=1.5,alpha=1,ellipse=TRUE,ellipse.alpha=0.2,theme=ggthemes::theme_base(base_family = "serif"),
                  save_pdf=FALSE,width =7, height=7){
  if(is.null(cols)){
    cols = pal_d3("category20",alpha=alpha)(20)
  }
  if(!is.null(BioM2_pathways_obj)){
    if(exact==FALSE){
      GO_Ancestor=NA
      data("GO_Ancestor",envir = environment())
      anno=GO_Ancestor
    }else{
      GO_Ancestor_exact=NA
      data("GO_Ancestor_exact",envir = environment())
      anno=GO_Ancestor_exact
    }
    if(type_text_table){
      Result=BioM2_pathways_obj$PathwaysResult
      Result=Result[1:10,]
      colnames(Result)=c('ID','Correlation','Pvalue','P-adjusted','Description')
      pic=ggtexttable(Result, rows = NULL,theme=text_table_theme)
      if(save_pdf){
        pic
        ggsave('PathwaysResult_Table.pdf',width =width, height =height)
        return(pic)
      }else{
        pic
        return(pic)
      }
    }else{
      pathways=BioM2_pathways_obj$PathwaysResult[,c('id','p.value')]
      colnames(pathways)=c('SNP','pvalue')
      new_anno=anno[anno$GO %in% pathways$SNP,]
      tb=table(new_anno$Ancestor)
      if(length(which(tb>20))>10){
        Names= names(tb[order(tb,decreasing=T)[1:8]])

      }else{
        Names=names(which(table(new_anno$Ancestor)>20))
      }
      new_anno=new_anno[which(new_anno$Ancestor %in% Names),]
      new_anno=new_anno[,c('GO','Ancestor')]
      colnames(new_anno)=c('SNP','Chromosome')
      Anno=merge(new_anno,pathways)
      a=names(which(table(Anno$Chromosome)>250))
      a2=names(which(table(Anno$Chromosome)<250))
      if(length(a)==0){
        Anno=Anno[Anno$Chromosome %in% a2,]
      }else{
        a2=names(which(table(Anno$Chromosome)<250))
        A=lapply(1:length(a),function(x){
          Anno2=Anno[Anno$Chromosome %in% a[x],]
          c1<- quantile(Anno2$pvalue, probs = 0.03)
          Anno3=Anno2[Anno2$pvalue < c1 ,]
          Anno4=Anno2[Anno2$pvalue > c1 ,]
          Anno4=Anno4[sample(1:nrow(Anno4),150,replace = F),]
          o=rbind(Anno3,Anno4)
        })
        A=do.call(rbind,A)
        A2=Anno[Anno$Chromosome %in% a2,]
        Anno=rbind(A,A2)
      }

      Anno$Position=sample(1:100000,nrow(Anno))
      Anno=Anno[,c(1,2,4,3)]
      if(save_pdf){
        pic=CMplot(Anno,plot.type="c",
                   threshold=c(0.001,0.05)/nrow(pathways),threshold.col=c('red','orange'),
                   multracks=FALSE, H=2,axis.cex=2,chr.den.col=NULL,col=cols,
                   r=2.5,lab.cex=1.7,
                   file.output=T,file='pdf',height=height, width=width)
        return(pic)

      }else{
        pic=CMplot(Anno,plot.type="c",
                   threshold=c(0.001,0.05)/nrow(pathways),threshold.col=c('red','orange'),
                   r=2,lab.cex=1.7,
                   multracks=FALSE, chr.den.col=NULL,H=2,axis.cex=1.7,col=cols,file.output=F)
        return(pic)
      }
    }


  }
  if(!is.null(FindParaModule_obj)){
    if(type_text_table){
      Result=FindParaModule_obj$TotalResult
      Result[,3:5]=round(Result[,3:5],2)
      pic=ggtexttable(Result, rows = NULL, theme=text_table_theme)
      if(save_pdf){
        pic
        ggsave('ParameterSelection_Table.pdf',width =width, height =height)
        return(pic)
      }else{
        pic
        return(pic)
      }
    }else{

      result=FindParaModule_obj$TotalResult
      result$minModuleSize=as.character(result$minModuleSize)
      pic=ggpubr::ggline(result,
                         size=size,
                         x = "mergeCutHeight",
                         y = "Mean_Fraction",
                         linetype = "minModuleSize",
                         shape = "minModuleSize",
                         color = "minModuleSize",
                         title = "Parameter Selection",
                         xlab = "mergeCutHeight",
                         ylab = "Mean_Fraction",
                         palette = cols)+theme
      if(save_pdf){
        pic
        ggsave('ParameterSelection.pdf',width =width, height =height)
        return(pic)
      }else{
        pic
        return(pic)
      }

    }

  }
  if(!is.null(ShowModule_obj)){
    if(type_text_table){
      NAME=names(ShowModule_obj)
      output=paste0(NAME,'_Table.pdf')
      p<-"
      Result=ShowModule_obj[[xxx]]
      colnames(Result)=c('GO','Description','Ancestor','AncestorGO')
      pic=ggtexttable(Result, rows = NULL, theme=text_table_theme)
      if(save_pdf){
         pic
         ggsave(output[xxx],width =width, height =height)
         return(pic)
      }else{
        pic
        return(pic)
      }
     "
      y=sapply(1:length(NAME),function(x) gsub('xxx',x,p))
      pic=eval(parse(text =y))
      return(pic)

    }else{
      NAME=names(ShowModule_obj)
      output=paste0(NAME,'_WordCloud.png')
      p<-"
    words=ShowModule_obj[[NAME[xxx]]]$Name
    engine <- worker()
    segment <- segment(words, engine)
    wordfreqs <- freq(segment)

    wordf <- wordfreqs[order(wordfreqs$freq,decreasing = T),]
    rm=c('of','in','by','for','via','process','regulation','lengthening','to')
    wordf=wordf[-which(wordf$char %in% rm),]
    #colors=rep('skyblue', nrow(wordf))
    colors=rep('darkseagreen', nrow(wordf))
    colors[1:5]='darkorange'
    my_graph <-wordcloud2(wordf,shape = 'circle',color = colors)

    if(save_pdf){
      my_graph
      saveWidget(my_graph,'tmp.html',selfcontained = F)
      webshot('tmp.html',output[xxx], delay =6)
      #return(my_graph)
    }else{
      my_graph
      return(my_graph)
    }
    "
      y=sapply(1:length(NAME),function(x) gsub('xxx',x,p))
      pic=eval(parse(text =y))
      return(pic)
    }

  }
  if(!is.null(PathwaysModule_obj)){
    if(type_text_table){
      Result=PathwaysModule_obj$DE_PathwaysModule[,-4]
      Result$Fraction=round(Result$Fraction,2)
      Result$cor=round(Result$cor,2)
      colnames(Result)=c('Modules','Num_Pathways','Fraction','P-adjusted','Correlation')
      pic=ggtexttable(Result, rows = NULL, theme=text_table_theme)
      if(save_pdf){
        pic
        ggsave('PathwaysModule_Table.pdf',width =width, height =height)
        return(pic)
      }else{
        pic
        return(pic)
      }

    }else if(volin){
      data=PathwaysModule_obj$Matrix[,PathwaysModule_obj$ModuleResult$ID]
      data=moduleEigengenes(data,PathwaysModule_obj$ModuleResult$cluster)$eigengenes
      data$label=PathwaysModule_obj$Matrix[,'label']
      data$label=ifelse(data$label==control_label,'Control','Case')
      colnames(data)[which(colnames(data)==paste0('ME',module))]='y'
      label=NA
      pic=ggstatsplot::ggbetweenstats(
        data=data,
        x = label,
        y = y,
        type = "nonparametric",
        p.adjust.method ='fdr'
      )+ labs(
        x = "Phenotype",
        y = "Module EigenPathways",
        #title = 'Distribution of Module Eigengenes across Phenotype',
        title = paste0('Module',module)
      ) +
        theme(
          # This is the new default font in the plot
          text = element_text(family = "serif", size = 8, color = "black"),
          plot.title = element_text(
            family = "serif",
            size = 20,
            face = "bold",
            color = "#2a475e"
          ),
          plot.subtitle = element_text(
            family = "serif",
            size = 15,
            face = "bold",
            color="#1b2838"
          ),
          plot.title.position = "plot", # slightly different from default
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12)
        )+
        theme(
          axis.ticks = element_blank(),
          axis.line = element_line(colour = "grey50"),
          panel.grid = element_line(color = "#b4aea9"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"),
          panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
        )
      if(save_pdf){
        pic
        ggsave(paste0('PathwaysModule_ME',module,'_VolinPlot.pdf'),width =width, height =height)
        return(pic)
      }else{
        pic
        return(pic)
      }
    }else{
      cluster=PathwaysModule_obj$ModuleResult
      cluster$cluster=paste0('ME',cluster$cluster)
      Result=PathwaysModule_obj$DE_PathwaysModule
      if(nrow(Result)>10){
        Result=Result[1:10,]
      }
      meta=cluster[cluster$cluster %in% Result$module,]
      data=as.data.frame(PathwaysModule_obj$Matrix)
      pdata=data[,meta$ID]

      test=as.data.frame(t(pdata))
      test$ID=rownames(test)
      test=merge(test,meta,by='ID')
      rownames(test)=test$ID
      test=test[,-1]
      test$cluster=as.factor(test$cluster)
      test_umap <- uwot::umap(test, n_neighbors = n_neighbors,spread=spread,min_dist = min_dist,
                              y = test$cluster, target_weight = target_weight)
      test_umap <- as.data.frame(test_umap)
      test_umap$Modules=test$cluster
      pic=ggpubr::ggscatter(test_umap,
                            x='V1',
                            y='V2',
                            size = size,
                            title ="Highly Biologically Explainable Differential Module",
                            subtitle="A UMAP visualization",
                            color = "Modules",
                            alpha = alpha,
                            ellipse = ellipse,
                            ellipse.alpha=ellipse.alpha,
                            ellipse.type="norm",
                            palette =cols,
                            xlab = "UMAP_1",
                            ylab = "UMAP_2",
      )+theme
      if(save_pdf){
        pic
        ggsave('PathwaysModule_UMAP.pdf',width =width, height =height)
        return(pic)
      }else{
        pic
        return(pic)
      }
    }
  }
}









#' Visualisation of significant pathway-level features
#'
#' @param BioM2_pathways_obj Results produced by BioM2(,target='pathways')
#' @param pathlistDB A list of pathways with pathway IDs and their corresponding genes ('entrezID' is used).
#' For details, please refer to ( data("GO2ALLEGS_BP") )
#' @param top Number of significant pathway-level features visualised
#' @param p.adjust.method p-value adjustment method.(holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY","fdr","none")
#' @param alpha The alpha transparency, a number in (0,1). Detail for scale_fill_viridis()
#' @param begin The (corrected) hue in (0,1) at which the color map begins. Detail for scale_fill_viridis().
#' @param end The (corrected) hue in (0,1) at which the color map ends. Detail for scale_fill_viridis()
#' @param option 	A character string indicating the color map option to use. Detail for scale_fill_viridis()
#' @param seq Interval of x-coordinate
#'
#' @return a ggplot2 object
#' @export
#' @import ggplot2
#' @import viridis
#' @importFrom stats p.adjust

PlotPathFearture=function(BioM2_pathways_obj=NULL,pathlistDB=NULL,top=10,p.adjust.method='none',begin=0.1,end=0.9,alpha=0.9,option='C',seq=1){

  data=BioM2_pathways_obj$PathwaysResult
  geneNum_pathways=sapply(1:length(pathlistDB),function(i) length(pathlistDB[[i]]))
  pathlistDB=pathlistDB[which(geneNum_pathways > 20 & geneNum_pathways < 200 )]
  data_top=data[1:top,]
  if(p.adjust.method=='none'){
    data_top$val=-log10(data_top$p.value)
  }else{
    data_top$val=-log10(p.adjust(data_top$p.value,method = p.adjust.method))
  }
  data_top=data_top[order(data_top$val),]
  data_top$size=sapply(1:nrow(data_top),function(x) length(pathlistDB[[data_top$id[x]]]))
  data_top$id=factor(data_top$id,levels = data_top$id)
  max=max(data_top$val)+1
  val=NA
  id=NA
  size=NA
  term=NA
  pic=ggplot(data_top) +
    geom_col(aes(val, id,fill=size), width = 0.6)+
    scale_fill_viridis(alpha=alpha,begin=begin,end=end,direction = -1,
                       name='size',option = option)+
    labs(
      x = '-log(p) ',
      y = NULL,
      title = paste0('Top ',top,' Pathway-Level Features'),
      shape='-log(P-value)'
    )+
    scale_x_continuous(
      limits = c(0, max),
      breaks = seq(0, max, by = seq),
      expand = c(0, 0),
      position = "top"
    ) +
    scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
      axis.ticks.length = unit(0, "mm"),
      #axis.title = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      axis.text.y = element_text(family = "serif", size = 12,face = 'bold',colour='black'),
      axis.text.x = element_text(family = "serif", size = 13,colour='black'),
      axis.title = element_text(size = 18,family = "serif",face = 'bold.italic',vjust = 5),
      plot.title = element_text(size = 20,family = "serif",face = 'bold'),
      legend.text = element_text(family = "serif",face = 'bold'),
      legend.title = element_text(size=13,family = "serif",face = 'bold'),
    )+
    geom_text(
      data = data_top,
      aes(0, y = id, label = term),
      hjust = 0,
      nudge_x = 0.5,
      colour = "white",
      family = "serif",
      size = 6
    )
  return(pic)
}




#' Visualisation Original features that make up the pathway
#'
#' @param data The input omics data
#' @param pathlistDB A list of pathways with pathway IDs and their corresponding genes ('entrezID' is used).
#' For details, please refer to ( data("GO2ALLEGS_BP") )
#' @param FeatureAnno The annotation data stored in a data.frame for probe mapping.
#' It must have at least two columns named 'ID' and 'entrezID'.
#' (For details, please refer to data( data("MethylAnno") )
#' @param PathNames A vector.A vector containing the names of pathways
#' @param p.adjust.method p-value adjustment method.(holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY","fdr","none")
#' @param save_pdf Whether to save images in PDF format
#' @param alpha The alpha transparency, a number in (0,1).
#' @param cols palette (vector of colour names)
#'
#' @return a plot object
#' @export
#' @import ggplot2
#' @import CMplot
#' @import ggsci
#' @importFrom stats p.adjust

PlotPathInner=function(data=NULL,pathlistDB=NULL,FeatureAnno=NULL,PathNames=NULL,
                       p.adjust.method='none',save_pdf=FALSE,alpha=1,cols=NULL){
  Result=list()
  if(is.null(cols)){
    cols = pal_d3("category20",alpha=alpha)(20)
  }
  featureAnno=FeatureAnno[FeatureAnno$ID %in% colnames(data),]
  for(i in 1:length(PathNames)){
    cpg_id=featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB[[PathNames[i]]])]
    cpg_data=data[,c('label',cpg_id)]
    cpg_0=cpg_data[which(cpg_data$label==unique(cpg_data$label)[1]),]
    cpg_1=cpg_data[which(cpg_data$label==unique(cpg_data$label)[2]),]
    pvalue=unlist(lapply(2:ncol(cpg_data),function(x) wilcox.test(cpg_0[,x],cpg_1[,x])$p.value))
    if(p.adjust.method=='none'){
      result=data.frame(SNP=cpg_id,Chromosome=rep(PathNames[i],length(cpg_id)),
                        Position=sample(1:100000,length(cpg_id)),pvalue=pvalue)
    }else{
      result=data.frame(SNP=cpg_id,Chromosome=rep(PathNames[i],length(cpg_id)),
                        Position=sample(1:100000,length(cpg_id)),pvalue=p.adjust(pvalue,method = p.adjust.method))
    }
    Result[[i]]=result
  }
  Result=do.call(rbind,Result)
  pic=CMplot(Result,plot.type="c",
             threshold=c(0.001,0.05),threshold.col=c('red','orange'),
             multracks=FALSE, H=2,axis.cex=2,chr.den.col=NULL,col=cols,
             r=2.5,lab.cex=2,outward = TRUE,signal.cex=2, signal.pch = 18,
             file.output=save_pdf,file='pdf',height=13, width=13)
  return(pic)
}




#' Correlalogram for Biological Differences Modules
#'
#' @param PathwaysModule_obj Results produced by PathwaysModule()
#' @param alpha The alpha transparency, a number in (0,1). Detail for scale_fill_viridis()
#' @param begin The (corrected) hue in (0,1) at which the color map begins. Detail for scale_fill_viridis().
#' @param end The (corrected) hue in (0,1) at which the color map ends. Detail for scale_fill_viridis()
#' @param option 	A character string indicating the color map option to use. Detail for scale_fill_viridis()
#' @param family calligraphic style
#' @return a ggplot object
#' @export
#' @import ggplot2
#' @importFrom  ggstatsplot ggcorrmat
#' @import viridis
#'
#'
#'

PlotCorModule=function(PathwaysModule_obj=NULL,
                       alpha=0.7,begin=0.2,end=0.9,option="C",family="serif"){
  colors=PathwaysModule_obj$ModuleResult$cluster
  names(colors)=PathwaysModule_obj$ModuleResult$ID
  ALL_eigengene=moduleEigengenes(PathwaysModule_obj$Matrix[,names(colors)],colors)$eigengenes
  data=ALL_eigengene[,PathwaysModule_obj$DE_PathwaysModule$module]

  pic=ggcorrmat(
    data  = data,
    type = "nonparametric",
    ggcorrplot.args = list(show.legend =T,pch.cex=10),
    title    = "Correlalogram for Biological Differences Modules",
    subtitle = " ",
    caption = " "
  )+theme(
    # This is the new default font in the plot
    text = element_text(family = family, size = 8, color = "black"),
    plot.title = element_text(
      family = family,
      size = 20,
      face = "bold",
      color = "black"
    ),
    plot.subtitle = element_text(
      family = family,
      size = 15,
      face = "bold",
      color="#1b2838"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 15)
  )+scale_fill_viridis(alpha=alpha,begin=begin,end=end,direction = -1,
                       name='correlation',option = option)+
    theme(legend.title =element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10,color = "black"))
  pic$labels$caption=NULL
  return(pic)
}






#' Network diagram of pathways-level features
#'
#' @param data The input omics data
#' @param pathlistDB A list of pathways with pathway IDs and their corresponding genes ('entrezID' is used).
#' For details, please refer to ( data("GO2ALLEGS_BP") )
#' @param FeatureAnno The annotation data stored in a data.frame for probe mapping.
#' It must have at least two columns named 'ID' and 'entrezID'.
#' (For details, please refer to data( data("MethylAnno") )
#' @param PathNames A vector.A vector containing the names of pathways
#' @param cutoff Threshold for correlation between features within a pathway
#' @param num    The first few internal features of each pathway that are most relevant to the phenotype
#' @param BioM2_pathways_obj Results produced by BioM2()
#'
#' @return a ggplot object
#' @export
#' @import ggplot2
#' @import ggnetwork
#' @import igraph
#' @import ggsci
#' @import ggforce
#'

PlotPathNet=function(data=NULL,BioM2_pathways_obj=NULL,FeatureAnno=NULL,pathlistDB=NULL,PathNames=NULL,
                     cutoff=0.2,num=20){
  featureAnno=FeatureAnno[FeatureAnno$ID %in% colnames(data),]
  sub=list()
  i=1
  for(i in 1:length(PathNames)){
    cpg_id=featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB[[PathNames[i]]])]
    cpg_data=data[,c('label',cpg_id)]
    COR=stats::cor(cpg_data$label,cpg_data[,-1])
    COR=ifelse(COR>0,COR,-COR)
    names(COR)=cpg_id
    n=names(COR)[order(COR,decreasing = T)][1:num]
    sub[[i]]=data[,n]
    names(sub)[i]=PathNames[i]
  }
  result=list()
  a=lapply(1:length(PathNames), function(x){
    data.frame(
      label=colnames(sub[[x]]),
      value=rep(names(sub)[x],length(sub[[x]]))
    )
  })
  result$vertices=do.call(rbind,a)
  comid=unique(result$vertices$label[duplicated(result$vertices$label)])
  result$vertices=result$vertices[!duplicated(result$vertices$label),]
  rownames(result$vertices)=result$vertices$label
  nonid=setdiff(result$vertices$label,comid)
  x=NA
  y=NA
  xend=NA
  yend=NA
  same.conf=NA
  conf=NA
  names(sub)=NULL
  dd=do.call(cbind,sub)
  dd=dd[,!duplicated(colnames(dd))]
  cor_matrix <- stats::cor(dd)
  upper_tri <- cor_matrix[upper.tri(cor_matrix)]
  n <- nrow(cor_matrix)
  upper_tri_matrix <- matrix(0, n, n)
  upper_tri_matrix[upper.tri(upper_tri_matrix)] <- upper_tri
  cor_matrix=upper_tri_matrix
  colnames(cor_matrix)=colnames(dd)
  rownames(cor_matrix)=colnames(dd)
  df <- data.frame(from = character(n^2), to = character(n^2), Correlation = numeric(n^2))
  count <- 1
  for (i in 1:n) {
    for (j in i:n) {
      df[count, "from"] <- rownames(cor_matrix)[i]
      df[count, "to"] <- rownames(cor_matrix)[j]
      df[count, "Correlation"] <- cor_matrix[i, j]
      count <- count + 1
    }
  }
  df$Correlation=ifelse(df$Correlation>0,df$Correlation,-df$Correlation)
  df=df[df$Correlation>0,]
  
  DF=df[df$Correlation> cutoff,]
  DF$same.conf=ifelse(result$vertices[DF$from,"value"]==result$vertices[DF$to,"value"],1,0)
  DF1=DF[DF$same.conf==1,]
  
  pname=BioM2_pathways_obj$PathwaysResult$id[1:10]
  cor_matrix <- stats::cor(BioM2_pathways_obj$PathwaysMatrix[,pname])
  upper_tri <- cor_matrix[upper.tri(cor_matrix)]
  n <- nrow(cor_matrix)
  upper_tri_matrix <- matrix(0, n, n)
  upper_tri_matrix[upper.tri(upper_tri_matrix)] <- upper_tri
  cor_matrix=upper_tri_matrix
  colnames(cor_matrix)=colnames(BioM2_pathways_obj$PathwaysMatrix[,pname])
  rownames(cor_matrix)=colnames(BioM2_pathways_obj$PathwaysMatrix[,pname])
  df <- data.frame(from = character(n^2), to = character(n^2), Correlation = numeric(n^2))
  count <- 1
  for (i in 1:n) {
    for (j in i:n) {
      df[count, "from"] <- rownames(cor_matrix)[i]
      df[count, "to"] <- rownames(cor_matrix)[j]
      df[count, "Correlation"] <- cor_matrix[i, j]
      count <- count + 1
    }
  }
  df$Correlation=ifelse(df$Correlation>0,df$Correlation,-df$Correlation)
  DF0=df[df$Correlation> 0,]
  DF0=DF0[DF0$Correlation>quantile(DF0$Correlation, probs = 0.75),]
  map=result$vertices[which(!duplicated(result$vertices$value)),]
  rownames(map)=map$value
  DF0$from=map[DF0$from,]$label
  DF0$to=map[DF0$to,]$label
  DF0$same.conf=rep(0,nrow(DF0))
  DF=rbind(DF1,DF0)
  result$edges=DF
  
  fb.igra=graph_from_data_frame(result$edges[,1:2],directed = FALSE)
  V(fb.igra)$conf=result$vertices[V(fb.igra)$name, "value"]
  E(fb.igra)$same.conf=result$edges$same.conf
  E(fb.igra)$lty=ifelse(E(fb.igra)$same.conf == 1, 1, 2)
  
  pic<-ggplot(ggnetwork(fb.igra), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(linetype= as.factor(same.conf)),
               #arrow = arrow(length = unit(6, "pt"), type = "closed") #if directed
               color = "grey50",
               curvature = 0.2,
               alpha=0.8,
               ncp=10,
               linewidth=0.7
    ) +
    geom_nodes(aes(color = conf),
               size = 7,
               alpha=0.5) +
    scale_color_brewer("Pathways",
                       palette = 'Paired') +
    scale_linetype_manual(values = c(2,1)) +
    guides(linetype =  "none") +
    theme_blank()+
    geom_nodes(aes(color = conf),
               size = 4)+labs(title = 'Network Diagram of TOP 10 Pathway-Level Features')+
    theme(legend.text = element_text(family = 'serif',face = 'bold.italic',color = 'grey15'),
          legend.title = element_text(family = 'serif',face = 'bold'),
          plot.title = element_text(family = 'serif',face = 'bold'))+
    geom_mark_ellipse(
      aes(fill=conf,label =conf),
      alpha = 0.2,
      show.legend = F
    )+scale_fill_brewer(palette = 'Paired')+xlim(-0.05,1.05)+ylim(-0.05,1.05)
  return(pic)
  
}