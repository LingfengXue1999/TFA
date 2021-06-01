scTFA_calculation=function(data,GRN,zscore=T,
                                              gene_weight_method=NULL,defined_weight_matrix=NULL,
                                              prefix=NULL,prefix_sep="_",verbose=T){
  # scTFA_calculation_zscore_v20200709_1
  # data: numeric matrix or data.frame, gene*cell, gene name will be toupper-ed, cell number should be more than 
  #       1 for z-score normalization (if zscore=T)
  # GRN: TF*gene, gene name will be toupper-ed 
  # zscore: whether to conduct z-score normalization for each gene, default: T
  # gene_weight_method: weight method for each gene in data. if "mean_expression", then weight is mean expression; 
  #                     if "defined", then weight matrix is defined_weight_matrix, weighted data is 
  #                     defined_weight_matrix*normalized_data
  # defined_weight_matrix: defined_weight_matrix (gene*gene) if gene_weight_method is "defined", rownames and 
  #                        colnames is the same as rownames of data
  # prefix: prefix add to TF names of TFA
  # prefix_sep: sep between prefix and names
  # return: TFA, TF*cell
  # verbose: print progress
  
  GA=data
  rownames(GA)=toupper(rownames(GA))
  
  if(!is.null(gene_weight_method)){
    if(gene_weight_method=="mean_expression"){
      weight_vector=apply(GA,1,mean)
      weight_matrix=diag(weight_vector)
      dimnames(weight_matrix) = list(names(weight_vector),names(weight_vector))
    }else if(gene_weight_method=="defined"){
      if(!is.null(defined_weight_matrix)){
        weight_matrix=defined_weight_matrix
      }else{
        stop("defined_weight_matrix should be provided when gene_weight_method is 'defined'")
      }
    }
  }else{
    weight_matrix=diag(1,nrow(GA),nrow(GA))
    dimnames(weight_matrix)=list(rownames(GA),rownames(GA))
  }
  
  if(zscore){
    zscore_function=function(x){
      if(sd(x)==0){
        return(x-mean(x))
      }else{
        return(scale(x,center = TRUE, scale = TRUE))
      }
    }
    GA_zscore=t(apply(GA,1,zscore_function))
    colnames(GA_zscore)=colnames(GA)
    GA=GA_zscore
  }
  GA=weight_matrix %*% GA
  
  colnames(GRN)=toupper(colnames(GRN))
  GRN=GRN[,intersect(rownames(GA),colnames(GRN))]
  GRN=GRN[rowSums(GRN)>0,]
  GRN_norm=t(apply(GRN,1,function(x) x/sum(x)))
  TFA=GRN_norm %*% GA[colnames(GRN_norm),]
  # # #TFA=TFA[rowSums(TFA)>0,,drop=F]
  if(!is.null(prefix)){
    rownames(TFA)=paste(prefix,rownames(TFA),sep=prefix_sep)
  }
  return(TFA)
}