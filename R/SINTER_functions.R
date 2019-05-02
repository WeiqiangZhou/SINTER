#' @import gam
#' @import matrixStats
#' @import preprocessCore
#' @import scran
#' @import ggplot2
#' @import parallel
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import irlba

#' @title Summarize single-cell ATAC-seq data according to ENCODE cluster features
#' @description This function is used for summarizing scATAC-seq data according to ENCODE clusters.
#' @param path Directory to bam files of scATAC-seq data
#' @param gr_ENCL GRanges object for ENCODE clusters
#' @param type Sequencing read type. "paired" or "single" end read.
#' @return
#'  \item{data}{Log2-transformed summarized data matrix.}
#'  \item{lib}{Total read counts (i.e., library size) for each single cell.}
#' @keywords summarize scATAC-seq
#' @examples
#' \dontrun{
#' scATAC_data <- ENCLfunc("/bam_file_path",gr_ENCL,type="paired")
#' }
#' @export
ENCLfunc <- function(path,gr_ENCL,type="paired") {

  allf <- list.files(path,pattern="\\.bam$")
  grf <- sapply(allf,function(f) {
    if (type=="paired") {
      GRanges(readGAlignmentPairs(paste0(path,"/",f)))
    } else if (type=="single") {
      GRanges(readGAlignments(paste0(path,"/",f)))
    }
  })
  libsize_org <- sapply(grf,length)
  libsize <- libsize_org/median(libsize_org)
  count <- sapply(grf,function(i) {
    countOverlaps(gr_ENCL,i)
  })
  data_mat <- log2(sweep(count,2,libsize,"/")+1)
  return(list(data=data_mat,lib=libsize_org))
}


#' @title Select variable gene from single-cell RNA-seq data
#' @description This function is used for selecting variable gene from scRNA-seq data.
#' @param expr Gene expression matrix of scRNA-seq (normalized counts).
#' @param RNA_train_all Gene expression matrix of bulk RNA-seq from ENCODE.
#' @param filter_th Threshold for filtering out genes with no expression in the given percentage of cells.
#' @param var_th Threshold for filtering out genes with low variability in the given percentage of genes.
#' @param dataplot Plot the mean and variance of all genes and the selected genes are highlighted.
#' @return
#'  \item{expr_select}{Log2-transformed gene expression matrix of scRNA-seq for the selected genes.}
#'  \item{RNA_train_select}{Gene expression matrix of bulk RNA-seq for the selected genes.}
#' @keywords variable gene scRNA-seq
#' @examples
#' \dontrun{
#' expr_data_select <- select_gene(expr_in,RNA_train_all,filter_th=0.1,var_th=0.2,dataplot=TRUE)
#' }
#' @export
select_gene <- function(expr,RNA_train_all,filter_th=0.1,var_th=0.2,dataplot=FALSE){

  retain_idx <- which(rowMeans(expr > 0) > filter_th)
  expr_filter <- expr[retain_idx,]

  row_mean <- rowMeans(expr_filter)
  row_var <- rowVars(as.matrix(expr_filter))
  data_fit <- data.frame(X=log2(row_mean+1),Y=log2(row_var+1))
  fit_model <- gam(Y~s(X),data=data_fit)
  hyper_var <- log2(row_var+1) - fit_model$fitted.values
  hyper_var_sort <- sort(hyper_var,decreasing=TRUE)
  hyper_var_th <- hyper_var_sort[round(length(hyper_var_sort)*var_th)]
  var_idx <- which(hyper_var > hyper_var_th)

  if(dataplot){
    plot(log2(row_mean+1),log2(row_var+1),pch=19)
    points(log2(row_mean+1),fit_model$fitted.values,pch=19,col="blue")
    points(log2(row_mean+1)[var_idx],log2(row_var+1)[var_idx],pch=19,col="red")
  }

  expr_select <- expr_filter[var_idx,]
  gene_names <- row.names(expr_select)
  gene_names <- sapply(gene_names,function(x) sub("\\..*","",x))

  train_gene <- row.names(RNA_train_all)
  train_gene <- sapply(train_gene,function(x) sub("\\..*","",x))

  match_idx <- match(gene_names,train_gene)
  print(paste0(length(which(!is.na(match_idx)))," genes are selected"))
  return(list(expr_select=log2(expr_select[which(!is.na(match_idx)),]+1),RNA_train_select=RNA_train_all[match_idx[which(!is.na(match_idx))],]))
}


#' @title Select variable ENCODE cluster features from scATAC-seq data
#' @description This function is used for selecting variable ENCODE cluster features from scATAC-seq data.
#' @param DHS_data Summarized ENCODE cluster features from scATAC-seq data.
#' @param lib_th Threshold for filtering out cells with low total number of reads.
#' @param var_th Threshold for filtering out ENCODE cluster features with low variability in the given percentage of features.
#' @param dataplot Plot the mean and variance of all ENCODE cluster features and the selected features are highlighted.
#' @return
#'  \item{data_select}{Selected ENCODE cluster features.}
#'  \item{select_idx}{Index for the selected ENCODE cluster features.}
#' @keywords variable feature scATAC-seq
#' @examples
#' \dontrun{
#' scATAC_data_select <- select_DHS_cluster(DHS_data,lib_th=1e3,var_th=0.2,dataplot=FALSE)
#' }
#' @export
select_DHS_cluster <- function(DHS_data,lib_th=1e3,var_th=0.2,dataplot=FALSE){

  filter_idx <- DHS_data$lib < lib_th
  DHS_data_filter <- DHS_data$data[,!filter_idx]

  DHS_data_sd <- DHS_data_filter - rowMeans(DHS_data_filter)
  colnames(DHS_data_sd) <- colnames(DHS_data_filter)

  row_mean <- rowMeans(DHS_data_filter)
  row_var <- rowVars(as.matrix(DHS_data_filter))
  data_fit <- data.frame(X=row_mean,Y=row_var)
  fit_model <- gam(Y~s(X),data=data_fit)
  hyper_var <- row_var - fit_model$fitted.values
  hyper_var_sort <- sort(hyper_var,decreasing=TRUE)
  hyper_var_th <- hyper_var_sort[round(length(hyper_var_sort)*var_th)]

  var_idx <- which(hyper_var > hyper_var_th)
  DHS_data_ex <- DHS_data_sd[var_idx,]

  if(dataplot){
    plot(row_mean,row_var,pch=19)
    points(row_mean,fit_model$fitted.values,pch=19,col="blue")
    points(row_mean[var_idx],row_var[var_idx],pch=19,col="red")
  }

  print(paste0(ncol(DHS_data_ex)," cells are retained and ",length(var_idx)," DHS clusters are selected"))
  return(list(data_select=DHS_data_ex,select_idx=var_idx))
}


#' @title Linear regression function based on sure independence screening
#' @description This function is used for building linear regression models based on the top N predictors that are most correlated with reponse.
#' @param DNase_train ENCODE cluster features from DNase-seq data for building the regression model.
#' @param RNA_train_mean Gene cluster mean from ENCODE RNA-seq data for building the regression model.
#' @param RNA_test_mean Gene cluster mean from scRNA-seq data for making predictions.
#' @param top_n Number of predictors used in the regression model.
#' @return
#'  \item{y_pre}{A vector of predicted ENCODE cluster features.}
#' @keywords prediction
#' @export
cluster_regression_topn <- function(DNase_train,RNA_train_mean,RNA_test_mean,top_n){
  y_cor <- cor(DNase_train,t(RNA_train_mean))
  y_cor[is.na(y_cor)] <- 0

  if(top_n > 1){
    max_idx <- sort(y_cor, decreasing=TRUE, index.return=TRUE)$ix[1:top_n]      #if the N>1, select the N most correlated predictiors
    data_train <- data.frame(y=DNase_train, t(RNA_train_mean[max_idx,]))
    data_test <- data.frame(t(RNA_test_mean[max_idx,]))
    fit <- lm(y~.,data_train)
    y_pre <- predict(fit,data_test)
  }else{
    max_idx <- which(y_cor == max(y_cor))[1]      #if N=1, select the most correlated predictor
    data_train <- data.frame(y=DNase_train, x=RNA_train_mean[max_idx,])
    data_test <- data.frame(x=RNA_test_mean[max_idx,])
    fit <- lm(y~x, data_train)
    y_pre <- predict(fit,data_test)
  }

  return(y_pre)
}


#' @title Prediction of ENCODE cluster features based on scRNA-seq data
#' @description This function is used for predicting ENCODE cluster features based on scRNA-seq data. The scRNA-seq data are first clustered into gene clusters and the cluster means are used as predictors.
#' @param expr_select Input gene expression data from scRNA-seq.
#' @param DNase_train ENCODE cluster features from DNase-seq data for building the regression model.
#' @param RNA_train Gene expression from ENCODE RNA-seq data for building the regression model.
#' @param num_predictor Number of predictors used in the prediction model.
#' @param cluster_scale The scale to determine the number of gene clusters. The number of gene clusters is obtained by [the number of genes]/[cluster_scale].
#' @param seed Set the seed in kmeans clustering for reproducible results.
#' @return
#'  \item{Y_pre}{A matrix of predicted ENCODE cluster features.}
#' @keywords prediction
#' @examples
#' \dontrun{
#' Y_pre <- pre_model(expr_select,DNase_train,RNA_train,num_predictor=10,cluster_scale=10,seed=12345)
#' }
#' @export
pre_model <- function(expr_select,DNase_train,RNA_train,num_predictor=10,cluster_scale=10,seed=12345){

  ##standardize each gene
  RNA_all_sd <- t(apply(cbind(expr_select,RNA_train),1,scale))
  RNA_train_sd <- RNA_all_sd[,-c(1:ncol(expr_select))]

  ##cluster genes based on RNA_train
  cluster_num <- round(nrow(RNA_train_sd)/cluster_scale)
  if(num_predictor >= cluster_num){
    num_predictor <- cluster_num - 1
  }
  print(paste0("Using ",cluster_num," clusters and ",num_predictor," predictors for prediction"))

  set.seed(seed)
  gene_cluster <- kmeans(RNA_train_sd, centers=cluster_num, nstart=10, iter.max=50)$cluster

  RNA_all_mean <- sapply(c(min(gene_cluster):max(gene_cluster)),function(x){
    if(length(which(gene_cluster == x)) == 1){
      RNA_all_sd[which(gene_cluster == x),]
    }
    else{
      colMeans(RNA_all_sd[which(gene_cluster == x),])
    }
  })

  RNA_all_mean <- t(RNA_all_mean)
  RNA_all_mean_norm <- normalize.quantiles(RNA_all_mean)
  colnames(RNA_all_mean_norm) <- c(colnames(expr_select),colnames(RNA_train))

  expr_select_norm <- RNA_all_mean_norm[,c(1:ncol(expr_select))]
  RNA_train_norm <- RNA_all_mean_norm[,-c(1:ncol(expr_select))]

  Y_pre <- matrix(data=NA,nrow=nrow(DNase_train),ncol=ncol(expr_select_norm))
  colnames(Y_pre) <- colnames(expr_select_norm)

  pb = txtProgressBar(min = 0, max = nrow(DNase_train), initial = 0, style = 3)
  for(i in 1:nrow(DNase_train)){
    setTxtProgressBar(pb,i)
    Y_pre[i,] <- cluster_regression_topn(DNase_train[i,],RNA_train_norm,expr_select_norm,num_predictor)
  }
  close(pb)

  return(Y_pre)

}


#' @title Evaluate the spacial discribution of the mixed single cells from different data types
#' @description This function is used for testing whether the single cells from different data types are mixed well.
#' For each cell in the input data, a fisher's extract test will be performed to test whether the ratio of input cells to reference cells in the given region is the same as the ratio of the total number of input cells to the total number of reference cells.
#' @param input_data Low dimensional representation of single cell from one data type as the input for matching (e.g., PCs from scRNA-seq data).
#' @param ref_data Low dimensional representation of single cell from another data type as the reference for matching (e.g., PCs from scATAC-seq data).
#' @param dist_scale Scale used to define the radius of the region for testing.
#' @param print_message Flag to print the radius used for the testing.
#' @return
#'  \item{input_count}{The number of cells in the input data for each test.}
#'  \item{ref_count}{The number of cells in the reference data for each test.}
#'  \item{pval}{P-values from fisher's extract tests.}
#' @keywords spacial test
#' @examples
#' \dontrun{
#' neighbor_test_p <- neighbor_test(input_data,ref_data,dist_scale=10)$pval
#' }
#' @export
neighbor_test <- function(input_data,ref_data,dist_scale=10,print_message=TRUE){

  dist_search <- max(max(dist(input_data)),max(dist(ref_data)))/dist_scale

  if(print_message){
    print(paste0("Using radius ",dist_search))
  }

  input_expect <- nrow(input_data)
  ref_expect <- nrow(ref_data)

  input_count <- rep(NA,nrow(input_data))
  ref_count <- rep(NA,nrow(input_data))

  dist_all <- as.matrix(dist(rbind(input_data,ref_data)))
  input_dist_all <- dist_all[c(1:nrow(input_data)),c(1:nrow(input_data))]
  ref_dist_all <- dist_all[c(1:nrow(input_data)),-c(1:nrow(input_data))]

  fisher_pval <- sapply(1:nrow(input_data),function(i){
    ##self included
    input_count[i] <- length(which(input_dist_all[i,] < dist_search))
    ref_count[i] <- length(which(ref_dist_all[i,] < dist_search))

    test_table <- matrix(c(input_count[i],ref_count[i],input_expect,ref_expect),
                  nrow = 2,dimnames = list(c("input", "ref"),c("observe", "expect")))
    fisher.test(test_table)$p.value
  })

  return(list(input_count=input_count,ref_count=ref_count,pval=fisher_pval))
}


#' @title Evaluate the spacial discribution of the mixed single cells from different data types while given the cell type labels
#' @description This function is used for evaluating whether single cells from the same cell type are mixed well.
#' For each cell in the input data, a fisher's extract test will be performed to test whether the ratio of input cells to reference cells in the given region is the same as the ratio of the total number of input cells to the total number of reference cells while only single cells from the same cell type are considered.
#' @param input_data Low dimensional representation of single cell from one data type as the input for matching (e.g., PCs from scRNA-seq data).
#' @param ref_data Low dimensional representation of single cell from another data type as the reference for matching (e.g., PCs from scATAC-seq data).
#' @param input_mem Cell type label for the input_data.
#' @param ref_mem Cell type label for the ref_data.
#' @param dist_scale Scale used to define the radius of the region for testing.
#' @param print_message Flag to print the radius used for the testing.
#' @return
#'  \item{test_stat}{P-values from fisher's extract tests.}
#' @keywords spacial test evaluation
#' @examples
#' \dontrun{
#' test_stat <- eval_neighbor_test(input_data,ref_data,input_mem,ref_mem,dist_scale=10,print_message=TRUE)
#' }
#' @export
eval_neighbor_test <- function(input_data,ref_data,input_mem,ref_mem,dist_scale=10,print_message=TRUE){

  dist_search <- max(max(dist(input_data)),max(dist(ref_data)))/dist_scale
  test_stat <- rep(NA,nrow(input_data))
  input_expect <- nrow(input_data)
  ref_expect <- nrow(ref_data)

  if(print_message){
    print(paste0("Using radius ",dist_search))
  }

  for(i in 1:nrow(input_data)){
    ref_dist <- sapply(1:nrow(ref_data),function(x){
      sqrt(sum((input_data[i,] - ref_data[x,]) ^ 2))})
    input_dist <- sapply(1:nrow(input_data),function(x){
      sqrt(sum((input_data[i,] - input_data[x,]) ^ 2))})
    ref_count <- length(intersect(which(ref_dist < dist_search),which(ref_mem == input_mem[i])))
    input_count <- length(intersect(which(input_dist < dist_search),which(input_mem == input_mem[i])))

    test_table <- matrix(c(input_count,ref_count,input_expect,ref_expect),
                  nrow = 2,dimnames = list(c("input", "ref"),c("observe", "expect")))
    test_stat[i] <- fisher.test(test_table)$p.value
  }
  return(test_stat)
}


#' @title Run MNN to correct for platform effects
#' @description This function is used run MNN to correct for platform effects between the experimental and predicted scATAC-seq data.
#' @param atac_data scATAC-seq data for matching.
#' @param pre_result Predicted scATAC-seq data based on scRNA-seq.
#' @param k The number of mutual nearest neighbor in MNN.
#' @param sigma The bandwidth of the Gaussian smoothing kernel used to compute the correction vector.
#' @param MNN_ref A flag to determine which data type is used as reference in MNN. Select from "scATAC" and "scRNA".
#' @return
#'  \item{data_combine}{The combined data matrix from all single cells.}
#' @keywords MNN
#' @examples
#' \dontrun{
#' data_combine <- run_MNN(atac_data,pre_result,k=param_opt$k_opt,sigma=param_opt$sigma_opt,MNN_ref="scATAC")
#' }
#' @export
run_MNN <- function(atac_data,pre_result,k,sigma,MNN_ref="scATAC",...){

  pre_result_sd <- pre_result - rowMeans(pre_result)
  colnames(pre_result_sd) <- colnames(pre_result)
  row.names(pre_result_sd) <- row.names(atac_data)
  
  if(MNN_ref == "scATAC"){
    data_MNN <- mnnCorrect(atac_data,pre_result_sd,k=k,sigma=sigma,...)
    data_combine <- cbind(data_MNN$corrected[[1]],data_MNN$corrected[[2]])
  }
  else{
    data_MNN <- mnnCorrect(pre_result_sd,atac_data,k=k,sigma=sigma,...)
    data_combine <- cbind(data_MNN$corrected[[2]],data_MNN$corrected[[1]])
  }
  
  colnames(data_combine) <- c(colnames(atac_data),colnames(pre_result_sd))
  return(data_combine)
}
