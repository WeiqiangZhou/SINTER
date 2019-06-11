#' @export
opt_fun <- function(param_est,k_in,ref_data,input_data,num_pc,dist_scale_in,fast=FALSE){
  
  MNN_result <- mnnCorrect(ref_data,input_data,k=k_in,sigma=param_est,svd.dim=0)
  data_combine <- cbind(MNN_result$corrected[[1]],MNN_result$corrected[[2]])
  data_pc <- prcomp_irlba(t(data_combine),n = num_pc, center = T, scale = T)$x
  if(fast==FALSE){
    neighbor_test_p <- neighbor_test(data_pc[-c(1:ncol(ref_data)),1:num_pc],data_pc[1:ncol(ref_data),1:num_pc],dist_scale=dist_scale_in,print_message=FALSE)$pval
    output <- mean(neighbor_test_p)
  }
  else{
    output <- neighbor_test_fast(data_pc[-c(1:ncol(ref_data)),1:num_pc],data_pc[1:ncol(ref_data),1:num_pc],dist_scale=dist_scale_in,print_message=FALSE)$pval
  }
  return(output)
}


#' @title Search for the optimal parameters for MNN
#' @description This function is used for searching for the optimal parameters used in MNN.
#' The goal is find parameters that maximize the average of p-values from the neighbor_test function.
#' @param data_ref Single cell data from one data type as the reference for matching (e.g., scATAC-seq data).
#' @param data_in Single cell data from another data type as the input for matching (e.g., scRNA-seq data).
#' @param k_range Searching space for k, the number of mutual nearest neighbor in MNN.
#' @param sigma_range Searching space for sigma, the bandwidth of the Gaussian smoothing kernel used to compute the correction vector.
#' @param dim Number of dimension used for matching the single cells. For example, the number of principal components.
#' @param dist_scale_in Scale used to define the radius of the region for testing.
#' @param tol_er The desired accuracy in function optimize.
#' @param ncore Number of CPU cores used for parallel processing. Use ncore = 1 to run the function without parallel processing.
#' @param fast A flag indicates whether or not to use a fast neighbor_test.
#' @return
#'  \item{k_opt}{The optimal value for k.}
#'  \item{sigma_opt}{The optimal value for sigma.}
#'  \item{max_obj}{The average p-value based on the optimal parameters.}
#' @keywords MNN optimization
#' @examples
#' \dontrun{
#' result_opt <- mnn_opt(data_ref,data_in,k_range=c(16:25),sigma_range=c(0.01,1),dim=3,dist_scale_in=10,tol_er=0.001,ncore=10)
#' }
#' @export
mnn_opt <- function(data_ref,data_in,k_range=c(16:25),sigma_range=c(0.01,1),dim=3,dist_scale_in=10,tol_er=0.001,ncore=10,fast=FALSE){
  
  if(ncore > 1){
    opt_result <- mclapply(k_range,function(x){
      print(paste0("Using k=",x," mutual nearest neighbors"))
      optimize(opt_fun, sigma_range, x, data_ref, data_in, dim, dist_scale_in, fast, maximum=TRUE, tol = tol_er)
    },mc.cores=ncore)
  }
  else{
    opt_result <- lapply(k_range,function(x){
      print(paste0("Using k=",x," mutual nearest neighbors"))
      optimize(opt_fun, sigma_range, x, data_ref, data_in, dim, dist_scale_in, fast, maximum=TRUE, tol = tol_er)
    })
  }
  
  opt_result_obj <- sapply(opt_result,function(x)x$objective)
  opt_result_sigma <- sapply(opt_result,function(x)x$maximum)
  max_idx <- which(opt_result_obj==max(opt_result_obj))[1]
  
  print(paste0("best k=",k_range[max_idx]))
  print(paste0("best sigma=",opt_result_sigma[max_idx]))
  print(paste0("max objective=",opt_result_obj[max_idx]))
  
  return(list(k_opt=k_range[max_idx],sigma_opt=opt_result_sigma[max_idx],max_obj=opt_result_obj[max_idx]))
}


#' @title Search for the optimal parameters for SINTER
#' @description This function is used for searching for the optimal parameters used in SINTER.
#' The goal is find parameters that maximize the average of p-values from the neighbor_test function.
#' @param atac_data scATAC-seq data for matching.
#' @param expr_data scRNA-seq data for matching.
#' @param DNase_train ENCODE cluster features from DNase-seq data for building the regression model.
#' @param RNA_train Gene expression from ENCODE RNA-seq data for building the regression model.
#' @param num_predictor Searching space for number of predictors used in the regression model.
#' @param cluster_scale Searching space for the scale to determine the number of gene clusters.
#' @param k_range Searching space for k, the number of mutual nearest neighbor in MNN if flag MNN_opt==TRUE.
#' @param sigma_range Searching space for sigma, the bandwidth of the Gaussian smoothing kernel used to compute the correction vector if flag MNN_opt==TRUE.
#' @param k_in Setting K, the number of mutual nearest neighbor in MNN if flag MNN_opt!=TRUE.
#' @param sigma_in Setting sigma, the bandwidth of the Gaussian smoothing kernel used to compute the correction vector if flag MNN_opt!=TRUE.
#' @param dim Number of dimension used for matching the single cells. For example, the number of principal components.
#' @param dist_scale_in Scale used to define the radius of the region for testing.
#' @param subsample A percentage value to determine whether the paramter searching should be done in a subset of cells instead of using all cells. Set subsample=FALSE to use all cells.
#' @param MNN_opt A flag to determine whether the parameters search should be performed for MNN.
#' @param fast A flag indicates whether or not to use a fast neighbor_test.
#' @param MNN_ref A flag to determine which data type is used as reference in MNN. Select from "scATAC" and "scRNA".
#' @param tol_er The desired accuracy in function optimize.
#' @param ncore Number of CPU cores used for parallel processing. Use ncore = 1 to run the function without parallel processing.
#' @param seed The seed used for subsampling if subsample!=FALSE.
#' @return
#'  \item{num_predictor_opt}{The optimal value for number of predictors.}
#'  \item{cluster_scale_opt}{The optimal value for cluster scale.}
#'  \item{k_opt}{The optimal value for k.}
#'  \item{sigma_opt}{The optimal value for sigma.}
#'  \item{max_obj}{The average p-value based on the optimal parameters.}
#' @keywords MNN optimization
#' @examples
#' \dontrun{
#' result_opt <- predict_opt(atac_data,expr_data,DNase_train,RNA_train,num_predictor=c(25,25,30),cluster_scale=c(10,20,50),k_range=c(20:29),sigma_range=c(0.01,1),
#' k_in=20,sigma_in=0.1,dim=3,dist_scale_in=10,subsample=FALSE,MNN_opt=TRUE,MNN_ref="scATAC",tol_er=0.001,ncore=10,seed=12345)
#' }
#' @export
predict_opt <- function(atac_data,expr_data,DNase_train,RNA_train,num_predictor=c(25,25,30),cluster_scale=c(10,20,50),k_range=c(20:29),sigma_range=c(0.01,1),
                        k_in=20,sigma_in=0.1,dim=3,dist_scale_in=10,subsample=FALSE,MNN_opt=TRUE,fast=FALSE,MNN_ref="scATAC",tol_er=0.001,ncore=10,seed=12345){
  
  if(subsample!=FALSE){
    set.seed(seed)
    sample_idx_input <- sample(c(1:ncol(expr_data)),round(ncol(expr_data)*subsample))
    expr_data <- expr_data[,sample_idx_input]
    
    sample_idx_ref <- sample(c(1:ncol(atac_data)),round(ncol(atac_data)*subsample))
    atac_data <- atac_data[,sample_idx_ref]
  }
  
  cluster_result <- lapply(cluster_scale,function(y){
    predictor_result <- lapply(num_predictor,function(x){
      pre_result <- pre_model(expr_data,DNase_train,RNA_train,num_predictor=x,cluster_scale=y)
      pre_result_sd <- pre_result - rowMeans(pre_result)
      colnames(pre_result_sd) <- colnames(pre_result)
      row.names(pre_result_sd) <- row.names(atac_data)
      
      if(MNN_opt == TRUE){
        if(MNN_ref == "scATAC"){
          return(unlist(mnn_opt(atac_data,pre_result_sd,k_range,sigma_range,dim,dist_scale_in,tol_er,ncore,fast)))
        }
        else{
          return(unlist(mnn_opt(pre_result_sd,atac_data,k_range,sigma_range,dim,dist_scale_in,tol_er,ncore,fast)))
        }
      }
      else{
        if(MNN_ref == "scATAC"){
          objective_value <- opt_fun(sigma_in,k_in,atac_data,pre_result_sd,dim,dist_scale_in,fast)
        }
        else{
          objective_value <- opt_fun(sigma_in,k_in,pre_result_sd,atac_data,dim,dist_scale_in,fast)
        }
        print(paste0("objective = ",objective_value))
        return(objective_value)
      }
    })
    
    if(MNN_opt == TRUE){
      predictor_result_combine <- do.call(rbind,predictor_result)
      max_idx <- which(predictor_result_combine[,3] == max(predictor_result_combine[,3]))
      return(c(num_predictor[max_idx],predictor_result_combine[max_idx,]))
    }
    else{
      predictor_result_combine <- do.call(c,predictor_result)
      max_idx <- which(predictor_result_combine == max(predictor_result_combine))
      return(c(num_predictor[max_idx],predictor_result_combine[max_idx]))
    }
  })
  
  if(MNN_opt == TRUE){
    cluster_result_combine <- do.call(rbind,cluster_result)
    max_idx <- which(cluster_result_combine[,4] == max(cluster_result_combine[,4]))
    
    print("Final results:")
    print(paste0("prediction: best number of predictor = ",as.integer(cluster_result_combine[max_idx,1])))
    print(paste0("prediction: best cluster scale = ",cluster_scale[max_idx]))
    print(paste0("mnn: best k = ",cluster_result_combine[max_idx,2]))
    print(paste0("mnn: best sigma = ",cluster_result_combine[max_idx,3]))
    print(paste0("max objective = ",cluster_result_combine[max_idx,4]))
    
    result_out <- data.frame(num_predictor_opt=as.integer(cluster_result_combine[max_idx,1]),cluster_scale_opt=cluster_scale[max_idx],
                             k_opt=cluster_result_combine[max_idx,2],sigma_opt=cluster_result_combine[max_idx,3],max_obj=cluster_result_combine[max_idx,4])
    row.names(result_out) <- NULL
    return(result_out)
  }
  else{
    cluster_result_combine <- do.call(rbind,cluster_result)
    max_idx <- which(cluster_result_combine[,2] == max(cluster_result_combine[,2]))
    
    print("Final results:")
    print(paste0("prediction: best number of predictor = ",as.integer(cluster_result_combine[max_idx,1])))
    print(paste0("prediction: best cluster scale = ",cluster_scale[max_idx]))
    print(paste0("mnn: using k = ",k_in))
    print(paste0("mnn: using sigma = ",sigma_in))
    print(paste0("max objective = ",cluster_result_combine[max_idx,2]))
    
    result_out <- data.frame(num_predictor_opt=as.integer(cluster_result_combine[max_idx,1]),cluster_scale_opt=cluster_scale[max_idx],
                             k_opt=k_in,sigma_opt=sigma_in,max_obj=cluster_result_combine[max_idx,2])
    row.names(result_out) <- NULL
    return(result_out)
  }
  
}
