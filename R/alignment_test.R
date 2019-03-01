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


#' @title Evaluate the spacial discribution of the mixed single cells from different data types. (Fast version of neighbor_test)
#' @description This function is used for testing whether the single cells from different data types are mixed well.
#' A fisher's extract test will be performed to test whether the average ratio of input cells to reference cells in the given region is the same as the ratio of the total number of input cells to the total number of reference cells.
#' @param input_data Low dimensional representation of single cell from one data type as the input for matching (e.g., PCs from scRNA-seq data).
#' @param ref_data Low dimensional representation of single cell from another data type as the reference for matching (e.g., PCs from scATAC-seq data).
#' @param dist_scale Scale used to define the radius of the region for testing.
#' @param print_message Flag to print the radius used for the testing.
#' @return
#'  \item{input_count}{The number of cells in the input data for each test.}
#'  \item{ref_count}{The number of cells in the reference data for each test.}
#'  \item{pval}{P-values from fisher's extract test.}
#' @keywords spacial test
#' @examples
#' \dontrun{
#' neighbor_test_p <- neighbor_test_fast(input_data,ref_data,dist_scale=10)$pval
#' }
#' @export
neighbor_test_fast <- function(input_data,ref_data,dist_scale=10,print_message=TRUE){
  
  dist_search <- max(max(dist(input_data)),max(dist(ref_data)))/dist_scale
  
  if(print_message){
    print(paste0("Using radius ",dist_search))
  }
  
  input_expect <- nrow(input_data)
  ref_expect <- nrow(ref_data)
  
  dist_all <- as.matrix(dist(rbind(input_data,ref_data)))
  input_dist_all <- dist_all[c(1:nrow(input_data)),c(1:nrow(input_data))]
  ref_dist_all <- dist_all[c(1:nrow(input_data)),-c(1:nrow(input_data))]
  
  input_count <- mean(rowSums(input_dist_all < dist_search))
  ref_count <- mean(rowSums(ref_dist_all < dist_search))
  
  test_table <- matrix(c(input_count,ref_count,input_expect,ref_expect),
                       nrow = 2,dimnames = list(c("input", "ref"),c("observe", "expect")))
  fisher_pval <- fisher.test(test_table)$p.value
  
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
