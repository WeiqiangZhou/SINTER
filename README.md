## SINTER:a tool for single-cell genomic data integration

### Overview
SINTER is a tool for integrating single-cell genomic data such as scRNA-seq and scATAC-seq.

### Installation
```
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("WeiqiangZhou/SINTER")
```

### How to use


### Example
```
##load SINTER
library(SINTER)

##load example training data
DNase_train <- DNase_hg19_cluster_2000_example
RNA_train <- RNA_hg19_example

##load example scRNA-seq and scATAC-seq data (contains H1, K562, and GM12878 cells)
expr_in <- expr_example
atac_in <- scatac_example

##select variable gene
expr_data_select <- select_gene(expr_in,RNA_train,var_th=0.2,dataplot=TRUE)
atac_data_select <- select_DHS_cluster(atac_in,lib_th=1e3,var_th=0.2,dataplot=TRUE)

atac_data_filter <- atac_data_select$data_select
atac_select_idx <- atac_data_select$select_idx

##obtain optimal parameters
param_opt <- predict_opt(atac_data_filter,expr_data_select$expr_select,DNase_train[atac_select_idx,],expr_data_select$RNA_train_select,
  num_predictor=c(20,25,30),cluster_scale=c(10,20,50),k_range=c(20:29),sigma_range=c(0.01,2),dim=3,dist_scale_in=10,
  subsample=FALSE,tol_er=0.001,ncore=10)

##perform prediction of chromatin accessibility based on gene expression using the optimal parameters
pre_result <- pre_model(expr_data_select$expr_select,DNase_train[atac_select_idx,],expr_data_select$RNA_train_select,
  num_predictor=param_opt$num_predictor_opt,cluster_scale=param_opt$cluster_scale_opt)

pre_result_sd <- pre_result - rowMeans(pre_result)
colnames(pre_result_sd) <- colnames(pre_result)
row.names(pre_result_sd) <- row.names(atac_data_filter)

##run MNN correction
data_MNN <- mnnCorrect(atac_data_filter,pre_result_sd,k=param_opt$k_opt,sigma=param_opt$sigma_opt)

##combine data and perform dimension reduction
data_combine <- cbind(data_MNN$corrected[[1]],data_MNN$corrected[[2]])
colnames(data_combine) <- c(colnames(atac_data_filter),colnames(pre_result_sd))

data_pc <- prcomp(t(data_combine),center = T, scale = T)$x

##generate a plot to showing the matching results
cell_idx <- rep(NA,ncol(data_combine))
cell_idx[grep("H1",colnames(data_combine))] <- "H1"
cell_idx[grep("GM12878",colnames(data_combine))] <- "GM12878"
cell_idx[grep("K562",colnames(data_combine))] <- "K562"

group_idx <- rep("scRNA",ncol(data_combine))
group_idx[grep("ATAC",colnames(data_combine))] <- "scATAC"
group_idx <- as.factor(group_idx)

plot_data <- data.frame(PC1=data_pc[,1],PC2=data_pc[,2],cell=cell_idx,group=group_idx)

p <- ggplot(plot_data,aes(x=PC1,y=PC2,color=cell,shape=group)) +
     geom_point(size=3)+
     scale_shape_manual(values=c(17,15))+
     theme_bw()+
     theme(axis.text.x = element_text(color="black",size=20),axis.text.y = element_text(color="black",size=20),
     axis.title.x = element_text(color="black",size=20),axis.title.y = element_text(color="black",size=20),
     legend.title=element_blank(),legend.text=element_text(size=20),legend.position="top")
print(p)
```
