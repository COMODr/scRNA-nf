#seurat_data <- newimport(SeuratObject)

###
newimport <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts
    
    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    
    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      
      message("This Seurat object doesn't provide any meta data");
      pd
    })
    
    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }
    
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    valid_data <- data[, row.names(pd)]
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
        
      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }
    
    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
      
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
    
  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")
    
    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 
      
    } else {
      mist_list <- list()
    }
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list
    
    monocle_cds@auxOrderingData$scran <- mist_list
    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(monocle_cds)
}



if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(! require(monocle)) BiocManager::install("monocle")
if(! require(BiocGenerics)) BiocManager::install("BiocGenerics")
if(! require(Seurat)) BiocManager::install("Seurat")
if(! require(dplyr)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)

option_list <- list(
  make_option(c("-C", "--CRdir"), type = "character", default = ".", action = "store", help = "The directory where the Seurat analysis results store"
  ),
  make_option(c("-S", "--sample"), type = "character", default = FALSE, action = "store", help = "Sample name"
  )
)
opt = parse_args(OptionParser(option_list = option_list, usage = "7-Trajectories.R --CRdir <CellRanger results path> --sample <Sample name>"))


##### 7.Trajectories
rslt <- opt$CRdir  
Sample_Name <- opt$sample 



###加载Seurat对象
Seurat_object <- list()
for (i in 1:length(Sample_Name)){
  Seurat_object[[i]] <- readRDS(paste0(rslt,"/",Sample_Name,"_seurat/",Sample_Name,".rds"))
}


if(!dir.exists(paste0(rslt, '/', Sample_Name, '_Trajectories'))){dir.create(paste0(rslt, '/', Sample_Name, '_Trajectories'), recursive = TRUE)}
HSMM <- newimport(Seurat_object[[i]], import_all = TRUE)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM, cores=8, relative_expr = TRUE)
HSMM <- detectGenes(HSMM, min_expr = 3 )
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
head(pData(HSMM))

##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
ordering_genes <- plot_ordering_genes(HSMM)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/ordering_genes.png'), plot = ordering_genes, width = 8, height = 5)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/ordering_genes.pdf'), plot = ordering_genes, width = 8, height = 5)


HSMM <- HSMM
#降维
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
#排序
HSMM <- orderCells(HSMM)
plot0 <- plot_cell_trajectory(HSMM)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/cell_trajectory.png'), plot = plot0, width = 8, height = 5)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/cell_trajectory.pdf'), plot = plot0, width = 8, height = 5)
#有点问题

#State轨迹分布图
plot1 <- plot_cell_trajectory(HSMM, color_by = "State") + facet_wrap(~State, nrow = 1)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/cell_trajectory_stat.png'), plot = plot1, width = 15, height = 5)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/cell_trajectory_stat.pdf'), plot = plot1, width = 15, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(HSMM, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 2)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/cell_trajectory_cluster.png'), plot = plot2, width = 20, height = 10)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/cell_trajectory_cluster.pdf'), plot = plot2, width = 20, height = 10)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/Pseudotime_cell_trajectory.png'), plot = plot3, width = 8, height = 5)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/Pseudotime_cell_trajectory.pdf'), plot = plot3, width = 8, height = 5)
##保存结果
write.csv(pData(HSMM), paste0(rslt, '/', Sample_Name, '_Trajectories/cell_Pseudotime.csv'))

lung_genes <- row.names(subset(fData(HSMM), gene_short_name %in% head(disp.genes,3)))
branched_pseudotime <- plot_genes_branched_pseudotime(HSMM[lung_genes,],
                                                      branch_point = 1,
                                                      color_by = "seurat_clusters",
                                                      ncol = 1)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/genes_branched_pseudotime.png'), plot = branched_pseudotime, width = 8, height = 5)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/genes_branched_pseudotime.pdf'), plot = branched_pseudotime, width = 8, height = 5)

genes_in_pseudotime <- plot_genes_in_pseudotime(HSMM[lung_genes,], color_by = "seurat_clusters")
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/genes_in_pseudotime.png'), plot = genes_in_pseudotime, width = 8, height = 5)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/genes_in_pseudotime.pdf'), plot = genes_in_pseudotime, width = 8, height = 5)

to_be_tested <- row.names(fData(HSMM))
cds_subset <- HSMM[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res <- diff_test_res[,c("gene_short_name", "pval", "qval")]
write.csv(diff_test_res, paste0(rslt, '/', Sample_Name, '_Trajectories/gene_related_to_branch.csv'))


# 选择myogenesis相关的基因集
all_markers <- FindAllMarkers(object = Seurat_object[[i]],
                              assay = 'RNA',#设置assay为RNA#
                              slot = 'counts', verbose = TRUE)
marker_genes <- row.names(subset(all_markers, cluster==0))
diff_test_res <- differentialGeneTest(HSMM[marker_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/branch_dependent_gene_heatmap.png'), 
       plot = plot_genes_branched_heatmap(HSMM[sig_gene_names,],
              branch_point = 1,
              num_clusters = 4,
              cores = 1,
              use_gene_short_name = T,
              show_rownames = T), width = 12, height = 12)
ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/branch_dependent_gene_heatmap.pdf'), 
       plot = plot_genes_branched_heatmap(HSMM[sig_gene_names,],
           branch_point = 1,
           num_clusters = 4,
           cores = 1,
           use_gene_short_name = T,
           show_rownames = T), width = 12, height = 12)

for (m in unique(all_markers$cluster)) {
  lung_genes <- all_markers %>%
    subset(cluster==m & p_val_adj<0.1) 
  marker_genes <- row.names(subset(fData(HSMM), gene_short_name %in% head(lung_genes$gene,50)))
  diff_test_res <- differentialGeneTest(HSMM[marker_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
  pseudotime_heatmap <- plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                                                num_clusters = 3,
                                                return_heatmap = TRUE,
                                                cluster_rows = TRUE,
                                                show_rownames = T)
  ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/', m, '_pseudotime_heatmap.png'), plot = pseudotime_heatmap, width = 6, height = 6)
  ggsave(paste0(rslt, '/', Sample_Name, '_Trajectories/', m, '_pseudotime_heatmap.pdf'), plot = pseudotime_heatmap, width = 6, height = 6)
}
saveRDS(object = HSMM,file = paste0(rslt,'/',Sample_Name,'_Trajectories/',Sample_Name, '.rds'))



