# Packages Library
if(! require(Seurat)) install.packages("Seurat", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(scater)) install.packages("scater", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(BiocParallel)) install.packages("BiocParallel", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(BiocNeighbors)) install.packages("BiocNeighbors", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(ggplot2)) install.packages("ggplot2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(reshape2)) install.packages("reshape2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(dplyr)) install.packages("dplyr", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(stringr)) install.packages("stringr", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(scales)) install.packages("scales", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(easyGgplot2)) install.packages("easyGgplot2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(openxlsx)) install.packages("openxlsx", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(xml2)) install.packages("xml2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(rvest)) install.packages("rvest", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)

option_list <- list(
  make_option(c("-C", "--CRdir"), type = "character", default = ".", action = "store", help = "The directory where the CellRanger in 2.summary and 4.1.Analysis_CellRanger results store"
  ),
  make_option(c("-S", "--sample"), type = "character", default = FALSE, action = "store", help = "Sample name"
  ),
  make_option(c("-x", "--x.low.cutoff"), type = "double", default = 0.125, action = "store", help = "The x.low.cutoff value"
  ),
  make_option(c("-X", "--x.high.cutoff"), type = "double", default = 5, action = "store", help = "The x.high.cutoff value"
  ),
  make_option(c("-y", "--y.cutoff"), type = "double", default = 0.1, action = "store", help = "The y.cutoff value"
  )
)
opt = parse_args(OptionParser(option_list = option_list, usage = "4_2-Analysis_Seurat.R --CRdir <CellRanger results path> --sample <Sample name>"))


#library(scran)
#library(monocle)

#单个样本分析

rslt <- opt$CRdir
Sample_Name <- opt$sample

x.low.cutoff <- opt$x.low.cutoff
x.high.cutoff <- opt$x.high.cutoff
y.cutoff <- opt$y.cutoff

###单个样本分析

if(!dir.exists(paste0(rslt, '/', Sample_Name, '_seurat'))){dir.create(paste0(rslt, '/', Sample_Name, '_seurat'))}
sample3 <- readRDS(file=paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'.rds'))

###PCA
if(!dir.exists(paste0(rslt, '/', Sample_Name, '_seurat/PCA'))){dir.create(paste0(rslt, '/', Sample_Name, '_seurat/PCA'))}
sample3 <- FindVariableFeatures(sample3, selection.method = "vst", nfeatures = 99999999)
high_varience_gene <- Seurat::VariableFeatures(sample3, max=10000)

cluster_anly <- paste0(rslt,'/',Sample_Name,'_reanalysis/analysis')
exp_dat <- read.csv(paste0(cluster_anly, '/diffexp/graphclust', '/differential_expression.csv'))
colnames(exp_dat)[1:2] <- c('Gene_ID',	'Gene_name')
high_varience_gene_dat <- subset(exp_dat[,c(1:2)], exp_dat$Gene_name %in% high_varience_gene)
write.csv(high_varience_gene_dat, paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_high_varience_gene.csv'), row.names = FALSE)


high_varience_gene_number <- length(high_varience_gene)
high_varience_gene_info <- cbind(x.low.cutoff, x.high.cutoff, y.cutoff, high_varience_gene_number)
write.csv(high_varience_gene_info, paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_high_varience_gene_information.csv'), row.names = FALSE)

PCA_embedding <- as.data.frame(sample3@reductions$pca@cell.embeddings)
PCA_embedding <- as.data.frame(cbind(Barcode=row.names(PCA_embedding), PCA_embedding))
write.csv(PCA_embedding, paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_pca_embeddings.csv'), row.names = FALSE)

PCA_loading <- as.data.frame(sample3@reductions$pca@feature.loadings)
PCA_loading <- as.data.frame(cbind(Features=row.names(PCA_loading), PCA_loading))
write.csv(PCA_loading, paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_pca_loadings.csv'), row.names = FALSE)

pca_stdev <- as.data.frame(sample3@reductions$pca@stdev)
PCA_sdev <- as.data.frame(cbind(PC=row.names(pca_stdev), stdev=pca_stdev[,1]))
write.csv(PCA_sdev, paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_pca_sdev.csv'), row.names = FALSE)

top10 <- head(high_varience_gene, 10)
plot1 <- VariableFeaturePlot(sample3) + RestoreLegend(position = "bottom")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + NoLegend()
disp_p <- plot1 + plot2
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_Dispersion.png'), 
       plot = disp_p, device = 'png', height = 200, width = 300, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_Dispersion.pdf'), 
       plot = disp_p, device = 'pdf', height = 200, width = 300, limitsize = FALSE, units = 'mm')

elo_p <- ElbowPlot(sample3, ndims = 50)
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_Elbow.png'), 
       plot = elo_p, device = 'png', height = 160, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_Elbow.pdf'), 
       plot = elo_p, device = 'pdf', height = 160, width = 200, limitsize = FALSE, units = 'mm')


elo_dat <- elo_p$data
library(ggpmisc)
library(ggtech)
elo_p50 <- ggplot(elo_dat, aes(dims, stdev)) + geom_point(size = 1)
elo_p50 <- elo_p50 + geom_smooth(formula = y ~ x, method = "loess", span=1, se=FALSE) 
elo_p50 <- elo_p50 + ylab('Standard Deviation of PC') + xlab('PC') 
elo_p50 <- elo_p50 + theme_bw()
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_PCA_sdev_fitted.png'), 
       plot = elo_p50, device = 'png', height = 160, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_PCA_sdev_fitted.pdf'), 
       plot = elo_p50, device = 'pdf', height = 160, width = 200, limitsize = FALSE, units = 'mm')


ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_PCA_Heatmap.png'), 
       plot = DimHeatmap(sample3, dims = 1:15, cells = 500, balanced = TRUE), device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_PCA_Heatmap.pdf'), 
       plot = DimHeatmap(sample3, dims = 1:15, cells = 500, balanced = TRUE), device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

pca_dmp <- DimPlot(sample3, reduction = "pca")
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_PCA.png'), 
       plot = pca_dmp, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/PCA/', Sample_Name, '_PCA.pdf'), 
       plot = pca_dmp, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

###Marker
if(!dir.exists(paste0(rslt, '/', Sample_Name, '_seurat/Marker'))){dir.create(paste0(rslt, '/', Sample_Name, '_seurat/Marker'))}
Avg_Exp <- data.frame(Seurat::AverageExpression(sample3))
write.csv(Avg_Exp, paste0(rslt, '/', Sample_Name, '_seurat/Marker/AverageExpression.xls'))

avg_cor <- as.data.frame(cor(Avg_Exp))
avg_corp <- pheatmap::pheatmap(avg_cor)
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_AverageExp.png'), 
       plot = avg_corp, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_AverageExp.pdf'), 
       plot = avg_corp, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

all_markers <- FindAllMarkers(object = sample3, test.use = 'wilcox', slot = 'counts', verbose = TRUE)
if (nrow(all_markers)>1){
  top20 <- all_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  hmp20 <- DoHeatmap(sample3, features = top20$gene) 
  ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Heatmap.png'), 
         plot = hmp20, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
  ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Heatmap.pdf'), 
         plot = hmp20, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')
}

annotation_col <- data.frame(seurat_clusters=sample3@meta.data$seurat_clusters)
rownames(annotation_col)<-row.names(sample3@meta.data)
cls <- stringr::str_replace(colnames(exp_dat)[3:ncol(exp_dat)], pattern = ".Mean.Counts|.Log2.fold.change|.Adjusted.p.value", '')
cls_all <- list()
cls <- unique(cls)
for(s in cls){
  cls_nm <- c(colnames(exp_dat)[1:2], colnames(exp_dat)[startsWith(colnames(exp_dat), paste0(s,'.'))])
  cls_nodat <- exp_dat[,cls_nm]
  cls_nodat <- subset(cls_nodat, cls_nodat[,5]<0.05 & abs(cls_nodat[,4])>0 & cls_nodat[,3]>1)  ###这里改动一下，将筛选条件从logFC>0，改成|logFC|>0
  colnames(cls_nodat) <-  c('Gene_ID',	'Gene_name', 'Mean.Counts','Log2.fold.change','Adjusted.p.value')
  if(nrow(cls_nodat)>1){
    cls_nodat$Log2.fold.change <- as.numeric(cls_nodat$Log2.fold.change)
    cls_dat <- cls_nodat %>% arrange(-Log2.fold.change) %>% top_n(n=4, wt=Log2.fold.change)
    cls_nodat$Cluster <- s
    cls_all[[s]] <- cls_nodat
  }
}
cls_all <- as.data.frame(bind_rows(cls_all))
mker_gene <- cls_all %>% group_by(Cluster) %>% top_n(n=3, wt=Log2.fold.change)
chosen <- mker_gene$Gene_name
assays_chosen <- as.data.frame(sample3@assays[["RNA"]]@counts)[chosen,]
assays_chosen <- assays_chosen[!startsWith(row.names(assays_chosen),'NA'),]
plot_phmp <- pheatmap::pheatmap(assays_chosen, 
                                show_colnames=FALSE, na.rm = TRUE,
                                annotation_col = annotation_col)
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_plotHeatmap.png'), 
       plot = plot_phmp, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_plotHeatmap.pdf'), 
       plot = plot_phmp, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

SingCellExp <- as.SingleCellExperiment(sample3)
exprs_mat <- SummarizedExperiment::assay(SingCellExp, "counts", withDimnames=FALSE)
ave_exprs <- matrixStats::rowSums2(as.matrix(exprs_mat))
chon <- order(ave_exprs, decreasing=TRUE)
chosen <- head(chon, 50)
sub_ave <- ave_exprs[chosen]
total_exprs <- sum(ave_exprs)
pct <- scales::percent(sum(sub_ave) / total_exprs, accuracy=0.1)
HighestExpPlt <- scater::plotHighestExprs(SingCellExp, exprs_values = "counts")+labs(x="% of total counts",y="Featrue",title=paste0("Top 50 account for ", pct, " of total"))
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_plotHighestExprs.png'), 
       plot = HighestExpPlt, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_plotHighestExprs.pdf'), 
       plot = HighestExpPlt, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

if (nrow(all_markers)>1){
  vln1 <- all_markers %>% dplyr::group_by(cluster) %>% top_n(n=4, wt=avg_log2FC)
  for (n in unique(all_markers$cluster)) {
    vln <- subset(vln1, vln1$cluster==n)
    Cluster_coExpression <- FeaturePlot(sample3, features = unlist(vln[c(1:2),'gene']), blend = TRUE)
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_coExpression.png'), 
           plot = Cluster_coExpression, device = 'png', height = 200, width = 250, limitsize = FALSE, units = 'mm')
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_coExpression.pdf'), 
           plot = Cluster_coExpression, device = 'pdf', height = 200, width = 250, limitsize = FALSE, units = 'mm')
    
    seclu_tsne <- FeaturePlot(sample3, features = vln$gene, reduction = 'tsne', ncol = 2)
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_tsne.png'), 
           plot = seclu_tsne, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_tsne.pdf'), 
           plot = seclu_tsne, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')
    
    seclu_umap <- FeaturePlot(sample3, features = vln$gene, reduction = 'umap', ncol = 2)
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_umap.png'), 
           plot = seclu_umap, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_umap.pdf'), 
           plot = seclu_umap, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')
    
    seclu_vln <- VlnPlot(sample3, features = vln$gene, pt.size = 0, ncol = 2)
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_violin.png'), 
           plot = seclu_vln, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
    ggsave(paste0(rslt, '/', Sample_Name, '_seurat/Marker/', Sample_Name, '_Cluster_', n, '_violin.pdf'), 
           plot = seclu_vln, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')
  }
}

###DIMENSION 
if(!dir.exists(paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION'))){dir.create(paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION'))}
dimention_cluster <- cbind(Barcode=row.names(sample3@meta.data), 
                           Cluster=sample3@meta.data[,'seurat_clusters'])
write.csv(dimention_cluster, paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_cluster.csv'), row.names = FALSE)

UMAP <- sample3@reductions[["umap"]]@cell.embeddings
UMAP <- as.data.frame(cbind(Barcode=row.names(UMAP),UMAP))
write.csv(UMAP, paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_UMAP.csv'), row.names = FALSE)

TSNE <- sample3@reductions[["tsne"]]@cell.embeddings
TSNE <- as.data.frame(cbind(Barcode=row.names(TSNE),TSNE))
write.csv(TSNE, paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_tSNE.csv'), row.names = FALSE)


dim_umap <- DimPlot(sample3, reduction = "umap")
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_cluster_UMAP.png'), 
       plot = dim_umap, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_cluster_UMAP.pdf'), 
       plot = dim_umap, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

dim_tsnep <- DimPlot(sample3, reduction = "tsne")
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_cluster_tSNE.png'), 
       plot = dim_tsnep, device = 'png', height = 200, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt, '/', Sample_Name, '_seurat/DIMENSION/', Sample_Name, '_cluster_tSNE.pdf'), 
       plot = dim_tsnep, device = 'pdf', height = 200, width = 200, limitsize = FALSE, units = 'mm')

###DIFF
if(!dir.exists(paste0(rslt, '/', Sample_Name, '_seurat/DIFF'))){dir.create(paste0(rslt, '/', Sample_Name, '_seurat/DIFF'))}
gzfile <- paste0(rslt,"/",Sample_Name,"_cellranger/raw_feature_bc_matrix/features.tsv.gz")
gf <- gzfile(gzfile,'r')
mt <- read.table(gf, header = F, sep = '\t')
close(gf)
if (nrow(all_markers)>1){
  Cluster_diff_sig <- all_markers
  Cluster_diff_sig <- left_join(Cluster_diff_sig, mt, by=c('gene'='V2'))
  Cluster_diff_sig <- Cluster_diff_sig[,c(8,7,6,2,1,5,3,4)]
  colnames(Cluster_diff_sig)[1:2] <- c('Gene', 'Gene_name')
  write.csv(Cluster_diff_sig, paste0(rslt, '/', Sample_Name, '_seurat/DIFF/Cluster_diff.csv'), row.names = FALSE)
  
  for (d in unique(Cluster_diff_sig$cluster)) {
    sub_cluster_diff <- subset(Cluster_diff_sig, Cluster_diff_sig$cluster==d)
    write.csv(sub_cluster_diff, paste0(rslt, '/', Sample_Name, '_seurat/DIFF/Cluster_', d, '_diff.csv'), row.names = FALSE)
    write.csv(subset(sub_cluster_diff, sub_cluster_diff$p_val_adj<0.05), paste0(rslt, '/', Sample_Name, '_seurat/DIFF/Cluster_', d, '_diff_significant.csv'), row.names = FALSE)
  }
}

###保存Seurat对象
saveRDS(object = sample3, file = paste0(rslt,'/',Sample_Name,'_seurat/',Sample_Name,".rds"))