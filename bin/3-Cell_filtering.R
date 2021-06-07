if(! require(Seurat)) install.packages("Seurat", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(ggplot2)) install.packages("ggplot2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)

option_list <- list(
  make_option(c("-C", "--CRdir"), type = "character", default = ".", action = "store", help = "The directory where the CellRanger results store"
  ),
  make_option(c("-S", "--sample"), type = "character", default = FALSE, action = "store", help = "The directory where the CellRanger results store"
  )
)
opt <-  parse_args(OptionParser(option_list = option_list, usage = "3-Cell_filtering --CRdir <CellRanger results path> --sample <Sample name>"))


########## All the paht
src <- opt$CRdir  # read source
rslt <- opt$CRdir # output path
Sample_Name <- opt$sample

mtr_summ <- read.csv(paste0(src, "/",Sample_Name,"_cellranger/metrics_summary.csv"))
mtr_summ <- as.data.frame(cbind(Sample_Name, mtr_summ))
#rmarkdown::paged_table(mtr_summ)


##单个样本

if(!dir.exists(paste0(rslt,'/',Sample_Name,'_cellfilter')))(dir.create(paste0(rslt,'/',Sample_Name,'_cellfilter'),recursive = TRUE))
raw.data <- Read10X(data.dir = paste0(src, "/",Sample_Name,"_cellranger/filtered_feature_bc_matrix"))
sample <- CreateSeuratObject(raw.data, project = Sample_Name, min.cells = 3)  # min.cells = 3

###写gene_bar文件
gene_bar <- as.data.frame(sample@assays[["RNA"]]@counts)
write.csv(gene_bar, paste0(rslt,'/',Sample_Name,'_gene_bar.csv'), row.names = FALSE)

#nGene_nUMI
nGene_nUMI <- sample@meta.data
nGene_nUMI <- as.data.frame(cbind(Barcode=row.names(nGene_nUMI), nGene_nUMI))
write.csv(nGene_nUMI, paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_nGene_nUMI.csv'), row.names = FALSE)

###样本过滤
sample[["percent.mito"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample[["percent.HB"]] <- PercentageFeatureSet(sample, pattern = "HB")
Vln_p <- VlnPlot(sample, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.HB"), ncol = 4)
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_nGene_nUMI_mito_HB.png'), plot = Vln_p, device = 'png', height = 150, width = 300, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_nGene_nUMI_mito_HB.pdf'), plot = Vln_p, device = 'pdf', height = 150, width = 300, limitsize = FALSE, units = 'mm')

##各个过滤指标相关性分析
plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
plot3 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.HB") + NoLegend()
relation_p <- plot1 + plot2 + plot3
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_nGene_nUMI_mito_HB_relation.png'), plot = relation_p, device = 'png', height = 150, width = 300, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_nGene_nUMI_mito_HB_relation.pdf'), plot = relation_p, device = 'pdf', height = 150, width = 300, limitsize = FALSE, units = 'mm')


##执行过滤
min.cells <- 3
low_nGene <- 200
high_nGene <- 2500
high_percent.mito <- 10
high_percent.HB <- 5
sample2 <- subset(sample, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mito<10 & percent.HB<5)
gene_filted <- as.data.frame(sample2@assays[["RNA"]]@counts)
gene_number <- raw.data@Dim[1]
cell_number <- mtr_summ$Estimated.Number.of.Cells
nGene_median <- mtr_summ$Median.Genes.per.Cell
gene_number_filtered <- sample2@assays$RNA@counts@Dim[1]
cell_number_filtered <- sample2@assays$RNA@counts@Dim[2]
nGene_median_filtered <- median(sample2@meta.data$nFeature_RNA)
filtered_information <- cbind.data.frame(Sample_Name, min.cells, low_nGene, high_nGene, high_percent.mito, high_percent.HB, gene_number, cell_number, nGene_median, gene_number_filtered, cell_number_filtered, nGene_median_filtered)
colnames(filtered_information) <- c('Sample_Name','min.cells','low.thresholds_nGene',
                                    'high.thresholds_nGene','high.thresholds_percent.mito',
                                    'high.thresholds_percent.HB','gene_number',
                                    'cell_number','nGene_median','gene_number_filtered',
                                    'cell_number_filtered','nGene_median_filtered'
)
write.csv(filtered_information, paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_filtered_information.csv'), row.names = FALSE)

##barcode_filtered
barcode_filtered <- as.data.frame(sample2@active.ident)
#colnames(barcode_filtered)[1] <- 'Barcode'
barcode_filtered <- as.data.frame(cbind(Barcode=row.names(barcode_filtered), barcode_filtered=barcode_filtered[,1]))
write.csv(barcode_filtered[,1], paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_barcode_filtered.csv'), row.names = FALSE)

###降维分群
#初始化
sample2 <- NormalizeData(sample2)
#找高变基因
x.low.cutoff <- 0.125
x.high.cutoff <- 5
y.cutoff <- 0.1
sample2 <- FindVariableFeatures(object = sample2, selection.method = "vst", mean.cutoff = c(x.low.cutoff,x.high.cutoff), dispersion.cutoff = c(y.cutoff,Inf))
#标准化
sample2 <- ScaleData(object = sample2, vars.to.regress = c("nUMI", "percent.mito"))
# PCA降维 
## 第一步，PCA分析
sample2 <- RunPCA(object = sample2, features = VariableFeatures(object = sample2))
## 第二步，评估最显著
sample2 <- JackStraw(object = sample2, num.replicate = 100)
sample2 <- ScoreJackStraw(sample2, dims = 1:20)
##comparing the distribution of p-values for each PC
#根据PCA结果找分群，要分两步，而Seurat V2版本在这里只有一步
sample2 <- FindNeighbors(sample2, reduction = "pca", dims = 1:10)
sample2 <- FindClusters(sample2, resolution = 0.4)

#tSNE降维
sample2 <- RunTSNE(sample2, dims.use = 1:10)
#UMAP降维
sample2 <- RunUMAP(sample2, dims = 1:10)
umap_p <- DimPlot(sample2, reduction = "umap") 
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_cluster_UMAP.png'), plot = umap_p, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_cluster_UMAP.pdf'), plot = umap_p, device = 'pdf', height = 150, width = 200, limitsize = FALSE, units = 'mm')

###检测Doublet
if(! require(DoubletFinder)) install.packages("DoubletFinder", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(parallel)) install.packages("parallel", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(spam)) install.packages("spam", repos="https://cloud.r-project.org/", dependencies = TRUE)

## pK Identification (no ground-truth) 
sweep.res.list_kidney <- paramSweep_v3(sample2, PCs = 1:10, sct = FALSE, num.cores = detectCores()-8)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
pk_plt <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric, group = 1)) + geom_line() + geom_point() 
pk_plt <- pk_plt + theme(plot.title=element_text(hjust=0.5),legend.title=element_text(face="bold")) 
pk_plt <- pk_plt + scale_fill_manual(values=c("red","green"))
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_pK_bcmvn.png'), plot = pk_plt, device = 'png', height = 150, width = 300, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_pK_bcmvn.pdf'), plot = pk_plt, device = 'pdf', height = 150, width = 300, limitsize = FALSE, units = 'mm')
write.csv(bcmvn_kidney, paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_pK_bcmvn.csv'), row.names = FALSE)

## Homotypic Doublet Proportion Estimate
annotations <- sample2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)    
#homotypic.prop <- 0.403973919489178
nExp_poi <- round(0.008*nrow(sample2@meta.data))  
#nExp_poi <- 9
## Assuming 0.8% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#nExp_poi.adj <- 5
mpK <- as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
#mpK <- 0.25
## Run DoubletFinder with varying classification stringencies 
#mypANN <- "DF.classifications_0.25_0.25_9"
#mypANN_adj <- " + labels(title = )"
sample3 <- doubletFinder_v3(sample2, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
meta_dbl <- as.data.frame(cbind(Barcode=row.names(sample3@meta.data), sample3@meta.data))
write.csv(meta_dbl, paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_meta_data_single_double.csv'), row.names = FALSE)
#nExp_poi <- 9
#sample3 <- doubletFinder_v3(sample3, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#detach("package:DoubletFinder", unload = TRUE)

###nGene_nUMI_filtered
meta_dbl2 <- as.data.frame(cbind(Barcode=row.names(sample3@meta.data), sample3@meta.data))
write.csv(meta_dbl2, paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_nGene_nUMI_filtered.csv'), row.names = FALSE)

###gene_filtered
gzfile <- paste0(src,"/",Sample_Name,"_cellranger/raw_feature_bc_matrix/features.tsv.gz")
gf <- gzfile(gzfile,'r')
mt <- read.table(gf, header = F, sep = '\t')
close(gf)
gene_filtered <- subset(mt[,c(1:2)], mt[,2] %in% row.names(sample3@assays$RNA@meta.features))
colnames(gene_filtered) <- c('Gene', 'Gene_name')
write.csv(gene_filtered, paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_gene_filtered.csv'), row.names = FALSE)

###single_doublet
single_doublt_plt <- DimPlot(sample3, reduction = "umap", group.by = paste0('DF.classifications_0.25_', mpK,'_',nExp_poi)) +labs(title=paste0('DF.classifications_0.25_', mpK,'_',nExp_poi))
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_cluster_UMAP_single_doublt.png'), plot = single_doublt_plt, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,'_cluster_UMAP_single_doublt.pdf'), plot = single_doublt_plt, device = 'pdf', height = 150, width = 200, limitsize = FALSE, units = 'mm')

###保存对象
saveRDS(object = sample3, file = paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,".rds"))
save(raw.data,sample,sample2,sample3, file = paste0(rslt,'/',Sample_Name,'_cellfilter/',Sample_Name,".Rdata"))

