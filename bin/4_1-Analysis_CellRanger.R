if(! require(reshape2)) install.packages("reshape2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(ggalluvial)) install.packages("ggalluvial", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(cowplot)) install.packages("cowplot", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(tidyverse)) install.packages("tidyverse", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(ggpubr)) install.packages("ggpubr", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(ggthemes)) install.packages("ggthemes", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(Seurat)) install.packages("Seurat", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)

option_list <- list(
  make_option(c("-C", "--Cell_filtering"), type = "character", default = ".", action = "store", help = "The directory where the Cell filtering results store"
  ),
  make_option(c("-R", "--reanalysis"), type = "character", default = ".", action = "store", help = "The directory where the Cell filtering results store"
  ),
  make_option(c("-S", "--sample"), type = "character", default = FALSE, action = "store", help = "The directory where the CellRanger results store"
  )
)
opt <-  parse_args(OptionParser(option_list = option_list, usage = "4_1-Analysis_CellRanger --Cell_filtering . --reanalysis . --sample <Sample name>"))


####4.1.Analysis_CellRanger
src <- opt$Cell_filtering  
rslt <- opt$reanalysis
Sample_Name <- opt$sample


###单个样本
if(!dir.exists(paste0(rslt, '/', Sample_Name, '_reanalysis'))) {dir.create(paste0(rslt, '/', Sample_Name, '_reanalysis'), recursive = TRUE)}
load(paste0(src, '/',Sample_Name,'_cellfilter/',Sample_Name,'.Rdata'))
#cellranger_reanalyze <- paste0(src, '/reanalyze/', Sample_Name, '_reanalyze/outs')

##graphclust
cluster_f <- paste0(rslt, '/', Sample_Name, '_cellranger/analysis/clustering/')
clu_file <- dir(cluster_f, '.csv$', recursive = TRUE, full.names = TRUE)
clu_nm <- stringr::str_split(clu_file, '/', simplify = TRUE)
clu_nm <- clu_nm[,ncol(clu_nm)-1]
clu_ls <- list()
for (c in 1:length(clu_nm)) {
    clu_df <- read.csv(clu_file[c])
    clu_df$group <- clu_nm[c]
    clu_ls[[clu_nm[c]]] <- clu_df
  }
clu_dat <- as.data.frame(dplyr::bind_rows(clu_ls))
clu_dat$Cluster <- as.factor(clu_dat$Cluster)
#levels(clu_dat$Cluster) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
nm <- c("kmeans_2_clusters", "kmeans_3_clusters","kmeans_4_clusters","kmeans_5_clusters",
        "kmeans_6_clusters","kmeans_7_clusters","kmeans_8_clusters","kmeans_9_clusters",
        "kmeans_10_clusters","graphclust")
x <- c('a','b','c','d','e','f','g','h','i','j')
lab_nm <- data.frame(nm, x)
clu_dat <- dplyr::left_join(clu_dat, lab_nm, by=c('group'='nm'))
ggalluvial <- ggplot(clu_dat, aes(x = factor(x), stratum = Cluster, alluvium =Barcode,
                                  y = 1, fill = Cluster, label = Cluster)) +
  scale_x_discrete(labels = nm, expand = c(0.05, 0.05)) +
  geom_flow() +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  theme_cowplot() + xlab('') + ylab('') +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45), 
        panel.grid.major =element_blank(), panel.grid.minor = element_blank())
ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/clustering/graphclust/',Sample_Name,'_ggalluvial.png'), plot = ggalluvial, device = 'png', height = 200, width = 400, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/clustering/graphclust/',Sample_Name,'_ggalluvial.pdf'), plot = ggalluvial, device = 'pdf', height = 200, width = 400, limitsize = FALSE, units = 'mm')

###ALL maps
##TSNE
tsne_dat <- read.csv(paste0(rslt,'/', Sample_Name, '_cellranger/analysis/tsne/2_components/projection.csv'))
tsne_dat <- dplyr::left_join(clu_dat, tsne_dat, by=c('Barcode'='Barcode'))
df2 <- subset(tsne_dat, !tsne_dat$group=='graphclust')
lvls <- c("kmeans_2_clusters","kmeans_3_clusters", "kmeans_4_clusters","kmeans_5_clusters","kmeans_6_clusters","kmeans_7_clusters","kmeans_8_clusters", "kmeans_9_clusters","kmeans_10_clusters")
df2$group <- factor(df2$group, lvls)
pltb2 <- ggscatter(df2, x = 'TSNE.1', y = 'TSNE.2', color = 'Cluster',size = 0.3) + theme_few(base_size = 14) + facet_wrap(~group)
ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/tsne/',Sample_Name,'_allkmeans_clustering.png'), plot = pltb2, device = 'png', height = 200, width = 250, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/tsne/',Sample_Name,'_allkmeans_clustering.pdf'), plot = pltb2, device = 'pdf', height = 200, width = 250, limitsize = FALSE, units = 'mm')

##UMAP
umap_dat <- read.csv(paste0(rslt,'/', Sample_Name, '_cellranger/analysis/umap/2_components/projection.csv'))
umap_dat <- dplyr::left_join(clu_dat, umap_dat, by=c('Barcode'='Barcode'))
umapdf2 <- subset(umap_dat, !umap_dat$group=='graphclust')
lvls <- c("kmeans_2_clusters","kmeans_3_clusters", "kmeans_4_clusters","kmeans_5_clusters","kmeans_6_clusters","kmeans_7_clusters","kmeans_8_clusters", "kmeans_9_clusters","kmeans_10_clusters")
umapdf2$group <- factor(umapdf2$group, lvls)
umappltb2 <- ggscatter(umapdf2, x = 'UMAP.1', y = 'UMAP.2', color = 'Cluster',size = 0.3) + theme_few(base_size = 14) + facet_wrap(~group)
ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/umap/',Sample_Name,'_allkmeans_clustering_umap.png'), plot = umappltb2, device = 'png', height = 200, width = 250, limitsize = FALSE, units = 'mm')
ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/umap/',Sample_Name,'_allkmeans_clustering_umap.pdf'), plot = umappltb2, device = 'pdf', height = 200, width = 250, limitsize = FALSE, units = 'mm')

##-------##
cluster_anly <- paste0(rslt,'/',Sample_Name,'_cellranger/analysis')
diffexps <- dir(paste0(cluster_anly, '/diffexp'))
for(d in diffexps){
  ### clustering  
  cluster_dat <- read.csv(paste0(cluster_f, '/',d, '/clusters.csv'))
  tsne_dat <- read.csv(paste0(rslt,'/', Sample_Name, '_cellranger/analysis/tsne/2_components/projection.csv'))
  tsne_dat_sub <- subset(tsne_dat, tsne_dat$Barcode %in% cluster_dat$Barcode)
  tsne_dat_sub <- left_join(tsne_dat_sub, cluster_dat, by=c('Barcode'='Barcode'))
  tsne_dat_sub$Cluster <- factor(tsne_dat_sub$Cluster)
  pltb1 <- ggscatter(tsne_dat_sub, x = 'TSNE.1', y = 'TSNE.2', color = 'Cluster', size = 1.5, 
                     title = unique(tsne_dat_sub$group)) + theme_few(base_size = 14)
  
  if(!dir.exists(paste0(rslt,'/', Sample_Name, '_reanalysis/analysis/tsne/',d))) dir.create(paste0(rslt,'/', Sample_Name, '_reanalysis/analysis/tsne/',d))
  cnm <- stringr::str_replace(d,'_clusters','')
  cnm <- stringr::str_replace(cnm,'graphclust','graphbased')
  ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/tsne/',d,'/',Sample_Name,'_',cnm,'_clustering.png'), plot = pltb1, device = 'png', height = 350, width = 400, limitsize = FALSE, units = 'mm')
  ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/tsne/',d,'/',Sample_Name,'_',cnm,'_clustering.pdf'), plot = pltb1, device = 'pdf', height = 350, width = 400, limitsize = FALSE, units = 'mm')
  
  ## UMAP
  umap_dat <- read.csv(paste0(rslt,'/', Sample_Name, '_cellranger/analysis/umap/2_components/projection.csv'))
  umap_dat_sub <- subset(umap_dat, umap_dat$Barcode %in% cluster_dat$Barcode)
  umap_dat_sub <- left_join(umap_dat_sub, cluster_dat, by=c('Barcode'='Barcode'))
  umap_dat_sub$Cluster <- factor(umap_dat_sub$Cluster)
  pltb1 <- ggscatter(umap_dat_sub, x = 'UMAP.1', y = 'UMAP.2', color = 'Cluster', size = 1.5, 
                     title = unique(umap_dat_sub$group)) + theme_few(base_size = 14)
  
  if(!dir.exists(paste0(rslt,'/', Sample_Name, '_reanalysis/analysis/umap/',d))) dir.create(paste0(rslt,'/', Sample_Name, '_reanalysis/analysis/umap/',d))
  ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/umap/',d,'/',Sample_Name,'_',cnm,'_clustering.png'), plot = pltb1, device = 'png', height = 350, width = 400, limitsize = FALSE, units = 'mm')
  ggsave(paste0(rslt,'/',Sample_Name,'_reanalysis/analysis/umap/',d,'/',Sample_Name,'_',cnm,'_clustering.pdf'), plot = pltb1, device = 'pdf', height = 350, width = 400, limitsize = FALSE, units = 'mm')
  
  
  
  ### diffexp  
  exp_dat <- read.csv(paste0(cluster_anly, '/diffexp/', d, '/differential_expression.csv'))
  colnames(exp_dat)[1:2] <- c('Gene_ID',	'Gene_name')
  cls <- stringr::str_replace(colnames(exp_dat)[3:ncol(exp_dat)], pattern = ".Mean.Counts|.Log2.fold.change|.Adjusted.p.value", '')
  cls_all <- list()
  cls <- unique(cls)
  for(s in cls){
      cls_nm <- c(colnames(exp_dat)[1:2], colnames(exp_dat)[startsWith(colnames(exp_dat), paste0(s,'.'))])
      cls_nodat <- exp_dat[,cls_nm]
      cls_nodat <- subset(cls_nodat, cls_nodat[,5]<0.05 & abs(cls_nodat[,4])>0 & cls_nodat[,3]>1)  ###这里改动一下，将筛选条件从logFC>0，改成|logFC|>0
      write.csv(subset(cls_nodat, cls_nodat[,5]<0.05), paste0(cluster_anly, '/diffexp/', d, '/', stringr::str_replace(s,'\\.','_'), '.csv'), row.names = FALSE)
      colnames(cls_nodat) <-  c('Gene_ID',	'Gene_name', 'Mean.Counts','Log2.fold.change','Adjusted.p.value')
      
      if(nrow(cls_nodat)>1){
        cls_nodat$Log2.fold.change <- as.numeric(cls_nodat$Log2.fold.change)
        cls_dat <- cls_nodat %>% arrange(-Log2.fold.change) %>% top_n(n=4, wt=Log2.fold.change)
        vln_p <- VlnPlot(sample3, features = cls_dat$`Gene_name`, pt.size = 0, ncol = 2) 
        ggsave(paste0(cluster_anly, '/diffexp/', d, '/', stringr::str_replace(s,'\\.','_'), '_gene_violin.png'), 
               plot = vln_p, device = 'png', height = 200, width = 300, limitsize = FALSE, units = 'mm')
        ggsave(paste0(cluster_anly, '/diffexp/', d, '/', stringr::str_replace(s,'\\.','_'), '_gene_violin.pdf'), 
               plot = vln_p, device = 'pdf', height = 200, width = 300, limitsize = FALSE, units = 'mm')
        
        cls_nodat$Cluster <- s
        cls_all[[s]] <- cls_nodat
      }
    }
  cls_all <- as.data.frame(bind_rows(cls_all))
  mker_gene <- cls_all %>% group_by(Cluster) %>% top_n(n=3, wt=Log2.fold.change)
  all_markers <- FindAllMarkers(object = sample3, test.use = 'DESeq2', slot = 'counts', verbose = TRUE)###放到d的循环前，可以节省时间
  if (nrow(all_markers)>1){
      top_gn <- dplyr::arrange(all_markers, all_markers$avg_logFC) %>% head(9)
      ftr_tsne <- FeaturePlot(sample3, features = top_gn$gene, reduction = 'tsne', ncol = 3)
      ggsave(paste0(cluster_anly, '/diffexp/', d, '/', Sample_Name, '_gene.tsne.png'), 
             plot = ftr_tsne, device = 'png', height = 250, width = 300, limitsize = FALSE, units = 'mm')
      ggsave(paste0(cluster_anly, '/diffexp/', d, '/', Sample_Name, '_gene.tsne.pdf'), 
             plot = ftr_tsne, device = 'pdf', height = 250, width = 300, limitsize = FALSE, units = 'mm')
      ftr_umap <- FeaturePlot(sample3, features = top_gn$gene, reduction = 'umap', ncol = 3)
      ggsave(paste0(cluster_anly, '/diffexp/', d, '/', Sample_Name, '_gene.umap.png'), 
             plot = ftr_umap, device = 'png', height = 250, width = 300, limitsize = FALSE, units = 'mm')
      ggsave(paste0(cluster_anly, '/diffexp/', d, '/', Sample_Name, '_gene.umap.pdf'), 
             plot = ftr_umap, device = 'pdf', height = 250, width = 300, limitsize = FALSE, units = 'mm')
    }
  chosen <- mker_gene$Gene_name
  rna_assays <- as.data.frame(sample3@assays$RNA@counts)[chosen,]
  annotation_col <- as.data.frame(cbind(Cluster=cluster_dat$Cluster))
  row.names(annotation_col) <- cluster_dat$Barcode
  annotation_col <- annotation_col %>% arrange(Cluster)
  anno_col <- as.data.frame(table(annotation_col$Cluster))
  annotation_col <- subset(annotation_col, annotation_col$Cluster %in% as.character(anno_col[which(anno_col$Freq>=10),1]))
  annotation_row <- mker_gene[,c('Gene_name','Cluster')]
  annotation_row$Cluster <- stringr::str_replace(annotation_row$Cluster,'Cluster.','')
  annotation_row <- annotation_row %>% arrange(Gene_name) %>% group_by(Gene_name) %>% top_n(n=1, wt=row_number()) %>% arrange(Cluster)
  annotation_row <- subset(annotation_row, annotation_row$Cluster %in% anno_col[which(anno_col$Freq>=10),1])
  gene_names <- annotation_row$Gene_name
  annotation_row <- as.data.frame(cbind(Cluster=annotation_row$Cluster))
  row.names(annotation_row) <- gene_names
  rna_assays <- subset(rna_assays, row.names(rna_assays) %in% row.names(annotation_row))
  annotation_col <- subset(annotation_col, row.names(annotation_col) %in% colnames(rna_assays))
  annotation_col <- subset(annotation_col, annotation_col$Cluster %in% annotation_row$Cluster)
  annotation_col$Cluster <- as.factor(annotation_col$Cluster)
  annotation_row$Cluster <- as.factor(annotation_row$Cluster)
  rna_assays <- rna_assays[,row.names(annotation_col)]
  
  p_htmp <- pheatmap::pheatmap(scale(rna_assays), 
                               show_colnames=FALSE, cluster_cols = F, cluster_rows = F,
                               annotation_col = annotation_col, annotation_row = annotation_row)
  ggsave(paste0(cluster_anly, '/diffexp/', d, '/', Sample_Name, '_',d,'_heatmap.png'), 
         plot = p_htmp, device = 'png', height = 300, width = 300, limitsize = FALSE, units = 'mm')
  ggsave(paste0(cluster_anly, '/diffexp/', d, '/', Sample_Name, '_',d,'_heatmap.pdf'), 
         plot = p_htmp, device = 'pdf', height = 300, width = 300, limitsize = FALSE, units = 'mm')
}





