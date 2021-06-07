rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(! require(clusterProfiler)) BiocManager::install("clusterProfiler", ask = FALSE, dependencies = TRUE)
if(! require(DOSE)) BiocManager::install("DOSE")
if(! require(pathview)) BiocManager::install("pathview")
if(! require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
if(! require(ggupset)) install.packages("ggupset", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(UpSetR)) install.packages("UpSetR", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(stringr)) install.packages("stringr", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(dplyr)) install.packages("dplyr", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(ggplot2)) install.packages("ggplot2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)

option_list <- list(
  make_option(c("-C", "--CRdir"), type = "character", default = ".", action = "store", help = "The directory where the Seurat results store"
  ),
  make_option(c("-S", "--sample"), type = "character", default = FALSE, action = "store", help = "Sample name"
  ),
  make_option(c("-O", "--Org"), type = "character", default = "org.Hs.eg.db", action = "store", help = "TF DB"
  ),
  make_option(c("-D", "--db"), type = "character", default = FALSE, action = "store", help = "TF DB"
  )
)
opt = parse_args(OptionParser(option_list = option_list, usage = "6-Enrichment.R --sample ${sample} <Sample name> --db <path of TF annotation>"))



enrichDotplot <- function(object, x = "geneRatio", color = "p.adjust",
                          showCategory=20, size=NULL, split = NULL,
                          font.size=12, title = "", orderBy="x", decreasing=TRUE) {
  
                                colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
                                if (x == "geneRatio" || x == "GeneRatio") {
                                  x <- "GeneRatio"
                                  if (is.null(size))
                                    size <- "Count"
                                } else if (x == "count" || x == "Count") {
                                  x <- "Count"
                                  if (is.null(size))
                                    size <- "GeneRatio"
                                } else if (is(x, "formula")) {
                                  x <- as.character(x)[2]
                                  if (is.null(size))
                                    size <- "Count"
                                } else {
                                  ## message("invalid x, setting to 'GeneRatio' by default")
                                  ## x <- "GeneRatio"
                                  ## size <- "Count"
                                  if (is.null(size))
                                    size  <- "Count"
                                }
                                
                                #nr <- unique(object@result[,x])
                                #if(length(nr)<20){
                                df <- object@result %>% arrange(p.adjust) %>% head(showCategory)
                                df$GeneRatio <- parse_ratio(df$GeneRatio)
                                #}else{
                                #  df <- fortify(object, showCategory = showCategory, split=split)
                                #}
                                ## already parsed in fortify
                                ## df$GeneRatio <- parse_ratio(df$GeneRatio)
                                
                                if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
                                  message('wrong orderBy parameter; set to default `orderBy = "x"`')
                                  orderBy <- "x"
                                }
                                
                                if (orderBy == "x") {
                                  df <- dplyr::mutate(df, x = eval(parse(text=x)))
                                }
                                
                                idx <- order(df[[orderBy]], decreasing = decreasing)
                                df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
                                ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
                                  geom_point() +
                                  scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
                                  ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
                                  ylab("Description") + ggtitle(title) + theme_dose(font.size) + scale_size(range=c(3, 8))
                                
}

#enrichDotplot(go, title = 'nm')



enrichBarplot <- function(height, x="Count", color='p.adjust', showCategory=20, font.size=12, title="", ...) {
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  object <- height
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  
  #nr <- unique(object@result[,x])
  #if(length(nr)<20){
  df <- object@result %>% arrange(p.adjust) %>% head(showCategory)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description)))
  #}else{
  #  df <- fortify(object, showCategory=showCategory, by=x, ...)
  #}
  
  
  if(colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
  } else {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = "Description")) +
      theme_dose(font.size) +
      theme(legend.position="none")
  }
  p + geom_bar(stat = "identity") + coord_flip() + ggtitle(title) + xlab("Description") + ylab("GeneRatio")
}


SC_strdb_network <- function(diff_genes, output_dir, species, nm, db){
  gzfile1 <- paste0(db, "/10090.protein.links.full.v11.0.txt.gz")
  gzfile2 <- paste0(db, "/9606.protein.links.full.v11.0.txt.gz")
  gzfile3 <- paste0(db, "/9986.protein.links.full.v11.0.txt.gz")
  gzfile4 <- paste0(db, "/7091.protein.links.full.v11.0.txt.gz")
  if(species==10090){
    gzfile <- gzfile1
  }else if(species==9606){
    gzfile <- gzfile2
  }else if(species==9986){
    gzfile <- gzfile3
  }else if(species==7091){
    gzfile <- gzfile4
  }
  gf <- gzfile(gzfile, 'rt')
  plinksfull <- read.table(gf, header = TRUE, quote = '', comment.char = '#')
  close(gf)
  
  if(length(diff_genes)>0){
    get_net <- function(ens_id){
      id <- paste0('%0d', ens_id)
      identifiers <- as.character(ens_id[1])
      for(i in 2:length(id)){
        identifiers <- paste0(identifiers, id[i])
      }
      
      required_score <- 400
      #species <- 10090
      stringdb_url <- paste0('https://string-db.org/api/tsv/network?identifiers=', identifiers, '&required_score=', required_score, '&species=', species)
      network_interac <- read.table(stringdb_url, header = TRUE)
      return(network_interac)
    }
    
    ends_id <- head(diff_genes,2100)
    m <- ceiling(length(ends_id)/100)
    if(m>2){
      aa <- seq(0, length(ends_id), length.out=m)
      aa <- round(aa)
      bb <- cut(1:length(ends_id), aa)
      cc <- unique(bb)
      
      sub_net <- list()
      p <- 1
      while(p<(m-1)) {
        for(n2 in cc[(p+1):(m-1)]){
          bb1 <- ends_id[bb==cc[p]]
          bb2 <- ends_id[bb==n2]
          ens_id <- c(bb1, bb2)
          get_res <- try(get_net(ens_id))
          if(inherits(get_res,"try-error") & p<(m-1)) { next }
          if(!class(get_res)=='try-error'){
            sub_net[[paste0(p,'+',n2)]] <- get_res
          }
          if(class(get_res)=='try-error'){print(paste0('error-', p,'+',which(cc==n2),'_',m-1))}else{print(paste0(p,'+',which(cc==n2),'_',m-1))}
        }
        p <- p + 1
      }
      network_interac <- dplyr::bind_rows(sub_net) %>% dplyr::distinct()
    }else{
      network_interac <- get_net(ends_id)
    }
    
    
    if(nrow(network_interac)>0){
      network_interac <- network_interac[,c("preferredName_A","preferredName_B","stringId_A","stringId_B","ncbiTaxonId","score",
                                            "nscore","fscore","pscore","ascore","escore","dscore","tscore")]
      network_interac$stringId_A <- paste0(network_interac$ncbiTaxonId,'.',network_interac$stringId_A)
      network_interac$stringId_B <- paste0(network_interac$ncbiTaxonId,'.',network_interac$stringId_B)
      
      homology <- plinksfull[,c('protein1','protein2','homology')]
      homology$homology <- homology$homology/100
      
      string_interactions <- dplyr::left_join(network_interac,homology,by=c('stringId_A'='protein1','stringId_B'='protein2'))
      string_interactions <- string_interactions[,c("preferredName_A",	"preferredName_B",	"stringId_A",	"stringId_B",	"nscore",	"fscore",	"pscore",	"homology",	"ascore",	"escore",	"dscore",	"tscore",	"score")]
      colnames(string_interactions) <- c("node1",	"node2",	"node1_external_id",	"node2_external_id",	"neighborhood_on_chromosome",	"gene_fusion",	"phylogenetic_cooccurrence",	"homology",	"coexpression",	"experimentally_determined_interaction",	"database_annotated",	"automated_textmining",	"combined_score")
      write.table(string_interactions, paste0(output_dir, '/', nm, '_ppi.tsv'), row.names = FALSE)
    }
  }
}













##### 6. Enrichment
rslt <- opt$CRdir  
Sample_Name <- opt$sample
#scr_path <- "/home/data/zgsdata/10X_scrna/bin"
dataset <- opt$Org
db <- opt$db

options(stringsAsFactors = FALSE)
grpcls_diff <- paste0(rslt, '/',Sample_Name,'_seurat/DIFF/')
csv_d <- dir(grpcls_diff, pattern = '_diff_significant.csv')
csv_d <- csv_d[startsWith(csv_d, 'Cluster_')]
csv_nm <- stringr::str_replace(csv_d, '_diff_significant.csv', '')

for (nm in csv_nm) {
  if (!is.null(nm)){
    diffexp_dat <- read.csv(paste0(grpcls_diff,'/', nm, '_diff_significant.csv'))
    colnames(diffexp_dat)[1] <- 'Gene_ID'
    gene_id <- diffexp_dat[,1]
    if(dataset=="org.Mm.eg.db"){require(org.Mm.eg.db)}
    if(dataset=="org.Hs.eg.db"){require(org.Hs.eg.db)}
    for (d in dev.list()) {
      dev.off(d)
    }
    
    #  GO
    if(!dir.exists(paste0(rslt, '/', Sample_Name, '_GO/', nm))){dir.create(paste0(rslt, '/', Sample_Name, '_GO/', nm), recursive = TRUE)}
    GOwd <- paste0(rslt, '/', Sample_Name, '_GO/', nm)
    eg <- bitr(gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb=dataset)  # id转换成ENTREZID
    genelist <- eg$ENTREZID
    go <- enrichGO(gene_id, OrgDb = dataset, ont='ALL',pAdjustMethod = 'BH',
                   #pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                   keyType = 'ENSEMBL')
    go_read <- enrichGO(eg$ENTREZID, OrgDb = dataset, ont='ALL',pAdjustMethod = 'BH', readable = TRUE,
                        #pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                        keyType = 'ENTREZID')
    go_genenm <- go_read@result[,c('ID','geneID')]
    colnames(go_genenm)[2] <- 'geneName'
    go_rslt <- left_join(go@result, go_genenm, by=c('ID'='ID'))
    go_rslt <- go_rslt[,c(1:9,11,10)]
    write.csv(go_rslt, paste0(GOwd, '/', nm, '_GOenrich.csv'), row.names = FALSE)
    write.csv(subset(go_rslt, go_rslt$p.adjust<0.05), paste0(GOwd, '/', nm, '_GOenrich_significant.csv'), row.names = FALSE)
    
    go@result <- subset(go@result, !go@result$ONTOLOGY=='NA')
    go@result$Description <- ifelse(nchar(go@result$Description)>55, paste0(substr(go@result$Description,1,55),'...'),go@result$Description)
    go_top10 <- go@result %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(-Count) %>%
      dplyr::mutate(row_number = row_number(-Count)) %>%
      subset(row_number<=20) %>%
      dplyr::arrange(ONTOLOGY,-Count)
    
    colors <- c(BP="#8DA1CB", CC="#FD8D62", MF="#66C3A5")
    GO_term_order <- factor(as.integer(rownames(go_top10)), labels = go_top10$Description)
    p <- ggplot(data=go_top10, aes(x=GO_term_order, y=Count, fill=ONTOLOGY))
    p <- p + geom_bar(stat='identity', width = 0.8) + theme_bw()
    p <- p + scale_colour_manual(values = colors, aesthetics = c("colour", "fill"))
    p <- p + xlab('GO term') + ylab('Num of Genes') + labs(title = paste0('The Most Enriched GO Terms\n', nm))
    p <- p + theme(axis.text.x = element_text(face = 'bold', color = 'gray50', angle = 70, vjust = 1, hjust = 1))
    goplt1 <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), plot.margin=unit(rep(3,4),'lines'))
    ggsave(paste0(nm, '_goBar_ONTO.png'), path = GOwd, plot = goplt1, limitsize = FALSE, width = 12, height = 8)
    ggsave(paste0(nm, '_goBar_ONTO.pdf'), path = GOwd,plot = goplt1, limitsize = FALSE, width = 12, height = 8)
    
    barplt1 <- enrichBarplot(go, title = paste0('The Most Enriched GO Terms\n',nm))
    ggsave(paste0(nm, '_GO_bar.png'), path = GOwd, plot = barplt1, limitsize = FALSE, width = 10, height = 8)
    ggsave(paste0(nm, '_GO_bar.pdf'), path = GOwd, plot = barplt1, limitsize = FALSE, width = 10, height = 8)
    
    goplt2 <- enrichDotplot(go, title = paste0('The Most Enriched GO Terms\n',nm))
    ggsave(paste0(nm, '_GO_dot.png'), path = GOwd, plot = goplt2, limitsize = FALSE, width = 10, height = 8)
    ggsave(paste0(nm, '_GO_dot.pdf'), path = GOwd, plot = goplt2, limitsize = FALSE, width = 10, height = 8)
    
    
    # KEGG 
    if(!dir.exists(paste0(rslt, '/', Sample_Name, '_KEGG/', nm, '/', nm))){dir.create(paste0(rslt, '/', Sample_Name, '_KEGG/', nm, '/', nm), recursive = TRUE)}
    KEGGwd <- paste0(rslt, '/', Sample_Name, '_KEGG/', nm, '/', nm)
    
    if(dataset=="org.Mm.eg.db") orgnm <- 'mmu'
    if(dataset=="org.Hs.eg.db") orgnm <- 'hsa'
    kegg <- enrichKEGG(eg$ENTREZID, organism = orgnm, keyType = 'kegg', pAdjustMethod = 'BH')
    keggrlst <- kegg@result
    keggnm_lst <- list()
    for (k in 1:nrow(keggrlst)) {
      gid <- keggrlst[k,1]
      gnid <-  keggrlst[k,8]
      gnid_spt <- stringr::str_split(gnid, pattern = '/', simplify = TRUE)[1,]
      esb_id <- bitr(gnid_spt, fromType = 'ENTREZID', toType = 'ENSEMBL', OrgDb = dataset)
      ENSEMBL_ID <- paste0(esb_id$ENSEMBL, collapse = '/')
      gnnm <- DOSE::EXTID2NAME(gnid_spt, keytype = 'ENTREZID', OrgDb = dataset)
      geneName <- paste0(gnnm, collapse = '/')
      keggID <- paste0(orgnm,':',paste0(gnid_spt, collapse = paste0('/',orgnm,':'), recycle0 = TRUE))
      kegg_gnnm <- as.data.frame(cbind(ID=gid, geneID=ENSEMBL_ID, geneName=geneName, keggID=keggID))
      keggnm_lst[[k]] <- kegg_gnnm
    }
    kegg_gnnm <- bind_rows(keggnm_lst)
    kegg_rlsts <- left_join(keggrlst[,-8], kegg_gnnm, by=c('ID'='ID'))
    kegg_rlsts <- kegg_rlsts[,c(1:7,9:11,8)]
    write.csv(kegg_rlsts, paste0(KEGGwd, '/', nm, '_KEGGenrich.csv'), row.names = FALSE)
    write.csv(subset(kegg_rlsts, kegg_rlsts$p.adjust<0.05), paste0(KEGGwd, '/', nm, '_KEGGenrich_significant.csv'), row.names = FALSE)
    
    keggplt1 <- enrichBarplot(kegg, title = paste0('The Most Enriched KEGG Terms\n',nm))
    ggsave(paste0(nm, '_kegg_bar.png'), path = KEGGwd, plot = keggplt1, limitsize = FALSE, width = 10, height = 8)
    ggsave(paste0(nm, '_kegg_bar.pdf'), path = KEGGwd, plot = keggplt1, limitsize = FALSE, width = 10, height = 8)
    
    keggplt2 <- enrichDotplot(kegg, title = paste0('The Most Enriched KEGG Terms\n',nm))
    ggsave(paste0(nm, '_kegg_dot.png'), path = KEGGwd, plot = keggplt2, limitsize = FALSE, width = 10, height = 8)
    ggsave(paste0(nm, '_kegg_dot.pdf'), path = KEGGwd, plot = keggplt2, limitsize = FALSE, width = 10, height = 8)
    
    
    ## kegg pathway
    degenes <- left_join(diffexp_dat[,c(1,4)], eg, by=c('Gene_ID'='ENSEMBL'))
    geneList <- arrange(degenes, desc(degenes[,2]))
    id_geneList <- geneList[,2]
    names(id_geneList) <- geneList$ENTREZID
    id_geneList <- sort(id_geneList, decreasing=TRUE)
    
    if(!dir.exists(paste0(rslt, '/', Sample_Name,'_KEGG/',nm, '/pathway_plot'))){dir.create(paste0(rslt, '/', Sample_Name,'_KEGG/',nm,'/pathway_plot'), recursive = TRUE)}
    keggpathway <- paste0(rslt, '/', Sample_Name, '_KEGG/',nm,'/pathway_plot/')
    path_id <- subset(keggrlst, keggrlst$p.adjust<0.05)[,1]
    for (pwid in path_id) {
      #pathview::download.kegg(pathway.id = '00190',kegg.dir = '/home/data/zgsdata/10X_scrna')
      
      out.pv <- try(pathview::pathview(gene.data = id_geneList, pathway.id = pwid, 
                             species = orgnm, cpd.idtype ="kegg", kegg.dir = keggpathway,
                             gene.idtype = "entrez", gene.annotpkg = NULL,min.nnodes = 3, kegg.native =TRUE)
                             , silent = TRUE)
      
      if(!'try-error' %in% class(out.pv)){
        if(class(out.pv$plot.data.gene)=="data.frame"){
        kegg_map <- unique(subset(out.pv$plot.data.gene, !out.pv$plot.data.gene$all.mapped=='')[,c(1:3)])
        colnames(kegg_map) <- c('pathway_gene_id',	'pathway_labels',	'all_mapped')
        trans_id <- bitr(kegg_map$pathway_gene_id, fromType = 'ENTREZID', toType = 'ENSEMBL', OrgDb = dataset)
        trans_id <- trans_id %>%
          group_by(ENTREZID) %>% 
          top_n(n=1, wt=row_number())
        kegg_map <- left_join(kegg_map, trans_id, by=c('pathway_gene_id'='ENTREZID'))
        for (m in 1:nrow(kegg_map)) {
          map_id <- stringr::str_split(kegg_map[m,3], ',', simplify = TRUE)[1,]
          map_etzid <- bitr(map_id, fromType = 'ENTREZID', toType = 'ENSEMBL', OrgDb = dataset)[,2]
          sig_dat <- diffexp_dat[which(diffexp_dat$Gene_ID %in% map_etzid), c(1,4)]
          sig_dat$sig <- case_when(sig_dat[,2] > 1 ~ 'up',
                                   sig_dat[,2] < -1 ~ 'down',
                                   TRUE ~ 'no')
          kegg_map$all_mapped_gene_id[m] <- paste0(paste0(sig_dat$Gene_ID, '(', sig_dat$sig,')-', map_id), collapse = ',')
        }
        kegg_maps <- left_join(kegg_map[,c(1,2,5,4)], diffexp_dat[,c(1,4)], by=c('ENSEMBL'='Gene_ID'))
        kegg_maps <- kegg_maps[,-4]
        write.csv(kegg_maps, paste0(keggpathway, '/', pwid, '.csv'), row.names = FALSE)
      }
      }
      Sys.sleep(1.5)
    }
    
    # ReactomePA
    if(!dir.exists(paste0(rslt, '/', Sample_Name, '_Reactome/', nm))){dir.create(paste0(rslt, '/', Sample_Name, '_Reactome/', nm), recursive = TRUE)}
    ReactomePAwd <- paste0(rslt, '/', Sample_Name, '_Reactome/', nm)
    
    if(dataset=="org.Mm.eg.db") org <- "mouse"
    if(dataset=="org.Hs.eg.db") org <- "human"
    reactome <- try(ReactomePA::enrichPathway(geneList$ENTREZID, readable=F,organism = org), silent = TRUE)
    if(!'try-error' %in% class(reactome) & nrow(reactome@result)>0){
      rctrlst <- reactome@result
      rctnm_lst <- list()
      for (k in 1:nrow(rctrlst)) {
        gid <- rctrlst[k,1]
        gnid <-  rctrlst[k,8]
        gnid_spt <- stringr::str_split(gnid, pattern = '/', simplify = TRUE)[1,]
        esb_id <- bitr(gnid_spt, fromType = 'ENTREZID', toType = 'ENSEMBL', OrgDb = dataset)
        ENSEMBL_ID <- paste0(esb_id$ENSEMBL, collapse = '/')
        gnnm <- DOSE::EXTID2NAME(gnid_spt, keytype = 'ENTREZID', OrgDb = dataset)
        geneName <- paste0(gnnm, collapse = '/')
        rctgID <- paste0(org,':',paste0(gnid_spt, collapse = paste0('/',org,':'), recycle0 = TRUE))
        rct_gnnm <- as.data.frame(cbind(ID=gid, geneID=ENSEMBL_ID, geneName=geneName, keggID=gnid))
        rctnm_lst[[k]] <- rct_gnnm
      }
      rct_gnnm <- bind_rows(rctnm_lst)
      rct_rlsts <- left_join(rctrlst[,-8], rct_gnnm, by=c('ID'='ID'))
      rct_rlsts <- rct_rlsts[,c(1:7,9:11,8)]
      write.csv(as.data.frame(reactome@result), paste0(ReactomePAwd, '/', nm, '_Reactome_enrich.csv'), row.names = FALSE)
      write.csv(subset(reactome@result, reactome@result$p.adjust<0.05), paste0(ReactomePAwd, '/', nm, '_Reactome_enrich_significant.csv'), row.names = FALSE)
      
      if(nrow(reactome@result)>1){
        reactome@result$Description <- ifelse(nchar(reactome@result$Description)>60, paste0(substr(reactome@result$Description,1,60),'...'),reactome@result$Description)
        reactomePlt1 <- enrichBarplot(reactome, title = paste0('The Most Enriched Reactome Terms\n',nm))
        ggsave(paste0(nm, '_Reactome_bar.png'), path = ReactomePAwd, plot = reactomePlt1, limitsize = FALSE, width = 10, height = 8)
        ggsave(paste0(nm, '_Reactome_bar.pdf'), path = ReactomePAwd, plot = reactomePlt1, limitsize = FALSE, width = 10, height = 8)
        
        reactomePlt2 <- enrichDotplot(reactome, title=paste0('The Most Enriched Reactome Terms\n',nm))
        ggsave(paste0(nm, '_Reactome_dot.png'), path = ReactomePAwd, plot = reactomePlt2, limitsize = FALSE, width = 10, height = 8)
        ggsave(paste0(nm, '_Reactome_dot.pdf'), path = ReactomePAwd, plot = reactomePlt2, limitsize = FALSE, width = 10, height = 8)
      }
    }
    
    # PPI
    if(!dir.exists(paste0(rslt, '/', Sample_Name, '_PPI/', nm))){dir.create(paste0(rslt, '/', Sample_Name, '_PPI/', nm), recursive = TRUE)}
    PPIwd <- paste0(rslt, '/', Sample_Name, '_PPI/', nm)
    
    diff_genes <- diffexp_dat$Gene_ID
    if(dataset=="org.Hs.eg.db"){
      species2 <- 9606  ##人9606,小鼠10090  
    }else if(dataset=="org.Ms.eg.db"){
      species2 <- 10090
    }else if(dataset=="Rabbit"){
      species2 <- 9986
    }else if(dataset=="Bomo"){
      species2 <- 7091
    }
    SC_strdb_network(diff_genes, output_dir = PPIwd, species = species2, nm, db)
  }
}


###TF
if(!dir.exists(paste0(rslt, '/', Sample_Name, '_TF'))){dir.create(paste0(rslt, '/', Sample_Name, '_TF'), recursive = TRUE)}
tf_list <- read.delim2(paste0(db,'/Homo_sapiens_TF.txt'))
Su_diffexp <- paste0(rslt,'/', Sample_Name, '_seurat/DIFF/')
Su_diff <- dir(Su_diffexp, pattern = 'Cluster_\\d+_diff.csv')
Su_diffnm <- stringr::str_replace(Su_diff, '_diff.csv', '')
for (t in 1:length(Su_diffnm)) {
  Su_diffdat <- read.csv(paste0(Su_diffexp,'/',Su_diff[t]))
  Su_diffdat <- left_join(Su_diffdat, tf_list[c(3,4)], by=c('Gene'='Ensembl'))
  Su_diffdat <- subset(Su_diffdat, !is.na(Su_diffdat$Family))
  Su_diffdat <- Su_diffdat[, c(1:2, 9, 4:8)]
  colnames(Su_diffdat)[c(2,3)] <- c('Gene_ID','TF_Family')
  if(!dir.exists(paste0(rslt, '/', Sample_Name, '_TF/', '/',Su_diffnm[t]))){dir.create(paste0(rslt, '/', Sample_Name, '_TF/', '/',Su_diffnm[t]))}
  write.csv(Su_diffdat, paste0(rslt, '/', Sample_Name, '_TF/', '/',Su_diffnm[t],'/',Su_diffnm[t],'_tf.csv'), row.names = FALSE)
} 

