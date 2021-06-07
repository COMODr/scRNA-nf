
options(stringsAsFactors = FALSE)

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
