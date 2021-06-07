# RawData数据质量评估
if(! require(stringr)) install.packages("stringr", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(jsonlite)) install.packages("jsonlite", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(xml2)) install.packages("xml2", repos="https://cloud.r-project.org/", dependencies = TRUE)
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = ".", action = "store", help = "This is the html files of clean data process results!")
  )
opt = parse_args(OptionParser(option_list = option_list, usage = "Rscript 1-QC.R -d <html and json files by clean data process results>"))

json_dir <- paste0(opt$dir, '/cleandata/')

json_fl <- dir(json_dir,pattern = 'json$')
Sample_Name <- str_replace(json_fl, '.fastp.json', '')

raw_lst <- list()
for(json in json_fl){
  jsnm <- str_replace(json,'.fastp.json','')
  jsondat <- fromJSON(paste0(json_dir,'/',json))
  raw_data <- jsondat$read2_before_filtering$total_reads
  raw_lst[json] <- raw_data
}

base_lst <- list()
for(json in json_fl){
  jsnm <- str_replace(json,'.json','')
  jsondat <- fromJSON(paste0(json_dir,'/',json))
  raw_data <- jsondat$read2_before_filtering$total_bases
  base_lst[json] <- raw_data
}

q20_lst <- list()
for(json in json_fl){
  jsnm <- str_replace(json,'.json','')
  jsondat <- fromJSON(paste0(json_dir,'/',json))
  raw_data <- jsondat$read2_before_filtering$q20_bases
  q20_lst[json] <- raw_data
}

q30_lst <- list()
for(json in json_fl){
  jsnm <- str_replace(json,'.json','')
  jsondat <- fromJSON(paste0(json_dir,'/',json))
  raw_data <- jsondat$read2_before_filtering$q30_bases
  q30_lst[json] <- raw_data
}

gc_lst <- list()
for(json in json_fl){
  jsnm <- stringr::str_replace(json, pattern = '.json', '.html')
  page = read_html(paste0(json_dir,'/',jsnm))
  text <- rvest::html_text(page)
  k <- gregexpr('GC\\(', text)[[1]][2]
  gc_lst[json] <- as.numeric(substr(text, k+3, k+7))
}

datatable <- data.frame(Sample_Name=character(), 
                        Raw_Reads=numeric(), Raw_Bases=numeric(),
                        Q20= numeric(),Q30= numeric(),
                        GC_Content= numeric(),stringsAsFactors=FALSE)   ###创建空的数据框
for(i in 1:length(json_fl)){
  Sample_Name <- Sample_Name[i]
  Raw_Reads <- sum(unlist(raw_lst[i]))
  Raw_Bases <- round(sum(unlist(base_lst[i]))*2/1000000000, 2)
  Q20 <- round(sum(unlist(q20_lst[i]))/sum(unlist(base_lst[i]))*100,2)
  Q30 <- round(sum(unlist(q30_lst[i]))/sum(unlist(base_lst[i]))*100,2)
  GC_Content <- round(sum(unlist(gc_lst[i]))/length(gc_lst[i]),2)
  datatable[i,] <- as.data.frame(cbind(Sample_Name,Raw_Reads,Raw_Bases,Q20,Q30,GC_Content))
}
colnames(datatable) <- c("Sample_Name", "Raw_Reads",	"Raw_Bases(G)",	"Q20(%)",	"Q30(%)",	"GC_Content(%)")

if(!dir.exists(paste0(opt$dir, '/1.QC/1.dataTable/'))) dir.create(paste0(opt$dir, '/1.QC/1.dataTable/'), recursive = TRUE)
write.table(datatable, paste0(opt$dir, '/1.QC/1.dataTable/datatable.txt'), quote=F,sep = '\t',row.names=F)














