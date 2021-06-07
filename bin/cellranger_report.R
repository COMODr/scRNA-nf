
if(! require(optparse)) install.packages("optparse", repos="https://cloud.r-project.org/", dependencies = TRUE)
option_list <- list(
  make_option(c("-o", "--outDir"), type = "character", default = ".", action = "store", help = "The directory where the Seurat analysis results store"
  ),
  make_option(c("-S", "--sample"), type = "character", default = FALSE, action = "store", help = "Sample name"
  ),
  make_option(c("-T", "--title"), type = "character", default = "项目名称标题", action = "store", help = "Sample name"
  )
)
opt = parse_args(OptionParser(option_list = option_list, usage = "cellranger_report.R --outDir <CellRanger results path> --sample <Sample name> --title <project title>"))


setwd(opt$outDir)
if(! require(rmarkdown)) install.packages("rmarkdown", repos="https://cloud.r-project.org/", dependencies = TRUE)
render(paste0(opt$outDir,'/scRNA-Cellranger.Rmd'), params = list(sample = opt$sample, project = opt$title), output_file = paste0(opt$outDir, '/', opt$sample, '-Cellranger.html'))
