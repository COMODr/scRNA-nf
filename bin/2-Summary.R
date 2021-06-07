
###单个样本
if(!dir.exists(paste0(rslt,'/2.Summary')))(dir.create(paste0(rslt,'/2.Summary')))
for (i in 1:length(Sample_Name)){
  if(!dir.exists(paste0(rslt, '/2.Summary/', Sample_Name[i])))
    dir.create(paste0(rslt, '/2.Summary/', Sample_Name[i]))
  cellranger_h <- paste0(CR_path, "/",Sample_Name[i],"/outs")
  if(dir.exists(paste0(rslt, '/2.Summary/', Sample_Name[i])))
    system(paste0('rm -rf ', rslt, '/2.Summary/', Sample_Name[i]))
  mv_p <- paste0('cp -rf ', cellranger_h, '/', ' ', rslt, '/2.Summary/')
  rnm_p <- paste0('mv ', rslt, '/2.Summary/outs ', rslt, '/2.Summary/', Sample_Name[i])
  system(paste0(mv_p, ' && ', rnm_p))
}

CELLRANGER.out.cellranger
mv ${sample}/outs ${sample}



