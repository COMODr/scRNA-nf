# Nextflow + Conda  
  
###    Nextflow<https://www.nextflow.io/docs/latest/index.html> is a bioinformatics workflow manager that enables the development of portable and reproducible workflows. It ###supports deploying workflows on a variety of execution platforms including local, HPC schedulers, AWS Batch, Google Genomics Pipelines, and Kubernetes. Additionally, it ###provides support for manage your workflow dependencies through built-in support for Conda, Docker, Singularity, and Modules.


1、安装nextflow & conda
a) 先装conda，方法自行百度
b) 安装nextflow，conda install -c bioconda nextflow

2、RUN
a) scRNA.nf --reads <'./10X_scrna/raw_data/*_{S1_L001_R1_001,S1_L001_R2_001}.fastq.gz'> --outdir <'./results'> --cpus 30，也可以修改nextflow.config配置文件后，直接scRNA.nf
b) from Github: nextflow run http://github.com/XXXX/scRNA-nf（未测试过）
c) -resume参数，跳过已执行process，继续后续process
d) 根据-env.yaml自动生成虚拟环境在work工作目录中
d) 执行会生成中间文件（work目录），产生大量冗余，消耗存储空间
e) 流程完成会自动发生邮件通知结果


3、nextflow.config

4、Process instructions

5、Bug unknown
a) 6.Enrichment/keggpathway/pathview()可能会download超时导致任务终断
b) 多样本整合，可能存在I/O冲突

6、Plan
a) 参数设置，需要最优算法自动设置
b) Seurat分析报告自动生成还未添加

```
#!/bin/bash 
set -e 

for x in *.nf; do 
  (
        printf "\n\n== Testing > $(basename $x) ==\n\n"  
        nextflow run $x 
  )
done
```
