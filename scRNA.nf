#!/usr/bin/env nextflow

log.info nfcoreHeader()

nextflow.enable.dsl = 2

if (params.reads){ raw_reads = Channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2 ).ifEmpty{ exit 1, "ERROR:Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
}

/*
--------------------------------------------------------------------------------
CLEAN FASTQ 
1.QC
--------------------------------------------------------------------------------
*/
process CLEAN_FASTQ_QC {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/cleandata", mode: 'copy'
	publishDir "${params.outdir}/1.QC/2.fastp_html", pattern: "*.html", mode: 'copy'
    input:
    tuple val(sample), path(raw_reads)
	
    output:
    tuple val(sample), path('*.fastq.gz')     , emit: cleanqs
    tuple val(sample), path('*.json')         , emit: json
    tuple val(sample), path('*.html')         , emit: html
    tuple val(sample), path('*.log')          , emit: log

    script:
    """
        fastp \\
            -i ${raw_reads[0]} \\
			-o ${sample}_clean_S1_L001_R1_001.fastq.gz \\
            -I ${raw_reads[1]} \\
            -O ${sample}_clean_S1_L001_R2_001.fastq.gz \\
			-n 15 -q 20 -u 50 \\
            --json ${sample}.fastp.json \\
            --html ${sample}.fastp.html \\
            --thread ${params.cpus} \\
            --detect_adapter_for_pe \\
            2> ${sample}.fastp.log 
    """

    }

/*
--------------------------------------------------------------------------------
2.Summary CellRanger
--------------------------------------------------------------------------------
*/
process CELLRANGER {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/2.Summary", pattern: "${sample}_cellranger", mode: 'copy', saveAs:{ filename -> "${sample}" }
    input:
    tuple val(sample), path(cleanqs)
	
    output:
	tuple val(sample), path("${sample}_cellranger") , emit: cellranger
    
    script:
    """
        $baseDir/${params.cellranger}/cellranger count --id=${sample} \
			--fastqs=. \
			--transcriptome=$baseDir/${params.refdata} \
		&& mkdir ${sample}_cellranger \
		&& mv ${sample}/outs/* ${sample}_cellranger 
    """
}

/*
--------------------------------------------------------------------------------
3.Cell_filtering
--------------------------------------------------------------------------------
*/
process CELL_FILTER {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/3.Cell_filtering", pattern: "${sample}_cellfilter", mode: 'copy', saveAs:{ filename -> "${sample}" }
    publishDir "${params.outdir}/2.Summary/${sample}", pattern: "${sample}_gene_bar.csv", mode: 'copy'
    input:
    tuple val(sample), path(cellranger)
	
    output:
    tuple val(sample), path("${sample}_gene_bar.csv"), emit: gene_bar
	tuple val(sample), path("${sample}_cellfilter") , emit: cellfilter
    
    script:
    """
		    Rscript $baseDir/bin/3-Cell_filtering.R --CRdir . --sample ${sample}
    """
}

/*
--------------------------------------------------------------------------------
4.1.Analysis_CellRanger
--------------------------------------------------------------------------------
*/
process CELLRANGER_REANALYSIS {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/4.1.Analysis_CellRanger", pattern: "${sample}_reanalysis", mode: 'copy', saveAs:{ filename -> "${sample}" }
    publishDir "${params.outdir}/4.1.Analysis_CellRanger", pattern: "${sample}.cloupe", mode: 'copy'
    publishDir "${params.outdir}/4.1.Analysis_CellRanger", pattern: "${sample}_web_summary.html", mode: 'copy'
    input:
    tuple val(sample), path(cellranger)
    tuple val(sample), path(cellfilter)
    
    output:
	tuple val(sample), path("${sample}_reanalysis") , emit: reanalysis

    script:
    """
        $baseDir/${params.cellranger}/cellranger reanalyze --id=${sample} \
			--matrix=${cellranger}/filtered_feature_bc_matrix.h5 \
        && mkdir ${sample}_reanalysis \
        && mv ${sample}/outs/analysis ./${sample}_reanalysis/ \
        && mv ${sample}/outs/cloupe.cloupe ./${sample}.cloupe \
        && mv ${sample}/outs/web_summary.html ./${sample}_web_summary.html \
        && Rscript $baseDir/bin/4_1-Analysis_CellRanger.R --sample ${sample} 
	"""
}

/*
--------------------------------------------------------------------------------
4.2.Analysis_Seurat
--------------------------------------------------------------------------------
*/
process SEURAT_ANALYSIS {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/4.2.Analysis_Seurat", pattern: "${sample}_seurat", mode: 'copy', saveAs:{ filename -> "${sample}" }
    input:
    tuple val(sample), path(cellranger)
    tuple val(sample), path(cellfilter)
    tuple val(sample), path(reanalysis)
    
    output:
    tuple val(sample), path("${sample}_seurat") , emit: seurat

    script:
    """
        Rscript $baseDir/bin/4_2-Analysis_Seurat.R --sample ${sample} 
	"""
}

/*
--------------------------------------------------------------------------------
5.*
6.Enrichment
--------------------------------------------------------------------------------
*/
process ENRICHMENT {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/6.Enrichment/GO", pattern: "${sample}_GO", mode: 'copy', saveAs:{ filename -> "${sample}" }
	publishDir "${params.outdir}/6.Enrichment/KEGG", pattern: "${sample}_KEGG", mode: 'copy', saveAs:{ filename -> "${sample}" }
	publishDir "${params.outdir}/6.Enrichment/Reactome", pattern: "${sample}_Reactome", mode: 'copy', saveAs:{ filename -> "${sample}" }
	publishDir "${params.outdir}/6.Enrichment/PPI", pattern: "${sample}_PPI", mode: 'copy', saveAs:{ filename -> "${sample}" }
	publishDir "${params.outdir}/6.Enrichment/TF", pattern: "${sample}_TF", mode: 'copy', saveAs:{ filename -> "${sample}" }
    input:
	tuple val(sample), path(cellfilter)
    tuple val(sample), path(seurat)

    
    output:
    tuple val(sample), path("${sample}_GO") , emit: GO
	tuple val(sample), path("${sample}_KEGG") , emit: KEGG
	tuple val(sample), path("${sample}_Reactome") , emit: Reactome
	tuple val(sample), path("${sample}_PPI") , emit: PPI
	tuple val(sample), path("${sample}_TF") , emit: TF
	
    script:
    """
        Rscript $baseDir/bin/6-Enrichment.R --sample ${sample} --db $baseDir/db
	"""
}

/*
--------------------------------------------------------------------------------
7.Trajectories
--------------------------------------------------------------------------------
*/
process TRAJECTORIES {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    publishDir "${params.outdir}/7.Trajectories", pattern: "${sample}_Trajectories", mode: 'copy', saveAs:{ filename -> "${sample}" }
    input:
    tuple val(sample), path(seurat)
    
    output:
    tuple val(sample), path("${sample}_Trajectories") , emit: trajectories
	
    script:
    """
        Rscript $baseDir/bin/7-Trajectories.R --sample ${sample} 
	"""
}

/*
--------------------------------------------------------------------------------
0.Report
--------------------------------------------------------------------------------
*/
process CELLRANGER_REPORT {
	conda "$baseDir/scrna-env.yaml"
	
    tag "$sample"
    
    input:
    path(cellranger_rmd)
    tuple val(sample), path(finall_out1)
	tuple val(sample), path(finall_out2)
    
    script:
    """
        Rscript -e "if(! require(rmarkdown)) install.packages('rmarkdown', repos='https://cloud.r-project.org/', dependencies = TRUE);rmarkdown::render(paste0('./scRNA-Cellranger.Rmd'),params = list(sample = '${sample}', project = 'project name'),output_file=paste0('${sample}-Cellranger.html'))"
        rm $outDir/scRNA-Cellranger.Rmd
    """
}


/*================================================================================
                      THIS IS THE START OF THE PIPELINE
--------------------------------------------------------------------------------
==================================================================================
*/
workflow ANALYSIS {
	take:
	raw_reads
	
	main:
    CLEAN_FASTQ_QC ( raw_reads )
	CELLRANGER ( CLEAN_FASTQ_QC.out.cleanqs )
	CELL_FILTER ( CELLRANGER.out.cellranger )
	CELLRANGER_REANALYSIS ( CELLRANGER.out.cellranger , CELL_FILTER.out.cellfilter )
	SEURAT_ANALYSIS ( CELLRANGER.out.cellranger , CELL_FILTER.out.cellfilter , CELLRANGER_REANALYSIS.out.reanalysis )
	ENRICHMENT ( CELL_FILTER.out.cellfilter , SEURAT_ANALYSIS.out.seurat )
	TRAJECTORIES ( SEURAT_ANALYSIS.out.seurat )
	
	emit:
	finall_out1 = ENRICHMENT.out.TF
	finall_out2 = TRAJECTORIES.out.trajectories
}

workflow REPORT {
	take:
	cellranger_rmd
	finall_out1
	finall_out2
	
	main:
    CELLRANGER_REPORT( cellranger_rmd , finall_out1 , finall_out2 )
}


if(params.outdir =~ '^.'){
    outDir = baseDir + '/' +params.outdir.replaceFirst('.', '')
}else{
    outDir = params.outdir
}
workflow {
    ANALYSIS ( raw_reads )
	
	"Rscript $baseDir/bin/1-QC_datatable.R --dir $outDir".execute().text
	"cp -rf $baseDir/bin/Loupe-Cell-Browser-3.1.0/ $outDir/4.1.Analysis_CellRanger/".execute().text
	"cp $baseDir/bin/scRNA-Cellranger.Rmd $outDir".execute().text
    "cp -rf $baseDir/images $outDir/".execute().text
	Channel.fromPath("${outDir}/*.Rmd").set { cellranger_rmd }
	REPORT ( cellranger_rmd , ANALYSIS.out.finall_out1.collect() , ANALYSIS.out.finall_out2.collect() )
}





/*================================================================================
                              PIPELINE COMPLETE
--------------------------------------------------------------------------------
==================================================================================
*/
workflow.onComplete {
	if(workflow.success){
		"find $workDir -name .report | xargs rm".execute().text
		"find $workDir -name stdout | xargs rm".execute().text
		"find $workDir -name .stdout | xargs rm".execute().text
		"find $workDir -name checks.out | xargs rm".execute().text
		"find $workDir -name .cache | xargs rm -rf".execute().text
		"find $workDir -name '.nextflow*' | xargs rm -rf".execute().text
		"find $workDir -name '.node-nextflow*' | xargs rm -rf".execute().text
		"find $workDir -name 'work' | xargs rm -rf".execute().text
	}

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: params.mailto, from: '1417188187@qq.com', subject: 'My pipeline execution', body: msg)
}













def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = "\033[0;30m";
    c_blue = "\033[0;34m";
    c_cyan = "\033[0;36m";
    c_dim = "\033[2m";
    c_green = "\033[0;32m";
    c_purple = "\033[0;35m";
    c_reset = "\033[0m";
    c_white = "\033[0;37m";
    c_yellow = "\033[0;33m";

    return """    -${c_dim}-----------------------------------------------------------------------------${c_reset}-
    
    ${c_blue}                         ${c_yellow}   .   ${c_green}                                                ${c_reset}
    ${c_blue}                         ${c_yellow}  @88> ${c_green}                                                ${c_reset}
    ${c_blue}       u.      .u    .   ${c_yellow}  %8P  ${c_green}                           u.    u.             ${c_reset}
    ${c_blue} ...ue888b   .d88B :@8c  ${c_yellow}   .   ${c_green}      uL          .u     x@88k u@88c.      .u   ${c_reset}
    ${c_blue} 888R Y888r ="8888f8888r ${c_yellow} .@88u ${c_green}  .ue888Nc..   ud8888.  ^"8888""8888"   ud8888. ${c_reset}
    ${c_blue} 888R I888>   4888>'88"  ${c_yellow}''888E`${c_green} d88E`"888E` :888'8888.   8888  888R  :888'8888.${c_reset}
    ${c_blue} 888R I888>   4888> '    ${c_yellow}  888E ${c_green} 888E  888E  d888 '88%"   8888  888R  d888 '88%"${c_reset}
    ${c_blue} 888R I888>   4888>      ${c_yellow}  888E ${c_green} 888E  888E  8888.+"      8888  888R  8888.+"   ${c_reset}
    ${c_blue}u8888cJ888   .d888L .+   ${c_yellow}  888E ${c_green} 888E  888E  8888L        8888  888R  8888L     ${c_reset}
    ${c_blue} "*888*P"    ^"8888*"    ${c_yellow}  888& ${c_green} 888& .888E  '8888c. .+  "*88*" 8888" '8888c. .+${c_reset}
    ${c_blue}   'Y"          "Y"      ${c_yellow}  R888"${c_green} *888" 888&   "88888%      ""   'Y"    "88888%  ${c_reset}
    ${c_blue}                         ${c_yellow}   ""  ${c_green}  `"   "888E    "YP'                     "YP'   ${c_reset}
    ${c_blue}                         ${c_yellow}       ${c_green} .dWi   `88E                                    ${c_reset}
    ${c_blue}                         ${c_yellow}       ${c_green} 4888~  J8%                                     ${c_reset}
    ${c_blue}                         ${c_yellow}       ${c_green}  ^"===*"`                                      ${c_reset}
    ${c_purple}         zgs0000@yeah.net         www.knorigene.com     ${c_reset}
    -${c_dim}-----------------------------------------------------------------------------${c_reset}-
    """.stripIndent()
}

















