process copyNumberVariants {
    input:
        path(inputBam)
        val(sample)

    output:
        val ${sample}_cnvpytor_100k.pdf, emit:cnv_plot
    
    publishDir("${out_dir}/cnv")

    script:
        """
        cnvpytor -root ${id}_CNV.pytor -rd ${inputBam} -j ${cnv_threads}
        cnvpytor -root ${id}_CNV.pytor -his 1000 10000 100000 -j ${cnv_threads} 
        cnvpytor -root ${id}_CNV.pytor -partition 1000 10000 100000 -j ${cnv_threads} # SLOW
        cnvpytor -root ${id}_CNV.pytor -call 1000 -j ${cnv_threads} > ${id}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 10000 -j ${cnv_threads} > ${id}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 100000 -j ${cnv_threads} > ${id}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${id}_CNV.pytor -plot manhattan 100000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o ${sample}_cnvpytor_100k.pdf
        """
}