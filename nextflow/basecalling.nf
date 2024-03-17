process basecalling {
    label 'GPU'
    publishDir "${params.out_dir}/basecalling", mode: 'copy', pattern: "*"

    input:
        path input
        path inputRef
	

    output:
        path "${params.out_dir}/basecalling/pass/", emit: bam_folder  //check out folder

    script:
	def r_args = params.reads ?: ''
	
        """
        dorado_basecall_server -c ${params.model_config} --port ${params.port} --device "cuda:all"  --log_path ${params.log_dir} --use_tcp &

        ont_basecaller_supervisor --input_path ${input} --save_path ${out_dir}/basecalling --config ${model_config} --bam_out  --align_ref ${ref} --compress_fastq --index --recursive --port ${port} --use_tcp --num_clients ${num_clients} --progress_stats_frequency 5
        """

}

process mergeBAM{
    publishDir "${params.out_dir}/bam", mode: 'copy', pattern: "*"

    input:
        path inputFolder
    
    output:
        path "${out_dir}/bam/${id}.sorted.index.bam", emit: inputBam

    script:
        """
        samtools merge -@ ${max_threads} -f ${out_dir}/bam/${id}.merged.sorted.bam ${out_dir}/basecalling/pass/*.bam 
        """
}

process index_bam {
    input:
        path(inputBam)
        val(index_threads)

    output:
        path "*.bai", emit: inputBai

    script:
        """
        samtools \
        index \
        -@${index_threads} \
        ${inputBam}
        """
}

process subsetBAM{
   //publishDir "${params.out_dir}/bam", mode: 'copy', pattern: "*"
    def r_args = params.reads ?: ''

    input:
        path inputBam
        path inputBai
        path panel
    
    output:
        path "${out_dir}/bam/${id}.merged.sorted.index.panel.bam", emit: subsetBam

    script:
        """
        samtools addreplacerg ${r_args} \
                -@10 -o ${out_dir}/bam/${id}.merged.sorted.index.bam \
                ${out_dir}/bam/${id}.merged.sorted.bam

        samtools index -@${max_threads} ${out_dir}/bam/${id}.merged.sorted.index.bam 
        
        bedtools intersect -a ${inputBam} -b ${panel} > ${out_dir}/bam/${id}.sorted.index.panel.bam

        samtools index -@${max_threads} ${out_dir}/bam/${id}.merged.sorted.index.panel.bam 
        """
}

// to re-write
process generatePOD5Summary {
    input:
    val inputPOD5s
    val runtime_attributes
    val pod5Docker

    output:
    file("summary.tsv") into summary

    script:
    """
    pod5 view --threads ${runtime_attributes.nThreads - 1} --include "read_id, channel" --output summary.tsv ${sep=" " inputPOD5s}
    """

    runtime:
    docker pod5Docker
    cpus runtime_attributes.nThreads
    memory "${runtime_attributes.gbRAM} GB"
    time 10.hours
}

//to re-write
process fast5ToPod5 {
    input:
    file inputFAST5
    val runtime_attributes
    val pod5Docker

    output:
    file("${inputFAST5.baseName}.pod5") into outputPOD5

    script:
    """
    pod5 convert fast5 --output ${inputFAST5.baseName}.pod5 ${inputFAST5}
    """
}

process splitPOD5ByChannel {
    input:
    val inputPOD5s
    val summaryFile
    val pod5Docker

    output:
    file "split_by_channel/*" into split_by_channel

    script:
    """
    set -e
    mkdir split_by_channel && \
    time pod5 subset --threads ${max_threads} --summary ${summaryFile} --columns channel --output split_by_channel ${sep=" " inputPOD5s}
    """
}