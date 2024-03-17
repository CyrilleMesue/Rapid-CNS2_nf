process check_bam_has_meth_tags {
    input:
        path(inputBam)
        val(threads)
    
    output:
        stdout emit: meth_check

    script:
        """
        samtools \
        view \
        -@${threads} \
        ${inputBam} \
        | grep -m 1 MM:Z
        """
}

process modkit_adjust_mods {
    input:
        path(inputBam)
        val(id)
        val(modkit_threads)

    publishDir("${outdir}/bam")

    output:
        path "*_modkit_merge.bam", emit: modkit_merged_bam
        path "*_modkit_merge.bam.bai", emit: modkit_merged_bai

    script:
        """
        modkit \
        adjust-mods \
        --convert h m \
        ${inputBam} \
        ${id}_modkit_merge.bam \
        --threads ${modkit_threads}

        samtools index ${id}_modkit_merge.bam
        """
}

process methylationCalls {
    memory {(1.GB * params.modkit_threads * task.attempt) + 3.GB}
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    publishDir("${out_dir}/mods/")
    input:
        path(modkit_merged_bam)
        path(ref)
        val options

    output:
        path "${id}*.bedmethyl", emit: bedmethyl_file

    script:
    """
    modkit pileup \\
        ${modkit_merged_bam} \\
        ${id}.mods.bedmethyl \\
        --ref ${ref} \\
        --threads ${modkit_threads}
    """
}

process check_mgmt_coverage {
    input:
        path(inputBam)
	    path(mgmtBed)
	    val(minimum_mgmt_cov)
        val(cov_threads)

    publishDir("${out_dir}/mgmt")

    output:
	val true
	path "*_cov.txt", emit: mgmt_avg_cov_file
    path "mgmt_cov.mosdepth.summary.txt"	
    stdout emit: mgmt_avg_cov

    script:
        """
        /mosdepth \
        -t ${cov_threads} \
        -n \
        --by ${mgmtBed} \
        mgmt_cov ${inputBam}
        
        cov="\$(grep "^chr10_region" mgmt_cov.mosdepth.summary.txt | awk '{ print \$4 }')"
        
        echo \${cov}
        if awk 'BEGIN{exit ARGV[1]>ARGV[2]}' "\$cov" ${minimum_mgmt_cov}
        then
            echo \${cov} > mgmt_below_thresh_cov.txt
        else
            echo \${cov} > mgmt_avg_cov.txt
        fi
	"""
}

process mgmtPromoter_methyartist {
    maxRetries 5
    errorStrategy { (task.attempt <= maxRetries) ? 'retry' : 'ignore' }

    input:
        path(modkit_merged_bam)
        path(modkit_merged_bai)
        path(ref)
	    val ready // mgmt_coverage has run

    publishDir("${out_dir}/mgmt")
    
    output:
        val true
        path "*.svg", emit: mgmt_plot optional true
    
    script:
        cov_file = file("${out_dir}/mgmt/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )    
            """
            methylartist \
            locus \
            -i chr10:129466536-129467536 \
            -b ${modkit_merged_bam} \
            --ref ${ref} \
            --motif CG \
            --mods m \
            --highlightpalette viridis \
            --samplepalette magma > ${id}_mgmt.svg

            """
        else

            """
            exit 1
            """
}

process mgmtPred {
    input:
        path(mgmtScript)
        file(mgmtBed)
        path(mgmtProbes)
        path(mgmtModel)
        val(id)

    output:
        val true

    script:
        cov_file = file("${out_dir}/mgmt/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )
            """
            Rscript ${mgmtScript} \
            --input ${mgmtBed} \
            --probes ${mgmtProbes} \
            --model ${mgmtModel} \
            --out_dir ${out_dir}/mgmt \
            --sample ${sample} \
            """
        else
            """
            """
}