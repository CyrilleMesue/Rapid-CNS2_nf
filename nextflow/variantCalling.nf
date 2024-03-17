process deepVariant {
    label 'GPU'
    publishDir "${out_dir}/snv/", mode: 'copy', pattern: "*"

    //stageInMode "copy"

    input:
    path inputBam
    path ref
    val pbDVMode
    val pbPATH
    val tmpDir

    output:
    path "${out_dir}/snv/${id}.dv.vcf", emit: dv_vcf

    script:
    """
    tar xzf ${ref.Name} && \
    mkdir -p ${tmpDir} && \
    time ${pbPATH} deepvariant \
    --tmp-dir ${tmpDir} \
    --in-bam ${inputBam} \
    --ref ${ref} \
    --out-variants ${id}.dv.vcf \
    --mode ${pbDVMode} \
    --run-partition --norealign-reads
    """

}

process recodeVCF {
    publishDir "${out_dir}/snv", mode: 'copy', pattern: "*"

    input:
        path dv_vcf

    output:
        path "${out_dir}/snv/${id}.dv.PASS.vcf.gz", emit: pass_vcf

    script:
    """
	vcftools --gzvcf ${OUT_DIR}/${SAMPLE}.vcf.gz --remove-filtered-all --recode --stdout | gzip -c > ${out_dir}/snv/${id}.dv.PASS.vcf.gz
    """

}


process convert2annovar{
    publishDir "${out_dir}/snv", mode: 'copy', pattern: "*"
    input:
        path(input)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        ${annovarPath}/convert2annovar.pl \
        -format vcf4 ${input} \
        -withfreq \
        -includeinfo \
        > ${id}_deepvariant_panel.avinput
        """
}

process table_annovar {
    publishDir "${out_dir}/snv", mode: 'copy', pattern: "*"
    input:
    input:
        path(annovar_input)
    
    output:
        path "*_multianno.csv", emit: dv_anno
      
    script:
        """
        /annovar/table_annovar.pl ${annovar_input} \
        /annovar/humandb/ \
        -buildver hg38 \
        -out ${id}_dv_panel \
        -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 \
        -operation gx,r,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo
        """
}

process filter_report {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    publishDir "${out_dir}/snv", mode: 'copy', pattern: "*"
    input:

    input:
        path(filter_report)
        path(dv_anno)
        val(id)

    output:
        val ${id}_dv_report.csv, emit: dv_report
    
    script:
        """
        Rscript ${filter_report} \
        --input ${dv_anno} \
        --output ${id}_dv_report.csv \
        --sample ${id}
        """        
}