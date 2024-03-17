process methylationClassification {
    input:
        path(methylationClassification)
        path(bedmethyl_file)
        val(id)
        val(out_dir)
        path(topProbes)
        path(trainingData)
        path(arrayFile)
        val(meth_threads)
        val(ready) // filter to just 5mC

    publishDir("${out_dir}/methylation_classification")

    output:
        val true

    script:
        """
        Rscript ${methylationClassification} \
        --sample ${ID} \
        --out_dir ${out_dir}/methylation_classification \
        --in_file ${bedmethyl_file} \
        --probes ${topProbes} \
        --training_data ${trainingData} \
        --array_file ${arrayFile} \
        --threads ${meth_threads}
        """
}
