#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
version = 3.0.0
/**
 ********************************** Rapid-CNS2 NextFlow ******************************************
 * 1 - Base calling, alligment and data preparation
 *  a. Base calling + alignment : ONT Dorado and minimap2
 *  b. Sorting, adding read group information and creating index : samtools
 *  c. Subsetting to target region
 * 2 - Long read variant calling and annotation
 *  a. Clara Parabricks Deepvariant 
 *  b. Annotation using ANNOVAR
 *  c. Filtering
 * 3 - Structural variants and annotation 
 *  a. Sniffles2
 *  b. Annotation using AnnotSV
 * 4 - Methylation analysis
 *  a. Methylation values using modkit
 *  b. Methylation classification using Rapid-CNS2 classifier
 *  c. MGMT promoter methylation status
 *  d. MGMT promoter region plot
 * 5 - Copy number variation calling using CNVpytor
 * 6 - Report rendering
 *******************************************************************************************
*/

   

params.input = null
params.ref = null
params.tmp_dir = "tempDir"
params.out_dir = "output"
params.log_dir = "logDir"
params.minimum_mgmt_cov = 5
params.model_config="dna_r10.4.1_e8.2_400bps_hac@v4.1.0"
params.port= 8887
params.num_gpu = 3
params.num_clients = par.num_gpu * 3


params.rParams = [] // Initialize an empty list to store -r parameters

// DeepVariant mode. 
params.pbDVMode = "ont"
params.pbPATH = "pbrun"

// set up and create an output directory
out_dir = file(params.out_dir)
out_dir.mkdir()

log.info """ \
        Inputs
        ================================================================
        sample id               : ${params.id}
        input                   : ${params.input}
        input_bam               : ${params.bam}
        outdir                  : ${params.out_dir}
        num_gpu                 : ${params.num_gpu}
        reference               : ${params.ref}
        annovar                 : ${params.annovar}
        threads                 : ${params.threads}
        panel                   : ${params.panel}
        min_mgmt_coverage       : ${params.minimum_mgmt_cov}
        mnp-flex                : ${params.mnp_flex}
        model_config

        To run with LSF, add -process.executor='lsf' to your nextflow command.

        To ALSO prepare input file for the MNP-Flex classifier, add the --mnp-flex flag. (Default behaviour is to not prepare the input file)
        """ 

params.help = null
params.test = null
params.reads = null

// Show help message
if (params.help) {
   log.info """
    Usage: nextflow run nextflow/main.nf  [options]

    Mandatory options:
        --input            Path to the directory containing Pod5/FAST5 file for Dorado basecalling and minimap2 alligement or BAM file if available
        --id                Sample identifier

    Options:
        --out_dir          Directory path to store all the outputs. [default : ${params.out_dir}]
        --ref              Path to hg19 reference file. [default : ${params.ref}]
        --tmp_dir          Directory to store temporary files. If it does not exists it will be created. [default : ${params.tmp_dir}]
        --basecalling      If data should be basecalled  // remove
        --sort_threads     Number of threads to use for samtools sort [default : ${params.sort_threads}]
        
        --sniffles_threads  Numbers of threads to run sniffles2 [default : ${params.sniffles_threads}]
        --cnv_threads Numbers of threads to run cnvpytor [default : ${params.cnv_threads}]
        --mod_threads Numbers of threads to run modkit [default : ${params.mod_threads}]
        --model_config     Basecalling model to be used [default : ${params.model_config}]
        --port             Port for basecall server [default : ${params.port}] 
        --reads            samtools addreplacerg -r option. It should be specified as this example :  --reads '-r "SM:GM24385" -r "ID:GM24385"'

    Example:
      nextflow run main.nf  --input "/input/pod5" --ref "/Ref/Homo_sapiens_assembly38.fasta" --id "Sample_XYZ"
    """
    exit 0
}

// Verify that the mandatory parameters are provided
if (params.basecalling & params.ref == null) error "The reference genome file is mandatory for basecalling. Please specify it with --ref"
if (params.input == null) error "The path to the input POD5 files or BAM file is mandatory, please specify it with --input"

include { basecalling } from './nextflow/basecalling.nf'
include { mergeBAM } from './nextflow/basecalling.nf'
include { variantCalling } from './nextflow/variantCalling.nf'
include { structuralVariants } from './nextflow/structuralVariants.nf'
include { methylationCalls } from './nextflow/methylationAnalysis.nf'
include { methylationClassification } from './nextflow/methylationAnalysis.nf'
include { mgmtPromoter } from './nextflow/methylationAnalysis.nf'
include { copyNumberVariants } from './nextflow/copyNumberVariants.nf'
include { reportRendering } from './nextflow/reportRendering.nf'
include { mnpFlex } from './nextflow/mnpFlex.nf'

workflow {
    Channel.fromPath(params.input, checkIfExists: true)
    .set {input}

    Channel.fromPath(params.ref, checkIfExists: true)
    .set {ref}

    Channel.from(params.id)
    .set {id}

    Channel.from(params.out_dir)
    .set {out_dir}

    Channel.from(params.model_config)
    .set {model_config}

    Channel.from(params.num_clients)
    .set {num_clients}

    Channel.from(params.port)
    .set {port}

    Channel.from(params.max_threads)
    .set {max_threads}

    Channel.from(params.pbDVMode)
    .set {pbDVMode}

    Channel.from(params.pbPATH)
    .set {pbPATH}

    Channel.from(params.annovarPath)
    .set {annovarPath}

    Channel.from(params.annovarDB)
    .set {annovarDB}

    Channel.from(params.annotsvAnnot)
    .set {annotsvAnnot}

    Channel.from(params.modkitThreads)
    .set {modkitThreads}

    Channel.from(params.cnvThreads)
    .set {cnvThreads}

    Channel.from(params.minimum_mgmt_cov)
    .set {minimum_mgmt_cov}

    Channel.fromPath(params.annotations, checkIfExists: true)
    .set {annotations}

    // Collect variables and scripts from bin

    Channel.fromPath("${projectDir}/data/NPHD_panel_hg38.bed", checkIfExists: true)
    .set {panel}

    Channel.fromPath("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmtBed}

    Channel.fromPath("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true)
    .set {mgmtProbes}

    Channel.fromPath("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {mgmtModel}

    Channel.fromPath("${projectDir}/scr/mgmt_pred.R", checkIfExists: true)
    .set {mgmtScript}

    Channel.fromPath("${projectDir}/scr/methylation_classification.R", checkIfExists: true)
    .set {methylationClassification}

    Channel.fromPath("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topProbes}
    
    Channel.fromPath("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingData}

    Channel.fromPath("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayFile}

    Channel.fromPath("${projectDir}/scr/filter_report.R", checkIfExists: true)
    .set {filter_report}

    Channel.fromPath("${projectDir}/scr/make_report.R", checkIfExists: true)
    .set {makereport}

    Channel.fromPath("${projectDir}/scr/Rapid_CNS2_report_UKHD_v3.0.0.Rmd", checkIfExists: true)
    .set {report_UKHD}

// Conditional step: check if basecalling is requested
    if (params.basecalling) {
        basecalling(input=input,
                inputRef=ref
            )

        mergeBAM(inputFolder=basecalling.out.bam_folder
            )
        set val(inputBam) from mergeBAM.out.inputBam
    } else {
        // If basecalling is not requested, directly use the input BAM file for variant calling
        set val(inputBam) from input
    }

    subsetBAM( inputBam=inputBam
            panel=panel
    )

    deepVariant( inputBam=subsetBAM.out.subsetBam,
            ref=ref,
            pbDVMode=pbDVMode,
            pbPATH=pbPATH,
            tmpDir=tmp_dir,
            annovarPath=annovarPath,
            annovarDB=annovarDB,
        )

    structuralVariants( inputBam=mergeBAM.out.bam_file,
            inputBai=mergeBAM.out.bai_file,
            inputRef=ref,
            annotsvAnnot=annotsvAnnot
        )
    
    methylationCalls( inputBam=mergeBAM.out.bam_file,
            inputBai=mergeBAM.out.bai_file,
            inputRef=ref,
            modkitThreads=modkitThreads
        )

    methylationClassification( bedmethyl_file=methylationCalls.out.bedmethyl_file,
            probes=topProbes,
            arrayFile=arrayFile,
            trainingData=trainingData
        )
    
    mgmtPromoter( bedmethyl_file=methylationCalls.out.bedmethyl_file,
            mgmtBed=mgmtBed,
            mgmtProbes=mgmtProbes,
            mgmtModel=mgmtModel
        )
    
    copyNumberVariants( inputBam=mergeBAM.out.bam_file,
            inputBai=mergeBAM.out.bai_file,
            inputRef=ref,
            cnvThreads=cnvThreads
        )
        // To do ***************************
    reportRendering( inputBam=basecalling.out.bam_file,
            inputBai=basecalling.out.bai_file,
            inputRef=ref
        )
    if ( params.mnpFlex) {
        mnpFlex( bedMethyl=methylationCalls.out.bedmethyl_file,
                mnpFlexBed=params.mnpFlexBed
            )
    }

}



workflow.onComplete {
    if(workflow.success) {
    println ( "The Rapid-CNS2 workflow is now complete!\n Your outpus are located in : " + params.out_dir )
    }
    else {
    println ( "Oops .. something went wrong, please look into the log file, and error messages into " + workDir )
    }
}