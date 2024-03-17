<div align="left">
<img alt="Screenshot 2024-03-17 at 10 41 16 PM" src="https://github.com/areebapatel/Rapid-CNS2_nf/assets/46373444/6ba243da-0dca-4f4e-9cea-a4df7b989ff6" width="200" height="200" style="float: left; margin-right: 10px;">
<h1 style="display: inline-block;">Rapid-CNS2 workflow</h1>
</div>

## Overview

The Rapid-CNS<sup>2</sup> nextflow pipeline is a bioinformatics workflow designed for comprehensive analysis of genomic and epigenomic data generated using adaptive sampling based sequencing of central nervous system (CNS) tumours. It performs tasks such as basecalling, variant calling, methylation analysis, structural variant calling, copy number variation calling, and provides a comprehensive molecular diagnostic-ready report.

This pipeline is implemented using Nextflow, allowing for easy execution and scalability on various compute environments, including local machines, clusters, and cloud platforms.

## Features

- Modular architecture for easy customization and extension.
- Supports both basecalling from raw ONT POD5s and analysis of pre-aligned BAM files.
- Accelerated variant calling with Clara Parabricks supported Deepvariant and Sniffles2
- Annotation and filtering of clinically relevant variants
- Includes methylation analysis with Rapid-CNS<sup>2</sup> classifier and MGMT promoter methylation status determination.
- Automated report generation for summarizing analysis results.
- Prepare input files for the MNP-Flex classifier (optional)

## Requirements

- Nextflow (version 3.0.0 or later)
- Conda, Docker or Singularity (optional, for containerized execution of tools)
- Required input data:
  - Raw ONT POD5 data (for basecalling) or pre-aligned BAM files
  - Reference genome file (hg38 required)

## Usage

1. Clone this repository:

    ```bash
    git clone https://github.com/areebapatel/Rapid-CNS2_nf.git
    ```

2. Edit the `nextflow.config` file to configure pipeline parameters according to your requirements.

3. Run the pipeline using Nextflow:

    ```bash
    nextflow run main.nf --input <input_directory> --id <sample_identifier> [--options]
    ```

    Replace `<input_directory>` with the path to the directory containing ONT POD5 data or pre-aligned BAM files, and `<sample_identifier>` with a unique identifier for the sample.

    Additional options can be specified to customize pipeline behavior. Use the `--help` option to view available options and their descriptions.

4. Monitor pipeline progress and access results in the specified output directory.

## Sequencing
This pipeline analyses CNS tumour data generated through Nanopore adaptive sampling using [ReadFish](https://github.com/LooseLab/readfish) or adaptive sampling on MinKNOW. It is compatible with data generated on MinION, GridION and PromethION


## Contributions
Contributions are welcome! If you encounter any issues, have suggestions for improvements, or would like to contribute new features, please open an issue or pull request on this repository.

## License

This project is licensed under the [MIT License](LICENSE).
