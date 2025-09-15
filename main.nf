#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/immunogenomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/immunogenomics
    Website: https://nf-co.re/immunogenomics
    Slack  : https://nfcore.slack.com/channels/immunogenomics
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IMMUNOGENOMICS  } from './workflows/immunogenomics'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_immunogenomics_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_immunogenomics_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_immunogenomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis workflow depending on type of input
//
workflow NFCORE_IMMUNOGENOMICS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run main analysis workflow
    //
    IMMUNOGENOMICS (
        samplesheet
    )

    emit:
    multiqc_report = IMMUNOGENOMICS.out.multiqc_report // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD TEST DATA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNLOAD_TEST_DATA {
    
    main:
    
    // Check if download_test_data script exists
    def download_script = file("$projectDir/bin/download_test_data.py")
    if (!download_script.exists()) {
        error("Test data download script not found at: $projectDir/bin/download_test_data.py")
    }
    
    // Set output directory for test data
    def test_data_dir = params.test_data_outdir ?: "assets/test_data"
    
    // Determine data type to download
    def data_type = params.test_data_type ?: "all"
    
    // Build download command
    def download_cmd = "python3 $download_script --outdir $test_data_dir --data-type $data_type"
    if (params.force_download) {
        download_cmd += " --force"
    }
    
    log.info """
    
    ================================================
    Downloading nf-core/immunogenomics Test Data
    ================================================
    
    Dataset: Texas Cancer Research Biobank (TCRB)
    Source: nf-core/test-datasets
    Data type: $data_type
    Output directory: $test_data_dir
    
    ⚠️  This data is for testing purposes only
    ⚠️  Do not attempt to re-identify participants
    
    ================================================
    """.stripIndent()
    
    // Execute download
    def result = download_cmd.execute()
    result.waitFor()
    
    if (result.exitValue() == 0) {
        log.info """
        
        ✅ Test data download completed successfully!
        
        To test the pipeline:
          WES:     nextflow run . --input $test_data_dir/samplesheet_wes.csv
          RNA-seq: nextflow run . --input $test_data_dir/samplesheet_rnaseq.csv
        
        Example commands:
          nextflow run . --input $test_data_dir/samplesheet_wes.csv --outdir results_wes -profile docker
          nextflow run . --input $test_data_dir/samplesheet_rnaseq.csv --outdir results_rnaseq -profile docker
        
        """.stripIndent()
    } else {
        error("Test data download failed. Check network connection and try again.")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // Handle special workflows
    //
    if (params.download_test_data) {
        DOWNLOAD_TEST_DATA()
        return
    }

    
    //
    // PIPELINE INITIALISATION
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_IMMUNOGENOMICS (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // PIPELINE COMPLETION
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_IMMUNOGENOMICS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/