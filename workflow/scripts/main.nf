#!/usr/bin/env nextflow

params.input_vcf = ''
params.output_dir = params.output_dir ?: "${workflow.launchDir}"
params.out_root = "${params.output_dir}/PS_scoring_" + new Date().format('yyyyMMdd-HHmmss')

include { SPLICEAI; VEP; PS } from './module/processes.nf'

workflow {
    Channel.fromPath(params.input_vcf)
        | map { it -> 
                tuple(it, "${it}.*i", "${params.reference}", "${params.annotation_gtf}") 
                }
        | SPLICEAI
        | VEP
        | PS
}

workflow.onComplete {
      println ""
      println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
      println "Pipeline completed at: $workflow.complete"
      println "Execution time       : $workflow.duration"
      println "Execution status     : ${ workflow.success ? 'OK' : 'failed' }"
      println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
}
