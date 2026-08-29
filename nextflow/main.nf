#!/usr/bin/env nextflow

/*
 * splicescope — Nextflow (DSL2) pipeline
 * Runs the same analysis as the CLI, but as portable, containerisable processes
 * that scale from a laptop to an HPC/cloud scheduler.
 *
 *   nextflow run nextflow/main.nf -profile test
 */

nextflow.enable.dsl = 2

params.outdir     = 'results'
params.genes      = 12
params.replicates = 6
params.sj_dir     = null   // provide real STAR SJ.out.tab dir to skip simulation
params.gtf        = null
params.groups     = null

process SIMULATE {
    tag 'simulate'
    publishDir "${params.outdir}/data", mode: 'copy'

    output:
    path 'demo_data', emit: data

    script:
    """
    splicescope simulate --outdir demo_data \\
        --genes ${params.genes} --replicates ${params.replicates}
    """
}

process SPLICESCOPE_RUN {
    tag 'run'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path data

    output:
    path 'out', emit: results

    script:
    """
    splicescope run \\
        --sj-dir ${data}/sj \\
        --gtf ${data}/annotation.gtf \\
        --groups ${data}/groups.tsv \\
        --outdir out
    """
}

workflow {
    if (params.sj_dir && params.gtf && params.groups) {
        // real-data mode: stage the provided inputs into one directory
        ch = Channel.fromPath(params.sj_dir, type: 'dir')
        SPLICESCOPE_RUN(ch)
    } else {
        // self-contained demo mode
        SIMULATE()
        SPLICESCOPE_RUN(SIMULATE.out.data)
    }
}
