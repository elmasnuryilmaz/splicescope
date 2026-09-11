#!/usr/bin/env nextflow

/*
 * splicescope — Nextflow (DSL2) pipeline
 * Runs the same analysis as the CLI, but as portable, containerisable processes
 * that scale from a laptop to an HPC/cloud scheduler.
 *
 *   nextflow run nextflow/main.nf -profile test          // self-contained demo
 *
 *   nextflow run nextflow/main.nf \
 *       --sj_dir  /path/to/SJ_out_tabs \
 *       --gtf     gencode.v47.annotation.gtf.gz \
 *       --groups  groups.tsv \
 *       --genome  GRCh38.primary_assembly.genome.fa   // optional; needs its .fai
 */

nextflow.enable.dsl = 2

params.outdir     = 'results'
params.genes      = 12
params.replicates = 6
params.sj_dir     = null   // provide real STAR SJ.out.tab dir to skip simulation
params.gtf        = null
params.groups     = null
params.genome     = null   // indexed FASTA; enables protein-consequence prediction
params.min_reads  = 10

process SIMULATE {
    tag 'simulate'
    publishDir "${params.outdir}/data", mode: 'copy'

    output:
    path 'demo_data/sj',             emit: sj
    path 'demo_data/annotation.gtf', emit: gtf
    path 'demo_data/groups.tsv',     emit: groups
    path 'demo_data/genome.fa*',     emit: genome   // the FASTA and its .fai

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
    path sj_dir
    path gtf
    path groups
    path genome_files   // empty, or the FASTA plus its .fai
    val  genome_arg     // '' when no genome was given

    output:
    path 'out', emit: results

    script:
    """
    splicescope run \\
        --sj-dir ${sj_dir} \\
        --gtf ${gtf} \\
        --groups ${groups} \\
        --min-reads ${params.min_reads} \\
        ${genome_arg} \\
        --outdir out
    """
}

workflow {
    if (params.sj_dir && params.gtf && params.groups) {
        // real-data mode: each input is staged separately, so they may live anywhere.
        // Previously a single channel was passed for all three and the script then
        // looked for <sj_dir>/annotation.gtf, which no real layout provides.
        sj_dir = Channel.fromPath(params.sj_dir, type: 'dir', checkIfExists: true)
        gtf    = Channel.fromPath(params.gtf,    checkIfExists: true)
        groups = Channel.fromPath(params.groups, checkIfExists: true)

        if (params.genome) {
            // the .fai must travel with the FASTA or GenomeFasta cannot address it
            genome_files = Channel.value(
                [file(params.genome, checkIfExists: true),
                 file("${params.genome}.fai", checkIfExists: true)]
            )
            genome_arg = Channel.value("--genome ${file(params.genome).name}")
        } else {
            genome_files = Channel.value([])
            genome_arg   = Channel.value('')
        }
    } else {
        // self-contained demo mode: the simulator writes a genome too, so the
        // consequence layer runs here as well.
        SIMULATE()
        sj_dir       = SIMULATE.out.sj
        gtf          = SIMULATE.out.gtf
        groups       = SIMULATE.out.groups
        genome_files = SIMULATE.out.genome
        genome_arg   = Channel.value('--genome genome.fa')
    }

    SPLICESCOPE_RUN(sj_dir, gtf, groups, genome_files, genome_arg)
}
