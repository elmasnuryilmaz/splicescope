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
params.install_splicescope = false   // set by the docker profile; see installCommand()
params.splicescope_wheel   = null    // docker profile only: install this wheel, not the release

// With the docker profile each task first installs the splicescope release that matches
// this pipeline's version. The install is part of the task script, so it runs inside the
// throwaway container, never in the user's own environment (a beforeScript would run on
// the host). Without the profile this is empty and splicescope must already be installed.
//
// --splicescope_wheel replaces the release with a wheel built from a checkout, staged into
// each task as an input so that any container engine can see it. CI uses it to run the
// code under review in the container: the release is the last published code, and the
// pull request that raises the version names one PyPI does not have yet.
def installCommand(wheel) {
    if (!params.install_splicescope) {
        return ''
    }
    def target = wheel ? "./${wheel}" : "'splicescope==${workflow.manifest.version}'"
    "pip install --quiet --no-cache-dir --disable-pip-version-check --root-user-action=ignore ${target}"
}

process SIMULATE {
    tag 'simulate'
    publishDir "${params.outdir}/data", mode: 'copy'

    input:
    path wheel          // empty, or the wheel to install; see installCommand()

    output:
    path 'demo_data/sj',             emit: sj
    path 'demo_data/annotation.gtf', emit: gtf
    path 'demo_data/groups.tsv',     emit: groups
    path 'demo_data/genome.fa*',     emit: genome   // the FASTA and its .fai

    script:
    """
    ${installCommand(wheel)}
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
    path wheel          // empty, or the wheel to install; see installCommand()

    output:
    path 'out', emit: results

    script:
    """
    ${installCommand(wheel)}
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
    // a wheel given without the docker profile would be silently ignored, and the run
    // would use whatever splicescope the environment has instead of the code asked for
    if (params.splicescope_wheel && !params.install_splicescope) {
        error "--splicescope_wheel is only used with -profile docker, which installs splicescope in each task"
    }
    wheel = params.splicescope_wheel
        ? Channel.value(file(params.splicescope_wheel, checkIfExists: true))
        : Channel.value([])

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
        SIMULATE(wheel)
        sj_dir       = SIMULATE.out.sj
        gtf          = SIMULATE.out.gtf
        groups       = SIMULATE.out.groups
        genome_files = SIMULATE.out.genome
        genome_arg   = Channel.value('--genome genome.fa')
    }

    SPLICESCOPE_RUN(sj_dir, gtf, groups, genome_files, genome_arg, wheel)
}
