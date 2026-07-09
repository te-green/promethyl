#!/usr/bin/env nextflow
/*
 * promethyl cohort pipeline
 *
 * Phase 1 (MODKIT, parallel per sample): modkit pileup + island aggregation.
 * Phase 2 (COHORT, single process):      cross-sample outlier analysis.
 */

nextflow.enable.dsl = 2

params.samples = null                     // YAML: samples, plus optional annotation/ref/include_bed/
                                           // region/min_coverage/mod_code/min_delta/z_threshold
params.outdir  = "results"
params.threads = 4

if (!params.samples) { error "Provide --samples cohort.yaml" }

// modkit pileup requires the BAM to already be indexed. Nextflow only stages files
// declared as process inputs into the task work dir, so the .bai/.csi sitting next
// to the BAM on the shared filesystem must be located and passed through explicitly.
def findBamIndex(bamPath) {
    for (candidate in [bamPath + '.bai', bamPath + '.csi',
                        bamPath.replaceAll(/\.bam$/, '.bai'), bamPath.replaceAll(/\.bam$/, '.csi')]) {
        def f = file(candidate)
        if (f.exists()) { return f }
    }
    error "No BAM index (.bai/.csi) found for ${bamPath} -- run `samtools index ${bamPath}`"
}

process MODKIT {
    tag "$sample_id"
    publishDir "${params.outdir}/modkit", mode: 'copy', pattern: "*_modkit.bed"
    cpus params.threads

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path annotation
    path ref
    path include_bed
    val region
    val min_coverage
    val mod_code

    output:
    tuple val(sample_id), path("${sample_id}.islands.tsv"), emit: islands
    path "*_modkit.bed"

    script:
    def ref_opt         = ref.name != 'NO_REF'         ? "--ref ${ref}"                 : ""
    def include_bed_opt = include_bed.name != 'NO_BED' ? "--include-bed ${include_bed}" : ""
    def region_opt      = region                        ? "--region ${region}"           : ""
    """
    run_sample.py \
        --id ${sample_id} \
        --bam ${bam} \
        --modkit-dir . \
        --annotation ${annotation} \
        --output ${sample_id}.islands.tsv \
        --threads ${task.cpus} \
        --min-coverage ${min_coverage} \
        --mod-code ${mod_code} \
        ${ref_opt} ${include_bed_opt} ${region_opt}
    """
}

process COHORT {
    publishDir params.outdir, mode: 'copy'

    input:
    val sample_args
    path annotation
    val min_delta
    val z_threshold

    output:
    path "cohort_methylation.tsv"

    script:
    """
    run_cohort.py \
        --samples ${sample_args.join(' ')} \
        --annotation ${annotation} \
        --output cohort_methylation.tsv \
        --min-delta ${min_delta} \
        --z-threshold ${z_threshold}
    """
}

workflow {
    cfg = new org.yaml.snakeyaml.Yaml().load(file(params.samples).text)

    if (!cfg.annotation) { error "Provide 'annotation' in ${params.samples}" }

    annotation   = file(cfg.annotation)
    ref          = cfg.ref          ? file(cfg.ref)         : file('NO_REF')
    include_bed  = cfg.include_bed  ? file(cfg.include_bed) : file('NO_BED')
    // '' rather than null: Nextflow forbids a null value on a process `val` input.
    region       = cfg.region       ?: ''
    min_coverage = cfg.min_coverage ?: 10
    mod_code     = cfg.mod_code     ?: "m"
    min_delta    = cfg.min_delta    ?: 0.3
    z_threshold  = cfg.z_threshold  ?: 2.0

    samples_ch = Channel
        .fromList(cfg.samples)
        .map { s -> tuple(s.id, file(s.bam), findBamIndex(s.bam)) }

    MODKIT(samples_ch, annotation, ref, include_bed, region, min_coverage, mod_code)

    // Collect all "id path" pairs into one list, then run the cohort process once.
    sample_args = MODKIT.out.islands
        .map { id, tsv -> "${id}=${tsv}" }
        .collect()

    COHORT(sample_args, annotation, min_delta, z_threshold)
}
