#!/usr/bin/env nextflow
/*
 * promethyl cohort pipeline
 *
 * Phase 1 (METHBAT, parallel per sample): MethBat pileup + island aggregation.
 * Phase 2 (COHORT, single process):       cross-sample outlier analysis.
 * Phase 3 (REPORT, single process):       render cohort TSV as HTML report.
 */

nextflow.enable.dsl = 2

params.samples = null
// YAML: samples, annotation, min_coverage, min_delta,
// z_threshold, dmr_bed, min_cpg_sites, panelapp_cache
params.outdir  = "results"
params.threads = 4

if (!params.samples) { error "Provide --samples cohort.yaml" }

// modkit pileup requires the BAM to already be indexed. Nextflow only stages files
// declared as process inputs into the task work dir, so the .bai/.csi sitting next
// to the BAM on the shared filesystem must be located and passed through explicitly.
def findBamIndex(bamPath) {
    def candidates = [bamPath + '.bai', bamPath + '.csi',
                       bamPath.replaceAll(/\.bam$/, '.bai'), bamPath.replaceAll(/\.bam$/, '.csi')]
    def found = candidates.collect { c -> file(c) }.find { f -> f.exists() }
    found ?: error("No BAM index (.bai/.csi) found for ${bamPath} -- run `samtools index ${bamPath}`")
}

process METHBAT {
    tag "$sample_id"

    publishDir "${params.outdir}/methbat",
        mode: 'copy',
        pattern: "*.5mC.bed.gz*"

    cpus params.threads

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path annotation
    path dmr_bed
    val min_coverage

    output:
    tuple val(sample_id), path("${sample_id}.islands.tsv"), emit: islands
    path "*.5mC.bed.gz"
    path "*.5mC.bed.gz.tbi"
    path "*.summary.json"

    script:
    """
    run_sample.py \
        --id ${sample_id} \
        --bam ${bam} \
        --methbat-dir . \
        --annotation ${annotation} \
        --output ${sample_id}.islands.tsv \
        --threads ${task.cpus} \
        --min-coverage ${min_coverage} \
        ${dmr_bed.name != 'NO_DMR' ? "--dmr-bed ${dmr_bed}" : ""}
    """
}


process COHORT {
    publishDir params.outdir, mode: 'copy'

    input:
    val sample_args
    val sample_meta_args
    path annotation
    val min_delta
    val z_threshold
    path dmr_bed
    val min_cpg_sites

    output:
    path "cohort_methylation.tsv"

    script:
    def sample_meta_opt = sample_meta_args ? "--sample-meta ${sample_meta_args.join(' ')}" : ""
    def dmr_bed_opt      = dmr_bed.name != 'NO_DMR' ? "--dmr-bed ${dmr_bed}" : ""
    """
    run_cohort.py \
        --samples ${sample_args.join(' ')} \
        ${sample_meta_opt} \
        --annotation ${annotation} \
        --output cohort_methylation.tsv \
        --min-delta ${min_delta} \
        --z-threshold ${z_threshold} \
        --min-cpg-sites ${min_cpg_sites} \
        ${dmr_bed_opt}
    """
}

process REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path cohort_tsv
    path panelapp_cache, stageAs: 'seed_cache.json'
    val has_panelapp_cache
    path report_script
    path panelapp_script

    output:
    path "report.html"
    path "panelapp_cache.json"

    script:
    // has_panelapp_cache is computed in the workflow block, against the *original* file
    // object -- not checked here via panelapp_cache.name, because inside `script:` that
    // variable refers to the staged/renamed local copy (always 'seed_cache.json' thanks to
    // stageAs below), which would make this check always true regardless of what was
    // actually passed in.
    //
    // report_script/panelapp_script are declared as real path inputs (not just referenced
    // via ${projectDir}) so their *contents* are part of this task's cache hash -- Nextflow
    // only hashes the literal script: block text plus its declared inputs, so a path merely
    // interpolated into a shell command (as this used to be) is invisible to -resume: editing
    // generate_report.py previously never invalidated the cache at all. Staged under their
    // own basenames (default stageAs) so `from panelapp import ...` inside generate_report.py
    // still resolves via same-directory sys.path, same as when it ran from the repo directly.
    """
    ${has_panelapp_cache ? "cp seed_cache.json panelapp_cache.json" : "echo '{}' > panelapp_cache.json"}
    python3 ${report_script} \
        --input ${cohort_tsv} \
        --output report.html \
        --panelapp-cache panelapp_cache.json \
        --outliers-only
    """
}


workflow {
    cfg = new org.yaml.snakeyaml.Yaml().load(file(params.samples).text)

    if (!cfg.annotation) {
        error "Provide 'annotation' in ${params.samples}"
    }

    annotation = file(cfg.annotation)

    min_coverage = cfg.min_coverage ?: 10
    min_delta = cfg.min_delta ?: 0.3
    z_threshold = cfg.z_threshold ?: 2.0

    dmr_bed = cfg.dmr_bed ? file(cfg.dmr_bed) : file('NO_DMR')
    min_cpg_sites = cfg.min_cpg_sites ?: 1

    panelapp_cache = (
        cfg.panelapp_cache &&
        file(cfg.panelapp_cache).exists()
    ) ? file(cfg.panelapp_cache) : file('NO_CACHE')

    has_panelapp_cache = panelapp_cache.name != 'NO_CACHE'

    samples_ch = Channel
        .fromList(cfg.samples)
        .map { s ->
            tuple(s.id, file(s.bam), findBamIndex(s.bam))
        }

    METHBAT(
        samples_ch,
        annotation,
        dmr_bed,
        min_coverage
    )

    sample_args = METHBAT.out.islands
        .map { id, tsv -> "${id}=${tsv}" }
        .collect()

    def hasSampleMeta = cfg.samples.any {
        s -> s.tissue || s.sex
    }

    sample_meta_args = hasSampleMeta
        ? cfg.samples.collect {
            s -> "${s.id}=${s.tissue ?: 'NA'}=${s.sex ?: 'NA'}"
        }
        : []

    COHORT(
        sample_args,
        sample_meta_args,
        annotation,
        min_delta,
        z_threshold,
        dmr_bed,
        min_cpg_sites
    )

    REPORT(
        COHORT.out,
        panelapp_cache,
        has_panelapp_cache,
        file("${projectDir}/src/generate_report.py"),
        file("${projectDir}/src/panelapp.py")
    )
}

