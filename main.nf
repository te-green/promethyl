#!/usr/bin/env nextflow
/*
 * promethyl cohort pipeline
 *
 * Phase 1 (MODKIT, parallel per sample): modkit pileup + island aggregation.
 * Phase 2 (COHORT, single process):      cross-sample outlier analysis.
 * Phase 3 (REPORT, single process):      render cohort TSV as HTML report.
 */

nextflow.enable.dsl = 2

params.samples = null                     // YAML: samples, plus optional annotation/ref/include_bed/
                                           // region/min_coverage/mod_code/min_delta/z_threshold/
                                           // dmr_bed/min_cpg_sites/panelapp_cache
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

// PacBio BAMs sometimes mix reads that carry MM/ML modification-call tags with reads
// that don't (e.g. low-pass reads the kinetics caller couldn't confidently call). modkit
// pileup can hang when it encounters that mix. This is a cheap check (seconds, tiny
// memory): sample the first 1000 records and report whether any lack MM, so the workflow
// can route only the affected samples through the expensive FIX_MM_TAGS process.
process CHECK_MM_TAGS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path(bam), path(bam_index), env(NEEDS_FIX), emit: checked

    script:
    """
    untagged=\$(samtools view ${bam} | head -n 1000 | awk -F'\\t' '{ok=0; for(i=12;i<=NF;i++) if (\$i ~ /^MM:/) {ok=1; break}; if (!ok) print}' | wc -l)

    NEEDS_FIX=false
    if [ "\$untagged" -gt 0 ]; then
        echo "[${sample_id}] \$untagged/1000 sampled reads missing MM tag -- needs fixing"
        NEEDS_FIX=true
    else
        echo "[${sample_id}] all sampled reads carry MM tag -- no fix needed"
    fi
    """
}

// Reads with neither MM nor ML are stamped with an empty "no modified bases called"
// MM/ML pair (C+m?;  /  B:C with no values -- note no trailing comma before ';', which
// would encode one empty delta entry instead of zero and leave MM/ML mismatched in
// length) so they still contribute canonical-base coverage instead of being dropped.
// Reads that already carry MM/ML pass through untouched. The two halves are then merged
// back together and re-sorted, which is the expensive part of this process -- give
// `samtools sort` most of the task's memory budget so it can use big in-memory chunks
// instead of spilling to disk.
process FIX_MM_TAGS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path("out/${sample_id}.bam"), path("out/${sample_id}.bam.bai"), emit: bam

    script:
    def sort_mem = (task.memory.toGiga() * 0.8 / task.cpus) as int
    """
    mkdir out
    samtools view -h -d MM -b -o with_mm.bam -U without_mm.bam ${bam}
    samtools view -h without_mm.bam \
        | awk -F'\\t' 'BEGIN{OFS="\\t"} /^@/{print; next} {print \$0, "MM:Z:C+m?;", "ML:B:C"}' \
        | samtools view -b -o tagged.bam -
    samtools cat with_mm.bam tagged.bam \
        | samtools sort -m ${sort_mem}G -@ ${task.cpus} -o out/${sample_id}.bam -
    samtools index out/${sample_id}.bam
    """
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

    output:
    path "report.html"
    path "panelapp_cache.json"

    script:
    // has_panelapp_cache is computed in the workflow block, against the *original* file
    // object -- not checked here via panelapp_cache.name, because inside `script:` that
    // variable refers to the staged/renamed local copy (always 'seed_cache.json' thanks to
    // stageAs below), which would make this check always true regardless of what was
    // actually passed in.
    """
    ${has_panelapp_cache ? "cp seed_cache.json panelapp_cache.json" : "echo '{}' > panelapp_cache.json"}
    python3 ${projectDir}/src/generate_report.py \
        --input ${cohort_tsv} \
        --output report.html \
        --panelapp-cache panelapp_cache.json
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
    dmr_bed       = cfg.dmr_bed       ? file(cfg.dmr_bed) : file('NO_DMR')
    min_cpg_sites = cfg.min_cpg_sites ?: 1
    // .exists() guards against a configured path that isn't there yet (e.g. a fresh run
    // directory pointing at a not-yet-created results/panelapp_cache.json from a future
    // run); computed here, before staging, since checking it inside the REPORT script
    // block would see the staged/renamed local filename instead of this original path.
    panelapp_cache     = (cfg.panelapp_cache && file(cfg.panelapp_cache).exists()) ? file(cfg.panelapp_cache) : file('NO_CACHE')
    has_panelapp_cache = panelapp_cache.name != 'NO_CACHE'

    samples_ch = Channel
        .fromList(cfg.samples)
        .map { s -> tuple(s.id, file(s.bam), findBamIndex(s.bam)) }

    CHECK_MM_TAGS(samples_ch)

    checked = CHECK_MM_TAGS.out.checked.branch { id, bam, bai, needs_fix ->
        needs_fix: needs_fix == 'true'
            return tuple(id, bam, bai)
        ok: true
            return tuple(id, bam, bai)
    }

    FIX_MM_TAGS(checked.needs_fix)

    modkit_input_ch = checked.ok.mix(FIX_MM_TAGS.out.bam)

    MODKIT(modkit_input_ch, annotation, ref, include_bed, region, min_coverage, mod_code)

    // Collect all "id path" pairs into one list, then run the cohort process once.
    sample_args = MODKIT.out.islands
        .map { id, tsv -> "${id}=${tsv}" }
        .collect()

    // tissue/sex come straight from the YAML (not the channel) since they're only
    // needed by the COHORT process, not MODKIT. Missing tissue/sex on a sample
    // becomes "NA" -- run_cohort.py excludes such samples from every comparison
    // group rather than guessing, so it's safe but that sample won't get tested.
    def hasSampleMeta = cfg.samples.any { s -> s.tissue || s.sex }
    sample_meta_args = hasSampleMeta
        ? cfg.samples.collect { s -> "${s.id}=${s.tissue ?: 'NA'}=${s.sex ?: 'NA'}" }
        : []

    COHORT(sample_args, sample_meta_args, annotation, min_delta, z_threshold, dmr_bed, min_cpg_sites)

    REPORT(COHORT.out, panelapp_cache, has_panelapp_cache)
}
