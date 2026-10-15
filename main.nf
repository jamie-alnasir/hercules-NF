#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/*
 * ============================================================================
 * Hercules-NF — Nextflow DSL2 Transcriptomic Motif Analysis
 * ============================================================================
 *
 * INPUT PARAMETERS
 *
 * Required:
 *   --sam                 Input aligned sequencing reads (SAM or SAM.gz).
 *   --gtf                 Genome annotation file (GTF).
 *
 * Optional:
 *   --outdir              Output directory (default: results).
 *   --feature             GTF feature type: CDS or exon (default: CDS).
 *   --k                   Motif length, 2–5 bases (default: 4).
 *   --shards              Number of parallel SAM partitions (default: 16).
 *   --primer_exclude      Read-start bases excluded from motif analysis
 *                         (default: 10; 0 disables exclusion).
 *   --min_mapq            Minimum mapping quality (default: 0).
 *   --exclude_duplicates  Exclude SAM-flagged duplicates (default: false).
 *   --min_pairs           Minimum pairs for Pearson correlation (default: 10).
 *   --outlier_count       Number of extreme correlations reported (default: 10).
 *   --save_intermediates  Retain intermediate SAM/SQLite files (default: false).
 *
 * EXECUTION
 *
 * Local:
 *   nextflow run main.nf -profile local \
 *       --sam reads.sam --gtf annotations.gtf
 *
 * SLURM:
 *   nextflow run main.nf -profile slurm,conda \
 *       --sam reads.sam --gtf annotations.gtf \
 *       --k 4 --shards 32 --outdir results -resume
 *
 * Nextflow options:
 *   -profile      Select execution profile (local, slurm, conda).
 *   -resume       Reuse completed tasks from previous executions.
 *   -w            Specify the shared working directory.
 *   -c            Load an additional configuration file.
 *
 * Phase I:  Parallel read assignment, motif counting and GC aggregation.
 * Phase II: Statistical correlations, GC analysis and report generation.
 *
 * Parameters and defaults are defined in nextflow.config.
 * ============================================================================
 */

/*
 * Hercules-NF: Chapter 6 motif-count analysis
 * Two scientific phases:
 *  - hercules_phase1.py: split, map, GC reduce, per-motif reduce
 *  - hercules_phase2.py: offline statistical interpretation and figures
 * All computing jobs are normal SLURM (or local) process tasks; no Spark/HDFS.
 */

// Divide the input SAM into independently processable read partitions.
// Each shard can subsequently be analysed as a separate SLURM task.
// The named `shards` output contains the generated SAM files.
process SPLIT_SAM {
    tag 'partition SAM'
    publishDir "${params.outdir}/intermediate/shards", mode: 'copy', enabled: params.save_intermediates

    input:
    path sam

    output:
    path 'shards/shard_*.sam', emit: shards

    script:
    """
    mkdir -p shards
    hercules_phase1.py split --sam '${sam}' \
        --partitions ${params.shards} --output-dir shards
    """
}

// Analyse a single SAM partition against the selected GTF annotations.
// Assign reads to features, count motif occurrences and accumulate GC values
// in a shard-specific SQLite database, retaining its shard identifier.
process MAP_SHARD {
    tag "${shard_id}"
    publishDir "${params.outdir}/intermediate/databases", mode: 'copy', enabled: params.save_intermediates

    input:
    tuple val(shard_id), path(shard), path(gtf)

    output:
    tuple val(shard_id), path("${shard_id}.sqlite"), emit: db

    script:
    def duplicateFlag = params.exclude_duplicates ? '--exclude-duplicates' : ''
    """
    hercules_phase1.py map \
        --sam '${shard}' --gtf '${gtf}' \
        --output '${shard_id}.sqlite' \
        --feature '${params.feature}' --k ${params.k} \
        --primer-exclude ${params.primer_exclude} \
        --min-mapq ${params.min_mapq} ${duplicateFlag}
    """
}

// Merge the feature-level GC statistics and mapping metrics from all shards.
// The `gc` and `metrics` emit names expose two separate output channels;
// Phase II uses the GC table to stratify motif-count correlations.
process GC_SUMMARY {
    tag 'feature mean read GC'
    publishDir "${params.outdir}/phase1", mode: 'copy'

    input:
    path dbs

    output:
    path 'GC-content.csv', emit: gc
    path 'mapping_metrics.json', emit: metrics

    script:
    """
    hercules_phase1.py gc-reduce \
        --dbs *.sqlite \
        --output GC-content.csv \
        --metrics mapping_metrics.json
    """
}

// Reduce counts for one nucleotide motif across all shard databases.
// Nextflow invokes this process independently for each of the 4^k motifs,
// emitting one feature-position count CSV per completed task.
process MOTIF_REDUCE {
    tag "motif ${motif}"
    publishDir "${params.outdir}/phase1/motifs", mode: 'copy'

    input:
    tuple val(motif), path(dbs)

    output:
    path "herc-final-${motif}.csv", emit: csv

    script:
    """
    hercules_phase1.py motif-reduce \
        --dbs *.sqlite \
        --motif '${motif}' \
        --output 'herc-final-${motif}.csv'
    """
}

// Run the non-distributed Phase II statistical interpretation program.
// Receive the feature-GC table and collected motif CSVs, then generate
// correlation statistics, comparative tests, figures and LaTeX tables.
process ANALYSE {
    tag 'Hercules Phase II'
    publishDir "${params.outdir}/phase2", mode: 'copy'

    input:
    tuple path(gc_file), path(motif_files)

    output:
    path 'analysis', emit: report

    script:
    """
    hercules_phase2.py \
        --counts-dir . --gc '${gc_file}' \
        --output-dir analysis \
        --k ${params.k} --min-pairs ${params.min_pairs} \
        --outlier-count ${params.outlier_count}
    """
}

// Assemble the two-phase pipeline and its task dependencies.
// Scatter read processing across shards, gather counts by motif, and wait
// for the GC and motif outputs before starting the final analysis.
workflow {
    if (!params.sam || !params.gtf) {
        error 'Specify both --sam and --gtf (see README.md)'
    }
    if ((params.k as int) < 2 || (params.k as int) > 5) {
        error '--k must be 2, 3, 4, or 5'
    }
    if ((params.shards as int) < 1) {
        error '--shards must be >= 1'
    }
    if (!(params.feature in ['CDS', 'exon'])) {
        error '--feature must be CDS or exon'
    }

    // Enumerate all A/C/G/T motifs of length k (4^k combinations).
    def bases = ['A', 'C', 'G', 'T']
    def motifs = ['']
    (1..(params.k as int)).each { n ->
        motifs = motifs.collectMany { prefix ->
            bases.collect { base -> prefix + base }
        }
    }

    def sam = file(params.sam, checkIfExists: true)
    def gtf = file(params.gtf, checkIfExists: true)

    // Launch the partitioning step with the resolved SAM input path.
    SPLIT_SAM(Channel.value(sam))

    // Expand the shard-file collection into individual (ID, SAM, GTF) tuples.
    // Each resulting tuple triggers one independent MAP_SHARD task.
    shard_inputs = SPLIT_SAM.out.shards
        .flatten()
        .map { shard -> tuple(shard.baseName, shard, gtf) }

    MAP_SHARD(shard_inputs)

    // collect waits for ALL shard jobs; read-only shard databases can then
    // be queried simultaneously by the GC and motif reduction tasks.
    databases = MAP_SHARD.out.db.map { shard_id, db -> db }.collect()

    // Reuse the collected shard databases to calculate overall GC summaries.
    GC_SUMMARY(databases)

    // Associate each motif with the shard-database collection for reduction.
    // Each emitted motif tuple schedules a separate MOTIF_REDUCE task.
    motif_inputs = databases
        .combine(Channel.fromList(motifs))
        .map { dbs, motif -> tuple(motif, dbs) }

    // Execute motif-specific reductions concurrently as resources permit.
    MOTIF_REDUCE(motif_inputs)

    // Synchronisation barrier: Phase II executes once all k-mer CSVs exist.
    analysis_inputs = GC_SUMMARY.out.gc
        .combine(MOTIF_REDUCE.out.csv.collect())

    ANALYSE(analysis_inputs)
}
