# Hercules-NF — motif-position analysis on Nextflow DSL2 / SLURM

A two-program, two-phase reimplementation of the **Hercules** transcriptomic
motif-bias analysis described in Chapter 6 of the 2018 Ph.D. thesis
*Bias detection in NGS data using sequence motifs in exons*.

**Status:** Executable Python 3 implementation with local integration tests.
This repository is a **optimised reimplementation**, **not** a numerically
validated reproduction of all results in the 2018 thesis. The Nextflow DSL2
workflow is provided but **was not run with Nextflow/SLURM** in the development
environment. Validate it on your own Nextflow-enabled installation before a
production analysis.

## Design

```text
                       main.nf (Nextflow DSL2)
                                 |
             SPLIT_SAM (1 task: SAM -> N partitions)
                                 |
            MAP_SHARD (N concurrent SLURM tasks)
     Python Phase I: GTF assignment + per-shard SQLite
                    /                         \
         GC_SUMMARY (1 task)       MOTIF_REDUCE (4^k tasks)
           GC-content.csv           herc-final-MOTIF.csv
                    \                         /
                     complete all required data
                                 |
                    ANALYSE (1 SLURM task)
    Python Phase II: correlations, GC comparisons, FDR, plots
```

Only two scientific Python executables are used:

* `bin/hercules_phase1.py`: `split`, `map`, `gc-reduce`, and `motif-reduce`.
* `bin/hercules_phase2.py`: per-motif/per-spacing correlations, GC-stratified
  comparisons, outliers, plots, and an analysis summary.

Each SLURM job is a normal single-node process; there are no Spark executors,
HDFS transfers, MPI dependencies, or Nextflow-internal sequence algorithms.
A read partition is processed **once** to count all `4^k` motifs, and the
per-motif reduction tasks query indexed shard databases in parallel.

## Requirements

* Nextflow >= 23.04, Java supported by that Nextflow release
* Python >= 3.10
* NumPy, SciPy and Matplotlib (versions in `environment.yml`)
* SLURM, for the `slurm` profile (not required for `local`)
* A **shared filesystem** for Nextflow work files and all SLURM compute nodes
* Plain SAM or gzip-compressed SAM, plus a matching GTF annotation

No BAM or CRAM support is included. To use BAM, convert to SAM first, e.g.
`samtools view -h alignments.bam > alignments.sam`. Large production runs may
benefit from a dedicated BAM reader rather than materialising an uncompressed
SAM file.

## Quick start — local

```bash
# Create the Python environment (or use -profile local,conda)
conda env create -f environment.yml
conda activate hercules-nf

pytest -q

nextflow run . -profile local \
  --sam data/example/reads.sam \
  --gtf data/example/features.gtf \
  --feature CDS --k 4 --shards 3 \
  --primer_exclude 0 --min_pairs 2 \
  --outdir results/example -resume
```

Because the included example is tiny, few, if any, Pearson coefficients
will be defined — there are insufficient non-constant count pairs. This is
expected and is **not** evidence of a workflow failure.

For a standalone, **Nextflow-independent** full integration smoke test:

```bash
bash tests/smoke_full.sh
```

## SLURM

```bash
# Run on a login/submission host with Nextflow, Java and sbatch available.
# All worker nodes must see the repo (including bin/) and the work directory.
nextflow run /shared/projects/Hercules-NF \
  -profile slurm,conda \
  -w /shared/scratch/hercules_nf_work \
  --sam /shared/data/sample.sam.gz \
  --gtf /shared/reference/annotation.gtf \
  --feature CDS --k 4 --shards 32 \
  --primer_exclude 10 \
  --outdir /shared/results/hercules/sample \
  -resume
```

If the cluster uses an existing Python environment, omit `,conda` and
activate the environment **for compute tasks**, e.g. through the site's
module configuration or a suitable SLURM profile. Activation only in the
Nextflow submit-shell may not propagate to compute nodes.

Site-specific SLURM partition, account, reservation, wall time and memory
limits can be configured in `conf/slurm-site.config` (see below):

```bash
nextflow run . -profile slurm -c conf/slurm-site.config \
  -w /shared/scratch/hercules_nf_work \
  --sam /shared/data/sample.sam.gz \
  --gtf /shared/reference/annotation.gtf
```

Resource requirements in `nextflow.config` are starting points, **not
measured production allocations**. Each `MAP_SHARD` reads the full GTF and
maintains an indexed SQLite DB; adjust resources and shard count for your
actual datasets. `-resume` allows successful Nextflow tasks to be reused.

## Configurable arguments

| Argument | Default | Description |
| --- | --- | --- |
| `--sam` | required | SAM or SAM.gz file |
| `--gtf` | required | Matching genomic features in GTF |
| `--feature` | `CDS` | GTF feature type: `CDS` or `exon` |
| `--k` | `4` | k-mer size (2–5 supported by workflow) |
| `--shards` | `16` | Number of parallel read partitions |
| `--primer_exclude` | `10` | Discard k-mer occurrences starting before this 0-based read offset; `0` disables |
| `--min_mapq` | `0` | Ignore reads with lower MAPQ |
| `--exclude_duplicates` | `false` | Drop SAM-flagged duplicate reads |
| `--min_pairs` | `10` | Minimum eligible motif-position pairs for Pearson r |
| `--outlier_count` | `10` | Highest/lowest motif count correlations to report |
| `--save_intermediates` | `false` | Publish large read shards and SQLite databases |
| `--outdir` | `results` | Publication directory |

Spacings and tolerances are currently the Chapter 6 defaults:
10/50 bp (±2 bp), 100/200 bp (±4 bp). The four mean-read-GC bins are
30–40, 40–50, 50–60 and 60–70 percent.

## Outputs

```text
<outdir>/
  phase1/
    GC-content.csv            feature,mean-read-GC (no header)
    mapping_metrics.json      filtering and assignment statistics
    motifs/
      herc-final-AAAA.csv     feature,zero-based-position,count (no header)
      ...
      herc-final-TTTT.csv
  phase2/
    analysis/
      correlations.csv        one record per motif / spacing / GC bin
      gc_comparisons.csv      Welch and Mann–Whitney U p-values + BH FDR
      outliers.csv            extreme defined correlations
      summary.json            methods and summary metadata
      tables.tex              insertable LaTeX tables
      plots/
        correlations_10bp.png
        exon_gc_10bp.png
        ...
```

The intermediate position/count files deliberately retain the **three-field
CSV convention** of Hercules (`feature,position,count`), allowing the original
graphing algorithms to be compared against the new implementation once adapted
to the new feature identifiers.


## Scientific design choices and methodological optimisations

Hercules-NF builds upon `Hercules.py`, `Hercules-graphit.py` and Chapter 6 of the Ph.D. thesis, preserving the **core scientific methodology** while enhancing computational efficiency, genomic data handling, statistical robustness and reproducibility. The original implementations are retained under `legacy/`.

1. **Feature annotation.** Defaults to `CDS`, with optional exon analysis. GTF features are automatically sorted and assigned chromosome-specific identifiers (`chromosome:start-end:strand`). Identical features are consolidated, and overlaps are resolved deterministically by maximum aligned coverage.

2. **Read assignment.** Supports comprehensive CIGAR processing, including splice-aware alignments. Unmapped, secondary, supplementary and QC-failed reads are excluded; duplicate filtering is optional. Motifs spanning alignment discontinuities are excluded.

3. **Read orientation.** Motifs are searched in the SAM `SEQ` field in reference-aligned orientation. Primer exclusion uses sequence-relative positions, without additional strand-orientation conversion.

4. **Motif positions.** Detects all valid adjacent and overlapping k-mers within annotated features, using zero-based, feature-relative coordinates.

5. **GC measurement.** Calculates mean read-level GC percentage per annotated feature, excluding non-ACGT bases. This preserves the original read-derived GC metric.

6. **Noise filtering.** Retains per-motif first-quartile filtering using `numpy.percentile(method='higher')`.

7. **Correlations.** Calculates Pearson correlations from feature-associated motif pairs, stratified by spacing and GC content. Insufficient or undefined correlations are represented as `null`, rather than zero.

8. **Statistics.** Uses Welch's t-test and Mann–Whitney U for independent comparisons of motif-GC classes against the 50% reference group. Benjamini–Hochberg FDR correction is applied across valid comparisons. Groups require at least two observations.

9. **Biological replicates.** Analyses one SAM dataset per execution. Multi-sample comparisons and biological replicate-level inference are not currently implemented.

10. **Memory and I/O.** Uses SQLite-based aggregation across parallel read partitions, with exact quartile computation in Phase II. Shared-storage requirements and memory allocation should be benchmarked for large datasets.

**Scientific reproducibility.** Hercules-NF maintains the original research objectives while incorporating methodological optimisations. Numerical differences from historical Ph.D. results should be independently validated before drawing new biological conclusions.


## Verification

```bash
python -m py_compile bin/*.py
pytest -q
bash tests/smoke_full.sh
```

The automated test suite checks same-coordinate features on different
chromosomes, CIGAR blocks, adjacent/overlapping motifs, primer exclusion,
shard aggregation, mean GC, Pearson correlation, FDR and report generation.

The complete Nextflow execution graph and a real SLURM allocation still need
site-level validation. The Python smoke test does **not** establish that
Nextflow was run successfully.

