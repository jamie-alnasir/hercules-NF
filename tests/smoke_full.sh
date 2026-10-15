#!/usr/bin/env bash
# Hercules-NF: end-to-end smoke test of both Python analysis phases.
# Runs the bundled synthetic SAM/GTF dataset through read partitioning,
# motif counting, GC aggregation and Phase II reporting without Nextflow/SLURM.
# Verifies all 256 four-mer outputs and eight figures; removes temporary data.

# Stop on errors, unset variables or failed pipeline commands.
set -euo pipefail

# Locate the repository and create an isolated workspace cleaned on exit.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# Phase I: split the example SAM into three independent read partitions.
python "$ROOT/bin/hercules_phase1.py" split \
  --sam "$ROOT/data/example/reads.sam" --partitions 3 --output-dir "$TMP/shards"

# Map each partition to CDS features and count its four-mer observations.
# Store intermediate motif counts and read-GC summaries in shard SQLite files.
for shard in "$TMP"/shards/*.sam; do
  stem="$(basename "$shard" .sam)"
  python "$ROOT/bin/hercules_phase1.py" map \
    --sam "$shard" --gtf "$ROOT/data/example/features.gtf" \
    --k 4 --feature CDS --primer-exclude 0 --output "$TMP/${stem}.sqlite"
done

# Merge feature-associated GC measurements and mapping statistics from shards.
python "$ROOT/bin/hercules_phase1.py" gc-reduce --dbs "$TMP"/*.sqlite \
  --output "$TMP/GC-content.csv" --metrics "$TMP/mapping_metrics.json"

# Reduce all 256 possible four-mers, merging their counts across shard databases.
# This smoke test uses one Python interpreter for speed; Nextflow schedules
# separate motif-reduction tasks when running the full workflow.
mkdir -p "$TMP/motifs"
python - "$ROOT" "$TMP" <<'PYREDUCE'
from argparse import Namespace
from itertools import product
from pathlib import Path
import sys

# Import the Phase I reducer from the repository's executable directory.
root, temp = map(Path, sys.argv[1:])
sys.path.insert(0, str(root / 'bin'))
from hercules_phase1 import motif_reduce

# Generate every A/C/G/T four-mer and write its aggregated position-count CSV.
dbs = sorted(temp.glob('*.sqlite'))
for letters in product('ACGT', repeat=4):
    motif = ''.join(letters)
    motif_reduce(Namespace(dbs=dbs, motif=motif,
                           output=str(temp / 'motifs' / f'herc-final-{motif}.csv')))
PYREDUCE

# Phase II: calculate motif-spacing correlations, GC comparisons and plots.
# Use a reduced minimum-pair threshold appropriate for the small test dataset.
python "$ROOT/bin/hercules_phase2.py" --counts-dir "$TMP/motifs" \
  --gc "$TMP/GC-content.csv" --output-dir "$TMP/analysis" --k 4 --min-pairs 2

# Confirm that Phase II processed every motif and produced the expected plots.
# Any failed assertion terminates the test with a non-zero exit status.
python - "$TMP" <<'PY'
from pathlib import Path
import json,sys

p=Path(sys.argv[1]); s=json.loads((p/'analysis/summary.json').read_text())
assert s['motifs_analysed'] == 256, s
assert len(list((p/'motifs').glob('herc-final-*.csv'))) == 256
assert len(list((p/'analysis/plots').glob('*.png'))) == 8
print('PASS: 256 motifs, per-feature GC, Phase II CSVs and 8 figures')
PY
