# Hercules-NF Phase I — pytest validation suite.
# Exercises feature lookup, CIGAR parsing, read partitioning, motif aggregation
# and primer exclusion using small, deterministic SAM/GTF fixtures.
# Run from the repository root using: python -m pytest tests/test_phase1.py -q

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin' / 'hercules_phase1.py'


# Execute the Phase I command-line script with the supplied arguments.
# Capture its output and fail immediately if the command exits unsuccessfully.
def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                          check=True, capture_output=True, text=True)


# Verify that equal coordinates on different chromosomes remain distinct.
# Also check unknown-chromosome handling and GC percentages with ambiguous bases.
def test_overlap_index_distinguishes_chromosomes():
    sys.path.insert(0, str(ROOT / 'bin'))
    from hercules_phase1 import FeatureIndex, matching_blocks, gc_percent
    index = FeatureIndex(ROOT / 'data/example/features.gtf', 'CDS')
    blocks = matching_blocks(1,'12M',12)
    assert index.best('chr1', blocks) == ('chr1',1,330,'+')
    assert index.best('chr2', blocks) == ('chr2',1,330,'+')
    assert index.best('unknown', blocks) is None
    assert gc_percent('GGAN') == pytest.approx(66.6666666666667)


# Check reference/query coordinate conversion for splices, clipping and insertions.
# Confirm adjacent match operators are merged and malformed CIGAR strings rejected.
def test_cigar_splice_and_insert():
    sys.path.insert(0, str(ROOT / 'bin'))
    from hercules_phase1 import matching_blocks
    assert matching_blocks(11,'4M10N4M',8) == [(11,14,0,4),(25,28,4,8)]
    assert matching_blocks(11,'2S4M1I3M',10) == [(11,14,2,6),(15,17,7,10)]
    assert matching_blocks(11,'2M2=2X',6) == [(11,16,0,6)]
    with pytest.raises(ValueError):
        matching_blocks(11,'4M1Z',5)


# Run the complete Phase I path on a small SAM/GTF calibration fixture.
# Verify shard reduction, chromosome-specific GC records and motif-position totals.
def test_full_phase1_counts_adjacent_motifs(tmp_path):
    reads = ROOT / 'data/example/reads.sam'
    gtf = ROOT / 'data/example/features.gtf'
    shard_dir = tmp_path / 'shards'
    run('split','--sam',reads,'--partitions',3,'--output-dir',shard_dir)
    dbs = []
    for file in sorted(shard_dir.glob('*.sam')):
        db = tmp_path / f'{file.stem}.sqlite'
        run('map','--sam',file,'--gtf',gtf,'--k',4,'--primer-exclude',0,
            '--feature','CDS','--output',db,'--flush-entries',2)
        dbs.append(db)
    assert len(dbs)==3
    out = tmp_path / 'GC-content.csv'
    metrics = tmp_path / 'metrics.json'
    run('gc-reduce','--dbs',*dbs,'--output',out,'--metrics',metrics)
    with out.open() as f:
        features = dict(csv.reader(f))
    assert set(features)=={'chr1:1-330:+','chr2:1-330:+'}
    output = tmp_path / 'herc-final-AGCT.csv'
    run('motif-reduce','--dbs',*dbs,'--motif','AGCT','--output',output)
    with output.open() as f:
        counts = {(a,int(b)):int(c) for a,b,c in csv.reader(f)}
    for feature in features:
        assert counts[(feature,0)] == 2
        assert counts[(feature,4)] == 2
        assert counts[(feature,8)] == 2
    assert len(counts) == 6


# Check that overlapping four-base motifs are both counted in a short read.
# Re-run with primer exclusion to verify that only eligible starts remain.
def test_overlapping_motifs_and_primer_exclusion(tmp_path):
    gtf = tmp_path/'x.gtf'
    gtf.write_text('chrX\ttest\tCDS\t1\t50\t.\t+\t0\tgene_id "g";\n')
    sam = tmp_path/'x.sam'
    sam.write_text('read\t0\tchrX\t1\t60\t5M\t*\t0\t0\tAAAAA\t*\n')
    raw = tmp_path/'raw.sqlite'
    run('map','--sam',sam,'--gtf',gtf,'--k',4,'--primer-exclude',0,'--output',raw)
    with sqlite3.connect(raw) as conn:
        records = conn.execute('SELECT feature,motif,pos,n FROM counts ORDER BY pos').fetchall()
    assert [r[2] for r in records] == [0,1]  # both overlapping windows
    filtered = tmp_path/'filtered.sqlite'
    run('map','--sam',sam,'--gtf',gtf,'--k',4,'--primer-exclude',1,'--output',filtered)
    with sqlite3.connect(filtered) as conn:
        rows = conn.execute('SELECT pos FROM counts').fetchall()
    assert rows == [(1,)]
