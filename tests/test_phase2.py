# Hercules-NF Phase II — pytest validation suite.
# Checks Pearson accumulation, Benjamini–Hochberg correction, spacing/GC
# stratification, report generation and handling of absent GC annotations.
# Run from the repository root using: python -m pytest tests/test_phase2.py -q

import csv
import json
import subprocess
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin' / 'hercules_phase2.py'
sys.path.insert(0, str(ROOT / 'bin'))
from hercules_phase2 import CorrelationAccumulator, bh_fdr, calculate_motif, group_comparisons


# Verify streaming Pearson correlation for perfectly associated count pairs.
# Confirm that too few pairs or zero variance produce an undefined result.
def test_pearson_and_undefined():
    x = CorrelationAccumulator()
    for a,b in [(2,4),(3,6),(4,8)]:
        x.add(a,b)
    assert x.pearson(2) == pytest.approx(1.0)
    assert x.pearson(4) is None
    c = CorrelationAccumulator()
    c.add(3,4)
    c.add(3,8)
    assert c.pearson(2) is None


# Check Benjamini–Hochberg false-discovery-rate adjusted p-values.
# Include the empty input case to verify that no corrections are fabricated.
def test_bh():
    assert bh_fdr([0.01,0.04,0.03]) == pytest.approx([.03,.04,.04])
    assert bh_fdr([]) == []


# Check first-quartile filtering and the pairing of motifs 10 bp apart.
# Verify that observations enter the correct feature-associated GC stratum.
def test_ten_bp_pair_and_gc_bin(tmp_path):
    file = tmp_path/'herc-final-AAAA.csv'
    file.write_text('chr1:1-100:+,0,1\nchr1:1-100:+,10,2\n'
                    'chr1:1-100:+,20,3\nchr1:1-100:+,30,4\n')
    rows = calculate_motif('AAAA',file,{'chr1:1-100:+':45.0},[10],{10:2},2)
    overall = next(x for x in rows if x['gc_bin']=='all')
    assert overall['quartile']==2
    assert overall['pairs']==2
    assert overall['pearson_r'] == pytest.approx(1.0)
    assert next(x for x in rows if x['gc_bin']=='40-50')['pairs']==2
    assert next(x for x in rows if x['gc_bin']=='30-40')['pairs']==0


# Run the Phase II command-line program on synthetic motif and GC tables.
# Verify its summary, statistical CSVs, LaTeX output and expected plot files.
def test_phase2_creates_csv_and_figures(tmp_path):
    folder = tmp_path/'counts'
    folder.mkdir()
    for name in ('AAAA','AGCT'):
        (folder/f'herc-final-{name}.csv').write_text(
            'chr1:1-100:+,0,1\nchr1:1-100:+,10,2\n'
            'chr1:1-100:+,20,3\nchr1:1-100:+,30,4\n')
    gc = tmp_path/'GC-content.csv'
    gc.write_text('chr1:1-100:+,45\n')
    out = tmp_path/'analysis'
    subprocess.run([sys.executable,str(SCRIPT),'--counts-dir',str(folder),
                    '--gc',str(gc),'--output-dir',str(out),'--k','4',
                    '--min-pairs','2'],check=True,capture_output=True)
    summary = json.loads((out/'summary.json').read_text())
    assert summary['motifs_analysed']==2
    assert (out/'correlations.csv').exists()
    assert (out/'gc_comparisons.csv').exists()
    assert (out/'tables.tex').exists()
    assert len(list((out/'plots').glob('*.png'))) == 8


# Ensure motif-count features must have corresponding GC annotations.
# Verify that missing feature metadata triggers a clear validation error.
def test_missing_feature_gc_fails(tmp_path):
    file = tmp_path/'herc-final-AAAA.csv'
    file.write_text('featureX,0,1\nfeatureX,10,2\nfeatureX,20,3\n')
    with pytest.raises(ValueError,match='Missing GC'):
        calculate_motif('AAAA',file,{},[10],{10:2},2)
