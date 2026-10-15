#!/usr/bin/env python3
"""Hercules-NF Phase II: Chapter 6-inspired k-mer count correlations and reports.

Consumes the original no-header 'feature,position,count' motif CSV contract and
'feature,mean_read_GC' table. Uses independent-sample
nonparametric tests instead of legacy arbitrary Wilcoxon pairing.
"""
from __future__ import annotations
import argparse
from bisect import bisect_left, bisect_right
from collections import defaultdict
import csv
from itertools import groupby
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GC_BINS = ((30, 40), (40, 50), (50, 60), (60, 70))
GC_CLASSES = (0, 25, 50, 75, 100)


# Accumulate motif-position pair statistics without retaining every pair.
# Track the sums required to calculate Pearson correlation efficiently.
class CorrelationAccumulator:
    """Exact integer sufficient statistics until calculation of Pearson r."""
    __slots__ = ('n', 'sx', 'sy', 'sxx', 'syy', 'sxy')

    # Initialise pair counts and the six running sufficient statistics.
    # Values are stored as integers until the final correlation is calculated.
    def __init__(self):
        self.n = self.sx = self.sy = self.sxx = self.syy = self.sxy = 0

    # Include one pair of observed motif counts in the running totals.
    # Update sample size, marginal sums, squared sums and cross-products.
    def add(self, x, y):
        self.n += 1
        self.sx += x
        self.sy += y
        self.sxx += x*x
        self.syy += y*y
        self.sxy += x*y

    # Compute Pearson r from the accumulated statistics once enough pairs exist.
    # Return None for inadequate data or zero variance in either variable.
    def pearson(self, min_pairs):
        if self.n < min_pairs or self.n < 2:
            return None
        vx = self.n*self.sxx-self.sx*self.sx
        vy = self.n*self.syy-self.sy*self.sy
        if vx <= 0 or vy <= 0:
            return None
        return float(self.n*self.sxy-self.sx*self.sy) / math.sqrt(vx*vy)


# Load feature-associated mean read-GC percentages from a two-column CSV.
# Index values by feature ID and reject malformed or duplicate entries.
def load_gc(file):
    result = {}
    with open(file, newline='') as source:
        for row in csv.reader(source):
            if not row:
                continue
            if len(row) != 2:
                raise ValueError('GC file must have exactly two columns')
            if row[0] in result:
                raise ValueError(f'Duplicate feature in GC file: {row[0]}')
            result[row[0]] = float(row[1])
    return result


# Stream motif counts as (feature ID, position, count) records from CSV.
# Convert numeric fields to integers and report invalid records.
def read_counts(file):
    with open(file, newline='') as source:
        for row in csv.reader(source):
            if not row:
                continue
            if len(row) != 3:
                raise ValueError(f'Expected feature,position,count in {file}: {row}')
            try:
                yield row[0], int(row[1]), int(row[2])
            except ValueError as error:
                raise ValueError(f'Invalid numeric count in {file}: {row}') from error


# Locate result files matching the requested k-mer length and alphabet.
# Return a motif-to-file mapping for subsequent statistical analysis.
def discover_files(folder, k):
    motif_files = sorted(Path(folder).glob('herc-final-*.csv'))
    result = {}
    for file in motif_files:
        motif = file.stem.removeprefix('herc-final-')
        if len(motif) != k or set(motif) - set('ACGT'):
            continue
        if motif in result:
            raise ValueError(f'Duplicate motif: {motif}')
        result[motif] = file
    if not result:
        raise ValueError(f'No herc-final-<motif>.csv (k={k}) files in {folder}')
    return result


# Calculate the first quartile of observed counts for one motif.
# Exclude unassigned features and retain the historical higher-percentile rule.
def first_quartile(file):
    # Historical method: numpy percentile with 'higher' interpolation.
    counts = [n for feature, pos, n in read_counts(file) if feature != '-1']
    if not counts:
        return None
    return float(np.percentile(np.asarray(counts, dtype=np.int64), 25, method='higher'))


# Pair feature-local motif counts at the requested spacings and tolerances.
# Apply quartile filtering, stratify by mean read-GC and calculate Pearson r.
def calculate_motif(motif, file, gc, spacings, tolerance, min_pairs):
    quartile = first_quartile(file)
    accum = {(d, label): CorrelationAccumulator()
             for d in spacings for label in ['all'] + [f'{lo}-{hi}' for lo, hi in GC_BINS]}
    if quartile is not None:
        # Phase I guarantees sorting by (feature,pos). groupby works only on sorted files.
        last_group = None
        for feature, group in groupby(read_counts(file), key=lambda rec: rec[0]):
            if last_group is not None and feature < last_group:
                raise ValueError(f'Count file not sorted by feature: {file}')
            last_group = feature
            accepted = sorted([(pos,n) for _,pos,n in group if n >= quartile])
            if len(accepted) < 2:
                continue
            positions = [p for p,n in accepted]
            gcv = gc.get(feature)
            if gcv is None:
                raise ValueError(f'Missing GC value for feature {feature}')
            bins = ['all'] + [f'{lo}-{hi}' for lo,hi in GC_BINS if lo <= gcv < hi]
            for i, (pos, x) in enumerate(accepted):
                for spacing in spacings:
                    tol = tolerance[spacing]
                    start = bisect_left(positions, pos+spacing-tol, lo=i+1)
                    stop = bisect_right(positions, pos+spacing+tol, lo=start)
                    for j in range(start, stop):
                        y = accepted[j][1]
                        for label in bins:
                            accum[spacing,label].add(x,y)
    gc_pct = 100.0*(motif.count('G')+motif.count('C'))/len(motif)
    rows = []
    for spacing in spacings:
        for label in ['all'] + [f'{lo}-{hi}' for lo,hi in GC_BINS]:
            acc = accum[spacing,label]
            rows.append({
                'motif':motif, 'motif_gc':gc_pct, 'quartile':quartile,
                'spacing':spacing, 'tolerance':tolerance[spacing],
                'gc_bin':label, 'pairs':acc.n,
                'pearson_r':acc.pearson(min_pairs)
            })
    return rows


# Apply the Benjamini–Hochberg false-discovery-rate adjustment.
# Return adjusted p-values in the same order as the supplied tests.
def bh_fdr(pvals):
    """Benjamini–Hochberg adjusted p-values for finite input values."""
    if not pvals:
        return []
    p = np.asarray(pvals, dtype=float)
    ix = np.argsort(p)
    sorted_p = p[ix]
    result = np.minimum.accumulate((sorted_p * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    ordered = np.empty_like(result)
    ordered[ix] = np.clip(result, 0, 1)
    return ordered.tolist()


# Compare motif-GC correlation distributions with the 50% GC reference class.
# Use Welch and Mann–Whitney tests within feature-GC bins and adjust p-values.
def group_comparisons(rows, spacings):
    """Compare motif-GC classes to 50% GC, per spacing and exon-GC stratum."""
    comparisons = []
    for spacing in spacings:
        for label in [f'{lo}-{hi}' for lo,hi in GC_BINS]:
            groups = defaultdict(list)
            for row in rows:
                if row['spacing']==spacing and row['gc_bin']==label and row['pearson_r'] is not None:
                    groups[row['motif_gc']].append(row['pearson_r'])
            control = groups[50.0]
            for other_gc in (0.0,25.0,75.0,100.0):
                vals = groups[other_gc]
                welch = mw = None
                if len(control) >= 2 and len(vals) >= 2:
                    with np.errstate(invalid='ignore', divide='ignore'):
                        result = stats.ttest_ind(control, vals, equal_var=False)
                    welch = float(result.pvalue) if np.isfinite(result.pvalue) else None
                    try:
                        result = stats.mannwhitneyu(control, vals, alternative='two-sided', method='auto')
                        mw = float(result.pvalue) if np.isfinite(result.pvalue) else None
                    except ValueError:
                        pass
                comparisons.append({
                    'spacing':spacing,'gc_bin':label, 'reference_motif_gc':50.0,
                    'comparison_motif_gc':other_gc,
                    'n_reference':len(control), 'n_comparison':len(vals),
                    'welch_p':welch,'mannwhitney_p':mw,
                    'welch_fdr':None,'mannwhitney_fdr':None
                })
    # One explicitly-defined testing family for every method, across all bins and spacings.
    for name, adjusted in (('welch_p','welch_fdr'), ('mannwhitney_p','mannwhitney_fdr')):
        indices = [i for i,r in enumerate(comparisons) if r[name] is not None]
        corrections = bh_fdr([comparisons[i][name] for i in indices])
        for i,v in zip(indices,corrections):
            comparisons[i][adjusted] = v
    return comparisons


# Rank the highest and lowest defined motif correlations for each spacing.
# Retain the requested number of motifs at each extreme for reporting.
def outliers(rows, spacings, number):
    result = []
    for spacing in spacings:
        candidates = [r for r in rows if r['spacing']==spacing and r['gc_bin']=='all'
                      and r['pearson_r'] is not None]
        asc = sorted(candidates, key=lambda r:(r['pearson_r'],r['motif']))
        for side, selected in [('lowest',asc[:number]), ('highest',list(reversed(asc[-number:])) )]:
            for rank, row in enumerate(selected,1):
                result.append({'spacing':spacing,'tail':side,'rank':rank,
                               'motif':row['motif'],'pearson_r':row['pearson_r'],
                               'pairs':row['pairs']})
    return result


# Write tabular analysis results using a fixed column order and CSV header.
# Ignore surplus dictionary fields not requested for the output table.
def write_csv(path, rows, columns):
    with open(path, 'w', newline='') as dest:
        writer = csv.DictWriter(dest, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


# Generate per-spacing box plots of correlations by motif and feature GC.
# Save headless Matplotlib figures, including when no valid values exist.
def plots(rows, spacings, destination):
    p = Path(destination) / 'plots'
    p.mkdir(parents=True,exist_ok=True)
    for spacing in spacings:
        fig, ax = plt.subplots(figsize=(9,5))
        by_gc = {gc:[r['pearson_r'] for r in rows if r['spacing']==spacing
                 and r['gc_bin']=='all' and r['motif_gc']==gc and r['pearson_r'] is not None]
                 for gc in GC_CLASSES}
        good = [(gc, values) for gc, values in by_gc.items() if values]
        if good:
            ax.boxplot([v for gc,v in good], tick_labels=[str(gc) for gc,v in good])
        else:
            ax.text(.5,.5,'No defined correlations',ha='center',va='center',transform=ax.transAxes)
        ax.set_xlabel('Motif GC content (%)')
        ax.set_ylabel('Pearson correlation (r)')
        ax.set_title(f'Intra-feature motif count correlation — {spacing} bp spacing')
        ax.set_ylim(-1.05,1.05)
        fig.tight_layout()
        fig.savefig(p / f'correlations_{spacing}bp.png',dpi=160)
        plt.close(fig)
        # Feature GC categories, pooled across motif sequences.
        fig, ax = plt.subplots(figsize=(9,5))
        bins = [f'{lo}-{hi}' for lo,hi in GC_BINS]
        vals = [[r['pearson_r'] for r in rows if r['spacing']==spacing
                 and r['gc_bin']==b and r['pearson_r'] is not None] for b in bins]
        useful = [(label,v) for label,v in zip(bins,vals) if v]
        if useful:
            ax.boxplot([v for label,v in useful],tick_labels=[label for label,v in useful])
        else:
            ax.text(.5,.5,'No defined correlations',ha='center',va='center',transform=ax.transAxes)
        ax.set_xlabel('Feature-associated mean read GC (%)')
        ax.set_ylabel('Pearson correlation (r)')
        ax.set_title(f'GC-stratified motif correlations — {spacing} bp spacing')
        ax.set_ylim(-1.05,1.05)
        fig.tight_layout()
        fig.savefig(p / f'exon_gc_{spacing}bp.png',dpi=160)
        plt.close(fig)


# Produce LaTeX tables of correlation extremes and GC-group comparisons.
# Write one insertable set of tables per sequence spacing.
def write_latex_tables(out, spacings, extreme, comparisons):
    """Write stand-alone insertable LaTeX tabular environments (no packages)."""
    # Format statistical values for LaTeX table cells.
    # Use a dash for missing values and concise numeric precision otherwise.
    def fmt(value):
        if value is None:
            return '--'
        if isinstance(value, (int, float)):
            return f'{value:.3g}'
        return str(value)

    lines = ['% Generated by Hercules-NF. Statistical methods are documented in README.md.']
    for spacing in spacings:
        lines.extend([
            r'\begin{table}[htbp]',
            r'\centering',
            r'\begin{tabular}{llllr}',
            r'\hline',
            r'Spacing & Tail & Rank & Motif & Pearson $r$ \\',
            r'\hline',
        ])
        for row in extreme:
            if row['spacing'] == spacing:
                lines.append(f"{spacing} bp & {row['tail']} & {row['rank']} & "
                             f"{row['motif']} & {fmt(row['pearson_r'])} " + r'\\')
        lines.extend([
            r'\hline',
            r'\end{tabular}',
            rf'\caption{{Extreme defined motif correlations at {spacing} bp spacing.}}',
            r'\end{table}',
            '',
            r'\begin{table}[htbp]',
            r'\centering',
            r'\begin{tabular}{llrrrr}',
            r'\hline',
            r'Feature GC & Motif GC & $n_{50}$ & $n_{x}$ & Welch FDR & MWU FDR \\',
            r'\hline',
        ])
        for row in comparisons:
            if row['spacing'] == spacing:
                lines.append(f"{row['gc_bin']}\\% & {fmt(row['comparison_motif_gc'])}\\% & "
                             f"{row['n_reference']} & {row['n_comparison']} & "
                             f"{fmt(row['welch_fdr'])} & {fmt(row['mannwhitney_fdr'])} " + r'\\')
        lines.extend([
            r'\hline',
            r'\end{tabular}',
            rf'\caption{{Comparisons with the 50\% motif-GC group at {spacing} bp; '
            r'Benjamini--Hochberg FDR-adjusted values. -- denotes insufficient data.}',
            r'\end{table}',
            '',
        ])
    (Path(out) / 'tables.tex').write_text('\n'.join(lines) + '\n')


# Parse command-line options and coordinate the Phase II analysis.
# Load counts and GC data, run statistics and write CSV, plot and LaTeX outputs.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts-dir',required=True)
    parser.add_argument('--gc',required=True)
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--k',type=int,default=4)
    parser.add_argument('--spacings',type=int,nargs='+',default=[10,50,100,200])
    parser.add_argument('--min-pairs',type=int,default=10)
    parser.add_argument('--outlier-count',type=int,default=10)
    args = parser.parse_args()
    if args.k < 1 or args.min_pairs < 2 or any(d<=0 for d in args.spacings):
        parser.error('--k >= 1, --min-pairs >= 2, and positive spacings required')
    out = Path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    gc = load_gc(args.gc)
    files = discover_files(args.counts_dir,args.k)
    # Historical settings: 10/50 bp have +/-2 bp; 100/200 have +/-4 bp.
    tolerance = {d:2 if d<=50 else 4 for d in args.spacings}
    correlations = []
    for motif,file in sorted(files.items()):
        correlations.extend(calculate_motif(motif,file,gc,args.spacings,tolerance,args.min_pairs))
    comparisons = group_comparisons(correlations,args.spacings)
    extreme = outliers(correlations,args.spacings,args.outlier_count)
    write_csv(out/'correlations.csv',correlations,
              ['motif','motif_gc','quartile','spacing','tolerance','gc_bin','pairs','pearson_r'])
    write_csv(out/'gc_comparisons.csv',comparisons,
              ['spacing','gc_bin','reference_motif_gc','comparison_motif_gc',
               'n_reference','n_comparison','welch_p','mannwhitney_p','welch_fdr','mannwhitney_fdr'])
    write_csv(out/'outliers.csv',extreme,['spacing','tail','rank','motif','pearson_r','pairs'])
    plots(correlations,args.spacings,out)
    write_latex_tables(out,args.spacings,extreme,comparisons)
    summary = {'motifs_analysed':len(files), 'spacings_bp':args.spacings,
               'tolerance_bp':tolerance,'min_pairs':args.min_pairs,
               'feature_gc_bands': [f'{lo}-{hi}' for lo,hi in GC_BINS],
               'correlation_rows':len(correlations),
               'defined_correlations':sum(r['pearson_r'] is not None for r in correlations),
               'test_family':'all valid tests per method across spacings and GC bins',
               'quartile_method':'numpy percentile(25,method=higher), per motif',
               'nonparametric_test':'Mann-Whitney U (independent samples; not legacy Wilcoxon)',
               'historical_equivalence':False}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,sort_keys=True))


if __name__=='__main__':
    main()
