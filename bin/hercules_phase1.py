#!/usr/bin/env python3
"""Hercules-NF Phase I: SAM partitioning, feature-indexed k-mer counting, reductions.

This is an optimised Python 3 implementation of the Chapter 6 measurement,

Positions are zero-based relative to the first genomic base of the selected feature. Counts are read observations of complete k-mers within one contiguous CIGAR match block.
"""
from __future__ import annotations
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import json
import os
from pathlib import Path
import re
import sqlite3
import sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
BASES = set('ACGT')


# Open plain-text or gzip-compressed files through the same interface.
# The compression format is selected from the filename suffix.
def open_text(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


# Partition SAM alignment records across a configurable number of shards.
# Skip SAM headers and distribute reads round-robin for parallel processing.
# Report the number of alignments assigned to the output files.
def split_sam(args):
    if args.partitions < 1:
        raise ValueError('--partitions must be >= 1')
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    files = [open(out / f'shard_{i:05d}.sam', 'w') for i in range(args.partitions)]
    n = 0
    try:
        with open_text(args.sam) as source:
            for line in source:
                if not line.strip() or line.startswith('@'):
                    continue
                files[n % args.partitions].write(line)
                if not line.endswith('\n'):
                    files[n % args.partitions].write('\n')
                n += 1
    finally:
        for f in files:
            f.close()
    print(json.dumps({'input_alignments': n, 'shards': args.partitions}))


# Index genomic annotations by chromosome for efficient feature assignment.
# Sorted intervals and cumulative endpoints support overlapping GTF features.
class FeatureIndex:
    """Chromosome-indexed GTF intervals with overlap-aware bounded reverse scan."""
    # Read the selected GTF feature type and discard duplicate intervals.
    # For each chromosome, build sorted starts and cumulative maximum ends
    # to support a bounded search for overlapping genomic features.
    def __init__(self, gtf, feature_type):
        grouped = defaultdict(set)
        with open_text(gtf) as source:
            for line in source:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 9:
                    raise ValueError('Malformed GTF line (requires 9 fields)')
                chrom, _, kind, s, e, _, strand, _, _ = fields[:9]
                if kind != feature_type:
                    continue
                start, end = int(s), int(e)
                if start < 1 or end < start:
                    raise ValueError(f'Bad GTF interval: {chrom}:{s}-{e}')
                grouped[chrom].add((start, end, strand))
        if not grouped:
            raise ValueError(f'No {feature_type!r} features in {gtf}')
        self.index = {}
        for chrom, intervals in grouped.items():
            ordered = sorted(intervals)
            starts = [x[0] for x in ordered]
            max_end = []
            max_seen = -1
            for _, end, _ in ordered:
                max_seen = max(max_seen, end)
                max_end.append(max_seen)
            self.index[chrom] = (ordered, starts, max_end)

    # Assign an alignment to the feature covering the most reference bases.
    # Search each matching CIGAR block for overlaps with indexed annotations;
    # resolve equal coverage deterministically, or return None if unmatched.
    def best(self, chrom, blocks):
        value = self.index.get(chrom)
        if value is None:
            return None
        intervals, starts, max_end = value
        covered = defaultdict(int)
        for ref_start, ref_end, _, _ in blocks:  # both inclusive
            i = bisect_right(starts, ref_end) - 1
            while i >= 0 and max_end[i] >= ref_start:
                start, end, strand = intervals[i]
                overlap = max(0, min(ref_end, end) - max(ref_start, start) + 1)
                if overlap:
                    covered[(start, end, strand)] += overlap
                i -= 1
        if not covered:
            return None
        # One feature per read, as in Hercules. Deterministic for ties.
        best = min(covered, key=lambda k: (-covered[k], k[0], k[1], k[2]))
        return (chrom, *best)


# Decode a SAM CIGAR string into contiguous aligned reference/query blocks.
# Merge adjacent matching operations but break at indels and skipped regions;
# validate that the CIGAR consumes the full query sequence.
def matching_blocks(pos, cigar, seq_length):
    """Returns (ref_start, ref_end, query_start, query_end) for contiguous matches.

    Reference coordinates are SAM/GTF 1-based inclusive; query indices are
    zero-based and half-open. Insertions/deletions/Ns break match blocks.
    """
    if cigar == '*':
        return []
    tokens = CIGAR_RE.findall(cigar)
    if not tokens or ''.join(n + op for n, op in tokens) != cigar:
        raise ValueError(f'Invalid CIGAR: {cigar}')
    ref, query = pos, 0
    result = []
    for length_str, op in tokens:
        length = int(length_str)
        if length < 1:
            raise ValueError('Zero-length CIGAR operation')
        if op in 'M=X':
            # Adjacent M/= /X operations remain one contiguous aligned run.
            # Only reference/query-consuming gaps actually break a k-mer.
            if result and result[-1][1] + 1 == ref and result[-1][3] == query:
                a, _, b, _ = result[-1]
                result[-1] = (a, ref + length - 1, b, query + length)
            else:
                result.append((ref, ref + length - 1, query, query + length))
            ref += length
            query += length
        elif op in 'IS':
            query += length
        elif op in 'DN':
            ref += length
        # H and P consume neither sequence nor reference
    if query != seq_length:
        raise ValueError(f'CIGAR consumes {query} query bases but SEQ length={seq_length}')
    return result


# Calculate the GC percentage of a sequencing read.
# Exclude ambiguous/non-ACGT characters from the denominator and return
# None if the sequence contains no canonical nucleotide bases.
def gc_percent(seq):
    valid = sum(c in BASES for c in seq)
    return 100.0 * sum(c in 'GC' for c in seq) / valid if valid else None


# Open an SQLite database for intermediate Phase I results.
# Initialise indexed motif-position counts, per-feature GC aggregates and
# alignment-processing metrics for later cross-shard reductions.
def database(path):
    conn = sqlite3.connect(path)
    conn.execute('PRAGMA journal_mode=DELETE')
    conn.execute('PRAGMA synchronous=NORMAL')
    conn.execute('PRAGMA temp_store=MEMORY')
    conn.executescript('''
        CREATE TABLE IF NOT EXISTS counts (
            feature TEXT NOT NULL, motif TEXT NOT NULL,
            pos INTEGER NOT NULL, n INTEGER NOT NULL,
            PRIMARY KEY (motif, feature, pos)
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS gc (
            feature TEXT PRIMARY KEY, gc_sum REAL NOT NULL, n INTEGER NOT NULL
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS metrics (key TEXT PRIMARY KEY, value INTEGER);
    ''')
    return conn


# Persist the in-memory motif and GC accumulators to SQLite in one transaction.
# Existing feature/position keys are incremented using SQL upserts, allowing
# bounded memory usage while scanning large read partitions.
def flush(conn, counts, gcs):
    with conn:
        conn.executemany('''
            INSERT INTO counts(feature,motif,pos,n) VALUES (?,?,?,?)
            ON CONFLICT(motif,feature,pos) DO UPDATE SET n=n+excluded.n
        ''', ((f,m,p,n) for (f,m,p),n in counts.items()))
        conn.executemany('''
            INSERT INTO gc(feature,gc_sum,n) VALUES (?,?,?)
            ON CONFLICT(feature) DO UPDATE SET
                gc_sum=gc_sum+excluded.gc_sum,n=n+excluded.n
        ''', ((f,v[0],v[1]) for f,v in gcs.items()))
    counts.clear()
    gcs.clear()


# Process one SAM shard, assigning aligned reads to annotated GTF features.
# Apply alignment filters, measure per-read GC and count every eligible k-mer
# at its zero-based feature-relative position; flush aggregates to SQLite.
def map_sam(args):
    index = FeatureIndex(args.gtf, args.feature)
    conn = database(args.output)
    counts = Counter()
    gcs = defaultdict(lambda: [0.0, 0])
    metrics = Counter()
    try:
        with open_text(args.sam) as source:
            for line in source:
                if not line.strip() or line.startswith('@'):
                    continue
                metrics['reads_seen'] += 1
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 11:
                    metrics['malformed'] += 1
                    continue
                try:
                    flag, rpos, mapq = int(parts[1]), int(parts[3]), int(parts[4])
                    if flag & 0x4 or flag & 0x100 or flag & 0x800 or flag & 0x200:
                        metrics['flag_filtered'] += 1
                        continue
                    if args.exclude_duplicates and flag & 0x400:
                        metrics['duplicate_filtered'] += 1
                        continue
                    if mapq < args.min_mapq:
                        metrics['mapq_filtered'] += 1
                        continue
                    chrom, seq = parts[2], parts[9].upper()
                    if seq == '*' or rpos < 1:
                        metrics['malformed'] += 1
                        continue
                    blocks = matching_blocks(rpos, parts[5], len(seq))
                    if not blocks:
                        continue
                except (ValueError, IndexError):
                    metrics['malformed'] += 1
                    continue
                feature = index.best(chrom, blocks)
                if feature is None:
                    metrics['unassigned'] += 1
                    continue
                chrom, start, end, strand = feature
                key = f'{chrom}:{start}-{end}:{strand}'
                metrics['assigned'] += 1
                gc = gc_percent(seq)
                if gc is not None:
                    gcs[key][0] += gc
                    gcs[key][1] += 1
                for ref_start, ref_end, q_start, q_end in blocks:
                    # Only complete motif windows wholly contained in feature.
                    low = max(0, start - ref_start)
                    high = min(q_end-q_start, end-ref_start+1)
                    for offset in range(low, high-args.k+1):
                        read_pos = q_start + offset
                        if read_pos < args.primer_exclude:
                            continue
                        motif = seq[read_pos:read_pos+args.k]
                        if len(motif) != args.k or any(b not in BASES for b in motif):
                            continue
                        motif_pos = ref_start + offset - start
                        counts[(key, motif, motif_pos)] += 1
                        metrics['motif_observations'] += 1
                if len(counts) >= args.flush_entries:
                    flush(conn, counts, gcs)
        flush(conn, counts, gcs)
        with conn:
            conn.executemany('INSERT OR REPLACE INTO metrics VALUES (?,?)', metrics.items())
    finally:
        conn.close()
    print(json.dumps(dict(metrics), sort_keys=True))


# Merge per-feature GC sums and read counts from all shard databases.
# Write the weighted mean of read-level GC percentages for each feature,
# alongside total mapping and filtering metrics in JSON format.
def gc_reduce(args):
    totals = defaultdict(lambda: [0.0, 0])
    metrics = Counter()
    for db in args.dbs:
        conn = sqlite3.connect(f'file:{Path(db).resolve()}?mode=ro', uri=True)
        try:
            for key, s, n in conn.execute('SELECT feature,gc_sum,n FROM gc'):
                totals[key][0] += s
                totals[key][1] += n
            metrics.update(dict(conn.execute('SELECT key,value FROM metrics')))
        finally:
            conn.close()
    with open(args.output, 'w', newline='') as dest:
        writer = csv.writer(dest)
        for key in sorted(totals):
            s, n = totals[key]
            writer.writerow((key, f'{s/n:.12g}'))
    with open(args.metrics, 'w') as dest:
        json.dump(dict(metrics), dest, indent=2, sort_keys=True)


# Combine the counts for one requested motif across all shard databases.
# Merge pre-sorted feature/position records without loading them all into
# memory, and write the final feature,position,count CSV file.
def motif_reduce(args):
    if not args.motif or any(b not in BASES for b in args.motif):
        raise ValueError('Motif must contain only A,C,G,T')
    # Each database is pre-indexed on (motif, feature, position).
    # Stream sorted data and heap-merge to avoid whole-genome memory loading.
    import heapq
    conns = []
    # Stream feature-position counts for the selected motif from one shard.
    # SQLite supplies records in the order required by the merge operation.
    def records(conn):
        yield from conn.execute('SELECT feature,pos,n FROM counts WHERE motif=? ORDER BY feature,pos', (args.motif,))
    try:
        streams = []
        for db in args.dbs:
            conn = sqlite3.connect(f'file:{Path(db).resolve()}?mode=ro', uri=True)
            conns.append(conn)
            streams.append(records(conn))
        merged = heapq.merge(*streams)
        with open(args.output, 'w', newline='') as dest:
            writer = csv.writer(dest)
            previous = None
            total = 0
            for feature, pos, count in merged:
                key = (feature, pos)
                if key != previous:
                    if previous is not None:
                        writer.writerow((*previous, total))
                    previous, total = key, count
                else:
                    total += count
            if previous is not None:
                writer.writerow((*previous, total))
    finally:
        for conn in conns:
            conn.close()


# Define the Phase I command-line interface and validate its parameters.
# Dispatch split, map, GC reduction or motif reduction to the corresponding
# routine, enabling Nextflow to invoke each operation independently.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    split = sub.add_parser('split')
    split.add_argument('--sam', required=True)
    split.add_argument('--partitions', type=int, default=16)
    split.add_argument('--output-dir', required=True)
    split.set_defaults(func=split_sam)
    map_p = sub.add_parser('map')
    map_p.add_argument('--sam', required=True)
    map_p.add_argument('--gtf', required=True)
    map_p.add_argument('--output', required=True)
    map_p.add_argument('--k', type=int, default=4)
    map_p.add_argument('--feature', choices=['exon', 'CDS'], default='CDS')
    map_p.add_argument('--primer-exclude', type=int, default=10)
    map_p.add_argument('--min-mapq', type=int, default=0)
    map_p.add_argument('--exclude-duplicates', action='store_true')
    map_p.add_argument('--flush-entries', type=int, default=150000)
    map_p.set_defaults(func=map_sam)
    gc_p = sub.add_parser('gc-reduce')
    gc_p.add_argument('--dbs', nargs='+', required=True)
    gc_p.add_argument('--output', required=True)
    gc_p.add_argument('--metrics', required=True)
    gc_p.set_defaults(func=gc_reduce)
    motif_p = sub.add_parser('motif-reduce')
    motif_p.add_argument('--dbs', nargs='+', required=True)
    motif_p.add_argument('--motif', required=True)
    motif_p.add_argument('--output', required=True)
    motif_p.set_defaults(func=motif_reduce)
    args = parser.parse_args()
    if getattr(args, 'k', 1) < 1 or getattr(args, 'primer_exclude', 0) < 0:
        parser.error('--k must be >= 1 and --primer-exclude >= 0')
    args.func(args)


if __name__ == '__main__':
    main()
