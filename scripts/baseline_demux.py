#!/usr/bin/env python3
# Copyright 2026 TechTune Media Group
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Readable Python reference authored for jseq-zero, with no profiling overhead.

Uses ordinary strings, sets and counters. Exact barcodes use a dictionary;
non-exact barcodes use a straightforward Hamming-distance comparison.
"""
import argparse
from collections import Counter
import csv
import gzip
from itertools import zip_longest
import time

ANCHOR = 'ACAGTCGCAATTCCGGTTAA'


def fastq_sequences(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        while handle.readline():
            sequence = handle.readline().rstrip('\r\n')
            separator = handle.readline()
            quality = handle.readline()
            if not separator or not quality:
                raise ValueError('Truncated FASTQ record')
            yield sequence


def find_barcode_offset(sequence):
    start = sequence.find(ANCHOR, 0, 40 + len(ANCHOR))
    if start != -1:
        return start + len(ANCHOR)
    for start in range(min(len(sequence) - 25, 40)):
        window = sequence[start:start + len(ANCHOR)]
        mismatches = sum(left != right for left, right in zip(window, ANCHOR))
        if mismatches <= 1:
            return start + len(ANCHOR)
    return -1


def match_barcode(barcode, library, exact_lookup):
    if barcode in exact_lookup:
        return exact_lookup[barcode]
    best_variant, best_distance, tied = None, 3, False
    for variant, reference in library:
        distance = sum(left != right for left, right in zip(barcode, reference))
        if distance < best_distance:
            best_variant, best_distance, tied = variant, distance, False
        elif distance == best_distance:
            tied = True
    return None if tied else best_variant


def run_reference(library_path, r1_path, r2_path, output_path):
    started = time.perf_counter()
    with open(library_path) as handle:
        library = [(row['variant_name'], row['barcode_25bp'])
                   for row in csv.DictReader(handle, delimiter='\t')]
    exact_lookup = {}
    for name, barcode in library:
        if len(barcode) != 25 or any(base not in 'ACGT' for base in barcode):
            raise ValueError('Reference barcodes require 25 uppercase A/C/G/T bases')
        exact_lookup.setdefault(barcode, name)
    seen_molecules = set()
    counts = Counter()
    total = no_anchor = unmatched = duplicates = 0
    for read1, read2 in zip_longest(fastq_sequences(r1_path), fastq_sequences(r2_path)):
        if read1 is None or read2 is None:
            raise ValueError('Paired FASTQ files have different record counts')
        total += 1
        cell_barcode, umi = read1[:16], read1[16:28]
        if len(read1) < 28 or 'N' in cell_barcode.upper() or 'N' in umi.upper():
            continue
        offset = find_barcode_offset(read2)
        if offset == -1:
            no_anchor += 1
            continue
        barcode = read2[offset:offset+25]
        if len(barcode) != 25 or 'N' in barcode.upper():
            unmatched += 1
            continue
        variant = match_barcode(barcode, library, exact_lookup)
        if variant is None:
            unmatched += 1
            continue
        molecule = (cell_barcode, umi, variant)
        if molecule in seen_molecules:
            duplicates += 1
        else:
            seen_molecules.add(molecule)
            counts[cell_barcode, variant] += 1
    with open(output_path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['cell_barcode', 'variant_name', 'umi_count'])
        for (cell, variant), count in sorted(counts.items()):
            writer.writerow([cell, variant, count])
    print(f'Total Pairs: {total:,} | Unique UMIs: {len(seen_molecules):,} | '
          f'PCR Dups: {duplicates:,} | No Flank: {no_anchor:,} | Unmatched: {unmatched:,}')
    print(f'Readable Python Elapsed: {time.perf_counter()-started:.2f} sec')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('library_tsv', 'r1_fastq', 'r2_fastq', 'out_tsv'):
        parser.add_argument(name)
    args = parser.parse_args()
    run_reference(args.library_tsv, args.r1_fastq, args.r2_fastq, args.out_tsv)
