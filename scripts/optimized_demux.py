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
"""Allocation-conscious Python reference authored for this project.

Binary input, exact lookup, a three-seed mismatch index, a bounded cache,
integer molecule identities and direct counters reduce avoidable work.
Python still allocates. Anchor matching retains the original reference rules.
"""
import argparse
from functools import lru_cache
import gzip
import time

EGFP_ANCHOR = b'ACAGTCGCAATTCCGGTTAA'
ANCHOR_LEN = len(EGFP_ANCHOR)
BASE_CODES = bytes.maketrans(b'ACTG', bytes((0, 1, 2, 3)))
BASE_LOW_BITS = int.from_bytes(b'\x01' * 25, 'big')
VARIANT_BITS = 20
VARIANT_MASK = (1 << VARIANT_BITS) - 1


def encode_barcode(barcode):
    # One byte per base lets C-backed operations encode without a Python loop.
    return int.from_bytes(barcode.translate(BASE_CODES), 'big')


class BarcodeIndex:
    def __init__(self, barcodes):
        self.encoded = tuple(encode_barcode(b) for b in barcodes)
        self.exact = {}
        self.seeds = ({}, {}, {})
        for index, barcode in enumerate(barcodes):
            if len(barcode) != 25 or any(b not in b'ACGT' for b in barcode):
                raise ValueError('Reference barcodes require 25 uppercase A/C/G/T bases')
            self.exact.setdefault(barcode, index)
            for table, seed in zip(self.seeds, (barcode[:9], barcode[9:17], barcode[17:])):
                table.setdefault(seed, []).append(index)
        self.approximate = lru_cache(maxsize=8192)(self._match_approximate)

    def match(self, barcode):
        exact = self.exact.get(barcode)
        return exact if exact is not None else self.approximate(barcode)

    def _match_approximate(self, barcode):
        query = encode_barcode(barcode)
        best, best_distance, tie = -1, 3, False
        # Two substitutions leave at least one of three disjoint seeds intact.
        buckets = (self.seeds[0].get(barcode[:9], ()),
                   self.seeds[1].get(barcode[9:17], ()),
                   self.seeds[2].get(barcode[17:], ()))
        for bucket in buckets:
            for candidate in bucket:
                if candidate == best:
                    continue
                xor = query ^ self.encoded[candidate]
                distance = ((xor | (xor >> 1)) & BASE_LOW_BITS).bit_count()
                if distance < best_distance:
                    best, best_distance, tie = candidate, distance, False
                elif distance == best_distance:
                    tie = True
        return -1 if tie else best


def find_barcode_offset(read2, max_search=40):
    exact = read2.find(EGFP_ANCHOR, 0, max_search + ANCHOR_LEN)
    if exact != -1:
        return exact + ANCHOR_LEN
    for start in range(min(len(read2)-25, max_search)):
        mismatches = 0
        for offset, base in enumerate(EGFP_ANCHOR):
            if read2[start+offset] != base:
                mismatches += 1
                if mismatches > 1:
                    break
        if mismatches <= 1:
            return start + ANCHOR_LEN
    return -1


def open_fastq(path):
    return gzip.open(path, 'rb') if str(path).endswith('.gz') else open(path, 'rb', buffering=1024*1024)


def run_baseline(lib_tsv, r1_path, r2_path, out_tsv):
    started = time.perf_counter()
    names, barcodes = [], []
    with open(lib_tsv, 'rb') as library:
        next(library)
        for line in library:
            if line.strip():
                name, barcode = line.rstrip(b'\r\n').split(b'\t', 1)
                names.append(name.decode('utf-8'))
                barcodes.append(barcode)
    if len(names) > VARIANT_MASK:
        raise ValueError('At most 1,048,575 variants are supported')
    index = BarcodeIndex(barcodes)
    seen_molecules = set()
    cell_variant_counts = {}
    total_pairs = flank_missed = unmatched = pcr_dups = 0
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2:
        while True:
            header1, header2 = f1.readline(), f2.readline()
            if not header1 or not header2:
                if bool(header1) != bool(header2):
                    raise ValueError('Paired FASTQ files have different record counts')
                break
            seq1, seq2 = f1.readline(), f2.readline()
            plus1, plus2 = f1.readline(), f2.readline()
            quality1, quality2 = f1.readline(), f2.readline()
            if not (plus1 and plus2 and quality1 and quality2):
                raise ValueError('Truncated FASTQ record')
            total_pairs += 1
            identity = seq1[:28]
            if len(identity) < 28 or b'N' in identity or b'n' in identity or b'\n' in identity or b'\r' in identity:
                continue
            barcode_start = find_barcode_offset(seq2)
            if barcode_start == -1:
                flank_missed += 1
                continue
            barcode = seq2[barcode_start:barcode_start+25]
            if len(barcode) != 25 or b'N' in barcode or b'n' in barcode or b'\n' in barcode or b'\r' in barcode:
                unmatched += 1
                continue
            variant = index.match(barcode)
            if variant < 0:
                unmatched += 1
                continue
            cb_umi = int.from_bytes(identity, 'big')
            molecule = (cb_umi << VARIANT_BITS) | variant
            if molecule in seen_molecules:
                pcr_dups += 1
                continue
            seen_molecules.add(molecule)
            cell_variant = ((cb_umi >> 96) << VARIANT_BITS) | variant
            cell_variant_counts[cell_variant] = cell_variant_counts.get(cell_variant, 0) + 1
    with open(out_tsv, 'w', buffering=1024*1024) as out:
        out.write('cell_barcode\tvariant_name\tumi_count\n')
        for key in sorted(cell_variant_counts):
            cell = (key >> VARIANT_BITS).to_bytes(16, 'big').decode('ascii')
            variant = names[key & VARIANT_MASK]
            out.write(f'{cell}\t{variant}\t{cell_variant_counts[key]}\n')
    print(f'Total Pairs: {total_pairs:,} | Unique UMIs: {len(seen_molecules):,} | '
          f'PCR Dups: {pcr_dups:,} | No Flank: {flank_missed:,} | Unmatched: {unmatched:,}')
    print(f'Python Elapsed: {time.perf_counter()-started:.2f} sec')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('library_tsv', 'r1_fastq', 'r2_fastq', 'out_tsv'):
        parser.add_argument(argument)
    args = parser.parse_args()
    run_baseline(args.library_tsv, args.r1_fastq, args.r2_fastq, args.out_tsv)
