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
import csv
import random
import tempfile
import unittest
from pathlib import Path
from optimized_demux import BarcodeIndex, EGFP_ANCHOR, run_baseline


def reference_match(query, barcodes):
    best, distance, tied = -1, 3, False
    for i, barcode in enumerate(barcodes):
        if query == barcode:
            return i
        d = sum(a != b for a, b in zip(query, barcode))
        if d < distance:
            best, distance, tied = i, d, False
        elif d == distance:
            tied = True
    return -1 if tied else best


class OptimizedReferenceTest(unittest.TestCase):
    def test_seed_index_matches_brute_force_oracle(self):
        rng = random.Random(42)
        refs = [bytes(rng.choices(b'ACGT', k=25)) for _ in range(150)]
        index = BarcodeIndex(refs)
        for ref in refs:
            for mutations in range(4):
                query = bytearray(ref)
                for pos in rng.sample(range(25), mutations):
                    query[pos] = rng.choice([b for b in b'ACGT' if b != query[pos]])
                query = bytes(query)
                self.assertEqual(reference_match(query, refs), index.match(query))

    def test_ambiguous_best_match_and_duplicate_seeds(self):
        index = BarcodeIndex([b'A'*25, b'C'+b'A'*24, b'G'*25])
        self.assertEqual(-1, index.match(b'T'+b'A'*24))
        self.assertEqual(0, index.match(b'A'*24+b'C'))
        self.assertEqual(0, BarcodeIndex([b'A'*25, b'A'*25]).match(b'A'*25))

    def test_direct_counts_preserve_cell_variant_and_umi_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            refs = [b'A'*25, b'C'*25]
            (p/'library.tsv').write_bytes(b'variant_name\tbarcode_25bp\nv0\t'+refs[0]+b'\nv1\t'+refs[1]+b'\n')
            cell = b'ACGT'*4
            molecules = [(cell+b'A'*12, refs[0]), (cell+b'A'*12, refs[0]),
                         (cell+b'C'*12, refs[0]), (cell+b'A'*12, refs[1])]
            for read in (1, 2):
                with (p/f'R{read}.fastq').open('wb') as out:
                    for i, (r1, barcode) in enumerate(molecules):
                        sequence = r1 if read == 1 else b'TTTGGCC'+EGFP_ANCHOR+barcode+b'A'*25
                        out.write(b'@'+str(i).encode()+b'\n'+sequence+b'\n+\n'+b'I'*len(sequence)+b'\n')
            run_baseline(p/'library.tsv', p/'R1.fastq', p/'R2.fastq', p/'out.tsv')
            with (p/'out.tsv').open() as handle:
                counts = {r['variant_name']: int(r['umi_count']) for r in csv.DictReader(handle, delimiter='\t')}
            self.assertEqual({'v0': 2, 'v1': 1}, counts)


if __name__ == '__main__':
    unittest.main()
