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
import tempfile
import unittest
from pathlib import Path
from baseline_demux import ANCHOR, match_barcode, run_reference


class BaselineTest(unittest.TestCase):
    def test_exact_and_one_or_two_mismatches(self):
        library=[('v0','A'*25),('v1','C'*25)]
        exact={barcode:name for name,barcode in library}
        for query in ('A'*25,'G'+'A'*24,'GG'+'A'*23):
            self.assertEqual('v0',match_barcode(query,library,exact))
        self.assertIsNone(match_barcode('GGG'+'A'*22,library,exact))

    def test_ambiguous_best_match_and_duplicate_barcodes(self):
        library=[('v0','A'*25),('v1','C'+'A'*24)]
        exact={barcode:name for name,barcode in library}
        self.assertIsNone(match_barcode('T'+'A'*24,library,exact))
        self.assertEqual('v0',match_barcode('A'*25,library+[('v2','A'*25)],exact))

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
                        sequence = r1 if read == 1 else b'TTTGGCC'+ANCHOR.encode()+barcode+b'A'*25
                        out.write(b'@'+str(i).encode()+b'\n'+sequence+b'\n+\n'+b'I'*len(sequence)+b'\n')
            run_reference(p/'library.tsv', p/'R1.fastq', p/'R2.fastq', p/'out.tsv')
            with (p/'out.tsv').open() as handle:
                counts = {r['variant_name']: int(r['umi_count']) for r in csv.DictReader(handle, delimiter='\t')}
            self.assertEqual({'v0': 2, 'v1': 1}, counts)


if __name__ == '__main__':
    unittest.main()
