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
"""Check the delivered native CLI, input failures and Python parity on a small fixture."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def counts(path):
    with path.open() as handle:
        return {(r['cell_barcode'],r['variant_name']):int(r['umi_count']) for r in csv.DictReader(handle,delimiter='\t')}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary',type=Path)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    binary=args.binary.resolve()
    project=Path(__file__).resolve().parents[1]
    checks=[]
    with tempfile.TemporaryDirectory(prefix='native-check-',dir=args.report.parent) as directory:
        root=Path(directory)
        def run(arguments,expect_success=True):
            result=subprocess.run([str(binary),*arguments],cwd=root,text=True,capture_output=True)
            if expect_success and result.returncode:
                raise RuntimeError(result.stderr)
            if not expect_success and result.returncode==0:
                raise AssertionError('Invalid CLI/input was accepted')
            return result
        env=os.environ.copy()
        env.pop('JAVA_HOME',None)
        env['PATH']='/nonexistent'
        help_result=subprocess.run([str(binary),'--help'],cwd=root,env=env,text=True,capture_output=True,check=True)
        assert 'amplicon-demux' in help_result.stdout
        checks.append('Native help works without Java or commands on PATH')
        run(['unknown-command'],False)
        checks.append('Unknown subcommand exits nonzero')
        run(['amplicon-demux'],False)
        checks.append('Missing demultiplexer arguments exit nonzero')
        run(['generate-benchmark','500','2','50'])
        checks.append('Native generator creates 1,000 paired reads')
        parameters=['library_barcodes.tsv','R1.fastq','R2.fastq']
        run(['-Xmx512m','amplicon-demux',*parameters,'native-counts.tsv'])
        expected=counts(root/'native-counts.tsv')
        for name in ('R1.fastq','R2.fastq'):
            with gzip.open(root/(name+'.gz'),'wb') as handle:handle.write((root/name).read_bytes())
        run(['-Xmx512m','amplicon-demux','library_barcodes.tsv','R1.fastq.gz','R2.fastq.gz','gzip-counts.tsv'])
        assert expected==counts(root/'gzip-counts.tsv')
        checks.append('Plain and gzip native count maps match')
        # An exact-match fixture has a known two-UMI result in both programs.
        cell='ACGT'*4
        (root/'exact-library.tsv').write_text('variant_name\tbarcode_25bp\nv0\t'+'A'*25+'\n')
        for read in (1,2):
            with (root/f'exact-R{read}.fastq').open('w') as handle:
                for i,umi in enumerate(('A'*12,'A'*12,'C'*12)):
                    sequence=cell+umi if read==1 else 'ACAGTCGCAATTCCGGTTAA'+'A'*25
                    handle.write(f'@{i}\n{sequence}\n+\n'+ 'I'*len(sequence)+'\n')
        exact_parameters=['exact-library.tsv','exact-R1.fastq','exact-R2.fastq']
        run(['amplicon-demux',*exact_parameters,'exact-native-counts.tsv'])
        subprocess.run([sys.executable,str(project/'scripts/baseline_demux.py'),*exact_parameters,'exact-python-counts.tsv'],cwd=root,check=True,capture_output=True,text=True)
        assert counts(root/'exact-native-counts.tsv')==counts(root/'exact-python-counts.tsv')=={(cell,'v0'):2}
        checks.append('Readable Python and native counts agree on a known exact-match fixture')
        (root/'truncated.fastq').write_text('@read\nACGT\n+\n')
        run(['amplicon-demux','library_barcodes.tsv','truncated.fastq','R2.fastq','bad-counts.tsv'],False)
        checks.append('Truncated FASTQ is rejected')
        (root/'empty.fastq').write_text('')
        run(['amplicon-demux','library_barcodes.tsv','R1.fastq','empty.fastq','unequal-counts.tsv'],False)
        checks.append('Unequal pair counts are rejected')
        run(['generate-benchmark','0','2','50'],False)
        checks.append('Invalid generator sizes are rejected')
    args.report.write_text(json.dumps({'passed':len(checks),'checks':checks},indent=2)+'\n')
    print(f'{len(checks)} native integration checks passed.')


if __name__=='__main__':
    main()
