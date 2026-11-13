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
"""One-off, reproducible Linux benchmark: readable Python and native GraalVM 25."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time

NAMES = ('readable-python', 'graalvm-native')


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def count_map(path):
    with path.open() as handle:
        return {(r['cell_barcode'],r['variant_name']):int(r['umi_count'])
                for r in csv.DictReader(handle,delimiter='\t')}


def capture(command):
    result = subprocess.run(command,capture_output=True,text=True)
    return (result.stdout or result.stderr).strip()


def hardware():
    cpu = {r['field'].rstrip(':'):r['data'] for r in json.loads(capture(['lscpu','-J']))['lscpu']}
    os_info = dict(line.split('=',1) for line in Path('/etc/os-release').read_text().splitlines() if '=' in line)
    memory_kib = int(re.search(r'MemTotal:\s+(\d+)',Path('/proc/meminfo').read_text())[1])
    affinity = sorted(os.sched_getaffinity(0))
    return {'os':os_info.get('PRETTY_NAME','').strip('"'), 'kernel':platform.release(),
            'architecture':platform.machine(), 'cpu_model':cpu.get('Model name'),
            'visible_cpus':os.cpu_count(), 'available_cpu_affinity':affinity,
            'benchmark_cpu_affinity':[affinity[0]], 'memory_gib':memory_kib/1024**2,
            'virtualization_detected':capture(['systemd-detect-virt']),
            'topology':{k:cpu.get(k) for k in ('Socket(s)','Core(s) per socket','Thread(s) per core','CPU max MHz','CPU min MHz','L1d cache','L1i cache','L2 cache','L3 cache')},
            'python':platform.python_version(), 'python_implementation':platform.python_implementation(),
            'load_average_before':os.getloadavg()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--build-info',type=Path,required=True)
    parser.add_argument('--read-pairs',type=int,required=True)
    parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args()
    project,data,out,binary=(p.resolve() for p in (args.project,args.data,args.out,args.binary))
    if args.repeats<1 or args.read_pairs<1:
        parser.error('positive repeats and read-pairs required')
    with binary.open('rb') as handle:
        binary_magic = handle.read(4)
    if binary_magic != b'\x7fELF':
        parser.error('this Linux comparison requires the native ELF executable')
    build=json.loads(args.build_info.read_text())
    if not re.search(r'native-image 25[. ]',build['graalvm_version']):
        parser.error('build metadata must identify GraalVM Native Image 25')
    if build['binary_sha256'] != sha256(binary):
        parser.error('native executable does not match build metadata')
    if (out/'results.json').exists():
        parser.error('choose a fresh output directory')
    out.mkdir(parents=True,exist_ok=True)
    files={name:{'sha256':sha256(data/name),'bytes':(data/name).stat().st_size}
           for name in ('library_barcodes.tsv','R1.fastq','R2.fastq')}
    with (data/'library_barcodes.tsv').open() as handle:
        variants=sum(1 for line in handle if line.strip())-1
    environment=hardware()
    source_files=['pom.xml','scripts/baseline_demux.py','scripts/benchmark_native.py']
    source_files += [str(p.relative_to(project)) for p in sorted((project/'src/main/java').rglob('*.java'))]
    report={'schema_version':2,'measured_at_utc':datetime.now(timezone.utc).isoformat(),
            'environment':environment,'build':build,
            'dataset':{'read_pairs':args.read_pairs,'library_variants':variants,'files':files},
            'source_sha256':{name:sha256(project/name) for name in source_files},
            'method':{'repeats':args.repeats,'implementations':list(NAMES),
                      'timing':'Full fresh-process wall time, including startup and input/output',
                      'memory':'Linux wait4 ru_maxrss in KiB, converted to MiB',
                      'native_heap_limit':'512m','allocation_profiling':False,
                      'run_order':'Rotated across repetitions',
                      'filesystem_cache':'Not flushed; all inputs read for hashing before measurements',
                      'cpu_affinity':'Only benchmark child processes pinned to one already-allowed CPU; no system settings changed',
                      'scope':'Both implementations authored in this repository; not a language-wide or biotech-lab comparison'},
            'runs':[]}
    outputs={}
    # Compare large count tables only after all timed children have exited.
    # A forked child can inherit the parent's resident-memory high-water mark.
    for repetition in range(args.repeats):
        order=NAMES[repetition%len(NAMES):]+NAMES[:repetition%len(NAMES)]
        for name in order:
            stem=f'{name}-{repetition+1}'
            output=out/(stem+'-counts.tsv')
            parameters=[str(data/'library_barcodes.tsv'),str(data/'R1.fastq'),str(data/'R2.fastq'),str(output)]
            if name=='graalvm-native':
                command=[str(binary),'-Xmx512m','amplicon-demux',*parameters]
            else:
                script='baseline_demux.py'
                command=[sys.executable,str(project/'scripts'/script),*parameters]
            print(f'Starting {name}, repetition {repetition+1}/{args.repeats}',flush=True)
            started=time.perf_counter()
            with (out/(stem+'.log')).open('w') as log:
                child=subprocess.Popen(command,cwd=project,stdout=log,stderr=subprocess.STDOUT,
                                       preexec_fn=lambda:os.sched_setaffinity(0,set(environment['benchmark_cpu_affinity'])))
                _,status,usage=os.wait4(child.pid,0)
                child.returncode=os.waitstatus_to_exitcode(status)
            elapsed=time.perf_counter()-started
            if child.returncode:
                raise RuntimeError(f'{stem} failed; inspect its log')
            text=(out/(stem+'.log')).read_text()
            telemetry={}
            for label,field in [('Total Pairs','read_pairs'),('Unique UMIs','unique_umis'),('PCR Dups','pcr_duplicates'),('No Flank','no_flank'),('Unmatched(?: Barcode)?','unmatched_barcode')]:
                match=re.search(label+r': ([\d,]+)',text)
                if not match:
                    raise RuntimeError(f'Missing {field} telemetry')
                telemetry[field]=int(match[1].replace(',',''))
            assert telemetry['read_pairs']==args.read_pairs
            output_hash=sha256(output)
            if name in outputs:
                assert output_hash==outputs[name]['sha256'], 'Output changed between repeats'
            outputs[name]={'path':output,'sha256':output_hash,'unique_umis':telemetry['unique_umis']}
            run={'implementation':name,'repeat':repetition+1,'wall_seconds':elapsed,
                 'peak_rss_mib':usage.ru_maxrss/1024,'user_cpu_seconds':usage.ru_utime,
                 'system_cpu_seconds':usage.ru_stime,'telemetry':telemetry,'output_sha256':output_hash}
            report['runs'].append(run)
            (out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
            print(f'{stem}: {elapsed:.3f}s, {run["peak_rss_mib"]:.1f} MiB, {telemetry["unique_umis"]:,} UMIs',flush=True)
    for name, output in outputs.items():
        observed=count_map(output['path'])
        assert sum(observed.values())==output['unique_umis']
        outputs[name]=observed
    report['summary']={}
    for name in NAMES:
        runs=[r for r in report['runs'] if r['implementation']==name]
        report['summary'][name]={'median_wall_seconds':statistics.median(r['wall_seconds'] for r in runs),
                                'min_wall_seconds':min(r['wall_seconds'] for r in runs),
                                'max_wall_seconds':max(r['wall_seconds'] for r in runs),
                                'median_peak_rss_mib':statistics.median(r['peak_rss_mib'] for r in runs),
                                'telemetry':runs[0]['telemetry']}
    python,native=outputs['readable-python'],outputs['graalvm-native']
    keys=python.keys()|native.keys()
    report['output_comparison']={'native_identical_to_python':python==native,
         'cell_variant_rows_differing':sum(python.get(k,0)!=native.get(k,0) for k in keys),
         'extra_native_counts':sum(max(0,native.get(k,0)-python.get(k,0)) for k in keys),
         'extra_python_counts':sum(max(0,python.get(k,0)-native.get(k,0)) for k in keys),
         'interpretation':'Different anchor matching; count differences are not a ground-truth accuracy measurement'}
    report['environment']['load_average_after']=os.getloadavg()
    (out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    with (out/'results.csv').open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['implementation','repeat','wall_seconds','peak_rss_mib','read_pairs','unique_umis','pcr_duplicates','no_flank','unmatched_barcode'])
        for r in report['runs']:
            t=r['telemetry']
            writer.writerow([r['implementation'],r['repeat'],r['wall_seconds'],r['peak_rss_mib'],t['read_pairs'],t['unique_umis'],t['pcr_duplicates'],t['no_flank'],t['unmatched_barcode']])
    print('All native and Python runs completed; output totals reconciled.',flush=True)


if __name__=='__main__':
    main()
