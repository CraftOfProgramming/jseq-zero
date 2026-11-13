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
"""Create shareable PNG, SVG and PDF comparison charts from measured results."""
import argparse
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    data = json.loads(args.results.read_text())
    if 'summary' not in data or len(data['runs']) != 2*data['method']['repeats']:
        parser.error('paired benchmark must be complete before plotting')
    output = args.out or args.results.parent
    output.mkdir(parents=True, exist_ok=True)
    baseline = data['summary']['readable-python']
    java = data['summary']['graalvm-native']
    speedup = baseline['median_wall_seconds']/java['median_wall_seconds']
    memory_reduction = 100*(1-java['median_peak_rss_mib']/baseline['median_peak_rss_mib'])
    ink, muted, accent, python_color, background = '#183044', '#526174', '#087F8C', '#7C8BA1', '#F5F8FA'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.spines.left': False, 'axes.edgecolor': '#CCD5DD',
                         'text.color': ink, 'axes.labelcolor': muted,
                         'xtick.color': muted, 'ytick.color': ink})
    fig = plt.figure(figsize=(16, 9), facecolor=background)
    fig.text(.06, .947, 'JSEQ ZERO  /  MEASURED PERFORMANCE', color=accent, weight='bold', size=12)
    pair_count = data['dataset']['read_pairs']
    title = f'{pair_count/1_000_000:g} million FASTQ read pairs' if pair_count >= 1_000_000 else f'{pair_count:,} FASTQ read pairs'
    fig.text(.06, .88, title, weight='bold', size=31)
    fig.text(.06, .837, 'Readable Python reference vs. jseq-zero GraalVM 25 native  |  Synthetic 50-variant library', size=14, color=muted)
    fig.add_artist(Line2D([.06,.94],[.793,.793],transform=fig.transFigure,color='#D9E2E8',linewidth=1))
    fig.text(.06, .742, 'ELAPSED TIME', weight='bold', size=12, color=muted)
    fig.text(.56, .742, 'PEAK PROCESS MEMORY', weight='bold', size=12, color=muted)
    fig.text(.06, .684, f'{speedup:.1f}× faster', weight='bold', size=29, color=accent)
    memory_title = f'{memory_reduction:.0f}% lower peak RSS' if memory_reduction >= 0 else f'{-memory_reduction:.0f}% higher peak RSS'
    fig.text(.56, .684, memory_title, weight='bold', size=29, color=accent if memory_reduction >= 0 else ink)
    repeats = data['method']['repeats']
    fig.text(.06,.64,f'Native Java speedup over our readable Python reference',color=muted,size=11)
    fig.text(.56,.64,'Same OS measurement for both processes',color=muted,size=11)

    for rectangle, field, xlabel in [([.20,.365,.29,.23], 'wall_seconds', 'Seconds · lower is better'),
                                     ([.70,.365,.23,.23], 'peak_rss_mib', 'MiB · lower is better')]:
        ax = fig.add_axes(rectangle, facecolor=background)
        values = [[r[field] for r in data['runs'] if r['implementation'] == name]
                  for name in ('readable-python', 'graalvm-native')]
        medians = [statistics.median(v) for v in values]
        upper = max(max(v) for v in values)
        for y, median, samples, color in zip((1,0),medians,values,(python_color,accent)):
            ax.barh(y,median,height=.38,color=color,zorder=3)
            ax.errorbar(median,y,xerr=[[median-min(samples)],[max(samples)-median]],fmt='none',ecolor=ink,capsize=5,lw=1.2,zorder=4)
            label = f'{median:.2f} s' if field == 'wall_seconds' else f'{median:.1f} MiB'
            ax.text(max(samples)+upper*.035,y,label,va='center',weight='bold',size=12)
        ax.set_yticks([1,0],['Readable Python\nreference','jseq-zero\nGraalVM 25 native'])
        ax.tick_params(axis='y',length=0,pad=12)
        ax.set_xlim(0,upper*1.42)
        ax.set_ylim(-.65,1.65)
        ax.set_xlabel(xlabel,labelpad=10,size=11)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.grid(axis='x',color='#DDE5EB',lw=.7,zorder=0)
        ax.set_axisbelow(True)

    fig.add_artist(Line2D([.06,.94],[.275,.275],transform=fig.transFigure,color='#D9E2E8',linewidth=1))
    fig.text(.06,.235,f'Counted UMIs: Python {baseline["telemetry"]["unique_umis"]:,}  |  Java {java["telemetry"]["unique_umis"]:,}. Outputs differ because anchor matching differs.',size=12,weight='bold')
    fig.text(.06,.19,'Both implementations were authored for this project. This is not a benchmark of biotech-lab software or clinical accuracy.',size=11,color=muted)
    env = data['environment']
    fig.text(.06,.149,f'Intel Core i7-12700  ·  RHEL 9.6  ·  CPython {env["python"]}  ·  Oracle GraalVM 25.0.3  ·  one CPU  ·  {data["measured_at_utc"][:10]}',size=11,color=muted)
    fig.text(.06,.108,f'Medians of {repeats} fresh processes; whiskers show min–max. Timing includes startup, imports and I/O. Allocation tracing is disabled.',size=10.5,color=muted)
    fig.text(.06,.071,'Filesystem cache was not flushed. Peak RSS includes the runtime; not Java heap use. 30.8 GiB system RAM. Raw results and methods accompany this chart.',size=10.5,color=muted)
    for extension in ('png','svg','pdf'):
        path = output/('benchmark-comparison.'+extension)
        fig.savefig(path,dpi=200,facecolor=background,metadata={'Creator':'jseq-zero benchmark plotting script'})
        if extension == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        print(path)
    plt.close(fig)


if __name__ == '__main__':
    main()
