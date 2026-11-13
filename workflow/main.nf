#!/usr/bin/env nextflow
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
nextflow.enable.dsl = 2

params.reads       = "${projectDir}/../data/*_R{1,2}.fastq.gz"
params.library_tsv = "${projectDir}/../data/library_barcodes.tsv"
params.jar_path    = "${projectDir}/../target/jseq-zero.jar"
params.scanpy_h5ad = ""
params.outdir      = "${projectDir}/../results"

process AMPLICON_FAST_DEMUX {
    tag "${sample_id}"
    publishDir "${params.outdir}/umi_counts", mode: 'copy'

    cpus 2
    memory '1 GB'

    input:
    tuple val(sample_id), path(reads)
    path library_tsv
    path jar_file

    output:
    tuple val(sample_id), path("${sample_id}_umi_counts.tsv"), emit: counts
    path "${sample_id}_telemetry.log", emit: logs

    script:
    def (r1, r2) = reads
    """
    java -Xmx768m -jar ${jar_file} amplicon-demux \
        ${library_tsv} \
        ${r1} \
        ${r2} \
        ${sample_id}_umi_counts.tsv > ${sample_id}_telemetry.log
    """
}

process JOIN_SCANPY_ANNDATA {
    tag "${sample_id}"
    publishDir "${params.outdir}/anndata_annotated", mode: 'copy'

    cpus 2
    memory '4 GB'

    input:
    tuple val(sample_id), path(counts_tsv)
    path h5ad_file

    output:
    path "${sample_id}_annotated.h5ad"

    script:
    """
    python3 -c '
    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad("${h5ad_file}")
    counts = pd.read_csv("${counts_tsv}", sep="\t")

    umi_matrix = counts.pivot(index="cell_barcode", columns="variant_name", values="umi_count").fillna(0)
    umi_matrix = umi_matrix.reindex(adata.obs_names, fill_value=0)
    adata.obsm["X_amplicon_counts"] = umi_matrix.to_numpy()
    adata.uns["amplicon_variants"] = list(umi_matrix.columns)

    adata.write_h5ad("${sample_id}_annotated.h5ad")
    '
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    library_ch    = Channel.fromPath(params.library_tsv, checkIfExists: true).first()
    jar_ch        = Channel.fromPath(params.jar_path, checkIfExists: true).first()

    AMPLICON_FAST_DEMUX(read_pairs_ch, library_ch, jar_ch)

    if (params.scanpy_h5ad) {
        h5ad_ch = Channel.fromPath(params.scanpy_h5ad, checkIfExists: true).first()
        JOIN_SCANPY_ANNDATA(AMPLICON_FAST_DEMUX.out.counts, h5ad_ch)
    }
}
