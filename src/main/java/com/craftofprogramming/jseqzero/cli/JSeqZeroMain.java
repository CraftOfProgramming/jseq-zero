/*
 * Copyright 2026 TechTune Media Group
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.craftofprogramming.jseqzero.cli;

import com.craftofprogramming.jseqzero.tools.AmpliconDemuxTool;
import com.craftofprogramming.jseqzero.tools.SyntheticFastqGenerator;

import java.util.Arrays;

/**
 * Unified Picard/fgbio-style subcommand dispatcher for jseq-zero.
 */
public final class JSeqZeroMain {

    private JSeqZeroMain() {
    }

    public static void main(final String[] args) throws Exception {
        if (args.length == 0 || "-h".equals(args[0]) || "--help".equals(args[0])) {
            printUsage();
            return;
        }

        final String subcommand = args[0];
        final String[] subArgs = Arrays.copyOfRange(args, 1, args.length);

        switch (subcommand) {
            case "amplicon-demux" -> AmpliconDemuxTool.runCli(subArgs);
            case "generate-benchmark" -> SyntheticFastqGenerator.runCli(subArgs);
            default -> {
                System.err.println("Unknown subcommand: " + subcommand);
                printUsage();
                System.exit(1);
            }
        }
    }

    private static void printUsage() {
        System.out.println(
                "jseq-zero v1.0.0 - Zero-Allocation JVM Primitives & Tools for High-Throughput Genomics\n\n"
                + "Usage:\n"
                + "  java -jar jseq-zero.jar <subcommand> [options]\n\n"
                + "Available Subcommands:\n"
                + "  amplicon-demux      Indel-tolerant paired-end FASTQ barcode demultiplexer & 10x UMI deduplicator\n"
                + "                      Args: <library.tsv> <R1.fastq[.gz]> <R2.fastq[.gz]> <out.tsv> [--anchor SEQ]\n\n"
                + "  generate-benchmark  Generate synthetic paired-end FASTQ dataset with PCR duplicates & indels\n"
                + "                      Args: [unique_molecules=1000000] [pcr_depth=5] [num_variants=50]\n");
    }
}
