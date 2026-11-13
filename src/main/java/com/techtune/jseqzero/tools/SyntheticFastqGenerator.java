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
package com.techtune.jseqzero.tools;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.SplittableRandom;

/**
 * Generates synthetic paired-end amplicon FASTQ benchmarks with configurable library size,
 * PCR duplication depth, upstream frame-shift indels, anchor deletions, and barcode substitutions.
 */
public final class SyntheticFastqGenerator {

    private static final byte[] BASES = {'A', 'C', 'G', 'T'};
    private static final byte[] EGFP_ANCHOR = "ACAGTCGCAATTCCGGTTAA".getBytes(StandardCharsets.US_ASCII);
    private static final byte[] POLY_A_TAIL = "AAAAAAAAAAAAAAAAAAAAAAAAA".getBytes(StandardCharsets.US_ASCII);

    private SyntheticFastqGenerator() {
    }

    public static void runCli(final String[] args) throws IOException {
        final int numUniqueMolecules = args.length > 0 ? Integer.parseInt(args[0]) : 1_000_000;
        final int pcrFactor = args.length > 1 ? Integer.parseInt(args[1]) : 5;
        final int numVariants = args.length > 2 ? Integer.parseInt(args[2]) : 50;
        generate(Path.of("."), numUniqueMolecules, pcrFactor, numVariants);
    }

    /** Generates deterministic synthetic inputs in the supplied directory. */
    public static void generate(final Path outputDirectory, final int numUniqueMolecules,
                                final int pcrFactor, final int numVariants) throws IOException {
        if (numUniqueMolecules < 1 || pcrFactor < 1 || numVariants < 1 || numVariants > 0xFFFFF) {
            throw new IllegalArgumentException("Molecules and PCR depth must be positive; variants must be in [1..1048575]");
        }
        final int totalReadPairs = Math.multiplyExact(numUniqueMolecules, pcrFactor);
        final int numCells = 5_000;

        final SplittableRandom rng = new SplittableRandom(42L);

        final byte[][] barcodeLibrary = new byte[numVariants][25];
        try (PrintWriter libOut = new PrintWriter(outputDirectory.resolve("library_barcodes.tsv").toFile())) {
            libOut.println("variant_name\tbarcode_25bp");
            for (int v = 0; v < numVariants; v++) {
                for (int i = 0; i < 25; i++) {
                    barcodeLibrary[v][i] = BASES[rng.nextInt(4)];
                }
                libOut.printf("AAV_Variant_%05d\t%s%n", v, new String(barcodeLibrary[v], StandardCharsets.US_ASCII));
            }
        }

        final byte[][] cellBarcodes = new byte[numCells][16];
        for (int c = 0; c < numCells; c++) {
            for (int i = 0; i < 16; i++) {
                cellBarcodes[c][i] = BASES[rng.nextInt(4)];
            }
        }

        final int[] molCellIdx = new int[numUniqueMolecules];
        final int[] molVariantIdx = new int[numUniqueMolecules];
        final byte[][] molUmi = new byte[numUniqueMolecules][12];

        for (int m = 0; m < numUniqueMolecules; m++) {
            molCellIdx[m] = rng.nextInt(numCells);
            molVariantIdx[m] = rng.nextInt(numVariants);
            for (int i = 0; i < 12; i++) {
                molUmi[m][i] = BASES[rng.nextInt(4)];
            }
        }

        System.out.printf(
                "Generating %,d paired-end reads (%,d unique molecules, %dx PCR depth, %,d library variants)...%n",
                totalReadPairs, numUniqueMolecules, pcrFactor, numVariants);

        try (OutputStream r1Out = new BufferedOutputStream(new FileOutputStream(outputDirectory.resolve("R1.fastq").toFile()), 1024 * 1024);
             OutputStream r2Out = new BufferedOutputStream(new FileOutputStream(outputDirectory.resolve("R2.fastq").toFile()), 1024 * 1024)) {

            final byte[] r1Qual = new byte[28];
            final byte[] r2Qual = new byte[100];
            Arrays.fill(r1Qual, (byte) 'I');
            Arrays.fill(r2Qual, (byte) 'I');

            final byte[] r2Buf = new byte[128];

            for (int r = 0; r < totalReadPairs; r++) {
                final int molIdx = rng.nextInt(numUniqueMolecules);
                final byte[] cb = cellBarcodes[molCellIdx[molIdx]];
                final byte[] umi = molUmi[molIdx];
                final byte[] bc = barcodeLibrary[molVariantIdx[molIdx]];

                writeFastqHeader(r1Out, r);
                r1Out.write(cb);
                r1Out.write(umi);
                r1Out.write('\n');
                r1Out.write('+');
                r1Out.write('\n');
                r1Out.write(r1Qual);
                r1Out.write('\n');

                int pos = 0;
                final int leaderLen = 4 + rng.nextInt(15);
                for (int i = 0; i < leaderLen; i++) {
                    r2Buf[pos++] = BASES[rng.nextInt(4)];
                }

                // 4% 1-bp deletion inside anchor, 4% 1-bp substitution inside anchor
                final int anchorNoise = rng.nextInt(100);
                final int dropIdx = (anchorNoise < 4) ? rng.nextInt(20) : -1;
                final int subIdx = (anchorNoise >= 4 && anchorNoise < 8) ? rng.nextInt(20) : -1;

                for (int i = 0; i < EGFP_ANCHOR.length; i++) {
                    if (i == dropIdx) {
                        continue;
                    }
                    byte b = EGFP_ANCHOR[i];
                    if (i == subIdx) {
                        b = mutateBase(b);
                    }
                    r2Buf[pos++] = b;
                }

                // 25% of reads carry 1-2 point mutations in the 25-bp barcode
                final int bcNoise = rng.nextInt(100);
                final int mut1 = (bcNoise < 25) ? rng.nextInt(25) : -1;
                final int mut2 = (bcNoise < 5) ? (mut1 + 1 + rng.nextInt(24)) % 25 : -1;

                for (int i = 0; i < 25; i++) {
                    byte b = bc[i];
                    if (i == mut1 || i == mut2) {
                        b = mutateBase(b);
                    }
                    r2Buf[pos++] = b;
                }

                System.arraycopy(POLY_A_TAIL, 0, r2Buf, pos, POLY_A_TAIL.length);
                pos += POLY_A_TAIL.length;

                writeFastqHeader(r2Out, r);
                r2Out.write(r2Buf, 0, pos);
                r2Out.write('\n');
                r2Out.write('+');
                r2Out.write('\n');
                r2Out.write(r2Qual, 0, pos);
                r2Out.write('\n');
            }
        }
        System.out.println("Completed: wrote library_barcodes.tsv, R1.fastq, and R2.fastq.");
    }

    private static byte mutateBase(final byte b) {
        return switch (b) {
            case 'A' -> 'C';
            case 'C' -> 'G';
            case 'G' -> 'T';
            default -> 'A';
        };
    }

    private static void writeFastqHeader(final OutputStream out, final int readId) throws IOException {
        out.write('@');
        out.write('R');
        final String s = Integer.toString(readId);
        for (int i = 0; i < s.length(); i++) {
            out.write(s.charAt(i));
        }
        out.write('\n');
    }
}
