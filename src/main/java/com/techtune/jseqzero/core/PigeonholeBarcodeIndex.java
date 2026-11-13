package com.techtune.jseqzero.core;

import java.util.Arrays;

/**
 * Compressed Sparse Row (CSR) 3-slice Pigeonhole Principle seed-and-extend index for 25-bp barcodes.
 *
 * <p>Partitions a 50-bit encoded 25-bp barcode into three disjoint slices (9 bp / 18 bits,
 * 8 bp / 16 bits, and 8 bp / 16 bits). Any read with {@code <= 2} substitutions must match
 * at least one slice with 0 errors, enabling {@code O(1)} expected matching across 100,000+ variants
 * with zero heap object overhead (~2.7 MB L3-cache-resident footprint).
 */
public final class PigeonholeBarcodeIndex {

    private static final int SLICE0_SPACE = 1 << 18; // 9 bp = 18 bits (262,144 buckets)
    private static final int SLICE1_SPACE = 1 << 16; // 8 bp = 16 bits (65,536 buckets)
    private static final int SLICE2_SPACE = 1 << 16; // 8 bp = 16 bits (65,536 buckets)

    private static final long MASK_18_BITS = SLICE0_SPACE - 1L;
    private static final long MASK_16_BITS = SLICE1_SPACE - 1L;

    private final long[] referenceBarcodes;
    private final int maxMismatches;

    private final int[] offsets0 = new int[SLICE0_SPACE + 1];
    private final int[] candidates0;

    private final int[] offsets1 = new int[SLICE1_SPACE + 1];
    private final int[] candidates1;

    private final int[] offsets2 = new int[SLICE2_SPACE + 1];
    private final int[] candidates2;

    public PigeonholeBarcodeIndex(final long[] referenceBarcodes, final int maxMismatches) {
        if (maxMismatches < 0 || maxMismatches > 2) {
            throw new IllegalArgumentException("3-slice pigeonhole index requires maxMismatches in [0..2]");
        }
        this.referenceBarcodes = referenceBarcodes;
        this.maxMismatches = maxMismatches;

        final int n = referenceBarcodes.length;
        this.candidates0 = new int[n];
        this.candidates1 = new int[n];
        this.candidates2 = new int[n];

        buildSliceIndex(referenceBarcodes, 0, MASK_18_BITS, offsets0, candidates0);
        buildSliceIndex(referenceBarcodes, 18, MASK_16_BITS, offsets1, candidates1);
        buildSliceIndex(referenceBarcodes, 34, MASK_16_BITS, offsets2, candidates2);
    }

    private static void buildSliceIndex(
            final long[] barcodes,
            final int shift,
            final long mask,
            final int[] offsets,
            final int[] candidates) {

        for (int i = 0; i < barcodes.length; i++) {
            final int key = (int) ((barcodes[i] >>> shift) & mask);
            offsets[key + 1]++;
        }
        for (int k = 0; k < offsets.length - 1; k++) {
            offsets[k + 1] += offsets[k];
        }
        final int[] writeCursors = Arrays.copyOf(offsets, offsets.length);
        for (int variantIdx = 0; variantIdx < barcodes.length; variantIdx++) {
            final int key = (int) ((barcodes[variantIdx] >>> shift) & mask);
            candidates[writeCursors[key]++] = variantIdx;
        }
    }

    /**
     * Matches an encoded 25-bp barcode against the indexed library within {@code maxMismatches}.
     *
     * @return unique best-matching variant index, or {@code -1} if unmatched or tied
     */
    public int match(final long readBarcode) {
        final int key0 = (int) (readBarcode & MASK_18_BITS);
        final int key1 = (int) ((readBarcode >>> 18) & MASK_16_BITS);
        final int key2 = (int) ((readBarcode >>> 34) & MASK_16_BITS);

        int bestVariant = -1;
        int bestDist = maxMismatches + 1;
        boolean tie = false;

        for (int i = offsets0[key0], end = offsets0[key0 + 1]; i < end; i++) {
            final int candidateIdx = candidates0[i];
            final int dist = Dna2Bit.hammingDistance(readBarcode, referenceBarcodes[candidateIdx]);
            if (dist < bestDist) {
                bestDist = dist;
                bestVariant = candidateIdx;
                tie = false;
                if (dist == 0) {
                    return candidateIdx;
                }
            } else if (dist == bestDist && candidateIdx != bestVariant) {
                tie = true;
            }
        }

        for (int i = offsets1[key1], end = offsets1[key1 + 1]; i < end; i++) {
            final int candidateIdx = candidates1[i];
            if (candidateIdx == bestVariant) {
                continue;
            }
            final int dist = Dna2Bit.hammingDistance(readBarcode, referenceBarcodes[candidateIdx]);
            if (dist < bestDist) {
                bestDist = dist;
                bestVariant = candidateIdx;
                tie = false;
                if (dist == 0) {
                    return candidateIdx;
                }
            } else if (dist == bestDist) {
                tie = true;
            }
        }

        for (int i = offsets2[key2], end = offsets2[key2 + 1]; i < end; i++) {
            final int candidateIdx = candidates2[i];
            if (candidateIdx == bestVariant) {
                continue;
            }
            final int dist = Dna2Bit.hammingDistance(readBarcode, referenceBarcodes[candidateIdx]);
            if (dist < bestDist) {
                bestDist = dist;
                bestVariant = candidateIdx;
                tie = false;
                if (dist == 0) {
                    return candidateIdx;
                }
            } else if (dist == bestDist) {
                tie = true;
            }
        }

        return tie ? -1 : bestVariant;
    }
}
