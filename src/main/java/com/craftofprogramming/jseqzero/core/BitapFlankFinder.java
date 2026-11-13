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
package com.craftofprogramming.jseqzero.core;

import java.nio.charset.StandardCharsets;

/**
 * Two-tier constant flanking/primer sequence detector for noisy amplicon reads.
 *
 * <p>Cascades a 6-cycle 2-bit sliding window (handling upstream frame-shift indels and substitutions)
 * into a branchless Wu-Manber / Bitap {@code k <= 2} Levenshtein bit-parallel automaton to recover
 * valid downstream barcodes even when PCR introduces insertions or deletions inside the anchor itself.
 */
public final class BitapFlankFinder {

    private final int anchorLength;
    private final long encodedAnchor;
    private final long anchorWindowMask;

    private final long[] bitapPatternMask = new long[128];
    private final long bitapMatchBit;

    public BitapFlankFinder(final String flankSequence) {
        this.anchorLength = flankSequence.length();
        if (anchorLength < 4 || anchorLength > 31) {
            throw new IllegalArgumentException("Anchor length must be between 4 and 31 bp");
        }

        final byte[] anchorBytes = flankSequence.getBytes(StandardCharsets.US_ASCII);
        this.encodedAnchor = Dna2Bit.encode(anchorBytes, 0, anchorLength);
        if (this.encodedAnchor < 0L) {
            throw new IllegalArgumentException("Anchor sequence must not contain ambiguous 'N' bases");
        }
        this.anchorWindowMask = (1L << (anchorLength * 2)) - 1L;

        for (int i = 0; i < 128; i++) {
            bitapPatternMask[i] = ~0L;
        }
        for (int i = 0; i < anchorLength; i++) {
            final byte b = anchorBytes[i];
            final long bit = ~(1L << i);
            bitapPatternMask[b] &= bit;
            bitapPatternMask[b | 0x20] &= bit; // case-insensitive ASCII
        }
        this.bitapMatchBit = 1L << (anchorLength - 1);
    }

    /**
     * Fast-path: Slides a 2-bit encoded window across {@code [minStartOffset .. maxStartOffset]}
     * allowing up to {@code maxSubstitutions} mismatches inside the anchor.
     *
     * @return exact starting byte offset of the downstream payload/barcode, or {@code -1} if not found
     */
    public int findDownstreamOffsetFast(
            final byte[] readBuffer,
            final int readLength,
            final int downstreamMinBytes,
            final int minStartOffset,
            final int maxStartOffset,
            final int maxSubstitutions) {

        final int maxSearchEnd = Math.min(readLength - downstreamMinBytes, maxStartOffset + anchorLength);
        if (minStartOffset + anchorLength > maxSearchEnd) {
            return -1;
        }

        long rollingWindow = 0L;
        int ambiguousRun = 0;
        for (int i = minStartOffset; i < minStartOffset + anchorLength - 1; i++) {
            final byte b = readBuffer[i];
            if (b == 'N' || b == 'n') {
                ambiguousRun = anchorLength;
            } else if (ambiguousRun > 0) {
                ambiguousRun--;
            }
            rollingWindow = (rollingWindow << 2) | ((b >> 1) & 3L);
        }

        int bestDownstreamOffset = -1;
        int bestDist = maxSubstitutions + 1;

        for (int i = minStartOffset + anchorLength - 1; i < maxSearchEnd; i++) {
            final byte b = readBuffer[i];
            if (b == 'N' || b == 'n') {
                ambiguousRun = anchorLength;
            } else if (ambiguousRun > 0) {
                ambiguousRun--;
            }
            rollingWindow = ((rollingWindow << 2) | ((b >> 1) & 3L)) & anchorWindowMask;

            if (ambiguousRun > 0) {
                continue;
            }

            final long xor = rollingWindow ^ encodedAnchor;
            final int dist = Long.bitCount((xor | (xor >>> 1)) & Dna2Bit.EVEN_BITS_MASK);
            if (dist == 0) {
                return i + 1;
            }
            if (dist < bestDist) {
                bestDist = dist;
                bestDownstreamOffset = i + 1;
            }
        }
        return bestDownstreamOffset;
    }

    /**
     * Fallback path: Branchless Wu-Manber / Bitap Levenshtein automaton ({@code k <= 2} edit distance)
     * handling insertions, deletions, and substitutions inside the flanking anchor.
     *
     * @return exact starting byte offset of the downstream payload/barcode, or {@code -1} if not found
     */
    public int findDownstreamOffsetBitap(
            final byte[] readBuffer,
            final int readLength,
            final int downstreamMinBytes,
            final int searchStart,
            final int searchEnd) {

        long r0 = ~0L;
        long r1 = ~0L << 1;
        long r2 = ~0L << 2;

        final int limit = Math.min(readLength - downstreamMinBytes, searchEnd);
        int bestDownstreamOffset = -1;
        int bestEditDist = 3;

        for (int i = searchStart; i < limit; i++) {
            final long charMask = bitapPatternMask[readBuffer[i] & 0x7F];

            final long oldR0 = r0;
            final long oldR1 = r1;

            r0 = (r0 << 1) | charMask;
            r1 = ((r1 << 1) | charMask) & (oldR0 & (r0 << 1)) & (oldR0 << 1);
            r2 = ((r2 << 1) | charMask) & (oldR1 & (r1 << 1)) & (oldR1 << 1);

            if ((r0 & bitapMatchBit) == 0L) {
                return i + 1;
            } else if ((r1 & bitapMatchBit) == 0L && bestEditDist > 1) {
                bestEditDist = 1;
                bestDownstreamOffset = i + 1;
            } else if ((r2 & bitapMatchBit) == 0L && bestEditDist > 2) {
                bestEditDist = 2;
                bestDownstreamOffset = i + 1;
            }
        }
        return bestDownstreamOffset;
    }
}
