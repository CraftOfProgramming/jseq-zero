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
package com.techtune.jseqzero.core;

/**
 * Branchless 2-bit DNA sequence encoding and register-level Hamming distance primitives.
 *
 * <p>Uses bits 1 and 2 of ASCII bytes {@code (b >> 1) & 3L} to map bases without lookup tables:
 * {@code A -> 0 (00)}, {@code C -> 1 (01)}, {@code T -> 2 (10)}, {@code G -> 3 (11)}.
 */
public final class Dna2Bit {

    /** Alternating bitmask (01010101...) isolating bit 0 of all 32 2-bit base slots in a 64-bit word. */
    public static final long EVEN_BITS_MASK = 0x5555555555555555L;

    private static final char[] TWO_BIT_TO_ASCII = {'A', 'C', 'T', 'G'};

    private Dna2Bit() {
    }

    /**
     * Encodes up to 31 DNA bases from a raw ASCII byte slice into a 64-bit primitive {@code long}.
     *
     * @param buffer raw ASCII byte array
     * @param offset starting byte index
     * @param length number of bases to encode (1..31)
     * @return non-negative 2-bit packed sequence, or {@code -1L} if an ambiguous base ('N') is present
     */
    public static long encode(final byte[] buffer, final int offset, final int length) {
        long encoded = 0L;
        final int end = offset + length;
        for (int i = offset; i < end; i++) {
            final byte b = buffer[i];
            if (b == 'N' || b == 'n') {
                return -1L;
            }
            encoded = (encoded << 2) | ((b >> 1) & 3L);
        }
        return encoded;
    }

    /**
     * Decodes a 2-bit packed DNA sequence back into a pre-allocated char buffer (zero allocation).
     */
    public static void decode(final long encoded, final int length, final char[] target, final int offset) {
        for (int i = length - 1; i >= 0; i--) {
            final int twoBits = (int) ((encoded >>> (i * 2)) & 3L);
            target[offset + (length - 1 - i)] = TWO_BIT_TO_ASCII[twoBits];
        }
    }

    /**
     * Computes exact base-level Hamming distance between two equal-length 2-bit encoded DNA sequences
     * in 4 CPU instructions (XOR, SHR, OR, AND, POPCNT).
     */
    public static int hammingDistance(final long seqA, final long seqB) {
        final long xor = seqA ^ seqB;
        final long collapsed = (xor | (xor >>> 1)) & EVEN_BITS_MASK;
        return Long.bitCount(collapsed);
    }

    /**
     * Packs a 16-bp Cell Barcode (32 bits) and 12-bp UMI (24 bits) from Read 1 into a 56-bit {@code long}.
     *
     * @return packed 56-bit value {@code [CB_32 | UMI_24]}, or {@code -1L} if any base is 'N'
     */
    public static long packCellBarcode16AndUmi12(
            final byte[] r1Buffer, final int cbOffset, final int umiOffset) {
        final long cb16 = encode(r1Buffer, cbOffset, 16);
        final long umi12 = encode(r1Buffer, umiOffset, 12);
        if ((cb16 | umi12) < 0L) {
            return -1L;
        }
        return (cb16 << 24) | umi12;
    }
}
