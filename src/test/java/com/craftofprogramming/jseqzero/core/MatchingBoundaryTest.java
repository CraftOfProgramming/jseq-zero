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

import org.junit.jupiter.api.Test;
import java.nio.charset.StandardCharsets;
import static org.junit.jupiter.api.Assertions.*;

class MatchingBoundaryTest {
    private static byte[] ascii(String value) { return value.getBytes(StandardCharsets.US_ASCII); }
    private static long encoded(String value) { return Dna2Bit.encode(ascii(value), 0, value.length()); }

    @Test void encodesLowercaseAndRejectsAmbiguousMolecules() {
        assertEquals(encoded("ACGT"), encoded("acgt"));
        assertEquals(-1, encoded("AN"));
        assertEquals(-1, encoded("an"));
        assertEquals(-1, Dna2Bit.packCellBarcode16AndUmi12(ascii("A".repeat(16)+"N"+"A".repeat(11)), 0, 16));
        String value = "ACGT".repeat(7)+"ACG";
        char[] result = new char[31];
        Dna2Bit.decode(encoded(value), 31, result, 0);
        assertEquals(value, new String(result));
    }

    @Test void rejectsAmbiguousBestMatchWithoutConfusingRepeatedSeedCandidates() {
        var index = new PigeonholeBarcodeIndex(new long[]{encoded("A".repeat(25)), encoded("C"+"A".repeat(24)), encoded("G".repeat(25))}, 2);
        assertEquals(-1, index.match(encoded("T"+"A".repeat(24))));
        assertEquals(0, index.match(encoded("A".repeat(24)+"C")));
        assertEquals(2, index.match(encoded("G".repeat(25))));
        assertEquals(-1, new PigeonholeBarcodeIndex(new long[0], 2).match(0));
        assertThrows(IllegalArgumentException.class, () -> new PigeonholeBarcodeIndex(new long[0], 3));
    }

    @Test void findsExactSubstitutedAndInsertedAnchors() {
        String anchor = "ACAGTCGCAATTCCGGTTAA";
        var finder = new BitapFlankFinder(anchor);
        byte[] exact = ascii("TTTGGCC"+anchor+"C".repeat(25));
        assertEquals(27, finder.findDownstreamOffsetFast(exact, exact.length, 25, 0, 40, 1));
        byte[] changed = ascii("TTTGGCC"+"TCAGTCGCAATTCCGGTTAA"+"C".repeat(25));
        assertEquals(27, finder.findDownstreamOffsetFast(changed, changed.length, 25, 0, 40, 1));
        byte[] inserted = ascii("TTTGGCC"+anchor.substring(0,8)+"G"+anchor.substring(8)+"C".repeat(25));
        assertEquals(28, finder.findDownstreamOffsetBitap(inserted, inserted.length, 25, 0, 40));
    }

    @Test void handlesMissingAnchorsAndInsufficientDownstreamPayload() {
        var finder = new BitapFlankFinder("ACAGTCGCAATTCCGGTTAA");
        byte[] missing = ascii("G".repeat(70));
        assertEquals(-1, finder.findDownstreamOffsetFast(missing, missing.length, 25, 0, 40, 1));
        assertEquals(-1, finder.findDownstreamOffsetBitap(missing, missing.length, 25, 0, 40));
        assertEquals(-1, finder.findDownstreamOffsetFast(ascii("ACGT"), 4, 25, 0, 40, 1));
        assertThrows(IllegalArgumentException.class, () -> new BitapFlankFinder("ACG"));
        assertThrows(IllegalArgumentException.class, () -> new BitapFlankFinder("ACGN"));
    }
}
