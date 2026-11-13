package com.techtune.jseqzero.core;

import org.junit.jupiter.api.Test;

import java.nio.charset.StandardCharsets;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class CorePrimitivesTest {

    @Test
    void encodesAndComputesHammingDistance() {
        final byte[] ref = "ACGTACGTACGTACGTACGTACGTA".getBytes(StandardCharsets.US_ASCII);
        final byte[] mut = "ACGTACGTCCGTACGTACGTACGTT".getBytes(StandardCharsets.US_ASCII);

        final long refEnc = Dna2Bit.encode(ref, 0, 25);
        final long mutEnc = Dna2Bit.encode(mut, 0, 25);

        assertEquals(0, Dna2Bit.hammingDistance(refEnc, refEnc));
        assertEquals(2, Dna2Bit.hammingDistance(refEnc, mutEnc));

        final char[] decoded = new char[25];
        Dna2Bit.decode(refEnc, 25, decoded, 0);
        assertEquals("ACGTACGTACGTACGTACGTACGTA", new String(decoded));
    }

    @Test
    void bitapRecoversBarcodeOffsetAcrossAnchorDeletion() {
        final BitapFlankFinder finder = new BitapFlankFinder("ACAGTCGCAATTCCGGTTAA");

        // 7-bp 5' leader + 19-bp anchor (1-bp deletion of 'A') + 25-bp barcode
        final String read2 = "TTTGGCC" + "ACAGTCGCATTCCGGTTAA" + "ACGTACGTACGTACGTACGTACGTA" + "AAAAA";
        final byte[] r2Bytes = read2.getBytes(StandardCharsets.US_ASCII);

        final int offset = finder.findDownstreamOffsetBitap(r2Bytes, r2Bytes.length, 25, 0, 40);
        assertEquals(26, offset);

        final long extracted = Dna2Bit.encode(r2Bytes, offset, 25);
        final long expected = Dna2Bit.encode("ACGTACGTACGTACGTACGTACGTA".getBytes(StandardCharsets.US_ASCII), 0, 25);
        assertEquals(expected, extracted);
    }

    @Test
    void pigeonholeIndexMatchesWithTwoSubstitutions() {
        final String[] seqs = {
                "AAAAAAAAAAAAAAAAAAAAAAAAA",
                "CCCCCCCCCCCCCCCCCCCCCCCCC",
                "ACGTACGTACGTACGTACGTACGTA",
                "TTTTTTTTTTTTTTTTTTTTTTTTT"
        };
        final long[] encoded = new long[seqs.length];
        for (int i = 0; i < seqs.length; i++) {
            encoded[i] = Dna2Bit.encode(seqs[i].getBytes(StandardCharsets.US_ASCII), 0, 25);
        }

        final PigeonholeBarcodeIndex index = new PigeonholeBarcodeIndex(encoded, 2);
        final long queryWith2Muts = Dna2Bit.encode(
                "ACGTACGTCCGTACGTACGTACGTT".getBytes(StandardCharsets.US_ASCII), 0, 25);

        assertEquals(2, index.match(queryWith2Muts));
        assertTrue(index.match(Dna2Bit.encode("GGGGGGGGGGGGGGGGGGGGGGGGG".getBytes(StandardCharsets.US_ASCII), 0, 25)) < 0);
    }
}
