package com.techtune.jseqzero.tools;

import com.techtune.jseqzero.core.BitapFlankFinder;
import com.techtune.jseqzero.core.Dna2Bit;
import com.techtune.jseqzero.core.FastqRecordReader;
import com.techtune.jseqzero.core.PigeonholeBarcodeIndex;
import org.agrona.collections.Long2LongCounterMap;
import org.agrona.collections.LongHashSet;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * High-throughput single-cell barcoded amplicon demultiplexer and UMI deduplicator.
 *
 * <p>Drop-in compatible with scAAVengr Stage 5 paired-end FASTQ workflows
 * (Read 1: 16-bp 10x Cell Barcode + 12-bp UMI; Read 2: upstream eGFP flank + 25-bp barcode),
 * scaling from 50 variants up to 1,048,575 combinatorial variants in flat primitive memory.
 */
public final class AmpliconDemuxTool {

    private static final String DEFAULT_EGFP_ANCHOR = "ACAGTCGCAATTCCGGTTAA";

    private final BitapFlankFinder flankFinder;
    private final PigeonholeBarcodeIndex barcodeIndex;
    private final String[] variantNames;

    private final LongHashSet seenMolecules;
    private final Long2LongCounterMap cellVariantUmiCounts;

    private long totalReadPairs = 0;
    private long flankNotFound = 0;
    private long barcodeUnmatched = 0;
    private long pcrDuplicates = 0;
    private long uniqueUmisCounted = 0;

    public AmpliconDemuxTool(
            final String anchorSequence,
            final String[] variantNames,
            final long[] referenceBarcodes,
            final int maxBarcodeMismatches,
            final int expectedUniqueMolecules) {

        if (variantNames.length > 0xFFFFF) {
            throw new IllegalArgumentException("20-bit variant packing supports up to 1,048,575 variants");
        }
        this.flankFinder = new BitapFlankFinder(anchorSequence);
        this.variantNames = variantNames;
        this.barcodeIndex = new PigeonholeBarcodeIndex(referenceBarcodes, maxBarcodeMismatches);

        this.seenMolecules = new LongHashSet(expectedUniqueMolecules);
        this.cellVariantUmiCounts = new Long2LongCounterMap(Math.max(1024, expectedUniqueMolecules / 4), 0.65f, 0);
    }

    public void processPairedFastq(final FastqRecordReader r1Reader, final FastqRecordReader r2Reader)
            throws IOException {
        final byte[] r1Buf = new byte[512];
        final byte[] r2Buf = new byte[512];

        while (true) {
            final int r1Len = r1Reader.nextSequence(r1Buf);
            final int r2Len = r2Reader.nextSequence(r2Buf);
            if (r1Len < 0 || r2Len < 0) {
                break;
            }
            totalReadPairs++;

            if (r1Len < 28) {
                continue;
            }
            final long packedCbUmi = Dna2Bit.packCellBarcode16AndUmi12(r1Buf, 0, 16);
            if (packedCbUmi < 0L) {
                continue;
            }

            int barcodeOffset = flankFinder.findDownstreamOffsetFast(r2Buf, r2Len, 25, 0, 40, 1);
            if (barcodeOffset < 0) {
                barcodeOffset = flankFinder.findDownstreamOffsetBitap(r2Buf, r2Len, 25, 0, 40);
                if (barcodeOffset < 0) {
                    flankNotFound++;
                    continue;
                }
            }

            if (barcodeOffset + 25 > r2Len) {
                flankNotFound++;
                continue;
            }

            final long readBarcode = Dna2Bit.encode(r2Buf, barcodeOffset, 25);
            if (readBarcode < 0L) {
                barcodeUnmatched++;
                continue;
            }

            final int variantIdx = barcodeIndex.match(readBarcode);
            if (variantIdx < 0) {
                barcodeUnmatched++;
                continue;
            }

            // 52-bit key: [32-bit Cell Barcode | 20-bit Variant ID] (up to 1,048,575 variants)
            final long cellAndVariantKey = ((packedCbUmi >>> 24) << 20) | (variantIdx & 0xFFFFFL);
            final long moleculeSignature = mix64(packedCbUmi, variantIdx);

            if (seenMolecules.add(moleculeSignature)) {
                uniqueUmisCounted++;
                cellVariantUmiCounts.getAndAdd(cellAndVariantKey, 1);
            } else {
                pcrDuplicates++;
            }
        }
    }

    public void writeOutputTsv(final PrintWriter out) {
        out.println("cell_barcode\tvariant_name\tumi_count");
        final char[] cbChars = new char[16];

        cellVariantUmiCounts.forEach((cellAndVariantKey, umiCount) -> {
            final long cb32 = cellAndVariantKey >>> 20;
            final int variantIdx = (int) (cellAndVariantKey & 0xFFFFFL);

            Dna2Bit.decode(cb32, 16, cbChars, 0);
            out.print(cbChars);
            out.print('\t');
            out.print(variantNames[variantIdx]);
            out.print('\t');
            out.println(umiCount);
        });
        out.flush();
    }

    public void printTelemetry(final double elapsedSeconds) {
        final long usedMemMb =
                (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / (1024 * 1024);
        System.out.printf(
                "Total Pairs: %,d | Unique UMIs: %,d | PCR Dups: %,d | No Flank: %,d | Unmatched Barcode: %,d%n",
                totalReadPairs, uniqueUmisCounted, pcrDuplicates, flankNotFound, barcodeUnmatched);
        System.out.printf(
                "jseq-zero Elapsed: %.2f sec (%,.0f pairs/sec) | Approx Heap Used: %d MB%n",
                elapsedSeconds, totalReadPairs / Math.max(0.001, elapsedSeconds), usedMemMb);
    }

    public long getUniqueUmisCounted() {
        return uniqueUmisCounted;
    }

    private static long mix64(final long packedCbUmi, final int variantIdx) {
        long h = packedCbUmi ^ ((long) variantIdx * 0x9E3779B97F4A7C15L);
        h = (h ^ (h >>> 30)) * 0xBF58476D1CE4E5B9L;
        h = (h ^ (h >>> 27)) * 0x94D049BB133111EBL;
        return h ^ (h >>> 31);
    }

    public static void runCli(final String[] args) throws IOException {
        if (args.length < 4) {
            System.err.println(
                    "Usage: jseq-zero amplicon-demux <library_barcodes.tsv> <R1.fastq[.gz]> <R2.fastq[.gz]> <out_counts.tsv> [--anchor SEQ]");
            System.exit(1);
        }

        String anchor = DEFAULT_EGFP_ANCHOR;
        for (int i = 4; i < args.length - 1; i++) {
            if ("--anchor".equals(args[i])) {
                anchor = args[i + 1];
            }
        }

        final long startNs = System.nanoTime();

        final List<String> namesList = new ArrayList<>();
        final List<Long> barcodesList = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
            String line = br.readLine(); // skip TSV header
            while ((line = br.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                final int tab = line.indexOf('\t');
                namesList.add(line.substring(0, tab));
                final byte[] bcBytes = line.substring(tab + 1).trim().getBytes(StandardCharsets.US_ASCII);
                barcodesList.add(Dna2Bit.encode(bcBytes, 0, 25));
            }
        }

        final String[] variantNames = namesList.toArray(new String[0]);
        final long[] refBarcodes = barcodesList.stream().mapToLong(Long::longValue).toArray();

        final AmpliconDemuxTool tool = new AmpliconDemuxTool(
                anchor,
                variantNames,
                refBarcodes,
                2,
                2_000_000
        );

        try (FastqRecordReader r1 = FastqRecordReader.open(args[1]);
             FastqRecordReader r2 = FastqRecordReader.open(args[2]);
             PrintWriter out = new PrintWriter(new BufferedOutputStream(new FileOutputStream(args[3]), 64 * 1024))) {
            tool.processPairedFastq(r1, r2);
            tool.writeOutputTsv(out);
        }

        final double elapsedSec = (System.nanoTime() - startNs) / 1_000_000_000.0;
        tool.printTelemetry(elapsedSec);
    }
}
