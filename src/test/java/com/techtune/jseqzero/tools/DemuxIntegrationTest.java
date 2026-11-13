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

import com.techtune.jseqzero.cli.JSeqZeroMain;
import com.techtune.jseqzero.core.Dna2Bit;
import com.techtune.jseqzero.core.FastqRecordReader;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;
import static org.junit.jupiter.api.Assertions.*;

class DemuxIntegrationTest {
    @TempDir Path directory;
    private static final String CELL = "ACGTACGTACGTACGT";
    private static final String ANCHOR = "ACAGTCGCAATTCCGGTTAA";

    private static String record(String sequence) {
        return "@r\n"+sequence+"\n+\n"+"I".repeat(sequence.length())+"\n";
    }
    private static FastqRecordReader reader(String text) {
        return new FastqRecordReader(new ByteArrayInputStream(text.getBytes(StandardCharsets.US_ASCII)));
    }
    private static long encoded(String value) {
        return Dna2Bit.encode(value.getBytes(StandardCharsets.US_ASCII), 0, 25);
    }
    private static AmpliconDemuxTool tool() {
        return new AmpliconDemuxTool(ANCHOR,new String[]{"v0","v1"},new long[]{encoded("A".repeat(25)),encoded("C".repeat(25))},2,32);
    }

    @Test void countsMoleculesByCellVariantAndUmiAndRejectsAmbiguousRead1() throws Exception {
        String r1 = record(CELL+"A".repeat(12)).repeat(2)
                +record(CELL+"C".repeat(12))+record(CELL+"A".repeat(12))
                +record("T".repeat(16)+"A".repeat(12))+record(CELL+"N"+"A".repeat(11));
        String a = record("TTTGGCC"+ANCHOR+"A".repeat(25));
        String c = record("TTTGGCC"+ANCHOR+"C".repeat(25));
        var tool = tool();
        try (var first=reader(r1);var second=reader(a.repeat(3)+c+a.repeat(2))) {
            tool.processPairedFastq(first,second);
        }
        assertEquals(4,tool.getUniqueUmisCounted());
        StringWriter output = new StringWriter();
        tool.writeOutputTsv(new PrintWriter(output));
        assertTrue(output.toString().contains(CELL+"\tv0\t2"));
        assertTrue(output.toString().contains(CELL+"\tv1\t1"));
        assertTrue(output.toString().contains("T".repeat(16)+"\tv0\t1"));
    }

    @Test void refusesUnequalPairCounts() throws Exception {
        try (var first=reader(record(CELL+"A".repeat(12)));var second=reader("")) {
            assertThrows(IOException.class, () -> tool().processPairedFastq(first,second));
        }
    }

    @Test void cliDispatchesGzipInputAndWritesExpectedTsv() throws Exception {
        Path library=directory.resolve("library.tsv"),r1=directory.resolve("R1.fastq.gz"),r2=directory.resolve("R2.fastq.gz"),out=directory.resolve("counts.tsv");
        Files.writeString(library,"variant_name\tbarcode_25bp\nv0\t"+"A".repeat(25)+"\n");
        try(var stream=new GZIPOutputStream(Files.newOutputStream(r1))) {stream.write(record(CELL+"A".repeat(12)).getBytes(StandardCharsets.US_ASCII));}
        try(var stream=new GZIPOutputStream(Files.newOutputStream(r2))) {stream.write(record("TTTGGCC"+ANCHOR+"A".repeat(25)).getBytes(StandardCharsets.US_ASCII));}
        JSeqZeroMain.main(new String[]{"amplicon-demux",library.toString(),r1.toString(),r2.toString(),out.toString(),"--anchor",ANCHOR});
        assertEquals("cell_barcode\tvariant_name\tumi_count\n"+CELL+"\tv0\t1\n",Files.readString(out));
    }

    @Test void helpExplainsBothSubcommands() throws Exception {
        PrintStream previous=System.out;
        ByteArrayOutputStream buffer=new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, StandardCharsets.UTF_8));
            JSeqZeroMain.main(new String[]{"--help"});
        } finally {System.setOut(previous);}
        assertTrue(buffer.toString(StandardCharsets.UTF_8).contains("amplicon-demux"));
        assertTrue(buffer.toString(StandardCharsets.UTF_8).contains("generate-benchmark"));
    }
}
