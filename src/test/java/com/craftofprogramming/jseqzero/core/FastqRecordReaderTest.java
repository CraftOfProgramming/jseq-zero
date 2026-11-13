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
import org.junit.jupiter.api.io.TempDir;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;
import static org.junit.jupiter.api.Assertions.*;

class FastqRecordReaderTest {
    @TempDir Path directory;

    private static FastqRecordReader reader(String text) {
        return new FastqRecordReader(new ByteArrayInputStream(text.getBytes(StandardCharsets.US_ASCII)));
    }

    @Test void readsLfAndCrLfRecordsIntoReusableBuffer() throws Exception {
        try (var input = reader("@r1\nACGT\n+\nIIII\n@r2\r\nTGCA\r\n+\r\nIIII\r\n")) {
            byte[] buffer = new byte[4];
            assertEquals(4, input.nextSequence(buffer));
            assertEquals("ACGT", new String(buffer, StandardCharsets.US_ASCII));
            assertEquals(4, input.nextSequence(buffer));
            assertEquals("TGCA", new String(buffer, StandardCharsets.US_ASCII));
            assertEquals(-1, input.nextSequence(buffer));
        }
    }

    @Test void acceptsFinalQualityLineWithoutTrailingNewline() throws Exception {
        try (var input = reader("@r\nACGT\n+\nIIII")) {
            assertEquals(4, input.nextSequence(new byte[4]));
            assertEquals(-1, input.nextSequence(new byte[4]));
        }
    }

    @Test void rejectsTruncatedRecordsInsteadOfReportingCleanEof() {
        for (String data : new String[]{"@r\n", "@r\nACGT\n", "@r\nACGT\n+\n"}) {
            try (var input = reader(data)) {
                assertThrows(IOException.class, () -> input.nextSequence(new byte[4]));
            } catch (IOException error) { fail(error); }
        }
    }

    @Test void rejectsSequencesLargerThanBufferInsteadOfTruncating() throws Exception {
        try (var input = reader("@r\nACGTA\n+\nIIIII\n")) {
            assertThrows(IOException.class, () -> input.nextSequence(new byte[4]));
        }
    }

    @Test void readsPlainAndGzipFilesIdentically() throws Exception {
        byte[] data = "@r\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII);
        Path plain = directory.resolve("reads.fastq"), gzip = directory.resolve("reads.fastq.gz");
        Files.write(plain, data);
        try (var out = new GZIPOutputStream(Files.newOutputStream(gzip))) { out.write(data); }
        for (Path path : new Path[]{plain, gzip}) {
            try (var input = FastqRecordReader.open(path.toString())) {
                byte[] buffer = new byte[4];
                assertEquals(4, input.nextSequence(buffer));
                assertArrayEquals("ACGT".getBytes(StandardCharsets.US_ASCII), buffer);
                assertEquals(-1, input.nextSequence(buffer));
            }
        }
    }

    @Test void readsRecordsAcrossInternalBufferBoundaries() throws Exception {
        String sequence = "ACGT".repeat(25);
        String record = "@read\n" + sequence + "\n+\n" + "I".repeat(100) + "\n";
        try (var input = reader(record.repeat(1000))) {
            byte[] buffer = new byte[100];
            for (int i = 0; i < 1000; i++) {
                assertEquals(100, input.nextSequence(buffer));
                assertEquals(sequence, new String(buffer, StandardCharsets.US_ASCII));
            }
            assertEquals(-1, input.nextSequence(buffer));
        }
    }
}
