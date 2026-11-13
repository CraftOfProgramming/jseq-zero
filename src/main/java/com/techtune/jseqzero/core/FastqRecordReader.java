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

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

/**
 * Zero-allocation streaming 4-line FASTQ record reader.
 * Copies only Line 2 (sequence bases) into a caller-supplied byte buffer while skipping
 * Line 1 (@header), Line 3 (+), and Line 4 (quality) in place without allocating Strings.
 */
public final class FastqRecordReader implements AutoCloseable {

    private final InputStream in;
    private final byte[] ioBuf = new byte[64 * 1024];
    private int pos = 0;
    private int limit = 0;

    public FastqRecordReader(final InputStream in) {
        this.in = in;
    }

    public static FastqRecordReader open(final String filePath) throws IOException {
        final InputStream raw = new BufferedInputStream(new FileInputStream(filePath), 64 * 1024);
        final InputStream stream = filePath.endsWith(".gz") ? new GZIPInputStream(raw, 64 * 1024) : raw;
        return new FastqRecordReader(stream);
    }

    /**
     * Reads the next 4-line FASTQ record, writing the sequence line into {@code targetBuffer}.
     *
     * @return sequence length in bytes, or {@code -1} on EOF
     */
    public int nextSequence(final byte[] targetBuffer) throws IOException {
        if (!skipLine()) {
            return -1; // Line 1: @Header
        }
        final int seqLen = readLineInto(targetBuffer); // Line 2: Sequence
        if (seqLen < 0) {
            throw new EOFException("FASTQ record is missing its sequence line");
        }
        if (!skipLine()) {
            throw new EOFException("FASTQ record is missing its separator line");
        }
        if (!skipLine()) {
            throw new EOFException("FASTQ record is missing its quality line");
        }
        return seqLen;
    }

    private boolean skipLine() throws IOException {
        boolean sawAnyByte = false;
        while (true) {
            if (pos >= limit && !fill()) {
                return sawAnyByte;
            }
            sawAnyByte = true;
            if (ioBuf[pos++] == '\n') {
                return true;
            }
        }
    }

    private int readLineInto(final byte[] target) throws IOException {
        int outIdx = 0;
        while (true) {
            if (pos >= limit && !fill()) {
                return outIdx == 0 ? -1 : outIdx;
            }
            final byte b = ioBuf[pos++];
            if (b == '\n') {
                return outIdx;
            }
            if (b != '\r') {
                if (outIdx == target.length) {
                    throw new IOException("FASTQ sequence exceeds buffer capacity of " + target.length + " bases");
                }
                target[outIdx++] = b;
            }
        }
    }

    private boolean fill() throws IOException {
        pos = 0;
        limit = in.read(ioBuf, 0, ioBuf.length);
        return limit > 0;
    }

    @Override
    public void close() throws IOException {
        in.close();
    }
}
