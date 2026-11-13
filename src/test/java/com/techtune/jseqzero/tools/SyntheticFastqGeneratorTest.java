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

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import java.nio.file.Files;
import java.nio.file.Path;
import static org.junit.jupiter.api.Assertions.*;

class SyntheticFastqGeneratorTest {
    @TempDir Path directory;

    @Test void generatesDeterministicPairedFastqWithMatchingQualityLengths() throws Exception {
        Path first=Files.createDirectory(directory.resolve("first")),second=Files.createDirectory(directory.resolve("second"));
        SyntheticFastqGenerator.generate(first,100,5,8);
        SyntheticFastqGenerator.generate(second,100,5,8);
        assertEquals(9,Files.readAllLines(first.resolve("library_barcodes.tsv")).size());
        for(String file:new String[]{"R1.fastq","R2.fastq"}) {
            assertArrayEquals(Files.readAllBytes(first.resolve(file)),Files.readAllBytes(second.resolve(file)));
            var lines=Files.readAllLines(first.resolve(file));
            assertEquals(500*4,lines.size());
            for(int i=0;i<lines.size();i+=4) {
                assertEquals("@R"+(i/4),lines.get(i));
                assertEquals("+",lines.get(i+2));
                assertEquals(lines.get(i+1).length(),lines.get(i+3).length());
                assertTrue(lines.get(i+1).matches("[ACGT]+"));
            }
        }
    }

    @Test void rejectsInvalidSizesBeforeGeneratingFiles() {
        assertThrows(IllegalArgumentException.class,()->SyntheticFastqGenerator.generate(directory,0,5,50));
        assertThrows(IllegalArgumentException.class,()->SyntheticFastqGenerator.generate(directory,10,0,50));
        assertThrows(IllegalArgumentException.class,()->SyntheticFastqGenerator.generate(directory,10,5,0));
        assertThrows(ArithmeticException.class,()->SyntheticFastqGenerator.generate(directory,Integer.MAX_VALUE,2,50));
        assertFalse(Files.exists(directory.resolve("R1.fastq")));
    }
}
