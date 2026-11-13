# jseq-zero

Zero-allocation, cache-friendly JVM primitives and CLI tools for high-throughput Next-Generation Sequencing (NGS) analysis.

`jseq-zero` applies mechanical-sympathy techniques—2-bit register encodings, `POPCNT` Hamming distance, branchless Wu-Manber/Bitap Levenshtein automata, Compressed Sparse Row (CSR) seed-and-extend indices, and open-addressed primitive hash tables (`Agrona`)—to reduce allocation pressure in sequencing workloads. Setup, dataset generation, collection growth, and output still allocate; zero allocation has not been established by profiling.

---

## Included Subcommands

* **`amplicon-demux`** — Indel-tolerant paired-end `FASTQ` barcode demultiplexer and 10x single-cell UMI deduplicator. Drop-in compatible with **scAAVengr Stage 5** (`scCAG-eGFP-Barcode-bghPolyA`), CRISPR pooled screens, and combinatorial capsid libraries up to **1,048,575 variants**.
* **`generate-benchmark`** — Synthetic paired-end `FASTQ` dataset generator injecting realistic PCR duplication depth, variable 5' frame-shifts, upstream anchor deletions, and barcode point mutations.

---

## Case Study Benchmark: Single-Cell AAV Barcoded Amplicon Screen (`amplicon-demux`)

Measured locally on macOS 14.4 (Apple Silicon), OpenJDK 25.0.1, with a 512 MB maximum heap on 2026-11-13. These are single-run results, not cross-machine guarantees.

| Measurement | Result |
| :--- | :--- |
| Dataset | 5,000,000 read pairs; 1,000,000 candidate molecules; mean PCR depth 5; 50 variants |
| Demultiplexer-reported elapsed time | 2.32 seconds |
| Whole Java process wall time | 2.459 seconds |
| Throughput | 2,154,031 read pairs/second |
| Heap use reported at completion | 44 MB (not peak RSS or a measured allocation rate) |
| Unique UMIs counted | 993,032 |
| PCR duplicates | 3,981,341 |
| Anchor not found | 0 |
| Unmatched barcode | 25,627 |

The generator samples candidate molecules with replacement, so the configured molecule count is not the number guaranteed to appear in the reads. The output TSV's UMI sum was checked against telemetry. These checks do not establish complete recovery of noisy reads or biological validation. Python speedups, 50,000-variant timing, and zero allocation have not been measured in this run.

---

## Core Architecture (`com.techtune.jseqzero.core`)

1. **`Dna2Bit`:** Branchless ASCII-to-2-bit DNA encoding `(b >> 1) & 3L` and 4-instruction `POPCNT` Hamming distance (`(xor | (xor >>> 1)) & 0x5555555555555555L`).
2. **`BitapFlankFinder`:** Cascades a 6-cycle 2-bit rolling window into a branchless Wu-Manber / Bitap Levenshtein (`k <= 2`) automaton, resolving exact downstream barcode offsets even when PCR amplification introduces indels inside the constant primer/flank.
3. **`PigeonholeBarcodeIndex`:** Flat CSR 3-slice seed-and-extend index (~2.7 MB L3-resident footprint) partitioning 50-bit barcodes into `18-bit + 16-bit + 16-bit` slices so any read with `<= 2` mismatches probes an average of ~3.4 candidates instead of scanning `100,000+` variants.
4. **`FastqRecordReader` & Primitive UMI Deduplication:** Streams uncompressed or `.gz` `FASTQ` records without allocating `String` objects, packing `[32-bit Cell Barcode | 20-bit Variant ID]` and 64-bit molecule signatures into flat `Agrona` primitive tables (`LongHashSet`, `Long2LongCounterMap`).

---

## Quick Start

Requires **JDK 25+** and **Maven 3.9+**:

```bash
# 1. Build zero-dependency executable JAR and run unit tests
mvn clean package

# 2. Generate 5,000,000 synthetic paired-end reads (1M unique molecules, 5x PCR depth, 50 variants)
java -jar target/jseq-zero.jar generate-benchmark 1000000 5 50

# 3. Run jseq-zero amplicon-demux
java -Xmx512m -jar target/jseq-zero.jar amplicon-demux library_barcodes.tsv R1.fastq R2.fastq jseq_counts.tsv

# 4. (Optional) Compare against baseline Python/Pandas implementation
python3 scripts/baseline_demux.py library_barcodes.tsv R1.fastq R2.fastq python_counts.tsv

# 5. (Optional) Test combinatorial library scaling (50,000 barcoded variants!)
java -jar target/jseq-zero.jar generate-benchmark 1000000 5 50000
java -Xmx512m -jar target/jseq-zero.jar amplicon-demux library_barcodes.tsv R1.fastq R2.fastq jseq_50k_counts.tsv
```

### Optional GraalVM Native Binary & Nextflow Integration

```bash
# Compile standalone native executable (when GraalVM JDK 25 is active)
mvn -Pnative package

# Run Nextflow DSL2 pipeline (with optional Scanpy .h5ad matrix injection)
nextflow run workflow/main.nf --reads 'data/*_R{1,2}.fastq.gz' --library_tsv data/library_barcodes.tsv
```

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
