# Native benchmark — 2026-11-13

This benchmark compares the repository's readable Python reference with its GraalVM 25 AOT executable. See [results.csv](results.csv) for every measured run and [results.json](results.json) for complete provenance, dataset hashes, source hashes, counts and hardware.

## Hardware and software

| Property | Observed value |
| --- | --- |
| Host role | Performance lab |
| CPU | 12th Gen Intel Core i7-12700 |
| Visible logical CPUs | 12; one socket, 12 reported cores, one thread per core |
| CPU availability to the session | 0, 1, 10, 11 |
| Benchmark child CPU | 0 for every implementation |
| Reported frequency range | 800–4900 MHz; governor and frequency settings unchanged |
| RAM | 30.787 GiB visible to Linux |
| L3 cache | 25 MiB |
| OS | Red Hat Enterprise Linux 9.6, x86-64 |
| Kernel | 5.14.0-570.12.1.el9_6.x86_64 |
| Virtualization detection | `systemd-detect-virt`: none |
| Python | CPython 3.11.16 |
| Native Image | Oracle GraalVM 25.0.3+9.1, serial GC, compressed references |
| Build tools | Maven 3.6.3; GCC 11.5.0 |

The build and measurements reused existing tools in a temporary directory. No packages were installed and no system settings were changed. Native Image was built with `--no-fallback -O3 -march=compatibility -J-Xmx8g --parallelism=4`. Compilation time is excluded from runtime measurements. The application requires no JVM at runtime.

## Method

1. Generate the synthetic dataset with the measured native executable: `generate-benchmark 1000000 5 50`, fixed seed 42. This produces 5,000,000 paired reads and a 50-variant barcode library.
2. Hash all input files before timing. This warms filesystem caches; caches are not flushed between runs.
3. Run each implementation three times as a fresh process, rotating which starts first. Pin only each new benchmark child to CPU 0. Use `-Xmx512m` for native Java; Python receives no artificial memory cap.
4. Measure full process elapsed time, including startup, imports and I/O. Use Linux `wait4` peak RSS, converting KiB to MiB. Neither allocation tracing nor a profiler is enabled.
5. Defer loading output tables until **after all timed runs**. This avoids inflating a forked child's initial RSS with large tables retained by the measurement parent. Check repeat output SHA-256 equality, TSV totals against telemetry, and differences between implementations.
6. Report medians and min–max ranges without selecting a favorable run. Both implementations are repository-authored, and their anchor-matching semantics differ.

The native executable SHA-256 is `8d2158e7905301cfa02e75bab955dc717e15172766c9b4256d1a474971db20c1`. Source and input hashes are in the JSON. A Java/native count difference is not evidence of better biological accuracy.

## Reproduce

From a source checkout, place the extracted native release directory alongside it as `../jseq-zero-native`. Use a fresh output directory for each measurement. The bundled `build-info.json` identifies the release binary; a different build needs matching build metadata.

```bash
mkdir benchmark-data
cd benchmark-data
../../jseq-zero-native/jseq-zero -Xmx512m generate-benchmark 1000000 5 50
cd ..
python3 scripts/benchmark_native.py \
  --binary ../jseq-zero-native/jseq-zero \
  --build-info ../jseq-zero-native/build-info.json \
  --data benchmark-data --out benchmark-results --read-pairs 5000000 --repeats 3

# Plotting tools are optional and separate from the measurement environment.
python3 -m venv .venv-benchmark
.venv-benchmark/bin/python -m pip install -r scripts/benchmark-requirements.txt
.venv-benchmark/bin/python scripts/plot_benchmark.py benchmark-results/results.json
```

PNG, editable SVG and PDF share the same data. Results depend on input composition, library size, CPU, runtime and cache state; rerun on representative workloads before making capacity decisions.
