#!/usr/bin/env python3
import sys
import time
import tracemalloc
import pandas as pd

EGFP_ANCHOR = "ACAGTCGCAATTCCGGTTAA"
ANCHOR_LEN = len(EGFP_ANCHOR)


def find_barcode_offset(read2: str, max_search: int = 40) -> int:
    idx = read2.find(EGFP_ANCHOR, 0, max_search + ANCHOR_LEN)
    if idx != -1:
        return idx + ANCHOR_LEN

    limit = min(len(read2) - 25, max_search)
    for i in range(limit):
        window = read2[i : i + ANCHOR_LEN]
        if sum(c1 != c2 for c1, c2 in zip(window, EGFP_ANCHOR)) <= 1:
            return i + ANCHOR_LEN
    return -1


def match_library(barcode: str, ref_barcodes: list[tuple[str, str]], max_dist: int = 2) -> str | None:
    best_name = None
    best_dist = max_dist + 1
    tie = False
    for name, ref_bc in ref_barcodes:
        if barcode == ref_bc:
            return name
        dist = sum(c1 != c2 for c1, c2 in zip(barcode, ref_bc))
        if dist < best_dist:
            best_dist = dist
            best_name = name
            tie = False
        elif dist == best_dist:
            tie = True
    return None if tie else best_name


def run_baseline(lib_tsv: str, r1_path: str, r2_path: str, out_tsv: str):
    t0 = time.perf_counter()
    tracemalloc.start()

    lib_df = pd.read_csv(lib_tsv, sep="\t")
    ref_barcodes = list(zip(lib_df["variant_name"], lib_df["barcode_25bp"]))

    seen_molecules: set[tuple[str, str, str]] = set()
    records: list[tuple[str, str]] = []

    total_pairs = 0
    flank_missed = 0
    unmatched = 0
    pcr_dups = 0

    with open(r1_path, "r") as f1, open(r2_path, "r") as f2:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            seq1 = f1.readline().rstrip("\n")
            f1.readline()
            f1.readline()

            f2.readline()
            seq2 = f2.readline().rstrip("\n")
            f2.readline()
            f2.readline()

            total_pairs += 1
            cb = seq1[:16]
            umi = seq1[16:28]
            if "N" in cb or "N" in umi:
                continue

            bc_start = find_barcode_offset(seq2)
            if bc_start == -1 or bc_start + 25 > len(seq2):
                flank_missed += 1
                continue

            read_bc = seq2[bc_start : bc_start + 25]
            if "N" in read_bc:
                unmatched += 1
                continue

            variant = match_library(read_bc, ref_barcodes, max_dist=2)
            if variant is None:
                unmatched += 1
                continue

            mol_key = (cb, umi, variant)
            if mol_key not in seen_molecules:
                seen_molecules.add(mol_key)
                records.append((cb, variant))
            else:
                pcr_dups += 1

    df = pd.DataFrame(records, columns=["cell_barcode", "variant_name"])
    counts_df = df.groupby(["cell_barcode", "variant_name"], as_index=False).size()
    counts_df.rename(columns={"size": "umi_count"}, inplace=True)
    counts_df.to_csv(out_tsv, sep="\t", index=False)

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed = time.perf_counter() - t0

    print(
        f"Total Pairs: {total_pairs:,} | Unique UMIs: {len(records):,} | "
        f"PCR Dups: {pcr_dups:,} | No Flank: {flank_missed:,} | Unmatched: {unmatched:,}"
    )
    print(f"Python/Pandas Elapsed: {elapsed:.2f} sec | Peak Python Traced RAM: {peak_bytes / (1024**2):.1f} MB")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python3 scripts/baseline_demux.py <library_barcodes.tsv> <R1.fastq> <R2.fastq> <out.tsv>")
        sys.exit(1)
    run_baseline(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
