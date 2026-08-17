# MethylExtractor vs. MethylDackel — Software Engineering Analysis

**Date:** 2026-06-11 (refreshed 2026-08-17)
**Scope:** Source-level review of MethylExtractor (`src/`, `include/`, build system, scripts, docs) and comparison against MethylDackel as the reference WGBS methylation extractor.
**Goal:** Confirm whether MethylExtractor performs a *similar job* to MethylDackel, and whether it does so with *better performance* and *better organization* from a Software Engineering standpoint.

---

## 1. Executive Summary

| Question | Verdict |
|----------|---------|
| Does it do a **similar job** to MethylDackel? | **Yes, for the common extraction path.** Per-cytosine mC/uC by context/strand, MAPQ/Phred filtering, Bismark `XG` strand detection, plus overlapping-mate coordinate clip (implemented). Several MethylDackel features remain unported (M-bias plots, bedGraph interchange). |
| Is the **output organization better**? | **Yes.** Per-chromosome / per-context HDF5 with a compact 12-byte record and Zstd compression is a genuinely better storage design than MethylDackel's text bedGraph/cytosine-report output, for large cohorts. |
| Is the **performance better**? | **Partially.** Intra-chromosome region parallelism plus memory-gated `--chrom-parallel` (default 2). Published “`-z 9` = 5× faster” claims are still unsubstantiated. Gate: ≤6 min wall on a sample that was ~12 min at `threads=10` + read-level. |
| Is the **engineering quality better**? | **Improved since the June 2026 review.** Overlap mate handling and extraction QC JSON contracts are in production; BAM→manifest integration tests exist under `tests/integration/`. Remaining gaps: chromosome parallelism, stronger automated coverage, doc hygiene. |

**Bottom line:** MethylExtractor is the MethylPipeline **linear / stock-pangenome** extraction tool. Its JSON exports feed **`sample.extraction_qc`** (`methylextractionqc`), not alignment QC (`methylalignmentqc`). For `pangenome_wgbs`, methylGrapher MethylCall replaces MethylExtractor but emits a compatible extraction manifest.

---

## 2. What the Tool Does (as implemented)

MethylExtractor is a multi-file C program built on HTSlib + HDF5:

- **`src/main.c`** — CLI parsing, chromosome-mapping load, reference + BAM header load, thread orchestration.
- **`src/bam_processing.c`** — per-region BAM iteration, strand inference, mC/uC counting, **overlapping-mate clip**, chromosome driver, JSON mapping loader.
- **`src/output_formats.c`** — coverage filtering/capping, TXT writer, HDF5 writer (compound type + Zstd/gzip), statistics JSON trigger.
- **`src/extraction_export.c`** — context QC sidecars + sample `extraction_manifest.json`.
- **`src/utils.c`** — context detection (CpG/CHG/CHH), trinucleotide encode/decode, site enumeration, buffer init.
- **`include/methyl_extractor.h`** — shared structs and constants.

Pipeline per chromosome:
1. Enumerate all candidate cytosine sites from the reference and pre-fill a dense `MethylRecord` buffer.
2. Split the site list into regions; each region thread tallies mC/uC into private dense arrays.
3. Merge private counts back into the per-site buffer.
4. Filter by minimum coverage (optionally cap coverage), then write HDF5 and/or TXT plus statistics JSON.

The compact on-disk record (12 bytes/site, Zstd HDF5, `{chr}-{ctx}.h5` with `--split`) is the strongest storage improvement over MethylDackel text outputs.

---

## 3. Functional Comparison with MethylDackel

| Capability | MethylDackel | MethylExtractor | Notes |
|------------|--------------|-----------------|-------|
| Input | Sorted+indexed BAM/CRAM via HTSlib | Sorted+indexed BAM via HTSlib | Equivalent backend. |
| Contexts | CpG, CHG, CHH | CpG (always), CHG (`-G`), CHH (`-H`) | Equivalent. |
| Per-strand reporting | Yes (+ `--mergeContext`) | Yes (always per-strand; no merge option) | No merge flag in ME. |
| Strand inference | OT/OB/CTOT/CTOB model | Bismark `XG` tag + flag heuristic fallback | See `getRealStrand`. |
| Read filters | MAPQ, dup, secondary, supplementary | MAPQ, dup, secondary, supplementary, QC-fail, `NH>1` | ME also drops multimappers via `NH`. |
| Base filter | Phred threshold | Phred threshold (`-p`) | Equivalent. |
| **Overlapping mate de-dup** | **Yes** | **Yes (coordinate clip)** | `bases_skipped_overlap_clip` in JSON. |
| M-bias analysis (`mbias`) | Yes | No | Not ported. |
| Read-level methylation | `--perRead` text | `--read-level` tile-pattern HDF5 sidecars | MethylPipeline contract. |
| Output formats | bedGraph, methylKit, cytosine report | HDF5, TXT, both; per-context split | Compact HDF5 focus. |
| Coverage cap | External | Native `-C` | ME adds this. |
| Stats summary | No native JSON | Context QC + extraction manifest | → `methylextractionqc`. |

---

## 4. Organization / Storage — where MethylExtractor is genuinely better

1. Compact 12-byte binary record vs multi-column ASCII.
2. Per-chromosome / per-context partitioning for parallel downstream reads.
3. Chunked Zstd HDF5 as a standard scientific container.
4. Extraction QC contract (`docs/extraction_qc_contract.md`) for pipeline guardrails.

---

## 5. Correctness & Performance Findings

### 5.1 Overlapping read pair handling — **Implemented (was High; now informational)**

**Status (2026-08-15):** Implemented in `process_region_direct` (`src/bam_processing.c`). For FR pairs on the same contig, the left mate yields reference positions `>= mpos` to the right mate; on equal starts, READ2 yields. Skipped bases increment `bases_skipped_overlap_clip`.

The June 2026 finding that this was missing is **obsolete**. Integration coverage: `tests/integration/test_bam_to_manifest.py`.

### 5.2 Chromosome-level threading — **Implemented (chrom-parallel pool)**

`--chrom-parallel K` (default 2, dropped to 1 when the two largest chromosomes exceed `--max-rss-gb`) runs up to K chromosomes in `process_chromosome` at once. `--threads T` is split as `max(1, T/K)` region workers per in-flight chromosome. Chromosomes are scheduled largest-first. Dense mC/uC arrays still scale with chromosome length; the RSS gate prevents stacking chr1+chr2 when the estimate exceeds the budget.

Phase elapsed-ms land in `{sample}.timing.json`. Intra-chromosome I/O: BAM index is shared; each region thread still has its own `samFile*` (htslib rule) with `hts_set_threads` for BGZF inflate. Header is not re-read per region.

### 5.3 Documented "`-z 9` → ~5× faster" optimization — **Doc/code drift**

Background multi-threaded zstd claims remain unsubstantiated; document actual Zstd filter levels only.

### 5.4 Silent uncompressed HDF5 risk — **Medium**

Missing Zstd plugin falls back to gzip (`H5Z_FLAG_MANDATORY` then gzip). Startup logs `HDF5_PLUGIN_PATH` and whether the Zstd filter is available. Set `HDF5_PLUGIN_PATH` (MethylPipeline worker install does this).

### 5.5 Per-region BAM handle — **Low (addressed)**

Each region thread still opens its own BAM handle (required for concurrent HTSlib readers). The index is loaded once per chromosome and shared read-only. `hts_set_threads` enables BGZF decompress workers on each handle.

### 5.6 Memory footprint — **Low/Medium**

Dense per-chromosome counting arrays; `--chrom-parallel` + `--max-rss-gb` bound peak RSS to in-flight chromosomes.

---

## 6. QC boundary (alignment vs extraction)

| Stage | Tooling | MethylExtractor role |
|-------|---------|----------------------|
| Alignment QC | MethylPipeline `methylalignmentqc` / `sample.methyl_qc` | **None** |
| Extraction QC | MethylPipeline `methylextractionqc` / `sample.extraction_qc` | **Producer** of manifests / context JSON |

Do not describe MethylExtractor JSON as “alignment QC exports.”

---

## 7. Testing

| Asset | Status |
|-------|--------|
| `tests/integration/test_bam_to_manifest.py` | BAM fixture → HDF5 + manifest schema assert |
| `tests/fixtures/toy_wgbs/` | Tiny FASTA/BAM/chrom_mapping |
| `tests/test_statistics_example.c` | Historical mock — prefer integration test |

```bash
python3 -m pytest tests/integration/test_bam_to_manifest.py -q
# or: METHYL_EXTRACTOR_BIN=/path/to/MethylExtractor pytest ...
```

---

## 8. MethylPipeline coupling

- Worker: `workers/methyl_worker/extract_runner.py` (`sample.methyl_extract`)
- Contract: `docs/extraction_qc_contract.md` ↔ `packages/methylextractionqc/`
- Cross-repo map: MethylPipeline `docs/architecture/sample-prep-tooling.md` (AB#703)
