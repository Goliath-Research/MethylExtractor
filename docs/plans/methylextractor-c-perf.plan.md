---
name: MethylExtractor C Perf
overview: Cut MethylExtractor wall time from ~12 minutes toward ≤6 minutes on the same production config (threads=10, read-level on) by instrumenting phases, fixing per-region BAM overhead, then adding memory-gated multi-chromosome parallelism—without changing extraction science or HDF5/QC contracts.
todos:
  - id: instrument-timing
    content: Add elapsed-ms phase timers + sample timing.json; capture before/after on a real ~12 min BAM
    status: completed
  - id: intra-chrom-io
    content: Share bam_hdr across region threads; hts_set_threads on per-region samFile; verify HDF5 Zstd plugin path
    status: completed
  - id: chrom-pool
    content: Implement --chrom-parallel + --max-rss-gb memory-gated chrom pool with largest-first schedule in main.c / bam_processing.c
    status: completed
  - id: pipeline-config
    content: Wire chrom_parallel/max_rss_gb through extract_runner + step schema/profile (config-not-code)
    status: completed
  - id: parity-docs
    content: Integration identity test chrom_parallel 1 vs 2; golden H5/manifest parity; update ANALYSIS/README and docs/plans
    status: completed
---

> **Status: IMPLEMENTED.** CLI `--chrom-parallel` / `--max-rss-gb`, `{sample}.timing.json`, HDF5 Zstd probe, and MethylPipeline `actionConfig.methyl_extract` knobs. Production ≤6 min gate is measured from `{sample}.timing.json` on a representative WGBS BAM (not a checked-in benchmark).

# MethylExtractor C performance (≤6 min gate)

## Context

Chromosomes previously ran **serially** in `src/main.c` while `-t/--threads` only split sites within one chrom. Production profiles pass `"threads": 10` and often enable `read_level`. Dense per-region `mC`/`uC` arrays are ~`8 × genomic_span` bytes, so multi-chrom is memory-gated.

**Success gate:** on one representative WGBS BAM that currently ~12 min with the production flags, wall ≤ **6 min**, bit-identical `{chrom}-{ctx}.h5` + manifest vs baseline. Toy identity: `tests/integration/test_bam_to_manifest.py::test_chrom_parallel_identity`.

## What landed

- Elapsed-ms phase logs + `{sample}.timing.json` (`methylextractor.timing`).
- Shared BAM index; skip per-region `sam_hdr_read`; `hts_set_threads` on each region `samFile*`; startup Zstd/`HDF5_PLUGIN_PATH` probe; HDF5/cJSON writes serialized with a recursive mutex.
- `--chrom-parallel` (default 2, dropped to 1 if two largest chroms exceed `--max-rss-gb`, default 32) with largest-first schedule; `--threads T` split as `max(1, T/K)` region workers per in-flight chrom.
- Pipeline: `MethylExtractStepConfig`, `extract_runner.py` flags, research profiles `chrom_parallel: 2` / `max_rss_gb: 32`.
