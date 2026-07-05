# MethylExtractor vs. MethylDackel — Software Engineering Analysis

**Date:** 2026-06-11
**Scope:** Source-level review of MethylExtractor (`src/`, `include/`, build system, scripts, docs) and comparison against MethylDackel as the reference WGBS methylation extractor.
**Goal:** Confirm whether MethylExtractor performs a *similar job* to MethylDackel, and whether it does so with *better performance* and *better organization* from a Software Engineering standpoint.

---

## 1. Executive Summary

| Question | Verdict |
|----------|---------|
| Does it do a **similar job** to MethylDackel? | **Yes, partially.** The core extraction (per-cytosine mC/uC counts by context and strand, with MAPQ/Phred filtering and Bismark `XG` strand detection) is functionally comparable for the common case. Several MethylDackel features are missing, and one headline feature (overlapping-mate de-duplication) is *claimed but not implemented*. |
| Is the **output organization better**? | **Yes.** Per-chromosome / per-context HDF5 with a compact 12-byte record and Zstd compression is a genuinely better storage design than MethylDackel's text bedGraph/cytosine-report output, for large cohorts. |
| Is the **performance better**? | **Unproven and partly undermined by bugs.** The intra-chromosome region parallelism is a sound idea, but the chromosome-level thread scheduler is effectively serial and contains undefined behavior, the `--threads` flag has almost no effect, and the documented "`-z 9` = 5× faster" optimization **does not exist in the code**. |
| Is the **engineering quality better**? | **Mixed.** Code is reasonably modular and readable, but there are correctness bugs, doc/code drift, broken helper scripts, dead tooling, no CI, no tests of substance, and inconsistent licensing. |

**Bottom line:** The *concept* and the *storage/organization layer* are a real improvement over MethylDackel. The *performance* and *correctness* claims in the documentation are **not currently substantiated by the code** and need to be either implemented or removed before the "better performance" claim can be confirmed.

---

## 2. What the Tool Does (as implemented)

MethylExtractor is a multi-file C program built on HTSlib + HDF5:

- **`src/main.c`** — CLI parsing, chromosome-mapping load, reference + BAM header load, thread orchestration.
- **`src/bam_processing.c`** — per-region BAM iteration, strand inference, mC/uC counting, chromosome driver, JSON mapping loader.
- **`src/output_formats.c`** — coverage filtering/capping, TXT writer, HDF5 writer (compound type + Zstd/gzip), statistics JSON.
- **`src/utils.c`** — context detection (CpG/CHG/CHH), trinucleotide encode/decode, site enumeration, buffer init.
- **`include/methyl_extractor.h`** — shared structs and constants.

Pipeline per chromosome:
1. Enumerate all candidate cytosine sites from the reference and pre-fill a dense `MethylRecord` buffer (position, context, strand, trinucleotide).
2. Split the site list into regions; each region thread re-opens the BAM, queries its coordinate window, and tallies methylated/unmethylated calls into private dense arrays.
3. Merge private counts back into the per-site buffer.
4. Filter by minimum coverage (optionally cap coverage), then write HDF5 and/or TXT plus a statistics JSON.

The compact on-disk record is a real strength:

```53:67:include/methyl_extractor.h
typedef struct
{
    unsigned tnc : 5;
    unsigned context : 2;
    unsigned strand : 1;
} tnc_bitfield_t;

typedef struct
{
    uint32_t pos;
    uint16_t mC;
    uint16_t uC;
    tnc_bitfield_t tnc;
    uint8_t _pad[1];
} MethylRecord;
```

12 bytes/site, columnar-compressible with Zstd, partitioned `{chr}.h5` (or `{chr}-{ctx}.h5` with `--split`). This is the part that delivers on the "better organization / highly-compressed HDF5" promise.

---

## 3. Functional Comparison with MethylDackel

| Capability | MethylDackel | MethylExtractor | Notes |
|------------|--------------|-----------------|-------|
| Input | Sorted+indexed BAM/CRAM via HTSlib | Sorted+indexed BAM via HTSlib | Equivalent backend. |
| Contexts | CpG, CHG, CHH | CpG (always), CHG (`-G`), CHH (`-H`) | Equivalent. |
| Per-strand reporting | Yes (+ `--mergeContext`) | Yes (always per-strand; no merge option) | ME records both strands of a CpG as separate sites; no merge flag. |
| Strand inference | OT/OB/CTOT/CTOB model | Bismark `XG` tag + flag heuristic fallback | See `getRealStrand`, `src/bam_processing.c`. Reasonable. |
| Read filters | MAPQ, dup, secondary, supplementary | MAPQ, dup, secondary, supplementary, QC-fail, `NH>1` | Comparable; ME also drops multimappers via `NH`. |
| Base filter | Phred threshold | Phred threshold (`-p`) | Equivalent. |
| **Overlapping mate de-dup** | **Yes** | **Claimed, but NOT implemented** | See §5.1 — major gap. |
| M-bias analysis (`mbias`) | Yes (plots + trimming) | No | MethylDackel feature not ported. |
| Read-level / per-read methylation | Yes (`--perRead`) | Yes (`--read-level` tile-pattern sidecars) | Sidecar `{chrom}-{ctx}.patterns.h5` per MethylPipeline contract; not per-read text dumps. |
| Output formats | bedGraph, methylKit, cytosine report | HDF5, TXT, both; per-context split | ME trades interchange formats for compact HDF5. |
| Coverage cap | No (external) | Yes (`-C`, mean-based downscale) | ME adds this natively. |
| Stats summary | No native JSON | Per-file stats JSON | ME adds this. |

**Conclusion for §3:** For the common "extract mC/uC per cytosine by context" task, MethylExtractor is functionally comparable to MethylDackel and adds useful native features (coverage capping, stats JSON, compact partitioned output, read-level pattern sidecars). It omits M-bias analysis and per-read text output, and — critically — does not implement the overlapping-mate handling it advertises.

---

## 4. Organization / Storage — where MethylExtractor is genuinely better

1. **Compact binary record** (12 bytes vs. multi-column ASCII text). For genome-wide CpG+CHG+CHH this is a large reduction before compression even applies.
2. **Per-chromosome / per-context partitioning** (`{chr}.h5`, `{chr}-{ctx}.h5`) enables selective, parallel downstream reads — much better than one large text file.
3. **Chunked + Zstd-compressed HDF5** is a standard, tool-agnostic scientific container; downstream consumers can memory-map and slice.
4. **Self-describing compound datatype** (`pos`, `mC`, `uC`, `tnc`) written via `H5Tcreate`/`H5Tinsert` in `flush_buffer` — good schema hygiene.

This storage/organization layer is the strongest justification for the project and is a real improvement over MethylDackel's text outputs for large-scale storage-critical cohorts.

---

## 5. Correctness & Performance Findings (Software Engineering)

### 5.1 Headline feature "overlapping read pair handling" is not implemented — **High severity**

README:

> "Proper handling of overlapping paired-end reads to prevent double-counting"

The counting loop in `process_region_direct` iterates each read independently and unconditionally increments counters; there is **no mate/overlap detection** anywhere in `src/`:

```113:132:src/bam_processing.c
                    int base = bam_seqi(seq, qpos);
                    uint32_t rel = (uint32_t)(pos - t->start_pos);

                    if (ref == 'C' && (strand == 1 || strand == 3))
                    {
                        if (base == 2)
                            r->counts.mC[rel]++; // C->C methylated
                        else if (base == 8)
                            r->counts.uC[rel]++; // C->T unmethylated
                    }
                    else if (ref == 'G' && (strand == 2 || strand == 4))
                    {
                        if (base == 4)
                            r->counts.mC[rel]++; // G->G methylated
                        else if (base == 1)
                            r->counts.uC[rel]++; // G->A unmethylated
                    }
```

A grep across `src/` for overlap/mate handling returns nothing. For paired-end WGBS with overlapping fragments, this **double-counts** the overlap, biasing per-site coverage and methylation level — the exact problem MethylDackel solves and that the README claims to solve. Either implement it or remove the claim.

### 5.2 Chromosome-level threading is effectively serial + undefined behavior — **High severity**

The orchestration loop creates one thread, then immediately joins *all* previously created threads on each iteration:

```419:443:src/main.c
        if (pthread_create(&threads[i], NULL, (void *(*)(void *))process_chromosome,
                           thread_args_copies[i]) != 0)
        {
            ...
        }
        active_threads++;
        for (int j = 0; j <= i; j++)
            if (pthread_join(threads[j], NULL) == 0)
                active_threads--;

        while (active_threads >= num_threads)
            for (int j = 0; j <= i; j++)
                if (pthread_join(threads[j], NULL) == 0)
                    active_threads--;
    }

    log_time("Joining threads...\n");
    for (int i = 0; i < valid_chr_count; i++)
        if (pthread_join(threads[i], NULL) == 0)
```

Problems:
- The inner `for (j = 0; j <= i; ...) pthread_join(threads[j], ...)` **re-joins already-joined threads** every iteration. Joining a thread more than once is **undefined behavior** (POSIX). The final join loop joins them yet again.
- Because each iteration joins thread `i` right after creating it, chromosomes are processed **one at a time** — the chromosome-level parallelism the design implies does not happen.
- `num_threads` (the user's `-t/--threads`) only gates a `while` loop whose condition is already false (active count is driven back to 0 by the join loop), so **`--threads` has essentially no effect** at the chromosome level.

The actual parallelism comes only from `process_chromosome`, which spawns its own region threads using `get_nprocs()` (capped at 32), independent of `--threads`:

```164:170:src/bam_processing.c
    int n_threads = get_nprocs();
    if (n_threads > 32)
        n_threads = 32;
    if (n_threads < 1)
        n_threads = 8;
```

Net effect: per-chromosome region parallelism works, but the documented top-level threading model is broken, contains UB, and ignores the user's thread setting.

### 5.3 Documented "`-z 9` → ~5× faster" optimization does not exist — **High severity (doc/code drift)**

README "Performance Optimization (New!)" promises that `-z 9` writes uncompressed then launches background multi-threaded `zstd` for ~5× speedup. `SUGGESTIONS.md` contains the proposed patch. **Neither is in the code.** `flush_buffer` only ever sets the HDF5 Zstd filter (levels 1–8) or gzip fallback; there is no level-9 branch, no `system("zstd ...")`, and `main.c`'s help text advertises only `0=none, 1-8=Zstd`:

```254:280:src/output_formats.c
        // Set compression (Zstd with gzip fallback)
        if (compression > 0)
        {
            // Levels 1–8: try Zstd first, fall back to gzip if not available
            unsigned int cd_values[1] = {(unsigned int)compression};
            status = H5Pset_filter(dcpl, ZSTD_FILTER, H5Z_FLAG_OPTIONAL, 1, cd_values);
            ...
        }
        else
        {
            // compression == 0 → truly raw
            log_time("Writing raw uncompressed HDF5\n");
        }
```

Passing `-z 9` just sets `cd_values = {9}`. The performance claim is unsubstantiated.

### 5.4 Silent "no compression" risk from `H5Z_FLAG_OPTIONAL` — **Medium severity**

The gzip-fallback logic keys off the return value of `H5Pset_filter`, but that call rarely fails for a known filter ID — actual unavailability of the Zstd plugin surfaces later at `H5Dwrite`. Because the filter is set `H5Z_FLAG_OPTIONAL`, a missing plugin causes HDF5 to **silently write uncompressed** rather than error or fall back to gzip. So if `HDF5_PLUGIN_PATH` is not configured, output may be uncompressed without any warning, and the gzip fallback never triggers. The fallback should verify filter availability (e.g. `H5Zfilter_avail(ZSTD_FILTER)`) before relying on it.

### 5.5 Repeated per-region index load — **Medium severity (performance)**

Every region thread re-opens the BAM and reloads the index:

```42:51:src/bam_processing.c
    samFile *in = sam_open(t->bam_file, "r");
    if (!in)
        return NULL;
    hts_idx_load = sam_index_load(in, t->bam_file);
```
(`sam_open` + `sam_index_load` + `sam_hdr_read` per region). With up to 32 regions per chromosome × every chromosome, this is substantial redundant I/O and parsing versus loading the index once and sharing it. Also, the pipeline assumes a coordinate-sorted, indexed BAM but `main.c` performs no check; indexing is left to the separate `sort_and_index.sh`.

### 5.6 Memory footprint — **Low/Medium severity**

Per chromosome, each region allocates two dense `uint32_t` arrays sized to its base-pair span (`r->counts.mC`/`uC`), plus the full-chromosome `MethylRecord` buffer and a `site_positions` array. For the largest human chromosomes this is on the order of ~2 GB of transient counting memory. It is bounded by chromosome size (as the README states) only because §5.2's bug serializes chromosomes; if the threading were "fixed" to run chromosomes concurrently, peak memory would multiply and could OOM. Counts could use `uint16_t` (matching the `uint16_t` output fields) or a sparse structure.

### 5.7 Region span `+10` heuristic — **Low severity (not a correctness bug)**

```186:189:src/bam_processing.c
        r->base.start_pos = site_positions[r->start_site_idx];
        r->base.end_pos = (r->end_site_idx < site_count)
                              ? site_positions[r->end_site_idx - 1] + 10
                              : targ->chr_len;
```
Adjacent regions can overlap by up to 10 bp, but per-position guards (`pos < start_pos || pos >= end_pos`) and the index-scoped merge ensure each site is attributed once, so output is correct. The overlap only wastes a little counting work. Worth a comment explaining why `+10` is safe.

---

## 6. Build, Tests, Tooling, Docs (Engineering Hygiene)

### 6.1 Broken helper script — **High severity**
`process_all_bams.sh` invokes flags that **do not exist** in the CLI and passes positional args in the wrong order:

```23:38:process_all_bams.sh
    /home/ubuntu/MethylExtractor/build/dynamic/MethylExtractor \
        --chrom-mapping "$CHROM_MAPPING" \
        --@ 16 \
        --o "$output_dir" \
        --hdf5-compression 9 \
        ...
        --no-cap-coverage \
        --split-context-files \
        "$REF_GENOME" \
        "$bam_file"
```
The real options are `-t/--threads`, `-o/--output-dir`, `-z/--compression`, `-c/--min-cov`, `-s/--split` (no `--@`, `--hdf5-compression`, `--no-cap-coverage`, `--split-context-files`), and the positional order is `<input.bam> <output_dir> [ref.fa]` — here `ref` and `bam` are swapped and `out_dir` is absent. This script cannot run successfully; it reflects an older CLI and was never updated. Doc/code drift in shipped tooling.

### 6.2 Stale / irrelevant tooling — **Medium severity**
- `tests/test_parquet.py` converts CSV→Parquet, but the tool emits TSV/`.txt` (not CSV) and Parquet support was explicitly removed (`output_formats.c` comment "Parquet output removed"). The test is dead.
- `src/cap_coverage.py` duplicates the C coverage-capping logic in Python — redundant, and not wired into anything.
- `pyproject.toml` + `poetry.lock` add Python packaging metadata to a C project with no Python package.

### 6.3 Tests — **Medium severity**
The only C "test" (`tests/test_statistics_example.c`) **copies** the struct and stats functions rather than linking against the real source, so it cannot catch regressions in the shipped code. There is no integration test against a known BAM, and no CI configuration anywhere in the repo.

### 6.4 Build system — **Low/Medium severity**
- `make all` runs `sudo apt-get update/install` as a side effect of building (`deps` is a prerequisite of `all`). Auto-`sudo` during a build is surprising and CI-hostile; deps should be opt-in.
- Hardcoded Debian/Ubuntu paths and `<sys/sysinfo.h>`/`get_nprocs` make it Linux-only (acceptable if that's the only target, but it limits portability and the README's general tone).

### 6.5 Documentation/licensing inconsistencies — **Low severity**
- License conflict: `README.md` says **MIT**, `pyproject.toml` says **"Copyright by Epimethyl Analytics"**, and there is **no `LICENSE` file**.
- README advertises features not present (overlap handling §5.1, `-z 9` optimization §5.3) and a `-z` range (`0–9`) that disagrees with the help text (`0–8`).
- Default `--threads` documented as 16, but the value is overridden to the host CPU count when left at default — and barely used (§5.2).

---

## 7. Recommendations (priority order)

1. **Implement overlapping-mate de-duplication** (or remove the claim). This is the single most important correctness item versus MethylDackel.
2. **Fix the chromosome-level thread scheduler** in `main.c`: remove the repeated `pthread_join` (UB), drive a real worker pool bounded by `--threads`, and make `--threads` actually govern parallelism. Then re-measure memory (§5.6).
3. **Reconcile compression docs with code**: either implement the `-z 9` external multi-threaded Zstd path or delete the "Performance Optimization (New!)" section and `SUGGESTIONS.md`.
4. **Make compression fallback robust**: check `H5Zfilter_avail(ZSTD_FILTER)` and warn/fall back to gzip explicitly instead of relying on `H5Pset_filter` return + `OPTIONAL`.
5. **Load the BAM index once** and share across region threads.
6. **Repair or delete shipped scripts/tests**: fix `process_all_bams.sh` to the real CLI; remove stale `test_parquet.py`, `cap_coverage.py` duplication, and unused Python packaging or document why they exist.
7. **Add a real test + CI**: a small fixture BAM + reference with expected mC/uC, linked against the actual sources, run in CI.
8. **Resolve licensing**: add a `LICENSE` file matching the intended terms and make `README`/`pyproject` agree.
9. **Benchmark to substantiate "better performance"**: publish wall-clock and output-size numbers vs. MethylDackel on the same dataset/hardware before asserting the performance advantage.

---

## 8. Final Assessment

- **Similar job:** Confirmed for the core extraction task, with notable feature gaps (overlap de-dup, M-bias, per-read text dumps). Read-level co-methylation sidecars (`--read-level`) are supported.
- **Better organization:** Confirmed. The compact, chunked, Zstd-compressed, per-chromosome/per-context HDF5 schema is a real, well-motivated improvement over MethylDackel's text outputs.
- **Better performance:** Not confirmed. The intra-chromosome parallel design is promising, but the top-level threading is bugged and serial, `--threads` is inert, and the advertised fast-compression path is absent. No benchmarks exist in-repo to support the claim.
- **Better software engineering overall:** Partially. Good modular layout and a strong storage concept are offset by correctness bugs, doc/code drift, broken scripts, weak tests, no CI, and licensing inconsistencies.

**Recommendation:** Treat the storage/organization layer as production-promising, but address §5.1, §5.2, and §5.3 before claiming functional parity with MethylDackel or a performance advantage. Once mate de-duplication is implemented, the threading is corrected, and benchmarks are published, MethylExtractor would have a credible case as a better-organized, competitive alternative.
