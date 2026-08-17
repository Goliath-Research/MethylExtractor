# Extraction QC export contract (MethylExtractor → MethylPipeline)

MethylExtractor exports extraction QC metrics at extraction time. **MethylPipeline** is responsible for reading H5 files (when needed), consuming these JSON exports, and evaluating Pass/Fail guardrails.

## Files produced

| File | Schema | When |
|------|--------|------|
| `{chrom}-{context}.json` | `methylextractor.context_qc` v1.0.0 | Per H5/TXT output (e.g. `1-CG.json`) |
| `{sample_id}.extraction_manifest.json` | `methylextractor.extraction_manifest` v1.0.0 | End of run, in output directory |
| `{sample_id}.timing.json` | `methylextractor.timing` v1.0.0 | Diagnostic phase elapsed-ms (not a QC guardrail) |

Chromosome names follow `chrom_mapping.json` (`"1"`..`"22"`, `"X"`, `"Y"`). Contexts: `CG` always; `CHG`/`CHH` when `-G`/`-H` are used with `-s`.

## Reference site counts

Each context JSON and manifest entry includes `sites_in_reference` — the number of cytosine sites enumerated from the reference for that chromosome and context. MethylPipeline can compute:

```
fraction_sites_covered = sites_passing_min_cov / sites_in_reference
```

without re-scanning the reference FASTA. Genome-wide fraction is in manifest `summary.cpg_fraction_sites_covered`.

## Pass/Fail (MethylPipeline)

Suggested blocking guardrails (configurable in MethylPipeline):

| Guardrail | Source field | Default |
|-----------|--------------|---------|
| CpG mean coverage | `summary.cpg_weighted_mean_coverage` | ≥ 10 |
| CHH methylation | `summary.chh_methylation_level` | ≤ 0.02 (if CHH extracted) |
| CHG methylation | `summary.chg_methylation_level` | ≤ 0.02 (if CHG extracted) |
| File completeness | `per_chromosome` keys vs `chrom_mapping.json` | all present |
| Chromosome uniformity | `per_chromosome.*.CG.mean_coverage` | min/median autosomal ≥ 0.5 |

CpG methylation level is **informational only** (tissue-dependent).

## Legacy compatibility

Per-file JSON retains top-level fields: `num_positions`, `total_methylated`, `total_unmethylated`, `avg_methylation_level`, `avg_coverage`.
