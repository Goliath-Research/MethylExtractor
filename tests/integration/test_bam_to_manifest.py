"""Integration: BAM fixture → HDF5 + extraction_manifest.json (production binary)."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "fixtures" / "toy_wgbs"
REPO_ROOT = Path(__file__).resolve().parents[2]


def _resolve_extractor_bin() -> Path:
    env = (os.environ.get("METHYL_EXTRACTOR_BIN") or "").strip()
    if env:
        p = Path(env)
        if p.is_file():
            return p
    for candidate in (
        REPO_ROOT / "build" / "dynamic" / "arm64" / "MethylExtractor",
        REPO_ROOT / "build" / "dynamic" / "x64" / "MethylExtractor",
        Path("/work/epimethyl/methyl-extractor-aarch64/bin/MethylExtractor"),
        Path("/work/epimethyl/methyl-extractor-amd64/bin/MethylExtractor"),
    ):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "MethylExtractor binary not found; set METHYL_EXTRACTOR_BIN or build with make"
    )


@pytest.fixture
def toy_run(tmp_path: Path) -> Path:
    if not (FIXTURE_DIR / "toy.bam").is_file():
        pytest.skip(f"missing fixture BAM under {FIXTURE_DIR}")
    work = tmp_path / "sample_toy"
    work.mkdir()
    for name in ("toy.fa", "toy.fa.fai", "toy.bam", "toy.bam.bai"):
        shutil.copy2(FIXTURE_DIR / name, work / name)
    mapping = json.loads((FIXTURE_DIR / "chrom_mapping.json").read_text(encoding="utf-8"))
    mapping["reference"] = str(work / "toy.fa")
    map_path = work / "chrom_mapping.json"
    map_path.write_text(json.dumps(mapping, indent=2) + "\n", encoding="utf-8")

    out_dir = work / "sample_toy"
    out_dir.mkdir()
    bin_path = _resolve_extractor_bin()
    cmd = [
        str(bin_path),
        "-t",
        "1",
        "-q",
        "0",
        "-p",
        "0",
        "-c",
        "1",
        "-s",
        "-f",
        "hdf5",
        "-m",
        str(map_path),
        "-o",
        str(out_dir),
        str(work / "toy.bam"),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise AssertionError(
            f"MethylExtractor failed ({proc.returncode}):\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return out_dir


def test_bam_to_manifest_schema(toy_run: Path) -> None:
    manifest_path = toy_run / "sample_toy.extraction_manifest.json"
    assert manifest_path.is_file(), f"missing manifest under {toy_run}"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    meta = manifest["metadata"]
    assert meta["schema_name"] == "methylextractor.extraction_manifest"
    assert meta["schema_version"] == "1.0.0"
    assert meta["sample_id"] == "sample_toy"
    assert "CG" in meta["contexts_extracted"]
    assert meta["filters"]["split_context_files"] in (1, True)

    summary = manifest["summary"]
    assert summary["chromosomes_processed"] >= 1
    assert summary["cpg_sites_passing_min_cov"] >= 1
    assert summary["cpg_sites_in_reference"] >= summary["cpg_sites_passing_min_cov"]
    assert summary["cpg_weighted_mean_coverage"] > 0

    assert "read_filtering" in manifest
    rf = manifest["read_filtering"]
    assert rf["reads_seen"] >= 1
    assert rf["reads_used"] >= 1
    # Overlap clip counter exists (may be zero for unpaired toy reads)
    assert "bases_skipped_overlap_clip" in rf

    per_chr = manifest["per_chromosome"]
    assert "1" in per_chr
    assert "CG" in per_chr["1"]
    assert per_chr["1"]["CG"]["num_positions"] >= 1

    h5 = toy_run / "1-CG.h5"
    ctx_json = toy_run / "1-CG.json"
    assert h5.is_file() and h5.stat().st_size > 0
    assert ctx_json.is_file()
    ctx = json.loads(ctx_json.read_text(encoding="utf-8"))
    assert ctx["metadata"]["schema_name"] == "methylextractor.context_qc"
    assert ctx["sites"]["sites_in_reference"] >= 1
    assert "read_filtering" in ctx


def test_context_qc_links_production_export_fields(toy_run: Path) -> None:
    """Sidecar retains legacy top-level fields consumed by older tooling."""
    ctx = json.loads((toy_run / "1-CG.json").read_text(encoding="utf-8"))
    for key in (
        "num_positions",
        "total_methylated",
        "total_unmethylated",
        "avg_methylation_level",
        "avg_coverage",
    ):
        assert key in ctx
