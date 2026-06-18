---
name: MethylExtractor release CI
overview: Per-arch Azure DevOps release and PR pipelines that build native binaries, package tarballs via MethylPipeline scripts, and publish to the methyl-extractor Universal Packages feed on git tags.
azure_devops:
  type: Task
  title: Register MethylExtractor release + PR pipelines; enable Universal Package publish per arch
  work_item_id: null  # set to AB# after creating Task in Azure DevOps Boards
  parent_epics:
    - repo: MethylPipeline
      plan: docs/plans/devops-ci-cd-release.plan.md
      task_id: me-release-ci
    - repo: MethylPipeline
      plan: docs/plans/production-gpu-worker-layout.plan.md
      task_id: extractor-ci
todos:
  - id: pr-pipeline
    content: Add ci/azure-pipelines-pr.yml (make + smoke on pull requests, no publish)
    status: completed
  - id: release-arm64
    content: Add ci/azure-pipelines-release-arm64.yml on pool build-arm64; tag v* and manual releaseVersion param
    status: completed
  - id: release-x64
    content: Add ci/azure-pipelines-release-x64.yml on Microsoft-hosted ubuntu-latest; tag v*
    status: completed
  - id: universal-publish
    content: az artifacts universal publish per arch to feed methyl-extractor (SemVer without leading zeros)
    status: completed
  - id: multi-repo-checkout
    content: Checkout MethylPipeline for scripts/package_methyl_extractor.sh and detect_platform.sh
    status: completed
  - id: ci-readme
    content: Document pipeline registration, ARM64 agent pool, and release flow in ci/README.md
    status: completed
isProject: false
---

# MethylExtractor release CI

## Scope

Build and publish **arch-specific tarballs** (`methyl-extractor-linux-aarch64`, `methyl-extractor-linux-amd64`) independently of MethylPipeline wheel releases. Orchestration (assemble + deploy to `/work/epimethyl`) is owned by MethylPipeline.

## Pipelines

| YAML | Pipeline name | Pool | Trigger |
|------|---------------|------|---------|
| `ci/azure-pipelines-pr.yml` | MethylExtractor-PR | `ubuntu-latest` | Pull requests |
| `ci/azure-pipelines-release-arm64.yml` | MethylExtractor-Release-ARM64 | `build-arm64` | Tags `v*`, manual |
| `ci/azure-pipelines-release-x64.yml` | MethylExtractor-Release-x64 | `ubuntu-latest` | Tags `v*`, manual |

Tag `v2026.6.1` runs **both** release pipelines. Each publishes one Universal Package version to feed **`methyl-extractor`**.

## Release flow (cross-repo)

1. Tag MethylExtractor `v2026.6.1` → both release pipelines publish.
2. Tag MethylPipeline when Python/worker code changes.
3. Run **Epimethyl-Release-Assemble** (MethylPipeline) with pinned ME + MP versions.
4. Approve **Epimethyl-Release-Deploy** to promote `/work/epimethyl/current`.

See MethylPipeline [`docs/deployment/production_release.md`](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline?path=/docs/deployment/production_release.md).

## Prerequisites

- Artifacts feed **`methyl-extractor`** — build service **Contributor**
- Multi-repo checkout permission for **Development/MethylPipeline**
- Self-hosted pool **`build-arm64`** with at least one online agent

Details: [`ci/README.md`](../../ci/README.md).

## Versioning

- Git tags: `v2026.6.1` (no leading zeros in minor, e.g. not `v2026.06.1`)
- Manual test runs from `main`: set pipeline parameter **`releaseVersion`** (e.g. `2026.6.1-test`)
