# Implementation plans (Azure DevOps traceability)

This repo owns **MethylExtractor** build and release CI. Cross-repo Epics and orchestration plans live in [MethylPipeline `docs/plans`](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline?path=/docs/plans).

## Work item mapping (this repo)

| Plan file | ADO type | Parent Epic (MethylPipeline plan) | Task id |
|-----------|----------|-----------------------------------|---------|
| [`methyl-extractor-release-ci.plan.md`](methyl-extractor-release-ci.plan.md) | **Task** | [DevOps CI/CD release](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline?path=/docs/plans/devops-ci-cd-release.plan.md) | `me-release-ci` |
| *(same plan)* | **Task** | [Production GPU worker layout](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline?path=/docs/plans/production-gpu-worker-layout.plan.md) | `extractor-ci` |

Create the **Task** in Azure DevOps under the appropriate Epic. Use the plan `todos[].content` for the task title. Mark **Closed** when merged.

Assemble and deploy pipelines are implemented in MethylPipeline only.

## Commit message convention

Link commits to work items (`AB#<id>`). Azure DevOps updates linked items when mention tracking is enabled.

```
Add ARM64 release pipeline with Universal Package publish (AB#1234)

Implements methyl-extractor-release-ci plan task release-pipelines.
Parent Epic: DevOps CI/CD release (AB#1200).
```

## Branch and PR workflow

1. Create or pick a **Task** under the Epic in Azure DevOps Boards.
2. Branch: `feature/AB1234-me-release-arm64` or `task/1234-universal-publish`.
3. Reference the plan file in the PR description.
4. PR title: `Enable Universal Package publish on ARM64 release (AB#1234)`.

## Related docs

- [`../../ci/README.md`](../../ci/README.md) — pipeline registration and ARM64 agent setup
- MethylPipeline [`docs/deployment/production_release.md`](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline?path=/docs/deployment/production_release.md) — full release workflow
