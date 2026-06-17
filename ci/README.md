# MethylExtractor Azure DevOps pipelines

Register in the **MethylExtractor** repository (Development project).

| YAML | Suggested pipeline name | Trigger |
|------|-------------------------|---------|
| [`azure-pipelines-pr.yml`](azure-pipelines-pr.yml) | MethylExtractor-PR | Pull requests |
| [`azure-pipelines-release.yml`](azure-pipelines-release.yml) | MethylExtractor-Release | Tags `v*` |

The release pipeline checks out **MethylPipeline** for `scripts/package_methyl_extractor.sh` only.

**Prerequisites:** agent pools `GPU-ARM64` and `GPU-x86_64`; Contributor on Artifacts feed `methyl-extractor`; multi-repo checkout permission for MethylPipeline.

Orchestration (assemble/deploy) lives in the [MethylPipeline](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline) repo `ci/` folder.
