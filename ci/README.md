# MethylExtractor Azure DevOps pipelines

Register in the **MethylExtractor** repository (Development project).

| YAML | Suggested pipeline name | Pool | Trigger |
|------|-------------------------|------|---------|
| [`azure-pipelines-pr.yml`](azure-pipelines-pr.yml) | MethylExtractor-PR | Microsoft-hosted `ubuntu-latest` | Pull requests |
| [`azure-pipelines-release-arm64.yml`](azure-pipelines-release-arm64.yml) | MethylExtractor-Release-ARM64 | Self-hosted `build-arm64` | Tags `v*` |
| [`azure-pipelines-release-x64.yml`](azure-pipelines-release-x64.yml) | MethylExtractor-Release-x64 | Microsoft-hosted `ubuntu-latest` | Tags `v*` |

Tag `v2026.6.1` on MethylExtractor runs **both** release pipelines. Each publishes its arch to the `methyl-extractor` Universal Packages feed (`methyl-extractor-linux-aarch64`, `methyl-extractor-linux-amd64`).

Release pipelines check out **MethylPipeline** for `scripts/package_methyl_extractor.sh` only.

## Prerequisites

1. **Artifacts feed** `methyl-extractor` — grant **Contributor** to the build service (`Project Collection Build Service` or `Development Build Service`).
2. **Multi-repo checkout** — allow this pipeline to check out `Development/MethylPipeline`.
3. **Agent pool `build-arm64`** — one self-hosted ARM64 Linux agent (small VM; no GPU).

Orchestration (assemble/deploy) lives in the [MethylPipeline](https://dev.azure.com/EpiMethyl/Development/_git/MethylPipeline) repo `ci/` folder.

## Register a pipeline

1. **Pipelines** → **New pipeline** → **Azure Repos Git** → **MethylExtractor**.
2. **Existing Azure Pipelines YAML file** → branch `main` → pick a path from the table above.
3. Save and rename to the suggested pipeline name.
4. Repeat for each YAML file (three pipelines total).

## ARM64 build VM (Azure CLI)

You can create the ARM64 Ubuntu VM from Windows, Linux, or [Azure Cloud Shell](https://shell.azure.com) — no portal required.

```bash
# Login and pick subscription (once per machine)
az login
az account set --subscription "<subscription-name-or-id>"

# Resource group (pick a region that offers Arm64 VMs, e.g. eastus, westus2)
az group create --name rg-methyl-build --location eastus

# ARM64 Ubuntu 24.04 — small compile-only size (no GPU)
az vm create \
  --resource-group rg-methyl-build \
  --name methyl-build-arm64 \
  --image Canonical:ubuntu-24_04-lts:server-arm64:latest \
  --size Standard_D2pls_v6 \
  --admin-username azureuser \
  --generate-ssh-keys \
  --public-ip-sku Standard

# Optional: open SSH from your IP only
MY_IP=$(curl -s https://api.ipify.org)
az vm open-port --resource-group rg-methyl-build --name methyl-build-arm64 --port 22 \
  --priority 1001 --source-address "$MY_IP/32"
```

SSH in and install build dependencies (or let the pipeline run `make deps` on first build):

```bash
az vm show -d --resource-group rg-methyl-build --name methyl-build-arm64 \
  --query publicIps -o tsv
# ssh azureuser@<public-ip>

sudo apt-get update
sudo apt-get install -y curl git jq
```

### Agent pool + pipeline agent

Create pool **build-arm64** once (no `az pipelines pool create` — use DevOps UI or REST):

- **UI:** Project **Development** → **Project settings** → **Agent pools** → **Add pool** → name `build-arm64`, type **Self-hosted**.
- **REST** (from any host with `az login` + PAT):

```bash
az devops invoke \
  --area distributedtask \
  --resource pools \
  --organization https://dev.azure.com/EpiMethyl \
  --project Development \
  --http-method POST \
  --api-version 7.1 \
  --in-file - <<'EOF'
{
  "name": "build-arm64",
  "poolType": "automation",
  "autoProvision": false
}
EOF
```

On the ARM64 VM, install the agent ([Linux agent docs](https://learn.microsoft.com/en-us/azure/devops/pipelines/agents/linux-agent?view=azure-devops)):

```bash
# On the VM — use a PAT with Agent Pools (read & manage)
ORG="https://dev.azure.com/EpiMethyl"
POOL="build-arm64"
PAT="<personal-access-token>"

mkdir -p ~/azagent && cd ~/azagent
curl -fsSL -u ":$PAT" \
  "$ORG/_apis/distribution/environments/agent-packages?api-version=7.1-preview.1" \
  | jq -r '.value[0].downloadUrl' \
  | xargs curl -fsSL -u ":$PAT" -o agent.tar.gz
tar xzf agent.tar.gz
./config.sh --unattended \
  --url "$ORG" \
  --auth pat \
  --token "$PAT" \
  --pool "$POOL" \
  --agent "$(hostname)" \
  --acceptTeeEula \
  --replace
sudo ./svc.sh install
sudo ./svc.sh start
```

Authorize the pool for YAML pipelines when Azure DevOps prompts on first run ([authorization](https://aka.ms/yamlauthz)).

### Cost control (optional)

```bash
# Stop VM when not building
az vm deallocate --resource-group rg-methyl-build --name methyl-build-arm64

# Start before a release tag
az vm start --resource-group rg-methyl-build --name methyl-build-arm64
```

## Release flow

1. Tag MethylExtractor `v2026.6.1` → **MethylExtractor-Release-ARM64** and **MethylExtractor-Release-x64** run in parallel.
2. Tag MethylPipeline → **MethylPipeline-Release**.
3. Run **Epimethyl-Release-Assemble** with both ME package versions pinned.
4. Approve **Epimethyl-Release-Deploy**.
