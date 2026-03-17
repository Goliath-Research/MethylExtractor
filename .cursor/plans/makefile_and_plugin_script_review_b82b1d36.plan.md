---
name: Makefile and plugin script review
overview: "HDF5 Zstd plugin is built under build/dynamic/<arch>/hdf5_zstd_plugin (with MethylExtractor); copy and HDF5_PLUGIN_PATH only on make install; arch check at install."
todos:
  - id: fix-script-so-name
    content: "Fix install_hdf5_zstd_plugin.sh – use libH5Zzstd.so everywhere; add build-only mode (path arg = build dir only, no install/.bashrc); no-args = current behavior."
    status: completed
  - id: makefile-plugin-build
    content: "Makefile – add HDF5_PLUGIN_BUILD = $(DYNAMIC_DIR)/hdf5_zstd_plugin and target that produces libH5Zzstd.so (script with path)."
    status: completed
  - id: makefile-install
    content: "Makefile – install depends on plugin .so; add arch check (file on binary and .so vs current machine); copy binary + plugin to /usr/local; append HDF5_PLUGIN_PATH to $$HOME/.bashrc."
    status: completed
  - id: makefile-clean
    content: "Makefile – clean removes plugin build dir(s), e.g. rm -rf build."
    status: completed
  - id: readme-docs
    content: "README – document plugin build in tree, install only on make install, and architecture check at install."
    status: completed
isProject: false
---

# Makefile and install_hdf5_zstd_plugin.sh review

## Relationship

- **Makefile**: Builds MethylExtractor and the HDF5 Zstd plugin (plugin built into the project’s build tree). `make install` copies the binary and the plugin .so to the system and configures `.bashrc`.
- **Plugin**: [HDF5 Zstandard filter](https://github.com/aparamon/HDF5Plugin-Zstandard). MethylExtractor uses it at runtime ([src/output_formats.c](src/output_formats.c): `H5Pset_filter(..., ZSTD_FILTER, ...)` with gzip fallback if the filter is not found).

---

## Plugin build layout (implemented)

Build the plugin **under the dynamic build dir, per architecture**: `HDF5_PLUGIN_BUILD = $(DYNAMIC_DIR)/hdf5_zstd_plugin` → `build/dynamic/$(ARCH_NAME)/hdf5_zstd_plugin/`, so it sits **alongside** the MethylExtractor binary (`build/dynamic/<arch>/MethylExtractor` and `build/dynamic/<arch>/hdf5_zstd_plugin/libH5Zzstd.so`). Copy to `/usr/local/hdf5/lib/plugin` only when the user runs **`make install`**. At install time, **verify** that the binary and plugin match the current machine’s architecture. No installation to the system during `make` or `make deps`.

- **Benefits**: All dynamic artifacts for an arch in one place; plugin is just another build artifact until install; `make install` is the only place that touches system paths and `.bashrc`; `make clean` (rm -rf build) removes everything. **Make** does not re-run the plugin recipe when `libH5Zzstd.so` already exists (target has no prerequisites).

**Correct .so name**: The plugin’s CMake uses `OUTPUT_NAME H5Zzstd`, so the built file is **`libH5Zzstd.so`** (not `libh5z_zstd.so`). Any script or Makefile that copies it must use this name.

---

## Makefile – deps, build plugin, install

**deps**  
- Only `apt-get install` (build-essential, libhdf5-dev, libhts-dev, …, cmake, git, hdf5-tools). No plugin build or install here.

**Build plugin into project tree (per-arch, implemented)**  
- `HDF5_PLUGIN_BUILD = $(DYNAMIC_DIR)/hdf5_zstd_plugin` → `build/dynamic/$(ARCH_NAME)/hdf5_zstd_plugin`. Plugin lives next to the binary.  
- Target `$(PLUGIN_SO)` (no prerequisites): runs `bash scripts/install_hdf5_zstd_plugin.sh $(HDF5_PLUGIN_BUILD)`. Make only runs it when the .so is missing.  
- `install` depends on `dynamic` and `$(PLUGIN_SO)`.

**install**  
- **Architecture check (required)**: Before copying anything, verify that the built binary and plugin match the current machine. Use `file $(DYNAMIC_DIR)/MethylExtractor` and `file $(HDF5_PLUGIN_BUILD)/libH5Zzstd.so` and require that they report the expected arch (e.g. output contains `aarch64` when `$(ARCH)` is `aarch64`, or `x86-64` when `$(ARCH)` is `x86_64`). If there is a mismatch, abort with a clear error (e.g. "Architecture mismatch: binary/plugin were built for X but this machine is Y. Run make on this machine first.") so you never copy an arm64 binary/plugin onto an x86_64 host or the reverse.
- Copy the binary to `/usr/local/bin/` and chmod (unchanged).
- Create `/usr/local/hdf5/lib/plugin` and copy `$(HDF5_PLUGIN_BUILD)/libH5Zzstd.so` there (e.g. `sudo cp ...`).
- Append `export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin` to the installing user’s `~/.bashrc` if not already present (use `$$HOME`).

---

## install_hdf5_zstd_plugin.sh (implemented)

- **Build-only mode** (first argument is a path): Use that dir as `WORKDIR`; after `cd "$WORKDIR"`, set `WORKDIR="$(pwd)"` so the final `cp libH5Zzstd.so "$WORKDIR/"` uses an absolute path and works from deep inside the tree (`.../HDF5Plugin-Zstandard/build`). Clone (if needed), build, copy .so to `$WORKDIR/`; do not install to `/usr` or touch `.bashrc`. No apt in this mode (deps from make).  
- **Standalone mode** (no args): Build in `$HOME/hdf5_plugins_build`, install to `/usr/local/hdf5/lib/plugin`, append `HDF5_PLUGIN_PATH` to `~/.bashrc`.  
- Use `libH5Zzstd.so` everywhere (CMake `OUTPUT_NAME` is H5Zzstd).

---

## To-dos (completed)

- [x] **fix-script-so-name**: [scripts/install_hdf5_zstd_plugin.sh](scripts/install_hdf5_zstd_plugin.sh) – `libH5Zzstd.so` everywhere; build-only mode (path arg); no-args = standalone install; `WORKDIR="$(pwd)"` after first cd for correct copy path.
- [x] **makefile-plugin-build**: Makefile – `HDF5_PLUGIN_BUILD = $(DYNAMIC_DIR)/hdf5_zstd_plugin`; target `$(PLUGIN_SO)` runs script with that path; deps include cmake, git, hdf5-tools.
- [x] **makefile-install**: Makefile – `install` depends on `dynamic` and `$(PLUGIN_SO)`; arch check via `file` on binary and .so; copy binary + plugin; append `HDF5_PLUGIN_PATH` to `$$HOME/.bashrc`.
- [x] **makefile-clean**: Makefile – `rm -rf build` (removes plugin under build/dynamic/$(ARCH_NAME)/).
- [x] **readme-docs**: README – plugin under `build/dynamic/<arch>/hdf5_zstd_plugin/`; optional “System install” section; arch check noted.

---

## Implemented summary

1. **Script** [scripts/install_hdf5_zstd_plugin.sh](scripts/install_hdf5_zstd_plugin.sh): `libH5Zzstd.so` everywhere; build-only when first arg is path (clone/build into that dir, set `WORKDIR="$(pwd)"` after first cd so copy works); standalone (no args) = build in $HOME, install to /usr, update .bashrc.
2. **Makefile – plugin**: `HDF5_PLUGIN_BUILD = $(DYNAMIC_DIR)/hdf5_zstd_plugin`; target `$(PLUGIN_SO)` with no prerequisites runs the script; Make skips when .so exists. deps add cmake, git, hdf5-tools.
3. **Makefile – install**: Depends on `dynamic` and `$(PLUGIN_SO)`. Arch check (`file` on binary and .so, require aarch64 or x86-64 to match); then copy binary, mkdir + cp plugin to /usr/local/hdf5/lib/plugin, append HDF5_PLUGIN_PATH to $$HOME/.bashrc if missing.
4. **clean**: `rm -rf build` (covers plugin under build/dynamic/$(ARCH_NAME)/).
5. **README**: Plugin path `build/dynamic/<arch>/hdf5_zstd_plugin/`; "System install (optional)" with make install and arch-check note.