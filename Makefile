# Makefile for MethylExtractor

# Detect architecture
ARCH := $(shell uname -m)
ifeq ($(ARCH),aarch64)
    ARCH_NAME := arm64
    HDF5_LIB_PATH := /usr/lib/aarch64-linux-gnu/hdf5/serial
else ifeq ($(ARCH),x86_64)
    ARCH_NAME := x64
    HDF5_LIB_PATH := /usr/lib/x86_64-linux-gnu/hdf5/serial
else
    $(error Unsupported architecture: $(ARCH))
endif

CC = gcc
CFLAGS = -O2 -Wall -Iinclude -Iinclude/cjson -I/usr/include/hdf5 -I/usr/include/htslib -I/usr/include
DEBUG_CFLAGS = -g -O0 -Wall -Iinclude -Iinclude/cjson -I/usr/include/hdf5 -I/usr/include/htslib -I/usr/include -DDEBUG
LDFLAGS = -Llib -L/usr/lib/$(ARCH_NAME)-linux-gnu

# Dependencies
HTSLIB_DIR = /usr
HDF5_DIR = /usr

# Libraries
HTSLIB_LIBS = -lhts
HDF5_LIBS = -lhdf5_serial -lhdf5_serial_hl

# Additional libraries for static linking
STATIC_LIBS = -lz -lm -ldl -lpthread -lbz2 -llzma -lcurl -lcrypto -lssl -lzstd

# Output directories
STATIC_DIR = build/static/$(ARCH_NAME)
DYNAMIC_DIR = build/dynamic/$(ARCH_NAME)
DEBUG_DIR = build/debug/$(ARCH_NAME)
HDF5_PLUGIN_BUILD = $(DYNAMIC_DIR)/hdf5_zstd_plugin
PLUGIN_SO = $(HDF5_PLUGIN_BUILD)/libH5Zzstd.so

# Source files
SRCS = src/main.c src/bam_processing.c src/output_formats.c src/extraction_export.c src/utils.c src/read_level.c src/mhap.c src/cjson/cJSON.c

# Targets
# Note: `deps` is intentionally NOT a prerequisite of `all` so that a plain
# `make` never runs `sudo apt-get`. Run `make deps` once to install system
# packages, then `make`.
all: dynamic

deps:
	@echo "Installing required dependencies..."
	sudo apt-get update
	sudo apt-get install -y \
		build-essential \
		cmake \
		git \
		libhdf5-dev \
		libhts-dev \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		libcurl4-gnutls-dev \
		libssl-dev \
		libzstd-dev \
		zstd \
		hdf5-tools
	@echo "Dependencies installed successfully"

# HDF5 Zstd plugin (built into project tree, installed only on make install).
# HDF5_PLUGIN_BUILD = build/dynamic/$(ARCH_NAME)/hdf5_zstd_plugin, so the plugin
# is always built into the directory matching this machine's architecture.
$(PLUGIN_SO):
	@echo "Building HDF5 Zstd plugin for $(ARCH) into $(HDF5_PLUGIN_BUILD)..."
	@mkdir -p $(HDF5_PLUGIN_BUILD)
	bash scripts/install_hdf5_zstd_plugin.sh $(HDF5_PLUGIN_BUILD)

# Convenience alias: `make plugin` builds the Zstd filter into the arch-correct
# tree and prints the HDF5_PLUGIN_PATH to export for in-tree (no-install) use.
plugin: $(PLUGIN_SO)
	@echo "HDF5 Zstd plugin ready: $(PLUGIN_SO)"
	@echo "To use it without 'make install', export:"
	@echo "  export HDF5_PLUGIN_PATH=$(abspath $(HDF5_PLUGIN_BUILD))"

# Note: Static linking is challenging due to missing static libraries for HTSlib and HDF5 dependencies
# (like libdeflate, rans, arith, fqz, tok3, and szip). Consider using dynamic linking instead.
# If static linking is required, you may need to build HTSlib and HDF5 from source with static libraries.
static: $(STATIC_DIR)/MethylExtractor

dynamic: $(DYNAMIC_DIR)/MethylExtractor

debug: $(DEBUG_DIR)/MethylExtractor

$(STATIC_DIR)/MethylExtractor: $(SRCS)
	mkdir -p $(STATIC_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -L$(HTSLIB_DIR)/lib -L$(HDF5_LIB_PATH) $(HTSLIB_LIBS) $(HDF5_LIBS) $(STATIC_LIBS) -static

$(DYNAMIC_DIR)/MethylExtractor: $(SRCS)
	mkdir -p $(DYNAMIC_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -L$(HTSLIB_DIR)/lib -L$(HDF5_LIB_PATH) $(HTSLIB_LIBS) $(HDF5_LIBS) -lpthread -lm

$(DEBUG_DIR)/MethylExtractor: $(SRCS)
	mkdir -p $(DEBUG_DIR)
	$(CC) $(DEBUG_CFLAGS) -o $@ $^ $(LDFLAGS) -L$(HTSLIB_DIR)/lib -L$(HDF5_LIB_PATH) $(HTSLIB_LIBS) $(HDF5_LIBS) -lpthread -lm

install: dynamic $(PLUGIN_SO)
	@echo "Checking architecture..."
	@bin_arch=$$(file -b $(DYNAMIC_DIR)/MethylExtractor); \
	plug_arch=$$(file -b $(PLUGIN_SO)); \
	case "$(ARCH)" in \
		aarch64) want="aarch64";; \
		x86_64) want="x86-64";; \
		*) want="";; \
	esac; \
	echo "$$bin_arch" | grep -q "$$want" || { echo "Architecture mismatch: binary is for $$bin_arch but this machine is $(ARCH). Run make on this machine first."; exit 1; }; \
	echo "$$plug_arch" | grep -q "$$want" || { echo "Architecture mismatch: plugin is for $$plug_arch but this machine is $(ARCH). Run make on this machine first."; exit 1; }
	@echo "Installing binary and HDF5 Zstd plugin..."
	sudo cp $(DYNAMIC_DIR)/MethylExtractor /usr/local/bin/
	sudo chmod +x /usr/local/bin/MethylExtractor
	sudo mkdir -p /usr/local/hdf5/lib/plugin
	sudo cp $(PLUGIN_SO) /usr/local/hdf5/lib/plugin/
	@if ! grep -q "HDF5_PLUGIN_PATH" $$HOME/.bashrc 2>/dev/null; then \
		echo "export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin" >> $$HOME/.bashrc; \
		echo "Appended HDF5_PLUGIN_PATH to $$HOME/.bashrc"; \
	fi

clean:
	rm -rf build
	@# Removes all build artifacts including the plugin under build/dynamic/$(ARCH_NAME)/

# Path of the dynamic binary for this machine, for CI and packaging scripts.
print-dynamic-bin:
	@echo $(DYNAMIC_DIR)/MethylExtractor

.PHONY: all deps plugin static dynamic debug install clean print-dynamic-bin
