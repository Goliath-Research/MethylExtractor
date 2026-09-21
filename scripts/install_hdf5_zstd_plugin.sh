#!/usr/bin/env bash

set -e

# When first argument is a path: build-only mode (build into that dir, no install, no .bashrc).
# When no argument: standalone mode (build in $HOME/hdf5_plugins_build, install to /usr, update .bashrc).
BUILD_ONLY_DIR="$1"

if [ -n "$BUILD_ONLY_DIR" ]; then
	echo "========================================="
	echo " Building HDF5 Zstandard plugin (build-only)"
	echo "========================================="
	WORKDIR="$BUILD_ONLY_DIR"
	INSTALL_TO_SYSTEM=false
else
	echo "========================================="
	echo " Installing HDF5 Zstandard plugin"
	echo "========================================="
	WORKDIR="$HOME/hdf5_plugins_build"
	INSTALL_TO_SYSTEM=true
fi

# ---------------------------------------------------
# 1. Install required packages (skip in build-only; caller ensures deps)
# ---------------------------------------------------

if [ "$INSTALL_TO_SYSTEM" = true ]; then
	echo "Installing dependencies..."
	sudo apt update
	sudo apt install -y \
		build-essential \
		cmake \
		git \
		libhdf5-dev \
		libzstd-dev \
		hdf5-tools
fi

# ---------------------------------------------------
# 2. Create working directory
# ---------------------------------------------------

mkdir -p "$WORKDIR"
cd "$WORKDIR"
# Use absolute path so copy works when we are deep inside the tree (build-only mode)
WORKDIR="$(pwd)"

# ---------------------------------------------------
# 3. Download Zstd HDF5 filter source
# ---------------------------------------------------

echo "Downloading HDF5 Zstd plugin..."

if [ ! -d "HDF5Plugin-Zstandard" ]; then
	git clone https://github.com/aparamon/HDF5Plugin-Zstandard.git
fi

cd HDF5Plugin-Zstandard

# ---------------------------------------------------
# 4. Build plugin
# ---------------------------------------------------

echo "Compiling plugin..."

mkdir -p build
cd build

# Upstream declares cmake_minimum_required(VERSION 2.8.10); CMake >= 4 rejects it.
# Older CMake ignores this unused cache variable.
cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
make -j$(nproc)

# Built library is libH5Zzstd.so (CMake OUTPUT_NAME H5Zzstd)
PLUGIN_SO="libH5Zzstd.so"
if [ ! -f "$PLUGIN_SO" ]; then
	echo "Error: $PLUGIN_SO not found after build" >&2
	exit 1
fi

if [ "$INSTALL_TO_SYSTEM" = false ]; then
	# Build-only: copy into the requested build dir root so Makefile can find it
	cp "$PLUGIN_SO" "$WORKDIR/"
	echo "========================================="
	echo " HDF5 Zstd plugin built: $WORKDIR/$PLUGIN_SO"
	echo "========================================="
	exit 0
fi

# ---------------------------------------------------
# 5. Install plugin (standalone mode only)
# ---------------------------------------------------

PLUGIN_DIR="/usr/local/hdf5/lib/plugin"

echo "Installing plugin into $PLUGIN_DIR"

sudo mkdir -p "$PLUGIN_DIR"
sudo cp "$PLUGIN_SO" "$PLUGIN_DIR/"

# ---------------------------------------------------
# 6. Configure plugin path
# ---------------------------------------------------

echo "Configuring HDF5_PLUGIN_PATH..."

if ! grep -q "HDF5_PLUGIN_PATH" ~/.bashrc ; then
	echo "export HDF5_PLUGIN_PATH=$PLUGIN_DIR" >> ~/.bashrc
fi

export HDF5_PLUGIN_PATH=$PLUGIN_DIR

# ---------------------------------------------------
# 7. Verify installation
# ---------------------------------------------------

echo "-----------------------------------------"
echo "Installed plugin files:"
ls -l "$PLUGIN_DIR"

echo ""
echo "HDF5 plugin path:"
echo "$HDF5_PLUGIN_PATH"

echo ""
echo "If you see $PLUGIN_SO above, installation succeeded."
echo ""
echo "Restart your shell or run:"
echo "source ~/.bashrc"

echo "========================================="
echo " HDF5 Zstd plugin installation complete"
echo "========================================="
