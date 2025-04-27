# Makefile for MethylExtractor

CC = gcc
CFLAGS = -O2 -Wall -Iinclude -I/usr/include/hdf5/serial -I/usr/include/htslib
DEBUG_CFLAGS = -g -O0 -Wall -Iinclude -I/usr/include/hdf5/serial -I/usr/include/htslib -DDEBUG
LDFLAGS = -Llib

# Dependencies
HTSLIB_DIR = /usr
HDF5_DIR = /usr

# Libraries
HTSLIB_LIBS = -lhts
HDF5_LIBS = -lhdf5 -lhdf5_hl

# Additional libraries for static linking
STATIC_LIBS = -lz -lm -ldl -lpthread -lbz2 -llzma -lcurl -lcrypto -lssl

# Output directories
STATIC_DIR = build/static
DYNAMIC_DIR = build/dynamic
DEBUG_DIR = build/debug

# Targets
all: dynamic

# Note: Static linking is challenging due to missing static libraries for HTSlib and HDF5 dependencies
# (like libdeflate, rans, arith, fqz, tok3, and szip). Consider using dynamic linking instead.
# If static linking is required, you may need to build HTSlib and HDF5 from source with static libraries.
static: $(STATIC_DIR)/MethylExtractor

dynamic: $(DYNAMIC_DIR)/MethylExtractor

debug: $(DEBUG_DIR)/MethylExtractor

$(STATIC_DIR)/MethylExtractor: src/MethylExtractor.c
	mkdir -p $(STATIC_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -L$(HTSLIB_DIR)/lib -L$(HDF5_DIR)/lib/x86_64-linux-gnu/hdf5/serial $(HTSLIB_LIBS) $(HDF5_LIBS) $(STATIC_LIBS) -static

$(DYNAMIC_DIR)/MethylExtractor: src/MethylExtractor.c
	mkdir -p $(DYNAMIC_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -L$(HTSLIB_DIR)/lib -L$(HDF5_DIR)/lib/x86_64-linux-gnu/hdf5/serial $(HTSLIB_LIBS) $(HDF5_LIBS) -lpthread

$(DEBUG_DIR)/MethylExtractor: src/MethylExtractor.c
	mkdir -p $(DEBUG_DIR)
	$(CC) $(DEBUG_CFLAGS) -o $@ $^ $(LDFLAGS) -L$(HTSLIB_DIR)/lib -L$(HDF5_DIR)/lib/x86_64-linux-gnu/hdf5/serial $(HTSLIB_LIBS) $(HDF5_LIBS) -lpthread

clean:
	rm -rf build

.PHONY: all static dynamic debug clean 