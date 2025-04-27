# Dockerfile for epiextractor

FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libhdf5-dev \
    libhts-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy source code
COPY . /app

# Build the application
RUN make all

# Set entrypoint
ENTRYPOINT ["./build/dynamic/epiextractor"]

# Default command
CMD ["--help"] 