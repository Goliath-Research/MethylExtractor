# Contributing to MethylExtractor

Thank you for your interest in contributing to MethylExtractor! We welcome contributions from the community to help improve this tool.

## Getting Started

1.  **Fork the repository** on GitHub.
2.  **Clone your fork** locally:
    ```bash
    git clone https://github.com/your-username/MethylExtractor.git
    cd MethylExtractor
    ```
3.  **Create a new branch** for your feature or bug fix:
    ```bash
    git checkout -b feature/my-new-feature
    ```

## Development Environment

-   **OS**: Linux (Debian/Ubuntu recommended)
-   **Compiler**: GCC
-   **Dependencies**: HTSlib, HDF5, Zstd

To install dependencies and build the project:
```bash
make deps
make
```

## Directory Structure

-   `src/`: Source code files (`.c`)
-   `include/`: Header files (`.h`)
-   `tests/`: Test scripts and data
-   `build/`: Output directory for binaries

## Coding Standards

-   Follow the existing code style (indentation, variable naming).
-   Keep functions focused and modular.
-   Add comments for complex logic.
-   Ensure all new code is covered by tests (if applicable).

## Submitting Changes

1.  **Commit your changes** with descriptive commit messages.
2.  **Push to your fork**:
    ```bash
    git push origin feature/my-new-feature
    ```
3.  **Open a Pull Request** on the main repository.
4.  Describe your changes clearly and reference any related issues.

## Reporting Issues

If you find a bug or have a feature request, please open an issue on GitHub. Provide as much detail as possible, including:
-   Steps to reproduce the issue
-   Expected vs. actual behavior
-   System information (OS, compiler version)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
