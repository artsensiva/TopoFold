# TopoFold: Enterprise Production Binary Wheel Distribution Guide

## Overview
This document specifies how to compile and package TopoFold into hardened, stripped, single-binary Python wheels (`.whl`) without distributing Rust or Python source code.

## Enterprise Release Profile Configuration
In `Cargo.toml`, the release profile is configured for maximum performance and intellectual property protection:

```toml
[profile.release]
opt-level = 3          # Aggressive compiler optimization
lto = "fat"            # Whole-program Link Time Optimization across all crates
codegen-units = 1      # Maximize cross-crate inlining and devirtualization
panic = "abort"        # Eliminates unwinding landing pads and stack trace symbols
strip = "symbols"      # Strips all debug symbols and internal function names
overflow-checks = false
```

## Building Stripped Binary Wheels with Maturin

### 1. Build Isolated Wheel
Run `maturin build` with `--release` and `--strip`:

```bash
# Build binary wheel for the current platform
maturin build --release --strip --out dist
```

This compiles `crates/topofold-python` and all workspace dependencies into a standalone binary extension (`topofold.abi3.so` or `topofold.cpython-*.so`), embedded within a `.whl` archive under `dist/`.

### 2. Verify Symbol Stripping
Verify that debug symbols and internal function identifiers have been removed:

```bash
# Extract the compiled shared library
unzip -p dist/topofold-*.whl "topofold/*.so" > /tmp/topofold_test.so 2>/dev/null || unzip -p dist/topofold-*.whl "*.so" > /tmp/topofold_test.so

# Check with nm / file
file /tmp/topofold_test.so
# Output: ELF 64-bit LSB shared object, ..., stripped

nm -D /tmp/topofold_test.so | grep -E "PyInit_topofold"
# Only the required PyO3 C-Python initialization symbol is exposed
```

### 3. Source-Free Wheel Verification
The resulting wheel contains **no Rust source code** (`.rs`), no internal build manifests (`Cargo.toml`), and only exposes the compiled native extension.

### 4. Direct Installation into Clean Environment
```bash
pip install dist/topofold-*.whl
```
