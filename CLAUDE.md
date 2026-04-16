# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Rust wrapper around the xpoa C++ library for partial-order alignment consensus sequence generation in DNA sequence analysis. The library builds a C++ static library during the Rust build process and links against it.

## Key Build Commands

```bash
# Build the project (compiles C++ code via build.rs)
cargo build

# Run tests
cargo test

# Run a specific test
cargo test test_poa_consensus

# Build for release
cargo build --release

# Run clippy linter
cargo clippy

# Check formatting
cargo fmt --check

# Generate documentation
cargo doc --open
```

## Development

### Build Process

The C++ code is built via a custom `build.rs` script:
1. Copies `cpp_code/` to `OUT_DIR` during build
2. Runs CMake to configure and build the C++ static library (`libxpoa.a`)
3. Links against the static library, `libm`, `libz`, and `libstdc++`
4. CMake generates a build directory in `OUT_DIR/xpoa/build/`

### Building from Scratch

To rebuild the C++ code:
```bash
rm -rf target
cargo build
```

### Regenerating Bindgen Header (for header changes)

```bash
bindgen PoaDraft.h \
  -o ../src/xpoa_sys.rs \
  --allowlist-function "PoaDraft.*" \
  -- \
  -x c++ \
  -std=c++17 \
  -I. \
  -I/usr/include/c++/13 \
  -I/usr/include/x86_64-linux-gnu/c++/13
```

## Architecture

### Rust Side Structure

```
src/
├── lib.rs           # Main public API, wrapper functions for C++ functions
├── xpoa.rs          # Rust types and traits for the Poa API (TSubread, PoaResult)
├── xpoa_sys.rs      # Bindgen-generated FFI bindings to C++ functions
└── dna.rs           # DNA sequence utilities (reverse_complement function)
```

**Key Rust Components:**
- `TSubread` trait: Abstraction for input reads with methods `get_seq()`, `get_cx()`, `is_fwd()`
- `PoaSetting`: C++ struct exposed via FFI with alignment parameters
- `PoaResult`: Wrapper around C++ `Result` with `seq()` and `n_passes()` methods
- `get_default_poa_setting()`: Returns default alignment configuration

**Public API Functions in `lib.rs`:**
- `poa_consensus()`: Basic consensus generation
- `poa_consensus_fwd()`: Forward strand only variant
- `poa_consensus_v2()`: Handles mixed orientation reads (reverses non-forward reads)
- `poa_consensus_with_fwd_string()`: Optimized for pre-verified forward strands

### C++ Side Structure

```
cpp_code/
├── PoaDraft.h              # C extern interface, defines Subread and Result structs
├── PoaConsensus.h          # Core consensus computation class
├── SparsePoa.h             # Sparse partial-order alignment
├── Settings.h              # PoaSetting definition
├── Sequence.h              # DNA sequence representation
├── edlib/                  # Edlib alignment library (local alignment)
│   ├── edlib.h
│   ├── AlignmentTools.h
│   └── ...
├── PoaGraph*.h,*.cpp       # Graph-based alignment structures
├── PoaAlignmentMatrix.*    # Alignment matrix computation
└── CMakeLists.txt          # C++ build configuration
```

**Key C++ Components:**
- `PoaDraftGen()`: Main consensus generation function (all orientations)
- `PoaDraftGenWithAllFwdStrand()`: Optimized when all reads are forward strand
- `PoaGraph`: Partial-order alignment graph structure
- `SparsePoa`: Sparse alignment implementation

## Key Concepts

### Orientation Flags
- `cx == 3`: Forward strand (cx3_cnt tracked for assertions)
- `cx` with reverse complement: Non-forward reads require reversal before processing

### Version Modes
- `version: 1`: Default PoaDraft behavior
- `version: 2`: RangeFinder2 mode (newer algorithm)

### Default Parameters
```rust
PoaSetting {
    min_identity: 0.82,     // Minimum sequence identity threshold
    match_score: 3,         // Match scoring
    mismatch_score: -5,     // Mismatch penalty
    insertion_score: -2,    // Insertion penalty
    deletion_score: -2,     // Deletion penalty
    ed_unify_strand: 1,     // Edlib unify strand flag
    version: 1,             // Algorithm version
}
```

## Testing

Tests are defined in `src/lib.rs`:
- `test_poa_consensus()`: Basic consensus with identical reads
- `test_poa_consensus_range_finder2()`: Tests version 2 (RangeFinder2)
- `test_poa_consensus_very_short()`: Edge cases with short sequences

## Files of Interest for Modifications

| File | Purpose |
|------|---------|
| `src/lib.rs` | Primary Rust API, consensus functions |
| `src/xpoa.rs` | Rust type abstractions (TSubread, PoaResult) |
| `src/dna.rs` | DNA reverse complement utilities |
| `cpp_code/PoaDraft.h` | C++ FFI interface definition |
| `build.rs` | Build script for C++ compilation |
| `cpp_code/PoaConsensus.h` | Core C++ consensus logic |
| `cpp_code/SparsePoa.h` | Sparse alignment algorithm |
