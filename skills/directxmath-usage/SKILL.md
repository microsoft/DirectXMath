---
name: directxmath-usage
description: >-
  Guide for integrating DirectXMath into new projects and an overview of the API, types, and classes.
  Use this skill when asked about how to use DirectXMath, integrate it into a project, or for an
  overview of available functionality.
license: MIT
metadata:
  author: chuckw
  version: "1.0"
---

# DirectXMath Usage

## Overview

DirectXMath is a header-only SIMD C++ linear algebra library for use in games, graphics engines, and other performance-sensitive applications. It provides vector, matrix, and quaternion math optimized for SSE/SSE2 (x86/x64) and ARM-NEON (ARM64) instruction sets.

- **Namespace**: `DirectX`
- **Minimum C++ Standard**: C++11 (with conditional C++14/17/20 features)
- **Minimum Compiler**: Visual C++ 2017 (`_MSC_VER >= 1910`), or equivalent Clang/GCC
- **Repository**: <https://github.com/microsoft/DirectXMath>
- **Documentation**: <https://github.com/microsoft/DirectXMath/wiki>, <https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal>
- **NuGet Packages**: `directxmath`
- **vcpkg Port**: `directxmath`

> An optional C++ wrapper, SimpleMath, is available that provides a more user-friendly interface to the DirectXMath API with fewer alignment and usage restrictions. It is part of the DirectX Tool Kit for DirectX 11 and DirectX 12, and is intended to mimic the XNA Game Studio framework math library.

## Integration Methods

### vcpkg manifest-mode (Recommended)

In your `vcpkg.json` file, add the following:

```json
{
  "$schema": "https://raw.githubusercontent.com/microsoft/vcpkg-tool/main/docs/vcpkg.schema.json",
  "dependencies": [
    "directxmath"
  ]
}
```

### vcpkg (classic)

```bash
vcpkg install directxmath
```

Then in your `CMakeLists.txt`:

```cmake
find_package(directxmath CONFIG REQUIRED)
target_link_libraries(YourTarget PRIVATE Microsoft::DirectXMath)
```

Features: `dx11` (Spherical Harmonics math library for DirectX 11), `dx12` (Spherical Harmonics math library for DirectX 12), `xdsp` (Digital Signal Processing library). Triplets: `x64-windows`, `x64-linux`, `arm64-windows`, etc.

### NuGet

Install the `directxmath` package from nuget.org.

### Windows SDK

The Windows SDK ships with DirectXMath, though it may not be the latest version. Include via:

```cpp
#include <DirectXMath.h>
```

## Basic Usage

```cpp
#include <DirectXMath.h>
using namespace DirectX;

// Load data from storage type to SIMD register
XMFLOAT3 position{ 1.0f, 2.0f, 3.0f };
XMVECTOR v = XMLoadFloat3(&position);

// Perform computation in SIMD
v = XMVector3Normalize(v);

// Store result back to memory
XMFLOAT3 result;
XMStoreFloat3(&result, v);
```

## Key Concepts

### Load / Compute / Store Pattern

DirectXMath separates storage types (for memory layout) from the computation type (`XMVECTOR`). Always:

1. **Load** from a storage type (`XMFLOAT3`, `XMFLOAT4`, etc.) into `XMVECTOR` or `XMMATRIX` using `XMLoad*` functions.
2. **Compute** using `XMVector*`, `XMMatrix*`, `XMQuaternion*`, or `XMColor*` functions.
3. **Store** results back to a storage type using `XMStore*` functions.

### Calling Conventions

Functions that accept `XMVECTOR` or `XMMATRIX` parameters use `XM_CALLCONV` and the parameter typedefs `FXMVECTOR`, `GXMVECTOR`, `HXMVECTOR`, `CXMVECTOR`, `FXMMATRIX`, and `CXMMATRIX` to maximize register usage across platforms.

### Alignment

`XMVECTOR` and `XMMATRIX` require 16-byte alignment. Storage types with an `A` suffix (e.g., `XMFLOAT4A`) are aligned; the unsuffixed variants (e.g., `XMFLOAT4`) are unaligned and safe for use in arbitrary data structures.

## API Reference

For detailed API signatures and function families, see the [reference overview](reference/overview.md).

The canonical API documentation is on [Microsoft Learn](https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal). All function signatures can be discovered in the public headers under the `Inc/` directory.

## Auxiliary Libraries

### SHMath (Spherical Harmonics)

The `SHMath/` directory provides spherical harmonic math functions for lighting computations. Unlike the core library, SHMath includes `.cpp` files that must be compiled.

See [reference overview - SHMath section](reference/overview.md#shmath-spherical-harmonics) for details.

### XDSP (Digital Signal Processing)

The `XDSP/` directory provides header-only FFT and DSP functions built on DirectXMath.

See [reference overview - XDSP section](reference/overview.md#xdsp-digital-signal-processing) for details.

## Platform Support

| Platform | SIMD Backend |
| --- | --- |
| Windows x86/x64 | SSE/SSE2 (baseline), up to AVX2 |
| Windows ARM64/ARM64EC | ARM-NEON |
| Linux x86_64 | SSE/SSE2 (baseline), up to AVX2 |
| Linux ARM64 | ARM-NEON |
| Any (fallback) | Pure C++ (`_XM_NO_INTRINSICS_`) |

## Additional Resources

- [Getting Started (Microsoft Learn)](https://learn.microsoft.com/windows/win32/dxmath/pg-xnamath-getting-started)
- [API Reference (Microsoft Learn)](https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal)
- [GitHub Wiki](https://github.com/microsoft/DirectXMath/wiki/)
- [Introductory Blog Post](https://walbourn.github.io/introducing-directxmath/)
