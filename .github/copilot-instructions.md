# GitHub Copilot Instructions

These instructions define how GitHub Copilot should assist with this project. The goal is to ensure consistent, high-quality code generation aligned with our conventions, stack, and best practices.

## Context

- **Project Type**: Math Library / DirectX / Direct3D
- **Project Name**: DirectXMath SIMD C++ linear algebra library
- **Language**: C++ (minimum C++11; C++14, C++17, and C++20 features used conditionally)
- **Framework / Libraries**: STL / CMake / CTest
- **Compiler Requirement**: Visual C++ 2017 or later (`_MSC_VER >= 1910`) when using MSVC

## Getting Started

- See the Getting Started guide on [Microsoft Learn](https://learn.microsoft.com/windows/win32/dxmath/pg-xnamath-getting-started).
- The recommended way to integrate *DirectXMath* into your project is by using the *vcpkg* Package Manager.
- You can make use of the nuget.org package **directxmath**.
- You can also use the library source code directly in your project or as a git submodule.
- The Windows SDK includes the DirectXMath library, although that version is not as up-to-date as other integration methods.

## General Guidelines

- **Code Style**: The project uses an .editorconfig file to enforce coding standards. Follow the rules defined in `.editorconfig` for indentation, line endings, and other formatting. Additional information can be found on the wiki at [Implementation](https://github.com/microsoft/DirectXMath/wiki/Implementation). The code requires C++11/C++14 features.
> Notable `.editorconfig` rules: C/C++ files use 4-space indentation, `crlf` line endings, and `latin1` charset — avoid non-ASCII characters in source files. HLSL files have separate indent/spacing rules defined in `.editorconfig`.
- **Documentation**: The project provides documentation on [Microsoft Learn](https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal) with additional wiki pages available on [GitHub](https://github.com/microsoft/DirectXMath/wiki/).
- **Error Handling**: The majority of functions have no error conditions and do not throw C++ exceptions which is why they are marked `noexcept`. A few functions have `bool` results to indicate success or failure.
- **Testing**: Unit tests for this project are implemented in this repository [Test Suite](https://github.com/walbourn/directxmathtest/) and can be run using CTest per the instructions at [Test Documentation](https://github.com/walbourn/directxmathtest/wiki).
- **Security**: This project uses secure coding practices from the Microsoft Secure Coding Guidelines, and is subject to the `SECURITY.md` file in the root of the repository.
- **Dependencies**: The project has minimal dependencies, primarily relying on compiler intrinsics. It is designed to be self-contained and portable across different platforms and toolsets.
- **Continuous Integration**: This project implements GitHub Actions for continuous integration, ensuring that all code changes are tested and validated before merging. This includes building the project for a number of configurations and toolsets, running unit tests, and static code analysis including GitHub super-linter, CodeQL, and MSVC Code Analysis.
- **Code of Conduct**: The project adheres to the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/). All contributors are expected to follow this code of conduct in all interactions related to the project.

## File Structure

```txt
.azuredevops/ # Azure DevOps pipeline configuration and policy files.
.github/      # GitHub Actions workflow files and linter configuration files.
.nuget/       # NuGet package configuration files.
build/        # Miscellaneous build files and scripts (CMake config templates, pkg-config, PowerShell helpers).
Inc/          # DirectXMath public and implementation files. The library is header-only, so all files are in this directory.
Extensions/   # Extensions to the DirectXMath library with standalone SSE-level specific functions for runtime selection of SIMD instruction set.
MatrixStack/  # D3DX9-like matrix stack implementation for DirectXMath.
SHMath/       # Spherical harmonic functions using DirectXMath.
Stereo3D/     # Stereo 3D projection matrix functions using DirectXMath created for HoloLens.
XDSP/         # Digital Signal Processing (DSP) functions using DirectXMath.
Tests/        # Tests are designed to be cloned from a separate repository at this location.
wiki/         # Local clone of the GitHub wiki documentation repository.
```

> The `Extensions` are not needed if building the library using `/arch:AVX` or `/arch:AVX2` which causes the DirectXMath library to build utilizing the additional SIMD instructions.

### Public Header Files (`Inc/`)

| File | Purpose |
| --- | --- |
| `DirectXMath.h` | Core library: types, constants, vector and matrix math functions |
| `DirectXMathConvert.inl` | Implementation of type conversion / load / store functions |
| `DirectXMathMatrix.inl` | Implementation of matrix math functions |
| `DirectXMathMisc.inl` | Implementation of miscellaneous functions (quaternion, color, plane, etc.) |
| `DirectXMathVector.inl` | Implementation of vector math functions |
| `DirectXCollision.h` / `.inl` | Bounding volumes: `BoundingSphere`, `BoundingBox`, `BoundingOrientedBox`, `BoundingFrustum` |
| `DirectXColors.h` | Predefined color constants (`DirectX::Colors::*` sRGB, `DirectX::ColorsLinear::*` linear) |
| `DirectXPackedVector.h` / `.inl` | GPU-friendly packed vector types (`DirectX::PackedVector::*`), e.g. `XMCOLOR`, `XMHALF4`, `XMBYTE4` |

### Extension Header Files (`Extensions/`)

Each extension header provides optimized overrides for a specific instruction set tier. They must be included **after** `DirectXMath.h` and live in their own sub-namespace:

| File | Namespace | Instruction set |
| --- | --- | --- |
| `DirectXMathSSE3.h` | `DirectX::SSE3` | SSE3 |
| `DirectXMathSSE4.h` | `DirectX::SSE4` | SSE4.1 |
| `DirectXMathAVX.h` | `DirectX::AVX` | AVX |
| `DirectXMathAVX2.h` | `DirectX::AVX2` | AVX2 |
| `DirectXMathFMA3.h` | `DirectX::FMA3` | FMA3 |
| `DirectXMathFMA4.h` | `DirectX::FMA4` | FMA4 (AMD) |
| `DirectXMathF16C.h` | `DirectX::F16C` | F16C (half-precision) |
| `DirectXMathBE.h` | `DirectX::BEMath` | Big-endian scalar fallback |

## Namespace Structure

All types and functions are in the `DirectX` namespace. Sub-namespaces are used for optional components:

| Namespace | Contents |
| --- | --- |
| `DirectX` | All core math types and functions |
| `DirectX::PackedVector` | GPU-friendly packed types from `DirectXPackedVector.h` |
| `DirectX::Colors` | sRGB named color constants from `DirectXColors.h` |
| `DirectX::ColorsLinear` | Linear-space color constants from `DirectXColors.h` |
| `DirectX::MathInternal` | Private implementation details (do not use directly) |
| `DirectX::SSE3`, `::SSE4`, `::AVX`, `::AVX2`, `::FMA3`, `::FMA4`, `::F16C`, `::BEMath` | Extension-specific optimized overrides |

## SIMD Intrinsic Macro Hierarchy

DirectXMath auto-selects the SIMD backend based on compiler/architecture. Macros cascade from highest to lowest:

```cpp
_XM_AVX2_INTRINSICS_     -> enables _XM_FMA3_INTRINSICS_ + _XM_F16C_INTRINSICS_
_XM_AVX_INTRINSICS_      -> enables _XM_SSE4_INTRINSICS_
_XM_SSE4_INTRINSICS_     -> enables _XM_SSE3_INTRINSICS_
_XM_SSE3_INTRINSICS_     -> enables _XM_SSE_INTRINSICS_
_XM_SSE_INTRINSICS_      (auto-enabled for x86/x64)
_XM_ARM_NEON_INTRINSICS_ (auto-enabled for ARM/ARM64)
_XM_NO_INTRINSICS_       (pure C++ fallback; force with this define)
```

Optional macros:
- `_XM_SVML_INTRINSICS_` — Intel Short Vector Math Library (auto-enabled with MSVC 2019+; opt-out with `_XM_DISABLE_INTEL_SVML_`)
- `_XM_NO_MOVNT_` — Disables non-temporal store instructions (`_mm_stream_ps`, etc.)
- `_XM_FAVOR_INTEL_` — Uses `_mm_permute_ps` (AVX) over `_mm_shuffle_ps` when available
- `_XM_NO_XMVECTOR_OVERLOADS_` — Disables `XMVECTOR` arithmetic operators (auto-set for GCC/Clang)

When writing multi-path implementations, follow this pattern:

```cpp
#if defined(_XM_NO_INTRINSICS_)
    // Pure C++ path
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // ARM NEON path
#elif defined(_XM_AVX2_INTRINSICS_)
    // AVX2 path
#elif defined(_XM_SSE4_INTRINSICS_)
    // SSE4.1 path
#elif defined(_XM_SSE_INTRINSICS_)
    // SSE/SSE2 path (baseline for x86/x64)
#endif
```

## Core Type System

### Computation Type

`XMVECTOR` is the fundamental 128-bit SIMD register type. It maps to `__m128` (SSE), `float32x4_t` (ARM NEON), or a plain struct (no-intrinsics). It must be **16-byte aligned**.

### Storage Types

Storage types hold data in memory. Use `XMLoad*` / `XMStore*` to move between storage and `XMVECTOR`.

| Type | Description |
| --- | --- |
| `XMFLOAT2` / `XMFLOAT2A` | 2D float vector (A = 16-byte aligned) |
| `XMFLOAT3` / `XMFLOAT3A` | 3D float vector |
| `XMFLOAT4` / `XMFLOAT4A` | 4D float vector |
| `XMINT2/3/4` | 2/3/4D signed int32 vector |
| `XMUINT2/3/4` | 2/3/4D unsigned int32 vector |
| `XMFLOAT3X3`, `XMFLOAT4X3`, `XMFLOAT4X4`, `XMFLOAT4X4A` | Matrix storage types |

### Constant Helper Types

Used for declaring compile-time SIMD constants (not for computation):

```cpp
XMVECTORF32  // float constants:    { { { 1.f, 0.f, 0.f, 1.f } } }
XMVECTORI32  // int32 constants:    { { { 1, 0, 0, 1 } } }
XMVECTORU32  // uint32 constants:   { { { 0xFFFFFFFF, 0, 0, 0 } } }
XMVECTORU8   // uint8 byte pattern: { { { 0x00, 0xFF, ... } } }
```

### Calling Convention Typedefs

XMVECTOR and XMMATRIX have special parameter-passing typedefs to maximize register usage. **Always use these instead of raw `XMVECTOR`/`XMMATRIX` for function parameters.**

| Typedef | Purpose |
| --- | --- |
| `FXMVECTOR` | 1st–3rd `XMVECTOR` parameters (in-register on x86/ARM/ARM64/vectorcall) |
| `GXMVECTOR` | 4th `XMVECTOR` parameter (in-register on ARM/ARM64/vectorcall) |
| `HXMVECTOR` | 5th–6th `XMVECTOR` parameters (in-register on ARM64/vectorcall) |
| `CXMVECTOR` | 7th+ `XMVECTOR` parameters (always by `const` reference) |
| `FXMMATRIX` | 1st `XMMATRIX` parameter (in-register on ARM64/vectorcall) |
| `CXMMATRIX` | 2nd+ `XMMATRIX` parameters (always by `const` reference) |

All functions that take or return `XMVECTOR`/`XMMATRIX` must use `XM_CALLCONV`:

```cpp
XMVECTOR XM_CALLCONV XMVectorHermite(
    FXMVECTOR Position0, FXMVECTOR Tangent0, FXMVECTOR Position1,
    GXMVECTOR Tangent1, float t) noexcept;

XMVECTOR XM_CALLCONV XMVector3Project(
    FXMVECTOR V,
    float ViewportX, float ViewportY, float ViewportWidth, float ViewportHeight,
    float ViewportMinZ, float ViewportMaxZ,
    FXMMATRIX Projection, CXMMATRIX View, CXMMATRIX World) noexcept;
```

### Load/Store Pattern

Always load storage types to `XMVECTOR` before computation, and store back afterward:

```cpp
XMFLOAT3 input{ 1.f, 2.f, 3.f };
XMVECTOR v = XMLoadFloat3(&input);        // storage → register
v = XMVector3Normalize(v);
XMFLOAT3 output;
XMStoreFloat3(&output, v);               // register → storage
```

Use the `A`-suffixed variants (`XMLoadFloat4A`, `XMStoreFloat4A`) only with 16-byte-aligned storage types (`XMFLOAT4A`, etc.).

## Naming Conventions

Source: `wiki/Implementation.md`

- **PascalCase**: class names, methods, free functions, enums
- **camelCase**: struct member variables
- **UPPERCASE**: preprocessor defines and nameless enums
- **`XM` prefix**: types (`XMVECTOR`, `XMFLOAT3`), utility macros (`XM_PI`, `XM_CALLCONV`)
- **`XMVector*/XMMatrix*/XMColor*/XMQuaternion*`**: function families
- **`XMLoad*/XMStore*`**: load/store conversion functions
- No [Hungarian notation](https://wikipedia.org/wiki/Hungarian_notation) except `p` for pointers and `sz` for strings

## SAL Annotations

The library uses SAL2 annotations on all function parameters. Include them on any new code:

```cpp
_In_                       // non-null input pointer
_Out_                      // non-null output pointer (written before read)
_Inout_                    // non-null input and output pointer
_In_reads_(Count)          // input array of Count elements
_Out_writes_(Count)        // output array of Count elements
_In_reads_bytes_(n)        // input pointer reading n bytes
```

Example:

```cpp
XMMATRIX XM_CALLCONV XMLoadFloat4x4(_In_ const XMFLOAT4X4* pSource) noexcept;
void XM_CALLCONV XMStoreFloat4x4(_Out_ XMFLOAT4X4* pDestination, _In_ FXMMATRIX M) noexcept;
```

SAL annotations compile to no-ops unless `/analyze` is used; they never affect code generation.

## Key Macros

| Macro | Description |
| --- | --- |
| `XM_CALLCONV` | `__vectorcall` (MSVC/clang-cl x86/x64), `__fastcall` (fallback), empty (GCC) |
| `XM_ALIGNED_STRUCT(n)` | Portable `alignas(n) struct` / `__declspec(align(n)) struct` |
| `XM_ALIGNED_DATA(n)` | Portable `alignas(n)` / `__attribute__((aligned(n)))` |
| `XM_DEPRECATED` | `[[deprecated]]` (C++14+), `__attribute__((deprecated))` (GCC), `__declspec(deprecated)` (MSVC) |
| `XM_CACHE_LINE_SIZE` | `64` (x86/x64), `128` (ARM/ARM64) |
| `XMGLOBALCONST` | Used for global color constants in `DirectXColors.h` |
| `XM_STREAM_PS` / `XM_SFENCE` | Non-temporal stores (controlled by `_XM_NO_MOVNT_`) |
| `XM_FMADD_PS` / `XM_FNMADD_PS` | FMA ops when `_XM_FMA3_INTRINSICS_` is active, else expanded to mul+add |
| `XM_PERMUTE_PS` | `_mm_permute_ps` (AVX) or `_mm_shuffle_ps` (SSE) |

## CMake Build Options

| CMake Option | Default | Description |
| --- | --- | --- |
| `BUILD_XDSP` | `OFF` | Build XDSP Digital Signal Processing extension |
| `BUILD_SHMATH` | `OFF` | Build Spherical Harmonics math extension |
| `BUILD_TESTING` | (CTest) | Enable CTest unit tests (requires `Tests/` submodule) |
| `BUILD_AVX_TEST` | `OFF` | Test preset: enable AVX instruction set |
| `BUILD_AVX2_TEST` | `OFF` | Test preset: enable AVX2 instruction set |
| `BUILD_F16C_TEST` | `OFF` | Test preset: enable F16C instruction set |
| `BUILD_NO_INTRINSICS` | `OFF` | Test preset: force `_XM_NO_INTRINSICS_` mode |

The CMake library exports as `Microsoft::DirectXMath` and requires **CMake 3.21+**.

## C++ Standard Feature Usage

The library uses conditional C++ standard feature guards:

```cpp
#if (__cplusplus >= 201402L)
    // C++14 features (e.g., [[deprecated]])
#endif

#if (__cplusplus >= 201703L)
    // C++17 features (e.g., alignas for XM_ALIGNED_STRUCT)
#endif

#if (__cplusplus >= 202002L)
    // C++20 features (e.g., spaceship <=> operator on XMFLOAT2/3/4)
    #include <compare>
    bool operator == (const XMFLOAT2&) const = default;
    auto operator <=> (const XMFLOAT2&) const = default;
#endif
```

The base requirement is **C++11**. Do not unconditionally use C++14/17/20 features without guards.

## File Header Convention

Every source file (`.h`, `.inl`, etc.) must begin with this block:

```cpp
//-------------------------------------------------------------------------------------
// {FileName}
//
// {One-line description}
//
// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
//
// https://go.microsoft.com/fwlink/?LinkID=615560
//-------------------------------------------------------------------------------------
```

Section separators within files use:
- Major sections: `//-------------------------------------------------------------------------------------`
- Subsections:   `//---------------------------------------------------------------------------------`

The project does **not** use Doxygen. API documentation is maintained exclusively on the GitHub wiki.
## References

- [Source git repository on GitHub](https://github.com/microsoft/DirectXMath.git)
- [DirectXMath wiki git repository on GitHub](https://github.com/microsoft/DirectXMath.wiki.git)
- [DirectXMath test suite git repository on GitHub](https://github.com/walbourn/directxmathtest.wiki.git).
- [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
- [Microsoft Secure Coding Guidelines](https://learn.microsoft.com/en-us/security/develop/secure-coding-guidelines)
- [CMake Documentation](https://cmake.org/documentation/)
- [VCPKG Documentation](https://learn.microsoft.com/vcpkg/)
- [DirectXMath Documentation](https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal)
- [Games for Windows and the DirectX SDK blog - March 2012](https://walbourn.github.io/introducing-directxmath/)

## No speculation

When creating documentation:

### Document Only What Exists

- Only document features, patterns, and decisions that are explicitly present in the source code.
- Only include configurations and requirements that are clearly specified.
- Do not make assumptions about implementation details.

### Handle Missing Information

- Ask the user questions to gather missing information.
- Document gaps in current implementation or specifications.
- List open questions that need to be addressed.

### Source Material

- Always cite the specific source file and line numbers for documented features.
- Link directly to relevant source code when possible.
- Indicate when information comes from requirements vs. implementation.

### Verification Process

- Review each documented item against source code whenever related to the task.
- Remove any speculative content.
- Ensure all documentation is verifiable against the current state of the codebase.

## Cross-platform Support Notes

- The code supports building for Windows and Linux.
- Portability and conformance of the code is validated by building with Visual C++, clang/LLVM for Windows, MinGW, and GCC for Linux.

### Platform and Compiler `#ifdef` Guards

Use these established guards — do not invent new ones:

| Guard | Purpose |
| --- | --- |
| `_WIN32` | Windows platform (desktop, UWP, Xbox) |
| `_MSC_VER` | MSVC-specific (and MSVC-like clang-cl) pragmas and warning suppression |
| `__clang__` | Clang/LLVM diagnostic suppressions |
| `__MINGW32__` | MinGW compatibility headers |
| `_M_ARM64` / `_M_X64` / `_M_IX86` | Architecture-specific code paths for MSVC (`#ifdef`) |
| `_M_ARM64EC` | ARM64EC ABI (ARM64 code with x64 interop using ARM-NEON) for MSVC |
| `__aarch64__` / `__x86_64__` / `__i386__` | Additional architecture-specific symbols for MinGW/GNUC (`#if`) |

> `_M_ARM`/ `__arm__` is legacy 32-bit ARM which is deprecated.

## Code Review Instructions

When reviewing code, focus on the following aspects:

- Adherence to coding standards defined in `.editorconfig` and on the [wiki](https://github.com/microsoft/DirectXMath/wiki/Implementation).
- Make coding recommendations based on the *C++ Core Guidelines*.
- Proper use of RAII and smart pointers.
- Correct error handling practices and C++ Exception safety.
- Clarity and maintainability of the code.
- Adequate comments where necessary.
- Public interfaces located in `Inc\*.h` should be clearly defined and documented on the Microsoft Docs pages.
- Optional functions are available in headers in the `SHMath`, `Stereo3D`, `MatrixStack`, and `XDSP` folders.
- Compliance with the project's architecture and design patterns.
- Ensure that all public functions and classes are covered by unit tests located on [GitHub](https://github.com/walbourn/directxmathtest.git) where applicable. Report any gaps in test coverage.
- Check for performance implications, especially in geometry processing algorithms.
- Provide brutally honest feedback on code quality, design, and potential improvements as needed.
