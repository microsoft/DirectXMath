# DirectXMath API Reference Overview

This document provides an overview of the DirectXMath API surface. All function signatures are defined in the public headers: `DirectXCollision.h`, `DirectXColors.h`, `DirectXMath.h`, `DirectXPackedVector.h`.

## Core Header: `DirectXMath.h`

The primary header defines all fundamental types, constants, and math functions.

### Fundamental Types

| Type | Description |
| --- | --- |
| `XMVECTOR` | 128-bit SIMD register type for computation (16-byte aligned) |
| `XMMATRIX` | 4x4 matrix of four `XMVECTOR` rows (16-byte aligned) |

### Storage Types

| Type | Description |
| --- | --- |
| `XMFLOAT2` / `XMFLOAT2A` | 2D float vector (A = 16-byte aligned) |
| `XMFLOAT3` / `XMFLOAT3A` | 3D float vector |
| `XMFLOAT4` / `XMFLOAT4A` | 4D float vector |
| `XMINT2` / `XMINT3` / `XMINT4` | 2/3/4D signed int32 vector |
| `XMUINT2` / `XMUINT3` / `XMUINT4` | 2/3/4D unsigned int32 vector |
| `XMFLOAT3X3` | 3x3 float matrix |
| `XMFLOAT3X4` / `XMFLOAT3X4A` | 3x4 float matrix |
| `XMFLOAT4X3` / `XMFLOAT4X3A` | 4x3 float matrix |
| `XMFLOAT4X4` / `XMFLOAT4X4A` | 4x4 float matrix |

### Constant Helper Types

| Type | Description |
| --- | --- |
| `XMVECTORF32` | Compile-time float vector constant |
| `XMVECTORI32` | Compile-time int32 vector constant |
| `XMVECTORU32` | Compile-time uint32 vector constant |
| `XMVECTORU8` | Compile-time byte vector constant |

### Calling Convention Typedefs

| Typedef | Usage |
| --- | --- |
| `FXMVECTOR` | 1st-3rd `XMVECTOR` parameters |
| `GXMVECTOR` | 4th `XMVECTOR` parameter |
| `HXMVECTOR` | 5th-6th `XMVECTOR` parameters |
| `CXMVECTOR` | 7th+ `XMVECTOR` parameters |
| `FXMMATRIX` | 1st `XMMATRIX` parameter |
| `CXMMATRIX` | 2nd+ `XMMATRIX` parameters |

### Function Families

#### Load and Store Functions

Convert between storage types and `XMVECTOR`/`XMMATRIX`:

- `XMLoadFloat2`, `XMLoadFloat2A`, `XMLoadFloat3`, `XMLoadFloat3A`, `XMLoadFloat4`, `XMLoadFloat4A`
- `XMLoadInt2`, `XMLoadInt3`, `XMLoadInt4`
- `XMLoadFloat3x3`, `XMLoadFloat3x4`, `XMLoadFloat3x4A`, `XMLoadFloat4x3`, `XMLoadFloat4x3A`, `XMLoadFloat4x4`, `XMLoadFloat4x4A`
- `XMStoreFloat2`, `XMStoreFloat2A`, `XMStoreFloat3`, `XMStoreFloat3A`, `XMStoreFloat4`, `XMStoreFloat4A`
- `XMStoreInt2`, `XMStoreInt3`, `XMStoreInt4`
- `XMStoreFloat3x3`, `XMStoreFloat3x4`, `XMStoreFloat3x4A`, `XMStoreFloat4x3`, `XMStoreFloat4x3A`, `XMStoreFloat4x4`, `XMStoreFloat4x4A`

#### Vector Functions (`XMVector*`)

Organized by component count (2D, 3D, 4D) and general operations:

- **General**: `XMVectorZero`, `XMVectorSet`, `XMVectorReplicate`, `XMVectorSplatX/Y/Z/W`, `XMVectorSwizzle`, `XMVectorPermute`, `XMVectorSelect`
- **Arithmetic**: `XMVectorAdd`, `XMVectorSubtract`, `XMVectorMultiply`, `XMVectorDivide`, `XMVectorScale`, `XMVectorNegate`, `XMVectorMultiplyAdd`, `XMVectorNegativeMultiplySubtract`
- **Comparison**: `XMVectorEqual`, `XMVectorGreater`, `XMVectorLess`, `XMVectorMin`, `XMVectorMax`, `XMVectorClamp`
- **Bitwise**: `XMVectorAndInt`, `XMVectorOrInt`, `XMVectorXorInt`, `XMVectorNotInt`
- **Transcendental**: `XMVectorSin`, `XMVectorCos`, `XMVectorSinCos`, `XMVectorTan`, `XMVectorSinH`, `XMVectorCosH`, `XMVectorATan`, `XMVectorATan2`
- **Exponential/Power**: `XMVectorExp2`, `XMVectorLog2`, `XMVectorPow`, `XMVectorSqrt`, `XMVectorReciprocalSqrt`

#### 2D Vector Functions (`XMVector2*`)

`XMVector2Length`, `XMVector2LengthSq`, `XMVector2Dot`, `XMVector2Cross`, `XMVector2Normalize`, `XMVector2Transform`, `XMVector2TransformCoord`, `XMVector2TransformNormal`

#### 3D Vector Functions (`XMVector3*`)

`XMVector3Length`, `XMVector3LengthSq`, `XMVector3Dot`, `XMVector3Cross`, `XMVector3Normalize`, `XMVector3Transform`, `XMVector3TransformCoord`, `XMVector3TransformNormal`, `XMVector3Project`, `XMVector3Unproject`, `XMVector3Rotate`, `XMVector3InverseRotate`, `XMVector3ComponentsFromNormal`

#### 4D Vector Functions (`XMVector4*`)

`XMVector4Length`, `XMVector4LengthSq`, `XMVector4Dot`, `XMVector4Cross`, `XMVector4Normalize`, `XMVector4Transform`

#### Matrix Functions (`XMMatrix*`)

- **Creation**: `XMMatrixIdentity`, `XMMatrixSet`, `XMMatrixTranslation`, `XMMatrixScaling`, `XMMatrixRotationX/Y/Z`, `XMMatrixRotationAxis`, `XMMatrixRotationQuaternion`, `XMMatrixRotationRollPitchYaw`
- **Operations**: `XMMatrixMultiply`, `XMMatrixTranspose`, `XMMatrixInverse`, `XMMatrixDeterminant`
- **Projection**: `XMMatrixPerspectiveFovLH/RH`, `XMMatrixPerspectiveLH/RH`, `XMMatrixOrthographicLH/RH`, `XMMatrixOrthographicOffCenterLH/RH`
- **View**: `XMMatrixLookAtLH/RH`, `XMMatrixLookToLH/RH`
- **Decomposition**: `XMMatrixDecompose`
- **Transforms**: `XMMatrixTransformation`, `XMMatrixTransformation2D`, `XMMatrixAffineTransformation`, `XMMatrixAffineTransformation2D`, `XMMatrixReflect`, `XMMatrixShadow`

#### Quaternion Functions (`XMQuaternion*`)

`XMQuaternionIdentity`, `XMQuaternionMultiply`, `XMQuaternionConjugate`, `XMQuaternionInverse`, `XMQuaternionNormalize`, `XMQuaternionDot`, `XMQuaternionLength`, `XMQuaternionSlerp`, `XMQuaternionRotationMatrix`, `XMQuaternionRotationAxis`, `XMQuaternionRotationRollPitchYaw`, `XMQuaternionToAxisAngle`

#### Color Functions (`XMColor*`)

`XMColorEqual`, `XMColorNegative`, `XMColorModulate`, `XMColorAdjustSaturation`, `XMColorAdjustContrast`, `XMColorRGBToHSL`, `XMColorHSLToRGB`, `XMColorRGBToHSV`, `XMColorHSVToRGB`, `XMColorRGBToXYZ`, `XMColorXYZToRGB`, `XMColorRGBToSRGB`, `XMColorSRGBToRGB`

#### Plane Functions (`XMPlane*`)

`XMPlaneNormalize`, `XMPlaneDot`, `XMPlaneDotCoord`, `XMPlaneDotNormal`, `XMPlaneFromPointNormal`, `XMPlaneFromPoints`, `XMPlaneIntersectLine`, `XMPlaneIntersectPlane`, `XMPlaneTransform`

### Utility Functions and Constants

- `XMConvertToRadians`, `XMConvertToDegrees`
- `XMScalarSin`, `XMScalarCos`, `XMScalarSinCos`, `XMScalarACos`, `XMScalarASin`
- `XMVerifyCPUSupport`
- Constants: `XM_PI`, `XM_2PI`, `XM_1DIVPI`, `XM_PIDIV2`, `XM_PIDIV4`

---

## Collision Header: `DirectXCollision.h`

Provides bounding volume types and intersection/containment tests.

### Bounding Volume Types

| Type | Description |
| --- | --- |
| `BoundingSphere` | Sphere defined by center and radius |
| `BoundingBox` | Axis-aligned bounding box (AABB) |
| `BoundingOrientedBox` | Oriented bounding box (OBB) |
| `BoundingFrustum` | View frustum |

### Common Methods on Bounding Volumes

Each bounding volume type provides:

- `Transform` — Transform the volume by a matrix or scale/rotation/translation
- `Contains` — Test containment of a point, triangle, sphere, box, or frustum
- `Intersects` — Test intersection with a ray, plane, triangle, sphere, box, or frustum
- `ContainedBy` — Test if contained by a set of planes (frustum culling)
- `CreateFromPoints` — Static factory from a point cloud
- `CreateMerged` — Static factory merging two volumes
- `GetCorners` — Retrieve the 8 corners (box/frustum)

### Intersection Test Functions

- `TriangleTests::Intersects` — Ray-triangle and triangle-triangle intersection
- `TriangleTests::ContainedBy` — Triangle vs. frustum planes

---

## Packed Vector Header: `DirectXPackedVector.h`

Provides GPU-friendly packed types in the `DirectX::PackedVector` namespace. Each type includes `XMLoad*` and `XMStore*` conversion functions.

### Packed Types

| Type | Format |
| --- | --- |
| `XMCOLOR` | BGRA 8-bit per channel (32-bit) |
| `XMHALF2` / `XMHALF4` | 16-bit half-precision float |
| `XMSHORTN2` / `XMSHORTN4` | Signed normalized 16-bit |
| `XMUSHORT2` / `XMUSHORT4` / `XMUSHORTN2` / `XMUSHORTN4` | Unsigned 16-bit |
| `XMBYTE2` / `XMBYTE4` / `XMBYTEN2` / `XMBYTEN4` | Signed 8-bit |
| `XMUBYTE2` / `XMUBYTE4` / `XMUBYTEN2` / `XMUBYTEN4` | Unsigned 8-bit |
| `XMU555` | 5-bit RGB + 1-bit alpha |
| `XMU565` | 5-6-5 RGB |
| `XMFLOAT3PK` | 11-11-10 packed float (R11G11B10_FLOAT) |
| `XMFLOAT3SE` | 9-9-9-5 shared exponent (R9G9B9E5) |
| `XMUNIBBLE4` | 4-4-4-4 unsigned |
| `XMUDECN4` / `XMUDEC4` | 10-10-10-2 unsigned |
| `XMDECN4` / `XMDEC4` | 10-10-10-2 signed (deprecated) |
| `XMXDEC4` / `XMXDECN4` | 10-10-10-2 extended (deprecated) |

---

## Colors Header: `DirectXColors.h`

Predefined color constants as `XMVECTORF32` values:

- `DirectX::Colors::*` — sRGB color constants (e.g., `Colors::Red`, `Colors::CornflowerBlue`)
- `DirectX::ColorsLinear::*` — Linear-space equivalents for use in linear rendering pipelines

> These are also known as ".NET colors"

---

## SHMath (Spherical Harmonics)

**Header**: `SHMath/DirectXSH.h`
**Namespace**: `DirectX`
**Compilation**: Requires compiling `DirectXSH.cpp` (and optionally `DirectXSHD3D11.cpp` / `DirectXSHD3D12.cpp`)

Spherical harmonics (SH) functions for real-time lighting. Supports orders 2 through 6.

### Constants

| Constant | Value |
| --- | --- |
| `XM_SH_MINORDER` | 2 |
| `XM_SH_MAXORDER` | 6 |

### Core SH Functions

| Function | Description |
| --- | --- |
| `XMSHEvalDirection` | Evaluate SH basis functions for a given direction |
| `XMSHRotate` | Rotate SH coefficients by a rotation matrix |
| `XMSHRotateZ` | Rotate SH coefficients around the Z axis |
| `XMSHAdd` | Add two sets of SH coefficients |
| `XMSHScale` | Scale SH coefficients by a scalar |
| `XMSHDot` | Dot product of two SH coefficient sets |
| `XMSHMultiply` | Multiply (convolve) two SH coefficient sets |
| `XMSHMultiply2` through `XMSHMultiply6` | Order-specific optimized multiply |

### Light Evaluation Functions

| Function | Description |
| --- | --- |
| `XMSHEvalDirectionalLight` | Project a directional light into SH coefficients |
| `XMSHEvalSphericalLight` | Project a spherical area light into SH coefficients |
| `XMSHEvalConeLight` | Project a cone light into SH coefficients |
| `XMSHEvalHemisphereLight` | Project a hemisphere light into SH coefficients |

### Cubemap Projection (Optional)

| Function | Description |
| --- | --- |
| `SHProjectCubeMap` (D3D11) | Project a D3D11 cubemap texture into SH (requires `DirectXSHD3D11.cpp`) |
| `SHProjectCubeMap` (D3D12) | Project D3D12 cubemap subresource data into SH (requires `DirectXSHD3D12.cpp`) |

---

## XDSP (Digital Signal Processing)

**Header**: `XDSP/XDSP.h`
**Namespace**: `XDSP`
**Compilation**: Header-only (all functions are `inline`)

Provides DirectXMath-based DSP functions primarily for audio FFT processing. All buffer parameters must be 16-byte aligned. Only single-precision floating-point is supported.

### Primary Functions

| Function | Description |
| --- | --- |
| `FFTInitializeUnityTable` | Initialize twiddle factor (unity) table for FFT |
| `FFT` | In-place radix-4 decimation-in-time FFT |
| `FFTUnswizzle` | Reorder FFT output from bit-reversed to natural order |
| `FFTPolar` | Convert complex FFT results to magnitude (polar form) |
| `FFTInterleaved` | Forward FFT on interleaved (real/imaginary alternating) data |
| `IFFTDeinterleaved` | Inverse FFT producing deinterleaved real/imaginary output |
| `Deinterleave` | Split interleaved real/imaginary data into separate buffers |
| `Interleave` | Combine separate real/imaginary buffers into interleaved format |

### Internal Butterfly Functions

These are used internally by the FFT routines:

| Function | Description |
| --- | --- |
| `vmulComplex` | Parallel multiplication of four complex numbers |
| `ButterflyDIT4_1` | Radix-4 DIT butterfly (single vector) |
| `ButterflyDIT4_4` | Radix-4 DIT butterfly (four parallel vectors) |
| `FFT4` | 4-point FFT |
| `FFT8` | 8-point FFT |
| `FFT16` | 16-point FFT |

### Usage Pattern

```cpp
#include <DirectXMath.h>
#include "XDSP.h"

// Typical FFT workflow:
// 1. Allocate 16-byte aligned buffers for real, imaginary, and unity table
// 2. Initialize unity table with FFTInitializeUnityTable
// 3. Perform FFT with FFT or FFTInterleaved
// 4. Unswizzle output with FFTUnswizzle
// 5. Optionally convert to polar with FFTPolar
```

---

## Extension Headers

Extension headers in `Extensions/` provide instruction-set-specific optimized overrides. Include them **after** `DirectXMath.h` and call functions from their respective namespaces:

| Header | Namespace | Instruction Set |
| --- | --- | --- |
| `DirectXMathSSE3.h` | `DirectX::SSE3` | SSE3 |
| `DirectXMathSSE4.h` | `DirectX::SSE4` | SSE4.1 |
| `DirectXMathAVX.h` | `DirectX::AVX` | AVX |
| `DirectXMathAVX2.h` | `DirectX::AVX2` | AVX2 |
| `DirectXMathFMA3.h` | `DirectX::FMA3` | FMA3 |
| `DirectXMathFMA4.h` | `DirectX::FMA4` | FMA4 (AMD) |
| `DirectXMathF16C.h` | `DirectX::F16C` | F16C (half-precision) |
| `DirectXMathBE.h` | `DirectX::BEMath` | Big-endian scalar fallback |

> **Note**: Extensions are not needed when building with `/arch:AVX` or `/arch:AVX2`, as the core library automatically uses the additional instructions.

---

## Further Reading

- [DirectXMath Programming Guide (Microsoft Learn)](https://learn.microsoft.com/windows/win32/dxmath/ovw-xnamath-progguide)
- [DirectXMath Reference (Microsoft Learn)](https://learn.microsoft.com/windows/win32/dxmath/ovw-xnamath-reference)
- [GitHub Wiki](https://github.com/microsoft/DirectXMath/wiki/)
