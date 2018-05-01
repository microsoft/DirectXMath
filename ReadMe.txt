-----------
DirectXMath
-----------

Copyright (c) Microsoft Corporation. All rights reserved.

February 2018

This package contains the DirectXMath library, an all inline SIMD C++ linear algebra library
for use in games and graphics apps

This code is designed to build with Visual Studio 2015 or 2017. It is recommended that you
make use of the latest updates (VS 2015 Update 3 or VS 2017 15.6 update).

These components are designed to work without requiring any content from the DirectX SDK. For details,
see "Where is the DirectX SDK?" <http://msdn.microsoft.com/en-us/library/ee663275.aspx>.

Inc\
    DirectXMath Files (in the DirectX C++ namespace)
        DirectXMath.h - Core library
        DirectXPackedVector.h - Load/Store functions and types for working with various compressed GPU formats
        DirectXColors.h - .NET-style Color defines in sRGB color space
        DirectXCollision.h - Bounding volume collision library

Extentions\
    Advanced instruction set variants for guarded codepaths
        DirectXMathSSE3.h - SSE3
        DirectXMathBE.h - Supplemental SSE3 (SSSE3)
        DirectXMathSSE4.h - SSE4.1
        DirectXMathAVX.h - Advanced Vector Extensions (AVX)
        DirectXMathAVX2.h - Advanced Vector Extensions 2 (AVX2)
        DirectXMathF16C.h - Half-precision conversions (F16C)
        DirectXMathFMA3.h - Fused multiply-accumulate (FMA3)
        DirectXMathFMA4.h - Fused multiply-accumulate (FMA4)

SHMath\
    Spherical Harmonics math functions
        DirectXSH.h - Header for SHMath functions
        DirectXSH.cpp, DirectXSHD3D11.cpp - Implementation

XDSP\
    XDSP.h - Digital Signal Processing helper functions

All content and source code for this package are subject to the terms of the MIT License.
<http://opensource.org/licenses/MIT>.

Documentation is available at <https://msdn.microsoft.com/en-us/library/windows/desktop/hh437833.aspx>.

For the latest version of DirectXMath, bug reports, etc. please visit the project site.
<https://github.com/Microsoft/DirectXMath>

This project has adopted the Microsoft Open Source Code of Conduct. For more information see the
Code of Conduct FAQ or contact opencode@microsoft.com with any additional questions or comments.

https://opensource.microsoft.com/codeofconduct/


---------------
RELEASE HISTORY
---------------

February 2018 (3.12)
    ARM64 use of fused multiply-accumulate intriniscs
    Conformance fix for XMConvertFloatToHalf
    Minor code cleanup

June 2017 (3.11)
    AVX optimization of XMMatrixMultiply and XMMatrixMultiplyTranspose
    AVX2 optimization for XMVectorSplatX
    FMA3 optimization of XMVectorMultiplyAdd and XMVectorNegativeMultiplySubtract (implied by /arch:AVX2)
    Conformance fixes to support compilation with Clang 3.7

January 2017 (3.10)
    Added XMVectorSum for horizontal adds
    ARMv8 intrinsics use for ARM64 platform (division, rounding, half-precision conversion)
    Added SSE3 codepaths using opt-in _XM_SSE3_INTRINSICS_
    XMVectorRound fix for no-intrinsics to match round to nearest (even)
    XMStoreFloat3SE fix when max channel isn't a perfect power of 2
    constexpr conformance fix and workaround for compiler bug in VS 2015 RTM
    Remove support for VS 2012 compilers
    Remove __vector4i deprecated type

June 2016 (3.09)
    Includes support for additional optimizations when built with /arch:AVX or /arch:AVX2
    Added use of constexpr for type constructors, XMConvertToRadians, and XMConvertToDegrees
    Marked __vector4i, XMXDEC4, XMDECN4, XMDEC4, and associated Load & Store functions as deprecated.
        These are vestiges of Xbox 360 support and will be removed in a future release
    Renamed parameter in XMMatrixPerspectiveFov* to reduce user confusion when relying on IntelliSense
    XMU565, XMUNIBBLE4 constructors take uint8_t instead of int8_t

May 2016
    DirectXMath 3.08 released under the MIT license

November 2015 (3.08)
    Added use of _mm_sfence for Stream methods
    Fixed bug with non-uniform scaling transforms for BoundingOrientedBox
    Added asserts for Near/FarZ in XMMatrix* methods
    Added use of =default for PODs with VS 2013/2015
    Additional SSE and ARM-NEON optimizations for PackedVector functions

April 2015 (3.07)
    Fix customer reported bugs in BoundingBox methods
    Fix customer reported bug in XMStoreFloat3SE  
    Fix customer reported bug in XMVectorATan2, XMVectorATan2Est  
    Fix customer reported bug in XMVectorRound 

October 2013 (3.06)
    Fixed load/store of XMFLOAT3SE to properly match the DXGI_FORMAT_R9G9B9E5_SHAREDEXP
    Added XMLoadUDecN4_XR and XMStoreUDecN4_XR to match DXGI_FORMAT_R10G10B10_XR_BIAS_A2_UNORM
    Added XMColorRGBToSRGB and XMColorSRGBToRGB to convert linear RGB <-> sRGB

July 2013 (3.05)
    Use x86/x64 __vectorcall calling-convention when available (XM_CALLCONV, HXMVECTOR, FXMMATRIX introduced)
    Fixed bug with XMVectorFloor and XMVectorCeiling when given whole odd numbers (i.e. 105.0)
    Improved XMVectorRound algorithm
    ARM-NEON optimizations for XMVectorExp2, XMVectorLog2, XMVectorExpE, and XMVectorLogE  
    ARM-NEON code paths use multiply-by-scalar intrinsics when supported
    Additional optimizations for ARM-NEON Stream functions
    Fixed potential warning C4723 using operator/ or operator/=

March 2013 (3.04)
    XMVectorExp2, XMVectorLog2, XMVectorExpE, and XMVectorLogE functions added to provide base-e support in addition to the existing base-2 support
    XMVectorExp and XMVectorLog are now aliases for XMVectorExp2 and XMVectorLog2  
    Additional optimizations for Stream functions
    XMVector3Cross now ensures w component is zero on ARM
    XMConvertHalfToFloat and XMConvertFloatToHalf  now use IEEE 754 standard float16 behavior for INF/QNAN
    Updated matrix version Transform for  BoundingOrientedBox  and  BoundingFrustum  to handle scaling

March 2012 (3.03)
    Breaking change: Removed union members from XMMATRIX type to make it a fully 'opaque' type
    Marked single-parameter C++ constructors for XMFLOAT2, XMFLOAT2A, XMFLOAT3, XMFLOAT3A, XMFLOAT4, and XMFLOAT4A explicit

February 2012 (3.02)
    ARM-NEON intrinsics (selected by default for the ARM platform)
    reworked XMVectorPermute, change of XM_PERMUTE_ defines, removal of XMVectorPermuteControl
    Addition of XM_SWIZZLE_ defines
    Optimizations for transcendental functions
    Template forms for permute, swizzle, shift-left, rotate-left, rotation-right, and insert
    Removal of deprecated types and functions
        (XM_CACHE_LINE_SIZE define, XMVectorExpEst, XMVectorLogEst, XMVectorPowEst, XMVectorSinHEs, XMVectorCosHEst, XMVectorTanHEst,
         XMVector2InBoundsR, XMVector3InBoundsR, XMVector4InBoundsR)
    Removed XM_STRICT_VECTOR4; XMVECTOR in NO-INTRINSICS always defined without .x, .y, .z, .w, .v, or .u
    Additional bounding types
    SAL fixes and improvements

September 2011 (3.00)
    Renamed and reorganized the headers
    Introduced C++ namespaces
    Removed the Xbox 360-specific GPU types
        (HENDN3, XMHEND3, XMUHENDN3, XMUHEND3, XMDHENN3, XMDHEN3,
         XMUDHENN3, XMUDHEN3, XMXICON4, XMXICO4, XMICON4, XMICO4, XMUICON4, XMUICO4 )
