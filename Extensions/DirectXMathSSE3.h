//-------------------------------------------------------------------------------------
// DirectXMathSSE3.h -- SSE3 extensions for SIMD C++ Math library
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//  
// Copyright (c) Microsoft Corporation. All rights reserved.
//
// http://go.microsoft.com/fwlink/?LinkID=615560
//-------------------------------------------------------------------------------------

#ifdef _MSC_VER
#pragma once
#endif

#ifdef _M_ARM
#error SSE3 not supported on ARM platform
#endif

#pragma warning(push)
#pragma warning(disable : 4987)
#include <intrin.h>
#pragma warning(pop)

#include <pmmintrin.h>

#include <DirectXMath.h>

namespace DirectX
{
#if (DIRECTXMATH_VERSION < 305) && !defined(XM_CALLCONV)
#define XM_CALLCONV __fastcall
typedef const DirectX::XMVECTOR& HXMVECTOR;
typedef const DirectX::XMMATRIX& FXMMATRIX;
#endif

namespace SSE3
{

inline bool XMVerifySSE3Support()
{
    // Should return true on AMD Athlon 64, AMD Phenom, and Intel Pentium 4 or later processors

    // See http://msdn.microsoft.com/en-us/library/hskdteyh.aspx
    int CPUInfo[4] = {-1};
    __cpuid( CPUInfo, 0 );

    if ( CPUInfo[0] < 1  )
        return false;

    __cpuid(CPUInfo, 1 );

    // We only check for SSE3 instruction set. SSSE3 instructions are not used.
    return ( (CPUInfo[2] & 0x1) != 0 );
}

inline XMVECTOR XM_CALLCONV XMVector2Dot
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    XMVECTOR vTemp = _mm_mul_ps(V1,V2);
    vTemp = _mm_hadd_ps(vTemp,vTemp);
    return _mm_shuffle_ps(vTemp,vTemp,_MM_SHUFFLE(0,0,0,0));
}

inline XMVECTOR XM_CALLCONV XMVector2LengthSq( FXMVECTOR V )
{
    return SSE3::XMVector2Dot(V, V);
}

inline XMVECTOR XM_CALLCONV XMVector3Dot
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    XMVECTOR vTemp = _mm_mul_ps(V1,V2);
    vTemp = _mm_and_ps( vTemp, g_XMMask3 );
    vTemp = _mm_hadd_ps(vTemp,vTemp);
    return _mm_hadd_ps(vTemp,vTemp);
}

inline XMVECTOR XM_CALLCONV XMVector3LengthSq( FXMVECTOR V )
{
    return SSE3::XMVector3Dot(V, V);
}

inline XMVECTOR XM_CALLCONV XMVector4Dot
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    XMVECTOR vTemp = _mm_mul_ps(V1,V2);
    vTemp = _mm_hadd_ps( vTemp, vTemp );
    return _mm_hadd_ps( vTemp, vTemp );
}

inline XMVECTOR XM_CALLCONV XMVector4LengthSq( FXMVECTOR V )
{
    return SSE3::XMVector4Dot(V, V);
}

inline XMVECTOR XM_CALLCONV XMVectorSwizzle_0022( FXMVECTOR V )
{
    return _mm_moveldup_ps(V);
}

inline XMVECTOR XM_CALLCONV XMVectorSwizzle_1133( FXMVECTOR V )
{
    return _mm_movehdup_ps(V);
}

}; // namespace SSE3

}; // namespace DirectX;