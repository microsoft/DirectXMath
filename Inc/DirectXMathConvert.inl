//-------------------------------------------------------------------------------------
// DirectXMathConvert.inl -- SIMD C++ Math library
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//  
// Copyright (c) Microsoft Corporation. All rights reserved.
//-------------------------------------------------------------------------------------

#ifdef _MSC_VER
#pragma once
#endif

/****************************************************************************
 *
 * Data conversion
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_)
// For VMX128, these routines are all defines in the main header

#pragma warning(push)
#pragma warning(disable:4701) // Prevent warnings about 'Result' potentially being used without having been initialized

inline XMVECTOR XMConvertVectorIntToFloat
(
    FXMVECTOR VInt,
    uint32_t     DivExponent
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t ElementIndex;
    float fScale;
    XMVECTOR Result;
    assert(DivExponent<32);
    fScale = 1.0f / (float)(1U << DivExponent);
    ElementIndex = 0;
    do {
        int32_t iTemp = (int32_t)VInt.vector4_u32[ElementIndex];
        Result.vector4_f32[ElementIndex] = ((float)iTemp) * fScale;
    } while (++ElementIndex<4);
    return Result;
#else // _XM_SSE_INTRINSICS_
    assert(DivExponent<32);
    // Convert to floats
    XMVECTOR vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&VInt)[0]);
    // Convert DivExponent into 1.0f/(1<<DivExponent)
    uint32_t uScale = 0x3F800000U - (DivExponent << 23);
    // Splat the scalar value
    __m128i vScale = _mm_set1_epi32(uScale);
    vResult = _mm_mul_ps(vResult,reinterpret_cast<const __m128 *>(&vScale)[0]);
    return vResult;
#endif
}

//------------------------------------------------------------------------------

inline XMVECTOR XMConvertVectorFloatToInt
(
    FXMVECTOR VFloat,
    uint32_t     MulExponent
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t ElementIndex;
    XMVECTOR Result;
    float fScale;
    assert(MulExponent<32);
    // Get the scalar factor.
    fScale = (float)(1U << MulExponent);
    ElementIndex = 0;
    do {
        int32_t iResult;
        float fTemp = VFloat.vector4_f32[ElementIndex]*fScale;
        if (fTemp <= -(65536.0f*32768.0f)) {
            iResult = (-0x7FFFFFFF)-1;
        } else if (fTemp > (65536.0f*32768.0f)-128.0f) {
            iResult = 0x7FFFFFFF;
        } else {
            iResult = (int32_t)fTemp;
        }
        Result.vector4_u32[ElementIndex] = (uint32_t)iResult;
    } while (++ElementIndex<4);
    return Result;
#else // _XM_SSE_INTRINSICS_
    assert(MulExponent<32);
    XMVECTOR vResult = _mm_set_ps1((float)(1U << MulExponent));
    vResult = _mm_mul_ps(vResult,VFloat);
    // In case of positive overflow, detect it
    XMVECTOR vOverflow = _mm_cmpgt_ps(vResult,g_XMMaxInt);
    // Float to int conversion
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // If there was positive overflow, set to 0x7FFFFFFF
    vResult = _mm_and_ps(vOverflow,g_XMAbsMask);
    vOverflow = _mm_andnot_ps(vOverflow,reinterpret_cast<const __m128 *>(&vResulti)[0]);
    vOverflow = _mm_or_ps(vOverflow,vResult);
    return vOverflow;
#endif
}

//------------------------------------------------------------------------------

inline XMVECTOR XMConvertVectorUIntToFloat
(
    FXMVECTOR VUInt,
    uint32_t      DivExponent
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t ElementIndex;
    float fScale;
    XMVECTOR Result;
    assert(DivExponent<32);
    fScale = 1.0f / (float)(1U << DivExponent);
    ElementIndex = 0;
    do {
        Result.vector4_f32[ElementIndex] = (float)VUInt.vector4_u32[ElementIndex] * fScale;
    } while (++ElementIndex<4);
    return Result;
#else // _XM_SSE_INTRINSICS_
    assert(DivExponent<32);
    // For the values that are higher than 0x7FFFFFFF, a fixup is needed
    // Determine which ones need the fix.
    XMVECTOR vMask = _mm_and_ps(VUInt,g_XMNegativeZero);
    // Force all values positive
    XMVECTOR vResult = _mm_xor_ps(VUInt,vMask);
    // Convert to floats
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Convert 0x80000000 -> 0xFFFFFFFF
    __m128i iMask = _mm_srai_epi32(reinterpret_cast<const __m128i *>(&vMask)[0],31);
    // For only the ones that are too big, add the fixup
    vMask = _mm_and_ps(reinterpret_cast<const __m128 *>(&iMask)[0],g_XMFixUnsigned);
    vResult = _mm_add_ps(vResult,vMask);
    // Convert DivExponent into 1.0f/(1<<DivExponent)
    uint32_t uScale = 0x3F800000U - (DivExponent << 23);
    // Splat
    iMask = _mm_set1_epi32(uScale);
    vResult = _mm_mul_ps(vResult,reinterpret_cast<const __m128 *>(&iMask)[0]);
    return vResult;
#endif
}

//------------------------------------------------------------------------------

inline XMVECTOR XMConvertVectorFloatToUInt
(
    FXMVECTOR VFloat,
    uint32_t      MulExponent
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t ElementIndex;
    XMVECTOR Result;
    float fScale;
    assert(MulExponent<32);
    // Get the scalar factor.
    fScale = (float)(1U << MulExponent);
    ElementIndex = 0;
    do {
        uint32_t uResult;
        float fTemp = VFloat.vector4_f32[ElementIndex]*fScale;
        if (fTemp <= 0.0f) {
            uResult = 0;
        } else if (fTemp >= (65536.0f*65536.0f)) {
            uResult = 0xFFFFFFFFU;
        } else {
            uResult = (uint32_t)fTemp;
        }
        Result.vector4_u32[ElementIndex] = uResult;
    } while (++ElementIndex<4);
    return Result;
#else // _XM_SSE_INTRINSICS_
    assert(MulExponent<32);
    XMVECTOR vResult = _mm_set_ps1(static_cast<float>(1U << MulExponent));
    vResult = _mm_mul_ps(vResult,VFloat);
    // Clamp to >=0
    vResult = _mm_max_ps(vResult,g_XMZero);
    // Any numbers that are too big, set to 0xFFFFFFFFU
    XMVECTOR vOverflow = _mm_cmpgt_ps(vResult,g_XMMaxUInt);
    XMVECTOR vValue = g_XMUnsignedFix;
    // Too large for a signed integer?
    XMVECTOR vMask = _mm_cmpge_ps(vResult,vValue);
    // Zero for number's lower than 0x80000000, 32768.0f*65536.0f otherwise
    vValue = _mm_and_ps(vValue,vMask);
    // Perform fixup only on numbers too large (Keeps low bit precision)
    vResult = _mm_sub_ps(vResult,vValue);
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Convert from signed to unsigned pnly if greater than 0x80000000
    vMask = _mm_and_ps(vMask,g_XMNegativeZero);
    vResult = _mm_xor_ps(reinterpret_cast<const __m128 *>(&vResulti)[0],vMask);
    // On those that are too large, set to 0xFFFFFFFF
    vResult = _mm_or_ps(vResult,vOverflow);
    return vResult;
#endif
}

#pragma warning(pop)

#endif // _XM_NO_INTRINSICS_ || _XM_SSE_INTRINSICS_

/****************************************************************************
 *
 * Vector and matrix load operations
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt(const uint32_t* pSource)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;
    assert(pSource);
    assert(((uintptr_t)pSource & 3) == 0);

    V.vector4_u32[0] = *pSource;
    V.vector4_u32[1] = 0;
    V.vector4_u32[2] = 0;
    V.vector4_u32[3] = 0;

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    assert(((uintptr_t)pSource & 3) == 0);

    return _mm_load_ss( (const float*)pSource );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat(const float* pSource)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;
    assert(pSource);
    assert(((uintptr_t)pSource & 3) == 0);

    V.vector4_f32[0] = *pSource;
    V.vector4_f32[1] = 0.f;
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    assert(((uintptr_t)pSource & 3) == 0);

    return _mm_load_ss( pSource );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt2
(
    const uint32_t* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_u32[0] = pSource[0];
    V.vector4_u32[1] = pSource[1];
    V.vector4_u32[2] = 0;
    V.vector4_u32[3] = 0;
    return V;
#elif defined(_XM_SSE_INTRINSICS_)

    assert(pSource);

    __m128 x = _mm_load_ss( (const float*)pSource );
    __m128 y = _mm_load_ss( (const float*)(pSource+1) );
    return _mm_unpacklo_ps( x, y );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt2A
(
    const uint32_t* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    V.vector4_u32[0] = pSource[0];
    V.vector4_u32[1] = pSource[1];
    V.vector4_u32[2] = 0;
    V.vector4_u32[3] = 0;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    __m128i V = _mm_loadl_epi64( (const __m128i*)pSource );
    return reinterpret_cast<__m128 *>(&V)[0];

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat2
(
    const XMFLOAT2* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    ((uint32_t *)(&V.vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->x))[0];
    ((uint32_t *)(&V.vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->y))[0];
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( &pSource->x );
    __m128 y = _mm_load_ss( &pSource->y );
    return _mm_unpacklo_ps( x, y );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat2A
(
    const XMFLOAT2A* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    V.vector4_f32[0] = pSource->x;
    V.vector4_f32[1] = pSource->y;
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    __m128i V = _mm_loadl_epi64( (const __m128i*)pSource );
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadSInt2
(
    const XMINT2* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x;
    V.vector4_f32[1] = (float)pSource->y;
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( (const float*)&pSource->x );
    __m128 y = _mm_load_ss( (const float*)&pSource->y );
    __m128 V = _mm_unpacklo_ps( x, y );
    return _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&V)[0]);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadUInt2
(
    const XMUINT2* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x;
    V.vector4_f32[1] = (float)pSource->y;
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( (const float*)&pSource->x );
    __m128 y = _mm_load_ss( (const float*)&pSource->y );
    __m128 V = _mm_unpacklo_ps( x, y );
    // For the values that are higher than 0x7FFFFFFF, a fixup is needed
    // Determine which ones need the fix.
    XMVECTOR vMask = _mm_and_ps(V,g_XMNegativeZero);
    // Force all values positive
    XMVECTOR vResult = _mm_xor_ps(V,vMask);
    // Convert to floats
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Convert 0x80000000 -> 0xFFFFFFFF
    __m128i iMask = _mm_srai_epi32(reinterpret_cast<const __m128i *>(&vMask)[0],31);
    // For only the ones that are too big, add the fixup
    vMask = _mm_and_ps(reinterpret_cast<const __m128 *>(&iMask)[0],g_XMFixUnsigned);
    vResult = _mm_add_ps(vResult,vMask);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt3
(
    const uint32_t* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_u32[0] = pSource[0];
    V.vector4_u32[1] = pSource[1];
    V.vector4_u32[2] = pSource[2];
    V.vector4_u32[3] = 0;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( (const float*)pSource );
    __m128 y = _mm_load_ss( (const float*)(pSource+1) );
    __m128 z = _mm_load_ss( (const float*)(pSource+2) );
    __m128 xy = _mm_unpacklo_ps( x, y );
    return _mm_movelh_ps( xy, z );

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt3A
(
    const uint32_t* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    V.vector4_u32[0] = pSource[0];
    V.vector4_u32[1] = pSource[1];
    V.vector4_u32[2] = pSource[2];
    V.vector4_u32[3] = 0;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    // Reads an extra integer
    __m128i V = _mm_load_si128( (const __m128i*)pSource );

    // TODO - Zero out .w

    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat3
(
    const XMFLOAT3* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    ((uint32_t *)(&V.vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->x))[0];
    ((uint32_t *)(&V.vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->y))[0];
    ((uint32_t *)(&V.vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->z))[0];
    V.vector4_f32[3] = 0.f;
    return V;
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( &pSource->x );
    __m128 y = _mm_load_ss( &pSource->y );
    __m128 z = _mm_load_ss( &pSource->z );
    __m128 xy = _mm_unpacklo_ps( x, y );
    return _mm_movelh_ps( xy, z );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat3A
(
    const XMFLOAT3A* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    V.vector4_f32[0] = pSource->x;
    V.vector4_f32[1] = pSource->y;
    V.vector4_f32[2] = pSource->z;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    // This reads 1 floats past the memory that should be ignored.
    return _mm_load_ps( &pSource->x );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadSInt3
(
    const XMINT3* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

#ifdef _XBOX_VER
    V = XMLoadInt3( (const uint32_t*)pSource );
    return XMConvertVectorIntToFloat( V, 0 );
#else
    V.vector4_f32[0] = (float)pSource->x;
    V.vector4_f32[1] = (float)pSource->y;
    V.vector4_f32[2] = (float)pSource->z;
    V.vector4_f32[3] = 0.f;
    return V;
#endif

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( (const float*)&pSource->x );
    __m128 y = _mm_load_ss( (const float*)&pSource->y );
    __m128 z = _mm_load_ss( (const float*)&pSource->z );
    __m128 xy = _mm_unpacklo_ps( x, y );
    __m128 V = _mm_movelh_ps( xy, z );
    return _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&V)[0]);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadUInt3
(
    const XMUINT3* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x;
    V.vector4_f32[1] = (float)pSource->y;
    V.vector4_f32[2] = (float)pSource->z;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    __m128 x = _mm_load_ss( (const float*)&pSource->x );
    __m128 y = _mm_load_ss( (const float*)&pSource->y );
    __m128 z = _mm_load_ss( (const float*)&pSource->z );
    __m128 xy = _mm_unpacklo_ps( x, y );
    __m128 V = _mm_movelh_ps( xy, z );
    // For the values that are higher than 0x7FFFFFFF, a fixup is needed
    // Determine which ones need the fix.
    XMVECTOR vMask = _mm_and_ps(V,g_XMNegativeZero);
    // Force all values positive
    XMVECTOR vResult = _mm_xor_ps(V,vMask);
    // Convert to floats
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Convert 0x80000000 -> 0xFFFFFFFF
    __m128i iMask = _mm_srai_epi32(reinterpret_cast<const __m128i *>(&vMask)[0],31);
    // For only the ones that are too big, add the fixup
    vMask = _mm_and_ps(reinterpret_cast<const __m128 *>(&iMask)[0],g_XMFixUnsigned);
    vResult = _mm_add_ps(vResult,vMask);
    return vResult; 

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt4
(
    const uint32_t* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_u32[0] = pSource[0];
    V.vector4_u32[1] = pSource[1];
    V.vector4_u32[2] = pSource[2];
    V.vector4_u32[3] = pSource[3];

    return V;

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pSource);

    __m128i V = _mm_loadu_si128( (const __m128i*)pSource );
    return reinterpret_cast<__m128 *>(&V)[0];

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadInt4A
(
    const uint32_t* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    V.vector4_u32[0] = pSource[0];
    V.vector4_u32[1] = pSource[1];
    V.vector4_u32[2] = pSource[2];
    V.vector4_u32[3] = pSource[3];
    return V;

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    __m128i V = _mm_load_si128( (const __m128i*)pSource );
    return reinterpret_cast<__m128 *>(&V)[0];

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat4
(
    const XMFLOAT4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    ((uint32_t *)(&V.vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->x))[0];
    ((uint32_t *)(&V.vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->y))[0];
    ((uint32_t *)(&V.vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->z))[0];
    ((uint32_t *)(&V.vector4_f32[3]))[0] = ((const uint32_t *)(&pSource->w))[0];
    return V;
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);

    return _mm_loadu_ps( &pSource->x );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadFloat4A
(
    const XMFLOAT4A* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    V.vector4_f32[0] = pSource->x;
    V.vector4_f32[1] = pSource->y;
    V.vector4_f32[2] = pSource->z;
    V.vector4_f32[3] = pSource->w;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    return _mm_load_ps( &pSource->x );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadSInt4
(
    const XMINT4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

#ifdef _XBOX_VER
    V = XMLoadInt4( (const uint32_t*)pSource );
    return XMConvertVectorIntToFloat( V, 0 );
#else
    V.vector4_f32[0] = (float)pSource->x;
    V.vector4_f32[1] = (float)pSource->y;
    V.vector4_f32[2] = (float)pSource->z;
    V.vector4_f32[3] = (float)pSource->w;
    return V;
#endif

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    __m128i V = _mm_loadu_si128( (const __m128i*)pSource );
    return _mm_cvtepi32_ps(V);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMLoadUInt4
(
    const XMUINT4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x;
    V.vector4_f32[1] = (float)pSource->y;
    V.vector4_f32[2] = (float)pSource->z;
    V.vector4_f32[3] = (float)pSource->w;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    __m128i V = _mm_loadu_si128( (const __m128i*)pSource );
    // For the values that are higher than 0x7FFFFFFF, a fixup is needed
    // Determine which ones need the fix.
    XMVECTOR vMask = _mm_and_ps(reinterpret_cast<const __m128 *>(&V)[0],g_XMNegativeZero);
    // Force all values positive
    XMVECTOR vResult = _mm_xor_ps(reinterpret_cast<const __m128 *>(&V)[0],vMask);
    // Convert to floats
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Convert 0x80000000 -> 0xFFFFFFFF
    __m128i iMask = _mm_srai_epi32(reinterpret_cast<const __m128i *>(&vMask)[0],31);
    // For only the ones that are too big, add the fixup
    vMask = _mm_and_ps(reinterpret_cast<const __m128 *>(&iMask)[0],g_XMFixUnsigned);
    vResult = _mm_add_ps(vResult,vMask);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMMATRIX XMLoadFloat3x3
(
    const XMFLOAT3X3* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;

    assert(pSource);

    M.r[0].vector4_f32[0] = pSource->m[0][0];
    M.r[0].vector4_f32[1] = pSource->m[0][1];
    M.r[0].vector4_f32[2] = pSource->m[0][2];
    M.r[0].vector4_f32[3] = 0.0f;

    M.r[1].vector4_f32[0] = pSource->m[1][0];
    M.r[1].vector4_f32[1] = pSource->m[1][1];
    M.r[1].vector4_f32[2] = pSource->m[1][2];
    M.r[1].vector4_f32[3] = 0.0f;

    M.r[2].vector4_f32[0] = pSource->m[2][0];
    M.r[2].vector4_f32[1] = pSource->m[2][1];
    M.r[2].vector4_f32[2] = pSource->m[2][2];
    M.r[2].vector4_f32[3] = 0.0f;

    M.r[3].vector4_f32[0] = 0.0f;
    M.r[3].vector4_f32[1] = 0.0f;
    M.r[3].vector4_f32[2] = 0.0f;
    M.r[3].vector4_f32[3] = 1.0f;

    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX M;
    XMVECTOR V1, V2, V3, Z, T1, T2, T3, T4, T5;

    Z = _mm_setzero_ps();

    assert(pSource);

    V1 = _mm_loadu_ps( &pSource->m[0][0] );
    V2 = _mm_loadu_ps( &pSource->m[1][1] );
    V3 = _mm_load_ss( &pSource->m[2][2] );

    T1 = _mm_unpackhi_ps( V1, Z );
    T2 = _mm_unpacklo_ps( V2, Z );
    T3 = _mm_shuffle_ps( V3, T2, _MM_SHUFFLE( 0, 1, 0, 0 ) );
    T4 = _mm_movehl_ps( T2, T3 );
    T5 = _mm_movehl_ps( Z, T1 );  

    M.r[0] = _mm_movelh_ps( V1, T1 );
    M.r[1] = _mm_add_ps( T4, T5 );
    M.r[2] = _mm_shuffle_ps( V2, V3, _MM_SHUFFLE(1, 0, 3, 2) );
    M.r[3] = g_XMIdentityR3;

    return M;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMMATRIX XMLoadFloat4x3
(
    const XMFLOAT4X3* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMMATRIX M;
    assert(pSource);

    ((uint32_t *)(&M.r[0].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[0][0]))[0];
    ((uint32_t *)(&M.r[0].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[0][1]))[0];
    ((uint32_t *)(&M.r[0].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[0][2]))[0];
    M.r[0].vector4_f32[3] = 0.0f;

    ((uint32_t *)(&M.r[1].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[1][0]))[0];
    ((uint32_t *)(&M.r[1].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[1][1]))[0];
    ((uint32_t *)(&M.r[1].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[1][2]))[0];
    M.r[1].vector4_f32[3] = 0.0f;

    ((uint32_t *)(&M.r[2].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[2][0]))[0];
    ((uint32_t *)(&M.r[2].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[2][1]))[0];
    ((uint32_t *)(&M.r[2].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[2][2]))[0];
    M.r[2].vector4_f32[3] = 0.0f;

    ((uint32_t *)(&M.r[3].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[3][0]))[0];
    ((uint32_t *)(&M.r[3].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[3][1]))[0];
    ((uint32_t *)(&M.r[3].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[3][2]))[0];
    M.r[3].vector4_f32[3] = 1.0f;

    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Use unaligned load instructions to 
    // load the 12 floats
    // vTemp1 = x1,y1,z1,x2
    XMVECTOR vTemp1 = _mm_loadu_ps(&pSource->m[0][0]);
    // vTemp2 = y2,z2,x3,y3
    XMVECTOR vTemp2 = _mm_loadu_ps(&pSource->m[1][1]);
    // vTemp4 = z3,x4,y4,z4
    XMVECTOR vTemp4 = _mm_loadu_ps(&pSource->m[2][2]);
    // vTemp3 = x3,y3,z3,z3
    XMVECTOR vTemp3 = _mm_shuffle_ps(vTemp2,vTemp4,_MM_SHUFFLE(0,0,3,2));
    // vTemp2 = y2,z2,x2,x2
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp1,_MM_SHUFFLE(3,3,1,0));
    // vTemp2 = x2,y2,z2,z2
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp2,_MM_SHUFFLE(1,1,0,2));
    // vTemp1 = x1,y1,z1,0
    vTemp1 = _mm_and_ps(vTemp1,g_XMMask3);
    // vTemp2 = x2,y2,z2,0
    vTemp2 = _mm_and_ps(vTemp2,g_XMMask3);
    // vTemp3 = x3,y3,z3,0
    vTemp3 = _mm_and_ps(vTemp3,g_XMMask3);
    // vTemp4i = x4,y4,z4,0
    __m128i vTemp4i = _mm_srli_si128(reinterpret_cast<const __m128i *>(&vTemp4)[0],32/8);
    // vTemp4i = x4,y4,z4,1.0f
    vTemp4i = _mm_or_si128(vTemp4i,g_XMIdentityR3);
    XMMATRIX M(vTemp1,
            vTemp2,
            vTemp3,
            reinterpret_cast<const __m128 *>(&vTemp4i)[0]);
    return M;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMMATRIX XMLoadFloat4x3A
(
    const XMFLOAT4X3A* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    M.r[0].vector4_f32[0] = pSource->m[0][0];
    M.r[0].vector4_f32[1] = pSource->m[0][1];
    M.r[0].vector4_f32[2] = pSource->m[0][2];
    M.r[0].vector4_f32[3] = 0.0f;

    M.r[1].vector4_f32[0] = pSource->m[1][0];
    M.r[1].vector4_f32[1] = pSource->m[1][1];
    M.r[1].vector4_f32[2] = pSource->m[1][2];
    M.r[1].vector4_f32[3] = 0.0f;

    M.r[2].vector4_f32[0] = pSource->m[2][0];
    M.r[2].vector4_f32[1] = pSource->m[2][1];
    M.r[2].vector4_f32[2] = pSource->m[2][2];
    M.r[2].vector4_f32[3] = 0.0f;

    M.r[3].vector4_f32[0] = pSource->m[3][0];
    M.r[3].vector4_f32[1] = pSource->m[3][1];
    M.r[3].vector4_f32[2] = pSource->m[3][2];
    M.r[3].vector4_f32[3] = 1.0f;

    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Use aligned load instructions to 
    // load the 12 floats
    // vTemp1 = x1,y1,z1,x2
    XMVECTOR vTemp1 = _mm_load_ps(&pSource->m[0][0]);
    // vTemp2 = y2,z2,x3,y3
    XMVECTOR vTemp2 = _mm_load_ps(&pSource->m[1][1]);
    // vTemp4 = z3,x4,y4,z4
    XMVECTOR vTemp4 = _mm_load_ps(&pSource->m[2][2]);
    // vTemp3 = x3,y3,z3,z3
    XMVECTOR vTemp3 = _mm_shuffle_ps(vTemp2,vTemp4,_MM_SHUFFLE(0,0,3,2));
    // vTemp2 = y2,z2,x2,x2
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp1,_MM_SHUFFLE(3,3,1,0));
    // vTemp2 = x2,y2,z2,z2
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp2,_MM_SHUFFLE(1,1,0,2));
    // vTemp1 = x1,y1,z1,0
    vTemp1 = _mm_and_ps(vTemp1,g_XMMask3);
    // vTemp2 = x2,y2,z2,0
    vTemp2 = _mm_and_ps(vTemp2,g_XMMask3);
    // vTemp3 = x3,y3,z3,0
    vTemp3 = _mm_and_ps(vTemp3,g_XMMask3);
    // vTemp4i = x4,y4,z4,0
    __m128i vTemp4i = _mm_srli_si128(reinterpret_cast<const __m128i *>(&vTemp4)[0],32/8);
    // vTemp4i = x4,y4,z4,1.0f
    vTemp4i = _mm_or_si128(vTemp4i,g_XMIdentityR3);
    XMMATRIX M(vTemp1,
            vTemp2,
            vTemp3,
            reinterpret_cast<const __m128 *>(&vTemp4i)[0]);
    return M;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMMATRIX XMLoadFloat4x4
(
    const XMFLOAT4X4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMMATRIX M;
    assert(pSource);

    ((uint32_t *)(&M.r[0].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[0][0]))[0];
    ((uint32_t *)(&M.r[0].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[0][1]))[0];
    ((uint32_t *)(&M.r[0].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[0][2]))[0];
    ((uint32_t *)(&M.r[0].vector4_f32[3]))[0] = ((const uint32_t *)(&pSource->m[0][3]))[0];

    ((uint32_t *)(&M.r[1].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[1][0]))[0];
    ((uint32_t *)(&M.r[1].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[1][1]))[0];
    ((uint32_t *)(&M.r[1].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[1][2]))[0];
    ((uint32_t *)(&M.r[1].vector4_f32[3]))[0] = ((const uint32_t *)(&pSource->m[1][3]))[0];

    ((uint32_t *)(&M.r[2].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[2][0]))[0];
    ((uint32_t *)(&M.r[2].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[2][1]))[0];
    ((uint32_t *)(&M.r[2].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[2][2]))[0];
    ((uint32_t *)(&M.r[2].vector4_f32[3]))[0] = ((const uint32_t *)(&pSource->m[2][3]))[0];

    ((uint32_t *)(&M.r[3].vector4_f32[0]))[0] = ((const uint32_t *)(&pSource->m[3][0]))[0];
    ((uint32_t *)(&M.r[3].vector4_f32[1]))[0] = ((const uint32_t *)(&pSource->m[3][1]))[0];
    ((uint32_t *)(&M.r[3].vector4_f32[2]))[0] = ((const uint32_t *)(&pSource->m[3][2]))[0];
    ((uint32_t *)(&M.r[3].vector4_f32[3]))[0] = ((const uint32_t *)(&pSource->m[3][3]))[0];

    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    XMMATRIX M;

    M.r[0] = _mm_loadu_ps( &pSource->_11 );
    M.r[1] = _mm_loadu_ps( &pSource->_21 );
    M.r[2] = _mm_loadu_ps( &pSource->_31 );
    M.r[3] = _mm_loadu_ps( &pSource->_41 );

    return M;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMMATRIX XMLoadFloat4x4A
(
    const XMFLOAT4X4A* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;

    assert(pSource);
    assert(((uintptr_t)pSource & 0xF) == 0);

    M.r[0].vector4_f32[0] = pSource->m[0][0];
    M.r[0].vector4_f32[1] = pSource->m[0][1];
    M.r[0].vector4_f32[2] = pSource->m[0][2];
    M.r[0].vector4_f32[3] = pSource->m[0][3];

    M.r[1].vector4_f32[0] = pSource->m[1][0];
    M.r[1].vector4_f32[1] = pSource->m[1][1];
    M.r[1].vector4_f32[2] = pSource->m[1][2];
    M.r[1].vector4_f32[3] = pSource->m[1][3];

    M.r[2].vector4_f32[0] = pSource->m[2][0];
    M.r[2].vector4_f32[1] = pSource->m[2][1];
    M.r[2].vector4_f32[2] = pSource->m[2][2];
    M.r[2].vector4_f32[3] = pSource->m[2][3];

    M.r[3].vector4_f32[0] = pSource->m[3][0];
    M.r[3].vector4_f32[1] = pSource->m[3][1];
    M.r[3].vector4_f32[2] = pSource->m[3][2];
    M.r[3].vector4_f32[3] = pSource->m[3][3];

    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX M;

    assert(pSource);

    M.r[0] = _mm_load_ps( &pSource->_11 );
    M.r[1] = _mm_load_ps( &pSource->_21 );
    M.r[2] = _mm_load_ps( &pSource->_31 );
    M.r[3] = _mm_load_ps( &pSource->_41 );

    return M;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * Vector and matrix store operations
 *
 ****************************************************************************/

inline void XMStoreInt
(
    uint32_t*    pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    *pDestination = XMVectorGetIntX( V );

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    _mm_store_ss( (float*)pDestination, V );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat
(
    float*    pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    *pDestination = XMVectorGetX( V );

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    _mm_store_ss( pDestination, V );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreInt2
(
    uint32_t*    pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination[0] = V.vector4_u32[0];
    pDestination[1] = V.vector4_u32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    XMVECTOR T = _mm_shuffle_ps( V, V, _MM_SHUFFLE( 1, 1, 1, 1 ) );
    _mm_store_ss( (float*)&pDestination[0], V );
    _mm_store_ss( (float*)&pDestination[1], T );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreInt2A
(
    uint32_t*    pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination[0] = V.vector4_u32[0];
    pDestination[1] = V.vector4_u32[1];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    _mm_storel_epi64( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&V)[0] );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat2
(
    XMFLOAT2* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination->x = V.vector4_f32[0];
    pDestination->y = V.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    XMVECTOR T = _mm_shuffle_ps( V, V, _MM_SHUFFLE( 1, 1, 1, 1 ) );
    _mm_store_ss( &pDestination->x, V );
    _mm_store_ss( &pDestination->y, T );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat2A
(
    XMFLOAT2A*   pDestination, 
    FXMVECTOR     V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination->x = V.vector4_f32[0];
    pDestination->y = V.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    _mm_storel_epi64( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&V)[0] );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreSInt2
(
    XMINT2* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination->x = (int32_t)V.vector4_f32[0];
    pDestination->y = (int32_t)V.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    // In case of positive overflow, detect it
    XMVECTOR vOverflow = _mm_cmpgt_ps(V,g_XMMaxInt);
    // Float to int conversion
    __m128i vResulti = _mm_cvttps_epi32(V);
    // If there was positive overflow, set to 0x7FFFFFFF
    XMVECTOR vResult = _mm_and_ps(vOverflow,g_XMAbsMask);
    vOverflow = _mm_andnot_ps(vOverflow,reinterpret_cast<const __m128 *>(&vResulti)[0]);
    vOverflow = _mm_or_ps(vOverflow,vResult);
    // Write two ints
    XMVECTOR T = _mm_shuffle_ps( vOverflow, vOverflow, _MM_SHUFFLE( 1, 1, 1, 1 ) );
    _mm_store_ss( (float*)&pDestination->x, vOverflow );
    _mm_store_ss( (float*)&pDestination->y, T );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreUInt2
(
    XMUINT2* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination->x = (uint32_t)V.vector4_f32[0];
    pDestination->y = (uint32_t)V.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    // Clamp to >=0
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    // Any numbers that are too big, set to 0xFFFFFFFFU
    XMVECTOR vOverflow = _mm_cmpgt_ps(vResult,g_XMMaxUInt);
    XMVECTOR vValue = g_XMUnsignedFix;
    // Too large for a signed integer?
    XMVECTOR vMask = _mm_cmpge_ps(vResult,vValue);
    // Zero for number's lower than 0x80000000, 32768.0f*65536.0f otherwise
    vValue = _mm_and_ps(vValue,vMask);
    // Perform fixup only on numbers too large (Keeps low bit precision)
    vResult = _mm_sub_ps(vResult,vValue);
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Convert from signed to unsigned pnly if greater than 0x80000000
    vMask = _mm_and_ps(vMask,g_XMNegativeZero);
    vResult = _mm_xor_ps(reinterpret_cast<const __m128 *>(&vResulti)[0],vMask);
    // On those that are too large, set to 0xFFFFFFFF
    vResult = _mm_or_ps(vResult,vOverflow);
    // Write two uints
    XMVECTOR T = _mm_shuffle_ps( vResult, vResult, _MM_SHUFFLE( 1, 1, 1, 1 ) );
    _mm_store_ss( (float*)&pDestination->x, vResult );
    _mm_store_ss( (float*)&pDestination->y, T );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreInt3
(
    uint32_t*    pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination[0] = V.vector4_u32[0];
    pDestination[1] = V.vector4_u32[1];
    pDestination[2] = V.vector4_u32[2];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    XMVECTOR T1 = _mm_shuffle_ps(V,V,_MM_SHUFFLE(1,1,1,1));
    XMVECTOR T2 = _mm_shuffle_ps(V,V,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss( (float*)pDestination, V );
    _mm_store_ss( (float*)&pDestination[1], T1 );
    _mm_store_ss( (float*)&pDestination[2], T2 );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreInt3A
(
    uint32_t*    pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination[0] = V.vector4_u32[0];
    pDestination[1] = V.vector4_u32[1];
    pDestination[2] = V.vector4_u32[2];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    XMVECTOR T = _mm_shuffle_ps(V,V,_MM_SHUFFLE(2,2,2,2));
    _mm_storel_epi64( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&V)[0] );
    _mm_store_ss( (float*)&pDestination[2], T );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat3
(
    XMFLOAT3* pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination->x = V.vector4_f32[0];
    pDestination->y = V.vector4_f32[1];
    pDestination->z = V.vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    XMVECTOR T1 = _mm_shuffle_ps(V,V,_MM_SHUFFLE(1,1,1,1));
    XMVECTOR T2 = _mm_shuffle_ps(V,V,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss( &pDestination->x, V );
    _mm_store_ss( &pDestination->y, T1 );
    _mm_store_ss( &pDestination->z, T2 );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat3A
(
    XMFLOAT3A*   pDestination, 
    FXMVECTOR     V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination->x = V.vector4_f32[0];
    pDestination->y = V.vector4_f32[1];
    pDestination->z = V.vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    XMVECTOR T = _mm_shuffle_ps(V,V,_MM_SHUFFLE(2,2,2,2));
    _mm_storel_epi64( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&V)[0] );
    _mm_store_ss( &pDestination->z, T );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreSInt3
(
    XMINT3* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination->x = (int32_t)V.vector4_f32[0];
    pDestination->y = (int32_t)V.vector4_f32[1];
    pDestination->z = (int32_t)V.vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    // In case of positive overflow, detect it
    XMVECTOR vOverflow = _mm_cmpgt_ps(V,g_XMMaxInt);
    // Float to int conversion
    __m128i vResulti = _mm_cvttps_epi32(V);
    // If there was positive overflow, set to 0x7FFFFFFF
    XMVECTOR vResult = _mm_and_ps(vOverflow,g_XMAbsMask);
    vOverflow = _mm_andnot_ps(vOverflow,reinterpret_cast<const __m128 *>(&vResulti)[0]);
    vOverflow = _mm_or_ps(vOverflow,vResult);
    // Write 3 uints
    XMVECTOR T1 = _mm_shuffle_ps(vOverflow,vOverflow,_MM_SHUFFLE(1,1,1,1));
    XMVECTOR T2 = _mm_shuffle_ps(vOverflow,vOverflow,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss( (float*)&pDestination->x, vOverflow );
    _mm_store_ss( (float*)&pDestination->y, T1 );
    _mm_store_ss( (float*)&pDestination->z, T2 );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreUInt3
(
    XMUINT3* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    pDestination->x = (uint32_t)V.vector4_f32[0];
    pDestination->y = (uint32_t)V.vector4_f32[1];
    pDestination->z = (uint32_t)V.vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 3) == 0);

    // Clamp to >=0
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    // Any numbers that are too big, set to 0xFFFFFFFFU
    XMVECTOR vOverflow = _mm_cmpgt_ps(vResult,g_XMMaxUInt);
    XMVECTOR vValue = g_XMUnsignedFix;
    // Too large for a signed integer?
    XMVECTOR vMask = _mm_cmpge_ps(vResult,vValue);
    // Zero for number's lower than 0x80000000, 32768.0f*65536.0f otherwise
    vValue = _mm_and_ps(vValue,vMask);
    // Perform fixup only on numbers too large (Keeps low bit precision)
    vResult = _mm_sub_ps(vResult,vValue);
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Convert from signed to unsigned pnly if greater than 0x80000000
    vMask = _mm_and_ps(vMask,g_XMNegativeZero);
    vResult = _mm_xor_ps(reinterpret_cast<const __m128 *>(&vResulti)[0],vMask);
    // On those that are too large, set to 0xFFFFFFFF
    vResult = _mm_or_ps(vResult,vOverflow);
    // Write 3 uints
    XMVECTOR T1 = _mm_shuffle_ps(vResult,vResult,_MM_SHUFFLE(1,1,1,1));
    XMVECTOR T2 = _mm_shuffle_ps(vResult,vResult,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss( (float*)&pDestination->x, vResult );
    _mm_store_ss( (float*)&pDestination->y, T1 );
    _mm_store_ss( (float*)&pDestination->z, T2 );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreInt4
(
    uint32_t*    pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination[0] = V.vector4_u32[0];
    pDestination[1] = V.vector4_u32[1];
    pDestination[2] = V.vector4_u32[2];
    pDestination[3] = V.vector4_u32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    
    _mm_storeu_si128( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&V)[0] );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreInt4A
(
    uint32_t*    pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination[0] = V.vector4_u32[0];
    pDestination[1] = V.vector4_u32[1];
    pDestination[2] = V.vector4_u32[2];
    pDestination[3] = V.vector4_u32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    _mm_store_si128( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&V)[0] );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline void XMStoreFloat4
(
    XMFLOAT4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    
    pDestination->x = V.vector4_f32[0];
    pDestination->y = V.vector4_f32[1];
    pDestination->z = V.vector4_f32[2];
    pDestination->w = V.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    
    _mm_storeu_ps( &pDestination->x, V );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat4A
(
    XMFLOAT4A*   pDestination, 
    FXMVECTOR     V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination->x = V.vector4_f32[0];
    pDestination->y = V.vector4_f32[1];
    pDestination->z = V.vector4_f32[2];
    pDestination->w = V.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    _mm_store_ps( &pDestination->x, V );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline void XMStoreSInt4
(
    XMINT4* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination->x = (int32_t)V.vector4_f32[0];
    pDestination->y = (int32_t)V.vector4_f32[1];
    pDestination->z = (int32_t)V.vector4_f32[2];
    pDestination->w = (int32_t)V.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    
    // In case of positive overflow, detect it
    XMVECTOR vOverflow = _mm_cmpgt_ps(V,g_XMMaxInt);
    // Float to int conversion
    __m128i vResulti = _mm_cvttps_epi32(V);
    // If there was positive overflow, set to 0x7FFFFFFF
    XMVECTOR vResult = _mm_and_ps(vOverflow,g_XMAbsMask);
    vOverflow = _mm_andnot_ps(vOverflow,reinterpret_cast<const __m128 *>(&vResulti)[0]);
    vOverflow = _mm_or_ps(vOverflow,vResult);
    _mm_storeu_si128( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&vOverflow)[0] );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreUInt4
(
    XMUINT4* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination->x = (uint32_t)V.vector4_f32[0];
    pDestination->y = (uint32_t)V.vector4_f32[1];
    pDestination->z = (uint32_t)V.vector4_f32[2];
    pDestination->w = (uint32_t)V.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);

    // Clamp to >=0
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    // Any numbers that are too big, set to 0xFFFFFFFFU
    XMVECTOR vOverflow = _mm_cmpgt_ps(vResult,g_XMMaxUInt);
    XMVECTOR vValue = g_XMUnsignedFix;
    // Too large for a signed integer?
    XMVECTOR vMask = _mm_cmpge_ps(vResult,vValue);
    // Zero for number's lower than 0x80000000, 32768.0f*65536.0f otherwise
    vValue = _mm_and_ps(vValue,vMask);
    // Perform fixup only on numbers too large (Keeps low bit precision)
    vResult = _mm_sub_ps(vResult,vValue);
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Convert from signed to unsigned pnly if greater than 0x80000000
    vMask = _mm_and_ps(vMask,g_XMNegativeZero);
    vResult = _mm_xor_ps(reinterpret_cast<const __m128 *>(&vResulti)[0],vMask);
    // On those that are too large, set to 0xFFFFFFFF
    vResult = _mm_or_ps(vResult,vOverflow);
    _mm_storeu_si128( (__m128i*)pDestination, reinterpret_cast<const __m128i *>(&vResult)[0] );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat3x3
(
    XMFLOAT3X3*	pDestination, 
    CXMMATRIX	M
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination->m[0][0] = M.r[0].vector4_f32[0];
    pDestination->m[0][1] = M.r[0].vector4_f32[1];
    pDestination->m[0][2] = M.r[0].vector4_f32[2];

    pDestination->m[1][0] = M.r[1].vector4_f32[0];
    pDestination->m[1][1] = M.r[1].vector4_f32[1];
    pDestination->m[1][2] = M.r[1].vector4_f32[2];

    pDestination->m[2][0] = M.r[2].vector4_f32[0];
    pDestination->m[2][1] = M.r[2].vector4_f32[1];
    pDestination->m[2][2] = M.r[2].vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    XMVECTOR vTemp1 = M.r[0];
    XMVECTOR vTemp2 = M.r[1];
    XMVECTOR vTemp3 = M.r[2];
    XMVECTOR vWork = _mm_shuffle_ps(vTemp1,vTemp2,_MM_SHUFFLE(0,0,2,2));
    vTemp1 = _mm_shuffle_ps(vTemp1,vWork,_MM_SHUFFLE(2,0,1,0));
    _mm_storeu_ps(&pDestination->m[0][0],vTemp1);
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp3,_MM_SHUFFLE(1,0,2,1));
    _mm_storeu_ps(&pDestination->m[1][1],vTemp2);
    vTemp3 = _mm_shuffle_ps(vTemp3,vTemp3,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss(&pDestination->m[2][2],vTemp3);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline void XMStoreFloat4x3
(
    XMFLOAT4X3* pDestination, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination->m[0][0] = M.r[0].vector4_f32[0];
    pDestination->m[0][1] = M.r[0].vector4_f32[1];
    pDestination->m[0][2] = M.r[0].vector4_f32[2];

    pDestination->m[1][0] = M.r[1].vector4_f32[0];
    pDestination->m[1][1] = M.r[1].vector4_f32[1];
    pDestination->m[1][2] = M.r[1].vector4_f32[2];

    pDestination->m[2][0] = M.r[2].vector4_f32[0];
    pDestination->m[2][1] = M.r[2].vector4_f32[1];
    pDestination->m[2][2] = M.r[2].vector4_f32[2];

    pDestination->m[3][0] = M.r[3].vector4_f32[0];
    pDestination->m[3][1] = M.r[3].vector4_f32[1];
    pDestination->m[3][2] = M.r[3].vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    XMVECTOR vTemp1 = M.r[0];
    XMVECTOR vTemp2 = M.r[1];
    XMVECTOR vTemp3 = M.r[2];
    XMVECTOR vTemp4 = M.r[3];
    XMVECTOR vTemp2x = _mm_shuffle_ps(vTemp2,vTemp3,_MM_SHUFFLE(1,0,2,1));
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp1,_MM_SHUFFLE(2,2,0,0));
    vTemp1 = _mm_shuffle_ps(vTemp1,vTemp2,_MM_SHUFFLE(0,2,1,0));
    vTemp3 = _mm_shuffle_ps(vTemp3,vTemp4,_MM_SHUFFLE(0,0,2,2));
    vTemp3 = _mm_shuffle_ps(vTemp3,vTemp4,_MM_SHUFFLE(2,1,2,0));
    _mm_storeu_ps(&pDestination->m[0][0],vTemp1);
    _mm_storeu_ps(&pDestination->m[1][1],vTemp2x);
    _mm_storeu_ps(&pDestination->m[2][2],vTemp3);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat4x3A
(
    XMFLOAT4X3A*	pDestination, 
    CXMMATRIX		M
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination->m[0][0] = M.r[0].vector4_f32[0];
    pDestination->m[0][1] = M.r[0].vector4_f32[1];
    pDestination->m[0][2] = M.r[0].vector4_f32[2];

    pDestination->m[1][0] = M.r[1].vector4_f32[0];
    pDestination->m[1][1] = M.r[1].vector4_f32[1];
    pDestination->m[1][2] = M.r[1].vector4_f32[2];

    pDestination->m[2][0] = M.r[2].vector4_f32[0];
    pDestination->m[2][1] = M.r[2].vector4_f32[1];
    pDestination->m[2][2] = M.r[2].vector4_f32[2];

    pDestination->m[3][0] = M.r[3].vector4_f32[0];
    pDestination->m[3][1] = M.r[3].vector4_f32[1];
    pDestination->m[3][2] = M.r[3].vector4_f32[2];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);
    // x1,y1,z1,w1
    XMVECTOR vTemp1 = M.r[0];
    // x2,y2,z2,w2
    XMVECTOR vTemp2 = M.r[1];
    // x3,y3,z3,w3
    XMVECTOR vTemp3 = M.r[2];
    // x4,y4,z4,w4
    XMVECTOR vTemp4 = M.r[3];
    // z1,z1,x2,y2
    XMVECTOR vTemp = _mm_shuffle_ps(vTemp1,vTemp2,_MM_SHUFFLE(1,0,2,2));
    // y2,z2,x3,y3 (Final)
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp3,_MM_SHUFFLE(1,0,2,1));
    // x1,y1,z1,x2 (Final)
    vTemp1 = _mm_shuffle_ps(vTemp1,vTemp,_MM_SHUFFLE(2,0,1,0));
    // z3,z3,x4,x4
    vTemp3 = _mm_shuffle_ps(vTemp3,vTemp4,_MM_SHUFFLE(0,0,2,2));
    // z3,x4,y4,z4 (Final)
    vTemp3 = _mm_shuffle_ps(vTemp3,vTemp4,_MM_SHUFFLE(2,1,2,0));
    // Store in 3 operations
    _mm_store_ps(&pDestination->m[0][0],vTemp1);
    _mm_store_ps(&pDestination->m[1][1],vTemp2);
    _mm_store_ps(&pDestination->m[2][2],vTemp3);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline void XMStoreFloat4x4
(
    XMFLOAT4X4* pDestination, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination->m[0][0] = M.r[0].vector4_f32[0];
    pDestination->m[0][1] = M.r[0].vector4_f32[1];
    pDestination->m[0][2] = M.r[0].vector4_f32[2];
    pDestination->m[0][3] = M.r[0].vector4_f32[3];

    pDestination->m[1][0] = M.r[1].vector4_f32[0];
    pDestination->m[1][1] = M.r[1].vector4_f32[1];
    pDestination->m[1][2] = M.r[1].vector4_f32[2];
    pDestination->m[1][3] = M.r[1].vector4_f32[3];

    pDestination->m[2][0] = M.r[2].vector4_f32[0];
    pDestination->m[2][1] = M.r[2].vector4_f32[1];
    pDestination->m[2][2] = M.r[2].vector4_f32[2];
    pDestination->m[2][3] = M.r[2].vector4_f32[3];

    pDestination->m[3][0] = M.r[3].vector4_f32[0];
    pDestination->m[3][1] = M.r[3].vector4_f32[1];
    pDestination->m[3][2] = M.r[3].vector4_f32[2];
    pDestination->m[3][3] = M.r[3].vector4_f32[3];
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);

    _mm_storeu_ps( &pDestination->_11, M.r[0] );
    _mm_storeu_ps( &pDestination->_21, M.r[1] );
    _mm_storeu_ps( &pDestination->_31, M.r[2] );
    _mm_storeu_ps( &pDestination->_41, M.r[3] );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void XMStoreFloat4x4A
(
    XMFLOAT4X4A*	pDestination, 
    CXMMATRIX		M
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);
    assert(((uintptr_t)pDestination & 0xF) == 0);

    pDestination->m[0][0] = M.r[0].vector4_f32[0];
    pDestination->m[0][1] = M.r[0].vector4_f32[1];
    pDestination->m[0][2] = M.r[0].vector4_f32[2];
    pDestination->m[0][3] = M.r[0].vector4_f32[3];

    pDestination->m[1][0] = M.r[1].vector4_f32[0];
    pDestination->m[1][1] = M.r[1].vector4_f32[1];
    pDestination->m[1][2] = M.r[1].vector4_f32[2];
    pDestination->m[1][3] = M.r[1].vector4_f32[3];

    pDestination->m[2][0] = M.r[2].vector4_f32[0];
    pDestination->m[2][1] = M.r[2].vector4_f32[1];
    pDestination->m[2][2] = M.r[2].vector4_f32[2];
    pDestination->m[2][3] = M.r[2].vector4_f32[3];

    pDestination->m[3][0] = M.r[3].vector4_f32[0];
    pDestination->m[3][1] = M.r[3].vector4_f32[1];
    pDestination->m[3][2] = M.r[3].vector4_f32[2];
    pDestination->m[3][3] = M.r[3].vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);

    _mm_store_ps( &pDestination->_11, M.r[0] );
    _mm_store_ps( &pDestination->_21, M.r[1] );
    _mm_store_ps( &pDestination->_31, M.r[2] );
    _mm_store_ps( &pDestination->_41, M.r[3] );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

