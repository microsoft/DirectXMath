//-------------------------------------------------------------------------------------
// DirectXMathVector.inl -- SIMD C++ Math library
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

#if defined(_XM_NO_INTRINSICS_)
#define XMISNAN(x)  ((*(uint32_t*)&(x) & 0x7F800000) == 0x7F800000 && (*(uint32_t*)&(x) & 0x7FFFFF) != 0)
#define XMISINF(x)  ((*(uint32_t*)&(x) & 0x7FFFFFFF) == 0x7F800000)
#endif

/****************************************************************************
 *
 * General Vector
 *
 ****************************************************************************/

//------------------------------------------------------------------------------
// Assignment operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Return a vector with all elements equaling zero
inline XMVECTOR XMVectorZero()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {0.0f,0.0f,0.0f,0.0f};
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32(0);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_setzero_ps();
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with four floating point values
inline XMVECTOR XMVectorSet
(
    float x, 
    float y, 
    float z, 
    float w
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = {x,y,z,w};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 V0 = vcreate_f32(((uint64_t)*(const uint32_t *)&x) | ((uint64_t)(*(const uint32_t *)&y) << 32));
    __n64 V1 = vcreate_f32(((uint64_t)*(const uint32_t *)&z) | ((uint64_t)(*(const uint32_t *)&w) << 32));
    return vcombine_f32(V0, V1);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_set_ps( w, z, y, x );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with four integer values
inline XMVECTOR XMVectorSetInt
(
    uint32_t x, 
    uint32_t y, 
    uint32_t z, 
    uint32_t w
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORU32 vResult = {x,y,z,w};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 V0 = vcreate_u32(((uint64_t)x) | ((uint64_t)y << 32));
    __n64 V1 = vcreate_u32(((uint64_t)z) | ((uint64_t)w << 32));
    return vcombine_u32(V0, V1);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_set_epi32( w, z, y, x );
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with a replicated floating point value
inline XMVECTOR XMVectorReplicate
(
    float Value
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
    XMVECTORF32 vResult = {Value,Value,Value,Value};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_f32( Value );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_set_ps1( Value );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with a replicated floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorReplicatePtr
(
    const float *pValue
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
    float Value = pValue[0];
    XMVECTORF32 vResult = {Value,Value,Value,Value};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_dup_f32( pValue );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_load_ps1( pValue );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with a replicated integer value
inline XMVECTOR XMVectorReplicateInt
(
    uint32_t Value
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
    XMVECTORU32 vResult = {Value,Value,Value,Value};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32( Value );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_set1_epi32( Value );
    return _mm_castsi128_ps(vTemp);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with a replicated integer value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorReplicateIntPtr
(
    const uint32_t *pValue
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
    uint32_t Value = pValue[0];
    XMVECTORU32 vResult = {Value,Value,Value,Value};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_dup_u32(pValue);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_load_ps1(reinterpret_cast<const float *>(pValue));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with all bits set (true mask)
inline XMVECTOR XMVectorTrueInt()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORU32 vResult = {0xFFFFFFFFU,0xFFFFFFFFU,0xFFFFFFFFU,0xFFFFFFFFU};
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_s32(-1);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_set1_epi32(-1);
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Initialize a vector with all bits clear (false mask)
inline XMVECTOR XMVectorFalseInt()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {0.0f,0.0f,0.0f,0.0f};
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32(0);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_setzero_ps();
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Replicate the x component of the vector
inline XMVECTOR XMVectorSplatX
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_f32[0] = 
    vResult.vector4_f32[1] = 
    vResult.vector4_f32[2] = 
    vResult.vector4_f32[3] = V.vector4_f32[0];
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32( vget_low_f32( V ), 0 );
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS( V, _MM_SHUFFLE(0, 0, 0, 0) );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Replicate the y component of the vector
inline XMVECTOR XMVectorSplatY
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_f32[0] = 
    vResult.vector4_f32[1] = 
    vResult.vector4_f32[2] = 
    vResult.vector4_f32[3] = V.vector4_f32[1];
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32( vget_low_f32( V ), 1 );
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS( V, _MM_SHUFFLE(1, 1, 1, 1) );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Replicate the z component of the vector
inline XMVECTOR XMVectorSplatZ
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_f32[0] = 
    vResult.vector4_f32[1] = 
    vResult.vector4_f32[2] = 
    vResult.vector4_f32[3] = V.vector4_f32[2];
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32( vget_high_f32( V ), 0 );
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS( V, _MM_SHUFFLE(2, 2, 2, 2) );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Replicate the w component of the vector
inline XMVECTOR XMVectorSplatW
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_f32[0] = 
    vResult.vector4_f32[1] = 
    vResult.vector4_f32[2] = 
    vResult.vector4_f32[3] = V.vector4_f32[3];
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32( vget_high_f32( V ), 1 );
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS( V, _MM_SHUFFLE(3, 3, 3, 3) );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return a vector of 1.0f,1.0f,1.0f,1.0f
inline XMVECTOR XMVectorSplatOne()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_f32[0] = 
    vResult.vector4_f32[1] = 
    vResult.vector4_f32[2] = 
    vResult.vector4_f32[3] = 1.0f;
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_f32(1.0f);
#elif defined(_XM_SSE_INTRINSICS_)
    return g_XMOne;
#else //  _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return a vector of INF,INF,INF,INF
inline XMVECTOR XMVectorSplatInfinity()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_u32[0] = 
    vResult.vector4_u32[1] = 
    vResult.vector4_u32[2] = 
    vResult.vector4_u32[3] = 0x7F800000;
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32(0x7F800000);
#elif defined(_XM_SSE_INTRINSICS_)
    return g_XMInfinity;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return a vector of Q_NAN,Q_NAN,Q_NAN,Q_NAN
inline XMVECTOR XMVectorSplatQNaN()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_u32[0] = 
    vResult.vector4_u32[1] = 
    vResult.vector4_u32[2] = 
    vResult.vector4_u32[3] = 0x7FC00000;
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32(0x7FC00000);
#elif defined(_XM_SSE_INTRINSICS_)
    return g_XMQNaN;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return a vector of 1.192092896e-7f,1.192092896e-7f,1.192092896e-7f,1.192092896e-7f
inline XMVECTOR XMVectorSplatEpsilon()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_u32[0] = 
    vResult.vector4_u32[1] = 
    vResult.vector4_u32[2] = 
    vResult.vector4_u32[3] = 0x34000000;
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32(0x34000000);
#elif defined(_XM_SSE_INTRINSICS_)
    return g_XMEpsilon;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return a vector of -0.0f (0x80000000),-0.0f,-0.0f,-0.0f
inline XMVECTOR XMVectorSplatSignMask()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult;
    vResult.vector4_u32[0] = 
    vResult.vector4_u32[1] = 
    vResult.vector4_u32[2] = 
    vResult.vector4_u32[3] = 0x80000000U;
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_n_u32(0x80000000U);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_set1_epi32( 0x80000000 );
    return reinterpret_cast<__m128*>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return a floating point value via an index. This is not a recommended
// function to use due to performance loss.
inline float XMVectorGetByIndex(FXMVECTOR V, size_t i)
{
    assert( i < 4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_f32[i];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return V.n128_f32[i];
#elif defined(_XM_SSE_INTRINSICS_)
    return V.m128_f32[i];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return the X component in an FPU register. 
inline float XMVectorGetX(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_f32[0];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_f32(V, 0);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cvtss_f32(V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Return the Y component in an FPU register. 
inline float XMVectorGetY(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_f32[1];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_f32(V, 1);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    return _mm_cvtss_f32(vTemp);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Return the Z component in an FPU register. 
inline float XMVectorGetZ(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_f32[2];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_f32(V, 2);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,2,2,2));
    return _mm_cvtss_f32(vTemp);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Return the W component in an FPU register. 
inline float XMVectorGetW(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_f32[3];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_f32(V, 3);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,3,3,3));
    return _mm_cvtss_f32(vTemp);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Store a component indexed by i into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XMVectorGetByIndexPtr(float *f, FXMVECTOR V, size_t i)
{
    assert( f != nullptr );
    assert( i <  4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    *f = V.vector4_f32[i];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    *f = V.n128_f32[i];
#elif defined(_XM_SSE_INTRINSICS_)
    *f = V.m128_f32[i];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Store the X component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XMVectorGetXPtr(float *x, FXMVECTOR V)
{
    assert( x != nullptr);
#if defined(_XM_NO_INTRINSICS_)
    *x = V.vector4_f32[0];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_f32(x,V,0);
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_store_ss(x,V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Store the Y component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XMVectorGetYPtr(float *y, FXMVECTOR V)
{
    assert( y != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *y = V.vector4_f32[1];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_f32(y,V,1);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    _mm_store_ss(y,vResult);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Store the Z component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XMVectorGetZPtr(float *z, FXMVECTOR V)
{
    assert( z != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *z = V.vector4_f32[2];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_f32(z,V,2);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss(z,vResult);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Store the W component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XMVectorGetWPtr(float *w, FXMVECTOR V)
{
    assert( w != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *w = V.vector4_f32[3];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_f32(w,V,3);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,3,3,3));
    _mm_store_ss(w,vResult);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Return an integer value via an index. This is not a recommended
// function to use due to performance loss.
inline uint32_t XMVectorGetIntByIndex(FXMVECTOR V, size_t i)
{
    assert( i < 4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_u32[i];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return V.n128_u32[i];
#elif defined(_XM_SSE_INTRINSICS_)
    return V.m128_u32[i];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Return the X component in an integer register. 
inline uint32_t XMVectorGetIntX(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_u32[0];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_u32(V, 0);
#elif defined(_XM_SSE_INTRINSICS_)
    return static_cast<uint32_t>(_mm_cvtsi128_si32(_mm_castps_si128(V)));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Return the Y component in an integer register. 
inline uint32_t XMVectorGetIntY(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_u32[1];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_u32(V, 1);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vResulti = _mm_shuffle_epi32(_mm_castps_si128(V),_MM_SHUFFLE(1,1,1,1));
    return static_cast<uint32_t>(_mm_cvtsi128_si32(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Return the Z component in an integer register. 
inline uint32_t XMVectorGetIntZ(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_u32[2];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_u32(V, 2);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vResulti = _mm_shuffle_epi32(_mm_castps_si128(V),_MM_SHUFFLE(2,2,2,2));
    return static_cast<uint32_t>(_mm_cvtsi128_si32(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Return the W component in an integer register. 
inline uint32_t XMVectorGetIntW(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V.vector4_u32[3];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vgetq_lane_u32(V, 3);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vResulti = _mm_shuffle_epi32(_mm_castps_si128(V),_MM_SHUFFLE(3,3,3,3));
    return static_cast<uint32_t>(_mm_cvtsi128_si32(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Store a component indexed by i into a 32 bit integer location in memory.
_Use_decl_annotations_
inline void XMVectorGetIntByIndexPtr(uint32_t *x, FXMVECTOR V, size_t i)
{
    assert( x != nullptr );
    assert( i <  4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    *x = V.vector4_u32[i];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    *x = V.n128_u32[i];
#elif defined(_XM_SSE_INTRINSICS_)
    *x = V.m128_u32[i];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Store the X component into a 32 bit integer location in memory.
_Use_decl_annotations_
inline void XMVectorGetIntXPtr(uint32_t *x, FXMVECTOR V)
{
    assert( x != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *x = V.vector4_u32[0];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_u32(x,V,0);
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_store_ss(reinterpret_cast<float *>(x),V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Store the Y component into a 32 bit integer location in memory.
_Use_decl_annotations_
inline void XMVectorGetIntYPtr(uint32_t *y, FXMVECTOR V)
{
    assert( y != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *y = V.vector4_u32[1];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_u32(y,V,1);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    _mm_store_ss(reinterpret_cast<float *>(y),vResult);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Store the Z component into a 32 bit integer locaCantion in memory.
_Use_decl_annotations_
inline void XMVectorGetIntZPtr(uint32_t *z, FXMVECTOR V)
{
    assert( z != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *z = V.vector4_u32[2];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_u32(z,V,2);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss(reinterpret_cast<float *>(z),vResult);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Store the W component into a 32 bit integer location in memory.
_Use_decl_annotations_
inline void XMVectorGetIntWPtr(uint32_t *w, FXMVECTOR V)
{
    assert( w != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    *w = V.vector4_u32[3];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    vst1q_lane_u32(w,V,3);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,3,3,3));
    _mm_store_ss(reinterpret_cast<float *>(w),vResult);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Set a single indexed floating point component
inline XMVECTOR XMVectorSetByIndex(FXMVECTOR V, float f, size_t i)
{
    assert( i < 4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U = V;
    U.vector4_f32[i] = f;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR U = V;
    U.n128_f32[i] = f;
    return U;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR U = V;
    U.m128_f32[i] = f;
    return U;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets the X component of a vector to a passed floating point value
inline XMVECTOR XMVectorSetX(FXMVECTOR V, float x)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = x;
    U.vector4_f32[1] = V.vector4_f32[1];
    U.vector4_f32[2] = V.vector4_f32[2];
    U.vector4_f32[3] = V.vector4_f32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_f32(x,V,0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_set_ss(x);
    vResult = _mm_move_ss(V,vResult);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Y component of a vector to a passed floating point value
inline XMVECTOR XMVectorSetY(FXMVECTOR V, float y)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = V.vector4_f32[0];
    U.vector4_f32[1] = y;
    U.vector4_f32[2] = V.vector4_f32[2];
    U.vector4_f32[3] = V.vector4_f32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_f32(y,V,1);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap y and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,2,0,1));
    // Convert input to vector
    XMVECTOR vTemp = _mm_set_ss(y);
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap y and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,2,0,1));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}
// Sets the Z component of a vector to a passed floating point value
inline XMVECTOR XMVectorSetZ(FXMVECTOR V, float z)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = V.vector4_f32[0];
    U.vector4_f32[1] = V.vector4_f32[1];
    U.vector4_f32[2] = z;
    U.vector4_f32[3] = V.vector4_f32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_f32(z,V,2);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap z and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,0,1,2));
    // Convert input to vector
    XMVECTOR vTemp = _mm_set_ss(z);
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap z and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,0,1,2));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the W component of a vector to a passed floating point value
inline XMVECTOR XMVectorSetW(FXMVECTOR V, float w)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = V.vector4_f32[0];
    U.vector4_f32[1] = V.vector4_f32[1];
    U.vector4_f32[2] = V.vector4_f32[2];
    U.vector4_f32[3] = w;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_f32(w,V,3);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap w and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,2,1,3));
    // Convert input to vector
    XMVECTOR vTemp = _mm_set_ss(w);
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap w and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,2,1,3));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets a component of a vector to a floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetByIndexPtr(FXMVECTOR V, const float *f, size_t i)
{
    assert( f != nullptr );
    assert( i < 4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U = V;
    U.vector4_f32[i] = *f;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR U = V;
    U.n128_f32[i] = *f;
    return U;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR U = V;
    U.m128_f32[i] = *f;
    return U;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets the X component of a vector to a floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetXPtr(FXMVECTOR V, const float *x)
{
    assert( x != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = *x;
    U.vector4_f32[1] = V.vector4_f32[1];
    U.vector4_f32[2] = V.vector4_f32[2];
    U.vector4_f32[3] = V.vector4_f32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_f32(x,V,0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_load_ss(x);
    vResult = _mm_move_ss(V,vResult);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Y component of a vector to a floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetYPtr(FXMVECTOR V, const float *y)
{
    assert( y != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = V.vector4_f32[0];
    U.vector4_f32[1] = *y;
    U.vector4_f32[2] = V.vector4_f32[2];
    U.vector4_f32[3] = V.vector4_f32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_f32(y,V,1);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap y and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,2,0,1));
    // Convert input to vector
    XMVECTOR vTemp = _mm_load_ss(y);
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap y and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,2,0,1));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Z component of a vector to a floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetZPtr(FXMVECTOR V, const float *z)
{
    assert( z != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = V.vector4_f32[0];
    U.vector4_f32[1] = V.vector4_f32[1];
    U.vector4_f32[2] = *z;
    U.vector4_f32[3] = V.vector4_f32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_f32(z,V,2);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap z and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,0,1,2));
    // Convert input to vector
    XMVECTOR vTemp = _mm_load_ss(z);
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap z and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,0,1,2));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the W component of a vector to a floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetWPtr(FXMVECTOR V, const float *w)
{
    assert( w != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_f32[0] = V.vector4_f32[0];
    U.vector4_f32[1] = V.vector4_f32[1];
    U.vector4_f32[2] = V.vector4_f32[2];
    U.vector4_f32[3] = *w;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_f32(w,V,3);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap w and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,2,1,3));
    // Convert input to vector
    XMVECTOR vTemp = _mm_load_ss(w);
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap w and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,2,1,3));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets a component of a vector to an integer passed by value
inline XMVECTOR XMVectorSetIntByIndex(FXMVECTOR V, uint32_t x, size_t i)
{
    assert( i < 4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U = V;
    U.vector4_u32[i] = x;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORU32 tmp;
    tmp.v = V;
    tmp.u[i] = x;
    return tmp;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTORU32 tmp;
    tmp.v = V;
    tmp.u[i] = x;
    return tmp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets the X component of a vector to an integer passed by value
inline XMVECTOR XMVectorSetIntX(FXMVECTOR V, uint32_t x)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = x;
    U.vector4_u32[1] = V.vector4_u32[1];
    U.vector4_u32[2] = V.vector4_u32[2];
    U.vector4_u32[3] = V.vector4_u32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_u32(x,V,0);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cvtsi32_si128(x);
    XMVECTOR vResult = _mm_move_ss(V,_mm_castsi128_ps(vTemp));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Y component of a vector to an integer passed by value
inline XMVECTOR XMVectorSetIntY(FXMVECTOR V, uint32_t y)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = V.vector4_u32[0];
    U.vector4_u32[1] = y;
    U.vector4_u32[2] = V.vector4_u32[2];
    U.vector4_u32[3] = V.vector4_u32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_u32(y,V,1);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap y and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,2,0,1));
    // Convert input to vector
    __m128i vTemp = _mm_cvtsi32_si128(y);
    // Replace the x component
    vResult = _mm_move_ss(vResult,_mm_castsi128_ps(vTemp));
    // Swap y and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,2,0,1));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Z component of a vector to an integer passed by value
inline XMVECTOR XMVectorSetIntZ(FXMVECTOR V, uint32_t z)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = V.vector4_u32[0];
    U.vector4_u32[1] = V.vector4_u32[1];
    U.vector4_u32[2] = z;
    U.vector4_u32[3] = V.vector4_u32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_u32(z,V,2);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap z and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,0,1,2));
    // Convert input to vector
    __m128i vTemp = _mm_cvtsi32_si128(z);
    // Replace the x component
    vResult = _mm_move_ss(vResult,_mm_castsi128_ps(vTemp));
    // Swap z and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,0,1,2));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the W component of a vector to an integer passed by value
inline XMVECTOR XMVectorSetIntW(FXMVECTOR V, uint32_t w)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = V.vector4_u32[0];
    U.vector4_u32[1] = V.vector4_u32[1];
    U.vector4_u32[2] = V.vector4_u32[2];
    U.vector4_u32[3] = w;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsetq_lane_u32(w,V,3);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap w and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,2,1,3));
    // Convert input to vector
    __m128i vTemp = _mm_cvtsi32_si128(w);
    // Replace the x component
    vResult = _mm_move_ss(vResult,_mm_castsi128_ps(vTemp));
    // Swap w and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,2,1,3));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets a component of a vector to an integer value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetIntByIndexPtr(FXMVECTOR V, const uint32_t *x, size_t i)
{
    assert( x != nullptr );
    assert( i < 4 );
    _Analysis_assume_( i < 4 );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U = V;
    U.vector4_u32[i] = *x;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORU32 tmp;
    tmp.v = V;
    tmp.u[i] = *x;
    return tmp;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTORU32 tmp;
    tmp.v = V;
    tmp.u[i] = *x;
    return tmp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

// Sets the X component of a vector to an integer value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetIntXPtr(FXMVECTOR V, const uint32_t *x)
{
    assert( x != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = *x;
    U.vector4_u32[1] = V.vector4_u32[1];
    U.vector4_u32[2] = V.vector4_u32[2];
    U.vector4_u32[3] = V.vector4_u32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_u32(x,V,0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_load_ss(reinterpret_cast<const float *>(x));
    XMVECTOR vResult = _mm_move_ss(V,vTemp);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Y component of a vector to an integer value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetIntYPtr(FXMVECTOR V, const uint32_t *y)
{
    assert( y != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = V.vector4_u32[0];
    U.vector4_u32[1] = *y;
    U.vector4_u32[2] = V.vector4_u32[2];
    U.vector4_u32[3] = V.vector4_u32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_u32(y,V,1);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap y and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,2,0,1));
    // Convert input to vector
    XMVECTOR vTemp = _mm_load_ss(reinterpret_cast<const float *>(y));
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap y and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,2,0,1));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the Z component of a vector to an integer value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetIntZPtr(FXMVECTOR V, const uint32_t *z)
{
    assert( z != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = V.vector4_u32[0];
    U.vector4_u32[1] = V.vector4_u32[1];
    U.vector4_u32[2] = *z;
    U.vector4_u32[3] = V.vector4_u32[3];
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_u32(z,V,2);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap z and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,0,1,2));
    // Convert input to vector
    XMVECTOR vTemp = _mm_load_ss(reinterpret_cast<const float *>(z));
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap z and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,0,1,2));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

// Sets the W component of a vector to an integer value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XMVectorSetIntWPtr(FXMVECTOR V, const uint32_t *w)
{
    assert( w != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR U;
    U.vector4_u32[0] = V.vector4_u32[0];
    U.vector4_u32[1] = V.vector4_u32[1];
    U.vector4_u32[2] = V.vector4_u32[2];
    U.vector4_u32[3] = *w;
    return U;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vld1q_lane_u32(w,V,3);
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap w and x
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,2,1,3));
    // Convert input to vector
    XMVECTOR vTemp = _mm_load_ss(reinterpret_cast<const float *>(w));
    // Replace the x component
    vResult = _mm_move_ss(vResult,vTemp);
    // Swap w and x again
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,2,1,3));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSwizzle
(
    FXMVECTOR V,
    uint32_t E0,
    uint32_t E1,
    uint32_t E2,
    uint32_t E3
)
{
    assert( (E0 < 4) && (E1 < 4) && (E2 < 4) && (E3 < 4) );
    _Analysis_assume_( (E0 < 4) && (E1 < 4) && (E2 < 4) && (E3 < 4) );
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result = { V.vector4_f32[E0],
                        V.vector4_f32[E1],
                        V.vector4_f32[E2],
                        V.vector4_f32[E3] };
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t ControlElement[ 4 ] =
    {
#ifdef _XM_LITTLEENDIAN_
        0x03020100, // XM_SWIZZLE_X
        0x07060504, // XM_SWIZZLE_Y
        0x0B0A0908, // XM_SWIZZLE_Z
        0x0F0E0D0C, // XM_SWIZZLE_W
#else
        0x00010203, // XM_SWIZZLE_X
        0x04050607, // XM_SWIZZLE_Y
        0x08090A0B, // XM_SWIZZLE_Z
        0x0C0D0E0F, // XM_SWIZZLE_W
#endif
    };

    int8x8x2_t tbl;
    tbl.val[0] = vget_low_f32(V);
    tbl.val[1] = vget_high_f32(V);

    __n64 idx = vcreate_u32( ((uint64_t)ControlElement[E0]) | (((uint64_t)ControlElement[E1]) << 32) );
    const __n64 rL = vtbl2_u8( tbl, idx );

    idx = vcreate_u32( ((uint64_t)ControlElement[E2]) | (((uint64_t)ControlElement[E3]) << 32) );
    const __n64 rH = vtbl2_u8( tbl, idx );

    return vcombine_f32( rL, rH );
#elif defined(_XM_VMX128_INTRINSICS_)
#else
    const uint32_t *aPtr = (const uint32_t* )(&V);

    XMVECTOR Result;
    uint32_t *pWork = (uint32_t*)(&Result);

    pWork[0] = aPtr[E0];
    pWork[1] = aPtr[E1];
    pWork[2] = aPtr[E2];
    pWork[3] = aPtr[E3];

    return Result;
#endif
}

//------------------------------------------------------------------------------
inline XMVECTOR XMVectorPermute
(
    FXMVECTOR V1,
    FXMVECTOR V2,
    uint32_t PermuteX,
    uint32_t PermuteY,
    uint32_t PermuteZ,
    uint32_t PermuteW
)
{
    assert( PermuteX <= 7 && PermuteY <= 7 && PermuteZ <= 7 && PermuteW <= 7 );
    _Analysis_assume_( PermuteX <= 7 && PermuteY <= 7 && PermuteZ <= 7 && PermuteW <= 7 );

#if defined(_XM_ARM_NEON_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const uint32_t ControlElement[ 8 ] =
    {
#ifdef _XM_LITTLEENDIAN_
        0x03020100, // XM_PERMUTE_0X
        0x07060504, // XM_PERMUTE_0Y
        0x0B0A0908, // XM_PERMUTE_0Z
        0x0F0E0D0C, // XM_PERMUTE_0W
        0x13121110, // XM_PERMUTE_1X
        0x17161514, // XM_PERMUTE_1Y
        0x1B1A1918, // XM_PERMUTE_1Z
        0x1F1E1D1C, // XM_PERMUTE_1W
#else
        0x00010203, // XM_PERMUTE_0X
        0x04050607, // XM_PERMUTE_0Y
        0x08090A0B, // XM_PERMUTE_0Z
        0x0C0D0E0F, // XM_PERMUTE_0W
        0x10111213, // XM_PERMUTE_1X
        0x14151617, // XM_PERMUTE_1Y
        0x18191A1B, // XM_PERMUTE_1Z
        0x1C1D1E1F, // XM_PERMUTE_1W
#endif
    };

    int8x8x4_t tbl;
    tbl.val[0] = vget_low_f32(V1);
    tbl.val[1] = vget_high_f32(V1);
    tbl.val[2] = vget_low_f32(V2);
    tbl.val[3] = vget_high_f32(V2);

    __n64 idx = vcreate_u32( ((uint64_t)ControlElement[PermuteX]) | (((uint64_t)ControlElement[PermuteY]) << 32) );
    const __n64 rL = vtbl4_u8( tbl, idx );

    idx = vcreate_u32( ((uint64_t)ControlElement[PermuteZ]) | (((uint64_t)ControlElement[PermuteW]) << 32) );
    const __n64 rH = vtbl4_u8( tbl, idx );

    return vcombine_f32( rL, rH );
#elif defined(_XM_VMX128_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
#else
 
    const uint32_t *aPtr[2];
    aPtr[0] = (const uint32_t* )(&V1);
    aPtr[1] = (const uint32_t* )(&V2);

    XMVECTOR Result;
    uint32_t *pWork = (uint32_t*)(&Result);

    const uint32_t i0 = PermuteX & 3;
    const uint32_t vi0 = PermuteX >> 2;
    pWork[0] = aPtr[vi0][i0];

    const uint32_t i1 = PermuteY & 3;
    const uint32_t vi1 = PermuteY >> 2;
    pWork[1] = aPtr[vi1][i1];

    const uint32_t i2 = PermuteZ & 3;
    const uint32_t vi2 = PermuteZ >> 2;
    pWork[2] = aPtr[vi2][i2];

    const uint32_t i3 = PermuteW & 3;
    const uint32_t vi3 = PermuteW >> 2;
    pWork[3] = aPtr[vi3][i3];

    return Result;
#endif
}

//------------------------------------------------------------------------------
// Define a control vector to be used in XMVectorSelect 
// operations.  The four integers specified in XMVectorSelectControl
// serve as indices to select between components in two vectors.
// The first index controls selection for the first component of 
// the vectors involved in a select operation, the second index 
// controls selection for the second component etc.  A value of
// zero for an index causes the corresponding component from the first 
// vector to be selected whereas a one causes the component from the
// second vector to be selected instead.

inline XMVECTOR XMVectorSelectControl
(
    uint32_t VectorIndex0, 
    uint32_t VectorIndex1, 
    uint32_t VectorIndex2, 
    uint32_t VectorIndex3
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    // x=Index0,y=Index1,z=Index2,w=Index3
    __m128i vTemp = _mm_set_epi32(VectorIndex3,VectorIndex2,VectorIndex1,VectorIndex0);
    // Any non-zero entries become 0xFFFFFFFF else 0
    vTemp = _mm_cmpgt_epi32(vTemp,g_XMZero);
    return reinterpret_cast<__m128 *>(&vTemp)[0];
#elif defined(_XM_ARM_NEON_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    __n64 V0 = vcreate_s32(((uint64_t)VectorIndex0) | ((uint64_t)VectorIndex1 << 32));
    __n64 V1 = vcreate_s32(((uint64_t)VectorIndex2) | ((uint64_t)VectorIndex3 << 32));
    __n128 vTemp = vcombine_s32(V0, V1);
    // Any non-zero entries become 0xFFFFFFFF else 0
    return vcgtq_s32(vTemp,g_XMZero);
#else
    XMVECTOR    ControlVector;
    const uint32_t  ControlElement[] =
                {
                    XM_SELECT_0,
                    XM_SELECT_1
                };

    assert(VectorIndex0 < 2);
    assert(VectorIndex1 < 2);
    assert(VectorIndex2 < 2);
    assert(VectorIndex3 < 2);
    _Analysis_assume_(VectorIndex0 < 2);
    _Analysis_assume_(VectorIndex1 < 2);
    _Analysis_assume_(VectorIndex2 < 2);
    _Analysis_assume_(VectorIndex3 < 2);

    ControlVector.vector4_u32[0] = ControlElement[VectorIndex0];
    ControlVector.vector4_u32[1] = ControlElement[VectorIndex1];
    ControlVector.vector4_u32[2] = ControlElement[VectorIndex2];
    ControlVector.vector4_u32[3] = ControlElement[VectorIndex3];

    return ControlVector;

#endif
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSelect
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR Control
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = (V1.vector4_u32[0] & ~Control.vector4_u32[0]) | (V2.vector4_u32[0] & Control.vector4_u32[0]);
    Result.vector4_u32[1] = (V1.vector4_u32[1] & ~Control.vector4_u32[1]) | (V2.vector4_u32[1] & Control.vector4_u32[1]);
    Result.vector4_u32[2] = (V1.vector4_u32[2] & ~Control.vector4_u32[2]) | (V2.vector4_u32[2] & Control.vector4_u32[2]);
    Result.vector4_u32[3] = (V1.vector4_u32[3] & ~Control.vector4_u32[3]) | (V2.vector4_u32[3] & Control.vector4_u32[3]);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vbslq_f32( Control, V2, V1 );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp1 = _mm_andnot_ps(Control,V1);
    XMVECTOR vTemp2 = _mm_and_ps(V2,Control);
    return _mm_or_ps(vTemp1,vTemp2);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMergeXY
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = V1.vector4_u32[0];
    Result.vector4_u32[1] = V2.vector4_u32[0];
    Result.vector4_u32[2] = V1.vector4_u32[1];
    Result.vector4_u32[3] = V2.vector4_u32[1];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vzipq_f32( V1, V2 ).val[0];
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_unpacklo_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMergeZW
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = V1.vector4_u32[2];
    Result.vector4_u32[1] = V2.vector4_u32[2];
    Result.vector4_u32[2] = V1.vector4_u32[3];
    Result.vector4_u32[3] = V2.vector4_u32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vzipq_f32( V1, V2 ).val[1];
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_unpackhi_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorShiftLeft(FXMVECTOR V1, FXMVECTOR V2, uint32_t Elements)
{
    assert( Elements < 4 );
    _Analysis_assume_( Elements < 4 );
    return XMVectorPermute(V1, V2, Elements, ((Elements) + 1), ((Elements) + 2), ((Elements) + 3));
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorRotateLeft(FXMVECTOR V, uint32_t Elements)
{
    assert( Elements < 4 );
    _Analysis_assume_( Elements < 4 );
    return XMVectorSwizzle( V, Elements & 3, (Elements + 1) & 3, (Elements + 2) & 3, (Elements + 3) & 3 );
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorRotateRight(FXMVECTOR V, uint32_t Elements)
{
    assert( Elements < 4 );
    _Analysis_assume_( Elements < 4 );
    return XMVectorSwizzle( V, (4 - (Elements)) & 3, (5 - (Elements)) & 3, (6 - (Elements)) & 3, (7 - (Elements)) & 3 );
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorInsert(FXMVECTOR VD, FXMVECTOR VS, uint32_t VSLeftRotateElements,
                                  uint32_t Select0, uint32_t Select1, uint32_t Select2, uint32_t Select3)
{
    XMVECTOR Control = XMVectorSelectControl(Select0&1, Select1&1, Select2&1, Select3&1);
    return XMVectorSelect( VD, XMVectorRotateLeft(VS, VSLeftRotateElements), Control );
}

//------------------------------------------------------------------------------
// Comparison operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_f32[0] == V2.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_f32[1] == V2.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_f32[2] == V2.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_f32[3] == V2.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vceqq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cmpeq_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMVECTOR XMVectorEqualR
(
    uint32_t*    pCR,
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    assert( pCR != nullptr );
#if defined(_XM_NO_INTRINSICS_)
    uint32_t ux = (V1.vector4_f32[0] == V2.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    uint32_t uy = (V1.vector4_f32[1] == V2.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    uint32_t uz = (V1.vector4_f32[2] == V2.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    uint32_t uw = (V1.vector4_f32[3] == V2.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
    uint32_t CR = 0;
    if (ux&uy&uz&uw)
    {
        // All elements are greater
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!(ux|uy|uz|uw))
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;

    XMVECTOR Control;
    Control.vector4_u32[0] = ux;
    Control.vector4_u32[1] = uy;
    Control.vector4_u32[2] = uz;
    Control.vector4_u32[3] = uw;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        // All elements are equal
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        // All elements are not equal
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
    uint32_t CR = 0;
    int iTest = _mm_movemask_ps(vTemp);
    if (iTest==0xf)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vTemp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Treat the components of the vectors as unsigned integers and
// compare individual bits between the two.  This is useful for
// comparing control vectors and result vectors returned from
// other comparison operations.

inline XMVECTOR XMVectorEqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_u32[0] == V2.vector4_u32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_u32[1] == V2.vector4_u32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_u32[2] == V2.vector4_u32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_u32[3] == V2.vector4_u32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vceqq_u32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_cmpeq_epi32( _mm_castps_si128(V1),_mm_castps_si128(V2) );
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMVECTOR XMVectorEqualIntR
(
    uint32_t*    pCR,
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    assert( pCR != nullptr );
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control = XMVectorEqualInt(V1, V2);

    *pCR = 0;
    if (XMVector4EqualInt(Control, XMVectorTrueInt()))
    {
        // All elements are equal
        *pCR |= XM_CRMASK_CR6TRUE;
    }
    else if (XMVector4EqualInt(Control, XMVectorFalseInt()))
    {
        // All elements are not equal
        *pCR |= XM_CRMASK_CR6FALSE;
    }
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        // All elements are equal
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        // All elements are not equal
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_cmpeq_epi32( _mm_castps_si128(V1),_mm_castps_si128(V2) );
    int iTemp = _mm_movemask_ps(reinterpret_cast<const __m128*>(&V)[0]);
    uint32_t CR = 0;
    if (iTemp==0x0F)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTemp)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorNearEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR Epsilon
)
{
#if defined(_XM_NO_INTRINSICS_)

    float fDeltax = V1.vector4_f32[0]-V2.vector4_f32[0];
    float fDeltay = V1.vector4_f32[1]-V2.vector4_f32[1];
    float fDeltaz = V1.vector4_f32[2]-V2.vector4_f32[2];
    float fDeltaw = V1.vector4_f32[3]-V2.vector4_f32[3];

    fDeltax = fabsf(fDeltax);
    fDeltay = fabsf(fDeltay);
    fDeltaz = fabsf(fDeltaz);
    fDeltaw = fabsf(fDeltaw);

    XMVECTOR Control;
    Control.vector4_u32[0] = (fDeltax <= Epsilon.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[1] = (fDeltay <= Epsilon.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[2] = (fDeltaz <= Epsilon.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[3] = (fDeltaw <= Epsilon.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR vDelta = vsubq_f32(V1,V2);
    return vacleq_f32( vDelta, Epsilon );
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the difference
    XMVECTOR vDelta = _mm_sub_ps(V1,V2);
    // Get the absolute value of the difference
    XMVECTOR vTemp = _mm_setzero_ps();
    vTemp = _mm_sub_ps(vTemp,vDelta);
    vTemp = _mm_max_ps(vTemp,vDelta);
    vTemp = _mm_cmple_ps(vTemp,Epsilon);
    return vTemp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorNotEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_f32[0] != V2.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_f32[1] != V2.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_f32[2] != V2.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_f32[3] != V2.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmvnq_u32(vceqq_f32(V1, V2));
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cmpneq_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorNotEqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_u32[0] != V2.vector4_u32[0]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[1] = (V1.vector4_u32[1] != V2.vector4_u32[1]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[2] = (V1.vector4_u32[2] != V2.vector4_u32[2]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[3] = (V1.vector4_u32[3] != V2.vector4_u32[3]) ? 0xFFFFFFFFU : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmvnq_u32(vceqq_u32(V1, V2));
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_cmpeq_epi32( _mm_castps_si128(V1),_mm_castps_si128(V2) );
    return _mm_xor_ps(reinterpret_cast<__m128 *>(&V)[0],g_XMNegOneMask);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorGreater
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_f32[0] > V2.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_f32[1] > V2.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_f32[2] > V2.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_f32[3] > V2.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vcgtq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cmpgt_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMVECTOR XMVectorGreaterR
(
    uint32_t*    pCR,
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    assert( pCR != nullptr );
#if defined(_XM_NO_INTRINSICS_)

    uint32_t ux = (V1.vector4_f32[0] > V2.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    uint32_t uy = (V1.vector4_f32[1] > V2.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    uint32_t uz = (V1.vector4_f32[2] > V2.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    uint32_t uw = (V1.vector4_f32[3] > V2.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
    uint32_t CR = 0;
    if (ux&uy&uz&uw)
    {
        // All elements are greater
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!(ux|uy|uz|uw))
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;

    XMVECTOR Control;
    Control.vector4_u32[0] = ux;
    Control.vector4_u32[1] = uy;
    Control.vector4_u32[2] = uz;
    Control.vector4_u32[3] = uw;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgtq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        // All elements are greater
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
    uint32_t CR = 0;
    int iTest = _mm_movemask_ps(vTemp);
    if (iTest==0xf)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vTemp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorGreaterOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_f32[0] >= V2.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_f32[1] >= V2.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_f32[2] >= V2.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_f32[3] >= V2.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vcgeq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cmpge_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMVECTOR XMVectorGreaterOrEqualR
(
    uint32_t*    pCR,
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    assert( pCR != nullptr );
#if defined(_XM_NO_INTRINSICS_)

    uint32_t ux = (V1.vector4_f32[0] >= V2.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    uint32_t uy = (V1.vector4_f32[1] >= V2.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    uint32_t uz = (V1.vector4_f32[2] >= V2.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    uint32_t uw = (V1.vector4_f32[3] >= V2.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
    uint32_t CR = 0;
    if (ux&uy&uz&uw)
    {
        // All elements are greater
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!(ux|uy|uz|uw))
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;

    XMVECTOR Control;
    Control.vector4_u32[0] = ux;
    Control.vector4_u32[1] = uy;
    Control.vector4_u32[2] = uz;
    Control.vector4_u32[3] = uw;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgeq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        // All elements are greater or equal
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        // All elements are not greater or equal
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    uint32_t CR = 0;
    int iTest = _mm_movemask_ps(vTemp);
    if (iTest==0xf)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        // All elements are not greater
        CR = XM_CRMASK_CR6FALSE;
    }
    *pCR = CR;
    return vTemp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorLess
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_f32[0] < V2.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_f32[1] < V2.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_f32[2] < V2.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_f32[3] < V2.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vcltq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cmplt_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorLessOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V1.vector4_f32[0] <= V2.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V1.vector4_f32[1] <= V2.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V1.vector4_f32[2] <= V2.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V1.vector4_f32[3] <= V2.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vcleq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cmple_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorInBounds
(
    FXMVECTOR V, 
    FXMVECTOR Bounds
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = (V.vector4_f32[0] <= Bounds.vector4_f32[0] && V.vector4_f32[0] >= -Bounds.vector4_f32[0]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[1] = (V.vector4_f32[1] <= Bounds.vector4_f32[1] && V.vector4_f32[1] >= -Bounds.vector4_f32[1]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[2] = (V.vector4_f32[2] <= Bounds.vector4_f32[2] && V.vector4_f32[2] >= -Bounds.vector4_f32[2]) ? 0xFFFFFFFF : 0;
    Control.vector4_u32[3] = (V.vector4_f32[3] <= Bounds.vector4_f32[3] && V.vector4_f32[3] >= -Bounds.vector4_f32[3]) ? 0xFFFFFFFF : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = vcleq_f32(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = vnegq_f32(Bounds);
    // Test if greater or equal (Reversed)
    vTemp2 = vcleq_f32(vTemp2,V);
    // Blend answers
    vTemp1 = vandq_u32(vTemp1,vTemp2);
    return vTemp1;
#elif defined(_XM_SSE_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = _mm_cmple_ps(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = _mm_mul_ps(Bounds,g_XMNegativeOne);
    // Test if greater or equal (Reversed)
    vTemp2 = _mm_cmple_ps(vTemp2,V);
    // Blend answers
    vTemp1 = _mm_and_ps(vTemp1,vTemp2);
    return vTemp1;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMVECTOR XMVectorInBoundsR
(
    uint32_t*    pCR,
    FXMVECTOR V, 
    FXMVECTOR Bounds
)
{
    assert( pCR != nullptr );
#if defined(_XM_NO_INTRINSICS_)

    uint32_t ux = (V.vector4_f32[0] <= Bounds.vector4_f32[0] && V.vector4_f32[0] >= -Bounds.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    uint32_t uy = (V.vector4_f32[1] <= Bounds.vector4_f32[1] && V.vector4_f32[1] >= -Bounds.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    uint32_t uz = (V.vector4_f32[2] <= Bounds.vector4_f32[2] && V.vector4_f32[2] >= -Bounds.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    uint32_t uw = (V.vector4_f32[3] <= Bounds.vector4_f32[3] && V.vector4_f32[3] >= -Bounds.vector4_f32[3]) ? 0xFFFFFFFFU : 0;

    uint32_t CR = 0;
    if (ux&uy&uz&uw)
    {
        // All elements are in bounds
        CR = XM_CRMASK_CR6BOUNDS;
    }
    *pCR = CR;

    XMVECTOR Control;
    Control.vector4_u32[0] = ux;
    Control.vector4_u32[1] = uy;
    Control.vector4_u32[2] = uz;
    Control.vector4_u32[3] = uw;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = vcleq_f32(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = vnegq_f32(Bounds);
    // Test if greater or equal (Reversed)
    vTemp2 = vcleq_f32(vTemp2,V);
    // Blend answers
    vTemp1 = vandq_u32(vTemp1,vTemp2);
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTemp1), vget_high_u8(vTemp1));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        // All elements are in bounds
        CR = XM_CRMASK_CR6BOUNDS;
    }
    *pCR = CR;
    return vTemp1;
#elif defined(_XM_SSE_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = _mm_cmple_ps(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = _mm_mul_ps(Bounds,g_XMNegativeOne);
    // Test if greater or equal (Reversed)
    vTemp2 = _mm_cmple_ps(vTemp2,V);
    // Blend answers
    vTemp1 = _mm_and_ps(vTemp1,vTemp2);

    uint32_t CR = 0;
    if (_mm_movemask_ps(vTemp1)==0xf) {
        // All elements are in bounds
        CR = XM_CRMASK_CR6BOUNDS;
    }
    *pCR = CR;
    return vTemp1;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorIsNaN
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = XMISNAN(V.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[1] = XMISNAN(V.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[2] = XMISNAN(V.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[3] = XMISNAN(V.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test against itself. NaN is always not equal
    __n128 vTempNan = vceqq_f32( V, V );
    // Flip results
    return vmvnq_u32( vTempNan );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test against itself. NaN is always not equal
    return _mm_cmpneq_ps(V,V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorIsInfinite
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Control;
    Control.vector4_u32[0] = XMISINF(V.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[1] = XMISINF(V.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[2] = XMISINF(V.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
    Control.vector4_u32[3] = XMISINF(V.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
    return Control;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Mask off the sign bit
    __n128 vTemp = vandq_u32(V,g_XMAbsMask);
    // Compare to infinity
    vTemp = vceqq_f32(vTemp,g_XMInfinity);
    // If any are infinity, the signs are true.
    return vTemp;
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bit
    __m128 vTemp = _mm_and_ps(V,g_XMAbsMask);
    // Compare to infinity
    vTemp = _mm_cmpeq_ps(vTemp,g_XMInfinity);
    // If any are infinity, the signs are true.
    return vTemp;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Rounding and clamping operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMin
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = (V1.vector4_f32[0] < V2.vector4_f32[0]) ? V1.vector4_f32[0] : V2.vector4_f32[0];
    Result.vector4_f32[1] = (V1.vector4_f32[1] < V2.vector4_f32[1]) ? V1.vector4_f32[1] : V2.vector4_f32[1];
    Result.vector4_f32[2] = (V1.vector4_f32[2] < V2.vector4_f32[2]) ? V1.vector4_f32[2] : V2.vector4_f32[2];
    Result.vector4_f32[3] = (V1.vector4_f32[3] < V2.vector4_f32[3]) ? V1.vector4_f32[3] : V2.vector4_f32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vminq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_min_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMax
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = (V1.vector4_f32[0] > V2.vector4_f32[0]) ? V1.vector4_f32[0] : V2.vector4_f32[0];
    Result.vector4_f32[1] = (V1.vector4_f32[1] > V2.vector4_f32[1]) ? V1.vector4_f32[1] : V2.vector4_f32[1];
    Result.vector4_f32[2] = (V1.vector4_f32[2] > V2.vector4_f32[2]) ? V1.vector4_f32[2] : V2.vector4_f32[2];
    Result.vector4_f32[3] = (V1.vector4_f32[3] > V2.vector4_f32[3]) ? V1.vector4_f32[3] : V2.vector4_f32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmaxq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_max_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorRound
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    const XMVECTOR Zero = XMVectorZero();
    const XMVECTOR BiasPos = XMVectorReplicate(0.5f);
    const XMVECTOR BiasNeg = XMVectorReplicate(-0.5f);

    XMVECTOR Bias = XMVectorLess(V, Zero);
    Bias = XMVectorSelect(BiasPos, BiasNeg, Bias);
    XMVECTOR Result = XMVectorAdd(V, Bias);
    Result = XMVectorTruncate(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vTest = vabsq_f32( V );
    vTest = vcltq_f32( vTest, g_XMNoFraction );

    __n128 Bias = vcltq_f32( V, vdupq_n_u32(0) );

    __n128 BiasPos = vdupq_n_f32( 0.5f );
    __n128 BiasNeg = vdupq_n_f32( -0.5f );
    Bias = vbslq_f32( Bias, BiasNeg, BiasPos );
    __n128 V0 = vaddq_f32( V, Bias );
    __n128 vInt = vcvtq_s32_f32( V0 );
    __n128 vResult = vcvtq_f32_s32( vInt );

    // All numbers less than 8388608 will use the round to int
    // All others, use the ORIGINAL value
    return vbslq_f32( vTest, vResult, V );
#elif defined(_XM_SSE_INTRINSICS_)
    // To handle NAN, INF and numbers greater than 8388608, use masking
    // Get the abs value
    __m128i vTest = _mm_and_si128(_mm_castps_si128(V),g_XMAbsMask);
    // Test for greater than 8388608 (All floats with NO fractionals, NAN and INF
    vTest = _mm_cmplt_epi32(vTest,g_XMNoFraction);
    // Convert to int and back to float for rounding
    __m128i vInt = _mm_cvtps_epi32(V);
    // Convert back to floats
    XMVECTOR vResult = _mm_cvtepi32_ps(vInt);
    // All numbers less than 8388608 will use the round to int
    vResult = _mm_and_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    // All others, use the ORIGINAL value
    vTest = _mm_andnot_si128(vTest,_mm_castps_si128(V));
    vResult = _mm_or_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorTruncate
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    uint32_t     i;

    // Avoid C4701
    Result.vector4_f32[0] = 0.0f;

    for (i = 0; i < 4; i++)
    {
        if (XMISNAN(V.vector4_f32[i]))
        {
            Result.vector4_u32[i] = 0x7FC00000;
        }
        else if (fabsf(V.vector4_f32[i]) < 8388608.0f)
        {
            Result.vector4_f32[i] = (float)((int32_t)V.vector4_f32[i]);
        }
        else
        {
            Result.vector4_f32[i] = V.vector4_f32[i];
        }
    }
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vTest = vabsq_f32( V );
    vTest = vcltq_f32( vTest, g_XMNoFraction );

    __n128 vInt = vcvtq_s32_f32( V );
    __n128 vResult = vcvtq_f32_s32( vInt );

    // All numbers less than 8388608 will use the round to int
    // All others, use the ORIGINAL value
    return vbslq_f32( vTest, vResult, V );
#elif defined(_XM_SSE_INTRINSICS_)
    // To handle NAN, INF and numbers greater than 8388608, use masking
    // Get the abs value
    __m128i vTest = _mm_and_si128(_mm_castps_si128(V),g_XMAbsMask);
    // Test for greater than 8388608 (All floats with NO fractionals, NAN and INF
    vTest = _mm_cmplt_epi32(vTest,g_XMNoFraction);
    // Convert to int and back to float for rounding with truncation
    __m128i vInt = _mm_cvttps_epi32(V);
    // Convert back to floats
    XMVECTOR vResult = _mm_cvtepi32_ps(vInt);
    // All numbers less than 8388608 will use the round to int
    vResult = _mm_and_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    // All others, use the ORIGINAL value
    vTest = _mm_andnot_si128(vTest,_mm_castps_si128(V));
    vResult = _mm_or_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorFloor
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR vResult = {
        floorf(V.vector4_f32[0]),
        floorf(V.vector4_f32[1]),
        floorf(V.vector4_f32[2]),
        floorf(V.vector4_f32[3])
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 V0 = vsubq_f32( V, vdupq_n_u32(0x3EFFFFA0) );
    return XMVectorRound(V0);
#elif defined(_XM_SSE_INTRINSICS_)
    // To handle NAN, INF and numbers greater than 8388608, use masking
    // Get the abs value
    __m128i vTest = _mm_and_si128(_mm_castps_si128(V),g_XMAbsMask);
    // Test for greater than 8388608 (All floats with NO fractionals, NAN and INF
    vTest = _mm_cmplt_epi32(vTest,g_XMNoFraction);
    // Convert to int and back to float for rounding
    XMVECTOR vResult = _mm_sub_ps(V,g_XMOneHalfMinusEpsilon);
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Convert back to floats
    vResult = _mm_cvtepi32_ps(vInt);
    // All numbers less than 8388608 will use the round to int
    vResult = _mm_and_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    // All others, use the ORIGINAL value
    vTest = _mm_andnot_si128(vTest,_mm_castps_si128(V));
    vResult = _mm_or_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorCeiling
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {
        ceilf(V.vector4_f32[0]),
        ceilf(V.vector4_f32[1]),
        ceilf(V.vector4_f32[2]),
        ceilf(V.vector4_f32[3])
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 V0 = vaddq_f32( V, vdupq_n_u32(0x3EFFFFA0) );
    return XMVectorRound(V0);
#elif defined(_XM_SSE_INTRINSICS_)
    // To handle NAN, INF and numbers greater than 8388608, use masking
    // Get the abs value
    __m128i vTest = _mm_and_si128(_mm_castps_si128(V),g_XMAbsMask);
    // Test for greater than 8388608 (All floats with NO fractionals, NAN and INF
    vTest = _mm_cmplt_epi32(vTest,g_XMNoFraction);
    // Convert to int and back to float for rounding
    XMVECTOR vResult = _mm_add_ps(V,g_XMOneHalfMinusEpsilon);
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Convert back to floats
    vResult = _mm_cvtepi32_ps(vInt);
    // All numbers less than 8388608 will use the round to int
    vResult = _mm_and_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    // All others, use the ORIGINAL value
    vTest = _mm_andnot_si128(vTest,_mm_castps_si128(V));
    vResult = _mm_or_ps(vResult,reinterpret_cast<const XMVECTOR *>(&vTest)[0]);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorClamp
(
    FXMVECTOR V, 
    FXMVECTOR Min, 
    FXMVECTOR Max
)
{
    assert(XMVector4LessOrEqual(Min, Max));

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVectorMax(Min, V);
    Result = XMVectorMin(Max, Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR vResult;
    vResult = vmaxq_f32(Min,V);
    vResult = vminq_f32(vResult,Max);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult;
    vResult = _mm_max_ps(Min,V);
    vResult = _mm_min_ps(vResult,Max);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSaturate
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    const XMVECTOR Zero = XMVectorZero();

    return XMVectorClamp(V, Zero, g_XMOne.v);

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Set <0 to 0
    XMVECTOR vResult = vmaxq_f32(V, vdupq_n_u32(0) );
    // Set>1 to 1
    return vminq_f32(vResult, vdupq_n_f32(1.0f) );
#elif defined(_XM_SSE_INTRINSICS_)
    // Set <0 to 0
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    // Set>1 to 1
    return _mm_min_ps(vResult,g_XMOne);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Bitwise logical operations
//------------------------------------------------------------------------------

inline XMVECTOR XMVectorAndInt
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = V1.vector4_u32[0] & V2.vector4_u32[0];
    Result.vector4_u32[1] = V1.vector4_u32[1] & V2.vector4_u32[1];
    Result.vector4_u32[2] = V1.vector4_u32[2] & V2.vector4_u32[2];
    Result.vector4_u32[3] = V1.vector4_u32[3] & V2.vector4_u32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vandq_u32(V1,V2);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_and_ps(V1,V2);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorAndCInt
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = V1.vector4_u32[0] & ~V2.vector4_u32[0];
    Result.vector4_u32[1] = V1.vector4_u32[1] & ~V2.vector4_u32[1];
    Result.vector4_u32[2] = V1.vector4_u32[2] & ~V2.vector4_u32[2];
    Result.vector4_u32[3] = V1.vector4_u32[3] & ~V2.vector4_u32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vbicq_u32(V1,V2);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_andnot_si128( _mm_castps_si128(V2), _mm_castps_si128(V1) );
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorOrInt
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = V1.vector4_u32[0] | V2.vector4_u32[0];
    Result.vector4_u32[1] = V1.vector4_u32[1] | V2.vector4_u32[1];
    Result.vector4_u32[2] = V1.vector4_u32[2] | V2.vector4_u32[2];
    Result.vector4_u32[3] = V1.vector4_u32[3] | V2.vector4_u32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vorrq_u32(V1,V2);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_or_si128( _mm_castps_si128(V1), _mm_castps_si128(V2) );
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorNorInt
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = ~(V1.vector4_u32[0] | V2.vector4_u32[0]);
    Result.vector4_u32[1] = ~(V1.vector4_u32[1] | V2.vector4_u32[1]);
    Result.vector4_u32[2] = ~(V1.vector4_u32[2] | V2.vector4_u32[2]);
    Result.vector4_u32[3] = ~(V1.vector4_u32[3] | V2.vector4_u32[3]);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 Result = vorrq_u32(V1,V2);
    return vbicq_u32(g_XMNegOneMask, Result);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i Result;
    Result = _mm_or_si128( _mm_castps_si128(V1), _mm_castps_si128(V2) );
    Result = _mm_andnot_si128( Result,g_XMNegOneMask);
    return reinterpret_cast<__m128 *>(&Result)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorXorInt
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_u32[0] = V1.vector4_u32[0] ^ V2.vector4_u32[0];
    Result.vector4_u32[1] = V1.vector4_u32[1] ^ V2.vector4_u32[1];
    Result.vector4_u32[2] = V1.vector4_u32[2] ^ V2.vector4_u32[2];
    Result.vector4_u32[3] = V1.vector4_u32[3] ^ V2.vector4_u32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return veorq_u32(V1,V2);
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i V = _mm_xor_si128( _mm_castps_si128(V1), _mm_castps_si128(V2) );
    return reinterpret_cast<__m128 *>(&V)[0];
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Computation operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorNegate
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = -V.vector4_f32[0];
    Result.vector4_f32[1] = -V.vector4_f32[1];
    Result.vector4_f32[2] = -V.vector4_f32[2];
    Result.vector4_f32[3] = -V.vector4_f32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vnegq_f32(V);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR Z;

    Z = _mm_setzero_ps();

    return _mm_sub_ps( Z, V );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorAdd
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = V1.vector4_f32[0] + V2.vector4_f32[0];
    Result.vector4_f32[1] = V1.vector4_f32[1] + V2.vector4_f32[1];
    Result.vector4_f32[2] = V1.vector4_f32[2] + V2.vector4_f32[2];
    Result.vector4_f32[3] = V1.vector4_f32[3] + V2.vector4_f32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vaddq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_add_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorAddAngles
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    const XMVECTOR Zero = XMVectorZero();

    // Add the given angles together.  If the range of V1 is such
    // that -Pi <= V1 < Pi and the range of V2 is such that
    // -2Pi <= V2 <= 2Pi, then the range of the resulting angle
    // will be -Pi <= Result < Pi.
    XMVECTOR Result = XMVectorAdd(V1, V2);

    XMVECTOR Mask = XMVectorLess(Result, g_XMNegativePi.v);
    XMVECTOR Offset = XMVectorSelect(Zero, g_XMTwoPi.v, Mask);

    Mask = XMVectorGreaterOrEqual(Result, g_XMPi.v);
    Offset = XMVectorSelect(Offset, g_XMNegativeTwoPi.v, Mask);

    Result = XMVectorAdd(Result, Offset);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Adjust the angles
    __n128 vResult = vaddq_f32(V1,V2);
    // Less than Pi?
    __n128 vOffset = vcltq_f32(vResult,g_XMNegativePi);
    vOffset = vandq_u32(vOffset,g_XMTwoPi);
    // Add 2Pi to all entries less than -Pi
    vResult = vaddq_f32(vResult,vOffset);
    // Greater than or equal to Pi?
    vOffset = vcgeq_f32(vResult,g_XMPi);
    vOffset = vandq_u32(vOffset,g_XMTwoPi);
    // Sub 2Pi to all entries greater than Pi
    vResult = vsubq_f32(vResult,vOffset);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Adjust the angles
    XMVECTOR vResult = _mm_add_ps(V1,V2);
    // Less than Pi?
    XMVECTOR vOffset = _mm_cmplt_ps(vResult,g_XMNegativePi);
    vOffset = _mm_and_ps(vOffset,g_XMTwoPi);
    // Add 2Pi to all entries less than -Pi
    vResult = _mm_add_ps(vResult,vOffset);
    // Greater than or equal to Pi?
    vOffset = _mm_cmpge_ps(vResult,g_XMPi);
    vOffset = _mm_and_ps(vOffset,g_XMTwoPi);
    // Sub 2Pi to all entries greater than Pi
    vResult = _mm_sub_ps(vResult,vOffset);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSubtract
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = V1.vector4_f32[0] - V2.vector4_f32[0];
    Result.vector4_f32[1] = V1.vector4_f32[1] - V2.vector4_f32[1];
    Result.vector4_f32[2] = V1.vector4_f32[2] - V2.vector4_f32[2];
    Result.vector4_f32[3] = V1.vector4_f32[3] - V2.vector4_f32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vsubq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_sub_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSubtractAngles
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    const XMVECTOR Zero = XMVectorZero();

    // Subtract the given angles.  If the range of V1 is such
    // that -Pi <= V1 < Pi and the range of V2 is such that
    // -2Pi <= V2 <= 2Pi, then the range of the resulting angle
    // will be -Pi <= Result < Pi.
    XMVECTOR Result = XMVectorSubtract(V1, V2);

    XMVECTOR Mask = XMVectorLess(Result, g_XMNegativePi.v);
    XMVECTOR Offset = XMVectorSelect(Zero, g_XMTwoPi.v, Mask);

    Mask = XMVectorGreaterOrEqual(Result, g_XMPi.v);
    Offset = XMVectorSelect(Offset, g_XMNegativeTwoPi.v, Mask);

    Result = XMVectorAdd(Result, Offset);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Adjust the angles
    __n128 vResult = vsubq_f32(V1,V2);
    // Less than Pi?
    __n128 vOffset = vcltq_f32(vResult,g_XMNegativePi);
    vOffset = vandq_u32(vOffset,g_XMTwoPi);
    // Add 2Pi to all entries less than -Pi
    vResult = vaddq_f32(vResult,vOffset);
    // Greater than or equal to Pi?
    vOffset = vcgeq_f32(vResult,g_XMPi);
    vOffset = vandq_u32(vOffset,g_XMTwoPi);
    // Sub 2Pi to all entries greater than Pi
    vResult = vsubq_f32(vResult,vOffset);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Adjust the angles
    XMVECTOR vResult = _mm_sub_ps(V1,V2);
    // Less than Pi?
    XMVECTOR vOffset = _mm_cmplt_ps(vResult,g_XMNegativePi);
    vOffset = _mm_and_ps(vOffset,g_XMTwoPi);
    // Add 2Pi to all entries less than -Pi
    vResult = _mm_add_ps(vResult,vOffset);
    // Greater than or equal to Pi?
    vOffset = _mm_cmpge_ps(vResult,g_XMPi);
    vOffset = _mm_and_ps(vOffset,g_XMTwoPi);
    // Sub 2Pi to all entries greater than Pi
    vResult = _mm_sub_ps(vResult,vOffset);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMultiply
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result = {
        V1.vector4_f32[0] * V2.vector4_f32[0],
        V1.vector4_f32[1] * V2.vector4_f32[1],
        V1.vector4_f32[2] * V2.vector4_f32[2],
        V1.vector4_f32[3] * V2.vector4_f32[3]
    };
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmulq_f32( V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_mul_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMultiplyAdd
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR V3
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {
        (V1.vector4_f32[0] * V2.vector4_f32[0]) + V3.vector4_f32[0],
        (V1.vector4_f32[1] * V2.vector4_f32[1]) + V3.vector4_f32[1],
        (V1.vector4_f32[2] * V2.vector4_f32[2]) + V3.vector4_f32[2],
        (V1.vector4_f32[3] * V2.vector4_f32[3]) + V3.vector4_f32[3]
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmlaq_f32( V3, V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_mul_ps( V1, V2 );
    return _mm_add_ps(vResult, V3 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorDivide
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = V1.vector4_f32[0] / V2.vector4_f32[0];
    Result.vector4_f32[1] = V1.vector4_f32[1] / V2.vector4_f32[1];
    Result.vector4_f32[2] = V1.vector4_f32[2] / V2.vector4_f32[2];
    Result.vector4_f32[3] = V1.vector4_f32[3] / V2.vector4_f32[3];
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // 2 iterations of Newton-Raphson refinement of reciprocal
    __n128 Reciprocal = vrecpeq_f32(V2);
    __n128 S = vrecpsq_f32( Reciprocal, V2 );
    Reciprocal = vmulq_f32( S, Reciprocal );
    S = vrecpsq_f32( Reciprocal, V2 );
    Reciprocal = vmulq_f32( S, Reciprocal );
    return vmulq_f32( V1, Reciprocal );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_div_ps( V1, V2 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorNegativeMultiplySubtract
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR V3
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR vResult = {
        V3.vector4_f32[0] - (V1.vector4_f32[0] * V2.vector4_f32[0]),
        V3.vector4_f32[1] - (V1.vector4_f32[1] * V2.vector4_f32[1]),
        V3.vector4_f32[2] - (V1.vector4_f32[2] * V2.vector4_f32[2]),
        V3.vector4_f32[3] - (V1.vector4_f32[3] * V2.vector4_f32[3])
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmlsq_f32( V3, V1, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR R = _mm_mul_ps( V1, V2 );
    return _mm_sub_ps( V3, R );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorScale
(
    FXMVECTOR V, 
    float    ScaleFactor
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {
        V.vector4_f32[0] * ScaleFactor,
        V.vector4_f32[1] * ScaleFactor,
        V.vector4_f32[2] * ScaleFactor,
        V.vector4_f32[3] * ScaleFactor
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vmulq_n_f32( V, ScaleFactor );
#elif defined(_XM_SSE_INTRINSICS_)
   XMVECTOR vResult = _mm_set_ps1(ScaleFactor);
   return _mm_mul_ps(vResult,V);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorReciprocalEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = 1.f / V.vector4_f32[0];
    Result.vector4_f32[1] = 1.f / V.vector4_f32[1];
    Result.vector4_f32[2] = 1.f / V.vector4_f32[2];
    Result.vector4_f32[3] = 1.f / V.vector4_f32[3];
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vrecpeq_f32(V);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_rcp_ps(V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorReciprocal
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = 1.f / V.vector4_f32[0];
    Result.vector4_f32[1] = 1.f / V.vector4_f32[1];
    Result.vector4_f32[2] = 1.f / V.vector4_f32[2];
    Result.vector4_f32[3] = 1.f / V.vector4_f32[3];
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // 2 iterations of Newton-Raphson refinement
    __n128 Reciprocal = vrecpeq_f32(V);
    __n128 S = vrecpsq_f32( Reciprocal, V );
    Reciprocal = vmulq_f32( S, Reciprocal );
    S = vrecpsq_f32( Reciprocal, V );
    return vmulq_f32( S, Reciprocal );
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_div_ps(g_XMOne,V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Return an estimated square root
inline XMVECTOR XMVectorSqrtEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = sqrtf( V.vector4_f32[0] );
    Result.vector4_f32[1] = sqrtf( V.vector4_f32[1] );
    Result.vector4_f32[2] = sqrtf( V.vector4_f32[2] );
    Result.vector4_f32[3] = sqrtf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // 1 iteration of Newton-Raphson refinment of sqrt
    __n128 S0 = vrsqrteq_f32(V);
    __n128 P0 = vmulq_f32( V, S0 );
    __n128 R0 = vrsqrtsq_f32( P0, S0 );
    __n128 S1 = vmulq_f32( S0, R0 );

    XMVECTOR VEqualsInfinity = XMVectorEqualInt(V, g_XMInfinity.v);
    XMVECTOR VEqualsZero = XMVectorEqual(V, vdupq_n_f32(0) );
    __n128 Result = vmulq_f32( V, S1 );
    XMVECTOR Select = XMVectorEqualInt(VEqualsInfinity, VEqualsZero);
    return XMVectorSelect(V, Result, Select);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_sqrt_ps(V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSqrt
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = sqrtf( V.vector4_f32[0] );
    Result.vector4_f32[1] = sqrtf( V.vector4_f32[1] );
    Result.vector4_f32[2] = sqrtf( V.vector4_f32[2] );
    Result.vector4_f32[3] = sqrtf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // 3 iterations of Newton-Raphson refinment of sqrt
    __n128 S0 = vrsqrteq_f32(V);
    __n128 P0 = vmulq_f32( V, S0 );
    __n128 R0 = vrsqrtsq_f32( P0, S0 );
    __n128 S1 = vmulq_f32( S0, R0 );
    __n128 P1 = vmulq_f32( V, S1 );
    __n128 R1 = vrsqrtsq_f32( P1, S1 );
    __n128 S2 = vmulq_f32( S1, R1 );
    __n128 P2 = vmulq_f32( V, S2 );
    __n128 R2 = vrsqrtsq_f32( P2, S2 );
    __n128 S3 = vmulq_f32( S2, R2 );

    XMVECTOR VEqualsInfinity = XMVectorEqualInt(V, g_XMInfinity.v);
    XMVECTOR VEqualsZero = XMVectorEqual(V, vdupq_n_f32(0) );
    __n128 Result = vmulq_f32( V, S3 );
    XMVECTOR Select = XMVectorEqualInt(VEqualsInfinity, VEqualsZero);
    return XMVectorSelect(V, Result, Select);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_sqrt_ps(V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorReciprocalSqrtEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = 1.f / sqrtf( V.vector4_f32[0] );
    Result.vector4_f32[1] = 1.f / sqrtf( V.vector4_f32[1] );
    Result.vector4_f32[2] = 1.f / sqrtf( V.vector4_f32[2] );
    Result.vector4_f32[3] = 1.f / sqrtf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vrsqrteq_f32(V);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_rsqrt_ps(V);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorReciprocalSqrt
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = 1.f / sqrtf( V.vector4_f32[0] );
    Result.vector4_f32[1] = 1.f / sqrtf( V.vector4_f32[1] );
    Result.vector4_f32[2] = 1.f / sqrtf( V.vector4_f32[2] );
    Result.vector4_f32[3] = 1.f / sqrtf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // 2 iterations of Newton-Raphson refinement of reciprocal
    __n128 S0 = vrsqrteq_f32(V);

    __n128 P0 = vmulq_f32( V, S0 );
    __n128 R0 = vrsqrtsq_f32( P0, S0 );

    __n128 S1 = vmulq_f32( S0, R0 );
    __n128 P1 = vmulq_f32( V, S1 );
    __n128 R1 = vrsqrtsq_f32( P1, S1 );

    return vmulq_f32( S1, R1 );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_sqrt_ps(V);
    vResult = _mm_div_ps(g_XMOne,vResult);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline XMVECTOR XMVectorExp
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = powf(2.0f, V.vector4_f32[0]);
    Result.vector4_f32[1] = powf(2.0f, V.vector4_f32[1]);
    Result.vector4_f32[2] = powf(2.0f, V.vector4_f32[2]);
    Result.vector4_f32[3] = powf(2.0f, V.vector4_f32[3]);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        powf(2.0f,vgetq_lane_f32(V, 0)),
        powf(2.0f,vgetq_lane_f32(V, 1)),
        powf(2.0f,vgetq_lane_f32(V, 2)),
        powf(2.0f,vgetq_lane_f32(V, 3))
    };
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    __declspec(align(16)) float a[4];
    _mm_store_ps( a, V );
    XMVECTOR vResult = _mm_setr_ps(
        powf(2.0f,a[0]),
        powf(2.0f,a[1]),
        powf(2.0f,a[2]),
        powf(2.0f,a[3]));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline XMVECTOR XMVectorLog
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    const float fScale = 1.4426950f; // (1.0f / logf(2.0f));

    XMVECTOR Result;
    Result.vector4_f32[0] = logf(V.vector4_f32[0])*fScale;
    Result.vector4_f32[1] = logf(V.vector4_f32[1])*fScale;
    Result.vector4_f32[2] = logf(V.vector4_f32[2])*fScale;
    Result.vector4_f32[3] = logf(V.vector4_f32[3])*fScale;
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR vScale = vdupq_n_f32(1.0f / logf(2.0f));
    XMVECTORF32 vResult = {
        logf(vgetq_lane_f32(V, 0)),
        logf(vgetq_lane_f32(V, 1)),
        logf(vgetq_lane_f32(V, 2)),
        logf(vgetq_lane_f32(V, 3))
    };
    return vmulq_f32( vResult, vScale );
#elif defined(_XM_SSE_INTRINSICS_)
    __declspec(align(16)) float a[4];
    _mm_store_ps( a, V );
    XMVECTOR vScale = _mm_set_ps1(1.0f / logf(2.0f));
    XMVECTOR vResult = _mm_setr_ps(
        logf(a[0]),
        logf(a[1]),
        logf(a[2]),
        logf(a[3]));
    vResult = _mm_mul_ps(vResult,vScale);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline XMVECTOR XMVectorPow
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = powf(V1.vector4_f32[0], V2.vector4_f32[0]);
    Result.vector4_f32[1] = powf(V1.vector4_f32[1], V2.vector4_f32[1]);
    Result.vector4_f32[2] = powf(V1.vector4_f32[2], V2.vector4_f32[2]);
    Result.vector4_f32[3] = powf(V1.vector4_f32[3], V2.vector4_f32[3]);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        powf(vgetq_lane_f32(V1, 0), vgetq_lane_f32(V2, 0)),
        powf(vgetq_lane_f32(V1, 1), vgetq_lane_f32(V2, 1)),
        powf(vgetq_lane_f32(V1, 2), vgetq_lane_f32(V2, 2)),
        powf(vgetq_lane_f32(V1, 3), vgetq_lane_f32(V2, 3))
    };
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    __declspec(align(16)) float a[4];
    __declspec(align(16)) float b[4];
    _mm_store_ps( a, V1 );
    _mm_store_ps( b, V2 );
    XMVECTOR vResult = _mm_setr_ps(
        powf(a[0],b[0]),
        powf(a[1],b[1]),
        powf(a[2],b[2]),
        powf(a[3],b[3]));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorAbs
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {
        fabsf(V.vector4_f32[0]),
        fabsf(V.vector4_f32[1]),
        fabsf(V.vector4_f32[2]),
        fabsf(V.vector4_f32[3])
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vabsq_f32( V );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_setzero_ps();
    vResult = _mm_sub_ps(vResult,V);
    vResult = _mm_max_ps(vResult,V);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMod
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    // V1 % V2 = V1 - V2 * truncate(V1 / V2)

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Quotient = XMVectorDivide(V1, V2);
    Quotient = XMVectorTruncate(Quotient);
    XMVECTOR Result = XMVectorNegativeMultiplySubtract(V2, Quotient, V1);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR vResult = XMVectorDivide(V1, V2);
    vResult = XMVectorTruncate(vResult);
    return vmlsq_f32( V1, vResult, V2 );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_div_ps(V1, V2);
    vResult = XMVectorTruncate(vResult);
    vResult = _mm_mul_ps(vResult,V2);
    vResult = _mm_sub_ps(V1,vResult);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorModAngles
(
    FXMVECTOR Angles
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;
    XMVECTOR Result;

    // Modulo the range of the given angles such that -XM_PI <= Angles < XM_PI
    V = XMVectorMultiply(Angles, g_XMReciprocalTwoPi.v);
    V = XMVectorRound(V);
    Result = XMVectorNegativeMultiplySubtract(g_XMTwoPi.v, V, Angles);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Modulo the range of the given angles such that -XM_PI <= Angles < XM_PI
    XMVECTOR vResult = vmulq_f32(Angles,g_XMReciprocalTwoPi);
    // Use the inline function due to complexity for rounding
    vResult = XMVectorRound(vResult);
    return vmlsq_f32( Angles, vResult, g_XMTwoPi );
#elif defined(_XM_SSE_INTRINSICS_)
    // Modulo the range of the given angles such that -XM_PI <= Angles < XM_PI
    XMVECTOR vResult = _mm_mul_ps(Angles,g_XMReciprocalTwoPi);
    // Use the inline function due to complexity for rounding
    vResult = XMVectorRound(vResult);
    vResult = _mm_mul_ps(vResult,g_XMTwoPi);
    vResult = _mm_sub_ps(Angles,vResult);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSin
(
    FXMVECTOR V
)
{
    // 11-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarSin( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarSin( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarSin( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarSin( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with sin(y) = sin(x).
    __n128 sign = vandq_u32(x, g_XMNegativeZero);
    __n128 c = vorrq_u32(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __n128 absx = vabsq_f32( x );
    __n128 rflx = vsubq_f32(c, x);
    __n128 comp = vcleq_f32(absx, g_XMHalfPi);
    x = vbslq_f32( comp, x, rflx );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation
    const XMVECTOR SC1 = g_XMSinCoefficients1;
    XMVECTOR Result = vdupq_lane_f32(vget_low_f32(SC1), 0);

    const XMVECTOR SC0 = g_XMSinCoefficients0;
    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(SC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_high_f32(SC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(SC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(SC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    Result = vmulq_f32(Result, x);
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with sin(y) = sin(x).
    __m128 sign = _mm_and_ps(x, g_XMNegativeZero);
    __m128 c = _mm_or_ps(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, g_XMHalfPi);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    x = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation
    const XMVECTOR SC1 = g_XMSinCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( SC1, _MM_SHUFFLE(0, 0, 0, 0) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    const XMVECTOR SC0 = g_XMSinCoefficients0;
    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(3, 3, 3, 3) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SC0,  _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);
    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, x);
    return Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorCos
(
    FXMVECTOR V
)
{
    // 10-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarCos( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarCos( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarCos( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarCos( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Map V to x in [-pi,pi].
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    __n128 sign = vandq_u32(x, g_XMNegativeZero);
    __n128 c = vorrq_u32(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __n128 absx = vabsq_f32( x );
    __n128 rflx = vsubq_f32(c, x);
    __n128 comp = vcleq_f32(absx, g_XMHalfPi);
    x = vbslq_f32( comp, x, rflx );
    sign = vbslq_f32( comp, g_XMOne, g_XMNegativeOne );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation
    const XMVECTOR CC1 = g_XMCosCoefficients1;
    XMVECTOR Result = vdupq_lane_f32(vget_low_f32(CC1), 0);

    const XMVECTOR CC0 = g_XMCosCoefficients0;
    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(CC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_high_f32(CC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(CC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(CC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    Result = vmulq_f32(Result, sign);
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    // Map V to x in [-pi,pi].
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    XMVECTOR sign = _mm_and_ps(x, g_XMNegativeZero);
    __m128 c = _mm_or_ps(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, g_XMHalfPi);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    x = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, g_XMOne);
    select1 = _mm_andnot_ps(comp, g_XMNegativeOne);
    sign = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation
    const XMVECTOR CC1 = g_XMCosCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( CC1, _MM_SHUFFLE(0, 0, 0, 0) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    const XMVECTOR CC0 = g_XMCosCoefficients0;
    vConstants = XM_PERMUTE_PS( CC0, _MM_SHUFFLE(3, 3, 3, 3) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CC0, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CC0, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CC0, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);
    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, sign);
    return Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline void XMVectorSinCos
(
    XMVECTOR* pSin, 
    XMVECTOR* pCos, 
    FXMVECTOR V
)
{
    assert(pSin != nullptr);
    assert(pCos != nullptr);

    // 11/10-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Sin;
    XMVECTOR Cos;

    XMScalarSinCos(&Sin.vector4_f32[0], &Cos.vector4_f32[0], V.vector4_f32[0]);
    XMScalarSinCos(&Sin.vector4_f32[1], &Cos.vector4_f32[1], V.vector4_f32[1]);
    XMScalarSinCos(&Sin.vector4_f32[2], &Cos.vector4_f32[2], V.vector4_f32[2]);
    XMScalarSinCos(&Sin.vector4_f32[3], &Cos.vector4_f32[3], V.vector4_f32[3]);

    *pSin = Sin;
    *pCos = Cos;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    __n128 sign = vandq_u32(x, g_XMNegativeZero);
    __n128 c = vorrq_u32(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __n128 absx = vabsq_f32( x );
    __n128 rflx = vsubq_f32(c, x);
    __n128 comp = vcleq_f32(absx, g_XMHalfPi);
    x = vbslq_f32( comp, x, rflx );
    sign = vbslq_f32( comp, g_XMOne, g_XMNegativeOne );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation for sine
    const XMVECTOR SC1 = g_XMSinCoefficients1;
    XMVECTOR Result = vdupq_lane_f32(vget_low_f32(SC1), 0);

    const XMVECTOR SC0 = g_XMSinCoefficients0;
    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(SC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_high_f32(SC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(SC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(SC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    *pSin = vmulq_f32(Result, x);

    // Compute polynomial approximation for cosine
    const XMVECTOR CC1 = g_XMCosCoefficients1;
    Result = vdupq_lane_f32(vget_low_f32(CC1), 0);

    const XMVECTOR CC0 = g_XMCosCoefficients0;
    vConstants = vdupq_lane_f32(vget_high_f32(CC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_high_f32(CC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(CC0), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(CC0), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    *pCos = vmulq_f32(Result, sign);
#elif defined(_XM_SSE_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with sin(y) = sin(x), cos(y) = sign*cos(x).
    XMVECTOR sign = _mm_and_ps(x, g_XMNegativeZero);
    __m128 c = _mm_or_ps(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, g_XMHalfPi);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    x = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, g_XMOne);
    select1 = _mm_andnot_ps(comp, g_XMNegativeOne);
    sign = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation of sine
    const XMVECTOR SC1 = g_XMSinCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( SC1, _MM_SHUFFLE(0, 0, 0, 0) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    const XMVECTOR SC0 = g_XMSinCoefficients0;
    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(3, 3, 3, 3) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SC0, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);
    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, x);
    *pSin = Result;

    // Compute polynomial approximation of cosine
    const XMVECTOR CC1 = g_XMCosCoefficients1;
    vConstants = XM_PERMUTE_PS( CC1, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_mul_ps(vConstants, x2);

    const XMVECTOR CC0 = g_XMCosCoefficients0;
    vConstants = XM_PERMUTE_PS( CC0, _MM_SHUFFLE(3, 3, 3, 3) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CC0,  _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CC0,  _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CC0, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);
    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, sign);
    *pCos = Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorTan
(
    FXMVECTOR V
)
{
    // Cody and Waite algorithm to compute tangent.

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = tanf( V.vector4_f32[0] );
    Result.vector4_f32[1] = tanf( V.vector4_f32[1] );
    Result.vector4_f32[2] = tanf( V.vector4_f32[2] );
    Result.vector4_f32[3] = tanf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_) 

    static const XMVECTORF32 TanCoefficients0 = {1.0f, -4.667168334e-1f, 2.566383229e-2f, -3.118153191e-4f};
    static const XMVECTORF32 TanCoefficients1 = {4.981943399e-7f, -1.333835001e-1f, 3.424887824e-3f, -1.786170734e-5f};
    static const XMVECTORF32 TanConstants = {1.570796371f, 6.077100628e-11f, 0.000244140625f, 0.63661977228f /*2 / Pi*/ };
    static const XMVECTORU32 Mask = {0x1, 0x1, 0x1, 0x1};

    XMVECTOR TwoDivPi = XMVectorSplatW(TanConstants.v);

    XMVECTOR Zero = XMVectorZero();

    XMVECTOR C0 = XMVectorSplatX(TanConstants.v);
    XMVECTOR C1 = XMVectorSplatY(TanConstants.v);
    XMVECTOR Epsilon = XMVectorSplatZ(TanConstants.v);

    XMVECTOR VA = XMVectorMultiply(V, TwoDivPi);

    VA = XMVectorRound(VA);

    XMVECTOR VC = XMVectorNegativeMultiplySubtract(VA, C0, V);

    XMVECTOR VB = XMVectorAbs(VA);

    VC = XMVectorNegativeMultiplySubtract(VA, C1, VC);

#if defined(_XM_ARM_NEON_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    VB = vcvtq_u32_f32( VB );
#elif defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    reinterpret_cast<__m128i *>(&VB)[0] = _mm_cvttps_epi32(VB);
#else
    for (size_t i = 0; i < 4; i++)
    {
        VB.vector4_u32[i] = (uint32_t)VB.vector4_f32[i];
    }
#endif

    XMVECTOR VC2 = XMVectorMultiply(VC, VC);

    XMVECTOR T7 = XMVectorSplatW(TanCoefficients1.v);
    XMVECTOR T6 = XMVectorSplatZ(TanCoefficients1.v);
    XMVECTOR T4 = XMVectorSplatX(TanCoefficients1.v);
    XMVECTOR T3 = XMVectorSplatW(TanCoefficients0.v);
    XMVECTOR T5 = XMVectorSplatY(TanCoefficients1.v);
    XMVECTOR T2 = XMVectorSplatZ(TanCoefficients0.v);
    XMVECTOR T1 = XMVectorSplatY(TanCoefficients0.v);
    XMVECTOR T0 = XMVectorSplatX(TanCoefficients0.v);

    XMVECTOR VBIsEven = XMVectorAndInt(VB, Mask.v);
    VBIsEven = XMVectorEqualInt(VBIsEven, Zero);

    XMVECTOR N = XMVectorMultiplyAdd(VC2, T7, T6);
    XMVECTOR D = XMVectorMultiplyAdd(VC2, T4, T3);
    N = XMVectorMultiplyAdd(VC2, N, T5);
    D = XMVectorMultiplyAdd(VC2, D, T2);
    N = XMVectorMultiply(VC2, N);
    D = XMVectorMultiplyAdd(VC2, D, T1);
    N = XMVectorMultiplyAdd(VC, N, VC);
    XMVECTOR VCNearZero = XMVectorInBounds(VC, Epsilon);
    D = XMVectorMultiplyAdd(VC2, D, T0);

    N = XMVectorSelect(N, VC, VCNearZero);
    D = XMVectorSelect(D, g_XMOne.v, VCNearZero);

    XMVECTOR R0 = XMVectorNegate(N);
    XMVECTOR R1 = XMVectorDivide(N,D);
    R0 = XMVectorDivide(D,R0);

    XMVECTOR VIsZero = XMVectorEqual(V, Zero);

    XMVECTOR Result = XMVectorSelect(R0, R1, VBIsEven);

    Result = XMVectorSelect(Result, Zero, VIsZero);

    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSinH
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = sinhf( V.vector4_f32[0] );
    Result.vector4_f32[1] = sinhf( V.vector4_f32[1] );
    Result.vector4_f32[2] = sinhf( V.vector4_f32[2] );
    Result.vector4_f32[3] = sinhf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Scale = {1.442695040888963f, 1.442695040888963f, 1.442695040888963f, 1.442695040888963f}; // 1.0f / ln(2.0f)

    XMVECTOR V1 = vmlaq_f32( g_XMNegativeOne.v, V, Scale.v );
    XMVECTOR V2 = vmlsq_f32( g_XMNegativeOne.v, V, Scale.v );
    XMVECTOR E1 = XMVectorExp(V1);
    XMVECTOR E2 = XMVectorExp(V2);

    return vsubq_f32(E1, E2);
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Scale = {1.442695040888963f, 1.442695040888963f, 1.442695040888963f, 1.442695040888963f}; // 1.0f / ln(2.0f)

    XMVECTOR V1 = _mm_mul_ps(V, Scale);
    V1 = _mm_add_ps(V1,g_XMNegativeOne);
    XMVECTOR V2 = _mm_mul_ps(V, Scale);
    V2 = _mm_sub_ps(g_XMNegativeOne,V2);
    XMVECTOR E1 = XMVectorExp(V1);
    XMVECTOR E2 = XMVectorExp(V2);

    return _mm_sub_ps(E1, E2);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorCosH
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = coshf( V.vector4_f32[0] );
    Result.vector4_f32[1] = coshf( V.vector4_f32[1] );
    Result.vector4_f32[2] = coshf( V.vector4_f32[2] );
    Result.vector4_f32[3] = coshf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Scale = {1.442695040888963f, 1.442695040888963f, 1.442695040888963f, 1.442695040888963f}; // 1.0f / ln(2.0f)

    XMVECTOR V1 = vmlaq_f32(g_XMNegativeOne.v, V, Scale.v);
    XMVECTOR V2 = vmlsq_f32(g_XMNegativeOne.v, V, Scale.v);
    XMVECTOR E1 = XMVectorExp(V1);
    XMVECTOR E2 = XMVectorExp(V2);
    return vaddq_f32(E1, E2);
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Scale = {1.442695040888963f, 1.442695040888963f, 1.442695040888963f, 1.442695040888963f}; // 1.0f / ln(2.0f)

    XMVECTOR V1 = _mm_mul_ps(V,Scale.v);
    V1 = _mm_add_ps(V1,g_XMNegativeOne.v);
    XMVECTOR V2 = _mm_mul_ps(V, Scale.v);
    V2 = _mm_sub_ps(g_XMNegativeOne.v,V2);
    XMVECTOR E1 = XMVectorExp(V1);
    XMVECTOR E2 = XMVectorExp(V2);
    return _mm_add_ps(E1, E2);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorTanH
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = tanhf( V.vector4_f32[0] );
    Result.vector4_f32[1] = tanhf( V.vector4_f32[1] );
    Result.vector4_f32[2] = tanhf( V.vector4_f32[2] );
    Result.vector4_f32[3] = tanhf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Scale = {2.8853900817779268f, 2.8853900817779268f, 2.8853900817779268f, 2.8853900817779268f}; // 2.0f / ln(2.0f)

    XMVECTOR E = vmulq_f32(V, Scale.v);
    E = XMVectorExp(E);
    E = vmlaq_f32( g_XMOneHalf.v, E, g_XMOneHalf.v );
    E = XMVectorReciprocal(E);
    return vsubq_f32(g_XMOne.v, E);
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Scale = {2.8853900817779268f, 2.8853900817779268f, 2.8853900817779268f, 2.8853900817779268f}; // 2.0f / ln(2.0f)

    XMVECTOR E = _mm_mul_ps(V, Scale.v);
    E = XMVectorExp(E);
    E = _mm_mul_ps(E,g_XMOneHalf.v);
    E = _mm_add_ps(E,g_XMOneHalf.v);
    E = _mm_div_ps(g_XMOne.v,E);
    return _mm_sub_ps(g_XMOne.v,E);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorASin
(
    FXMVECTOR V
)
{
    // 7-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarASin( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarASin( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarASin( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarASin( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 nonnegative = vcgeq_f32(V, g_XMZero);
    __n128 x = vabsq_f32(V);

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __n128 oneMValue = vsubq_f32(g_XMOne, x);
    __n128 clampOneMValue = vmaxq_f32(g_XMZero, oneMValue);
    __n128 root = XMVectorSqrt(clampOneMValue);

    // Compute polynomial approximation
    const XMVECTOR AC1 = g_XMArcCoefficients1;
    __n128 t0 = vdupq_lane_f32(vget_high_f32(AC1), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(AC1), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC1), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC1), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    const XMVECTOR AC0 = g_XMArcCoefficients0;
    vConstants = vdupq_lane_f32(vget_high_f32(AC0), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_high_f32(AC0), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC0), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC0), 0);
    t0 = vmlaq_f32( vConstants, t0, x );
    t0 = vmulq_f32(t0, root);

    __n128 t1 = vsubq_f32(g_XMPi, t0);
    t0 = vbslq_f32( nonnegative, t0, t1 );
    t0 = vsubq_f32(g_XMHalfPi, t0);
    return t0;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128 nonnegative = _mm_cmpge_ps(V, g_XMZero);
    __m128 mvalue = _mm_sub_ps(g_XMZero, V);
    __m128 x = _mm_max_ps(V, mvalue);  // |V|

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __m128 oneMValue = _mm_sub_ps(g_XMOne, x);
    __m128 clampOneMValue = _mm_max_ps(g_XMZero, oneMValue);
    __m128 root = _mm_sqrt_ps(clampOneMValue);  // sqrt(1-|V|)

    // Compute polynomial approximation
    const XMVECTOR AC1 = g_XMArcCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 t0 = _mm_mul_ps(vConstants, x);

    vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(2, 2, 2, 2) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(1, 1, 1, 1) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(0, 0, 0, 0) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    const XMVECTOR AC0 = g_XMArcCoefficients0;
    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(3, 3, 3, 3) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC0,_MM_SHUFFLE(2, 2, 2, 2) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(1, 1, 1, 1) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(0, 0, 0, 0) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, root);

    __m128 t1 = _mm_sub_ps(g_XMPi, t0);
    t0 = _mm_and_ps(nonnegative, t0);
    t1 = _mm_andnot_ps(nonnegative, t1);
    t0 = _mm_or_ps(t0, t1);
    t0 = _mm_sub_ps(g_XMHalfPi, t0);
    return t0;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorACos
(
    FXMVECTOR V
)
{
    // 7-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarACos( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarACos( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarACos( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarACos( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 nonnegative = vcgeq_f32(V, g_XMZero);
    __n128 x = vabsq_f32(V);

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __n128 oneMValue = vsubq_f32(g_XMOne, x);
    __n128 clampOneMValue = vmaxq_f32(g_XMZero, oneMValue);
    __n128 root = XMVectorSqrt(clampOneMValue);

    // Compute polynomial approximation
    const XMVECTOR AC1 = g_XMArcCoefficients1;
    __n128 t0 = vdupq_lane_f32(vget_high_f32(AC1), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(AC1), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC1), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC1), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    const XMVECTOR AC0 = g_XMArcCoefficients0;
    vConstants = vdupq_lane_f32(vget_high_f32(AC0), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_high_f32(AC0), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC0), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AC0), 0);
    t0 = vmlaq_f32( vConstants, t0, x );
    t0 = vmulq_f32(t0, root);

    __n128 t1 = vsubq_f32(g_XMPi, t0);
    t0 = vbslq_f32( nonnegative, t0, t1 );
    return t0;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128 nonnegative = _mm_cmpge_ps(V, g_XMZero);
    __m128 mvalue = _mm_sub_ps(g_XMZero, V);
    __m128 x = _mm_max_ps(V, mvalue);  // |V|

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __m128 oneMValue = _mm_sub_ps(g_XMOne, x);
    __m128 clampOneMValue = _mm_max_ps(g_XMZero, oneMValue);
    __m128 root = _mm_sqrt_ps(clampOneMValue);  // sqrt(1-|V|)

    // Compute polynomial approximation
    const XMVECTOR AC1 = g_XMArcCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 t0 = _mm_mul_ps(vConstants, x);

    vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(2, 2, 2, 2) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(1, 1, 1, 1) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC1, _MM_SHUFFLE(0, 0, 0, 0) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    const XMVECTOR AC0 = g_XMArcCoefficients0;
    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(3, 3, 3, 3) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(2, 2, 2, 2) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(1, 1, 1, 1) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AC0, _MM_SHUFFLE(0, 0, 0, 0) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, root);

    __m128 t1 = _mm_sub_ps(g_XMPi, t0);
    t0 = _mm_and_ps(nonnegative, t0);
    t1 = _mm_andnot_ps(nonnegative, t1);
    t0 = _mm_or_ps(t0, t1);
    return t0;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorATan
(
    FXMVECTOR V
)
{
    // 17-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = atanf( V.vector4_f32[0] );
    Result.vector4_f32[1] = atanf( V.vector4_f32[1] );
    Result.vector4_f32[2] = atanf( V.vector4_f32[2] );
    Result.vector4_f32[3] = atanf( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 absV = vabsq_f32(V);
    __n128 invV = XMVectorReciprocal(V);
    __n128 comp = vcgtq_f32(V, g_XMOne);
    __n128 sign = vbslq_f32(comp, g_XMOne, g_XMNegativeOne);
    comp = vcleq_f32(absV, g_XMOne);
    sign = vbslq_f32(comp, g_XMZero, sign);
    __n128 x = vbslq_f32(comp, V, invV);

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation
    const XMVECTOR TC1 = g_XMATanCoefficients1;
    __n128 Result = vdupq_lane_f32(vget_high_f32(TC1), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(TC1), 0);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_low_f32(TC1), 1);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_low_f32(TC1), 0);
    Result = vmlaq_f32( vConstants, Result, x2 );

    const XMVECTOR TC0 = g_XMATanCoefficients0;
    vConstants = vdupq_lane_f32(vget_high_f32(TC0), 1);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_high_f32(TC0), 0);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_low_f32(TC0), 1);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_low_f32(TC0), 0);
    Result = vmlaq_f32( vConstants, Result, x2 );

    Result = vmlaq_f32( g_XMOne, Result, x2 );
    Result = vmulq_f32( Result, x );

    __n128 result1 = vmulq_f32(sign, g_XMHalfPi);
    result1 = vsubq_f32(result1, Result);

    comp = vceqq_f32(sign, g_XMZero);
    Result = vbslq_f32( comp, Result, result1 );
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128 absV = XMVectorAbs(V);
    __m128 invV = _mm_div_ps(g_XMOne, V);
    __m128 comp = _mm_cmpgt_ps(V, g_XMOne);
    __m128 select0 = _mm_and_ps(comp, g_XMOne);
    __m128 select1 = _mm_andnot_ps(comp, g_XMNegativeOne);
    __m128 sign = _mm_or_ps(select0, select1);
    comp = _mm_cmple_ps(absV, g_XMOne);
    select0 = _mm_and_ps(comp, g_XMZero);
    select1 = _mm_andnot_ps(comp, sign);
    sign = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, V);
    select1 = _mm_andnot_ps(comp, invV);
    __m128 x = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation
    const XMVECTOR TC1 = g_XMATanCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( TC1, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    vConstants = XM_PERMUTE_PS( TC1, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( TC1, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( TC1, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    const XMVECTOR TC0 = g_XMATanCoefficients0;
    vConstants = XM_PERMUTE_PS( TC0, _MM_SHUFFLE(3, 3, 3, 3) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( TC0, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( TC0, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( TC0, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);
    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, x);
    __m128 result1 = _mm_mul_ps(sign, g_XMHalfPi);
    result1 = _mm_sub_ps(result1, Result);

    comp = _mm_cmpeq_ps(sign, g_XMZero);
    select0 = _mm_and_ps(comp, Result);
    select1 = _mm_andnot_ps(comp, result1);
    Result = _mm_or_ps(select0, select1);
    return Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorATan2
(
    FXMVECTOR Y, 
    FXMVECTOR X
)
{
    // Return the inverse tangent of Y / X in the range of -Pi to Pi with the following exceptions:

    //     Y == 0 and X is Negative         -> Pi with the sign of Y
    //     y == 0 and x is positive         -> 0 with the sign of y
    //     Y != 0 and X == 0                -> Pi / 2 with the sign of Y
    //     Y != 0 and X is Negative         -> atan(y/x) + (PI with the sign of Y)
    //     X == -Infinity and Finite Y      -> Pi with the sign of Y
    //     X == +Infinity and Finite Y      -> 0 with the sign of Y
    //     Y == Infinity and X is Finite    -> Pi / 2 with the sign of Y
    //     Y == Infinity and X == -Infinity -> 3Pi / 4 with the sign of Y
    //     Y == Infinity and X == +Infinity -> Pi / 4 with the sign of Y

    static const XMVECTORF32 ATan2Constants = {XM_PI, XM_PIDIV2, XM_PIDIV4, XM_PI * 3.0f / 4.0f};

    XMVECTOR Zero = XMVectorZero();
    XMVECTOR ATanResultValid = XMVectorTrueInt();

    XMVECTOR Pi = XMVectorSplatX(ATan2Constants);
    XMVECTOR PiOverTwo = XMVectorSplatY(ATan2Constants);
    XMVECTOR PiOverFour = XMVectorSplatZ(ATan2Constants);
    XMVECTOR ThreePiOverFour = XMVectorSplatW(ATan2Constants);

    XMVECTOR YEqualsZero = XMVectorEqual(Y, Zero);
    XMVECTOR XEqualsZero = XMVectorEqual(X, Zero);
    XMVECTOR XIsPositive = XMVectorAndInt(X, g_XMNegativeZero.v);
    XIsPositive = XMVectorEqualInt(XIsPositive, Zero);
    XMVECTOR YEqualsInfinity = XMVectorIsInfinite(Y);
    XMVECTOR XEqualsInfinity = XMVectorIsInfinite(X);

    XMVECTOR YSign = XMVectorAndInt(Y, g_XMNegativeZero.v);
    Pi = XMVectorOrInt(Pi, YSign);
    PiOverTwo = XMVectorOrInt(PiOverTwo, YSign);
    PiOverFour = XMVectorOrInt(PiOverFour, YSign);
    ThreePiOverFour = XMVectorOrInt(ThreePiOverFour, YSign);

    XMVECTOR R1 = XMVectorSelect(Pi, YSign, XIsPositive);
    XMVECTOR R2 = XMVectorSelect(ATanResultValid, PiOverTwo, XEqualsZero);
    XMVECTOR R3 = XMVectorSelect(R2, R1, YEqualsZero);
    XMVECTOR R4 = XMVectorSelect(ThreePiOverFour, PiOverFour, XIsPositive);
    XMVECTOR R5 = XMVectorSelect(PiOverTwo, R4, XEqualsInfinity);
    XMVECTOR Result = XMVectorSelect(R3, R5, YEqualsInfinity);
    ATanResultValid = XMVectorEqualInt(Result, ATanResultValid);

    XMVECTOR V = XMVectorDivide(Y, X);

    XMVECTOR R0 = XMVectorATan(V);

    R1 = XMVectorSelect( Pi, Zero, XIsPositive );
    R2 = XMVectorAdd(R0, R1);

    return XMVectorSelect(Result, R2, ATanResultValid);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorSinEst
(
    FXMVECTOR V
)
{
    // 7-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarSinEst( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarSinEst( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarSinEst( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarSinEst( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with sin(y) = sin(x).
    __n128 sign = vandq_u32(x, g_XMNegativeZero);
    __n128 c = vorrq_u32(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __n128 absx = vabsq_f32( x );
    __n128 rflx = vsubq_f32(c, x);
    __n128 comp = vcleq_f32(absx, g_XMHalfPi);
    x = vbslq_f32( comp, x, rflx );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation
    const XMVECTOR SEC = g_XMSinCoefficients1;
    XMVECTOR Result = vdupq_lane_f32(vget_high_f32(SEC), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(SEC), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(SEC), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    Result = vmulq_f32(Result, x);
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with sin(y) = sin(x).
    __m128 sign = _mm_and_ps(x, g_XMNegativeZero);
    __m128 c = _mm_or_ps(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, g_XMHalfPi);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    x = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation
    const XMVECTOR SEC = g_XMSinCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( SEC, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    vConstants = XM_PERMUTE_PS( SEC, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SEC, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, x);
    return Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorCosEst
(
    FXMVECTOR V
)
{
    // 6-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarCosEst( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarCosEst( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarCosEst( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarCosEst( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Map V to x in [-pi,pi].
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    __n128 sign = vandq_u32(x, g_XMNegativeZero);
    __n128 c = vorrq_u32(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __n128 absx = vabsq_f32( x );
    __n128 rflx = vsubq_f32(c, x);
    __n128 comp = vcleq_f32(absx, g_XMHalfPi);
    x = vbslq_f32( comp, x, rflx );
    sign = vbslq_f32( comp, g_XMOne, g_XMNegativeOne );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation
    const XMVECTOR CEC = g_XMCosCoefficients1;
    XMVECTOR Result = vdupq_lane_f32(vget_high_f32(CEC), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(CEC), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(CEC), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    Result = vmulq_f32(Result, sign);
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    // Map V to x in [-pi,pi].
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    XMVECTOR sign = _mm_and_ps(x, g_XMNegativeZero);
    __m128 c = _mm_or_ps(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, g_XMHalfPi);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    x = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, g_XMOne);
    select1 = _mm_andnot_ps(comp, g_XMNegativeOne);
    sign = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation
    const XMVECTOR CEC = g_XMCosCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( CEC, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    vConstants = XM_PERMUTE_PS( CEC, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CEC, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, sign);
    return Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline void XMVectorSinCosEst
(
    XMVECTOR* pSin, 
    XMVECTOR* pCos, 
    FXMVECTOR  V
)
{
    assert(pSin != nullptr);
    assert(pCos != nullptr);

    // 7/6-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Sin;
    XMVECTOR Cos;

    XMScalarSinCosEst(&Sin.vector4_f32[0], &Cos.vector4_f32[0], V.vector4_f32[0]);
    XMScalarSinCosEst(&Sin.vector4_f32[1], &Cos.vector4_f32[1], V.vector4_f32[1]);
    XMScalarSinCosEst(&Sin.vector4_f32[2], &Cos.vector4_f32[2], V.vector4_f32[2]);
    XMScalarSinCosEst(&Sin.vector4_f32[3], &Cos.vector4_f32[3], V.vector4_f32[3]);

    *pSin = Sin;
    *pCos = Cos;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    __n128 sign = vandq_u32(x, g_XMNegativeZero);
    __n128 c = vorrq_u32(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __n128 absx = vabsq_f32( x );
    __n128 rflx = vsubq_f32(c, x);
    __n128 comp = vcleq_f32(absx, g_XMHalfPi);
    x = vbslq_f32( comp, x, rflx );
    sign = vbslq_f32( comp, g_XMOne, g_XMNegativeOne );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation for sine
    const XMVECTOR SEC = g_XMSinCoefficients1;
    XMVECTOR Result = vdupq_lane_f32(vget_high_f32(SEC), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(SEC), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(SEC), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    *pSin = vmulq_f32(Result, x);

    // Compute polynomial approximation
    const XMVECTOR CEC = g_XMCosCoefficients1;
    Result = vdupq_lane_f32(vget_high_f32(CEC), 1);

    vConstants = vdupq_lane_f32(vget_high_f32(CEC), 0);
    Result = vmlaq_f32(vConstants, Result, x2);

    vConstants = vdupq_lane_f32(vget_low_f32(CEC), 1);
    Result = vmlaq_f32(vConstants, Result, x2);

    Result = vmlaq_f32(g_XMOne, Result, x2);
    *pCos = vmulq_f32(Result, sign);
#elif defined(_XM_SSE_INTRINSICS_)
    // Force the value within the bounds of pi
    XMVECTOR x = XMVectorModAngles(V);

    // Map in [-pi/2,pi/2] with sin(y) = sin(x), cos(y) = sign*cos(x).
    XMVECTOR sign = _mm_and_ps(x, g_XMNegativeZero);
    __m128 c = _mm_or_ps(g_XMPi, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, g_XMHalfPi);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    x = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, g_XMOne);
    select1 = _mm_andnot_ps(comp, g_XMNegativeOne);
    sign = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation for sine
    const XMVECTOR SEC = g_XMSinCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( SEC, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    vConstants = XM_PERMUTE_PS( SEC, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( SEC, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, x);
    *pSin = Result;

    // Compute polynomial approximation for cosine
    const XMVECTOR CEC = g_XMCosCoefficients1;
    vConstants = XM_PERMUTE_PS( CEC, _MM_SHUFFLE(3, 3, 3, 3) );
    Result = _mm_mul_ps(vConstants, x2);

    vConstants = XM_PERMUTE_PS( CEC, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( CEC, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    Result = _mm_add_ps(Result, g_XMOne);
    Result = _mm_mul_ps(Result, sign);
    *pCos = Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorTanEst
(
    FXMVECTOR V
)
{
    XMVECTOR OneOverPi = XMVectorSplatW(g_XMTanEstCoefficients.v);

    XMVECTOR V1 = XMVectorMultiply(V, OneOverPi);
    V1 = XMVectorRound(V1);

    V1 = XMVectorNegativeMultiplySubtract(g_XMPi.v, V1, V);

    XMVECTOR T0 = XMVectorSplatX(g_XMTanEstCoefficients.v);
    XMVECTOR T1 = XMVectorSplatY(g_XMTanEstCoefficients.v);
    XMVECTOR T2 = XMVectorSplatZ(g_XMTanEstCoefficients.v);

    XMVECTOR V2T2 = XMVectorNegativeMultiplySubtract(V1, V1, T2);
    XMVECTOR V2 = XMVectorMultiply(V1, V1);
    XMVECTOR V1T0 = XMVectorMultiply(V1, T0);
    XMVECTOR V1T1 = XMVectorMultiply(V1, T1);

    XMVECTOR D = XMVectorReciprocalEst(V2T2);
    XMVECTOR N = XMVectorMultiplyAdd(V2, V1T1, V1T0);

    return XMVectorMultiply(N, D);
}


//------------------------------------------------------------------------------

inline XMVECTOR XMVectorASinEst
(
    FXMVECTOR V
)
{
    // 3-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarASinEst( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarASinEst( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarASinEst( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarASinEst( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 nonnegative = vcgeq_f32(V, g_XMZero);
    __n128 x = vabsq_f32(V);

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __n128 oneMValue = vsubq_f32(g_XMOne, x);
    __n128 clampOneMValue = vmaxq_f32(g_XMZero, oneMValue);
    __n128 root = XMVectorSqrt(clampOneMValue);

    // Compute polynomial approximation
    const XMVECTOR AEC = g_XMArcEstCoefficients;
    __n128 t0 = vdupq_lane_f32(vget_high_f32(AEC), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(AEC), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AEC), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AEC), 0);
    t0 = vmlaq_f32( vConstants, t0, x );
    t0 = vmulq_f32(t0, root);

    __n128 t1 = vsubq_f32(g_XMPi, t0);
    t0 = vbslq_f32( nonnegative, t0, t1 );
    t0 = vsubq_f32(g_XMHalfPi, t0);
    return t0;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128 nonnegative = _mm_cmpge_ps(V, g_XMZero);
    __m128 mvalue = _mm_sub_ps(g_XMZero, V);
    __m128 x = _mm_max_ps(V, mvalue);  // |V|

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __m128 oneMValue = _mm_sub_ps(g_XMOne, x);
    __m128 clampOneMValue = _mm_max_ps(g_XMZero, oneMValue);
    __m128 root = _mm_sqrt_ps(clampOneMValue);  // sqrt(1-|V|)

    // Compute polynomial approximation
    const XMVECTOR AEC = g_XMArcEstCoefficients;
    XMVECTOR vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 t0 = _mm_mul_ps(vConstants, x);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(2, 2, 2, 2) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(1, 1, 1, 1) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(0, 0, 0, 0) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, root);

    __m128 t1 = _mm_sub_ps(g_XMPi, t0);
    t0 = _mm_and_ps(nonnegative, t0);
    t1 = _mm_andnot_ps(nonnegative, t1);
    t0 = _mm_or_ps(t0, t1);
    t0 = _mm_sub_ps(g_XMHalfPi, t0);
    return t0;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorACosEst
(
    FXMVECTOR V
)
{
    // 3-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = XMScalarACosEst( V.vector4_f32[0] );
    Result.vector4_f32[1] = XMScalarACosEst( V.vector4_f32[1] );
    Result.vector4_f32[2] = XMScalarACosEst( V.vector4_f32[2] );
    Result.vector4_f32[3] = XMScalarACosEst( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 nonnegative = vcgeq_f32(V, g_XMZero);
    __n128 x = vabsq_f32(V);

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __n128 oneMValue = vsubq_f32(g_XMOne, x);
    __n128 clampOneMValue = vmaxq_f32(g_XMZero, oneMValue);
    __n128 root = XMVectorSqrt(clampOneMValue);

    // Compute polynomial approximation
    const XMVECTOR AEC = g_XMArcEstCoefficients;
    __n128 t0 = vdupq_lane_f32(vget_high_f32(AEC), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(AEC), 0);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AEC), 1);
    t0 = vmlaq_f32( vConstants, t0, x );

    vConstants = vdupq_lane_f32(vget_low_f32(AEC), 0);
    t0 = vmlaq_f32( vConstants, t0, x );
    t0 = vmulq_f32(t0, root);

    __n128 t1 = vsubq_f32(g_XMPi, t0);
    t0 = vbslq_f32( nonnegative, t0, t1 );
    return t0;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128 nonnegative = _mm_cmpge_ps(V, g_XMZero);
    __m128 mvalue = _mm_sub_ps(g_XMZero, V);
    __m128 x = _mm_max_ps(V, mvalue);  // |V|

    // Compute (1-|V|), clamp to zero to avoid sqrt of negative number.
    __m128 oneMValue = _mm_sub_ps(g_XMOne, x);
    __m128 clampOneMValue = _mm_max_ps(g_XMZero, oneMValue);
    __m128 root = _mm_sqrt_ps(clampOneMValue);  // sqrt(1-|V|)

    // Compute polynomial approximation
    const XMVECTOR AEC = g_XMArcEstCoefficients;
    XMVECTOR vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 t0 = _mm_mul_ps(vConstants, x);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(2, 2, 2, 2) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(1, 1, 1, 1) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, x);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(0, 0, 0, 0) );
    t0 = _mm_add_ps(t0, vConstants);
    t0 = _mm_mul_ps(t0, root);

    __m128 t1 = _mm_sub_ps(g_XMPi, t0);
    t0 = _mm_and_ps(nonnegative, t0);
    t1 = _mm_andnot_ps(nonnegative, t1);
    t0 = _mm_or_ps(t0, t1);
    return t0;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

namespace Internal
{

inline float XMScalarATanEst
(
    float Value
)
{
    float y, sign;
    if (fabsf(Value) <= 1.0f)
    {
        y = Value;
        sign = 0.0f;
    }
    else if (Value > 1.0f)
    {
        y = 1.0f / Value;
        sign = 1.0f;
    }
    else
    {
        y = 1.0f / Value;
        sign = -1.0f;
    }

    // 9-degree minimax approximation
    float y2 = y*y;
    float poly = ((((0.0208351f*y2-0.085133f)*y2+0.180141f)*y2-0.3302995f)*y2+0.999866f)*y;

    return (sign == 0.0f ? poly : sign*XM_PIDIV2 - poly);
}

};  // namespace Internal

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorATanEst
(
    FXMVECTOR V
)
{
    // 9-degree minimax approximation

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result.vector4_f32[0] = Internal::XMScalarATanEst( V.vector4_f32[0] );
    Result.vector4_f32[1] = Internal::XMScalarATanEst( V.vector4_f32[1] );
    Result.vector4_f32[2] = Internal::XMScalarATanEst( V.vector4_f32[2] );
    Result.vector4_f32[3] = Internal::XMScalarATanEst( V.vector4_f32[3] );
    return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 absV = vabsq_f32(V);
    __n128 invV = XMVectorReciprocalEst(V);
    __n128 comp = vcgtq_f32(V, g_XMOne);
    __n128 sign = vbslq_f32(comp, g_XMOne, g_XMNegativeOne );
    comp = vcleq_f32(absV, g_XMOne);
    sign = vbslq_f32(comp, g_XMZero, sign );
    __n128 x = vbslq_f32(comp, V, invV );

    __n128 x2 = vmulq_f32(x, x);

    // Compute polynomial approximation
    const XMVECTOR AEC = g_XMATanEstCoefficients1;
    __n128 Result = vdupq_lane_f32(vget_high_f32(AEC), 1);

    XMVECTOR vConstants = vdupq_lane_f32(vget_high_f32(AEC), 0);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_low_f32(AEC), 1);
    Result = vmlaq_f32( vConstants, Result, x2 );

    vConstants = vdupq_lane_f32(vget_low_f32( AEC), 0);
    Result = vmlaq_f32( vConstants, Result, x2 );

    // ATanEstCoefficients0 is already splatted
    Result = vmlaq_f32( g_XMATanEstCoefficients0, Result, x2 );
    Result = vmulq_f32( Result, x );

    float32x4_t result1 = vmulq_f32(sign, g_XMHalfPi);
    result1 = vsubq_f32(result1, Result);

    comp = vceqq_f32(sign, g_XMZero);
    Result = vbslq_f32( comp, Result, result1 );
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128 absV = XMVectorAbs(V);
    __m128 invV = _mm_div_ps(g_XMOne, V);
    __m128 comp = _mm_cmpgt_ps(V, g_XMOne);
    __m128 select0 = _mm_and_ps(comp, g_XMOne);
    __m128 select1 = _mm_andnot_ps(comp, g_XMNegativeOne);
    __m128 sign = _mm_or_ps(select0, select1);
    comp = _mm_cmple_ps(absV, g_XMOne);
    select0 = _mm_and_ps(comp, g_XMZero);
    select1 = _mm_andnot_ps(comp, sign);
    sign = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, V);
    select1 = _mm_andnot_ps(comp, invV);
    __m128 x = _mm_or_ps(select0, select1);

    __m128 x2 = _mm_mul_ps(x, x);

    // Compute polynomial approximation
    const XMVECTOR AEC = g_XMATanEstCoefficients1;
    XMVECTOR vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(3, 3, 3, 3) );
    __m128 Result = _mm_mul_ps(vConstants, x2);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(2, 2, 2, 2) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(1, 1, 1, 1) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    vConstants = XM_PERMUTE_PS( AEC, _MM_SHUFFLE(0, 0, 0, 0) );
    Result = _mm_add_ps(Result, vConstants);
    Result = _mm_mul_ps(Result, x2);

    // ATanEstCoefficients0 is already splatted
    Result = _mm_add_ps(Result, g_XMATanEstCoefficients0);
    Result = _mm_mul_ps(Result, x);
    __m128 result1 = _mm_mul_ps(sign, g_XMHalfPi);
    result1 = _mm_sub_ps(result1, Result);

    comp = _mm_cmpeq_ps(sign, g_XMZero);
    select0 = _mm_and_ps(comp, Result);
    select1 = _mm_andnot_ps(comp, result1);
    Result = _mm_or_ps(select0, select1);
    return Result;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorATan2Est
(
    FXMVECTOR Y, 
    FXMVECTOR X
)
{
    static const XMVECTORF32 ATan2Constants = {XM_PI, XM_PIDIV2, XM_PIDIV4, 2.3561944905f /* Pi*3/4 */};

    const XMVECTOR Zero = XMVectorZero();
    XMVECTOR ATanResultValid = XMVectorTrueInt();

    XMVECTOR Pi = XMVectorSplatX(ATan2Constants);
    XMVECTOR PiOverTwo = XMVectorSplatY(ATan2Constants);
    XMVECTOR PiOverFour = XMVectorSplatZ(ATan2Constants);
    XMVECTOR ThreePiOverFour = XMVectorSplatW(ATan2Constants);

    XMVECTOR YEqualsZero = XMVectorEqual(Y, Zero);
    XMVECTOR XEqualsZero = XMVectorEqual(X, Zero);
    XMVECTOR XIsPositive = XMVectorAndInt(X, g_XMNegativeZero.v);
    XIsPositive = XMVectorEqualInt(XIsPositive, Zero);
    XMVECTOR YEqualsInfinity = XMVectorIsInfinite(Y);
    XMVECTOR XEqualsInfinity = XMVectorIsInfinite(X);

    XMVECTOR YSign = XMVectorAndInt(Y, g_XMNegativeZero.v);
    Pi = XMVectorOrInt(Pi, YSign);
    PiOverTwo = XMVectorOrInt(PiOverTwo, YSign);
    PiOverFour = XMVectorOrInt(PiOverFour, YSign);
    ThreePiOverFour = XMVectorOrInt(ThreePiOverFour, YSign);

    XMVECTOR R1 = XMVectorSelect(Pi, YSign, XIsPositive);
    XMVECTOR R2 = XMVectorSelect(ATanResultValid, PiOverTwo, XEqualsZero);
    XMVECTOR R3 = XMVectorSelect(R2, R1, YEqualsZero);
    XMVECTOR R4 = XMVectorSelect(ThreePiOverFour, PiOverFour, XIsPositive);
    XMVECTOR R5 = XMVectorSelect(PiOverTwo, R4, XEqualsInfinity);
    XMVECTOR Result = XMVectorSelect(R3, R5, YEqualsInfinity);
    ATanResultValid = XMVectorEqualInt(Result, ATanResultValid);

    XMVECTOR Reciprocal = XMVectorReciprocalEst(X);
    XMVECTOR V = XMVectorMultiply(Y, Reciprocal);
    XMVECTOR R0 = XMVectorATanEst(V);

    R1 = XMVectorSelect( Pi, Zero, XIsPositive );
    R2 = XMVectorAdd(R0, R1);

    Result = XMVectorSelect(Result, R2, ATanResultValid);

    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorLerp
(
    FXMVECTOR V0, 
    FXMVECTOR V1, 
    float    t
)
{
    // V0 + t * (V1 - V0)

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Scale = XMVectorReplicate(t);
    XMVECTOR Length = XMVectorSubtract(V1, V0);
    return XMVectorMultiplyAdd(Length, Scale, V0);

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR L = vsubq_f32( V1, V0 );
    return vmlaq_n_f32( V0, L, t );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR L = _mm_sub_ps( V1, V0 );
    XMVECTOR S = _mm_set_ps1( t );
    XMVECTOR Result = _mm_mul_ps( L, S );
    return _mm_add_ps( Result, V0 );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorLerpV
(
    FXMVECTOR V0, 
    FXMVECTOR V1, 
    FXMVECTOR T
)
{
    // V0 + T * (V1 - V0)

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Length = XMVectorSubtract(V1, V0);
    return XMVectorMultiplyAdd(Length, T, V0);

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR L = vsubq_f32( V1, V0 );
    return vmlaq_f32( V0, L, T );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR Length = _mm_sub_ps( V1, V0 );
    XMVECTOR Result = _mm_mul_ps( Length, T );
    return _mm_add_ps( Result, V0 );
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorHermite
(
    FXMVECTOR Position0, 
    FXMVECTOR Tangent0, 
    FXMVECTOR Position1, 
    GXMVECTOR Tangent1, 
    float    t
)
{
    // Result = (2 * t^3 - 3 * t^2 + 1) * Position0 +
    //          (t^3 - 2 * t^2 + t) * Tangent0 +
    //          (-2 * t^3 + 3 * t^2) * Position1 +
    //          (t^3 - t^2) * Tangent1

#if defined(_XM_NO_INTRINSICS_)

    float t2 = t * t;
    float t3 = t * t2;

    XMVECTOR P0 = XMVectorReplicate(2.0f * t3 - 3.0f * t2 + 1.0f);
    XMVECTOR T0 = XMVectorReplicate(t3 - 2.0f * t2 + t);
    XMVECTOR P1 = XMVectorReplicate(-2.0f * t3 + 3.0f * t2);
    XMVECTOR T1 = XMVectorReplicate(t3 - t2);

    XMVECTOR Result = XMVectorMultiply(P0, Position0);
    Result = XMVectorMultiplyAdd(T0, Tangent0, Result);
    Result = XMVectorMultiplyAdd(P1, Position1, Result);
    Result = XMVectorMultiplyAdd(T1, Tangent1, Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    float t2 = t * t;
    float t3 = t * t2;

    XMVECTOR P0 = vdupq_n_f32(2.0f * t3 - 3.0f * t2 + 1.0f);
    XMVECTOR T0 = vdupq_n_f32(t3 - 2.0f * t2 + t);
    XMVECTOR P1 = vdupq_n_f32(-2.0f * t3 + 3.0f * t2);
    XMVECTOR T1 = vdupq_n_f32(t3 - t2);

    XMVECTOR vResult = vmulq_f32(P0, Position0);
    vResult = vmlaq_f32( vResult, T0, Tangent0 );
    vResult = vmlaq_f32( vResult, P1, Position1 );
    vResult = vmlaq_f32( vResult, T1, Tangent1 );
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    float t2 = t * t;
    float t3 = t * t2;

    XMVECTOR P0 = _mm_set_ps1(2.0f * t3 - 3.0f * t2 + 1.0f);
    XMVECTOR T0 = _mm_set_ps1(t3 - 2.0f * t2 + t);
    XMVECTOR P1 = _mm_set_ps1(-2.0f * t3 + 3.0f * t2);
    XMVECTOR T1 = _mm_set_ps1(t3 - t2);

    XMVECTOR vResult = _mm_mul_ps(P0, Position0);
    XMVECTOR vTemp = _mm_mul_ps(T0, Tangent0);
    vResult = _mm_add_ps(vResult,vTemp);
    vTemp = _mm_mul_ps(P1, Position1);
    vResult = _mm_add_ps(vResult,vTemp);
    vTemp = _mm_mul_ps(T1, Tangent1);
    vResult = _mm_add_ps(vResult,vTemp);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorHermiteV
(
    FXMVECTOR Position0, 
    FXMVECTOR Tangent0, 
    FXMVECTOR Position1, 
    GXMVECTOR Tangent1, 
    CXMVECTOR T
)
{
    // Result = (2 * t^3 - 3 * t^2 + 1) * Position0 +
    //          (t^3 - 2 * t^2 + t) * Tangent0 +
    //          (-2 * t^3 + 3 * t^2) * Position1 +
    //          (t^3 - t^2) * Tangent1

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR T2 = XMVectorMultiply(T, T);
    XMVECTOR T3 = XMVectorMultiply(T , T2);

    XMVECTOR P0 = XMVectorReplicate(2.0f * T3.vector4_f32[0] - 3.0f * T2.vector4_f32[0] + 1.0f);
    XMVECTOR T0 = XMVectorReplicate(T3.vector4_f32[1] - 2.0f * T2.vector4_f32[1] + T.vector4_f32[1]);
    XMVECTOR P1 = XMVectorReplicate(-2.0f * T3.vector4_f32[2] + 3.0f * T2.vector4_f32[2]);
    XMVECTOR T1 = XMVectorReplicate(T3.vector4_f32[3] - T2.vector4_f32[3]);

    XMVECTOR Result = XMVectorMultiply(P0, Position0);
    Result = XMVectorMultiplyAdd(T0, Tangent0, Result);
    Result = XMVectorMultiplyAdd(P1, Position1, Result);
    Result = XMVectorMultiplyAdd(T1, Tangent1, Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 CatMulT2 = {-3.0f,-2.0f,3.0f,-1.0f};
    static const XMVECTORF32 CatMulT3 = {2.0f,1.0f,-2.0f,1.0f};

    XMVECTOR T2 = vmulq_f32(T,T);
    XMVECTOR T3 = vmulq_f32(T,T2);
    // Mul by the constants against t^2
    T2 = vmulq_f32(T2,CatMulT2);
    // Mul by the constants against t^3
    T3 = vmlaq_f32(T2, T3, CatMulT3 );
    // T3 now has the pre-result.
    // I need to add t.y only
    T2 = vandq_u32(T,g_XMMaskY);
    T3 = vaddq_f32(T3,T2);
    // Add 1.0f to x
    T3 = vaddq_f32(T3,g_XMIdentityR0);
    // Now, I have the constants created
    // Mul the x constant to Position0
    XMVECTOR vResult = vdupq_lane_f32( vget_low_f32( T3 ), 0 ); // T3[0]
    vResult = vmulq_f32(vResult,Position0);
    // Mul the y constant to Tangent0
    T2 = vdupq_lane_f32( vget_low_f32( T3 ), 1 ); // T3[1]
    vResult = vmlaq_f32(vResult, T2, Tangent0 );
    // Mul the z constant to Position1
    T2 = vdupq_lane_f32( vget_high_f32( T3 ), 0 ); // T3[2]
    vResult = vmlaq_f32(vResult, T2, Position1 );
    // Mul the w constant to Tangent1
    T3 = vdupq_lane_f32( vget_high_f32( T3 ), 1 ); // T3[3]
    vResult = vmlaq_f32(vResult, T3, Tangent1 );
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 CatMulT2 = {-3.0f,-2.0f,3.0f,-1.0f};
    static const XMVECTORF32 CatMulT3 = {2.0f,1.0f,-2.0f,1.0f};

    XMVECTOR T2 = _mm_mul_ps(T,T);
    XMVECTOR T3 = _mm_mul_ps(T,T2);
    // Mul by the constants against t^2
    T2 = _mm_mul_ps(T2,CatMulT2);
    // Mul by the constants against t^3
    T3 = _mm_mul_ps(T3,CatMulT3);
    // T3 now has the pre-result.
    T3 = _mm_add_ps(T3,T2);
    // I need to add t.y only
    T2 = _mm_and_ps(T,g_XMMaskY);
    T3 = _mm_add_ps(T3,T2);
    // Add 1.0f to x
    T3 = _mm_add_ps(T3,g_XMIdentityR0);
    // Now, I have the constants created
    // Mul the x constant to Position0
    XMVECTOR vResult = XM_PERMUTE_PS(T3,_MM_SHUFFLE(0,0,0,0));
    vResult = _mm_mul_ps(vResult,Position0);
    // Mul the y constant to Tangent0
    T2 = XM_PERMUTE_PS(T3,_MM_SHUFFLE(1,1,1,1));
    T2 = _mm_mul_ps(T2,Tangent0);
    vResult = _mm_add_ps(vResult,T2);
    // Mul the z constant to Position1
    T2 = XM_PERMUTE_PS(T3,_MM_SHUFFLE(2,2,2,2));
    T2 = _mm_mul_ps(T2,Position1);
    vResult = _mm_add_ps(vResult,T2);
    // Mul the w constant to Tangent1
    T3 = XM_PERMUTE_PS(T3,_MM_SHUFFLE(3,3,3,3));
    T3 = _mm_mul_ps(T3,Tangent1);
    vResult = _mm_add_ps(vResult,T3);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorCatmullRom
(
    FXMVECTOR Position0, 
    FXMVECTOR Position1, 
    FXMVECTOR Position2, 
    GXMVECTOR Position3, 
    float    t
)
{
    // Result = ((-t^3 + 2 * t^2 - t) * Position0 +
    //           (3 * t^3 - 5 * t^2 + 2) * Position1 +
    //           (-3 * t^3 + 4 * t^2 + t) * Position2 +
    //           (t^3 - t^2) * Position3) * 0.5

#if defined(_XM_NO_INTRINSICS_)

    float t2 = t * t;
    float t3 = t * t2;

    XMVECTOR P0 = XMVectorReplicate((-t3 + 2.0f * t2 - t) * 0.5f);
    XMVECTOR P1 = XMVectorReplicate((3.0f * t3 - 5.0f * t2 + 2.0f) * 0.5f);
    XMVECTOR P2 = XMVectorReplicate((-3.0f * t3 + 4.0f * t2 + t) * 0.5f);
    XMVECTOR P3 = XMVectorReplicate((t3 - t2) * 0.5f);

    XMVECTOR Result = XMVectorMultiply(P0, Position0);
    Result = XMVectorMultiplyAdd(P1, Position1, Result);
    Result = XMVectorMultiplyAdd(P2, Position2, Result);
    Result = XMVectorMultiplyAdd(P3, Position3, Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    float t2 = t * t;
    float t3 = t * t2;

    XMVECTOR P0 = vdupq_n_f32((-t3 + 2.0f * t2 - t) * 0.5f);
    XMVECTOR P1 = vdupq_n_f32((3.0f * t3 - 5.0f * t2 + 2.0f) * 0.5f);
    XMVECTOR P2 = vdupq_n_f32((-3.0f * t3 + 4.0f * t2 + t) * 0.5f);
    XMVECTOR P3 = vdupq_n_f32((t3 - t2) * 0.5f);

    P1 = vmulq_f32(P1, Position1);
    P0 = vmlaq_f32(P1, P0, Position0);
    P3 = vmulq_f32(P3, Position3);
    P2 = vmlaq_f32(P3, P2, Position2);
    P0 = vaddq_f32(P0,P2);
    return P0;
#elif defined(_XM_SSE_INTRINSICS_)
    float t2 = t * t;
    float t3 = t * t2;

    XMVECTOR P0 = _mm_set_ps1((-t3 + 2.0f * t2 - t) * 0.5f);
    XMVECTOR P1 = _mm_set_ps1((3.0f * t3 - 5.0f * t2 + 2.0f) * 0.5f);
    XMVECTOR P2 = _mm_set_ps1((-3.0f * t3 + 4.0f * t2 + t) * 0.5f);
    XMVECTOR P3 = _mm_set_ps1((t3 - t2) * 0.5f);

    P0 = _mm_mul_ps(P0, Position0);
    P1 = _mm_mul_ps(P1, Position1);
    P2 = _mm_mul_ps(P2, Position2);
    P3 = _mm_mul_ps(P3, Position3);
    P0 = _mm_add_ps(P0,P1);
    P2 = _mm_add_ps(P2,P3);
    P0 = _mm_add_ps(P0,P2);
    return P0;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorCatmullRomV
(
    FXMVECTOR Position0, 
    FXMVECTOR Position1, 
    FXMVECTOR Position2, 
    GXMVECTOR Position3, 
    CXMVECTOR T
)
{
#if defined(_XM_NO_INTRINSICS_)
    float fx = T.vector4_f32[0];
    float fy = T.vector4_f32[1];
    float fz = T.vector4_f32[2];
    float fw = T.vector4_f32[3];
    XMVECTOR vResult = {
        0.5f*((-fx*fx*fx+2*fx*fx-fx)*Position0.vector4_f32[0]+
        (3*fx*fx*fx-5*fx*fx+2)*Position1.vector4_f32[0]+
        (-3*fx*fx*fx+4*fx*fx+fx)*Position2.vector4_f32[0]+
        (fx*fx*fx-fx*fx)*Position3.vector4_f32[0]),
        0.5f*((-fy*fy*fy+2*fy*fy-fy)*Position0.vector4_f32[1]+
        (3*fy*fy*fy-5*fy*fy+2)*Position1.vector4_f32[1]+
        (-3*fy*fy*fy+4*fy*fy+fy)*Position2.vector4_f32[1]+
        (fy*fy*fy-fy*fy)*Position3.vector4_f32[1]),
        0.5f*((-fz*fz*fz+2*fz*fz-fz)*Position0.vector4_f32[2]+
        (3*fz*fz*fz-5*fz*fz+2)*Position1.vector4_f32[2]+
        (-3*fz*fz*fz+4*fz*fz+fz)*Position2.vector4_f32[2]+
        (fz*fz*fz-fz*fz)*Position3.vector4_f32[2]),
        0.5f*((-fw*fw*fw+2*fw*fw-fw)*Position0.vector4_f32[3]+
        (3*fw*fw*fw-5*fw*fw+2)*Position1.vector4_f32[3]+
        (-3*fw*fw*fw+4*fw*fw+fw)*Position2.vector4_f32[3]+
        (fw*fw*fw-fw*fw)*Position3.vector4_f32[3])
    };
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Catmul2 = {2.0f,2.0f,2.0f,2.0f};
    static const XMVECTORF32 Catmul3 = {3.0f,3.0f,3.0f,3.0f};
    static const XMVECTORF32 Catmul4 = {4.0f,4.0f,4.0f,4.0f};
    static const XMVECTORF32 Catmul5 = {5.0f,5.0f,5.0f,5.0f};
    // Cache T^2 and T^3
    XMVECTOR T2 = vmulq_f32(T,T);
    XMVECTOR T3 = vmulq_f32(T,T2);
    // Perform the Position0 term
    XMVECTOR vResult = vaddq_f32(T2,T2);
    vResult = vsubq_f32(vResult,T);
    vResult = vsubq_f32(vResult,T3);
    vResult = vmulq_f32(vResult,Position0);
    // Perform the Position1 term and add
    XMVECTOR vTemp = vmulq_f32(T3,Catmul3);
    vTemp = vmlsq_f32(vTemp, T2, Catmul5);
    vTemp = vaddq_f32(vTemp,Catmul2);
    vResult = vmlaq_f32(vResult, vTemp, Position1);
    // Perform the Position2 term and add
    vTemp = vmulq_f32(T2,Catmul4);
    vTemp = vmlsq_f32(vTemp, T3, Catmul3);
    vTemp = vaddq_f32(vTemp,T);
    vResult = vmlaq_f32(vResult, vTemp, Position2);
    // Position3 is the last term
    T3 = vsubq_f32(T3,T2);
    vResult = vmlaq_f32(vResult, T3, Position3);
    // Multiply by 0.5f and exit
    vResult = vmulq_f32(vResult,g_XMOneHalf);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Catmul2 = {2.0f,2.0f,2.0f,2.0f};
    static const XMVECTORF32 Catmul3 = {3.0f,3.0f,3.0f,3.0f};
    static const XMVECTORF32 Catmul4 = {4.0f,4.0f,4.0f,4.0f};
    static const XMVECTORF32 Catmul5 = {5.0f,5.0f,5.0f,5.0f};
    // Cache T^2 and T^3
    XMVECTOR T2 = _mm_mul_ps(T,T);
    XMVECTOR T3 = _mm_mul_ps(T,T2);
    // Perform the Position0 term
    XMVECTOR vResult = _mm_add_ps(T2,T2);
    vResult = _mm_sub_ps(vResult,T);
    vResult = _mm_sub_ps(vResult,T3);
    vResult = _mm_mul_ps(vResult,Position0);
    // Perform the Position1 term and add
    XMVECTOR vTemp = _mm_mul_ps(T3,Catmul3);
    XMVECTOR vTemp2 = _mm_mul_ps(T2,Catmul5);
    vTemp = _mm_sub_ps(vTemp,vTemp2);
    vTemp = _mm_add_ps(vTemp,Catmul2);
    vTemp = _mm_mul_ps(vTemp,Position1);
    vResult = _mm_add_ps(vResult,vTemp);
    // Perform the Position2 term and add
    vTemp = _mm_mul_ps(T2,Catmul4);
    vTemp2 = _mm_mul_ps(T3,Catmul3);
    vTemp = _mm_sub_ps(vTemp,vTemp2);
    vTemp = _mm_add_ps(vTemp,T);
    vTemp = _mm_mul_ps(vTemp,Position2);
    vResult = _mm_add_ps(vResult,vTemp);
    // Position3 is the last term
    T3 = _mm_sub_ps(T3,T2);
    T3 = _mm_mul_ps(T3,Position3);
    vResult = _mm_add_ps(vResult,T3);
    // Multiply by 0.5f and exit
    vResult = _mm_mul_ps(vResult,g_XMOneHalf);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorBaryCentric
(
    FXMVECTOR Position0, 
    FXMVECTOR Position1, 
    FXMVECTOR Position2, 
    float    f, 
    float    g
)
{
    // Result = Position0 + f * (Position1 - Position0) + g * (Position2 - Position0)

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR P10 = XMVectorSubtract(Position1, Position0);
    XMVECTOR ScaleF = XMVectorReplicate(f);

    XMVECTOR P20 = XMVectorSubtract(Position2, Position0);
    XMVECTOR ScaleG = XMVectorReplicate(g);

    XMVECTOR Result = XMVectorMultiplyAdd(P10, ScaleF, Position0);
    Result = XMVectorMultiplyAdd(P20, ScaleG, Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR R1 = vsubq_f32(Position1,Position0);
    XMVECTOR SF = vdupq_n_f32(f);
    XMVECTOR R2 = vsubq_f32(Position2,Position0);
    XMVECTOR SG = vdupq_n_f32(g);
    R1 = vmlaq_f32( Position0, R1, SF);
    return vmlaq_f32( R1, R2, SG );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR R1 = _mm_sub_ps(Position1,Position0);
    XMVECTOR SF = _mm_set_ps1(f);
    XMVECTOR R2 = _mm_sub_ps(Position2,Position0);
    XMVECTOR SG = _mm_set_ps1(g);
    R1 = _mm_mul_ps(R1,SF);
    R2 = _mm_mul_ps(R2,SG);
    R1 = _mm_add_ps(R1,Position0);
    R1 = _mm_add_ps(R1,R2);
    return R1;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorBaryCentricV
(
    FXMVECTOR Position0, 
    FXMVECTOR Position1, 
    FXMVECTOR Position2, 
    GXMVECTOR F, 
    CXMVECTOR G
)
{
    // Result = Position0 + f * (Position1 - Position0) + g * (Position2 - Position0)

#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR P10 = XMVectorSubtract(Position1, Position0);
    XMVECTOR P20 = XMVectorSubtract(Position2, Position0);

    XMVECTOR Result = XMVectorMultiplyAdd(P10, F, Position0);
    Result = XMVectorMultiplyAdd(P20, G, Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR R1 = vsubq_f32(Position1,Position0);
    XMVECTOR R2 = vsubq_f32(Position2,Position0);
    R1 = vmlaq_f32( Position0, R1, F );
    return vmlaq_f32( R1, R2, G);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR R1 = _mm_sub_ps(Position1,Position0);
    XMVECTOR R2 = _mm_sub_ps(Position2,Position0);
    R1 = _mm_mul_ps(R1,F);
    R2 = _mm_mul_ps(R2,G);
    R1 = _mm_add_ps(R1,Position0);
    R1 = _mm_add_ps(R1,R2);
    return R1;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * 2D Vector
 *
 ****************************************************************************/

//------------------------------------------------------------------------------
// Comparison operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline bool XMVector2Equal
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] == V2.vector4_f32[0]) && (V1.vector4_f32[1] == V2.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vceq_f32( vget_low_f32(V1), vget_low_f32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
// z and w are don't care
    return (((_mm_movemask_ps(vTemp)&3)==3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector2EqualR(V1, V2));
#endif
}


//------------------------------------------------------------------------------

inline uint32_t XMVector2EqualR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    uint32_t CR = 0;
    if ((V1.vector4_f32[0] == V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] == V2.vector4_f32[1]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] != V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] != V2.vector4_f32[1]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vceq_f32( vget_low_f32(V1), vget_low_f32(V2) );
    uint64_t r = vget_lane_u64( vTemp, 0 );
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
// z and w are don't care
    int iTest = _mm_movemask_ps(vTemp)&3;
    uint32_t CR = 0;
    if (iTest==3)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector2EqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_u32[0] == V2.vector4_u32[0]) && (V1.vector4_u32[1] == V2.vector4_u32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vceq_u32( vget_low_u32(V1), vget_low_u32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    return (((_mm_movemask_ps(_mm_castsi128_ps(vTemp))&3)==3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector2EqualIntR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector2EqualIntR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    uint32_t CR = 0;
    if ((V1.vector4_u32[0] == V2.vector4_u32[0]) && 
        (V1.vector4_u32[1] == V2.vector4_u32[1]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_u32[0] != V2.vector4_u32[0]) && 
        (V1.vector4_u32[1] != V2.vector4_u32[1]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vceq_u32( vget_low_u32(V1), vget_low_u32(V2) );
    uint64_t r = vget_lane_u64( vTemp, 0 );
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    int iTest = _mm_movemask_ps(_mm_castsi128_ps(vTemp))&3;
    uint32_t CR = 0;
    if (iTest==3)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector2NearEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR Epsilon
)
{
#if defined(_XM_NO_INTRINSICS_)
    float dx = fabsf(V1.vector4_f32[0]-V2.vector4_f32[0]);
    float dy = fabsf(V1.vector4_f32[1]-V2.vector4_f32[1]);
    return ((dx <= Epsilon.vector4_f32[0]) &&
            (dy <= Epsilon.vector4_f32[1]));
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vDelta = vsub_f32(vget_low_u32(V1), vget_low_u32(V2));
    __n64 vTemp = vacle_f32( vDelta, vget_low_u32(Epsilon) );
    uint64_t r = vget_lane_u64( vTemp, 0 );
    return ( r == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the difference
    XMVECTOR vDelta = _mm_sub_ps(V1,V2);
    // Get the absolute value of the difference
    XMVECTOR vTemp = _mm_setzero_ps();
    vTemp = _mm_sub_ps(vTemp,vDelta);
    vTemp = _mm_max_ps(vTemp,vDelta);
    vTemp = _mm_cmple_ps(vTemp,Epsilon);
    // z and w are don't care
    return (((_mm_movemask_ps(vTemp)&3)==0x3) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector2NotEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] != V2.vector4_f32[0]) || (V1.vector4_f32[1] != V2.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vceq_f32( vget_low_f32(V1), vget_low_f32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) != 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
// z and w are don't care
    return (((_mm_movemask_ps(vTemp)&3)!=3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAnyFalse(XMVector2EqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector2NotEqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_u32[0] != V2.vector4_u32[0]) || (V1.vector4_u32[1] != V2.vector4_u32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vceq_u32( vget_low_u32(V1), vget_low_u32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) != 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    return (((_mm_movemask_ps(_mm_castsi128_ps(vTemp))&3)!=3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAnyFalse(XMVector2EqualIntR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector2Greater
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] > V2.vector4_f32[0]) && (V1.vector4_f32[1] > V2.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vcgt_f32( vget_low_f32(V1), vget_low_f32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
// z and w are don't care
    return (((_mm_movemask_ps(vTemp)&3)==3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector2GreaterR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector2GreaterR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    uint32_t CR = 0;
    if ((V1.vector4_f32[0] > V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] > V2.vector4_f32[1]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] <= V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] <= V2.vector4_f32[1]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vcgt_f32( vget_low_f32(V1), vget_low_f32(V2) );
    uint64_t r = vget_lane_u64( vTemp, 0 );
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
    int iTest = _mm_movemask_ps(vTemp)&3;
    uint32_t CR = 0;
    if (iTest==3)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector2GreaterOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] >= V2.vector4_f32[0]) && (V1.vector4_f32[1] >= V2.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vcge_f32( vget_low_f32(V1), vget_low_f32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&3)==3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector2GreaterOrEqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector2GreaterOrEqualR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    uint32_t CR = 0;
    if ((V1.vector4_f32[0] >= V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] >= V2.vector4_f32[1]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] < V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] < V2.vector4_f32[1]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vcge_f32( vget_low_f32(V1), vget_low_f32(V2) );
    uint64_t r = vget_lane_u64( vTemp, 0 );
    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    int iTest = _mm_movemask_ps(vTemp)&3;
    uint32_t CR = 0;
    if (iTest == 3)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector2Less
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] < V2.vector4_f32[0]) && (V1.vector4_f32[1] < V2.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vclt_f32( vget_low_f32(V1), vget_low_f32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmplt_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&3)==3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector2GreaterR(V2, V1));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector2LessOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vTemp = vcle_f32( vget_low_f32(V1), vget_low_f32(V2) );
    return ( vget_lane_u64( vTemp, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmple_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&3)==3) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector2GreaterOrEqualR(V2, V1));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector2InBounds
(
    FXMVECTOR V, 
    FXMVECTOR Bounds
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V.vector4_f32[0] <= Bounds.vector4_f32[0] && V.vector4_f32[0] >= -Bounds.vector4_f32[0]) && 
        (V.vector4_f32[1] <= Bounds.vector4_f32[1] && V.vector4_f32[1] >= -Bounds.vector4_f32[1])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    __n64 B = vget_low_f32( Bounds );
    // Test if less than or equal
    __n64 vTemp1 = vcle_f32(VL,B);
    // Negate the bounds
    __n64 vTemp2 = vneg_f32(B);
    // Test if greater or equal (Reversed)
    vTemp2 = vcle_f32(vTemp2,VL);
    // Blend answers
    vTemp1 = vand_u32(vTemp1,vTemp2);
    // x and y in bounds?
    return ( vget_lane_u64( vTemp1, 0 ) == 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = _mm_cmple_ps(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = _mm_mul_ps(Bounds,g_XMNegativeOne);
    // Test if greater or equal (Reversed)
    vTemp2 = _mm_cmple_ps(vTemp2,V);
    // Blend answers
    vTemp1 = _mm_and_ps(vTemp1,vTemp2);
    // x and y in bounds? (z and w are don't care)
    return (((_mm_movemask_ps(vTemp1)&0x3)==0x3) != 0);
#else // _XM_VMX128_INTRINSICS_   
    return XMComparisonAllInBounds(XMVector2InBoundsR(V, Bounds));
#endif
}


//------------------------------------------------------------------------------

inline bool XMVector2IsNaN
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (XMISNAN(V.vector4_f32[0]) ||
            XMISNAN(V.vector4_f32[1]));
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    // Test against itself. NaN is always not equal
    __n64 vTempNan = vceq_f32( VL, VL );
    // If x or y are NaN, the mask is zero
    return ( vget_lane_u64( vTempNan, 0 ) != 0xFFFFFFFFFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test against itself. NaN is always not equal
    XMVECTOR vTempNan = _mm_cmpneq_ps(V,V);
    // If x or y are NaN, the mask is non-zero
    return ((_mm_movemask_ps(vTempNan)&3) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector2IsInfinite
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    return (XMISINF(V.vector4_f32[0]) ||
            XMISINF(V.vector4_f32[1]));
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Mask off the sign bit
    __n64 vTemp = vand_u32( vget_low_f32( V ) , vget_low_f32( g_XMAbsMask ) );
    // Compare to infinity
    vTemp = vceq_f32(vTemp, vget_low_f32( g_XMInfinity) );
    // If any are infinity, the signs are true.
    return vget_lane_u64( vTemp, 0 ) != 0;
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bit
    __m128 vTemp = _mm_and_ps(V,g_XMAbsMask);
    // Compare to infinity
    vTemp = _mm_cmpeq_ps(vTemp,g_XMInfinity);
    // If x or z are infinity, the signs are true.
    return ((_mm_movemask_ps(vTemp)&3) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Computation operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Dot
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] =
    Result.vector4_f32[1] =
    Result.vector4_f32[2] =
    Result.vector4_f32[3] = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Perform the dot product on x and y
    __n64 vTemp = vmul_f32( vget_low_f32(V1), vget_low_f32(V2) );
    vTemp = vpadd_f32( vTemp, vTemp );
    return vcombine_f32( vTemp, vTemp );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y
    XMVECTOR vLengthSq = _mm_mul_ps(V1,V2);
    // vTemp has y splatted
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    // x+y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Cross
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    // [ V1.x*V2.y - V1.y*V2.x, V1.x*V2.y - V1.y*V2.x ]

#if defined(_XM_NO_INTRINSICS_)
    float fCross = (V1.vector4_f32[0] * V2.vector4_f32[1]) - (V1.vector4_f32[1] * V2.vector4_f32[0]);
    XMVECTOR vResult = { 
        fCross,
        fCross,
        fCross,
        fCross
    };
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Negate = { 1.f, -1.f, 0, 0 };

    __n64 vTemp = vmul_f32( vget_low_f32( V1 ), vrev64_f32( vget_low_f32( V2 ) ) );
    vTemp = vmul_f32( vTemp, vget_low_f32( Negate ) );
    vTemp = vpadd_f32( vTemp, vTemp );
    return vcombine_f32( vTemp, vTemp );
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap x and y
    XMVECTOR vResult = XM_PERMUTE_PS(V2,_MM_SHUFFLE(0,1,0,1));
    // Perform the muls
    vResult = _mm_mul_ps(vResult,V1);
    // Splat y
    XMVECTOR vTemp = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(1,1,1,1));
    // Sub the values
    vResult = _mm_sub_ss(vResult,vTemp);
    // Splat the cross product
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,0,0,0));
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2LengthSq
(
    FXMVECTOR V
)
{
    return XMVector2Dot(V, V);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2ReciprocalLengthEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector2LengthSq(V);
    Result = XMVectorReciprocalSqrtEst(Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32(V);
    // Dot2
    __n64 vTemp = vmul_f32( VL, VL );
    vTemp = vpadd_f32( vTemp, vTemp );
    // Reciprocal sqrt (estimate)
    vTemp = vrsqrte_f32( vTemp );
    return vcombine_f32( vTemp, vTemp );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has y splatted
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    // x+y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = _mm_rsqrt_ss(vLengthSq);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2ReciprocalLength
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector2LengthSq(V);
    Result = XMVectorReciprocalSqrt(Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32(V);
    // Dot2
    __n64 vTemp = vmul_f32( VL, VL );
    vTemp = vpadd_f32( vTemp, vTemp );
    // Reciprocal sqrt
    __n64  S0 = vrsqrte_f32(vTemp);
    __n64  P0 = vmul_f32( vTemp, S0 );
    __n64  R0 = vrsqrts_f32( P0, S0 );
    __n64  S1 = vmul_f32( S0, R0 );
    __n64  P1 = vmul_f32( vTemp, S1 );
    __n64  R1 = vrsqrts_f32( P1, S1 );
    __n64 Result = vmul_f32( S1, R1 );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has y splatted
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    // x+y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = _mm_sqrt_ss(vLengthSq);
    vLengthSq = _mm_div_ss(g_XMOne,vLengthSq);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2LengthEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector2LengthSq(V);
    Result = XMVectorSqrtEst(Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32(V);
    // Dot2
    __n64 vTemp = vmul_f32( VL, VL );
    vTemp = vpadd_f32( vTemp, vTemp );
    const __n64 zero = vdup_n_u32(0);
    __n64 VEqualsZero = vceq_f32( vTemp, zero );
    // Sqrt (estimate)
    __n64 Result = vrsqrte_f32( vTemp );
    Result = vmul_f32( vTemp, Result );
    Result = vbsl_f32( VEqualsZero, zero, Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has y splatted
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    // x+y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = _mm_sqrt_ss(vLengthSq);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Length
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector2LengthSq(V);
    Result = XMVectorSqrt(Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32(V);
    // Dot2
    __n64 vTemp = vmul_f32( VL, VL );
    vTemp = vpadd_f32( vTemp, vTemp );
    const __n64 zero = vdup_n_u32(0);
    __n64 VEqualsZero = vceq_f32( vTemp, zero );
    // Sqrt
    __n64 S0 = vrsqrte_f32( vTemp );
    __n64 P0 = vmul_f32( vTemp, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( vTemp, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    __n64 Result = vmul_f32( S1, R1 );
    Result = vmul_f32( vTemp, Result );
    Result = vbsl_f32( VEqualsZero, zero, Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has y splatted
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    // x+y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// XMVector2NormalizeEst uses a reciprocal estimate and
// returns QNaN on zero and infinite vectors.

inline XMVECTOR XMVector2NormalizeEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector2ReciprocalLength(V);
    Result = XMVectorMultiply(V, Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32(V);
    // Dot2
    __n64 vTemp = vmul_f32( VL, VL );
    vTemp = vpadd_f32( vTemp, vTemp );
    // Reciprocal sqrt (estimate)
    vTemp = vrsqrte_f32( vTemp );
    // Normalize
    __n64 Result = vmul_f32( VL, vTemp );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has y splatted
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    // x+y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = _mm_rsqrt_ss(vLengthSq);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    vLengthSq = _mm_mul_ps(vLengthSq,V);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Normalize
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR vResult = XMVector2Length( V );
    float fLength = vResult.vector4_f32[0];

    // Prevent divide by zero
    if (fLength > 0) {
        fLength = 1.0f/fLength;
    }
    
    vResult.vector4_f32[0] = V.vector4_f32[0]*fLength;
    vResult.vector4_f32[1] = V.vector4_f32[1]*fLength;
    vResult.vector4_f32[2] = V.vector4_f32[2]*fLength;
    vResult.vector4_f32[3] = V.vector4_f32[3]*fLength;
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32(V);
    // Dot2
    __n64 vTemp = vmul_f32( VL, VL );
    vTemp = vpadd_f32( vTemp, vTemp );
    __n64 VEqualsZero = vceq_f32( vTemp, vdup_n_u32(0) );
    __n64 VEqualsInf = vceq_f32( vTemp, vget_low_f32(g_XMInfinity) );
    // Reciprocal sqrt (2 iterations of Newton-Raphson)
    __n64 S0 = vrsqrte_f32( vTemp );
    __n64 P0 = vmul_f32( vTemp, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( vTemp, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    vTemp = vmul_f32( S1, R1 );
    // Normalize
    __n64 Result = vmul_f32( VL, vTemp );
    Result = vbsl_f32( VEqualsZero, vdup_n_f32(0), Result );
    Result = vbsl_f32( VEqualsInf, vget_low_f32(g_XMQNaN), Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x and y only
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,1,1,1));
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    // Prepare for the division
    XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
    // Create zero with a single instruction
    XMVECTOR vZeroMask = _mm_setzero_ps();
    // Test for a divide by zero (Must be FP to detect -0.0)
    vZeroMask = _mm_cmpneq_ps(vZeroMask,vResult);
    // Failsafe on zero (Or epsilon) length planes
    // If the length is infinity, set the elements to zero
    vLengthSq = _mm_cmpneq_ps(vLengthSq,g_XMInfinity);
    // Reciprocal mul to perform the normalization
    vResult = _mm_div_ps(V,vResult);
    // Any that are infinity, set to zero
    vResult = _mm_and_ps(vResult,vZeroMask);
    // Select qnan or result based on infinite length
    XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq,g_XMQNaN);
    XMVECTOR vTemp2 = _mm_and_ps(vResult,vLengthSq);
    vResult = _mm_or_ps(vTemp1,vTemp2);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2ClampLength
(
    FXMVECTOR V, 
    float    LengthMin, 
    float    LengthMax
)
{
    XMVECTOR ClampMax = XMVectorReplicate(LengthMax);
    XMVECTOR ClampMin = XMVectorReplicate(LengthMin);
    return XMVector2ClampLengthV(V, ClampMin, ClampMax);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2ClampLengthV
(
    FXMVECTOR V, 
    FXMVECTOR LengthMin, 
    FXMVECTOR LengthMax
)
{
    assert((XMVectorGetY(LengthMin) == XMVectorGetX(LengthMin)));
    assert((XMVectorGetY(LengthMax) == XMVectorGetX(LengthMax)));
    assert(XMVector2GreaterOrEqual(LengthMin, g_XMZero));
    assert(XMVector2GreaterOrEqual(LengthMax, g_XMZero));
    assert(XMVector2GreaterOrEqual(LengthMax, LengthMin));

    XMVECTOR LengthSq = XMVector2LengthSq(V);

    const XMVECTOR Zero = XMVectorZero();

    XMVECTOR RcpLength = XMVectorReciprocalSqrt(LengthSq);

    XMVECTOR InfiniteLength = XMVectorEqualInt(LengthSq, g_XMInfinity.v);
    XMVECTOR ZeroLength = XMVectorEqual(LengthSq, Zero);

    XMVECTOR Length = XMVectorMultiply(LengthSq, RcpLength);

    XMVECTOR Normal = XMVectorMultiply(V, RcpLength);

    XMVECTOR Select = XMVectorEqualInt(InfiniteLength, ZeroLength);
    Length = XMVectorSelect(LengthSq, Length, Select);
    Normal = XMVectorSelect(LengthSq, Normal, Select);

    XMVECTOR ControlMax = XMVectorGreater(Length, LengthMax);
    XMVECTOR ControlMin = XMVectorLess(Length, LengthMin);

    XMVECTOR ClampLength = XMVectorSelect(Length, LengthMax, ControlMax);
    ClampLength = XMVectorSelect(ClampLength, LengthMin, ControlMin);

    XMVECTOR Result = XMVectorMultiply(Normal, ClampLength);

    // Preserve the original vector (with no precision loss) if the length falls within the given range
    XMVECTOR Control = XMVectorEqualInt(ControlMax, ControlMin);
    Result = XMVectorSelect(Result, V, Control);

    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Reflect
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal
)
{
    // Result = Incident - (2 * dot(Incident, Normal)) * Normal

    XMVECTOR Result;
    Result = XMVector2Dot(Incident, Normal);
    Result = XMVectorAdd(Result, Result);
    Result = XMVectorNegativeMultiplySubtract(Result, Normal, Incident);
    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Refract
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal, 
    float    RefractionIndex
)
{
    XMVECTOR Index = XMVectorReplicate(RefractionIndex);
    return XMVector2RefractV(Incident, Normal, Index);
}

//------------------------------------------------------------------------------

// Return the refraction of a 2D vector
inline XMVECTOR XMVector2RefractV
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal, 
    FXMVECTOR RefractionIndex
)
{
    // Result = RefractionIndex * Incident - Normal * (RefractionIndex * dot(Incident, Normal) + 
    // sqrt(1 - RefractionIndex * RefractionIndex * (1 - dot(Incident, Normal) * dot(Incident, Normal))))

#if defined(_XM_NO_INTRINSICS_)

    float IDotN = (Incident.vector4_f32[0]*Normal.vector4_f32[0])+(Incident.vector4_f32[1]*Normal.vector4_f32[1]);
    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    float RY = 1.0f-(IDotN*IDotN);
    float RX = 1.0f-(RY*RefractionIndex.vector4_f32[0]*RefractionIndex.vector4_f32[0]);
    RY = 1.0f-(RY*RefractionIndex.vector4_f32[1]*RefractionIndex.vector4_f32[1]);
    if (RX>=0.0f) {
        RX = (RefractionIndex.vector4_f32[0]*Incident.vector4_f32[0])-(Normal.vector4_f32[0]*((RefractionIndex.vector4_f32[0]*IDotN)+sqrtf(RX)));
    } else {
        RX = 0.0f;
    }
    if (RY>=0.0f) {
        RY = (RefractionIndex.vector4_f32[1]*Incident.vector4_f32[1])-(Normal.vector4_f32[1]*((RefractionIndex.vector4_f32[1]*IDotN)+sqrtf(RY)));
    } else {
        RY = 0.0f;
    }

    XMVECTOR vResult;
    vResult.vector4_f32[0] = RX;
    vResult.vector4_f32[1] = RY;
    vResult.vector4_f32[2] = 0.0f;   
    vResult.vector4_f32[3] = 0.0f;
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 IL = vget_low_f32( Incident );
    __n64 NL = vget_low_f32( Normal );
    __n64 RIL = vget_low_f32( RefractionIndex );
    // Get the 2D Dot product of Incident-Normal
    __n64 vTemp = vmul_f32(IL, NL);
    __n64 IDotN = vpadd_f32( vTemp, vTemp );
    // vTemp = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    vTemp = vmls_f32( vget_low_f32( g_XMOne ), IDotN, IDotN);
    vTemp = vmul_f32(vTemp,RIL);
    vTemp = vmls_f32(vget_low_f32( g_XMOne ), vTemp, RIL );
    // If any terms are <=0, sqrt() will fail, punt to zero
    __n64 vMask = vcgt_f32(vTemp, vget_low_f32(g_XMZero) );
    // Sqrt(vTemp)
    __n64 S0 = vrsqrte_f32(vTemp);
    __n64 P0 = vmul_f32( vTemp, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( vTemp, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    __n64 S2 = vmul_f32( S1, R1 );
    vTemp = vmul_f32( vTemp, S2 );
    // R = RefractionIndex * IDotN + sqrt(R)
    vTemp = vmla_f32( vTemp, RIL, IDotN );
    // Result = RefractionIndex * Incident - Normal * R
    __n64 vResult = vmul_f32(RIL,IL);
    vResult = vmls_f32( vResult, vTemp, NL );
    vResult = vand_u32(vResult,vMask);
    return vcombine_f32(vResult, vResult);
#elif defined(_XM_SSE_INTRINSICS_)
    // Result = RefractionIndex * Incident - Normal * (RefractionIndex * dot(Incident, Normal) + 
    // sqrt(1 - RefractionIndex * RefractionIndex * (1 - dot(Incident, Normal) * dot(Incident, Normal))))
    // Get the 2D Dot product of Incident-Normal
    XMVECTOR IDotN = XMVector2Dot(Incident, Normal);
    // vTemp = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    XMVECTOR vTemp = _mm_mul_ps(IDotN,IDotN);
    vTemp = _mm_sub_ps(g_XMOne,vTemp);
    vTemp = _mm_mul_ps(vTemp,RefractionIndex);
    vTemp = _mm_mul_ps(vTemp,RefractionIndex);
    vTemp = _mm_sub_ps(g_XMOne,vTemp);
    // If any terms are <=0, sqrt() will fail, punt to zero
    XMVECTOR vMask = _mm_cmpgt_ps(vTemp,g_XMZero);
    // R = RefractionIndex * IDotN + sqrt(R)
    vTemp = _mm_sqrt_ps(vTemp);
    XMVECTOR vResult = _mm_mul_ps(RefractionIndex,IDotN);
    vTemp = _mm_add_ps(vTemp,vResult);
    // Result = RefractionIndex * Incident - Normal * R
    vResult = _mm_mul_ps(RefractionIndex,Incident);
    vTemp = _mm_mul_ps(vTemp,Normal);
    vResult = _mm_sub_ps(vResult,vTemp);
    vResult = _mm_and_ps(vResult,vMask);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Orthogonal
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = -V.vector4_f32[1];
    Result.vector4_f32[1] = V.vector4_f32[0];
    Result.vector4_f32[2] = 0.f;
    Result.vector4_f32[3] = 0.f;
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Negate = { -1.f, 1.f, 0, 0 };
    const __n64 zero = vdup_n_f32(0);

    __n64 VL = vget_low_f32( V );
    __n64 Result = vmul_f32( vrev64_f32( VL ), vget_low_f32( Negate ) );
    return vcombine_f32( Result, zero );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,2,0,1));
    vResult = _mm_mul_ps(vResult,g_XMNegateX);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2AngleBetweenNormalsEst
(
    FXMVECTOR N1, 
    FXMVECTOR N2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_) 

    XMVECTOR Result = XMVector2Dot(N1, N2);
    Result = XMVectorClamp(Result, g_XMNegativeOne.v, g_XMOne.v);
    Result = XMVectorACosEst(Result);
    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2AngleBetweenNormals
(
    FXMVECTOR N1, 
    FXMVECTOR N2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Result = XMVector2Dot(N1, N2);
    Result = XMVectorClamp(Result, g_XMNegativeOne, g_XMOne);
    Result = XMVectorACos(Result);
    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2AngleBetweenVectors
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR L1 = XMVector2ReciprocalLength(V1);
    XMVECTOR L2 = XMVector2ReciprocalLength(V2);

    XMVECTOR Dot = XMVector2Dot(V1, V2);

    L1 = XMVectorMultiply(L1, L2);

    XMVECTOR CosAngle = XMVectorMultiply(Dot, L1);
    CosAngle = XMVectorClamp(CosAngle, g_XMNegativeOne.v, g_XMOne.v);

    return XMVectorACos(CosAngle);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2LinePointDistance
(
    FXMVECTOR LinePoint1, 
    FXMVECTOR LinePoint2, 
    FXMVECTOR Point
)
{
    // Given a vector PointVector from LinePoint1 to Point and a vector
    // LineVector from LinePoint1 to LinePoint2, the scaled distance 
    // PointProjectionScale from LinePoint1 to the perpendicular projection
    // of PointVector onto the line is defined as:
    //
    //     PointProjectionScale = dot(PointVector, LineVector) / LengthSq(LineVector)

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR PointVector = XMVectorSubtract(Point, LinePoint1);
    XMVECTOR LineVector = XMVectorSubtract(LinePoint2, LinePoint1);

    XMVECTOR LengthSq = XMVector2LengthSq(LineVector);

    XMVECTOR PointProjectionScale = XMVector2Dot(PointVector, LineVector);
    PointProjectionScale = XMVectorDivide(PointProjectionScale, LengthSq);

    XMVECTOR DistanceVector = XMVectorMultiply(LineVector, PointProjectionScale);
    DistanceVector = XMVectorSubtract(PointVector, DistanceVector);

    return XMVector2Length(DistanceVector);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2IntersectLine
(
    FXMVECTOR Line1Point1, 
    FXMVECTOR Line1Point2, 
    FXMVECTOR Line2Point1, 
    GXMVECTOR Line2Point2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR V1 = XMVectorSubtract(Line1Point2, Line1Point1);
    XMVECTOR V2 = XMVectorSubtract(Line2Point2, Line2Point1);
    XMVECTOR V3 = XMVectorSubtract(Line1Point1, Line2Point1);

    XMVECTOR C1 = XMVector2Cross(V1, V2);
    XMVECTOR C2 = XMVector2Cross(V2, V3);

    XMVECTOR Result;
    const XMVECTOR Zero = XMVectorZero();
    if (XMVector2NearEqual(C1, Zero, g_XMEpsilon.v))
    {
        if (XMVector2NearEqual(C2, Zero, g_XMEpsilon.v))
        {
            // Coincident
            Result = g_XMInfinity.v;
        }
        else
        {
            // Parallel
            Result = g_XMQNaN.v;
        }
    }
    else
    {
        // Intersection point = Line1Point1 + V1 * (C2 / C1)
        XMVECTOR Scale = XMVectorReciprocal(C1);
        Scale = XMVectorMultiply(C2, Scale);
        Result = XMVectorMultiplyAdd(V1, Scale, Line1Point1);
    }

    return Result;

#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR V1 = _mm_sub_ps(Line1Point2, Line1Point1);
    XMVECTOR V2 = _mm_sub_ps(Line2Point2, Line2Point1);
    XMVECTOR V3 = _mm_sub_ps(Line1Point1, Line2Point1);
    // Generate the cross products
    XMVECTOR C1 = XMVector2Cross(V1, V2);
    XMVECTOR C2 = XMVector2Cross(V2, V3);
    // If C1 is not close to epsilon, use the calculated value
    XMVECTOR vResultMask = _mm_setzero_ps();
    vResultMask = _mm_sub_ps(vResultMask,C1);
    vResultMask = _mm_max_ps(vResultMask,C1);
    // 0xFFFFFFFF if the calculated value is to be used
    vResultMask = _mm_cmpgt_ps(vResultMask,g_XMEpsilon);
    // If C1 is close to epsilon, which fail type is it? INFINITY or NAN?
    XMVECTOR vFailMask = _mm_setzero_ps();
    vFailMask = _mm_sub_ps(vFailMask,C2);
    vFailMask = _mm_max_ps(vFailMask,C2);
    vFailMask = _mm_cmple_ps(vFailMask,g_XMEpsilon);
    XMVECTOR vFail = _mm_and_ps(vFailMask,g_XMInfinity);
    vFailMask = _mm_andnot_ps(vFailMask,g_XMQNaN);
    // vFail is NAN or INF
    vFail = _mm_or_ps(vFail,vFailMask);
    // Intersection point = Line1Point1 + V1 * (C2 / C1)
    XMVECTOR vResult = _mm_div_ps(C2,C1);
    vResult = _mm_mul_ps(vResult,V1);
    vResult = _mm_add_ps(vResult,Line1Point1);
    // Use result, or failure value
    vResult = _mm_and_ps(vResult,vResultMask);
    vResultMask = _mm_andnot_ps(vResultMask,vFail);
    vResult = _mm_or_ps(vResult,vResultMask);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2Transform
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiplyAdd(Y, M.r[1], M.r[3]);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    __n128 Y = vdupq_lane_f32( VL, 1 );
    __n128 Result = vmlaq_f32( M.r[3], Y, M.r[1] );
    __n128 X = vdupq_lane_f32( VL, 0 );
    return vmlaq_f32( Result, X, M.r[0] );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,0,0,0));
    vResult = _mm_mul_ps(vResult,M.r[0]);
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    vTemp = _mm_mul_ps(vTemp,M.r[1]);
    vResult = _mm_add_ps(vResult,vTemp);
    vResult = _mm_add_ps(vResult,M.r[3]);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT4* XMVector2TransformStream
(
    XMFLOAT4*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT2* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t* pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];
    const XMVECTOR row3 = M.r[3];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat2((const XMFLOAT2*)pInputVector);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiplyAdd(Y, row1, row3);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMStoreFloat4((XMFLOAT4*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline XMVECTOR XMVector2TransformCoord
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiplyAdd(Y, M.r[1], M.r[3]);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    XMVECTOR W = XMVectorSplatW(Result);
    return XMVectorDivide( Result, W );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT2* XMVector2TransformCoordStream
(
    XMFLOAT2*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT2* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t*    pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];
    const XMVECTOR row3 = M.r[3];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat2((const XMFLOAT2*)pInputVector);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiplyAdd(Y, row1, row3);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMVECTOR W = XMVectorSplatW(Result);

        Result = XMVectorDivide(Result, W);

        XMStoreFloat2((XMFLOAT2*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector2TransformNormal
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiply(Y, M.r[1]);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    __n128 Y = vdupq_lane_f32( VL, 1 );
    __n128 Result = vmulq_f32( Y, M.r[1] );
    __n128 X = vdupq_lane_f32( VL, 0 );
    return vmlaq_f32( Result, X, M.r[0] );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,0,0,0));
    vResult = _mm_mul_ps(vResult,M.r[0]);
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    vTemp = _mm_mul_ps(vTemp,M.r[1]);
    vResult = _mm_add_ps(vResult,vTemp);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT2* XMVector2TransformNormalStream
(
    XMFLOAT2*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT2* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t*    pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat2((const XMFLOAT2*)pInputVector);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiply(Y, row1);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMStoreFloat2((XMFLOAT2*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * 3D Vector
 *
 ****************************************************************************/

//------------------------------------------------------------------------------
// Comparison operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline bool XMVector3Equal
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] == V2.vector4_f32[0]) && (V1.vector4_f32[1] == V2.vector4_f32[1]) && (V1.vector4_f32[2] == V2.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&7)==7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector3EqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector3EqualR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t CR = 0;
    if ((V1.vector4_f32[0] == V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] == V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] == V2.vector4_f32[2]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] != V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] != V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] != V2.vector4_f32[2]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU;

    uint32_t CR = 0;
    if ( r == 0xFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
    int iTest = _mm_movemask_ps(vTemp)&7;
    uint32_t CR = 0;
    if (iTest==7)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector3EqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_u32[0] == V2.vector4_u32[0]) && (V1.vector4_u32[1] == V2.vector4_u32[1]) && (V1.vector4_u32[2] == V2.vector4_u32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    return (((_mm_movemask_ps(_mm_castsi128_ps(vTemp))&7)==7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector3EqualIntR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector3EqualIntR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t CR = 0;
    if ((V1.vector4_u32[0] == V2.vector4_u32[0]) && 
        (V1.vector4_u32[1] == V2.vector4_u32[1]) &&
        (V1.vector4_u32[2] == V2.vector4_u32[2]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_u32[0] != V2.vector4_u32[0]) && 
        (V1.vector4_u32[1] != V2.vector4_u32[1]) &&
        (V1.vector4_u32[2] != V2.vector4_u32[2]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU;

    uint32_t CR = 0;
    if ( r == 0xFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    int iTemp = _mm_movemask_ps(_mm_castsi128_ps(vTemp))&7;
    uint32_t CR = 0;
    if (iTemp==7)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTemp)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector3NearEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR Epsilon
)
{
#if defined(_XM_NO_INTRINSICS_)
    float dx, dy, dz;

    dx = fabsf(V1.vector4_f32[0]-V2.vector4_f32[0]);
    dy = fabsf(V1.vector4_f32[1]-V2.vector4_f32[1]);
    dz = fabsf(V1.vector4_f32[2]-V2.vector4_f32[2]);
    return (((dx <= Epsilon.vector4_f32[0]) &&
            (dy <= Epsilon.vector4_f32[1]) &&
            (dz <= Epsilon.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vDelta = vsubq_f32( V1, V2 );
    __n128 vResult = vacleq_f32( vDelta, Epsilon );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the difference
    XMVECTOR vDelta = _mm_sub_ps(V1,V2);
    // Get the absolute value of the difference
    XMVECTOR vTemp = _mm_setzero_ps();
    vTemp = _mm_sub_ps(vTemp,vDelta);
    vTemp = _mm_max_ps(vTemp,vDelta);
    vTemp = _mm_cmple_ps(vTemp,Epsilon);
    // w is don't care
    return (((_mm_movemask_ps(vTemp)&7)==0x7) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector3NotEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] != V2.vector4_f32[0]) || (V1.vector4_f32[1] != V2.vector4_f32[1]) || (V1.vector4_f32[2] != V2.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1)  & 0xFFFFFFU) != 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&7)!=7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAnyFalse(XMVector3EqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector3NotEqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_u32[0] != V2.vector4_u32[0]) || (V1.vector4_u32[1] != V2.vector4_u32[1]) || (V1.vector4_u32[2] != V2.vector4_u32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1)  & 0xFFFFFFU) != 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    return (((_mm_movemask_ps(_mm_castsi128_ps(vTemp))&7)!=7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAnyFalse(XMVector3EqualIntR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector3Greater
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] > V2.vector4_f32[0]) && (V1.vector4_f32[1] > V2.vector4_f32[1]) && (V1.vector4_f32[2] > V2.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgtq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&7)==7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector3GreaterR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector3GreaterR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t CR = 0;
    if ((V1.vector4_f32[0] > V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] > V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] > V2.vector4_f32[2]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] <= V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] <= V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] <= V2.vector4_f32[2]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgtq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU;

    uint32_t CR = 0;
    if ( r == 0xFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
    uint32_t CR = 0;
    int iTest = _mm_movemask_ps(vTemp)&7;
    if (iTest==7) 
    {
        CR =  XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector3GreaterOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] >= V2.vector4_f32[0]) && (V1.vector4_f32[1] >= V2.vector4_f32[1]) && (V1.vector4_f32[2] >= V2.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgeq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&7)==7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector3GreaterOrEqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector3GreaterOrEqualR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    uint32_t CR = 0;
    if ((V1.vector4_f32[0] >= V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] >= V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] >= V2.vector4_f32[2]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] < V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] < V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] < V2.vector4_f32[2]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgeq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU;

    uint32_t CR = 0;
    if ( r == 0xFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    uint32_t CR = 0;
    int iTest = _mm_movemask_ps(vTemp)&7;
    if (iTest==7) 
    {
        CR =  XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector3Less
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] < V2.vector4_f32[0]) && (V1.vector4_f32[1] < V2.vector4_f32[1]) && (V1.vector4_f32[2] < V2.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcltq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmplt_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&7)==7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector3GreaterR(V2, V1));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector3LessOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1]) && (V1.vector4_f32[2] <= V2.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcleq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmple_ps(V1,V2);
    return (((_mm_movemask_ps(vTemp)&7)==7) != 0);
#else // _XM_VMX128_INTRINSICS_
    return XMComparisonAllTrue(XMVector3GreaterOrEqualR(V2, V1));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector3InBounds
(
    FXMVECTOR V, 
    FXMVECTOR Bounds
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V.vector4_f32[0] <= Bounds.vector4_f32[0] && V.vector4_f32[0] >= -Bounds.vector4_f32[0]) && 
        (V.vector4_f32[1] <= Bounds.vector4_f32[1] && V.vector4_f32[1] >= -Bounds.vector4_f32[1]) &&
        (V.vector4_f32[2] <= Bounds.vector4_f32[2] && V.vector4_f32[2] >= -Bounds.vector4_f32[2])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test if less than or equal
    __n128 vTemp1 = vcleq_f32(V,Bounds);
    // Negate the bounds
    __n128 vTemp2 = vnegq_f32(Bounds);
    // Test if greater or equal (Reversed)
    vTemp2 = vcleq_f32(vTemp2,V);
    // Blend answers
    vTemp1 = vandq_u32(vTemp1,vTemp2);
    // in bounds?
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTemp1), vget_high_u8(vTemp1));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) == 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = _mm_cmple_ps(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = _mm_mul_ps(Bounds,g_XMNegativeOne);
    // Test if greater or equal (Reversed)
    vTemp2 = _mm_cmple_ps(vTemp2,V);
    // Blend answers
    vTemp1 = _mm_and_ps(vTemp1,vTemp2);
    // x,y and z in bounds? (w is don't care)
    return (((_mm_movemask_ps(vTemp1)&0x7)==0x7) != 0);
#else
    return XMComparisonAllInBounds(XMVector3InBoundsR(V, Bounds));
#endif
}


//------------------------------------------------------------------------------

inline bool XMVector3IsNaN
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    return (XMISNAN(V.vector4_f32[0]) ||
            XMISNAN(V.vector4_f32[1]) ||
            XMISNAN(V.vector4_f32[2]));

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test against itself. NaN is always not equal
    __n128 vTempNan = vceqq_f32( V, V );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTempNan), vget_high_u8(vTempNan));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    // If x or y or z are NaN, the mask is zero
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) != 0xFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test against itself. NaN is always not equal
    XMVECTOR vTempNan = _mm_cmpneq_ps(V,V);
    // If x or y or z are NaN, the mask is non-zero
    return ((_mm_movemask_ps(vTempNan)&7) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector3IsInfinite
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (XMISINF(V.vector4_f32[0]) ||
            XMISINF(V.vector4_f32[1]) ||
            XMISINF(V.vector4_f32[2]));
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Mask off the sign bit
    __n128 vTempInf = vandq_u32( V, g_XMAbsMask );
    // Compare to infinity
    vTempInf = vceqq_f32(vTempInf, g_XMInfinity );
    // If any are infinity, the signs are true.
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTempInf), vget_high_u8(vTempInf));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( (vget_lane_u32(vTemp.val[1], 1) & 0xFFFFFFU) != 0 );
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bit
    __m128 vTemp = _mm_and_ps(V,g_XMAbsMask);
    // Compare to infinity
    vTemp = _mm_cmpeq_ps(vTemp,g_XMInfinity);
    // If x,y or z are infinity, the signs are true.
    return ((_mm_movemask_ps(vTemp)&7) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Computation operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Dot
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    float fValue = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1] + V1.vector4_f32[2] * V2.vector4_f32[2];
    XMVECTOR vResult = {
        fValue,
        fValue,
        fValue,
        fValue
    };            
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vTemp = vmulq_f32( V1, V2 );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    return vcombine_f32( v1, v1 );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product
    XMVECTOR vDot = _mm_mul_ps(V1,V2);
    // x=Dot.vector4_f32[1], y=Dot.vector4_f32[2]
    XMVECTOR vTemp = XM_PERMUTE_PS(vDot,_MM_SHUFFLE(2,1,2,1));
    // Result.vector4_f32[0] = x+y
    vDot = _mm_add_ss(vDot,vTemp);
    // x=Dot.vector4_f32[2]
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // Result.vector4_f32[0] = (x+y)+z
    vDot = _mm_add_ss(vDot,vTemp);
    // Splat x
    return XM_PERMUTE_PS(vDot,_MM_SHUFFLE(0,0,0,0));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Cross
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
    // [ V1.y*V2.z - V1.z*V2.y, V1.z*V2.x - V1.x*V2.z, V1.x*V2.y - V1.y*V2.x ]

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR vResult = {
        (V1.vector4_f32[1] * V2.vector4_f32[2]) - (V1.vector4_f32[2] * V2.vector4_f32[1]),
        (V1.vector4_f32[2] * V2.vector4_f32[0]) - (V1.vector4_f32[0] * V2.vector4_f32[2]),
        (V1.vector4_f32[0] * V2.vector4_f32[1]) - (V1.vector4_f32[1] * V2.vector4_f32[0]),
        0.0f
    };
    return vResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 v1xy = vget_low_f32(V1);
    __n64 v2xy = vget_low_f32(V2);

    __n64 v1yx = vrev64_f32( v1xy );
    __n64 v2yx = vrev64_f32( v2xy );

    __n64 v1zz = vdup_lane_f32( vget_high_f32(V1), 0 );
    __n64 v2zz = vdup_lane_f32( vget_high_f32(V2), 0 );

    __n128 vResult = vmulq_f32( vcombine_f32(v1yx,v1xy), vcombine_f32(v2zz,v2yx) );
    vResult = vmlsq_f32( vResult, vcombine_f32(v1zz,v1yx), vcombine_f32(v2yx,v2xy) );
    return veorq_u32( vResult, g_XMFlipY );
#elif defined(_XM_SSE_INTRINSICS_)
    // y1,z1,x1,w1
    XMVECTOR vTemp1 = XM_PERMUTE_PS(V1,_MM_SHUFFLE(3,0,2,1));
    // z2,x2,y2,w2
    XMVECTOR vTemp2 = XM_PERMUTE_PS(V2,_MM_SHUFFLE(3,1,0,2));
    // Perform the left operation
    XMVECTOR vResult = _mm_mul_ps(vTemp1,vTemp2);
    // z1,x1,y1,w1
    vTemp1 = XM_PERMUTE_PS(vTemp1,_MM_SHUFFLE(3,0,2,1));
    // y2,z2,x2,w2
    vTemp2 = XM_PERMUTE_PS(vTemp2,_MM_SHUFFLE(3,1,0,2));
    // Perform the right operation
    vTemp1 = _mm_mul_ps(vTemp1,vTemp2);
    // Subract the right from left, and return answer
    vResult = _mm_sub_ps(vResult,vTemp1);
    // Set w to zero
    return _mm_and_ps(vResult,g_XMMask3);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3LengthSq
(
    FXMVECTOR V
)
{
    return XMVector3Dot(V, V);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3ReciprocalLengthEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector3LengthSq(V);
    Result = XMVectorReciprocalSqrtEst(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot3
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    // Reciprocal sqrt (estimate)
    v2 = vrsqrte_f32( v1 );
    return vcombine_f32(v2, v2);
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and y
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,2,1,2));
    // x+z, y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    // y,y,y,y
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // x+z+y,??,??,??
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    // Splat the length squared
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    // Get the reciprocal
    vLengthSq = _mm_rsqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3ReciprocalLength
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector3LengthSq(V);
    Result = XMVectorReciprocalSqrt(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot3
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    // Reciprocal sqrt
    __n64  S0 = vrsqrte_f32(v1);
    __n64  P0 = vmul_f32( v1, S0 );
    __n64  R0 = vrsqrts_f32( P0, S0 );
    __n64  S1 = vmul_f32( S0, R0 );
    __n64  P1 = vmul_f32( v1, S1 );
    __n64  R1 = vrsqrts_f32( P1, S1 );
    __n64 Result = vmul_f32( S1, R1 );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
     // Perform the dot product
    XMVECTOR vDot = _mm_mul_ps(V,V);
    // x=Dot.y, y=Dot.z
    XMVECTOR vTemp = XM_PERMUTE_PS(vDot,_MM_SHUFFLE(2,1,2,1));
    // Result.x = x+y
    vDot = _mm_add_ss(vDot,vTemp);
    // x=Dot.z
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // Result.x = (x+y)+z
    vDot = _mm_add_ss(vDot,vTemp);
    // Splat x
    vDot = XM_PERMUTE_PS(vDot,_MM_SHUFFLE(0,0,0,0));
    // Get the reciprocal
    vDot = _mm_sqrt_ps(vDot);
    // Get the reciprocal
    vDot = _mm_div_ps(g_XMOne,vDot);
    return vDot;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3LengthEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector3LengthSq(V);
    Result = XMVectorSqrtEst(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot3
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    const __n64 zero = vdup_n_u32(0);
    __n64 VEqualsZero = vceq_f32( v1, zero );
    // Sqrt (estimate)
    __n64 Result = vrsqrte_f32( v1 );
    Result = vmul_f32( v1, Result );
    Result = vbsl_f32( VEqualsZero, zero, Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and y
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,2,1,2));
    // x+z, y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    // y,y,y,y
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // x+z+y,??,??,??
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    // Splat the length squared
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    // Get the length
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Length
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector3LengthSq(V);
    Result = XMVectorSqrt(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot3
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    const __n64 zero = vdup_n_u32(0);
    __n64 VEqualsZero = vceq_f32( v1, zero );
    // Sqrt
    __n64 S0 = vrsqrte_f32( v1 );
    __n64 P0 = vmul_f32( v1, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( v1, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    __n64 Result = vmul_f32( S1, R1 );
    Result = vmul_f32( v1, Result );
    Result = vbsl_f32( VEqualsZero, zero, Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and y
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,2,1,2));
    // x+z, y
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    // y,y,y,y
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // x+z+y,??,??,??
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    // Splat the length squared
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    // Get the length
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// XMVector3NormalizeEst uses a reciprocal estimate and
// returns QNaN on zero and infinite vectors.

inline XMVECTOR XMVector3NormalizeEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector3ReciprocalLength(V);
    Result = XMVectorMultiply(V, Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot3
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    // Reciprocal sqrt (estimate)
    v2 = vrsqrte_f32( v1 );
    // Normalize
    return vmulq_f32( V, vcombine_f32(v2,v2) );
#elif defined(_XM_SSE_INTRINSICS_)
     // Perform the dot product
    XMVECTOR vDot = _mm_mul_ps(V,V);
    // x=Dot.y, y=Dot.z
    XMVECTOR vTemp = XM_PERMUTE_PS(vDot,_MM_SHUFFLE(2,1,2,1));
    // Result.x = x+y
    vDot = _mm_add_ss(vDot,vTemp);
    // x=Dot.z
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // Result.x = (x+y)+z
    vDot = _mm_add_ss(vDot,vTemp);
    // Splat x
    vDot = XM_PERMUTE_PS(vDot,_MM_SHUFFLE(0,0,0,0));
    // Get the reciprocal
    vDot = _mm_rsqrt_ps(vDot);
    // Perform the normalization
    vDot = _mm_mul_ps(vDot,V);
    return vDot;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Normalize
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    float fLength;
    XMVECTOR vResult;

    vResult = XMVector3Length( V );
    fLength = vResult.vector4_f32[0];

    // Prevent divide by zero
    if (fLength > 0) {
        fLength = 1.0f/fLength;
    }
    
    vResult.vector4_f32[0] = V.vector4_f32[0]*fLength;
    vResult.vector4_f32[1] = V.vector4_f32[1]*fLength;
    vResult.vector4_f32[2] = V.vector4_f32[2]*fLength;
    vResult.vector4_f32[3] = V.vector4_f32[3]*fLength;
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot3
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vdup_lane_f32( v2, 0 );
    v1 = vadd_f32( v1, v2 );
    __n64 VEqualsZero = vceq_f32( v1, vdup_n_u32(0) );
    __n64 VEqualsInf = vceq_f32( v1, vget_low_f32(g_XMInfinity) );
    // Reciprocal sqrt (2 iterations of Newton-Raphson)
    __n64 S0 = vrsqrte_f32( v1 );
    __n64 P0 = vmul_f32( v1, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( v1, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    v2 = vmul_f32( S1, R1 );
    // Normalize
    __n128 vResult = vmulq_f32( V, vcombine_f32(v2,v2) );
    vResult = vbslq_f32( vcombine_f32(VEqualsZero,VEqualsZero), vdupq_n_f32(0), vResult );
    return vbslq_f32( vcombine_f32(VEqualsInf,VEqualsInf), g_XMQNaN, vResult );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z only
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,1,2,1));
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    vLengthSq = _mm_add_ss(vLengthSq,vTemp);
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(0,0,0,0));
    // Prepare for the division
    XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
    // Create zero with a single instruction
    XMVECTOR vZeroMask = _mm_setzero_ps();
    // Test for a divide by zero (Must be FP to detect -0.0)
    vZeroMask = _mm_cmpneq_ps(vZeroMask,vResult);
    // Failsafe on zero (Or epsilon) length planes
    // If the length is infinity, set the elements to zero
    vLengthSq = _mm_cmpneq_ps(vLengthSq,g_XMInfinity);
    // Divide to perform the normalization
    vResult = _mm_div_ps(V,vResult);
    // Any that are infinity, set to zero
    vResult = _mm_and_ps(vResult,vZeroMask);
    // Select qnan or result based on infinite length
    XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq,g_XMQNaN);
    XMVECTOR vTemp2 = _mm_and_ps(vResult,vLengthSq);
    vResult = _mm_or_ps(vTemp1,vTemp2);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3ClampLength
(
    FXMVECTOR V, 
    float    LengthMin, 
    float    LengthMax
)
{
    XMVECTOR ClampMax = XMVectorReplicate(LengthMax);
    XMVECTOR ClampMin = XMVectorReplicate(LengthMin);

    return XMVector3ClampLengthV(V, ClampMin, ClampMax);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3ClampLengthV
(
    FXMVECTOR V, 
    FXMVECTOR LengthMin, 
    FXMVECTOR LengthMax
)
{
    assert((XMVectorGetY(LengthMin) == XMVectorGetX(LengthMin)) && (XMVectorGetZ(LengthMin) == XMVectorGetX(LengthMin)));
    assert((XMVectorGetY(LengthMax) == XMVectorGetX(LengthMax)) && (XMVectorGetZ(LengthMax) == XMVectorGetX(LengthMax)));
    assert(XMVector3GreaterOrEqual(LengthMin, XMVectorZero()));
    assert(XMVector3GreaterOrEqual(LengthMax, XMVectorZero()));
    assert(XMVector3GreaterOrEqual(LengthMax, LengthMin));

    XMVECTOR LengthSq = XMVector3LengthSq(V);

    const XMVECTOR Zero = XMVectorZero();

    XMVECTOR RcpLength = XMVectorReciprocalSqrt(LengthSq);

    XMVECTOR InfiniteLength = XMVectorEqualInt(LengthSq, g_XMInfinity.v);
    XMVECTOR ZeroLength = XMVectorEqual(LengthSq, Zero);

    XMVECTOR Normal = XMVectorMultiply(V, RcpLength);

    XMVECTOR Length = XMVectorMultiply(LengthSq, RcpLength);

    XMVECTOR Select = XMVectorEqualInt(InfiniteLength, ZeroLength);
    Length = XMVectorSelect(LengthSq, Length, Select);
    Normal = XMVectorSelect(LengthSq, Normal, Select);

    XMVECTOR ControlMax = XMVectorGreater(Length, LengthMax);
    XMVECTOR ControlMin = XMVectorLess(Length, LengthMin);

    XMVECTOR ClampLength = XMVectorSelect(Length, LengthMax, ControlMax);
    ClampLength = XMVectorSelect(ClampLength, LengthMin, ControlMin);

    XMVECTOR Result = XMVectorMultiply(Normal, ClampLength);

    // Preserve the original vector (with no precision loss) if the length falls within the given range
    XMVECTOR Control = XMVectorEqualInt(ControlMax, ControlMin);
    Result = XMVectorSelect(Result, V, Control);

    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Reflect
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal
)
{
    // Result = Incident - (2 * dot(Incident, Normal)) * Normal

    XMVECTOR Result = XMVector3Dot(Incident, Normal);
    Result = XMVectorAdd(Result, Result);
    Result = XMVectorNegativeMultiplySubtract(Result, Normal, Incident);

    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Refract
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal, 
    float    RefractionIndex
)
{
    XMVECTOR Index = XMVectorReplicate(RefractionIndex);
    return XMVector3RefractV(Incident, Normal, Index);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3RefractV
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal, 
    FXMVECTOR RefractionIndex
)
{
    // Result = RefractionIndex * Incident - Normal * (RefractionIndex * dot(Incident, Normal) + 
    // sqrt(1 - RefractionIndex * RefractionIndex * (1 - dot(Incident, Normal) * dot(Incident, Normal))))

#if defined(_XM_NO_INTRINSICS_)

    const XMVECTOR  Zero = XMVectorZero();

    XMVECTOR IDotN = XMVector3Dot(Incident, Normal);

    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    XMVECTOR R = XMVectorNegativeMultiplySubtract(IDotN, IDotN, g_XMOne.v);
    R = XMVectorMultiply(R, RefractionIndex);
    R = XMVectorNegativeMultiplySubtract(R, RefractionIndex, g_XMOne.v);

    if (XMVector4LessOrEqual(R, Zero))
    {
        // Total internal reflection
        return Zero;
    }
    else
    {
        // R = RefractionIndex * IDotN + sqrt(R)
        R = XMVectorSqrt(R);
        R = XMVectorMultiplyAdd(RefractionIndex, IDotN, R);

        // Result = RefractionIndex * Incident - Normal * R
        XMVECTOR Result = XMVectorMultiply(RefractionIndex, Incident);
        Result = XMVectorNegativeMultiplySubtract(Normal, R, Result);

        return Result;
    }

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR IDotN = XMVector3Dot(Incident,Normal);

    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    __n128 R = vmlsq_f32( g_XMOne, IDotN, IDotN);
    R = vmulq_f32(R, RefractionIndex);
    R = vmlsq_f32(g_XMOne, R, RefractionIndex );

    __n128 vResult = vcleq_f32(R,g_XMZero);
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    if ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU )
    {
        // Total internal reflection
        vResult = g_XMZero;
    }
    else
    {
        // Sqrt(R)
        __n128 S0 = vrsqrteq_f32(R);
        __n128 P0 = vmulq_f32( R, S0 );
        __n128 R0 = vrsqrtsq_f32( P0, S0 );
        __n128 S1 = vmulq_f32( S0, R0 );
        __n128 P1 = vmulq_f32( R, S1 );
        __n128 R1 = vrsqrtsq_f32( P1, S1 );
        __n128 S2 = vmulq_f32( S1, R1 );
        R = vmulq_f32( R, S2 );
        // R = RefractionIndex * IDotN + sqrt(R)
        R = vmlaq_f32( R, RefractionIndex, IDotN );
        // Result = RefractionIndex * Incident - Normal * R
        vResult = vmulq_f32(RefractionIndex, Incident);
        vResult = vmlsq_f32( vResult, R, Normal );
    }
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Result = RefractionIndex * Incident - Normal * (RefractionIndex * dot(Incident, Normal) + 
    // sqrt(1 - RefractionIndex * RefractionIndex * (1 - dot(Incident, Normal) * dot(Incident, Normal))))
    XMVECTOR IDotN = XMVector3Dot(Incident, Normal);
    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    XMVECTOR R = _mm_mul_ps(IDotN, IDotN);
    R = _mm_sub_ps(g_XMOne,R);
    R = _mm_mul_ps(R, RefractionIndex);
    R = _mm_mul_ps(R, RefractionIndex);
    R = _mm_sub_ps(g_XMOne,R);

    XMVECTOR vResult = _mm_cmple_ps(R,g_XMZero);
    if (_mm_movemask_ps(vResult)==0x0f)
    {
        // Total internal reflection
        vResult = g_XMZero;
    }
    else
    {
        // R = RefractionIndex * IDotN + sqrt(R)
        R = _mm_sqrt_ps(R);
        vResult = _mm_mul_ps(RefractionIndex,IDotN);
        R = _mm_add_ps(R,vResult);
        // Result = RefractionIndex * Incident - Normal * R
        vResult = _mm_mul_ps(RefractionIndex, Incident);
        R = _mm_mul_ps(R,Normal);
        vResult = _mm_sub_ps(vResult,R);
    }
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Orthogonal
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Zero = XMVectorZero();
    XMVECTOR Z = XMVectorSplatZ(V);
    XMVECTOR YZYY = XMVectorSwizzle<XM_SWIZZLE_Y, XM_SWIZZLE_Z, XM_SWIZZLE_Y, XM_SWIZZLE_Y>(V);

    XMVECTOR NegativeV = XMVectorSubtract(Zero, V);

    XMVECTOR ZIsNegative = XMVectorLess(Z, Zero);
    XMVECTOR YZYYIsNegative = XMVectorLess(YZYY, Zero);

    XMVECTOR S = XMVectorAdd(YZYY, Z);
    XMVECTOR D = XMVectorSubtract(YZYY, Z);

    XMVECTOR Select = XMVectorEqualInt(ZIsNegative, YZYYIsNegative);

    XMVECTOR R0 = XMVectorPermute<XM_PERMUTE_1X, XM_PERMUTE_0X, XM_PERMUTE_0X, XM_PERMUTE_0X>(NegativeV, S);
    XMVECTOR R1 = XMVectorPermute<XM_PERMUTE_1X, XM_PERMUTE_0X, XM_PERMUTE_0X, XM_PERMUTE_0X>(V, D);

    return XMVectorSelect(R1, R0, Select);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3AngleBetweenNormalsEst
(
    FXMVECTOR N1, 
    FXMVECTOR N2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Result = XMVector3Dot(N1, N2);
    Result = XMVectorClamp(Result, g_XMNegativeOne.v, g_XMOne.v);
    Result = XMVectorACosEst(Result);
    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3AngleBetweenNormals
(
    FXMVECTOR N1, 
    FXMVECTOR N2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Result = XMVector3Dot(N1, N2);
    Result = XMVectorClamp(Result, g_XMNegativeOne.v, g_XMOne.v);
    Result = XMVectorACos(Result);
    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3AngleBetweenVectors
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR L1 = XMVector3ReciprocalLength(V1);
    XMVECTOR L2 = XMVector3ReciprocalLength(V2);

    XMVECTOR Dot = XMVector3Dot(V1, V2);

    L1 = XMVectorMultiply(L1, L2);

    XMVECTOR CosAngle = XMVectorMultiply(Dot, L1);
    CosAngle = XMVectorClamp(CosAngle, g_XMNegativeOne.v, g_XMOne.v);

    return XMVectorACos(CosAngle);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3LinePointDistance
(
    FXMVECTOR LinePoint1, 
    FXMVECTOR LinePoint2, 
    FXMVECTOR Point
)
{
    // Given a vector PointVector from LinePoint1 to Point and a vector
    // LineVector from LinePoint1 to LinePoint2, the scaled distance 
    // PointProjectionScale from LinePoint1 to the perpendicular projection
    // of PointVector onto the line is defined as:
    //
    //     PointProjectionScale = dot(PointVector, LineVector) / LengthSq(LineVector)

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR PointVector = XMVectorSubtract(Point, LinePoint1);
    XMVECTOR LineVector = XMVectorSubtract(LinePoint2, LinePoint1);

    XMVECTOR LengthSq = XMVector3LengthSq(LineVector);

    XMVECTOR PointProjectionScale = XMVector3Dot(PointVector, LineVector);
    PointProjectionScale = XMVectorDivide(PointProjectionScale, LengthSq);

    XMVECTOR DistanceVector = XMVectorMultiply(LineVector, PointProjectionScale);
    DistanceVector = XMVectorSubtract(PointVector, DistanceVector);

    return XMVector3Length(DistanceVector);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline void XMVector3ComponentsFromNormal
(
    XMVECTOR* pParallel, 
    XMVECTOR* pPerpendicular, 
    FXMVECTOR  V, 
    FXMVECTOR  Normal
)
{
    assert(pParallel != nullptr);
    assert(pPerpendicular != nullptr);

    XMVECTOR Scale = XMVector3Dot(V, Normal);

    XMVECTOR Parallel = XMVectorMultiply(Normal, Scale);

    *pParallel = Parallel;
    *pPerpendicular = XMVectorSubtract(V, Parallel);
}

//------------------------------------------------------------------------------
// Transform a vector using a rotation expressed as a unit quaternion

inline XMVECTOR XMVector3Rotate
(
    FXMVECTOR V, 
    FXMVECTOR RotationQuaternion
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR A = XMVectorSelect(g_XMSelect1110.v, V, g_XMSelect1110.v);
    XMVECTOR Q = XMQuaternionConjugate(RotationQuaternion);
    XMVECTOR Result = XMQuaternionMultiply(Q, A);
    return XMQuaternionMultiply(Result, RotationQuaternion);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Transform a vector using the inverse of a rotation expressed as a unit quaternion

inline XMVECTOR XMVector3InverseRotate
(
    FXMVECTOR V, 
    FXMVECTOR RotationQuaternion
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR A = XMVectorSelect(g_XMSelect1110.v, V, g_XMSelect1110.v);
    XMVECTOR Result = XMQuaternionMultiply(RotationQuaternion, A);
    XMVECTOR Q = XMQuaternionConjugate(RotationQuaternion);
    return XMQuaternionMultiply(Result, Q);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Transform
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Z = XMVectorSplatZ(V);
    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiplyAdd(Z, M.r[2], M.r[3]);
    Result = XMVectorMultiplyAdd(Y, M.r[1], Result);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    XMVECTOR vResult = vdupq_lane_f32( VL, 0 ); // X
    XMVECTOR vTemp = vdupq_lane_f32( VL, 1 ); // Y
    vResult = vmlaq_f32( M.r[3], vResult, M.r[0] );
    vResult = vmlaq_f32( vResult, vTemp, M.r[1] );
    vTemp = vdupq_lane_f32( vget_high_f32( V ), 0 ); // Z
    return vmlaq_f32( vResult, vTemp, M.r[2] );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,0,0,0));
    vResult = _mm_mul_ps(vResult,M.r[0]);
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    vTemp = _mm_mul_ps(vTemp,M.r[1]);
    vResult = _mm_add_ps(vResult,vTemp);
    vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,2,2,2));
    vTemp = _mm_mul_ps(vTemp,M.r[2]);
    vResult = _mm_add_ps(vResult,vTemp);
    vResult = _mm_add_ps(vResult,M.r[3]);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT4* XMVector3TransformStream
(
    XMFLOAT4*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT3* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t* pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];
    const XMVECTOR row2 = M.r[2];
    const XMVECTOR row3 = M.r[3];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat3((const XMFLOAT3*)pInputVector);
        XMVECTOR Z = XMVectorSplatZ(V);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiplyAdd(Z, row2, row3);
        Result = XMVectorMultiplyAdd(Y, row1, Result);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMStoreFloat4((XMFLOAT4*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}


//------------------------------------------------------------------------------

inline XMVECTOR XMVector3TransformCoord
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Z = XMVectorSplatZ(V);
    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiplyAdd(Z, M.r[2], M.r[3]);
    Result = XMVectorMultiplyAdd(Y, M.r[1], Result);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    XMVECTOR W = XMVectorSplatW(Result);
    return XMVectorDivide( Result, W );

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT3* XMVector3TransformCoordStream
(
    XMFLOAT3*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT3* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t*    pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];
    const XMVECTOR row2 = M.r[2];
    const XMVECTOR row3 = M.r[3];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat3((const XMFLOAT3*)pInputVector);
        XMVECTOR Z = XMVectorSplatZ(V);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiplyAdd(Z, row2, row3);
        Result = XMVectorMultiplyAdd(Y, row1, Result);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMVECTOR W = XMVectorSplatW(Result);

        Result = XMVectorDivide(Result, W);

        XMStoreFloat3((XMFLOAT3*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3TransformNormal
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Z = XMVectorSplatZ(V);
    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiply(Z, M.r[2]);
    Result = XMVectorMultiplyAdd(Y, M.r[1], Result);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    XMVECTOR vResult = vdupq_lane_f32( VL, 0 ); // X
    XMVECTOR vTemp = vdupq_lane_f32( VL, 1 ); // Y
    vResult = vmulq_f32( vResult, M.r[0] );
    vResult = vmlaq_f32( vResult, vTemp, M.r[1] );
    vTemp = vdupq_lane_f32( vget_high_f32( V ), 0 ); // Z
    return vmlaq_f32( vResult, vTemp, M.r[2] );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,0,0,0));
    vResult = _mm_mul_ps(vResult,M.r[0]);
    XMVECTOR vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    vTemp = _mm_mul_ps(vTemp,M.r[1]);
    vResult = _mm_add_ps(vResult,vTemp);
    vTemp = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,2,2,2));
    vTemp = _mm_mul_ps(vTemp,M.r[2]);
    vResult = _mm_add_ps(vResult,vTemp);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT3* XMVector3TransformNormalStream
(
    XMFLOAT3*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT3* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t* pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];
    const XMVECTOR row2 = M.r[2];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat3((const XMFLOAT3*)pInputVector);
        XMVECTOR Z = XMVectorSplatZ(V);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiply(Z, row2);
        Result = XMVectorMultiplyAdd(Y, row1, Result);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMStoreFloat3((XMFLOAT3*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Project
(
    FXMVECTOR V, 
    float    ViewportX, 
    float    ViewportY, 
    float    ViewportWidth, 
    float    ViewportHeight, 
    float    ViewportMinZ, 
    float    ViewportMaxZ, 
    CXMMATRIX Projection, 
    CXMMATRIX View, 
    CXMMATRIX World
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    const float HalfViewportWidth = ViewportWidth * 0.5f;
    const float HalfViewportHeight = ViewportHeight * 0.5f;

    XMVECTOR Scale = XMVectorSet(HalfViewportWidth, -HalfViewportHeight, ViewportMaxZ - ViewportMinZ, 0.0f);
    XMVECTOR Offset = XMVectorSet(ViewportX + HalfViewportWidth, ViewportY + HalfViewportHeight, ViewportMinZ, 0.0f);

    XMMATRIX Transform = XMMatrixMultiply(World, View);
    Transform = XMMatrixMultiply(Transform, Projection);

    XMVECTOR Result = XMVector3TransformCoord(V, Transform);

    Result = XMVectorMultiplyAdd(Result, Scale, Offset);

    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT3* XMVector3ProjectStream
(
    XMFLOAT3*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT3* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    float           ViewportX, 
    float           ViewportY, 
    float           ViewportWidth, 
    float           ViewportHeight, 
    float           ViewportMinZ, 
    float           ViewportMaxZ, 
    CXMMATRIX     Projection, 
    CXMMATRIX     View, 
    CXMMATRIX     World
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)

    const float HalfViewportWidth = ViewportWidth * 0.5f;
    const float HalfViewportHeight = ViewportHeight * 0.5f;

    XMVECTOR Scale = XMVectorSet(HalfViewportWidth, -HalfViewportHeight, ViewportMaxZ - ViewportMinZ, 1.0f);
    XMVECTOR Offset = XMVectorSet(ViewportX + HalfViewportWidth, ViewportY + HalfViewportHeight, ViewportMinZ, 0.0f);

    XMMATRIX Transform = XMMatrixMultiply(World, View);
    Transform = XMMatrixMultiply(Transform, Projection);

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t* pOutputVector = (uint8_t*)pOutputStream;

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat3((const XMFLOAT3*)pInputVector);

        XMVECTOR Result = XMVector3TransformCoord(V, Transform);
        Result = XMVectorMultiplyAdd(Result, Scale, Offset);

        XMStoreFloat3((XMFLOAT3*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Unproject
(
    FXMVECTOR V, 
    float     ViewportX, 
    float     ViewportY, 
    float     ViewportWidth, 
    float     ViewportHeight, 
    float     ViewportMinZ, 
    float     ViewportMaxZ, 
    CXMMATRIX Projection, 
    CXMMATRIX View, 
    CXMMATRIX World
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 D = { -1.0f, 1.0f, 0.0f, 0.0f };

    XMVECTOR Scale = XMVectorSet(ViewportWidth * 0.5f, -ViewportHeight * 0.5f, ViewportMaxZ - ViewportMinZ, 1.0f);
    Scale = XMVectorReciprocal(Scale);

    XMVECTOR Offset = XMVectorSet(-ViewportX, -ViewportY, -ViewportMinZ, 0.0f);
    Offset = XMVectorMultiplyAdd(Scale, Offset, D.v);

    XMMATRIX Transform = XMMatrixMultiply(World, View);
    Transform = XMMatrixMultiply(Transform, Projection);
    Transform = XMMatrixInverse(nullptr, Transform);

    XMVECTOR Result = XMVectorMultiplyAdd(V, Scale, Offset);

    return XMVector3TransformCoord(Result, Transform);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

_Use_decl_annotations_
inline XMFLOAT3* XMVector3UnprojectStream
(
    XMFLOAT3*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT3* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    float           ViewportX, 
    float           ViewportY, 
    float           ViewportWidth, 
    float           ViewportHeight, 
    float           ViewportMinZ, 
    float           ViewportMaxZ, 
    CXMMATRIX       Projection, 
    CXMMATRIX       View, 
    CXMMATRIX       World)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 D = { -1.0f, 1.0f, 0.0f, 0.0f };

    XMVECTOR Scale = XMVectorSet(ViewportWidth * 0.5f, -ViewportHeight * 0.5f, ViewportMaxZ - ViewportMinZ, 1.0f);
    Scale = XMVectorReciprocal(Scale);

    XMVECTOR Offset = XMVectorSet(-ViewportX, -ViewportY, -ViewportMinZ, 0.0f);
    Offset = XMVectorMultiplyAdd(Scale, Offset, D.v);

    XMMATRIX Transform = XMMatrixMultiply(World, View);
    Transform = XMMatrixMultiply(Transform, Projection);
    Transform = XMMatrixInverse(nullptr, Transform);

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t* pOutputVector = (uint8_t*)pOutputStream;

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat3((const XMFLOAT3*)pInputVector);

        XMVECTOR Result = XMVectorMultiplyAdd(V, Scale, Offset);

        Result = XMVector3TransformCoord(Result, Transform);

        XMStoreFloat3((XMFLOAT3*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * 4D Vector
 *
 ****************************************************************************/

//------------------------------------------------------------------------------
// Comparison operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline bool XMVector4Equal
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] == V2.vector4_f32[0]) && (V1.vector4_f32[1] == V2.vector4_f32[1]) && (V1.vector4_f32[2] == V2.vector4_f32[2]) && (V1.vector4_f32[3] == V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
    return ((_mm_movemask_ps(vTemp)==0x0f) != 0);
#else
    return XMComparisonAllTrue(XMVector4EqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector4EqualR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    uint32_t CR = 0;

    if ((V1.vector4_f32[0] == V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] == V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] == V2.vector4_f32[2]) &&
        (V1.vector4_f32[3] == V2.vector4_f32[3]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] != V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] != V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] != V2.vector4_f32[2]) &&
        (V1.vector4_f32[3] != V2.vector4_f32[3]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);

    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(V1,V2);
    int iTest = _mm_movemask_ps(vTemp);
    uint32_t CR = 0;
    if (iTest==0xf)     // All equal?
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (iTest==0)  // All not equal?
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector4EqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_u32[0] == V2.vector4_u32[0]) && (V1.vector4_u32[1] == V2.vector4_u32[1]) && (V1.vector4_u32[2] == V2.vector4_u32[2]) && (V1.vector4_u32[3] == V2.vector4_u32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    return ((_mm_movemask_ps(_mm_castsi128_ps(vTemp))==0xf) != 0);
#else
    return XMComparisonAllTrue(XMVector4EqualIntR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector4EqualIntR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t CR = 0;
    if (V1.vector4_u32[0] == V2.vector4_u32[0] && 
        V1.vector4_u32[1] == V2.vector4_u32[1] &&
        V1.vector4_u32[2] == V2.vector4_u32[2] &&
        V1.vector4_u32[3] == V2.vector4_u32[3])
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (V1.vector4_u32[0] != V2.vector4_u32[0] && 
        V1.vector4_u32[1] != V2.vector4_u32[1] &&
        V1.vector4_u32[2] != V2.vector4_u32[2] &&
        V1.vector4_u32[3] != V2.vector4_u32[3])
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);

    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    int iTest = _mm_movemask_ps(_mm_castsi128_ps(vTemp));
    uint32_t CR = 0;
    if (iTest==0xf)     // All equal?
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (iTest==0)  // All not equal?
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

inline bool XMVector4NearEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR Epsilon
)
{
#if defined(_XM_NO_INTRINSICS_)
    float dx, dy, dz, dw;

    dx = fabsf(V1.vector4_f32[0]-V2.vector4_f32[0]);
    dy = fabsf(V1.vector4_f32[1]-V2.vector4_f32[1]);
    dz = fabsf(V1.vector4_f32[2]-V2.vector4_f32[2]);
    dw = fabsf(V1.vector4_f32[3]-V2.vector4_f32[3]);
    return (((dx <= Epsilon.vector4_f32[0]) &&
            (dy <= Epsilon.vector4_f32[1]) &&
            (dz <= Epsilon.vector4_f32[2]) &&
            (dw <= Epsilon.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vDelta = vsubq_f32( V1, V2 );
    __n128 vResult = vacleq_f32( vDelta, Epsilon );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the difference
    XMVECTOR vDelta = _mm_sub_ps(V1,V2);
    // Get the absolute value of the difference
    XMVECTOR vTemp = _mm_setzero_ps();
    vTemp = _mm_sub_ps(vTemp,vDelta);
    vTemp = _mm_max_ps(vTemp,vDelta);
    vTemp = _mm_cmple_ps(vTemp,Epsilon);
    return ((_mm_movemask_ps(vTemp)==0xf) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector4NotEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] != V2.vector4_f32[0]) || (V1.vector4_f32[1] != V2.vector4_f32[1]) || (V1.vector4_f32[2] != V2.vector4_f32[2]) || (V1.vector4_f32[3] != V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) != 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpneq_ps(V1,V2);
    return ((_mm_movemask_ps(vTemp)) != 0);
#else
    return XMComparisonAnyFalse(XMVector4EqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector4NotEqualInt
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_u32[0] != V2.vector4_u32[0]) || (V1.vector4_u32[1] != V2.vector4_u32[1]) || (V1.vector4_u32[2] != V2.vector4_u32[2]) || (V1.vector4_u32[3] != V2.vector4_u32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vceqq_u32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) != 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    __m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1),_mm_castps_si128(V2));
    return ((_mm_movemask_ps(_mm_castsi128_ps(vTemp))!=0xF) != 0);
#else
    return XMComparisonAnyFalse(XMVector4EqualIntR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector4Greater
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] > V2.vector4_f32[0]) && (V1.vector4_f32[1] > V2.vector4_f32[1]) && (V1.vector4_f32[2] > V2.vector4_f32[2]) && (V1.vector4_f32[3] > V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgtq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
    return ((_mm_movemask_ps(vTemp)==0x0f) != 0);
#else
    return XMComparisonAllTrue(XMVector4GreaterR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector4GreaterR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t CR = 0;
    if (V1.vector4_f32[0] > V2.vector4_f32[0] && 
        V1.vector4_f32[1] > V2.vector4_f32[1] &&
        V1.vector4_f32[2] > V2.vector4_f32[2] &&
        V1.vector4_f32[3] > V2.vector4_f32[3])
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (V1.vector4_f32[0] <= V2.vector4_f32[0] && 
        V1.vector4_f32[1] <= V2.vector4_f32[1] &&
        V1.vector4_f32[2] <= V2.vector4_f32[2] &&
        V1.vector4_f32[3] <= V2.vector4_f32[3])
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgtq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);

    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    uint32_t CR = 0;
    XMVECTOR vTemp = _mm_cmpgt_ps(V1,V2);
    int iTest = _mm_movemask_ps(vTemp);
    if (iTest==0xf) {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector4GreaterOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] >= V2.vector4_f32[0]) && (V1.vector4_f32[1] >= V2.vector4_f32[1]) && (V1.vector4_f32[2] >= V2.vector4_f32[2]) && (V1.vector4_f32[3] >= V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgeq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    return ((_mm_movemask_ps(vTemp)==0x0f) != 0);
#else
    return XMComparisonAllTrue(XMVector4GreaterOrEqualR(V1, V2));
#endif
}

//------------------------------------------------------------------------------

inline uint32_t XMVector4GreaterOrEqualR
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    uint32_t CR = 0;
    if ((V1.vector4_f32[0] >= V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] >= V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] >= V2.vector4_f32[2]) &&
        (V1.vector4_f32[3] >= V2.vector4_f32[3]))
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ((V1.vector4_f32[0] < V2.vector4_f32[0]) && 
        (V1.vector4_f32[1] < V2.vector4_f32[1]) &&
        (V1.vector4_f32[2] < V2.vector4_f32[2]) &&
        (V1.vector4_f32[3] < V2.vector4_f32[3]))
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcgeq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    uint32_t r = vget_lane_u32(vTemp.val[1], 1);

    uint32_t CR = 0;
    if ( r == 0xFFFFFFFFU )
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if ( !r )
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#elif defined(_XM_SSE_INTRINSICS_)
    uint32_t CR = 0;
    XMVECTOR vTemp = _mm_cmpge_ps(V1,V2);
    int iTest = _mm_movemask_ps(vTemp);
    if (iTest==0x0f)
    {
        CR = XM_CRMASK_CR6TRUE;
    }
    else if (!iTest)
    {
        CR = XM_CRMASK_CR6FALSE;
    }
    return CR;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector4Less
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] < V2.vector4_f32[0]) && (V1.vector4_f32[1] < V2.vector4_f32[1]) && (V1.vector4_f32[2] < V2.vector4_f32[2]) && (V1.vector4_f32[3] < V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcltq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmplt_ps(V1,V2);
    return ((_mm_movemask_ps(vTemp)==0x0f) != 0);
#else
    return XMComparisonAllTrue(XMVector4GreaterR(V2, V1));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector4LessOrEqual
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1]) && (V1.vector4_f32[2] <= V2.vector4_f32[2]) && (V1.vector4_f32[3] <= V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vcleq_f32( V1, V2 );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmple_ps(V1,V2);
    return ((_mm_movemask_ps(vTemp)==0x0f) != 0);
#else
    return XMComparisonAllTrue(XMVector4GreaterOrEqualR(V2, V1));
#endif
}

//------------------------------------------------------------------------------

inline bool XMVector4InBounds
(
    FXMVECTOR V, 
    FXMVECTOR Bounds
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V.vector4_f32[0] <= Bounds.vector4_f32[0] && V.vector4_f32[0] >= -Bounds.vector4_f32[0]) && 
        (V.vector4_f32[1] <= Bounds.vector4_f32[1] && V.vector4_f32[1] >= -Bounds.vector4_f32[1]) &&
        (V.vector4_f32[2] <= Bounds.vector4_f32[2] && V.vector4_f32[2] >= -Bounds.vector4_f32[2]) &&
        (V.vector4_f32[3] <= Bounds.vector4_f32[3] && V.vector4_f32[3] >= -Bounds.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test if less than or equal
    __n128 vTemp1 = vcleq_f32(V,Bounds);
    // Negate the bounds
    __n128 vTemp2 = vnegq_f32(Bounds);
    // Test if greater or equal (Reversed)
    vTemp2 = vcleq_f32(vTemp2,V);
    // Blend answers
    vTemp1 = vandq_u32(vTemp1,vTemp2);
    // in bounds?
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTemp1), vget_high_u8(vTemp1));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test if less than or equal
    XMVECTOR vTemp1 = _mm_cmple_ps(V,Bounds);
    // Negate the bounds
    XMVECTOR vTemp2 = _mm_mul_ps(Bounds,g_XMNegativeOne);
    // Test if greater or equal (Reversed)
    vTemp2 = _mm_cmple_ps(vTemp2,V);
    // Blend answers
    vTemp1 = _mm_and_ps(vTemp1,vTemp2);
    // All in bounds?
    return ((_mm_movemask_ps(vTemp1)==0x0f) != 0);
#else
    return XMComparisonAllInBounds(XMVector4InBoundsR(V, Bounds));
#endif
}


//------------------------------------------------------------------------------

inline bool XMVector4IsNaN
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (XMISNAN(V.vector4_f32[0]) ||
            XMISNAN(V.vector4_f32[1]) ||
            XMISNAN(V.vector4_f32[2]) ||
            XMISNAN(V.vector4_f32[3]));
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Test against itself. NaN is always not equal
    __n128 vTempNan = vceqq_f32( V, V );
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTempNan), vget_high_u8(vTempNan));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    // If any are NaN, the mask is zero
    return ( vget_lane_u32(vTemp.val[1], 1) != 0xFFFFFFFFU );
#elif defined(_XM_SSE_INTRINSICS_)
    // Test against itself. NaN is always not equal
    XMVECTOR vTempNan = _mm_cmpneq_ps(V,V);
    // If any are NaN, the mask is non-zero
    return (_mm_movemask_ps(vTempNan)!=0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline bool XMVector4IsInfinite
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    return (XMISINF(V.vector4_f32[0]) ||
            XMISINF(V.vector4_f32[1]) ||
            XMISINF(V.vector4_f32[2]) ||
            XMISINF(V.vector4_f32[3]));

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Mask off the sign bit
    __n128 vTempInf = vandq_u32( V, g_XMAbsMask );
    // Compare to infinity
    vTempInf = vceqq_f32(vTempInf, g_XMInfinity );
    // If any are infinity, the signs are true.
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vTempInf), vget_high_u8(vTempInf));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    return ( vget_lane_u32(vTemp.val[1], 1) != 0 );
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bit
    XMVECTOR vTemp = _mm_and_ps(V,g_XMAbsMask);
    // Compare to infinity
    vTemp = _mm_cmpeq_ps(vTemp,g_XMInfinity);
    // If any are infinity, the signs are true.
    return (_mm_movemask_ps(vTemp) != 0);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// Computation operations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Dot
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] =
    Result.vector4_f32[1] =
    Result.vector4_f32[2] =
    Result.vector4_f32[3] = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1] + V1.vector4_f32[2] * V2.vector4_f32[2] + V1.vector4_f32[3] * V2.vector4_f32[3];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vTemp = vmulq_f32( V1, V2 );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    return vcombine_f32( v1, v1 );
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp2 = V2;
    XMVECTOR vTemp = _mm_mul_ps(V1,vTemp2);
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp,_MM_SHUFFLE(1,0,0,0)); // Copy X to the Z position and Y to the W position
    vTemp2 = _mm_add_ps(vTemp2,vTemp);          // Add Z = X+Z; W = Y+W;
    vTemp = _mm_shuffle_ps(vTemp,vTemp2,_MM_SHUFFLE(0,3,0,0));  // Copy W to the Z position
    vTemp = _mm_add_ps(vTemp,vTemp2);           // Add Z and W together
    return XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(2,2,2,2));    // Splat Z and return
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Cross
(
    FXMVECTOR V1, 
    FXMVECTOR V2, 
    FXMVECTOR V3
)
{
    // [ ((v2.z*v3.w-v2.w*v3.z)*v1.y)-((v2.y*v3.w-v2.w*v3.y)*v1.z)+((v2.y*v3.z-v2.z*v3.y)*v1.w),
    //   ((v2.w*v3.z-v2.z*v3.w)*v1.x)-((v2.w*v3.x-v2.x*v3.w)*v1.z)+((v2.z*v3.x-v2.x*v3.z)*v1.w),
    //   ((v2.y*v3.w-v2.w*v3.y)*v1.x)-((v2.x*v3.w-v2.w*v3.x)*v1.y)+((v2.x*v3.y-v2.y*v3.x)*v1.w),
    //   ((v2.z*v3.y-v2.y*v3.z)*v1.x)-((v2.z*v3.x-v2.x*v3.z)*v1.y)+((v2.y*v3.x-v2.x*v3.y)*v1.z) ]

#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;   

    Result.vector4_f32[0] = (((V2.vector4_f32[2]*V3.vector4_f32[3])-(V2.vector4_f32[3]*V3.vector4_f32[2]))*V1.vector4_f32[1])-(((V2.vector4_f32[1]*V3.vector4_f32[3])-(V2.vector4_f32[3]*V3.vector4_f32[1]))*V1.vector4_f32[2])+(((V2.vector4_f32[1]*V3.vector4_f32[2])-(V2.vector4_f32[2]*V3.vector4_f32[1]))*V1.vector4_f32[3]);
    Result.vector4_f32[1] = (((V2.vector4_f32[3]*V3.vector4_f32[2])-(V2.vector4_f32[2]*V3.vector4_f32[3]))*V1.vector4_f32[0])-(((V2.vector4_f32[3]*V3.vector4_f32[0])-(V2.vector4_f32[0]*V3.vector4_f32[3]))*V1.vector4_f32[2])+(((V2.vector4_f32[2]*V3.vector4_f32[0])-(V2.vector4_f32[0]*V3.vector4_f32[2]))*V1.vector4_f32[3]);
    Result.vector4_f32[2] = (((V2.vector4_f32[1]*V3.vector4_f32[3])-(V2.vector4_f32[3]*V3.vector4_f32[1]))*V1.vector4_f32[0])-(((V2.vector4_f32[0]*V3.vector4_f32[3])-(V2.vector4_f32[3]*V3.vector4_f32[0]))*V1.vector4_f32[1])+(((V2.vector4_f32[0]*V3.vector4_f32[1])-(V2.vector4_f32[1]*V3.vector4_f32[0]))*V1.vector4_f32[3]);
    Result.vector4_f32[3] = (((V2.vector4_f32[2]*V3.vector4_f32[1])-(V2.vector4_f32[1]*V3.vector4_f32[2]))*V1.vector4_f32[0])-(((V2.vector4_f32[2]*V3.vector4_f32[0])-(V2.vector4_f32[0]*V3.vector4_f32[2]))*V1.vector4_f32[1])+(((V2.vector4_f32[1]*V3.vector4_f32[0])-(V2.vector4_f32[0]*V3.vector4_f32[1]))*V1.vector4_f32[2]);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    const __n64 select = vget_low_f32( g_XMMaskX );

    // Term1: V2zwyz * V3wzwy
    const __n64 v2xy = vget_low_f32(V2);
    const __n64 v2zw = vget_high_f32(V2);
    const __n64 v2yx = vrev64_f32(v2xy);
    const __n64 v2wz = vrev64_f32(v2zw);
    const __n64 v2yz = vbsl_f32( select, v2yx, v2wz );

    const __n64 v3zw = vget_high_f32(V3);
    const __n64 v3wz = vrev64_f32(v3zw);
    const __n64 v3xy = vget_low_f32(V3);
    const __n64 v3wy = vbsl_f32( select, v3wz, v3xy );

    __n128 vTemp1 = vcombine_f32(v2zw,v2yz);
    __n128 vTemp2 = vcombine_f32(v3wz,v3wy);
    __n128 vResult = vmulq_f32( vTemp1, vTemp2 );

    // - V2wzwy * V3zwyz
    const __n64 v2wy = vbsl_f32( select, v2wz, v2xy );

    const __n64 v3yx = vrev64_f32(v3xy);
    const __n64 v3yz = vbsl_f32( select, v3yx, v3wz );

    vTemp1 = vcombine_f32(v2wz,v2wy);
    vTemp2 = vcombine_f32(v3zw,v3yz);
    vResult = vmlsq_f32( vResult, vTemp1, vTemp2 );

    // term1 * V1yxxx
    const __n64 v1xy = vget_low_f32(V1);
    const __n64 v1yx = vrev64_f32(v1xy);

    vTemp1 = vcombine_f32( v1yx, vdup_lane_f32( v1yx, 1 ) );
    vResult = vmulq_f32( vResult, vTemp1 );

    // Term2: V2ywxz * V3wxwx
    const __n64 v2yw = vrev64_f32(v2wy);
    const __n64 v2xz = vbsl_f32( select, v2xy, v2wz );

    const __n64 v3wx = vbsl_f32( select, v3wz, v3yx );

    vTemp1 = vcombine_f32(v2yw,v2xz);
    vTemp2 = vcombine_f32(v3wx,v3wx);
    __n128 vTerm = vmulq_f32( vTemp1, vTemp2 );

    // - V2wxwx * V3ywxz
    const __n64 v2wx = vbsl_f32( select, v2wz, v2yx );

    const __n64 v3yw = vrev64_f32(v3wy);
    const __n64 v3xz = vbsl_f32( select, v3xy, v3wz );

    vTemp1 = vcombine_f32(v2wx,v2wx);
    vTemp2 = vcombine_f32(v3yw,v3xz);
    vTerm = vmlsq_f32( vTerm, vTemp1, vTemp2 );

    // vResult - term2 * V1zzyy
    const __n64 v1zw = vget_high_f32(V1);

    vTemp1 = vcombine_f32( vdup_lane_f32(v1zw, 0), vdup_lane_f32(v1yx, 0) );
    vResult = vmlsq_f32( vResult, vTerm, vTemp1 );

    // Term3: V2yzxy * V3zxyx
    const __n64 v3zx = vrev64_f32(v3xz);

    vTemp1 = vcombine_f32(v2yz,v2xy);
    vTemp2 = vcombine_f32(v3zx,v3yx);
    vTerm = vmulq_f32( vTemp1, vTemp2 );

    // - V2zxyx * V3yzxy
    const __n64 v2zx = vrev64_f32(v2xz);

    vTemp1 = vcombine_f32(v2zx,v2yx);
    vTemp2 = vcombine_f32(v3yz,v3xy);
    vTerm = vmlsq_f32( vTerm, vTemp1, vTemp2 );

    // vResult + term3 * V1wwwz
    const __n64 v1wz = vrev64_f32(v1zw);

    vTemp1 = vcombine_f32( vdup_lane_f32( v1wz, 0 ), v1wz );
    return vmlaq_f32( vResult, vTerm, vTemp1 );
#elif defined(_XM_SSE_INTRINSICS_)
    // V2zwyz * V3wzwy
    XMVECTOR vResult = XM_PERMUTE_PS(V2,_MM_SHUFFLE(2,1,3,2));
    XMVECTOR vTemp3 = XM_PERMUTE_PS(V3,_MM_SHUFFLE(1,3,2,3));
    vResult = _mm_mul_ps(vResult,vTemp3);
    // - V2wzwy * V3zwyz
    XMVECTOR vTemp2 = XM_PERMUTE_PS(V2,_MM_SHUFFLE(1,3,2,3));
    vTemp3 = XM_PERMUTE_PS(vTemp3,_MM_SHUFFLE(1,3,0,1));
    vTemp2 = _mm_mul_ps(vTemp2,vTemp3);
    vResult = _mm_sub_ps(vResult,vTemp2);
    // term1 * V1yxxx
    XMVECTOR vTemp1 = XM_PERMUTE_PS(V1,_MM_SHUFFLE(0,0,0,1));
    vResult = _mm_mul_ps(vResult,vTemp1);

    // V2ywxz * V3wxwx
    vTemp2 = XM_PERMUTE_PS(V2,_MM_SHUFFLE(2,0,3,1));
    vTemp3 = XM_PERMUTE_PS(V3,_MM_SHUFFLE(0,3,0,3));
    vTemp3 = _mm_mul_ps(vTemp3,vTemp2);
    // - V2wxwx * V3ywxz
    vTemp2 = XM_PERMUTE_PS(vTemp2,_MM_SHUFFLE(2,1,2,1));
    vTemp1 = XM_PERMUTE_PS(V3,_MM_SHUFFLE(2,0,3,1));
    vTemp2 = _mm_mul_ps(vTemp2,vTemp1);
    vTemp3 = _mm_sub_ps(vTemp3,vTemp2);
    // vResult - temp * V1zzyy
    vTemp1 = XM_PERMUTE_PS(V1,_MM_SHUFFLE(1,1,2,2));
    vTemp1 = _mm_mul_ps(vTemp1,vTemp3);
    vResult = _mm_sub_ps(vResult,vTemp1);

    // V2yzxy * V3zxyx
    vTemp2 = XM_PERMUTE_PS(V2,_MM_SHUFFLE(1,0,2,1));
    vTemp3 = XM_PERMUTE_PS(V3,_MM_SHUFFLE(0,1,0,2));
    vTemp3 = _mm_mul_ps(vTemp3,vTemp2);
    // - V2zxyx * V3yzxy
    vTemp2 = XM_PERMUTE_PS(vTemp2,_MM_SHUFFLE(2,0,2,1));
    vTemp1 = XM_PERMUTE_PS(V3,_MM_SHUFFLE(1,0,2,1));
    vTemp1 = _mm_mul_ps(vTemp1,vTemp2);
    vTemp3 = _mm_sub_ps(vTemp3,vTemp1);
    // vResult + term * V1wwwz
    vTemp1 = XM_PERMUTE_PS(V1,_MM_SHUFFLE(2,3,3,3));
    vTemp3 = _mm_mul_ps(vTemp3,vTemp1);
    vResult = _mm_add_ps(vResult,vTemp3);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4LengthSq
(
    FXMVECTOR V
)
{
    return XMVector4Dot(V, V);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4ReciprocalLengthEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector4LengthSq(V);
    Result = XMVectorReciprocalSqrtEst(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot4
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    // Reciprocal sqrt (estimate)
    v2 = vrsqrte_f32( v1 );
    return vcombine_f32(v2, v2);
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(3,2,3,2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,0,0,0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp,vLengthSq,_MM_SHUFFLE(3,3,0,0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,2,2,2));
    // Get the reciprocal
    vLengthSq = _mm_rsqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4ReciprocalLength
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector4LengthSq(V);
    Result = XMVectorReciprocalSqrt(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot4
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    // Reciprocal sqrt
    __n64  S0 = vrsqrte_f32(v1);
    __n64  P0 = vmul_f32( v1, S0 );
    __n64  R0 = vrsqrts_f32( P0, S0 );
    __n64  S1 = vmul_f32( S0, R0 );
    __n64  P1 = vmul_f32( v1, S1 );
    __n64  R1 = vrsqrts_f32( P1, S1 );
    __n64 Result = vmul_f32( S1, R1 );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(3,2,3,2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,0,0,0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp,vLengthSq,_MM_SHUFFLE(3,3,0,0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,2,2,2));
    // Get the reciprocal
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    // Accurate!
    vLengthSq = _mm_div_ps(g_XMOne,vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4LengthEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector4LengthSq(V);
    Result = XMVectorSqrtEst(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot4
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    const __n64 zero = vdup_n_u32(0);
    __n64 VEqualsZero = vceq_f32( v1, zero );
    // Sqrt (estimate)
    __n64 Result = vrsqrte_f32( v1 );
    Result = vmul_f32( v1, Result );
    Result = vbsl_f32( VEqualsZero, zero, Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(3,2,3,2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,0,0,0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp,vLengthSq,_MM_SHUFFLE(3,3,0,0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,2,2,2));
    // Prepare for the division
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Length
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_) 

    XMVECTOR Result;

    Result = XMVector4LengthSq(V);
    Result = XMVectorSqrt(Result);

    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot4
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    const __n64 zero = vdup_n_u32(0);
    __n64 VEqualsZero = vceq_f32( v1, zero );
    // Sqrt
    __n64 S0 = vrsqrte_f32( v1 );
    __n64 P0 = vmul_f32( v1, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( v1, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    __n64 Result = vmul_f32( S1, R1 );
    Result = vmul_f32( v1, Result );
    Result = vbsl_f32( VEqualsZero, zero, Result );
    return vcombine_f32( Result, Result );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(3,2,3,2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,0,0,0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp,vLengthSq,_MM_SHUFFLE(3,3,0,0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,2,2,2));
    // Prepare for the division
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
// XMVector4NormalizeEst uses a reciprocal estimate and
// returns QNaN on zero and infinite vectors.

inline XMVECTOR XMVector4NormalizeEst
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector4ReciprocalLength(V);
    Result = XMVectorMultiply(V, Result);
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot4
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    // Reciprocal sqrt (estimate)
    v2 = vrsqrte_f32( v1 );
    // Normalize
    return vmulq_f32( V, vcombine_f32(v2,v2) );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(3,2,3,2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,0,0,0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp,vLengthSq,_MM_SHUFFLE(3,3,0,0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,2,2,2));
    // Get the reciprocal
    XMVECTOR vResult = _mm_rsqrt_ps(vLengthSq);
    // Reciprocal mul to perform the normalization
    vResult = _mm_mul_ps(vResult,V);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Normalize
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)
    float fLength;
    XMVECTOR vResult;

    vResult = XMVector4Length( V );
    fLength = vResult.vector4_f32[0];

    // Prevent divide by zero
    if (fLength > 0) {
        fLength = 1.0f/fLength;
    }
    
    vResult.vector4_f32[0] = V.vector4_f32[0]*fLength;
    vResult.vector4_f32[1] = V.vector4_f32[1]*fLength;
    vResult.vector4_f32[2] = V.vector4_f32[2]*fLength;
    vResult.vector4_f32[3] = V.vector4_f32[3]*fLength;
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    // Dot4
    __n128 vTemp = vmulq_f32( V, V );
    __n64 v1 = vget_low_f32( vTemp );
    __n64 v2 = vget_high_f32( vTemp );
    v1 = vpadd_f32( v1, v1 );
    v2 = vpadd_f32( v2, v2 );
    v1 = vadd_f32( v1, v2 );
    __n64 VEqualsZero = vceq_f32( v1, vdup_n_u32(0) );
    __n64 VEqualsInf = vceq_f32( v1, vget_low_f32(g_XMInfinity) );
    // Reciprocal sqrt (2 iterations of Newton-Raphson)
    __n64 S0 = vrsqrte_f32( v1 );
    __n64 P0 = vmul_f32( v1, S0 );
    __n64 R0 = vrsqrts_f32( P0, S0 );
    __n64 S1 = vmul_f32( S0, R0 );
    __n64 P1 = vmul_f32( v1, S1 );
    __n64 R1 = vrsqrts_f32( P1, S1 );
    v2 = vmul_f32( S1, R1 );
    // Normalize
    __n128 vResult = vmulq_f32( V, vcombine_f32(v2,v2) );
    vResult = vbslq_f32( vcombine_f32(VEqualsZero,VEqualsZero), vdupq_n_f32(0), vResult );
    return vbslq_f32( vcombine_f32(VEqualsInf,VEqualsInf), g_XMQNaN, vResult );
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(V,V);
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(3,2,3,2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(1,0,0,0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp,vLengthSq,_MM_SHUFFLE(3,3,0,0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq,vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq,_MM_SHUFFLE(2,2,2,2));
    // Prepare for the division
    XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
    // Create zero with a single instruction
    XMVECTOR vZeroMask = _mm_setzero_ps();
    // Test for a divide by zero (Must be FP to detect -0.0)
    vZeroMask = _mm_cmpneq_ps(vZeroMask,vResult);
    // Failsafe on zero (Or epsilon) length planes
    // If the length is infinity, set the elements to zero
    vLengthSq = _mm_cmpneq_ps(vLengthSq,g_XMInfinity);
    // Divide to perform the normalization
    vResult = _mm_div_ps(V,vResult);
    // Any that are infinity, set to zero
    vResult = _mm_and_ps(vResult,vZeroMask);
    // Select qnan or result based on infinite length
    XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq,g_XMQNaN);
    XMVECTOR vTemp2 = _mm_and_ps(vResult,vLengthSq);
    vResult = _mm_or_ps(vTemp1,vTemp2);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4ClampLength
(
    FXMVECTOR V, 
    float    LengthMin, 
    float    LengthMax
)
{
    XMVECTOR ClampMax = XMVectorReplicate(LengthMax);
    XMVECTOR ClampMin = XMVectorReplicate(LengthMin);

    return XMVector4ClampLengthV(V, ClampMin, ClampMax);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4ClampLengthV
(
    FXMVECTOR V, 
    FXMVECTOR LengthMin, 
    FXMVECTOR LengthMax
)
{
    assert((XMVectorGetY(LengthMin) == XMVectorGetX(LengthMin)) && (XMVectorGetZ(LengthMin) == XMVectorGetX(LengthMin)) && (XMVectorGetW(LengthMin) == XMVectorGetX(LengthMin)));
    assert((XMVectorGetY(LengthMax) == XMVectorGetX(LengthMax)) && (XMVectorGetZ(LengthMax) == XMVectorGetX(LengthMax)) && (XMVectorGetW(LengthMax) == XMVectorGetX(LengthMax)));
    assert(XMVector4GreaterOrEqual(LengthMin, XMVectorZero()));
    assert(XMVector4GreaterOrEqual(LengthMax, XMVectorZero()));
    assert(XMVector4GreaterOrEqual(LengthMax, LengthMin));

    XMVECTOR LengthSq = XMVector4LengthSq(V);

    const XMVECTOR Zero = XMVectorZero();

    XMVECTOR RcpLength = XMVectorReciprocalSqrt(LengthSq);

    XMVECTOR InfiniteLength = XMVectorEqualInt(LengthSq, g_XMInfinity.v);
    XMVECTOR ZeroLength = XMVectorEqual(LengthSq, Zero);

    XMVECTOR Normal = XMVectorMultiply(V, RcpLength);

    XMVECTOR Length = XMVectorMultiply(LengthSq, RcpLength);

    XMVECTOR Select = XMVectorEqualInt(InfiniteLength, ZeroLength);
    Length = XMVectorSelect(LengthSq, Length, Select);
    Normal = XMVectorSelect(LengthSq, Normal, Select);

    XMVECTOR ControlMax = XMVectorGreater(Length, LengthMax);
    XMVECTOR ControlMin = XMVectorLess(Length, LengthMin);

    XMVECTOR ClampLength = XMVectorSelect(Length, LengthMax, ControlMax);
    ClampLength = XMVectorSelect(ClampLength, LengthMin, ControlMin);

    XMVECTOR Result = XMVectorMultiply(Normal, ClampLength);

    // Preserve the original vector (with no precision loss) if the length falls within the given range
    XMVECTOR Control = XMVectorEqualInt(ControlMax, ControlMin);
    Result = XMVectorSelect(Result, V, Control);

    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Reflect
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal
)
{
    // Result = Incident - (2 * dot(Incident, Normal)) * Normal

    XMVECTOR Result = XMVector4Dot(Incident, Normal);
    Result = XMVectorAdd(Result, Result);
    Result = XMVectorNegativeMultiplySubtract(Result, Normal, Incident);

    return Result;
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Refract
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal, 
    float    RefractionIndex
)
{
    XMVECTOR Index = XMVectorReplicate(RefractionIndex);
    return XMVector4RefractV(Incident, Normal, Index);
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4RefractV
(
    FXMVECTOR Incident, 
    FXMVECTOR Normal, 
    FXMVECTOR RefractionIndex
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR        IDotN;
    XMVECTOR        R;
    const XMVECTOR  Zero = XMVectorZero();

    // Result = RefractionIndex * Incident - Normal * (RefractionIndex * dot(Incident, Normal) + 
    // sqrt(1 - RefractionIndex * RefractionIndex * (1 - dot(Incident, Normal) * dot(Incident, Normal))))

    IDotN = XMVector4Dot(Incident, Normal);

    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    R = XMVectorNegativeMultiplySubtract(IDotN, IDotN, g_XMOne.v);
    R = XMVectorMultiply(R, RefractionIndex);
    R = XMVectorNegativeMultiplySubtract(R, RefractionIndex, g_XMOne.v);

    if (XMVector4LessOrEqual(R, Zero))
    {
        // Total internal reflection
        return Zero;
    }
    else
    {
        XMVECTOR Result;

        // R = RefractionIndex * IDotN + sqrt(R)
        R = XMVectorSqrt(R);
        R = XMVectorMultiplyAdd(RefractionIndex, IDotN, R);

        // Result = RefractionIndex * Incident - Normal * R
        Result = XMVectorMultiply(RefractionIndex, Incident);
        Result = XMVectorNegativeMultiplySubtract(Normal, R, Result);

        return Result;
    }

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTOR IDotN = XMVector4Dot(Incident,Normal);

    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    __n128 R = vmlsq_f32( g_XMOne, IDotN, IDotN);
    R = vmulq_f32(R, RefractionIndex);
    R = vmlsq_f32(g_XMOne, R, RefractionIndex );

    __n128 vResult = vcleq_f32(R,g_XMZero);
    int8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
    vTemp = vzip_u16(vTemp.val[0], vTemp.val[1]);
    if ( vget_lane_u32(vTemp.val[1], 1) == 0xFFFFFFFFU )
    {
        // Total internal reflection
        vResult = g_XMZero;
    }
    else
    {
        // Sqrt(R)
        __n128 S0 = vrsqrteq_f32(R);
        __n128 P0 = vmulq_f32( R, S0 );
        __n128 R0 = vrsqrtsq_f32( P0, S0 );
        __n128 S1 = vmulq_f32( S0, R0 );
        __n128 P1 = vmulq_f32( R, S1 );
        __n128 R1 = vrsqrtsq_f32( P1, S1 );
        __n128 S2 = vmulq_f32( S1, R1 );
        R = vmulq_f32( R, S2 );
        // R = RefractionIndex * IDotN + sqrt(R)
        R = vmlaq_f32( R, RefractionIndex, IDotN );
        // Result = RefractionIndex * Incident - Normal * R
        vResult = vmulq_f32(RefractionIndex, Incident);
        vResult = vmlsq_f32( vResult, R, Normal );
    }
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR IDotN = XMVector4Dot(Incident,Normal);

    // R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
    XMVECTOR R = _mm_mul_ps(IDotN,IDotN);
    R = _mm_sub_ps(g_XMOne,R);
    R = _mm_mul_ps(R, RefractionIndex);
    R = _mm_mul_ps(R, RefractionIndex);
    R = _mm_sub_ps(g_XMOne,R);

    XMVECTOR vResult = _mm_cmple_ps(R,g_XMZero);
    if (_mm_movemask_ps(vResult)==0x0f)
    {
        // Total internal reflection
        vResult = g_XMZero;
    }
    else
    {
        // R = RefractionIndex * IDotN + sqrt(R)
        R = _mm_sqrt_ps(R);
        vResult = _mm_mul_ps(RefractionIndex, IDotN);
        R = _mm_add_ps(R,vResult);
        // Result = RefractionIndex * Incident - Normal * R
        vResult = _mm_mul_ps(RefractionIndex, Incident);
        R = _mm_mul_ps(R,Normal);
        vResult = _mm_sub_ps(vResult,R);
    }
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Orthogonal
(
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result.vector4_f32[0] = V.vector4_f32[2];
    Result.vector4_f32[1] = V.vector4_f32[3];
    Result.vector4_f32[2] = -V.vector4_f32[0];
    Result.vector4_f32[3] = -V.vector4_f32[1];
    return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Negate = { 1.f, 1.f, -1.f, -1.f };

    __n128 Result = vcombine_f32( vget_high_f32( V ), vget_low_f32( V ) );
    return vmulq_f32( Result, Negate );
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 FlipZW = {1.0f,1.0f,-1.0f,-1.0f};
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,0,3,2));
    vResult = _mm_mul_ps(vResult,FlipZW);
    return vResult;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4AngleBetweenNormalsEst
(
    FXMVECTOR N1, 
    FXMVECTOR N2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Result = XMVector4Dot(N1, N2);
    Result = XMVectorClamp(Result, g_XMNegativeOne.v, g_XMOne.v);
    Result = XMVectorACosEst(Result);
    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4AngleBetweenNormals
(
    FXMVECTOR N1, 
    FXMVECTOR N2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR Result = XMVector4Dot(N1, N2);
    Result = XMVectorClamp(Result, g_XMNegativeOne.v, g_XMOne.v);
    Result = XMVectorACos(Result);
    return Result;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4AngleBetweenVectors
(
    FXMVECTOR V1, 
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    XMVECTOR L1 = XMVector4ReciprocalLength(V1);
    XMVECTOR L2 = XMVector4ReciprocalLength(V2);

    XMVECTOR Dot = XMVector4Dot(V1, V2);

    L1 = XMVectorMultiply(L1, L2);

    XMVECTOR CosAngle = XMVectorMultiply(Dot, L1);
    CosAngle = XMVectorClamp(CosAngle, g_XMNegativeOne.v, g_XMOne.v);

    return XMVectorACos(CosAngle);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector4Transform
(
    FXMVECTOR V, 
    CXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)
    float fX = (M.m[0][0]*V.vector4_f32[0])+(M.m[1][0]*V.vector4_f32[1])+(M.m[2][0]*V.vector4_f32[2])+(M.m[3][0]*V.vector4_f32[3]);
    float fY = (M.m[0][1]*V.vector4_f32[0])+(M.m[1][1]*V.vector4_f32[1])+(M.m[2][1]*V.vector4_f32[2])+(M.m[3][1]*V.vector4_f32[3]);
    float fZ = (M.m[0][2]*V.vector4_f32[0])+(M.m[1][2]*V.vector4_f32[1])+(M.m[2][2]*V.vector4_f32[2])+(M.m[3][2]*V.vector4_f32[3]);
    float fW = (M.m[0][3]*V.vector4_f32[0])+(M.m[1][3]*V.vector4_f32[1])+(M.m[2][3]*V.vector4_f32[2])+(M.m[3][3]*V.vector4_f32[3]);
    XMVECTOR vResult = {
        fX,
        fY,
        fZ,
        fW
    };
    return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 VL = vget_low_f32( V );
    XMVECTOR vTemp1 = vdupq_lane_f32( VL, 0 ); // X
    XMVECTOR vTemp2 = vdupq_lane_f32( VL, 1 ); // Y
    XMVECTOR vResult = vmulq_f32( vTemp1, M.r[0] );
    vResult = vmlaq_f32( vResult, vTemp2, M.r[1] );
    __n64 VH = vget_high_f32( V );
    vTemp1 = vdupq_lane_f32( VH, 0 ); // Z
    vTemp2 = vdupq_lane_f32( VH, 1 ); // W
    vResult = vmlaq_f32( vResult, vTemp1, M.r[2] );
    return vmlaq_f32( vResult, vTemp2, M.r[3] );
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat x,y,z and w
    XMVECTOR vTempX = XM_PERMUTE_PS(V,_MM_SHUFFLE(0,0,0,0));
    XMVECTOR vTempY = XM_PERMUTE_PS(V,_MM_SHUFFLE(1,1,1,1));
    XMVECTOR vTempZ = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,2,2,2));
    XMVECTOR vTempW = XM_PERMUTE_PS(V,_MM_SHUFFLE(3,3,3,3));
    // Mul by the matrix
    vTempX = _mm_mul_ps(vTempX,M.r[0]);
    vTempY = _mm_mul_ps(vTempY,M.r[1]);
    vTempZ = _mm_mul_ps(vTempZ,M.r[2]);
    vTempW = _mm_mul_ps(vTempW,M.r[3]);
    // Add them all together
    vTempX = _mm_add_ps(vTempX,vTempY);
    vTempZ = _mm_add_ps(vTempZ,vTempW);
    vTempX = _mm_add_ps(vTempX,vTempZ);
    return vTempX;
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMFLOAT4* XMVector4TransformStream
(
    XMFLOAT4*       pOutputStream, 
    size_t          OutputStride, 
    const XMFLOAT4* pInputStream, 
    size_t          InputStride, 
    size_t          VectorCount, 
    CXMMATRIX       M
)
{
    assert(pOutputStream != nullptr);
    assert(pInputStream != nullptr);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS) || defined(_XM_ARM_NEON_INTRINSICS_)

    const uint8_t* pInputVector = (const uint8_t*)pInputStream;
    uint8_t* pOutputVector = (uint8_t*)pOutputStream;

    const XMVECTOR row0 = M.r[0];
    const XMVECTOR row1 = M.r[1];
    const XMVECTOR row2 = M.r[2];
    const XMVECTOR row3 = M.r[3];

    for (size_t i = 0; i < VectorCount; i++)
    {
        XMVECTOR V = XMLoadFloat4((const XMFLOAT4*)pInputVector);
        XMVECTOR W = XMVectorSplatW(V);
        XMVECTOR Z = XMVectorSplatZ(V);
        XMVECTOR Y = XMVectorSplatY(V);
        XMVECTOR X = XMVectorSplatX(V);

        XMVECTOR Result = XMVectorMultiply(W, row3);
        Result = XMVectorMultiplyAdd(Z, row2, Result);
        Result = XMVectorMultiplyAdd(Y, row1, Result);
        Result = XMVectorMultiplyAdd(X, row0, Result);

        XMStoreFloat4((XMFLOAT4*)pOutputVector, Result);

        pInputVector += InputStride; 
        pOutputVector += OutputStride;
    }

    return pOutputStream;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * XMVECTOR operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline XMVECTOR operator+ (FXMVECTOR V)
{
    return V;
}

//------------------------------------------------------------------------------

inline XMVECTOR operator- (FXMVECTOR V)
{
    return XMVectorNegate(V);
}

//------------------------------------------------------------------------------

inline XMVECTOR& operator+=
(
    XMVECTOR&       V1,
    FXMVECTOR       V2
)
{
    V1 = XMVectorAdd(V1, V2);
    return V1;
}

//------------------------------------------------------------------------------

inline XMVECTOR& operator-=
(
    XMVECTOR&       V1,
    FXMVECTOR       V2
)
{
    V1 = XMVectorSubtract(V1, V2);
    return V1;
}

//------------------------------------------------------------------------------

inline XMVECTOR& operator*=
(
    XMVECTOR&       V1,
    FXMVECTOR       V2
)
{
    V1 = XMVectorMultiply(V1, V2);
    return V1;
}

//------------------------------------------------------------------------------

inline XMVECTOR& operator/=
(
    XMVECTOR&       V1,
    FXMVECTOR       V2
)
{
    V1 = XMVectorDivide(V1,V2);
    return V1;
}

//------------------------------------------------------------------------------

inline XMVECTOR& operator*=
(
    XMVECTOR&   V,
    const float S
)
{
    V = XMVectorScale(V, S);
    return V;
}

//------------------------------------------------------------------------------

inline XMVECTOR& operator/=
(
    XMVECTOR&   V,
    const float S
)
{
    assert( S != 0.0f );
    V = XMVectorScale(V, 1.0f / S);
    return V;
}

//------------------------------------------------------------------------------

inline XMVECTOR operator+
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    return XMVectorAdd(V1, V2);
}

//------------------------------------------------------------------------------

inline XMVECTOR operator-
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    return XMVectorSubtract(V1, V2);
}

//------------------------------------------------------------------------------

inline XMVECTOR operator*
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    return XMVectorMultiply(V1, V2);
}

//------------------------------------------------------------------------------

inline XMVECTOR operator/
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    return XMVectorDivide(V1,V2);
}

//------------------------------------------------------------------------------

inline XMVECTOR operator*
(
    FXMVECTOR      V,
    const float    S
)
{
    return XMVectorScale(V, S);
}

//------------------------------------------------------------------------------

inline XMVECTOR operator/
(
    FXMVECTOR      V,
    const float    S
)
{
    assert( S != 0.0f );
    return XMVectorScale(V, 1.0f / S);
}

//------------------------------------------------------------------------------

inline XMVECTOR operator*
(
    float           S,
    FXMVECTOR  	    V
)
{
    return XMVectorScale(V, S);
}

#if defined(_XM_NO_INTRINSICS_)
#undef XMISNAN
#undef XMISINF
#endif


