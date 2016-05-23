//-------------------------------------------------------------------------------------
// DirectXPackedVector.inl -- SIMD C++ Math library
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

#define XM3_PACK_FACTOR                  (float)(1 << 22)
#define XM3_UNPACK_FACTOR_UNSIGNED       (float)(1 << 23)
#define XM3_UNPACK_FACTOR_SIGNED         XM3_PACK_FACTOR

#define XM3_UNPACK_UNSIGNEDN_OFFSET(BitsX, BitsY, BitsZ, BitsW) \
                                        {-XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsX)) - 1), \
                                         -XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsY)) - 1), \
                                         -XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsZ)) - 1), \
                                         -XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsW)) - 1)}

#define XM3_UNPACK_UNSIGNEDN_SCALE(BitsX, BitsY, BitsZ, BitsW) \
                                        {XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsX)) - 1), \
                                         XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsY)) - 1), \
                                         XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsZ)) - 1), \
                                         XM3_UNPACK_FACTOR_UNSIGNED / (float)((1 << (BitsW)) - 1)}

#define XM3_UNPACK_SIGNEDN_SCALE(BitsX, BitsY, BitsZ, BitsW) \
                                        {-XM3_UNPACK_FACTOR_SIGNED / (float)((1 << ((BitsX) - 1)) - 1), \
                                         -XM3_UNPACK_FACTOR_SIGNED / (float)((1 << ((BitsY) - 1)) - 1), \
                                         -XM3_UNPACK_FACTOR_SIGNED / (float)((1 << ((BitsZ) - 1)) - 1), \
                                         -XM3_UNPACK_FACTOR_SIGNED / (float)((1 << ((BitsW) - 1)) - 1)}

#define XM3_PACK_UNSIGNEDN_SCALE(BitsX, BitsY, BitsZ, BitsW) \
                                        {-(float)((1 << (BitsX)) - 1) / XM3_PACK_FACTOR, \
                                         -(float)((1 << (BitsY)) - 1) / XM3_PACK_FACTOR, \
                                         -(float)((1 << (BitsZ)) - 1) / XM3_PACK_FACTOR, \
                                         -(float)((1 << (BitsW)) - 1) / XM3_PACK_FACTOR}

#define XM3_PACK_SIGNEDN_SCALE(BitsX, BitsY, BitsZ, BitsW) \
                                        {-(float)((1 << ((BitsX) - 1)) - 1) / XM3_PACK_FACTOR, \
                                         -(float)((1 << ((BitsY) - 1)) - 1) / XM3_PACK_FACTOR, \
                                         -(float)((1 << ((BitsZ) - 1)) - 1) / XM3_PACK_FACTOR, \
                                         -(float)((1 << ((BitsW) - 1)) - 1) / XM3_PACK_FACTOR}

#define XM3_PACK_OFFSET                  XMVectorSplatConstant(3, 0)

/****************************************************************************
 *
 * Data conversion
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline float PackedVector::XMConvertHalfToFloat
(
    HALF Value
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_)

    uint32_t Mantissa;
    uint32_t Exponent;
    uint32_t Result;

    Mantissa = (uint32_t)(Value & 0x03FF);

    if ((Value & 0x7C00) != 0)  // The value is normalized
    {
        Exponent = (uint32_t)((Value >> 10) & 0x1F);
    }
    else if (Mantissa != 0)     // The value is denormalized
    {
        // Normalize the value in the resulting float
        Exponent = 1;

        do
        {
            Exponent--;
            Mantissa <<= 1;
        } while ((Mantissa & 0x0400) == 0);

        Mantissa &= 0x03FF;
    }
    else                        // The value is zero
    {
        Exponent = (uint32_t)-112;
    }

    Result = ((Value & 0x8000) << 16) | // Sign
             ((Exponent + 112) << 23) | // Exponent
             (Mantissa << 13);          // Mantissa

    return *(float*)&Result;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif
}

//------------------------------------------------------------------------------

inline float* PackedVector::XMConvertHalfToFloatStream
(
    float*      pOutputStream, 
    size_t      OutputStride, 
    const HALF* pInputStream, 
    size_t      InputStride, 
    size_t      HalfCount
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_)

    const uint8_t* pHalf = (const uint8_t*)pInputStream;
    uint8_t* pFloat = (uint8_t*)pOutputStream;

    assert(pOutputStream);
    assert(pInputStream);

    for (size_t i = 0; i < HalfCount; i++)
    {
        *(float*)pFloat = XMConvertHalfToFloat(*(const HALF*)pHalf);
        pHalf += InputStride;
        pFloat += OutputStride; 
    }

    return pOutputStream;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline PackedVector::HALF PackedVector::XMConvertFloatToHalf
(
    float Value
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_)
    uint32_t Result;

    uint32_t IValue = ((uint32_t *)(&Value))[0];
    uint32_t Sign = (IValue & 0x80000000U) >> 16U;
    IValue = IValue & 0x7FFFFFFFU;      // Hack off the sign

    if (IValue > 0x47FFEFFFU)
    {
        // The number is too large to be represented as a half.  Saturate to infinity.
        Result = 0x7FFFU;
    }
    else
    {
        if (IValue < 0x38800000U)
        {
            // The number is too small to be represented as a normalized half.
            // Convert it to a denormalized value.
            uint32_t Shift = 113U - (IValue >> 23U);
            IValue = (0x800000U | (IValue & 0x7FFFFFU)) >> Shift;
        }
        else
        {
            // Rebias the exponent to represent the value as a normalized half.
            IValue += 0xC8000000U;
        }

        Result = ((IValue + 0x0FFFU + ((IValue >> 13U) & 1U)) >> 13U)&0x7FFFU; 
    }
    return (HALF)(Result|Sign);

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif
}

//------------------------------------------------------------------------------

inline PackedVector::HALF* PackedVector::XMConvertFloatToHalfStream
(
    HALF*        pOutputStream, 
    size_t       OutputStride, 
    const float* pInputStream, 
    size_t       InputStride, 
    size_t       FloatCount
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_)

    uint8_t* pFloat = (uint8_t*)pInputStream;
    uint8_t* pHalf = (uint8_t*)pOutputStream;

    assert(pOutputStream);
    assert(pInputStream);

    for (size_t i = 0; i < FloatCount; i++)
    {
        *(HALF*)pHalf = XMConvertFloatToHalf(*(float*)pFloat);
        pFloat += InputStride; 
        pHalf += OutputStride;
    }
    return pOutputStream;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * Vector and matrix load operations
 *
 ****************************************************************************/

inline XMVECTOR PackedVector::XMLoadColor
(
    const XMCOLOR* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    assert(pSource);
    {
    // int32_t -> Float conversions are done in one instruction.
    // uint32_t -> Float calls a runtime function. Keep in int32_t
    int32_t iColor = (int32_t)(pSource->c);
    XMVECTOR vColor = {
        (float)((iColor >> 16) & 0xFF) * (1.0f/255.0f),
        (float)((iColor >> 8) & 0xFF) * (1.0f/255.0f),
        (float)(iColor & 0xFF) * (1.0f/255.0f),
        (float)((iColor >> 24) & 0xFF) * (1.0f/255.0f)
    };
    return vColor;
    }
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Splat the color in all four entries
    __m128i vInt = _mm_set1_epi32(pSource->c);
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vInt = _mm_and_si128(vInt,g_XMMaskA8R8G8B8);
    // a is unsigned! Flip the bit to convert the order to signed
    vInt = _mm_xor_si128(vInt,g_XMFlipA8R8G8B8);
    // Convert to floating point numbers
    XMVECTOR vTemp = _mm_cvtepi32_ps(vInt);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMFixAA8R8G8B8);
    // Convert 0-255 to 0.0f-1.0f
    return _mm_mul_ps(vTemp,g_XMNormalizeA8R8G8B8);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadHalf2
(
    const XMHALF2* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    assert(pSource);
    {
    XMVECTOR vResult = {
        XMConvertHalfToFloat(pSource->x),
        XMConvertHalfToFloat(pSource->y),
        0.0f,
        0.0f
    };
    return vResult;
    }
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    XMVECTOR vResult = {
        XMConvertHalfToFloat(pSource->x),
        XMConvertHalfToFloat(pSource->y),
        0.0f,
        0.0f
    };
    return vResult;

#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadShortN2
(
    const XMSHORTN2* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    assert(pSource);
    {
    XMVECTOR vResult = {
        (pSource->x == -32768) ? -1.f : ((float)pSource->x * (1.0f/32767.0f)),
        (pSource->y == -32768) ? -1.f : ((float)pSource->y * (1.0f/32767.0f)),
        0.0f,
        0.0f
    };
    return vResult;
    }

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // x needs to be sign extended
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // x - 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16);
    // Convert -1.0f - 1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMNormalizeX16Y16);
    // Clamp result (for case of -32768)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadShort2
(
    const XMSHORT2* pSource
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
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // x needs to be sign extended
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // x - 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16);
    // Y is 65536 too large
    return _mm_mul_ps(vTemp,g_XMFixupY16);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUShortN2
(
    const XMUSHORTN2* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x / 65535.0f;
    V.vector4_f32[1] = (float)pSource->y / 65535.0f;
    V.vector4_f32[2] = 0.f;
    V.vector4_f32[3] = 0.f;
    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 FixupY16 = {1.0f/65535.0f,1.0f/(65535.0f*65536.0f),0.0f,0.0f};
    static const XMVECTORF32 FixaddY16 = {0,32768.0f*65536.0f,0,0};
    assert(pSource);
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // y needs to be sign flipped
    vTemp = _mm_xor_ps(vTemp,g_XMFlipY);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // y + 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,FixaddY16);
    // Y is 65536 times too large
    vTemp = _mm_mul_ps(vTemp,FixupY16);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUShort2
(
    const XMUSHORT2* pSource
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
    static const XMVECTORF32 FixaddY16 = {0,32768.0f,0,0};
    assert(pSource);
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // y needs to be sign flipped
    vTemp = _mm_xor_ps(vTemp,g_XMFlipY);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // Y is 65536 times too large
    vTemp = _mm_mul_ps(vTemp,g_XMFixupY16);
    // y + 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,FixaddY16);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadByteN2
(
    const XMBYTEN2* pSource
)
{
    assert(pSource);
    {
    XMVECTOR vResult = {
        (pSource->x == -128) ? -1.f : ((float)pSource->x * (1.0f/127.0f)),
        (pSource->y == -128) ? -1.f : ((float)pSource->y * (1.0f/127.0f)),
        0.0f,
        0.0f
    };
    return vResult;
    }
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadByte2
(
    const XMBYTE2* pSource
)
{
    assert(pSource);
    {
    XMVECTOR vResult = {
        (float)pSource->x,
        (float)pSource->y,
        0.0f,
        0.0f
    };
    return vResult;
    }
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUByteN2
(
    const XMUBYTEN2* pSource
)
{
    assert(pSource);
    {
    XMVECTOR vResult = {
        (float)pSource->x * (1.0f/255.0f),
        (float)pSource->y * (1.0f/255.0f),
        0.0f,
        0.0f
    };
    return vResult;
    }
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUByte2
(
    const XMUBYTE2* pSource
)
{
    assert(pSource);
    {
    XMVECTOR vResult = {
        (float)pSource->x,
        (float)pSource->y,
        0.0f,
        0.0f
    };
    return vResult;
    }
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadU565
(
    const XMU565* pSource
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORI32 U565And = {0x1F,0x3F<<5,0x1F<<11,0};
    static const XMVECTORF32 U565Mul = {1.0f,1.0f/32.0f,1.0f/2048.f,0};
    assert(pSource);
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,U565And);
    // Convert to float
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Normalize x, y, and z
    vResult = _mm_mul_ps(vResult,U565Mul);
    return vResult;
#else
    XMVECTOR          V;
    uint32_t              Element;

    assert(pSource);

    Element = pSource->v & 0x1F;
    V.vector4_f32[0] = (float)Element;
    Element = (pSource->v >> 5) & 0x3F;
    V.vector4_f32[1] = (float)Element;
    Element = (pSource->v >> 11) & 0x1F;
    V.vector4_f32[2] = (float)Element;
    V.vector4_f32[3] = 0.f;
    return V;
#endif // !_XM_SSE_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadFloat3PK
(
    const XMFLOAT3PK* pSource
)
{
    __declspec(align(16)) uint32_t Result[4];
    uint32_t Mantissa;
    uint32_t Exponent;

    assert(pSource);

    // X Channel (6-bit mantissa)
    Mantissa = pSource->xm;

    if ( pSource->xe == 0x1f ) // INF or NAN
    {
        Result[0] = 0x7f800000 | (pSource->xm << 17);
    }
    else
    {
        if ( pSource->xe != 0 ) // The value is normalized
        {
            Exponent = pSource->xe;
        }
        else if (Mantissa != 0) // The value is denormalized
        {
            // Normalize the value in the resulting float
            Exponent = 1;
    
            do
            {
                Exponent--;
                Mantissa <<= 1;
            } while ((Mantissa & 0x40) == 0);
    
            Mantissa &= 0x3F;
        }
        else // The value is zero
        {
            Exponent = (uint32_t)-112;
        }
    
        Result[0] = ((Exponent + 112) << 23) | (Mantissa << 17);
    }

    // Y Channel (6-bit mantissa)
    Mantissa = pSource->ym;

    if ( pSource->ye == 0x1f ) // INF or NAN
    {
        Result[1] = 0x7f800000 | (pSource->ym << 17);
    }
    else
    {
        if ( pSource->ye != 0 ) // The value is normalized
        {
            Exponent = pSource->ye;
        }
        else if (Mantissa != 0) // The value is denormalized
        {
            // Normalize the value in the resulting float
            Exponent = 1;
    
            do
            {
                Exponent--;
                Mantissa <<= 1;
            } while ((Mantissa & 0x40) == 0);
    
            Mantissa &= 0x3F;
        }
        else // The value is zero
        {
            Exponent = (uint32_t)-112;
        }
    
        Result[1] = ((Exponent + 112) << 23) | (Mantissa << 17);
    }

    // Z Channel (5-bit mantissa)
    Mantissa = pSource->zm;

    if ( pSource->ze == 0x1f ) // INF or NAN
    {
        Result[2] = 0x7f800000 | (pSource->zm << 17);
    }
    else
    {
        if ( pSource->ze != 0 ) // The value is normalized
        {
            Exponent = pSource->ze;
        }
        else if (Mantissa != 0) // The value is denormalized
        {
            // Normalize the value in the resulting float
            Exponent = 1;
    
            do
            {
                Exponent--;
                Mantissa <<= 1;
            } while ((Mantissa & 0x20) == 0);
    
            Mantissa &= 0x1F;
        }
        else // The value is zero
        {
            Exponent = (uint32_t)-112;
        }

        Result[2] = ((Exponent + 112) << 23) | (Mantissa << 18);
    }

    return XMLoadFloat3A( (const XMFLOAT3A*)&Result );
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadFloat3SE
(
    const XMFLOAT3SE* pSource
)
{
    __declspec(align(16)) uint32_t Result[4];
    uint32_t Mantissa;
    uint32_t Exponent, ExpBits;

    assert(pSource);

    if ( pSource->e == 0x1f ) // INF or NAN
    {
        Result[0] = 0x7f800000 | (pSource->xm << 14);
        Result[1] = 0x7f800000 | (pSource->ym << 14);
        Result[2] = 0x7f800000 | (pSource->zm << 14);
    }
    else if ( pSource->e != 0 ) // The values are all normalized
    {
        Exponent = pSource->e;

        ExpBits = (Exponent + 112) << 23;

        Mantissa = pSource->xm;
        Result[0] = ExpBits | (Mantissa << 14);

        Mantissa = pSource->ym;
        Result[1] = ExpBits | (Mantissa << 14);

        Mantissa = pSource->zm;
        Result[2] = ExpBits | (Mantissa << 14);
    }
    else
    {
        // X Channel
        Mantissa = pSource->xm;

        if (Mantissa != 0) // The value is denormalized
        {
            // Normalize the value in the resulting float
            Exponent = 1;

            do
            {
                Exponent--;
                Mantissa <<= 1;
            } while ((Mantissa & 0x200) == 0);

            Mantissa &= 0x1FF;
        }
        else // The value is zero
        {
            Exponent = (uint32_t)-112;
        }

        Result[0] = ((Exponent + 112) << 23) | (Mantissa << 14);

        // Y Channel
        Mantissa = pSource->ym;

        if (Mantissa != 0) // The value is denormalized
        {
            // Normalize the value in the resulting float
            Exponent = 1;

            do
            {
                Exponent--;
                Mantissa <<= 1;
            } while ((Mantissa & 0x200) == 0);

            Mantissa &= 0x1FF;
        }
        else // The value is zero
        {
            Exponent = (uint32_t)-112;
        }

        Result[1] = ((Exponent + 112) << 23) | (Mantissa << 14);

        // Z Channel
        Mantissa = pSource->zm;

        if (Mantissa != 0) // The value is denormalized
        {
            // Normalize the value in the resulting float
            Exponent = 1;

            do
            {
                Exponent--;
                Mantissa <<= 1;
            } while ((Mantissa & 0x200) == 0);

            Mantissa &= 0x1FF;
        }
        else // The value is zero
        {
            Exponent = (uint32_t)-112;
        }

        Result[2] = ((Exponent + 112) << 23) | (Mantissa << 14);
    }

    return XMLoadFloat3A( (const XMFLOAT3A*)&Result );
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadHalf4
(
    const XMHALF4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    assert(pSource);
    {
    XMVECTOR vResult = {
        XMConvertHalfToFloat(pSource->x),
        XMConvertHalfToFloat(pSource->y),
        XMConvertHalfToFloat(pSource->z),
        XMConvertHalfToFloat(pSource->w)
    };
    return vResult;
    }
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    XMVECTOR vResult = {
        XMConvertHalfToFloat(pSource->x),
        XMConvertHalfToFloat(pSource->y),
        XMConvertHalfToFloat(pSource->z),
        XMConvertHalfToFloat(pSource->w)
    };
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadShortN4
(
    const XMSHORTN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    assert(pSource);
    {
    XMVECTOR vResult = {
        (pSource->x == -32768) ? -1.f : ((float)pSource->x * (1.0f/32767.0f)),
        (pSource->y == -32768) ? -1.f : ((float)pSource->y * (1.0f/32767.0f)),
        (pSource->z == -32768) ? -1.f : ((float)pSource->z * (1.0f/32767.0f)),
        (pSource->w == -32768) ? -1.f : ((float)pSource->w * (1.0f/32767.0f))
    };
    return vResult;
    }
#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(reinterpret_cast<const __m128 *>(&vIntd)[0],g_XMMaskX16Y16Z16W16);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16Z16W16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // x and z - 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16Z16W16);
    // Convert to -1.0f - 1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMNormalizeX16Y16Z16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    vTemp = _mm_shuffle_ps(vTemp,vTemp,_MM_SHUFFLE(3,1,2,0));
    // Clamp result (for case of -32768)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadShort4
(
    const XMSHORT4* pSource
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
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(reinterpret_cast<const __m128 *>(&vIntd)[0],g_XMMaskX16Y16Z16W16);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16Z16W16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // x and z - 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16Z16W16);
    // Fix y and w because they are 65536 too large
    vTemp = _mm_mul_ps(vTemp,g_XMFixupY16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    return _mm_shuffle_ps(vTemp,vTemp,_MM_SHUFFLE(3,1,2,0));
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUShortN4
(
    const XMUSHORTN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x / 65535.0f;
    V.vector4_f32[1] = (float)pSource->y / 65535.0f;
    V.vector4_f32[2] = (float)pSource->z / 65535.0f;
    V.vector4_f32[3] = (float)pSource->w / 65535.0f;

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    static const XMVECTORF32 FixupY16W16 = {1.0f/65535.0f,1.0f/65535.0f,1.0f/(65535.0f*65536.0f),1.0f/(65535.0f*65536.0f)};
    static const XMVECTORF32 FixaddY16W16  = {0,0,32768.0f*65536.0f,32768.0f*65536.0f};
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(reinterpret_cast<const __m128 *>(&vIntd)[0],g_XMMaskX16Y16Z16W16);
    // y and w are signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipZW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // y and w + 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,FixaddY16W16);
    // Fix y and w because they are 65536 too large
    vTemp = _mm_mul_ps(vTemp,FixupY16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    return _mm_shuffle_ps(vTemp,vTemp,_MM_SHUFFLE(3,1,2,0));
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUShort4
(
    const XMUSHORT4* pSource
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
    static const XMVECTORF32 FixaddY16W16  = {0,0,32768.0f,32768.0f};
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(reinterpret_cast<const __m128 *>(&vIntd)[0],g_XMMaskX16Y16Z16W16);
    // y and w are signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipZW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // Fix y and w because they are 65536 too large
    vTemp = _mm_mul_ps(vTemp,g_XMFixupY16W16);
    // y and w + 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,FixaddY16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    return _mm_shuffle_ps(vTemp,vTemp,_MM_SHUFFLE(3,1,2,0));
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadXDecN4
(
    const XMXDECN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR V;
    uint32_t Element;
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};

    assert(pSource);
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);

    Element = pSource->v & 0x3FF;
    V.vector4_f32[0] = (float)(int16_t)(Element | SignExtend[Element >> 9]) / 511.0f;
    Element = (pSource->v >> 10) & 0x3FF;
    V.vector4_f32[1] = (float)(int16_t)(Element | SignExtend[Element >> 9]) / 511.0f;
    Element = (pSource->v >> 20) & 0x3FF;
    V.vector4_f32[2] = (float)(int16_t)(Element | SignExtend[Element >> 9]) / 511.0f;
    V.vector4_f32[3] = (float)(pSource->v >> 30) / 3.0f;

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Splat the color in all four entries
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskA2B10G10R10);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipA2B10G10R10);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMFixAA2B10G10R10);
    // Convert 0-255 to 0.0f-1.0f
    return _mm_mul_ps(vTemp,g_XMNormalizeA2B10G10R10);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadXDec4
(
    const XMXDEC4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR          V;
    uint32_t              Element;
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};

    assert(pSource);
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);

    Element = pSource->v & 0x3FF;
    V.vector4_f32[0] = (float)(int16_t)(Element | SignExtend[Element >> 9]);
    Element = (pSource->v >> 10) & 0x3FF;
    V.vector4_f32[1] = (float)(int16_t)(Element | SignExtend[Element >> 9]);
    Element = (pSource->v >> 20) & 0x3FF;
    V.vector4_f32[2] = (float)(int16_t)(Element | SignExtend[Element >> 9]);
    V.vector4_f32[3] = (float)(pSource->v >> 30);

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);
    static const XMVECTORI32 XDec4Xor = {0x200, 0x200<<10, 0x200<<20, 0x80000000};
    static const XMVECTORF32 XDec4Add = {-512.0f,-512.0f*1024.0f,-512.0f*1024.0f*1024.0f,32768*65536.0f};
    assert(pSource);
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,XDec4Xor);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,XDec4Add);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMMulDec4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUDecN4
(
    const XMUDECN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR          V;
    uint32_t              Element;

    assert(pSource);

    Element = pSource->v & 0x3FF;
    V.vector4_f32[0] = (float)Element / 1023.0f;
    Element = (pSource->v >> 10) & 0x3FF;
    V.vector4_f32[1] = (float)Element / 1023.0f;
    Element = (pSource->v >> 20) & 0x3FF;
    V.vector4_f32[2] = (float)Element / 1023.0f;
    V.vector4_f32[3] = (float)(pSource->v >> 30) / 3.0f;

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    static const XMVECTORF32 UDecN4Mul = {1.0f/1023.0f,1.0f/(1023.0f*1024.0f),1.0f/(1023.0f*1024.0f*1024.0f),1.0f/(3.0f*1024.0f*1024.0f*1024.0f)};
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,UDecN4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUDec4
(
    const XMUDEC4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR          V;
    uint32_t              Element;

    assert(pSource);

    Element = pSource->v & 0x3FF;
    V.vector4_f32[0] = (float)Element;
    Element = (pSource->v >> 10) & 0x3FF;
    V.vector4_f32[1] = (float)Element;
    Element = (pSource->v >> 20) & 0x3FF;
    V.vector4_f32[2] = (float)Element;
    V.vector4_f32[3] = (float)(pSource->v >> 30);

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMMulDec4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadDecN4
(
    const XMDECN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR          V;
    uint32_t              Element;
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};
    static const uint32_t SignExtendW[] = {0x00000000, 0xFFFFFFFC};

    assert(pSource);
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);
    assert(((pSource->v >> 30) & 0x3) != 0x2);

    Element = pSource->v & 0x3FF;
    V.vector4_f32[0] = (float)(int16_t)(Element | SignExtend[Element >> 9]) / 511.0f;
    Element = (pSource->v >> 10) & 0x3FF;
    V.vector4_f32[1] = (float)(int16_t)(Element | SignExtend[Element >> 9]) / 511.0f;
    Element = (pSource->v >> 20) & 0x3FF;
    V.vector4_f32[2] = (float)(int16_t)(Element | SignExtend[Element >> 9]) / 511.0f;
    Element = pSource->v >> 30;
    V.vector4_f32[3] = (float)(int16_t)(Element | SignExtendW[Element >> 1]);

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pSource);
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);
    assert(((pSource->v >> 30) & 0x3) != 0x2);
    static const XMVECTORF32 DecN4Mul = {1.0f/511.0f,1.0f/(511.0f*1024.0f),1.0f/(511.0f*1024.0f*1024.0f),1.0f/(1024.0f*1024.0f*1024.0f)};
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorDec4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,DecN4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadDec4
(
    const XMDEC4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR          V;
    uint32_t              Element;
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};
    static const uint32_t SignExtendW[] = {0x00000000, 0xFFFFFFFC};

    assert(pSource);
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);
    assert(((pSource->v >> 30) & 0x3) != 0x2);

    Element = pSource->v & 0x3FF;
    V.vector4_f32[0] = (float)(int16_t)(Element | SignExtend[Element >> 9]);
    Element = (pSource->v >> 10) & 0x3FF;
    V.vector4_f32[1] = (float)(int16_t)(Element | SignExtend[Element >> 9]);
    Element = (pSource->v >> 20) & 0x3FF;
    V.vector4_f32[2] = (float)(int16_t)(Element | SignExtend[Element >> 9]);
    Element = pSource->v >> 30;
    V.vector4_f32[3] = (float)(int16_t)(Element | SignExtendW[Element >> 1]);

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    assert((pSource->v & 0x3FF) != 0x200);
    assert(((pSource->v >> 10) & 0x3FF) != 0x200);
    assert(((pSource->v >> 20) & 0x3FF) != 0x200);
    assert(((pSource->v >> 30) & 0x3) != 0x2);
    assert(pSource);
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorDec4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMMulDec4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUByteN4
(
    const XMUBYTEN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_f32[0] = (float)pSource->x / 255.0f;
    V.vector4_f32[1] = (float)pSource->y / 255.0f;
    V.vector4_f32[2] = (float)pSource->z / 255.0f;
    V.vector4_f32[3] = (float)pSource->w / 255.0f;

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadUByteN4Mul = {1.0f/255.0f,1.0f/(255.0f*256.0f),1.0f/(255.0f*65536.0f),1.0f/(255.0f*65536.0f*256.0f)};
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // w is signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // w + 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadUByteN4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUByte4
(
    const XMUBYTE4* pSource
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
    static const XMVECTORF32 LoadUByte4Mul = {1.0f,1.0f/256.0f,1.0f/65536.0f,1.0f/(65536.0f*256.0f)};
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // w is signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // w + 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadUByte4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadByteN4
(
    const XMBYTEN4* pSource
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR V;

    assert(pSource);

    V.vector4_f32[0] = (pSource->x == -128) ? -1.f : ((float)pSource->x / 127.0f);
    V.vector4_f32[1] = (pSource->y == -128) ? -1.f : ((float)pSource->y / 127.0f);
    V.vector4_f32[2] = (pSource->z == -128) ? -1.f : ((float)pSource->z / 127.0f);
    V.vector4_f32[3] = (pSource->w == -128) ? -1.f : ((float)pSource->w / 127.0f);

    return V;

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadByteN4Mul = {1.0f/127.0f,1.0f/(127.0f*256.0f),1.0f/(127.0f*65536.0f),1.0f/(127.0f*65536.0f*256.0f)};
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // x,y and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorByte4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // x, y and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddByte4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadByteN4Mul);
    // Clamp result (for case of -128)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadByte4
(
    const XMBYTE4* pSource
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
    static const XMVECTORF32 LoadByte4Mul = {1.0f,1.0f/256.0f,1.0f/65536.0f,1.0f/(65536.0f*256.0f)};
    assert(pSource);
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // x,y and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorByte4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vTemp)[0]);
    // x, y and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddByte4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadByte4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadUNibble4
(
     const XMUNIBBLE4* pSource
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORI32 UNibble4And = {0xF,0xF0,0xF00,0xF000};
    static const XMVECTORF32 UNibble4Mul = {1.0f,1.0f/16.f,1.0f/256.f,1.0f/4096.f};
    assert(pSource);
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,UNibble4And);
    // Convert to float
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Normalize x, y, and z
    vResult = _mm_mul_ps(vResult,UNibble4Mul);
    return vResult;
#else
    XMVECTOR          V;
    uint32_t              Element;

    assert(pSource);

    Element = pSource->v & 0xF;
    V.vector4_f32[0] = (float)Element;
    Element = (pSource->v >> 4) & 0xF;
    V.vector4_f32[1] = (float)Element;
    Element = (pSource->v >> 8) & 0xF;
    V.vector4_f32[2] = (float)Element;
    Element = (pSource->v >> 12) & 0xF;
    V.vector4_f32[3] = (float)Element;

    return V;
#endif // !_XM_SSE_INTRISICS_
}

//------------------------------------------------------------------------------

inline XMVECTOR PackedVector::XMLoadU555
(
     const XMU555* pSource
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORI32 U555And = {0x1F,0x1F<<5,0x1F<<10,0x8000};
    static const XMVECTORF32 U555Mul = {1.0f,1.0f/32.f,1.0f/1024.f,1.0f/32768.f};
    assert(pSource);
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,U555And);
    // Convert to float
    vResult = _mm_cvtepi32_ps(reinterpret_cast<const __m128i *>(&vResult)[0]);
    // Normalize x, y, and z
    vResult = _mm_mul_ps(vResult,U555Mul);
    return vResult;
#else
    XMVECTOR          V;
    uint32_t              Element;

    assert(pSource);

    Element = pSource->v & 0x1F;
    V.vector4_f32[0] = (float)Element;
    Element = (pSource->v >> 5) & 0x1F;
    V.vector4_f32[1] = (float)Element;
    Element = (pSource->v >> 10) & 0x1F;
    V.vector4_f32[2] = (float)Element;
    Element = (pSource->v >> 15) & 0x1;
    V.vector4_f32[3] = (float)Element;

    return V;
#endif // !_XM_SSE_INTRISICS_
}


/****************************************************************************
 *
 * Vector and matrix store operations
 *
 ****************************************************************************/

inline void PackedVector::XMStoreColor
(
    XMCOLOR* pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {255.0f, 255.0f, 255.0f, 255.0f};

    assert(pDestination);

    N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    pDestination->c = ((uint32_t)N.vector4_f32[3] << 24) |
                      ((uint32_t)N.vector4_f32[0] << 16) |
                      ((uint32_t)N.vector4_f32[1] <<  8) |
                      ((uint32_t)N.vector4_f32[2]);

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32  Scale = {255.0f,255.0f,255.0f,255.0f};
    // Set <0 to 0
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    // Set>1 to 1
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Convert to 0-255
    vResult = _mm_mul_ps(vResult,Scale);
    // Shuffle RGBA to ARGB
    vResult = _mm_shuffle_ps(vResult,vResult,_MM_SHUFFLE(3,0,1,2));
    // Convert to int 
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Mash to shorts
    vInt = _mm_packs_epi32(vInt,vInt);
    // Mash to bytes
    vInt = _mm_packus_epi16(vInt,vInt);
    // Store the color
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->c),reinterpret_cast<__m128 *>(&vInt)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreHalf2
(
    XMHALF2* pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    assert(pDestination);

    pDestination->x = XMConvertFloatToHalf(V.vector4_f32[0]);
    pDestination->y = XMConvertFloatToHalf(V.vector4_f32[1]);

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    pDestination->x = XMConvertFloatToHalf(XMVectorGetX(V));
    pDestination->y = XMConvertFloatToHalf(XMVectorGetY(V));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreShortN2
(
    XMSHORTN2* pDestination, 
    FXMVECTOR   V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR N;
    static const XMVECTORF32  Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    assert(pDestination);

    N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    vResult = _mm_mul_ps(vResult,Scale);
    __m128i vResulti = _mm_cvtps_epi32(vResult);
    vResulti = _mm_packs_epi32(vResulti,vResulti);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->x),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreShort2
(
    XMSHORT2* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTOR  Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32 Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,Min);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Pack the ints into shorts
    vInt = _mm_packs_epi32(vInt,vInt);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->x),reinterpret_cast<const __m128 *>(&vInt)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUShortN2
(
    XMUSHORTN2* pDestination, 
    FXMVECTOR    V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    assert(pDestination);

    N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);
    N = XMVectorTruncate(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 Scale = {65535.0f, 65535.0f, 65535.0f, 65535.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    vResult = _mm_mul_ps(vResult,Scale);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Since the SSE pack instruction clamps using signed rules,
    // manually extract the values to store them to memory
    pDestination->x = static_cast<int16_t>(_mm_extract_epi16(vInt,0));
    pDestination->y = static_cast<int16_t>(_mm_extract_epi16(vInt,2));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUShort2
(
    XMUSHORT2* pDestination, 
    FXMVECTOR   V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32  Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Since the SSE pack instruction clamps using signed rules,
    // manually extract the values to store them to memory
    pDestination->x = static_cast<int16_t>(_mm_extract_epi16(vInt,0));
    pDestination->y = static_cast<int16_t>(_mm_extract_epi16(vInt,2));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreByteN2
(
    XMBYTEN2* pDestination, 
    FXMVECTOR   V
)
{
    XMVECTOR N;
    XMFLOAT4A tmp;
    static const XMVECTORF32  Scale = {127.0f, 127.0f, 127.0f, 127.0f};

    assert(pDestination);

    N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int8_t)tmp.x;
    pDestination->y = (int8_t)tmp.y;
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreByte2
(
    XMBYTE2* pDestination, 
    FXMVECTOR  V
)
{
    XMVECTOR               N;
    XMFLOAT4A              tmp;
    static const XMVECTOR  Min = {-127.0f, -127.0f, -127.0f, -127.0f};
    static const XMVECTOR  Max = {127.0f, 127.0f, 127.0f, 127.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int8_t)tmp.x;
    pDestination->y = (int8_t)tmp.y;
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUByteN2
(
    XMUBYTEN2* pDestination, 
    FXMVECTOR    V
)
{
    XMVECTOR               N;
    XMFLOAT4A              tmp;
    static const XMVECTORF32  Scale = {255.0f, 255.0f, 255.0f, 255.0f};

    assert(pDestination);

    N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);
    N = XMVectorTruncate(N);

    XMStoreFloat4A( &tmp, N );

    pDestination->x = (uint8_t)tmp.x;
    pDestination->y = (uint8_t)tmp.y;
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUByte2
(
    XMUBYTE2* pDestination, 
    FXMVECTOR   V
)
{
    XMVECTOR               N;
    static const XMVECTOR  Max = {255.0f, 255.0f, 255.0f, 255.0f};
    XMFLOAT4A              tmp;

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    XMStoreFloat4A( &tmp, N );

    pDestination->x = (uint8_t)tmp.x;
    pDestination->y = (uint8_t)tmp.y;
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreU565
(
    XMU565* pDestination,
    FXMVECTOR V
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32  Max = {31.0f, 63.0f, 31.0f, 0.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // No SSE operations will write to 16-bit values, so we have to extract them manually
    uint16_t x = static_cast<uint16_t>(_mm_extract_epi16(vInt,0));
    uint16_t y = static_cast<uint16_t>(_mm_extract_epi16(vInt,2));
    uint16_t z = static_cast<uint16_t>(_mm_extract_epi16(vInt,4));
    pDestination->v = ((z & 0x1F) << 11) |
                      ((y & 0x3F) << 5) |
                      ((x & 0x1F));
#else
    XMVECTOR               N;
    static const XMVECTORF32  Max = {31.0f, 63.0f, 31.0f, 0.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max.v);
    N = XMVectorRound(N);

    pDestination->v = (((uint16_t)N.vector4_f32[2] & 0x1F) << 11) |
                      (((uint16_t)N.vector4_f32[1] & 0x3F) << 5) |
                      (((uint16_t)N.vector4_f32[0] & 0x1F));
#endif !_XM_SSE_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreFloat3PK
(
    XMFLOAT3PK* pDestination,
    FXMVECTOR V
)
{
    __declspec(align(16)) uint32_t IValue[4];
    uint32_t I, Sign, j;
    uint32_t Result[3];

    assert(pDestination);

    XMStoreFloat3A( (XMFLOAT3A*)&IValue, V );

    // X & Y Channels (5-bit exponent, 6-bit mantissa)
    for(j=0; j < 2; ++j)
    {
        Sign = IValue[j] & 0x80000000;
        I = IValue[j] & 0x7FFFFFFF;

        if ((I & 0x7F800000) == 0x7F800000)
        {
            // INF or NAN
            Result[j] = 0x7c0;
            if (( I & 0x7FFFFF ) != 0)
            {
                Result[j] = 0x7c0 | (((I>>17)|(I>11)|(I>>6)|(I))&0x3f);
            }
            else if ( Sign )
            {
                // -INF is clamped to 0 since 3PK is positive only
                Result[j] = 0;
            }
        }
        else if ( Sign )
        {
            // 3PK is positive only, so clamp to zero
            Result[j] = 0;
        }
        else if (I > 0x477E0000U)
        {
            // The number is too large to be represented as a float11, set to max
            Result[j] = 0x7BF;
        }
        else
        {
            if (I < 0x38800000U)
            {
                // The number is too small to be represented as a normalized float11
                // Convert it to a denormalized value.
                uint32_t Shift = 113U - (I >> 23U);
                I = (0x800000U | (I & 0x7FFFFFU)) >> Shift;
            }
            else
            {
                // Rebias the exponent to represent the value as a normalized float11
                I += 0xC8000000U;
            }
     
            Result[j] = ((I + 0xFFFFU + ((I >> 17U) & 1U)) >> 17U)&0x7ffU;
        }
    }

    // Z Channel (5-bit exponent, 5-bit mantissa)
    Sign = IValue[2] & 0x80000000;
    I = IValue[2] & 0x7FFFFFFF;

    if ((I & 0x7F800000) == 0x7F800000)
    {
        // INF or NAN
        Result[2] = 0x3e0;
        if ( I & 0x7FFFFF )
        {
            Result[2] = 0x3e0 | (((I>>18)|(I>13)|(I>>3)|(I))&0x1f);
        }
        else if ( Sign )
        {
            // -INF is clamped to 0 since 3PK is positive only
            Result[2] = 0;
        }
    }
    else if ( Sign )
    {
        // 3PK is positive only, so clamp to zero
        Result[2] = 0;
    }
    else if (I > 0x477C0000U)
    {
        // The number is too large to be represented as a float10, set to max
        Result[2] = 0x3df;
    }
    else
    {
        if (I < 0x38800000U)
        {
            // The number is too small to be represented as a normalized float10
            // Convert it to a denormalized value.
            uint32_t Shift = 113U - (I >> 23U);
            I = (0x800000U | (I & 0x7FFFFFU)) >> Shift;
        }
        else
        {
            // Rebias the exponent to represent the value as a normalized float10
            I += 0xC8000000U;
        }
     
        Result[2] = ((I + 0x1FFFFU + ((I >> 18U) & 1U)) >> 18U)&0x3ffU;
    }

    // Pack Result into memory
    pDestination->v = (Result[0] & 0x7ff)
                      | ( (Result[1] & 0x7ff) << 11 )
                      | ( (Result[2] & 0x3ff) << 22 );
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreFloat3SE
(
    XMFLOAT3SE* pDestination,
    FXMVECTOR V
)
{
    __declspec(align(16)) uint32_t IValue[4];
    uint32_t I, Sign, j, T;
    uint32_t Frac[3];
    uint32_t Exp[3];
    

    assert(pDestination);

    XMStoreFloat3A( (XMFLOAT3A*)&IValue, V );

    // X, Y, Z Channels (5-bit exponent, 9-bit mantissa)
    for(j=0; j < 3; ++j)
    {
        Sign = IValue[j] & 0x80000000;
        I = IValue[j] & 0x7FFFFFFF;

        if ((I & 0x7F800000) == 0x7F800000)
        {
            // INF or NAN
            Exp[j] = 0x1f;
            if (( I & 0x7FFFFF ) != 0)
            {
                Frac[j] = ((I>>14)|(I>5)|(I))&0x1ff;
            }
            else if ( Sign )
            {
                // -INF is clamped to 0 since 3SE is positive only
                Exp[j] = Frac[j] = 0;
            }
        }
        else if ( Sign )
        {
            // 3SE is positive only, so clamp to zero
            Exp[j] = Frac[j] = 0;
        }
        else if (I > 0x477FC000U)
        {
            // The number is too large, set to max
            Exp[j] = 0x1e;
            Frac[j] = 0x1ff;
        }
        else
        {
            if (I < 0x38800000U)
            {
                // The number is too small to be represented as a normalized float11
                // Convert it to a denormalized value.
                uint32_t Shift = 113U - (I >> 23U);
                I = (0x800000U | (I & 0x7FFFFFU)) >> Shift;
            }
            else
            {
                // Rebias the exponent to represent the value as a normalized float11
                I += 0xC8000000U;
            }
     
            T = ((I + 0x1FFFU + ((I >> 14U) & 1U)) >> 14U)&0x3fffU;

            Exp[j] = (T & 0x3E00) >> 9;
            Frac[j] = T & 0x1ff;
        }
    }

    // Adjust to a shared exponent
    T = XMMax( Exp[0], XMMax( Exp[1], Exp[2] ) );

    Frac[0] = Frac[0] >> (T - Exp[0]);
    Frac[1] = Frac[1] >> (T - Exp[1]);
    Frac[2] = Frac[2] >> (T - Exp[2]);

    // Store packed into memory
    pDestination->xm = Frac[0];
    pDestination->ym = Frac[1];
    pDestination->zm = Frac[2];
    pDestination->e = T;
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreHalf4
(
    XMHALF4* pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_) 

    assert(pDestination);

    pDestination->x = XMConvertFloatToHalf(V.vector4_f32[0]);
    pDestination->y = XMConvertFloatToHalf(V.vector4_f32[1]);
    pDestination->z = XMConvertFloatToHalf(V.vector4_f32[2]);
    pDestination->w = XMConvertFloatToHalf(V.vector4_f32[3]);

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    pDestination->x = XMConvertFloatToHalf(XMVectorGetX(V));
    pDestination->y = XMConvertFloatToHalf(XMVectorGetY(V));
    pDestination->z = XMConvertFloatToHalf(XMVectorGetZ(V));
    pDestination->w = XMConvertFloatToHalf(XMVectorGetW(V));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreShortN4
(
    XMSHORTN4* pDestination, 
    FXMVECTOR   V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    assert(pDestination);

    N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];
    pDestination->z = (int16_t)N.vector4_f32[2];
    pDestination->w = (int16_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    vResult = _mm_mul_ps(vResult,Scale);
    __m128i vResulti = _mm_cvtps_epi32(vResult);
    vResulti = _mm_packs_epi32(vResulti,vResulti);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->x),reinterpret_cast<const __m128d *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreShort4
(
    XMSHORT4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTOR  Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];
    pDestination->z = (int16_t)N.vector4_f32[2];
    pDestination->w = (int16_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32  Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,Min);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Pack the ints into shorts
    vInt = _mm_packs_epi32(vInt,vInt);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->x),reinterpret_cast<const __m128d *>(&vInt)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUShortN4
(
    XMUSHORTN4* pDestination, 
    FXMVECTOR    V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    assert(pDestination);

    N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);
    N = XMVectorTruncate(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];
    pDestination->z = (int16_t)N.vector4_f32[2];
    pDestination->w = (int16_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 Scale = {65535.0f, 65535.0f, 65535.0f, 65535.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    vResult = _mm_mul_ps(vResult,Scale);
    // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Since the SSE pack instruction clamps using signed rules,
    // manually extract the values to store them to memory
    pDestination->x = static_cast<int16_t>(_mm_extract_epi16(vInt,0));
    pDestination->y = static_cast<int16_t>(_mm_extract_epi16(vInt,2));
    pDestination->z = static_cast<int16_t>(_mm_extract_epi16(vInt,4));
    pDestination->w = static_cast<int16_t>(_mm_extract_epi16(vInt,6));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUShort4
(
    XMUSHORT4* pDestination, 
    FXMVECTOR   V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    pDestination->x = (int16_t)N.vector4_f32[0];
    pDestination->y = (int16_t)N.vector4_f32[1];
    pDestination->z = (int16_t)N.vector4_f32[2];
    pDestination->w = (int16_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32  Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Since the SSE pack instruction clamps using signed rules,
    // manually extract the values to store them to memory
    pDestination->x = static_cast<int16_t>(_mm_extract_epi16(vInt,0));
    pDestination->y = static_cast<int16_t>(_mm_extract_epi16(vInt,2));
    pDestination->z = static_cast<int16_t>(_mm_extract_epi16(vInt,4));
    pDestination->w = static_cast<int16_t>(_mm_extract_epi16(vInt,6));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreXDecN4
(
    XMXDECN4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Min = {-1.0f, -1.0f, -1.0f, 0.0f};
    static const XMVECTORF32  Scale = {511.0f, 511.0f, 511.0f, 3.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    pDestination->v = ((uint32_t)N.vector4_f32[3] << 30) |
                       (((int32_t)N.vector4_f32[2] & 0x3FF) << 20) |
                       (((int32_t)N.vector4_f32[1] & 0x3FF) << 10) |
                       (((int32_t)N.vector4_f32[0] & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Min = {-1.0f, -1.0f, -1.0f, 0.0f};
    static const XMVECTORF32 Scale = {511.0f, 511.0f*1024.0f, 511.0f*1048576.0f,3.0f*536870912.0f};
    static const XMVECTORI32 ScaleMask = {0x3FF,0x3FF<<10,0x3FF<<20,0x3<<29};
    assert(pDestination);
    XMVECTOR vResult = _mm_max_ps(V,Min);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,Scale);
    // Convert to int (W is unsigned)
    __m128i vResulti = _mm_cvtps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,ScaleMask);
    // To fix W, add itself to shift it up to <<30 instead of <<29
    __m128i vResultw = _mm_and_si128(vResulti,g_XMMaskW);
    vResulti = _mm_add_epi32(vResulti,vResultw);
    // Do a horizontal or of all 4 entries
    vResult = _mm_shuffle_ps(reinterpret_cast<const __m128 *>(&vResulti)[0],reinterpret_cast<const __m128 *>(&vResulti)[0],_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,reinterpret_cast<const __m128i *>(&vResult)[0]);
    vResult = _mm_shuffle_ps(vResult,vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,reinterpret_cast<const __m128i *>(&vResult)[0]);
    vResult = _mm_shuffle_ps(vResult,vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,reinterpret_cast<const __m128i *>(&vResult)[0]);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreXDec4
(
    XMXDEC4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Min = {-511.0f, -511.0f, -511.0f, 0.0f};
    static const XMVECTOR  Max = {511.0f, 511.0f, 511.0f, 3.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min, Max);

    pDestination->v = ((uint32_t)N.vector4_f32[3] << 30) |
                       (((int32_t)N.vector4_f32[2] & 0x3FF) << 20) |
                       (((int32_t)N.vector4_f32[1] & 0x3FF) << 10) |
                       (((int32_t)N.vector4_f32[0] & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 MinXDec4 = {-511.0f,-511.0f,-511.0f, 0.0f};
    static const XMVECTORF32 MaxXDec4 = { 511.0f, 511.0f, 511.0f, 3.0f};
    static const XMVECTORF32 ScaleXDec4 = {1.0f,1024.0f/2.0f,1024.0f*1024.0f,1024.0f*1024.0f*1024.0f/2.0f};
    static const XMVECTORI32 MaskXDec4= {0x3FF,0x3FF<<(10-1),0x3FF<<20,0x3<<(30-1)};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,MinXDec4);
    vResult = _mm_min_ps(vResult,MaxXDec4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleXDec4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskXDec4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // Perform a single bit left shift on y|w
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUDecN4
(
    XMUDECN4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {1023.0f, 1023.0f, 1023.0f, 3.0f};

    assert(pDestination);

    N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);

    pDestination->v = ((uint32_t)N.vector4_f32[3] << 30) |
                       (((uint32_t)N.vector4_f32[2] & 0x3FF) << 20) |
                       (((uint32_t)N.vector4_f32[1] & 0x3FF) << 10) |
                       (((uint32_t)N.vector4_f32[0] & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 ScaleUDecN4 = {1023.0f,1023.0f*1024.0f*0.5f,1023.0f*1024.0f*1024.0f,3.0f*1024.0f*1024.0f*1024.0f*0.5f};
    static const XMVECTORI32 MaskUDecN4= {0x3FF,0x3FF<<(10-1),0x3FF<<20,0x3<<(30-1)};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUDecN4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUDecN4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // Perform a left shift by one bit on y|w
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUDec4
(
    XMUDEC4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Max = {1023.0f, 1023.0f, 1023.0f, 3.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max);

    pDestination->v = ((uint32_t)N.vector4_f32[3] << 30) |
                       (((uint32_t)N.vector4_f32[2] & 0x3FF) << 20) |
                       (((uint32_t)N.vector4_f32[1] & 0x3FF) << 10) |
                       (((uint32_t)N.vector4_f32[0] & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 MaxUDec4 = { 1023.0f, 1023.0f, 1023.0f, 3.0f};
    static const XMVECTORF32 ScaleUDec4 = {1.0f,1024.0f/2.0f,1024.0f*1024.0f,1024.0f*1024.0f*1024.0f/2.0f};
    static const XMVECTORI32 MaskUDec4= {0x3FF,0x3FF<<(10-1),0x3FF<<20,0x3<<(30-1)};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,MaxUDec4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUDec4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUDec4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // Perform a left shift by one bit on y|w
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreDecN4
(
    XMDECN4* pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {511.0f, 511.0f, 511.0f, 1.0f};

    assert(pDestination);

    N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);

    pDestination->v = ((int32_t)N.vector4_f32[3] << 30) |
                       (((int32_t)N.vector4_f32[2] & 0x3FF) << 20) |
                       (((int32_t)N.vector4_f32[1] & 0x3FF) << 10) |
                       (((int32_t)N.vector4_f32[0] & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 ScaleDecN4 = {511.0f,511.0f*1024.0f,511.0f*1024.0f*1024.0f,1.0f*1024.0f*1024.0f*1024.0f};
    static const XMVECTORI32 MaskDecN4= {0x3FF,0x3FF<<10,0x3FF<<20,0x3<<30};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleDecN4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskDecN4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreDec4
(
    XMDEC4*  pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Min = {-511.0f, -511.0f, -511.0f, -1.0f};
    static const XMVECTOR  Max = {511.0f, 511.0f, 511.0f, 1.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min, Max);

    pDestination->v = ((int32_t)N.vector4_f32[3] << 30) |
                       (((int32_t)N.vector4_f32[2] & 0x3FF) << 20) |
                       (((int32_t)N.vector4_f32[1] & 0x3FF) << 10) |
                       (((int32_t)N.vector4_f32[0] & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 MinDec4 = {-511.0f,-511.0f,-511.0f,-1.0f};
    static const XMVECTORF32 MaxDec4 = { 511.0f, 511.0f, 511.0f, 1.0f};
    static const XMVECTORF32 ScaleDec4 = {1.0f,1024.0f,1024.0f*1024.0f,1024.0f*1024.0f*1024.0f};
    static const XMVECTORI32 MaskDec4= {0x3FF,0x3FF<<10,0x3FF<<20,0x3<<30};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,MinDec4);
    vResult = _mm_min_ps(vResult,MaxDec4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleDec4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskDec4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUByteN4
(
    XMUBYTEN4* pDestination, 
    FXMVECTOR V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {255.0f, 255.0f, 255.0f, 255.0f};

    assert(pDestination);

    N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    pDestination->x = (uint8_t)N.vector4_f32[0];
    pDestination->y = (uint8_t)N.vector4_f32[1];
    pDestination->z = (uint8_t)N.vector4_f32[2];
    pDestination->w = (uint8_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 ScaleUByteN4 = {255.0f,255.0f*256.0f*0.5f,255.0f*256.0f*256.0f,255.0f*256.0f*256.0f*256.0f*0.5f};
    static const XMVECTORI32 MaskUByteN4 = {0xFF,0xFF<<(8-1),0xFF<<16,0xFF<<(24-1)};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUByteN4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUByteN4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // Perform a single bit left shift to fix y|w 
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUByte4
(
    XMUBYTE4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Max = {255.0f, 255.0f, 255.0f, 255.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    pDestination->x = (uint8_t)N.vector4_f32[0];
    pDestination->y = (uint8_t)N.vector4_f32[1];
    pDestination->z = (uint8_t)N.vector4_f32[2];
    pDestination->w = (uint8_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 MaxUByte4 = { 255.0f, 255.0f, 255.0f, 255.0f};
    static const XMVECTORF32 ScaleUByte4 = {1.0f,256.0f*0.5f,256.0f*256.0f,256.0f*256.0f*256.0f*0.5f};
    static const XMVECTORI32 MaskUByte4 = {0xFF,0xFF<<(8-1),0xFF<<16,0xFF<<(24-1)};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,MaxUByte4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUByte4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUByte4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // Perform a single bit left shift to fix y|w 
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreByteN4
(
    XMBYTEN4* pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTORF32  Scale = {127.0f, 127.0f, 127.0f, 127.0f};

    assert(pDestination);

    N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(V, Scale.v);
    N = XMVectorRound(N);

    pDestination->x = (int8_t)N.vector4_f32[0];
    pDestination->y = (int8_t)N.vector4_f32[1];
    pDestination->z = (int8_t)N.vector4_f32[2];
    pDestination->w = (int8_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 ScaleByteN4 = {127.0f,127.0f*256.0f,127.0f*256.0f*256.0f,127.0f*256.0f*256.0f*256.0f};
    static const XMVECTORI32 MaskByteN4 = {0xFF,0xFF<<8,0xFF<<16,0xFF<<24};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleByteN4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskByteN4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreByte4
(
    XMBYTE4*  pDestination, 
    FXMVECTOR  V
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR               N;
    static const XMVECTOR  Min = {-127.0f, -127.0f, -127.0f, -127.0f};
    static const XMVECTOR  Max = {127.0f, 127.0f, 127.0f, 127.0f};

    assert(pDestination);

    N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    pDestination->x = (int8_t)N.vector4_f32[0];
    pDestination->y = (int8_t)N.vector4_f32[1];
    pDestination->z = (int8_t)N.vector4_f32[2];
    pDestination->w = (int8_t)N.vector4_f32[3];

#elif defined(_XM_SSE_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32 MinByte4 = {-127.0f,-127.0f,-127.0f,-127.0f};
    static const XMVECTORF32 MaxByte4 = { 127.0f, 127.0f, 127.0f, 127.0f};
    static const XMVECTORF32 ScaleByte4 = {1.0f,256.0f,256.0f*256.0f,256.0f*256.0f*256.0f};
    static const XMVECTORI32 MaskByte4 = {0xFF,0xFF<<8,0xFF<<16,0xFF<<24};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,MinByte4);
    vResult = _mm_min_ps(vResult,MaxByte4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleByte4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskByte4);
    // Do a horizontal or of 4 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(3,2,3,2));
    // x = x|z, y = y|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(1,1,1,1));
    // i = x|y|z|w
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),reinterpret_cast<const __m128 *>(&vResulti)[0]);
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreUNibble4
(
     XMUNIBBLE4* pDestination,
     FXMVECTOR V
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32  Max = {15.0f,15.0f,15.0f,15.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // No SSE operations will write to 16-bit values, so we have to extract them manually
    uint16_t x = static_cast<uint16_t>(_mm_extract_epi16(vInt,0));
    uint16_t y = static_cast<uint16_t>(_mm_extract_epi16(vInt,2));
    uint16_t z = static_cast<uint16_t>(_mm_extract_epi16(vInt,4));
    uint16_t w = static_cast<uint16_t>(_mm_extract_epi16(vInt,6));
    pDestination->v = ((w & 0xF) << 12) |
                      ((z & 0xF) << 8) |
                      ((y & 0xF) << 4) |
                      ((x & 0xF));
#else
    XMVECTOR               N;
    static const XMVECTORF32  Max = {15.0f,15.0f,15.0f,15.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max.v);
    N = XMVectorRound(N);

    pDestination->v = (((uint16_t)N.vector4_f32[3] & 0xF) << 12) |
                      (((uint16_t)N.vector4_f32[2] & 0xF) << 8) |
                      (((uint16_t)N.vector4_f32[1] & 0xF) << 4) |
                      (((uint16_t)N.vector4_f32[0] & 0xF));
#endif !_XM_SSE_INTRINSICS_
}

//------------------------------------------------------------------------------

inline void PackedVector::XMStoreU555
(
     XMU555* pDestination,
     FXMVECTOR V
)
{
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    assert(pDestination);
    static const XMVECTORF32  Max = {31.0f, 31.0f, 31.0f, 1.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // No SSE operations will write to 16-bit values, so we have to extract them manually
    uint16_t x = static_cast<uint16_t>(_mm_extract_epi16(vInt,0));
    uint16_t y = static_cast<uint16_t>(_mm_extract_epi16(vInt,2));
    uint16_t z = static_cast<uint16_t>(_mm_extract_epi16(vInt,4));
    uint16_t w = static_cast<uint16_t>(_mm_extract_epi16(vInt,6));
    pDestination->v = ((w) ? 0x8000 : 0) |
                      ((z & 0x1F) << 10) |
                      ((y & 0x1F) << 5) |
                      ((x & 0x1F));
#else
    XMVECTOR               N;
    static const XMVECTORF32  Max = {31.0f, 31.0f, 31.0f, 1.0f};

    assert(pDestination);

    N = XMVectorClamp(V, XMVectorZero(), Max.v);
    N = XMVectorRound(N);

    pDestination->v = ((N.vector4_f32[3] > 0.f) ? 0x8000 : 0) |
                      (((uint16_t)N.vector4_f32[2] & 0x1F) << 10) |
                      (((uint16_t)N.vector4_f32[1] & 0x1F) << 5) |
                      (((uint16_t)N.vector4_f32[0] & 0x1F));
#endif !_XM_SSE_INTRINSICS_
}


/****************************************************************************
 *
 * XMCOLOR operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMCOLOR::XMCOLOR
(
    float _r,
    float _g,
    float _b,
    float _a
)
{
    XMStoreColor(this, XMVectorSet(_r, _g, _b, _a));
}

//------------------------------------------------------------------------------

inline PackedVector::XMCOLOR::XMCOLOR
(
    const float* pArray
)
{
    XMStoreColor(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMHALF2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMHALF2::XMHALF2
(
    float _x,
    float _y
)
{
    x = XMConvertFloatToHalf(_x);
    y = XMConvertFloatToHalf(_y);
}

//------------------------------------------------------------------------------

inline PackedVector::XMHALF2::XMHALF2
(
    const float* pArray
)
{
    assert( pArray != nullptr );
    x = XMConvertFloatToHalf(pArray[0]);
    y = XMConvertFloatToHalf(pArray[1]);
}

/****************************************************************************
 *
 * XMSHORTN2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMSHORTN2::XMSHORTN2
(
    float _x,
    float _y
)
{
    XMStoreShortN2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMSHORTN2::XMSHORTN2
(
    const float* pArray
)
{
    XMStoreShortN2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMSHORT2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMSHORT2::XMSHORT2
(
    float _x,
    float _y
)
{
    XMStoreShort2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMSHORT2::XMSHORT2
(
    const float* pArray
)
{
    XMStoreShort2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMUSHORTN2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORTN2::XMUSHORTN2
(
    float _x,
    float _y
)
{
    XMStoreUShortN2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORTN2::XMUSHORTN2
(
    const float* pArray
)
{
    XMStoreUShortN2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMUSHORT2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORT2::XMUSHORT2
(
    float _x,
    float _y
)
{
    XMStoreUShort2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORT2::XMUSHORT2
(
    const float* pArray
)
{
    XMStoreUShort2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMBYTEN2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMBYTEN2::XMBYTEN2
(
    float _x,
    float _y
)
{
    XMStoreByteN2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMBYTEN2::XMBYTEN2
(
    const float* pArray
)
{
    XMStoreByteN2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMBYTE2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMBYTE2::XMBYTE2
(
    float _x,
    float _y
)
{
    XMStoreByte2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMBYTE2::XMBYTE2
(
    const float* pArray
)
{
    XMStoreByte2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMUBYTEN2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTEN2::XMUBYTEN2
(
    float _x,
    float _y
)
{
    XMStoreUByteN2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTEN2::XMUBYTEN2
(
    const float* pArray
)
{
    XMStoreUByteN2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMUBYTE2 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTE2::XMUBYTE2
(
    float _x,
    float _y
)
{
    XMStoreUByte2(this, XMVectorSet(_x, _y, 0.0f, 0.0f));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTE2::XMUBYTE2
(
    const float* pArray
)
{
    XMStoreUByte2(this, XMLoadFloat2((const XMFLOAT2*)pArray));
}

/****************************************************************************
 *
 * XMU565 operators
 *
 ****************************************************************************/

inline PackedVector::XMU565::XMU565
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreU565(this, XMVectorSet( _x, _y, _z, 0.0f ));
}

inline PackedVector::XMU565::XMU565
(
    const float *pArray
)
{
    XMStoreU565(this, XMLoadFloat3((const XMFLOAT3*)pArray ));
}

/****************************************************************************
 *
 * XMFLOAT3PK operators
 *
 ****************************************************************************/

inline PackedVector::XMFLOAT3PK::XMFLOAT3PK
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreFloat3PK(this, XMVectorSet( _x, _y, _z, 0.0f ));
}

inline PackedVector::XMFLOAT3PK::XMFLOAT3PK
(
    const float *pArray
)
{
    XMStoreFloat3PK(this, XMLoadFloat3((const XMFLOAT3*)pArray ));
}

/****************************************************************************
 *
 * XMFLOAT3SE operators
 *
 ****************************************************************************/

inline PackedVector::XMFLOAT3SE::XMFLOAT3SE
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreFloat3SE(this, XMVectorSet( _x, _y, _z, 0.0f ));
}

inline PackedVector::XMFLOAT3SE::XMFLOAT3SE
(
    const float *pArray
)
{
    XMStoreFloat3SE(this, XMLoadFloat3((const XMFLOAT3*)pArray ));
}

/****************************************************************************
 *
 * XMHALF4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMHALF4::XMHALF4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    x = XMConvertFloatToHalf(_x);
    y = XMConvertFloatToHalf(_y);
    z = XMConvertFloatToHalf(_z);
    w = XMConvertFloatToHalf(_w);
}

//------------------------------------------------------------------------------

inline PackedVector::XMHALF4::XMHALF4
(
    const float* pArray
)
{
    XMConvertFloatToHalfStream(&x, sizeof(HALF), pArray, sizeof(float), 4);
}

/****************************************************************************
 *
 * XMSHORTN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMSHORTN4::XMSHORTN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreShortN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMSHORTN4::XMSHORTN4
(
    const float* pArray
)
{
    XMStoreShortN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMSHORT4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMSHORT4::XMSHORT4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreShort4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMSHORT4::XMSHORT4
(
    const float* pArray
)
{
    XMStoreShort4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUSHORTN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORTN4::XMUSHORTN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUShortN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORTN4::XMUSHORTN4
(
    const float* pArray
)
{
    XMStoreUShortN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUSHORT4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORT4::XMUSHORT4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUShort4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUSHORT4::XMUSHORT4
(
    const float* pArray
)
{
    XMStoreUShort4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMXDECN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMXDECN4::XMXDECN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreXDecN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMXDECN4::XMXDECN4
(
    const float* pArray
)
{
    XMStoreXDecN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMXDEC4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMXDEC4::XMXDEC4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreXDec4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMXDEC4::XMXDEC4
(
    const float* pArray
)
{
    XMStoreXDec4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMDECN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMDECN4::XMDECN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreDecN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMDECN4::XMDECN4
(
    const float* pArray
)
{
    XMStoreDecN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMDEC4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMDEC4::XMDEC4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreDec4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMDEC4::XMDEC4
(
    const float* pArray
)
{
    XMStoreDec4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUDECN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUDECN4::XMUDECN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUDecN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUDECN4::XMUDECN4
(
    const float* pArray
)
{
    XMStoreUDecN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUDEC4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUDEC4::XMUDEC4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUDec4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUDEC4::XMUDEC4
(
    const float* pArray
)
{
    XMStoreUDec4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMBYTEN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMBYTEN4::XMBYTEN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreByteN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMBYTEN4::XMBYTEN4
(
    const float* pArray
)
{
    XMStoreByteN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMBYTE4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMBYTE4::XMBYTE4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreByte4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMBYTE4::XMBYTE4
(
    const float* pArray
)
{
    XMStoreByte4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUBYTEN4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTEN4::XMUBYTEN4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUByteN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTEN4::XMUBYTEN4
(
    const float* pArray
)
{
    XMStoreUByteN4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUBYTE4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTE4::XMUBYTE4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUByte4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUBYTE4::XMUBYTE4
(
    const float* pArray
)
{
    XMStoreUByte4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMUNIBBLE4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMUNIBBLE4::XMUNIBBLE4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUNibble4(this, XMVectorSet( _x, _y, _z, _w ));
}

//------------------------------------------------------------------------------

inline PackedVector::XMUNIBBLE4::XMUNIBBLE4
(
    const float *pArray
)
{
    XMStoreUNibble4(this, XMLoadFloat4((const XMFLOAT4*)pArray));
}

/****************************************************************************
 *
 * XMU555 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::XMU555::XMU555
(
    float _x,
    float _y,
    float _z,
    bool _w
)
{
    XMStoreU555(this, XMVectorSet(_x, _y, _z, ((_w) ? 1.0f : 0.0f) ));
}

//------------------------------------------------------------------------------

inline PackedVector::XMU555::XMU555
(
    const float *pArray,
    bool _w
)
{
    XMVECTOR V = XMLoadFloat3((const XMFLOAT3*)pArray);
    XMStoreU555(this, XMVectorSetW(V, ((_w) ? 1.0f : 0.0f) ));
}


#undef XM3_PACK_FACTOR
#undef XM3_UNPACK_FACTOR_UNSIGNED
#undef XM3_UNPACK_FACTOR_SIGNED
#undef XM3_UNPACK_UNSIGNEDN_OFFSET
#undef XM3_UNPACK_UNSIGNEDN_SCALE
#undef XM3_UNPACK_SIGNEDN_SCALE
#undef XM3_PACK_UNSIGNEDN_SCALE
#undef XM3_PACK_SIGNEDN_SCALE
#undef XM3_PACK_OFFSET

