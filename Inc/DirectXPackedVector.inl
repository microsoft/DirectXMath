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
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    uint32_t Mantissa = (uint32_t)(Value & 0x03FF);

    uint32_t Exponent = (Value & 0x7C00);
    if ( Exponent == 0x7C00 ) // INF/NAN
    {
        Exponent = (uint32_t)143;
    }
    else if (Exponent != 0)  // The value is normalized
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

    uint32_t Result = ((Value & 0x8000) << 16) | // Sign
                      ((Exponent + 112) << 23) | // Exponent
                      (Mantissa << 13);          // Mantissa

    return reinterpret_cast<float*>(&Result)[0];
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline float* PackedVector::XMConvertHalfToFloatStream
(
    float*      pOutputStream, 
    size_t      OutputStride, 
    const HALF* pInputStream, 
    size_t      InputStride, 
    size_t      HalfCount
)
{
    assert(pOutputStream);
    assert(pInputStream);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)

    const uint8_t* pHalf = reinterpret_cast<const uint8_t*>(pInputStream);
    uint8_t* pFloat = reinterpret_cast<uint8_t*>(pOutputStream);

    for (size_t i = 0; i < HalfCount; i++)
    {
        *reinterpret_cast<float*>(pFloat) = XMConvertHalfToFloat(reinterpret_cast<const HALF*>(pHalf)[0]);
        pHalf += InputStride;
        pFloat += OutputStride; 
    }

    return pOutputStream;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------

inline PackedVector::HALF PackedVector::XMConvertFloatToHalf
(
    float Value
)
{
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    uint32_t Result;

    uint32_t IValue = reinterpret_cast<uint32_t *>(&Value)[0];
    uint32_t Sign = (IValue & 0x80000000U) >> 16U;
    IValue = IValue & 0x7FFFFFFFU;      // Hack off the sign

    if (IValue > 0x477FE000U)
    {
        // The number is too large to be represented as a half.  Saturate to infinity.
        if (((IValue & 0x7F800000) == 0x7F800000) && ((IValue & 0x7FFFFF ) != 0))
        {
            Result = 0x7FFF; // NAN
        }
        else
        {
            Result = 0x7C00U; // INF
        }
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
_Use_decl_annotations_
inline PackedVector::HALF* PackedVector::XMConvertFloatToHalfStream
(
    HALF* pOutputStream, 
    size_t       OutputStride, 
    const float* pInputStream, 
    size_t       InputStride, 
    size_t       FloatCount
)
{
    assert(pOutputStream);
    assert(pInputStream);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_) || defined(XM_NO_MISALIGNED_VECTOR_ACCESS)

    const uint8_t* pFloat = reinterpret_cast<const uint8_t*>(pInputStream);
    uint8_t* pHalf = reinterpret_cast<uint8_t*>(pOutputStream);

    for (size_t i = 0; i < FloatCount; i++)
    {
        *reinterpret_cast<HALF*>(pHalf) = XMConvertFloatToHalf(reinterpret_cast<const float*>(pFloat)[0]);
        pFloat += InputStride; 
        pHalf += OutputStride;
    }
    return pOutputStream;

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

/****************************************************************************
 *
 * Vector and matrix load operations
 *
 ****************************************************************************/
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadColor
(
    const XMCOLOR* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    // int32_t -> Float conversions are done in one instruction.
    // uint32_t -> Float calls a runtime function. Keep in int32_t
    int32_t iColor = (int32_t)(pSource->c);
    XMVECTORF32 vColor = {
        (float)((iColor >> 16) & 0xFF) * (1.0f/255.0f),
        (float)((iColor >> 8) & 0xFF) * (1.0f/255.0f),
        (float)(iColor & 0xFF) * (1.0f/255.0f),
        (float)((iColor >> 24) & 0xFF) * (1.0f/255.0f)
    };
    return vColor.v;
#elif defined(_XM_SSE_INTRINSICS_)
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
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadHalf2
(
    const XMHALF2* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        XMConvertHalfToFloat(pSource->x),
        XMConvertHalfToFloat(pSource->y),
        0.0f,
        0.0f
    };
    return vResult.v;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadShortN2
(
    const XMSHORTN2* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (pSource->x == -32768) ? -1.f : ((float)pSource->x * (1.0f/32767.0f)),
        (pSource->y == -32768) ? -1.f : ((float)pSource->y * (1.0f/32767.0f)),
        0.0f,
        0.0f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // x needs to be sign extended
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
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
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadShort2
(
    const XMSHORT2* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        0.f,
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // x needs to be sign extended
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x - 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16);
    // Y is 65536 too large
    return _mm_mul_ps(vTemp,g_XMFixupY16);
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUShortN2
(
    const XMUSHORTN2* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x / 65535.0f,
        (float)pSource->y / 65535.0f,
        0.f,
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 FixupY16 = {1.0f/65535.0f,1.0f/(65535.0f*65536.0f),0.0f,0.0f};
    static const XMVECTORF32 FixaddY16 = {0,32768.0f*65536.0f,0,0};
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // y needs to be sign flipped
    vTemp = _mm_xor_ps(vTemp,g_XMFlipY);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // y + 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,FixaddY16);
    // Y is 65536 times too large
    vTemp = _mm_mul_ps(vTemp,FixupY16);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUShort2
(
    const XMUSHORT2* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        0.f,
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 FixaddY16 = {0,32768.0f,0,0};
    // Splat the two shorts in all four entries (WORD alignment okay,
    // DWORD alignment preferred)
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
    vTemp = _mm_and_ps(vTemp,g_XMMaskX16Y16);
    // y needs to be sign flipped
    vTemp = _mm_xor_ps(vTemp,g_XMFlipY);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // Y is 65536 times too large
    vTemp = _mm_mul_ps(vTemp,g_XMFixupY16);
    // y + 0x8000 to undo the signed order.
    vTemp = _mm_add_ps(vTemp,FixaddY16);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadByteN2
(
    const XMBYTEN2* pSource
)
{
    assert(pSource);
    XMVECTORF32 vResult = {
        (pSource->x == -128) ? -1.f : ((float)pSource->x * (1.0f/127.0f)),
        (pSource->y == -128) ? -1.f : ((float)pSource->y * (1.0f/127.0f)),
        0.0f,
        0.0f
    };
    return vResult.v;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadByte2
(
    const XMBYTE2* pSource
)
{
    assert(pSource);
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        0.0f,
        0.0f
    };
    return vResult.v;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUByteN2
(
    const XMUBYTEN2* pSource
)
{
    assert(pSource);
    XMVECTORF32 vResult = {
        (float)pSource->x * (1.0f/255.0f),
        (float)pSource->y * (1.0f/255.0f),
        0.0f,
        0.0f
    };
    return vResult.v;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUByte2
(
    const XMUBYTE2* pSource
)
{
    assert(pSource);
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        0.0f,
        0.0f
    };
    return vResult.v;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadU565
(
    const XMU565* pSource
)
{
    assert(pSource);
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORI32 U565And = {0x1F,0x3F<<5,0x1F<<11,0};
    static const XMVECTORF32 U565Mul = {1.0f,1.0f/32.0f,1.0f/2048.f,0};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,U565And);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Normalize x, y, and z
    vResult = _mm_mul_ps(vResult,U565Mul);
    return vResult;
#else
    XMVECTORF32 vResult = {
        float(pSource->v & 0x1F),
        float((pSource->v >> 5) & 0x3F),
        float((pSource->v >> 11) & 0x1F),
        0.f,
    };
    return vResult.v;
#endif // !_XM_SSE_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadFloat3PK
(
    const XMFLOAT3PK* pSource
)
{
    assert(pSource);

    __declspec(align(16)) uint32_t Result[4];
    uint32_t Mantissa;
    uint32_t Exponent;

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

    return XMLoadFloat3A( reinterpret_cast<const XMFLOAT3A*>(&Result) );
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadFloat3SE
(
    const XMFLOAT3SE* pSource
)
{
    assert(pSource);

    __declspec(align(16)) uint32_t Result[4];
    uint32_t Mantissa;
    uint32_t Exponent, ExpBits;

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

    return XMLoadFloat3A( reinterpret_cast<const XMFLOAT3A*>(&Result) );
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadHalf4
(
    const XMHALF4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_) 
    XMVECTORF32 vResult = {
        XMConvertHalfToFloat(pSource->x),
        XMConvertHalfToFloat(pSource->y),
        XMConvertHalfToFloat(pSource->z),
        XMConvertHalfToFloat(pSource->w)
    };
    return vResult.v;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadShortN4
(
    const XMSHORTN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = {
        (pSource->x == -32768) ? -1.f : ((float)pSource->x * (1.0f/32767.0f)),
        (pSource->y == -32768) ? -1.f : ((float)pSource->y * (1.0f/32767.0f)),
        (pSource->z == -32768) ? -1.f : ((float)pSource->z * (1.0f/32767.0f)),
        (pSource->w == -32768) ? -1.f : ((float)pSource->w * (1.0f/32767.0f))
    };
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vInt = vld1_s16( (const int16_t*)pSource );
    __n128 V = vmovl_s16( vInt );
    V = vcvtq_f32_s32( V );
    const __n128 Scale = vdupq_n_f32( 1.0f/32767.0f );
    V = vmulq_f32( V, Scale );
    return vmaxq_f32( V, g_XMNegativeOne );
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd),g_XMMaskX16Y16Z16W16);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16Z16W16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16Z16W16);
    // Convert to -1.0f - 1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMNormalizeX16Y16Z16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(3,1,2,0));
    // Clamp result (for case of -32768)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadShort4
(
    const XMSHORT4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        (float)pSource->z,
        (float)pSource->w
    };
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vInt = vld1_s16( (const int16_t*)pSource );
    __n128 V = vmovl_s16( vInt );
    return vcvtq_f32_s32( V );
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd),g_XMMaskX16Y16Z16W16);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipX16Y16Z16W16);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMFixX16Y16Z16W16);
    // Fix y and w because they are 65536 too large
    vTemp = _mm_mul_ps(vTemp,g_XMFixupY16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    return XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(3,1,2,0));
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUShortN4
(
    const XMUSHORTN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x / 65535.0f,
        (float)pSource->y / 65535.0f,
        (float)pSource->z / 65535.0f,
        (float)pSource->w / 65535.0f
    };
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vInt = vld1_u16( (const uint16_t*)pSource );
    __n128 V = vmovl_u16( vInt );
    V = vcvtq_f32_u32( V );
    const __n128 Scale = vdupq_n_f32( 1.0f/65535.0f );
    return vmulq_f32( V, Scale );
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 FixupY16W16 = {1.0f/65535.0f,1.0f/65535.0f,1.0f/(65535.0f*65536.0f),1.0f/(65535.0f*65536.0f)};
    static const XMVECTORF32 FixaddY16W16  = {0,0,32768.0f*65536.0f,32768.0f*65536.0f};
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd),g_XMMaskX16Y16Z16W16);
    // y and w are signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipZW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // y and w + 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,FixaddY16W16);
    // Fix y and w because they are 65536 too large
    vTemp = _mm_mul_ps(vTemp,FixupY16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    return XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(3,1,2,0));
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUShort4
(
    const XMUSHORT4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        (float)pSource->z,
        (float)pSource->w
    };
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n64 vInt = vld1_u16( (const uint16_t*)pSource );
    __n128 V = vmovl_u16( vInt );
    return vcvtq_f32_u32( V );
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 FixaddY16W16  = {0,0,32768.0f,32768.0f};
    // Splat the color in all four entries (x,z,y,w)
    __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double *>(&pSource->x));
    // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
    __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd),g_XMMaskX16Y16Z16W16);
    // y and w are signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipZW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // Fix y and w because they are 65536 too large
    vTemp = _mm_mul_ps(vTemp,g_XMFixupY16W16);
    // y and w + 0x8000 to complete the conversion
    vTemp = _mm_add_ps(vTemp,FixaddY16W16);
    // Very important! The entries are x,z,y,w, flip it to x,y,z,w
    return XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(3,1,2,0));
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadXDecN4
(
    const XMXDECN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
    uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

    XMVECTORF32 vResult = {
        (ElementX == 0x200) ? -1.f : ((float)(int16_t)(ElementX | SignExtend[ElementX >> 9]) / 511.0f),
        (ElementY == 0x200) ? -1.f : ((float)(int16_t)(ElementY | SignExtend[ElementY >> 9]) / 511.0f),
        (ElementZ == 0x200) ? -1.f : ((float)(int16_t)(ElementZ | SignExtend[ElementZ >> 9]) / 511.0f),
        (float)(pSource->v >> 30) / 3.0f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the color in all four entries
    __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskA2B10G10R10);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipA2B10G10R10);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMFixAA2B10G10R10);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMNormalizeA2B10G10R10);
    // Clamp result (for case of -512)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadXDec4
(
    const XMXDEC4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
    uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

    XMVECTORF32 vResult = {
        (float)(int16_t)(ElementX | SignExtend[ElementX >> 9]),
        (float)(int16_t)(ElementY | SignExtend[ElementY >> 9]),
        (float)(int16_t)(ElementZ | SignExtend[ElementZ >> 9]),
        (float)(pSource->v >> 30)
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORI32 XDec4Xor = {0x200, 0x200<<10, 0x200<<20, 0x80000000};
    static const XMVECTORF32 XDec4Add = {-512.0f,-512.0f*1024.0f,-512.0f*1024.0f*1024.0f,32768*65536.0f};
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,XDec4Xor);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,XDec4Add);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMMulDec4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUDecN4
(
    const XMUDECN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
    uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

    XMVECTORF32 vResult = {
        (float)ElementX / 1023.0f,
        (float)ElementY / 1023.0f,
        (float)ElementZ / 1023.0f,
        (float)(pSource->v >> 30) / 3.0f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 UDecN4Mul = {1.0f/1023.0f,1.0f/(1023.0f*1024.0f),1.0f/(1023.0f*1024.0f*1024.0f),1.0f/(3.0f*1024.0f*1024.0f*1024.0f)};
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,UDecN4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUDec4
(
    const XMUDEC4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
    uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

    XMVECTORF32 vResult = {
        (float)ElementX,
        (float)ElementY,
        (float)ElementZ,
        (float)(pSource->v >> 30)
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMMulDec4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadDecN4
(
    const XMDECN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};
    static const uint32_t SignExtendW[] = {0x00000000, 0xFFFFFFFC};

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
    uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;
    uint32_t ElementW = pSource->v >> 30;

    XMVECTORF32 vResult = {
        (ElementX == 0x200) ? -1.f : ((float)(int16_t)(ElementX | SignExtend[ElementX >> 9]) / 511.0f),
        (ElementY == 0x200) ? -1.f : ((float)(int16_t)(ElementY | SignExtend[ElementY >> 9]) / 511.0f),
        (ElementZ == 0x200) ? -1.f : ((float)(int16_t)(ElementZ | SignExtend[ElementZ >> 9]) / 511.0f),
        (ElementW == 0x2)   ? -1.f : ((float)(int16_t)(ElementW | SignExtendW[(ElementW >> 1) & 1]))
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 DecN4Mul = {1.0f/511.0f,1.0f/(511.0f*1024.0f),1.0f/(511.0f*1024.0f*1024.0f),1.0f/(1024.0f*1024.0f*1024.0f)};
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorDec4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,DecN4Mul);
    // Clamp result (for case of -512/-1)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadDec4
(
    const XMDEC4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFFFFC00};
    static const uint32_t SignExtendW[] = {0x00000000, 0xFFFFFFFC};

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
    uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;
    uint32_t ElementW = pSource->v >> 30;

    XMVECTORF32 vResult = {
        (float)(int16_t)(ElementX | SignExtend[ElementX >> 9]),
        (float)(int16_t)(ElementY | SignExtend[ElementY >> 9]),
        (float)(int16_t)(ElementZ | SignExtend[ElementZ >> 9]),
        (float)(int16_t)(ElementW | SignExtendW[ElementW >> 1])
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Splat the color in all four entries
    XMVECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskDec4);
    // a is unsigned! Flip the bit to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorDec4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // RGB + 0, A + 0x80000000.f to undo the signed order.
    vTemp = _mm_add_ps(vTemp,g_XMAddDec4);
    // Convert 0-255 to 0.0f-1.0f
    vTemp = _mm_mul_ps(vTemp,g_XMMulDec4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUByteN4
(
    const XMUBYTEN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x / 255.0f,
        (float)pSource->y / 255.0f,
        (float)pSource->z / 255.0f,
        (float)pSource->w / 255.0f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadUByteN4Mul = {1.0f/255.0f,1.0f/(255.0f*256.0f),1.0f/(255.0f*65536.0f),1.0f/(255.0f*65536.0f*256.0f)};
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // w is signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // w + 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadUByteN4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUByte4
(
    const XMUBYTE4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        (float)pSource->z,
        (float)pSource->w
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadUByte4Mul = {1.0f,1.0f/256.0f,1.0f/65536.0f,1.0f/(65536.0f*256.0f)};
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // w is signed! Flip the bits to convert the order to unsigned
    vTemp = _mm_xor_ps(vTemp,g_XMFlipW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // w + 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddUDec4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadUByte4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadByteN4
(
    const XMBYTEN4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (pSource->x == -128) ? -1.f : ((float)pSource->x / 127.0f),
        (pSource->y == -128) ? -1.f : ((float)pSource->y / 127.0f),
        (pSource->z == -128) ? -1.f : ((float)pSource->z / 127.0f),
        (pSource->w == -128) ? -1.f : ((float)pSource->w / 127.0f)
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadByteN4Mul = {1.0f/127.0f,1.0f/(127.0f*256.0f),1.0f/(127.0f*65536.0f),1.0f/(127.0f*65536.0f*256.0f)};
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // x,y and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorByte4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
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
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadByte4
(
    const XMBYTE4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)pSource->x,
        (float)pSource->y,
        (float)pSource->z,
        (float)pSource->w
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadByte4Mul = {1.0f,1.0f/256.0f,1.0f/65536.0f,1.0f/(65536.0f*256.0f)};
    // Splat the color in all four entries (x,z,y,w)
    XMVECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float *>(&pSource->x));
    // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
    vTemp = _mm_and_ps(vTemp,g_XMMaskByte4);
    // x,y and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorByte4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x, y and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddByte4);
    // Fix y, z and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadByte4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadUNibble4
(
     const XMUNIBBLE4* pSource
)
{
    assert(pSource);
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORI32 UNibble4And = {0xF,0xF0,0xF00,0xF000};
    static const XMVECTORF32 UNibble4Mul = {1.0f,1.0f/16.f,1.0f/256.f,1.0f/4096.f};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,UNibble4And);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Normalize x, y, and z
    vResult = _mm_mul_ps(vResult,UNibble4Mul);
    return vResult;
#else
    XMVECTORF32 vResult = {
        float(pSource->v & 0xF),
        float((pSource->v >> 4) & 0xF),
        float((pSource->v >> 8) & 0xF),
        float((pSource->v >> 12) & 0xF)
    };
    return vResult.v;
#endif // !_XM_SSE_INTRISICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::XMLoadU555
(
     const XMU555* pSource
)
{
    assert(pSource);
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORI32 U555And = {0x1F,0x1F<<5,0x1F<<10,0x8000};
    static const XMVECTORF32 U555Mul = {1.0f,1.0f/32.f,1.0f/1024.f,1.0f/32768.f};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,U555And);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Normalize x, y, and z
    vResult = _mm_mul_ps(vResult,U555Mul);
    return vResult;
#else
    XMVECTORF32 vResult = {
        float(pSource->v & 0x1F),
        float((pSource->v >> 5) & 0x1F),
        float((pSource->v >> 10) & 0x1F),
        float((pSource->v >> 15) & 0x1)
    };
    return vResult.v;
#endif // !_XM_SSE_INTRISICS_
}

///begin_xbox360
////////////////////////////////////////////////////////////////////////////////
// PackedVector::Xbox
////////////////////////////////////////////////////////////////////////////////
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadUHenDN3
(
    const XMUHENDN3* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        float(pSource->v & 0x7FF) / 2047.0f,
        float((pSource->v >> 11) & 0x7FF) / 2047.0f,
        float((pSource->v >> 22) & 0x3FF) / 1023.0f,
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 UHenDN3Mul = {1.0f/2047.0f,1.0f/(2047.0f*2048.0f),1.0f/(1023.0f*2048.0f*2048.0f),0};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskHenD3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMFlipZ);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddUHenD3);
    // Normalize x,y and z to -1.0f-1.0f
    vResult = _mm_mul_ps(vResult,UHenDN3Mul);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadUHenD3
(
    const XMUHEND3* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        float(pSource->v & 0x7FF),
        float((pSource->v >> 11) & 0x7FF),
        float((pSource->v >> 22) & 0x3FF),
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskHenD3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMFlipZ);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddUHenD3);
    // Normalize x and y to -1024-1023.0f and z to -512-511.0f
    vResult = _mm_mul_ps(vResult,g_XMMulHenD3);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadHenDN3
(
    const XMHENDN3* pSource
)
{
    assert(pSource);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtendXY[] = {0x00000000, 0xFFFFF800};
    static const uint32_t SignExtendZ[] = {0x00000000, 0xFFFFFC00};

    uint32_t ElementX = pSource->v & 0x7FF;
    uint32_t ElementY = (pSource->v >> 11) & 0x7FF;
    uint32_t ElementZ = (pSource->v >> 22) & 0x3FF;

    XMVECTORF32 vResult = {
        (ElementX == 0x400) ? -1.f : ((float)(int16_t)(ElementX | SignExtendXY[ElementX >> 10]) / 1023.0f),
        (ElementY == 0x400) ? -1.f : ((float)(int16_t)(ElementY | SignExtendXY[ElementY >> 10]) / 1023.0f),
        (ElementZ == 0x200) ? -1.f : ((float)(int16_t)(ElementZ | SignExtendZ[ElementZ >> 9]) / 511.0f),
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 HenDN3Mul = {1.0f/1023.0f,1.0f/(1023.0f*2048.0f),1.0f/(511.0f*2048.0f*2048.0f),0};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskHenD3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMXorHenD3);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddHenD3);
    // Normalize x,y and z to -1.0f-1.0f
    vResult = _mm_mul_ps(vResult,HenDN3Mul);
    // Clamp result (for case of -1024/-512)
    return _mm_max_ps( vResult, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadHenD3
(
    const XMHEND3* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtendXY[] = {0x00000000, 0xFFFFF800};
    static const uint32_t SignExtendZ[] = {0x00000000, 0xFFFFFC00};

    uint32_t ElementX = pSource->v & 0x7FF;
    uint32_t ElementY = (pSource->v >> 11) & 0x7FF;
    uint32_t ElementZ = (pSource->v >> 22) & 0x3FF;

    XMVECTORF32 vResult = {
        (float)(int16_t)(ElementX | SignExtendXY[ElementX >> 10]),
        (float)(int16_t)(ElementY | SignExtendXY[ElementY >> 10]),
        (float)(int16_t)(ElementZ | SignExtendZ[ElementZ >> 9]),
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskHenD3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMXorHenD3);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddHenD3);
    // Normalize x and y to -1024-1023.0f and z to -512-511.0f
    vResult = _mm_mul_ps(vResult,g_XMMulHenD3);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadUDHenN3
(
    const XMUDHENN3* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x7FF;
    uint32_t ElementZ = (pSource->v >> 21) & 0x7FF;

    XMVECTORF32 vResult = {
        (float)ElementX / 1023.0f,
        (float)ElementY / 2047.0f,
        (float)ElementZ / 2047.0f,
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 UDHenN3Mul = {1.0f/1023.0f,1.0f/(2047.0f*1024.0f),1.0f/(2047.0f*1024.0f*2048.0f),0};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskDHen3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMFlipZ);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddUHenD3);
    // Normalize x,y and z to -1.0f-1.0f
    vResult = _mm_mul_ps(vResult,UDHenN3Mul);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadUDHen3
(
    const XMUDHEN3* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x7FF;
    uint32_t ElementZ = (pSource->v >> 21) & 0x7FF;

    XMVECTORF32 vResult = {
        (float)ElementX,
        (float)ElementY,
        (float)ElementZ,
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskDHen3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMFlipZ);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddUHenD3);
    // Normalize x to 0-1023.0f and y and z to 0-2047.0f
    vResult = _mm_mul_ps(vResult,g_XMMulDHen3);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadDHenN3
(
    const XMDHENN3* pSource
)
{
    assert(pSource);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtendX[] = {0x00000000, 0xFFFFFC00};
    static const uint32_t SignExtendYZ[] = {0x00000000, 0xFFFFF800};

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x7FF;
    uint32_t ElementZ = (pSource->v >> 21) & 0x7FF;

    XMVECTORF32 vResult = {
        (ElementX == 0x200) ? -1.f : ((float)(int16_t)(ElementX | SignExtendX[ElementX >> 9]) / 511.0f),
        (ElementY == 0x400) ? -1.f : ((float)(int16_t)(ElementY | SignExtendYZ[ElementY >> 10]) / 1023.0f),
        (ElementZ == 0x400) ? -1.f : ((float)(int16_t)(ElementZ | SignExtendYZ[ElementZ >> 10]) / 1023.0f),
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 DHenN3Mul = {1.0f/511.0f,1.0f/(1023.0f*1024.0f),1.0f/(1023.0f*1024.0f*2048.0f),0};
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskDHen3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMXorDHen3);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddDHen3);
    // Normalize x,y and z to -1.0f-1.0f
    vResult = _mm_mul_ps(vResult,DHenN3Mul);
    // Clamp result (for case of -512/-1024)
    return _mm_max_ps( vResult, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadDHen3
(
    const XMDHEN3* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtendX[] = {0x00000000, 0xFFFFFC00};
    static const uint32_t SignExtendYZ[] = {0x00000000, 0xFFFFF800};

    uint32_t ElementX = pSource->v & 0x3FF;
    uint32_t ElementY = (pSource->v >> 10) & 0x7FF;
    uint32_t ElementZ = (pSource->v >> 21) & 0x7FF;

    XMVECTORF32 vResult = {
        (float)(int16_t)(ElementX | SignExtendX[ElementX >> 9]),
        (float)(int16_t)(ElementY | SignExtendYZ[ElementY >> 10]),
        (float)(int16_t)(ElementZ | SignExtendYZ[ElementZ >> 10]),
        0.f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Get the 32 bit value and splat it
    XMVECTOR vResult = _mm_load_ps1(reinterpret_cast<const float *>(&pSource->v));
    // Mask off x, y and z
    vResult = _mm_and_ps(vResult,g_XMMaskDHen3);
    // Convert x and y to unsigned
    vResult = _mm_xor_ps(vResult,g_XMXorDHen3);
    // Convert to float
    vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
    // Convert x and y back to signed
    vResult = _mm_add_ps(vResult,g_XMAddDHen3);
    // Normalize x to -210-511.0f and y and z to -1024-1023.0f
    vResult = _mm_mul_ps(vResult,g_XMMulDHen3);
    return vResult;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadXIcoN4
(
    const XMXICON4* pSource
)
{
    assert(pSource);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFF00000};

    uint32_t ElementX = (uint32_t)(pSource->v & 0xFFFFF);
    uint32_t ElementY = (uint32_t)((pSource->v >> 20) & 0xFFFFF);
    uint32_t ElementZ = (uint32_t)((pSource->v >> 40) & 0xFFFFF);

    XMVECTORF32 vResult = {
        (ElementX == 0x80000) ? -1.f : ((float)(int32_t)(ElementX | SignExtend[ElementX >> 19]) / 524287.0f),
        (ElementY == 0x80000) ? -1.f : ((float)(int32_t)(ElementY | SignExtend[ElementY >> 19]) / 524287.0f),
        (ElementZ == 0x80000) ? -1.f : ((float)(int32_t)(ElementZ | SignExtend[ElementZ >> 19]) / 524287.0f),
        (float)(pSource->v >> 60) / 15.0f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadXIcoN4Mul = {1.0f/524287.0f,1.0f/(524287.0f*4096.0f),1.0f/524287.0f,1.0f/(15.0f*4096.0f*65536.0f)};
    // Grab the 64 bit structure
    __m128d vResultd = _mm_load_sd(reinterpret_cast<const double *>(&pSource->v));
    // By shifting down 8 bits, y and z are in seperate 32 bit elements
    __m128i vResulti = _mm_srli_si128(_mm_castpd_si128(vResultd),8/8);
    // vResultd has x and w, vResulti has y and z, merge into one as x,w,y,z
    XMVECTOR vTemp = _mm_shuffle_ps(_mm_castsi128_ps(vResultd),_mm_castsi128_ps(vResulti),_MM_SHUFFLE(1,0,1,0));
    // Fix the entries to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,3,2,0));
    // Mask x,y,z and w
    vTemp = _mm_and_ps(vTemp,g_XMMaskIco4);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorXIco4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddXIco4);
    // Fix y and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadXIcoN4Mul);
    // Clamp result (for case of -524288)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadXIco4
(
    const XMXICO4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFF00000};

    uint32_t ElementX = (uint32_t)(pSource->v & 0xFFFFF);
    uint32_t ElementY = (uint32_t)((pSource->v >> 20) & 0xFFFFF);
    uint32_t ElementZ = (uint32_t)((pSource->v >> 40) & 0xFFFFF);

    XMVECTORF32 vResult = {
        (float)(int32_t)(ElementX | SignExtend[ElementX >> 19]),
        (float)(int32_t)(ElementY | SignExtend[ElementY >> 19]),
        (float)(int32_t)(ElementZ | SignExtend[ElementZ >> 19]),
        (float)(pSource->v >> 60)
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Grab the 64 bit structure
    __m128d vResultd = _mm_load_sd(reinterpret_cast<const double *>(&pSource->v));
    // By shifting down 8 bits, y and z are in seperate 32 bit elements
    __m128i vResulti = _mm_srli_si128(_mm_castpd_si128(vResultd),8/8);
    // vResultd has x and w, vResulti has y and z, merge into one as x,w,y,z
    XMVECTOR vTemp = _mm_shuffle_ps(_mm_castsi128_ps(vResultd),_mm_castsi128_ps(vResulti),_MM_SHUFFLE(1,0,1,0));
    // Fix the entries to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,3,2,0));
    // Mask x,y,z and w
    vTemp = _mm_and_ps(vTemp,g_XMMaskIco4);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorXIco4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddXIco4);
    // Fix y and w because they are too large
    vTemp = _mm_mul_ps(vTemp,g_XMMulIco4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadUIcoN4
(
    const XMUICON4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)(pSource->v & 0xFFFFF) / 1048575.0f,
        (float)((pSource->v >> 20) & 0xFFFFF) / 1048575.0f,
        (float)((pSource->v >> 40) & 0xFFFFF) / 1048575.0f,
        (float)(pSource->v >> 60) / 15.0f
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadUIcoN4Mul = {1.0f/1048575.0f,1.0f/(1048575.0f*4096.0f),1.0f/1048575.0f,1.0f/(15.0f*4096.0f*65536.0f)};
    // Grab the 64 bit structure
    __m128d vResultd = _mm_load_sd(reinterpret_cast<const double *>(&pSource->v));
    // By shifting down 8 bits, y and z are in seperate 32 bit elements
    __m128i vResulti = _mm_srli_si128(_mm_castpd_si128(vResultd),8/8);
    // vResultd has x and w, vResulti has y and z, merge into one as x,w,y,z
    XMVECTOR vTemp = _mm_shuffle_ps(_mm_castsi128_ps(vResultd),_mm_castsi128_ps(vResulti),_MM_SHUFFLE(1,0,1,0));
    // Fix the entries to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,3,2,0));
    // Mask x,y,z and w
    vTemp = _mm_and_ps(vTemp,g_XMMaskIco4);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipYW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddUIco4);
    // Fix y and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadUIcoN4Mul);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadUIco4
(
    const XMUICO4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    XMVECTORF32 vResult = {
        (float)(pSource->v & 0xFFFFF),
        (float)((pSource->v >> 20) & 0xFFFFF),
        (float)((pSource->v >> 40) & 0xFFFFF),
        (float)(pSource->v >> 60)
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Grab the 64 bit structure
    __m128d vResultd = _mm_load_sd(reinterpret_cast<const double *>(&pSource->v));
    // By shifting down 8 bits, y and z are in seperate 32 bit elements
    __m128i vResulti = _mm_srli_si128(_mm_castpd_si128(vResultd),8/8);
    // vResultd has x and w, vResulti has y and z, merge into one as x,w,y,z
    XMVECTOR vTemp = _mm_shuffle_ps(_mm_castsi128_ps(vResultd),_mm_castsi128_ps(vResulti),_MM_SHUFFLE(1,0,1,0));
    // Fix the entries to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,3,2,0));
    // Mask x,y,z and w
    vTemp = _mm_and_ps(vTemp,g_XMMaskIco4);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMFlipYW);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddUIco4);
    // Fix y and w because they are too large
    vTemp = _mm_mul_ps(vTemp,g_XMMulIco4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadIcoN4
(
    const XMICON4* pSource
)
{
    assert(pSource);

#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFF00000};
    static const uint32_t SignExtendW[] = {0x00000000, 0xFFFFFFF0};

    uint32_t ElementX = (uint32_t)(pSource->v & 0xFFFFF);
    uint32_t ElementY = (uint32_t)((pSource->v >> 20) & 0xFFFFF);
    uint32_t ElementZ = (uint32_t)((pSource->v >> 40) & 0xFFFFF);
    uint32_t ElementW = (uint32_t)(pSource->v >> 60);

    XMVECTORF32 vResult = {
        (ElementX == 0x80000) ? -1.f : ((float)(int32_t)(ElementX | SignExtend[ElementX >> 19]) / 524287.0f),
        (ElementY == 0x80000) ? -1.f : ((float)(int32_t)(ElementY | SignExtend[ElementY >> 19]) / 524287.0f),
        (ElementZ == 0x80000) ? -1.f : ((float)(int32_t)(ElementZ | SignExtend[ElementZ >> 19]) / 524287.0f),
        (ElementW == 0x80) ? -1.f : ((float)(int32_t)(ElementW | SignExtendW[ElementW >> 3]) / 7.0f)
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 LoadIcoN4Mul = {1.0f/524287.0f,1.0f/(524287.0f*4096.0f),1.0f/524287.0f,1.0f/(7.0f*4096.0f*65536.0f)};
    // Grab the 64 bit structure
    __m128d vResultd = _mm_load_sd(reinterpret_cast<const double *>(&pSource->v));
    // By shifting down 8 bits, y and z are in seperate 32 bit elements
    __m128i vResulti = _mm_srli_si128(_mm_castpd_si128(vResultd),8/8);
    // vResultd has x and w, vResulti has y and z, merge into one as x,w,y,z
    XMVECTOR vTemp = _mm_shuffle_ps(_mm_castsi128_ps(vResultd),_mm_castsi128_ps(vResulti),_MM_SHUFFLE(1,0,1,0));
    // Fix the entries to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,3,2,0));
    // Mask x,y,z and w
    vTemp = _mm_and_ps(vTemp,g_XMMaskIco4);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorIco4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddIco4);
    // Fix y and w because they are too large
    vTemp = _mm_mul_ps(vTemp,LoadIcoN4Mul);
    // Clamp result (for case of -524288/-8)
    return _mm_max_ps( vTemp, g_XMNegativeOne );
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline XMVECTOR PackedVector::Xbox::XMLoadIco4
(
    const XMICO4* pSource
)
{
    assert(pSource);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
    static const uint32_t SignExtend[] = {0x00000000, 0xFFF00000};
    static const uint32_t SignExtendW[] = {0x00000000, 0xFFFFFFF0};

    uint32_t ElementX = (uint32_t)(pSource->v & 0xFFFFF);
    uint32_t ElementY = (uint32_t)((pSource->v >> 20) & 0xFFFFF);
    uint32_t ElementZ = (uint32_t)((pSource->v >> 40) & 0xFFFFF);
    uint32_t ElementW = (uint32_t)(pSource->v >> 60);

    XMVECTORF32 vResult = {
        (float)(int32_t)(ElementX | SignExtend[ElementX >> 19]),
        (float)(int32_t)(ElementY | SignExtend[ElementY >> 19]),
        (float)(int32_t)(ElementZ | SignExtend[ElementZ >> 19]),
        (float)(int32_t)(ElementW | SignExtendW[ElementW >> 3])
    };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Grab the 64 bit structure
    __m128d vResultd = _mm_load_sd(reinterpret_cast<const double *>(&pSource->v));
    // By shifting down 8 bits, y and z are in seperate 32 bit elements
    __m128i vResulti = _mm_srli_si128(_mm_castpd_si128(vResultd),8/8);
    // vResultd has x and w, vResulti has y and z, merge into one as x,w,y,z
    XMVECTOR vTemp = _mm_shuffle_ps(_mm_castsi128_ps(vResultd),_mm_castsi128_ps(vResulti),_MM_SHUFFLE(1,0,1,0));
    // Fix the entries to x,y,z,w
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,3,2,0));
    // Mask x,y,z and w
    vTemp = _mm_and_ps(vTemp,g_XMMaskIco4);
    // x and z are unsigned! Flip the bits to convert the order to signed
    vTemp = _mm_xor_ps(vTemp,g_XMXorIco4);
    // Convert to floating point numbers
    vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
    // x and z - 0x80 to complete the conversion
    vTemp = _mm_add_ps(vTemp,g_XMAddIco4);
    // Fix y and w because they are too large
    vTemp = _mm_mul_ps(vTemp,g_XMMulIco4);
    return vTemp;
#elif defined(XM_NO_MISALIGNED_VECTOR_ACCESS)
#endif // _XM_VMX128_INTRINSICS_
}
///end_xbox360


/****************************************************************************
 *
 * Vector and matrix store operations
 *
 ****************************************************************************/
_Use_decl_annotations_
inline void PackedVector::XMStoreColor
(
    XMCOLOR* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {255.0f, 255.0f, 255.0f, 255.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->c = ((uint32_t)tmp.w << 24) |
                      ((uint32_t)tmp.x << 16) |
                      ((uint32_t)tmp.y <<  8) |
                      ((uint32_t)tmp.z);

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32  Scale = {255.0f,255.0f,255.0f,255.0f};
    // Set <0 to 0
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    // Set>1 to 1
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Convert to 0-255
    vResult = _mm_mul_ps(vResult,Scale);
    // Shuffle RGBA to ARGB
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,0,1,2));
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
_Use_decl_annotations_
inline void PackedVector::XMStoreHalf2
(
    XMHALF2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    pDestination->x = XMConvertFloatToHalf(XMVectorGetX(V));
    pDestination->y = XMConvertFloatToHalf(XMVectorGetY(V));

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreShortN2
(
    XMSHORTN2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    vResult = _mm_mul_ps(vResult,Scale);
    __m128i vResulti = _mm_cvtps_epi32(vResult);
    vResulti = _mm_packs_epi32(vResulti,vResulti);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->x),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreShort2
(
    XMSHORT2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32 Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32 Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,Min);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Pack the ints into shorts
    vInt = _mm_packs_epi32(vInt,vInt);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->x),_mm_castsi128_ps(vInt));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUShortN2
(
    XMUSHORTN2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);
    N = XMVectorTruncate(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;

#elif defined(_XM_SSE_INTRINSICS_)
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
_Use_decl_annotations_
inline void PackedVector::XMStoreUShort2
(
    XMUSHORT2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;

#elif defined(_XM_SSE_INTRINSICS_)
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
_Use_decl_annotations_
inline void PackedVector::XMStoreByteN2
(
    XMBYTEN2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);

    static const XMVECTORF32  Scale = {127.0f, 127.0f, 127.0f, 127.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int8_t)tmp.x;
    pDestination->y = (int8_t)tmp.y;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreByte2
(
    XMBYTE2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);

    static const XMVECTORF32 Min = {-127.0f, -127.0f, -127.0f, -127.0f};
    static const XMVECTORF32 Max = {127.0f, 127.0f, 127.0f, 127.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (int8_t)tmp.x;
    pDestination->y = (int8_t)tmp.y;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUByteN2
(
    XMUBYTEN2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);

    static const XMVECTORF32  Scale = {255.0f, 255.0f, 255.0f, 255.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);
    N = XMVectorTruncate(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (uint8_t)tmp.x;
    pDestination->y = (uint8_t)tmp.y;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUByte2
(
    XMUBYTE2* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);

    static const XMVECTORF32 Max = {255.0f, 255.0f, 255.0f, 255.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->x = (uint8_t)tmp.x;
    pDestination->y = (uint8_t)tmp.y;
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreU565
(
    XMU565* pDestination,
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
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
    static const XMVECTORF32  Max = {31.0f, 63.0f, 31.0f, 0.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A( &tmp, N );

    pDestination->v = (((uint16_t)tmp.z & 0x1F) << 11) |
                      (((uint16_t)tmp.y & 0x3F) << 5) |
                      (((uint16_t)tmp.x & 0x1F));
#endif !_XM_SSE_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreFloat3PK
(
    XMFLOAT3PK* pDestination,
    FXMVECTOR V
)
{
    assert(pDestination);

    __declspec(align(16)) uint32_t IValue[4];
    XMStoreFloat3A( reinterpret_cast<XMFLOAT3A*>(&IValue), V );

    uint32_t Result[3];

    // X & Y Channels (5-bit exponent, 6-bit mantissa)
    for(uint32_t j=0; j < 2; ++j)
    {
        uint32_t Sign = IValue[j] & 0x80000000;
        uint32_t I = IValue[j] & 0x7FFFFFFF;

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
    uint32_t Sign = IValue[2] & 0x80000000;
    uint32_t I = IValue[2] & 0x7FFFFFFF;

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
_Use_decl_annotations_
inline void PackedVector::XMStoreFloat3SE
(
    XMFLOAT3SE* pDestination,
    FXMVECTOR V
)
{
    assert(pDestination);

    __declspec(align(16)) uint32_t IValue[4];
    XMStoreFloat3A( reinterpret_cast<XMFLOAT3A*>(&IValue), V );

    uint32_t Exp[3];
    uint32_t Frac[3];

    // X, Y, Z Channels (5-bit exponent, 9-bit mantissa)
    for(uint32_t j=0; j < 3; ++j)
    {
        uint32_t Sign = IValue[j] & 0x80000000;
        uint32_t I = IValue[j] & 0x7FFFFFFF;

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
     
            uint32_t T = ((I + 0x1FFFU + ((I >> 14U) & 1U)) >> 14U)&0x3fffU;

            Exp[j] = (T & 0x3E00) >> 9;
            Frac[j] = T & 0x1ff;
        }
    }

    // Adjust to a shared exponent
    uint32_t T = XMMax( Exp[0], XMMax( Exp[1], Exp[2] ) );

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
_Use_decl_annotations_
inline void PackedVector::XMStoreHalf4
(
    XMHALF4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_SSE_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)
 
    XMFLOAT4A t;
    XMStoreFloat4A(&t, V );

    pDestination->x = XMConvertFloatToHalf(t.x);
    pDestination->y = XMConvertFloatToHalf(t.y);
    pDestination->z = XMConvertFloatToHalf(t.z);
    pDestination->w = XMConvertFloatToHalf(t.w);

#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreShortN4
(
    XMSHORTN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)

    static const XMVECTORF32  Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;
    pDestination->z = (int16_t)tmp.z;
    pDestination->w = (int16_t)tmp.w;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vmaxq_f32( V, g_XMNegativeOne );
    vResult = vminq_f32( vResult, g_XMOne );
    const __n128 Scale = vdupq_n_f32( 32767.0f );
    vResult = vmulq_f32( vResult, Scale );
    vResult = vcvtq_s32_f32( vResult );
    __n64 vInt = vmovn_s32( vResult );
    vst1_s16( (int16_t*)pDestination, vInt );
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Scale = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    vResult = _mm_mul_ps(vResult,Scale);
    __m128i vResulti = _mm_cvtps_epi32(vResult);
    vResulti = _mm_packs_epi32(vResulti,vResulti);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->x),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreShort4
(
    XMSHORT4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)

    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32 Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;
    pDestination->z = (int16_t)tmp.z;
    pDestination->w = (int16_t)tmp.w;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32 Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};

    __n128 vResult = vmaxq_f32( V, Min );
    vResult = vminq_f32( vResult, Max );
    vResult = vcvtq_s32_f32( vResult );
    __n64 vInt = vmovn_s32( vResult );
    vst1_s16( (int16_t*)pDestination, vInt );
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Min = {-32767.0f, -32767.0f, -32767.0f, -32767.0f};
    static const XMVECTORF32 Max = {32767.0f, 32767.0f, 32767.0f, 32767.0f};
    // Bounds check
    XMVECTOR vResult = _mm_max_ps(V,Min);
    vResult = _mm_min_ps(vResult,Max);
     // Convert to int with rounding
    __m128i vInt = _mm_cvtps_epi32(vResult);
    // Pack the ints into shorts
    vInt = _mm_packs_epi32(vInt,vInt);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->x),_mm_castsi128_pd(vInt));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUShortN4
(
    XMUSHORTN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)

    static const XMVECTORF32  Scale = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);
    N = XMVectorTruncate(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;
    pDestination->z = (int16_t)tmp.z;
    pDestination->w = (int16_t)tmp.w;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    __n128 vResult = vmaxq_f32( V, g_XMZero );
    vResult = vminq_f32( vResult, g_XMOne );
    const __n128 Scale = vdupq_n_f32( 65535.0f );
    vResult = vmulq_f32( vResult, Scale );
    vResult = vcvtq_u32_f32( vResult );
    __n64 vInt = vmovn_u32( vResult );
    vst1_u16( (uint16_t*)pDestination, vInt );
#elif defined(_XM_SSE_INTRINSICS_)
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
_Use_decl_annotations_
inline void PackedVector::XMStoreUShort4
(
    XMUSHORT4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)

    static const XMVECTORF32 Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (int16_t)tmp.x;
    pDestination->y = (int16_t)tmp.y;
    pDestination->z = (int16_t)tmp.z;
    pDestination->w = (int16_t)tmp.w;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
    static const XMVECTORF32 Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};

    __n128 vResult = vmaxq_f32( V, g_XMZero );
    vResult = vminq_f32( vResult, Max );
    vResult = vcvtq_u32_f32( vResult );
    __n64 vInt = vmovn_u32( vResult );
    vst1_u16( (uint16_t*)pDestination, vInt );
#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Max = {65535.0f, 65535.0f, 65535.0f, 65535.0f};
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
_Use_decl_annotations_
inline void PackedVector::XMStoreXDecN4
(
    XMXDECN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Min = {-1.0f, -1.0f, -1.0f, 0.0f};
    static const XMVECTORF32  Scale = {511.0f, 511.0f, 511.0f, 3.0f};

    XMVECTOR N = XMVectorClamp(V, Min.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint32_t)tmp.w << 30) |
                       (((int32_t)tmp.z & 0x3FF) << 20) |
                       (((int32_t)tmp.y & 0x3FF) << 10) |
                       (((int32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 Min = {-1.0f, -1.0f, -1.0f, 0.0f};
    static const XMVECTORF32 Scale = {511.0f, 511.0f*1024.0f, 511.0f*1048576.0f,3.0f*536870912.0f};
    static const XMVECTORI32 ScaleMask = {0x3FF,0x3FF<<10,0x3FF<<20,0x3<<29};
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
    vResult = XM_PERMUTE_PS(_mm_castsi128_ps(vResulti),_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreXDec4
(
    XMXDEC4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-511.0f, -511.0f, -511.0f, 0.0f};
    static const XMVECTORF32 Max = {511.0f, 511.0f, 511.0f, 3.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint32_t)tmp.w << 30) |
                       (((int32_t)tmp.z & 0x3FF) << 20) |
                       (((int32_t)tmp.y & 0x3FF) << 10) |
                       (((int32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUDecN4
(
    XMUDECN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {1023.0f, 1023.0f, 1023.0f, 3.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint32_t)tmp.w << 30) |
                       (((uint32_t)tmp.z & 0x3FF) << 20) |
                       (((uint32_t)tmp.y & 0x3FF) << 10) |
                       (((uint32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUDec4
(
    XMUDEC4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Max = {1023.0f, 1023.0f, 1023.0f, 3.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint32_t)tmp.w << 30) |
                       (((uint32_t)tmp.z & 0x3FF) << 20) |
                       (((uint32_t)tmp.y & 0x3FF) << 10) |
                       (((uint32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreDecN4
(
    XMDECN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {511.0f, 511.0f, 511.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((int32_t)tmp.w << 30) |
                       (((int32_t)tmp.z & 0x3FF) << 20) |
                       (((int32_t)tmp.y & 0x3FF) << 10) |
                       (((int32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreDec4
(
    XMDEC4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-511.0f, -511.0f, -511.0f, -1.0f};
    static const XMVECTORF32 Max = {511.0f, 511.0f, 511.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((int32_t)tmp.w << 30) |
                       (((int32_t)tmp.z & 0x3FF) << 20) |
                       (((int32_t)tmp.y & 0x3FF) << 10) |
                       (((int32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUByteN4
(
    XMUBYTEN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {255.0f, 255.0f, 255.0f, 255.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (uint8_t)tmp.x;
    pDestination->y = (uint8_t)tmp.y;
    pDestination->z = (uint8_t)tmp.z;
    pDestination->w = (uint8_t)tmp.w;

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUByte4
(
    XMUBYTE4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Max = {255.0f, 255.0f, 255.0f, 255.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (uint8_t)tmp.x;
    pDestination->y = (uint8_t)tmp.y;
    pDestination->z = (uint8_t)tmp.z;
    pDestination->w = (uint8_t)tmp.w;

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreByteN4
(
    XMBYTEN4* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {127.0f, 127.0f, 127.0f, 127.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(V, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (int8_t)tmp.x;
    pDestination->y = (int8_t)tmp.y;
    pDestination->z = (int8_t)tmp.z;
    pDestination->w = (int8_t)tmp.w;

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreByte4
(
    XMBYTE4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-127.0f, -127.0f, -127.0f, -127.0f};
    static const XMVECTORF32 Max = {127.0f, 127.0f, 127.0f, 127.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->x = (int8_t)tmp.x;
    pDestination->y = (int8_t)tmp.y;
    pDestination->z = (int8_t)tmp.z;
    pDestination->w = (int8_t)tmp.w;

#elif defined(_XM_SSE_INTRINSICS_)
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
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreUNibble4
(
     XMUNIBBLE4* pDestination,
     FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
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
    static const XMVECTORF32  Max = {15.0f,15.0f,15.0f,15.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((uint16_t)tmp.w & 0xF) << 12) |
                      (((uint16_t)tmp.z & 0xF) << 8) |
                      (((uint16_t)tmp.y & 0xF) << 4) |
                      (((uint16_t)tmp.x & 0xF));
#endif !_XM_SSE_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::XMStoreU555
(
     XMU555* pDestination,
     FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
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
    static const XMVECTORF32  Max = {31.0f, 31.0f, 31.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((tmp.w > 0.f) ? 0x8000 : 0) |
                      (((uint16_t)tmp.z & 0x1F) << 10) |
                      (((uint16_t)tmp.y & 0x1F) << 5) |
                      (((uint16_t)tmp.x & 0x1F));
#endif !_XM_SSE_INTRINSICS_
}

///begin_xbox360
////////////////////////////////////////////////////////////////////////////////
// PackedVector::Xbox
////////////////////////////////////////////////////////////////////////////////
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreUHenDN3
(
    XMUHENDN3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {2047.0f, 2047.0f, 1023.0f, 0.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((uint32_t)tmp.z & 0x3FF) << 22) |
                      (((uint32_t)tmp.y & 0x7FF) << 11) |
                      (((uint32_t)tmp.x & 0x7FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 ScaleUHenDN3 = {2047.0f, 2047.0f*2048.0f,1023.0f*(2048.0f*2048.0f)/2.0f,1.0f};
    static const XMVECTORI32 MaskUHenDN3 = {0x7FF,0x7FF<<11,0x3FF<<(22-1),0};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUHenDN3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUHenDN3);
    // Do a horizontal or of 3 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(0,3,2,1));
    // i = x|y
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti2,_MM_SHUFFLE(0,3,2,1));
    // Add Z to itself to perform a single bit left shift
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreUHenD3
(
    XMUHEND3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Max = {2047.0f, 2047.0f, 1023.0f, 0.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((uint32_t)tmp.z & 0x3FF) << 22) |
                      (((uint32_t)tmp.y & 0x7FF) << 11) |
                      (((uint32_t)tmp.x & 0x7FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 MaxUHenD3 = { 2047.0f, 2047.0f, 1023.0f, 1.0f};
    static const XMVECTORF32 ScaleUHenD3 = {1.0f, 2048.0f,(2048.0f*2048.0f)/2.0f,1.0f};
    static const XMVECTORI32 MaskUHenD3 = {0x7FF,0x7FF<<11,0x3FF<<(22-1),0};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,MaxUHenD3);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUHenD3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUHenD3);
    // Do a horizontal or of 3 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(0,3,2,1));
    // i = x|y
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti2,_MM_SHUFFLE(0,3,2,1));
    // Add Z to itself to perform a single bit left shift
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreHenDN3
(
    XMHENDN3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {1023.0f, 1023.0f, 511.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((int32_t)tmp.z & 0x3FF) << 22) |
                      (((int32_t)tmp.y & 0x7FF) << 11) |
                      (((int32_t)tmp.x & 0x7FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 ScaleHenDN3 = {1023.0f, 1023.0f*2048.0f,511.0f*(2048.0f*2048.0f),1.0f};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleHenDN3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,g_XMMaskHenD3);
    // Do a horizontal or of all 4 entries
    vResult = XM_PERMUTE_PS(_mm_castsi128_ps(vResulti),_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreHenD3
(
    XMHEND3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-1023.0f, -1023.0f, -511.0f, -1.0f};
    static const XMVECTORF32 Max = {1023.0f, 1023.0f, 511.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((int32_t)tmp.z & 0x3FF) << 22) |
                      (((int32_t)tmp.y & 0x7FF) << 11) |
                      (((int32_t)tmp.x & 0x7FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 MinHenD3 = {-1023.0f,-1023.0f,-511.0f,-1.0f};
    static const XMVECTORF32 MaxHenD3 = { 1023.0f, 1023.0f, 511.0f, 1.0f};
    static const XMVECTORF32 ScaleHenD3 = {1.0f, 2048.0f,(2048.0f*2048.0f),1.0f};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,MinHenD3);
    vResult = _mm_min_ps(vResult,MaxHenD3);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleHenD3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,g_XMMaskHenD3);
    // Do a horizontal or of all 4 entries
    vResult = XM_PERMUTE_PS(_mm_castsi128_ps(vResulti),_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreUDHenN3
(
    XMUDHENN3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {1023.0f, 2047.0f, 2047.0f, 0.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiply(N, Scale.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((uint32_t)tmp.z & 0x7FF) << 21) |
                      (((uint32_t)tmp.y & 0x7FF) << 10) |
                      (((uint32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 ScaleUDHenN3 = {1023.0f,2047.0f*1024.0f,2047.0f*(1024.0f*2048.0f)/2.0f,1.0f};
    static const XMVECTORI32 MaskUDHenN3 = {0x3FF,0x7FF<<10,0x7FF<<(21-1),0};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUDHenN3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUDHenN3);
    // Do a horizontal or of 3 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(0,3,2,1));
    // i = x|y
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti2,_MM_SHUFFLE(0,3,2,1));
    // Add Z to itself to perform a single bit left shift
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreUDHen3
(
    XMUDHEN3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Max = {1023.0f, 2047.0f, 2047.0f, 0.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((uint32_t)tmp.z & 0x7FF) << 21) |
                      (((uint32_t)tmp.y & 0x7FF) << 10) |
                      (((uint32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 MaxUDHen3 = { 1023.0f, 2047.0f, 2047.0f, 1.0f};
    static const XMVECTORF32 ScaleUDHen3 = {1.0f, 1024.0f,(1024.0f*2048.0f)/2.0f,1.0f};
    static const XMVECTORI32 MaskUDHen3 = {0x3FF,0x7FF<<10,0x7FF<<(21-1),0};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMZero);
    vResult = _mm_min_ps(vResult,MaxUDHen3);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUDHen3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUDHen3);
    // Do a horizontal or of 3 entries
    __m128i vResulti2 = _mm_shuffle_epi32(vResulti,_MM_SHUFFLE(0,3,2,1));
    // i = x|y
    vResulti = _mm_or_si128(vResulti,vResulti2);
    // Move Z to the x position
    vResulti2 = _mm_shuffle_epi32(vResulti2,_MM_SHUFFLE(0,3,2,1));
    // Add Z to itself to perform a single bit left shift
    vResulti2 = _mm_add_epi32(vResulti2,vResulti2);
    // i = x|y|z
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreDHenN3
(
    XMDHENN3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {511.0f, 1023.0f, 1023.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((int32_t)tmp.z & 0x7FF) << 21) |
                      (((int32_t)tmp.y & 0x7FF) << 10) |
                      (((int32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 ScaleDHenN3 = {511.0f, 1023.0f*1024.0f,1023.0f*(1024.0f*2048.0f),1.0f};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleDHenN3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,g_XMMaskDHen3);
    // Do a horizontal or of all 4 entries
    vResult = XM_PERMUTE_PS(_mm_castsi128_ps(vResulti),_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreDHen3
(
    XMDHEN3* pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-511.0f, -1023.0f, -1023.0f, -1.0f};
    static const XMVECTORF32 Max = {511.0f, 1023.0f, 1023.0f, 1.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = (((int32_t)tmp.z & 0x7FF) << 21) |
                      (((int32_t)tmp.y & 0x7FF) << 10) |
                      (((int32_t)tmp.x & 0x3FF));

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 MinDHen3 = {-511.0f,-1023.0f,-1023.0f,-1.0f};
    static const XMVECTORF32 MaxDHen3 = { 511.0f, 1023.0f, 1023.0f, 1.0f};
    static const XMVECTORF32 ScaleDHen3 = {1.0f, 1024.0f,(1024.0f*2048.0f),1.0f};
    // Clamp to bounds
    XMVECTOR vResult = _mm_max_ps(V,MinDHen3);
    vResult = _mm_min_ps(vResult,MaxDHen3);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleDHen3);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,g_XMMaskDHen3);
    // Do a horizontal or of all 4 entries
    vResult = XM_PERMUTE_PS(_mm_castsi128_ps(vResulti),_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(0,3,2,1));
    vResulti = _mm_or_si128(vResulti,_mm_castps_si128(vResult));
    _mm_store_ss(reinterpret_cast<float *>(&pDestination->v),_mm_castsi128_ps(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreXIcoN4
(
    XMXICON4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Min = {-1.0f, -1.0f, -1.0f, 0.0f};
    static const XMVECTORF32  Scale = {524287.0f, 524287.0f, 524287.0f, 15.0f};

    XMVECTOR N = XMVectorClamp(V, Min.v, g_XMOne.v);
    N = XMVectorMultiply(N, Scale.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint64_t)tmp.w << 60) |
                       (((int64_t)tmp.z & 0xFFFFF) << 40) |
                       (((int64_t)tmp.y & 0xFFFFF) << 20) |
                       (((int64_t)tmp.x & 0xFFFFF));

#elif defined(_XM_SSE_INTRINSICS_)
    // Note: Masks are x,w,y and z
    static const XMVECTORF32 MinXIcoN4 = {-1.0f, 0.0f,-1.0f,-1.0f};
    static const XMVECTORF32 ScaleXIcoN4 = {524287.0f,15.0f*4096.0f*65536.0f*0.5f,524287.0f*4096.0f,524287.0f};
    static const XMVECTORI32 MaskXIcoN4 = {0xFFFFF,0xF<<((60-32)-1),0xFFFFF000,0xFFFFF};

    // Clamp to bounds
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,1,3,0));
    vResult = _mm_max_ps(vResult,MinXIcoN4);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleXIcoN4);
    // Convert to integer (w is unsigned)
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off unused bits
    vResulti = _mm_and_si128(vResulti,MaskXIcoN4);
    // Isolate Y
    __m128i vResulti2 = _mm_and_si128(vResulti,g_XMMaskY);
    // Double Y (Really W) to fixup for unsigned conversion
    vResulti = _mm_add_epi32(vResulti,vResulti2);
    // Shift y and z to straddle the 32-bit boundary
    vResulti2 = _mm_srli_si128(vResulti,(64+12)/8);
    // Shift it into place
    vResulti2 = _mm_slli_si128(vResulti2,20/8);
    // i = x|y<<20|z<<40|w<<60
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->v),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreXIco4
(
    XMXICO4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-524287.0f, -524287.0f, -524287.0f, 0.0f};
    static const XMVECTORF32 Max = {524287.0f, 524287.0f, 524287.0f, 15.0f};

    XMVECTOR N = XMVectorClamp(V, Min.v, Max.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint64_t)tmp.w << 60) |
                       (((int64_t)tmp.z & 0xFFFFF) << 40) |
                       (((int64_t)tmp.y & 0xFFFFF) << 20) |
                       (((int64_t)tmp.x & 0xFFFFF));

#elif defined(_XM_SSE_INTRINSICS_)
    // Note: Masks are x,w,y and z
    static const XMVECTORF32 MinXIco4 = {-524287.0f, 0.0f,-524287.0f,-524287.0f};
    static const XMVECTORF32 MaxXIco4 = { 524287.0f,15.0f, 524287.0f, 524287.0f};
    static const XMVECTORF32 ScaleXIco4 = {1.0f,4096.0f*65536.0f*0.5f,4096.0f,1.0f};
    static const XMVECTORI32 MaskXIco4 = {0xFFFFF,0xF<<((60-1)-32),0xFFFFF000,0xFFFFF};
    // Clamp to bounds
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,1,3,0));
    vResult = _mm_max_ps(vResult,MinXIco4);
    vResult = _mm_min_ps(vResult,MaxXIco4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleXIco4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskXIco4);
    // Isolate Y
    __m128i vResulti2 = _mm_and_si128(vResulti,g_XMMaskY);
    // Double Y (Really W) to fixup for unsigned conversion
    vResulti = _mm_add_epi32(vResulti,vResulti2);
    // Shift y and z to straddle the 32-bit boundary
    vResulti2 = _mm_srli_si128(vResulti,(64+12)/8);
    // Shift it into place
    vResulti2 = _mm_slli_si128(vResulti2,20/8);
    // i = x|y<<20|z<<40|w<<60
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->v),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreUIcoN4
(
    XMUICON4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Scale = {1048575.0f, 1048575.0f, 1048575.0f, 15.0f};

    XMVECTOR N = XMVectorSaturate(V);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMOneHalf.v);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint64_t)tmp.w << 60) |
                       (((uint64_t)tmp.z & 0xFFFFF) << 40) |
                       (((uint64_t)tmp.y & 0xFFFFF) << 20) |
                       (((uint64_t)tmp.x & 0xFFFFF));

#elif defined(_XM_SSE_INTRINSICS_)
    // Note: Masks are x,w,y and z
    static const XMVECTORF32 ScaleUIcoN4 = {1048575.0f,15.0f*4096.0f*65536.0f,1048575.0f*4096.0f,1048575.0f};
    static const XMVECTORI32 MaskUIcoN4 = {0xFFFFF,0xF<<(60-32),0xFFFFF000,0xFFFFF};
    static const XMVECTORF32 AddUIcoN4 = {0.0f,-32768.0f*65536.0f,-32768.0f*65536.0f,0.0f};
    // Clamp to bounds
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,1,3,0));
    vResult = _mm_max_ps(vResult,g_XMZero);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUIcoN4);
    // Adjust for unsigned entries
    vResult = _mm_add_ps(vResult,AddUIcoN4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Fix the signs on the unsigned entries
    vResulti = _mm_xor_si128(vResulti,g_XMFlipYZ);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUIcoN4);
    // Shift y and z to straddle the 32-bit boundary
    __m128i vResulti2 = _mm_srli_si128(vResulti,(64+12)/8);
    // Shift it into place
    vResulti2 = _mm_slli_si128(vResulti2,20/8);
    // i = x|y<<20|z<<40|w<<60
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->v),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_

}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreUIco4
(
    XMUICO4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Max = {1048575.0f, 1048575.0f, 1048575.0f, 15.0f};

    XMVECTOR N = XMVectorClamp(V, XMVectorZero(), Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint64_t)tmp.w << 60) |
                       (((uint64_t)tmp.z & 0xFFFFF) << 40) |
                       (((uint64_t)tmp.y & 0xFFFFF) << 20) |
                       (((uint64_t)tmp.x & 0xFFFFF));

#elif defined(_XM_SSE_INTRINSICS_)
    // Note: Masks are x,w,y and z
    static const XMVECTORF32 MaxUIco4 = { 1048575.0f, 15.0f, 1048575.0f, 1048575.0f};
    static const XMVECTORF32 ScaleUIco4 = {1.0f,4096.0f*65536.0f,4096.0f,1.0f};
    static const XMVECTORI32 MaskUIco4 = {0xFFFFF,0xF<<(60-32),0xFFFFF000,0xFFFFF};
    static const XMVECTORF32 AddUIco4 = {0.0f,-32768.0f*65536.0f,-32768.0f*65536.0f,0.0f};
    // Clamp to bounds
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,1,3,0));
    vResult = _mm_max_ps(vResult,g_XMZero);
    vResult = _mm_min_ps(vResult,MaxUIco4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleUIco4);
    vResult = _mm_add_ps(vResult,AddUIco4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    vResulti = _mm_xor_si128(vResulti,g_XMFlipYZ);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskUIco4);
    // Shift y and z to straddle the 32-bit boundary
    __m128i vResulti2 = _mm_srli_si128(vResulti,(64+12)/8);
    // Shift it into place
    vResulti2 = _mm_slli_si128(vResulti2,20/8);
    // i = x|y<<20|z<<40|w<<60
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->v),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreIcoN4
(
    XMICON4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32  Scale = {524287.0f, 524287.0f, 524287.0f, 7.0f};

    XMVECTOR N = XMVectorClamp(V, g_XMNegativeOne.v, g_XMOne.v);
    N = XMVectorMultiplyAdd(N, Scale.v, g_XMNegativeZero.v);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((uint64_t)tmp.w << 60) |
                       (((uint64_t)tmp.z & 0xFFFFF) << 40) |
                       (((uint64_t)tmp.y & 0xFFFFF) << 20) |
                       (((uint64_t)tmp.x & 0xFFFFF));

#elif defined(_XM_SSE_INTRINSICS_)
    // Note: Masks are x,w,y and z
    static const XMVECTORF32 ScaleIcoN4 = {524287.0f,7.0f*4096.0f*65536.0f,524287.0f*4096.0f,524287.0f};
    static const XMVECTORI32 MaskIcoN4 = {0xFFFFF,0xF<<(60-32),0xFFFFF000,0xFFFFF};
    // Clamp to bounds
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,1,3,0));
    vResult = _mm_max_ps(vResult,g_XMNegativeOne);
    vResult = _mm_min_ps(vResult,g_XMOne);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleIcoN4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskIcoN4);
    // Shift y and z to straddle the 32-bit boundary
    __m128i vResulti2 = _mm_srli_si128(vResulti,(64+12)/8);
    // Shift it into place
    vResulti2 = _mm_slli_si128(vResulti2,20/8);
    // i = x|y<<20|z<<40|w<<60
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->v),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline void PackedVector::Xbox::XMStoreIco4
(
    XMICO4*  pDestination, 
    FXMVECTOR V
)
{
    assert(pDestination);
#if defined(_XM_NO_INTRINSICS_) || defined(_XM_ARM_NEON_INTRINSICS_)

    static const XMVECTORF32 Min = {-524287.0f, -524287.0f, -524287.0f, -7.0f};
    static const XMVECTORF32 Max = {524287.0f, 524287.0f, 524287.0f, 7.0f};

    XMVECTOR N = XMVectorClamp(V, Min, Max);
    N = XMVectorRound(N);

    XMFLOAT4A tmp;
    XMStoreFloat4A(&tmp, N );

    pDestination->v = ((int64_t)tmp.w << 60) |
                       (((int64_t)tmp.z & 0xFFFFF) << 40) |
                       (((int64_t)tmp.y & 0xFFFFF) << 20) |
                       (((int64_t)tmp.x & 0xFFFFF));

#elif defined(_XM_SSE_INTRINSICS_)
    // Note: Masks are x,w,y and z
    static const XMVECTORF32 MinIco4 = {-524287.0f,-7.0f,-524287.0f,-524287.0f};
    static const XMVECTORF32 MaxIco4 = { 524287.0f, 7.0f, 524287.0f, 524287.0f};
    static const XMVECTORF32 ScaleIco4 = {1.0f,4096.0f*65536.0f,4096.0f,1.0f};
    static const XMVECTORI32 MaskIco4 = {0xFFFFF,0xF<<(60-32),0xFFFFF000,0xFFFFF};
    // Clamp to bounds
    XMVECTOR vResult = XM_PERMUTE_PS(V,_MM_SHUFFLE(2,1,3,0));
    vResult = _mm_max_ps(vResult,MinIco4);
    vResult = _mm_min_ps(vResult,MaxIco4);
    // Scale by multiplication
    vResult = _mm_mul_ps(vResult,ScaleIco4);
    // Convert to int
    __m128i vResulti = _mm_cvttps_epi32(vResult);
    // Mask off any fraction
    vResulti = _mm_and_si128(vResulti,MaskIco4);
    // Shift y and z to straddle the 32-bit boundary
    __m128i vResulti2 = _mm_srli_si128(vResulti,(64+12)/8);
    // Shift it into place
    vResulti2 = _mm_slli_si128(vResulti2,20/8);
    // i = x|y<<20|z<<40|w<<60
    vResulti = _mm_or_si128(vResulti,vResulti2);
    _mm_store_sd(reinterpret_cast<double *>(&pDestination->v),_mm_castsi128_pd(vResulti));
#else // _XM_VMX128_INTRINSICS_
#endif // _XM_VMX128_INTRINSICS_
}
///end_xbox360


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
_Use_decl_annotations_
inline PackedVector::XMCOLOR::XMCOLOR
(
    const float* pArray
)
{
    XMStoreColor(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
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
_Use_decl_annotations_
inline PackedVector::XMSHORTN2::XMSHORTN2
(
    const float* pArray
)
{
    XMStoreShortN2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMSHORT2::XMSHORT2
(
    const float* pArray
)
{
    XMStoreShort2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUSHORTN2::XMUSHORTN2
(
    const float* pArray
)
{
    XMStoreUShortN2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUSHORT2::XMUSHORT2
(
    const float* pArray
)
{
    XMStoreUShort2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMBYTEN2::XMBYTEN2
(
    const float* pArray
)
{
    XMStoreByteN2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMBYTE2::XMBYTE2
(
    const float* pArray
)
{
    XMStoreByte2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUBYTEN2::XMUBYTEN2
(
    const float* pArray
)
{
    XMStoreUByteN2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUBYTE2::XMUBYTE2
(
    const float* pArray
)
{
    XMStoreUByte2(this, XMLoadFloat2(reinterpret_cast<const XMFLOAT2*>(pArray)));
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

_Use_decl_annotations_
inline PackedVector::XMU565::XMU565
(
    const float *pArray
)
{
    XMStoreU565(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
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

_Use_decl_annotations_
inline PackedVector::XMFLOAT3PK::XMFLOAT3PK
(
    const float *pArray
)
{
    XMStoreFloat3PK(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
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

_Use_decl_annotations_
inline PackedVector::XMFLOAT3SE::XMFLOAT3SE
(
    const float *pArray
)
{
    XMStoreFloat3SE(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
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

_Use_decl_annotations_
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
_Use_decl_annotations_
inline PackedVector::XMSHORTN4::XMSHORTN4
(
    const float* pArray
)
{
    XMStoreShortN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMSHORT4::XMSHORT4
(
    const float* pArray
)
{
    XMStoreShort4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUSHORTN4::XMUSHORTN4
(
    const float* pArray
)
{
    XMStoreUShortN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUSHORT4::XMUSHORT4
(
    const float* pArray
)
{
    XMStoreUShort4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMXDECN4::XMXDECN4
(
    const float* pArray
)
{
    XMStoreXDecN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMXDEC4::XMXDEC4
(
    const float* pArray
)
{
    XMStoreXDec4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMDECN4::XMDECN4
(
    const float* pArray
)
{
    XMStoreDecN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMDEC4::XMDEC4
(
    const float* pArray
)
{
    XMStoreDec4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUDECN4::XMUDECN4
(
    const float* pArray
)
{
    XMStoreUDecN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUDEC4::XMUDEC4
(
    const float* pArray
)
{
    XMStoreUDec4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMBYTEN4::XMBYTEN4
(
    const float* pArray
)
{
    XMStoreByteN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMBYTE4::XMBYTE4
(
    const float* pArray
)
{
    XMStoreByte4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUBYTEN4::XMUBYTEN4
(
    const float* pArray
)
{
    XMStoreUByteN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUBYTE4::XMUBYTE4
(
    const float* pArray
)
{
    XMStoreUByte4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMUNIBBLE4::XMUNIBBLE4
(
    const float *pArray
)
{
    XMStoreUNibble4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
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
_Use_decl_annotations_
inline PackedVector::XMU555::XMU555
(
    const float *pArray,
    bool _w
)
{
    XMVECTOR V = XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray));
    XMStoreU555(this, XMVectorSetW(V, ((_w) ? 1.0f : 0.0f) ));
}

///begin_xbox360
/****************************************************************************
 *
 * XMHENDN3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMHENDN3::XMHENDN3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreHenDN3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMHENDN3::XMHENDN3
(
    const float* pArray
)
{
    XMStoreHenDN3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMHEND3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMHEND3::XMHEND3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreHenD3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMHEND3::XMHEND3
(
    const float* pArray
)
{
    XMStoreHenD3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMUHENDN3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMUHENDN3::XMUHENDN3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreUHenDN3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMUHENDN3::XMUHENDN3
(
    const float* pArray
)
{
    XMStoreUHenDN3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMUHEND3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMUHEND3::XMUHEND3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreUHenD3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMUHEND3::XMUHEND3
(
    const float* pArray
)
{
    XMStoreUHenD3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMDHENN3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMDHENN3::XMDHENN3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreDHenN3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMDHENN3::XMDHENN3
(
    const float* pArray
)
{
    XMStoreDHenN3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMDHEN3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMDHEN3::XMDHEN3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreDHen3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMDHEN3::XMDHEN3
(
    const float* pArray
)
{
    XMStoreDHen3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMUDHENN3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMUDHENN3::XMUDHENN3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreUDHenN3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMUDHENN3::XMUDHENN3
(
    const float* pArray
)
{
    XMStoreUDHenN3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMUDHEN3 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMUDHEN3::XMUDHEN3
(
    float _x,
    float _y,
    float _z
)
{
    XMStoreUDHen3(this, XMVectorSet(_x, _y, _z, 0.0f));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMUDHEN3::XMUDHEN3
(
    const float* pArray
)
{
    XMStoreUDHen3(this, XMLoadFloat3(reinterpret_cast<const XMFLOAT3*>(pArray)));
}

/****************************************************************************
 *
 * XMXICON4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMXICON4::XMXICON4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreXIcoN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMXICON4::XMXICON4
(
    const float* pArray
)
{
    XMStoreXIcoN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
}

/****************************************************************************
 *
 * XMXICO4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMXICO4::XMXICO4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreXIco4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMXICO4::XMXICO4
(
    const float* pArray
)
{
    XMStoreXIco4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
}

/****************************************************************************
 *
 * XMICON4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMICON4::XMICON4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreIcoN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMICON4::XMICON4
(
    const float* pArray
)
{
    XMStoreIcoN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
}

/****************************************************************************
 *
 * XMICO4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMICO4::XMICO4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreIco4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMICO4::XMICO4
(
    const float* pArray
)
{
    XMStoreIco4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
}

/****************************************************************************
 *
 * XMUICON4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMUICON4::XMUICON4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUIcoN4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMUICON4::XMUICON4
(
    const float* pArray
)
{
    XMStoreUIcoN4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
}

/****************************************************************************
 *
 * XMUICO4 operators
 *
 ****************************************************************************/

//------------------------------------------------------------------------------

inline PackedVector::Xbox::XMUICO4::XMUICO4
(
    float _x,
    float _y,
    float _z,
    float _w
)
{
    XMStoreUIco4(this, XMVectorSet(_x, _y, _z, _w));
}

//------------------------------------------------------------------------------
_Use_decl_annotations_
inline PackedVector::Xbox::XMUICO4::XMUICO4
(
    const float* pArray
)
{
    XMStoreUIco4(this, XMLoadFloat4(reinterpret_cast<const XMFLOAT4*>(pArray)));
}
///end_xbox360

