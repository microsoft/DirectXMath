//--------------------------------------------------------------------------------------
// File: XDSP.h
//
// DirectXMath based Digital Signal Processing (DSP) functions for audio,
// primarily Fast Fourier Transform (FFT)
//
// All buffer parameters must be 16-byte aligned
//
// All FFT functions support only single-precision floating-point audio
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//
// http://go.microsoft.com/fwlink/?LinkID=615557
//--------------------------------------------------------------------------------------

#pragma once

#include <assert.h>
#include <directxmath.h>

#pragma warning(push)
#pragma warning(disable : 4005 4668)
#include <stdint.h>
#pragma warning(pop)

#pragma warning(push)
#pragma warning(disable: 4328 4640 6001 6262)

namespace XDSP
{
    #if (DIRECTXMATH_VERSION < 305) && !defined(XM_CALLCONV)
    #define XM_CALLCONV __fastcall
    typedef const DirectX::XMVECTOR& HXMVECTOR;
    typedef const DirectX::XMMATRIX& FXMMATRIX;
    #endif

    typedef DirectX::XMVECTOR XMVECTOR;
    typedef DirectX::FXMVECTOR FXMVECTOR;
    typedef DirectX::GXMVECTOR GXMVECTOR;
    typedef DirectX::CXMVECTOR CXMVECTOR;

    inline bool ISPOWEROF2(size_t n) { return ( ((n)&((n)-1)) == 0 && (n) != 0 ); }

    // Parallel multiplication of four complex numbers, assuming real and imaginary values are stored in separate vectors.
    __forceinline void XM_CALLCONV vmulComplex (_Out_ XMVECTOR& rResult, _Out_ XMVECTOR& iResult,
                                                _In_ FXMVECTOR r1, _In_ FXMVECTOR i1, _In_ FXMVECTOR r2, _In_ GXMVECTOR i2)
    {
        using namespace DirectX;
        // (r1, i1) * (r2, i2) = (r1r2 - i1i2, r1i2 + r2i1)
		XMVECTOR vr1r2 = XMVectorMultiply(r1, r2);
		XMVECTOR vr1i2 = XMVectorMultiply(r1, i2);
		rResult = XMVectorNegativeMultiplySubtract(i1, i2, vr1r2); // real: (r1*r2 - i1*i2)
		iResult = XMVectorMultiplyAdd(r2, i1, vr1i2); // imaginary: (r1*i2 + r2*i1)
	}

    __forceinline void XM_CALLCONV vmulComplex (_Inout_ XMVECTOR& r1, _Inout_ XMVECTOR& i1, _In_ FXMVECTOR r2, _In_ FXMVECTOR i2)
    {
        using namespace DirectX;
        // (r1, i1) * (r2, i2) = (r1r2 - i1i2, r1i2 + r2i1)
		XMVECTOR vr1r2 = XMVectorMultiply(r1, r2);
		XMVECTOR vr1i2 = XMVectorMultiply(r1, i2);
		r1 = XMVectorNegativeMultiplySubtract(i1, i2, vr1r2); // real: (r1*r2 - i1*i2)
		i1 = XMVectorMultiplyAdd(r2, i1, vr1i2); // imaginary: (r1*i2 + r2*i1)
	}

    //----------------------------------------------------------------------------------
    // Radix-4 decimation-in-time FFT butterfly.
    // This version assumes that all four elements of the butterfly are
    // adjacent in a single vector.
    //
    // Compute the product of the complex input vector and the
    // 4-element DFT matrix:
    //     | 1  1  1  1 |    | (r1X,i1X) |
    //     | 1 -j -1  j |    | (r1Y,i1Y) |
    //     | 1 -1  1 -1 |    | (r1Z,i1Z) |
    //     | 1  j -1 -j |    | (r1W,i1W) |
    //
    // This matrix can be decomposed into two simpler ones to reduce the
    // number of additions needed. The decomposed matrices look like this:
    //     | 1  0  1  0 |    | 1  0  1  0 |
    //     | 0  1  0 -j |    | 1  0 -1  0 |
    //     | 1  0 -1  0 |    | 0  1  0  1 |
    //     | 0  1  0  j |    | 0  1  0 -1 |
    //
    // Combine as follows:
    //          | 1  0  1  0 |   | (r1X,i1X) |         | (r1X + r1Z, i1X + i1Z) |
    // Temp   = | 1  0 -1  0 | * | (r1Y,i1Y) |       = | (r1X - r1Z, i1X - i1Z) |
    //          | 0  1  0  1 |   | (r1Z,i1Z) |         | (r1Y + r1W, i1Y + i1W) |
    //          | 0  1  0 -1 |   | (r1W,i1W) |         | (r1Y - r1W, i1Y - i1W) |
    //
    //          | 1  0  1  0 |   | (rTempX,iTempX) |   | (rTempX + rTempZ, iTempX + iTempZ) |
    // Result = | 0  1  0 -j | * | (rTempY,iTempY) | = | (rTempY + iTempW, iTempY - rTempW) |
    //          | 1  0 -1  0 |   | (rTempZ,iTempZ) |   | (rTempX - rTempZ, iTempX - iTempZ) |
    //          | 0  1  0  j |   | (rTempW,iTempW) |   | (rTempY - iTempW, iTempY + rTempW) |
    //----------------------------------------------------------------------------------
    __forceinline void ButterflyDIT4_1 (_Inout_ XMVECTOR& r1, _Inout_ XMVECTOR& i1)
    {
        using namespace DirectX;

        // sign constants for radix-4 butterflies
        const static XMVECTORF32 vDFT4SignBits1 = { 1.0f, -1.0f,  1.0f, -1.0f };
        const static XMVECTORF32 vDFT4SignBits2 = { 1.0f,  1.0f, -1.0f, -1.0f };
        const static XMVECTORF32 vDFT4SignBits3 = { 1.0f, -1.0f, -1.0f,  1.0f };

        // calculating Temp
        // [r1X| r1X|r1Y| r1Y] + [r1Z|-r1Z|r1W|-r1W]
        // [i1X| i1X|i1Y| i1Y] + [i1Z|-i1Z|i1W|-i1W]
        XMVECTOR r1L = XMVectorSwizzle<0,0,1,1>( r1 );
        XMVECTOR r1H = XMVectorSwizzle<2,2,3,3>( r1 );

        XMVECTOR i1L = XMVectorSwizzle<0,0,1,1>( i1 );
        XMVECTOR i1H = XMVectorSwizzle<2,2,3,3>( i1 );

        XMVECTOR rTemp = XMVectorMultiplyAdd( r1H, vDFT4SignBits1, r1L );  
        XMVECTOR iTemp = XMVectorMultiplyAdd( i1H, vDFT4SignBits1, i1L ); 

        // calculating Result
        XMVECTOR rZrWiZiW = XMVectorPermute<2,3,6,7>(rTemp,iTemp);  // [rTempZ|rTempW|iTempZ|iTempW]
        XMVECTOR rZiWrZiW = XMVectorSwizzle<0,3,0,3>(rZrWiZiW);     // [rTempZ|iTempW|rTempZ|iTempW]
        XMVECTOR iZrWiZrW = XMVectorSwizzle<2,1,2,1>(rZrWiZiW);     // [rTempZ|iTempW|rTempZ|iTempW]

        // [rTempX| rTempY| rTempX| rTempY] + [rTempZ| iTempW|-rTempZ|-iTempW]
        // [iTempX| iTempY| iTempX| iTempY] + // [iTempZ|-rTempW|-iTempZ| rTempW]
        XMVECTOR rTempL = XMVectorSwizzle<0,1,0,1>(rTemp);
        XMVECTOR iTempL = XMVectorSwizzle<0,1,0,1>(iTemp);

        r1 = XMVectorMultiplyAdd( rZiWrZiW, vDFT4SignBits2, rTempL );
        i1 = XMVectorMultiplyAdd( iZrWiZrW, vDFT4SignBits3, iTempL );                
    }

    //----------------------------------------------------------------------------------
    // Radix-4 decimation-in-time FFT butterfly.
    // This version assumes that elements of the butterfly are
    // in different vectors, so that each vector in the input
    // contains elements from four different butterflies.
    // The four separate butterflies are processed in parallel.
    //
    // The calculations here are the same as the ones in the single-vector
    // radix-4 DFT, but instead of being done on a single vector (X,Y,Z,W)
    // they are done in parallel on sixteen independent complex values.
    // There is no interdependence between the vector elements:
    // | 1  0  1  0 |    | (rIn0,iIn0) |               | (rIn0 + rIn2, iIn0 + iIn2) |
    // | 1  0 -1  0 | *  | (rIn1,iIn1) |  =   Temp   = | (rIn0 - rIn2, iIn0 - iIn2) |
    // | 0  1  0  1 |    | (rIn2,iIn2) |               | (rIn1 + rIn3, iIn1 + iIn3) |
    // | 0  1  0 -1 |    | (rIn3,iIn3) |               | (rIn1 - rIn3, iIn1 - iIn3) |
    //
    //          | 1  0  1  0 |   | (rTemp0,iTemp0) |   | (rTemp0 + rTemp2, iTemp0 + iTemp2) |
    // Result = | 0  1  0 -j | * | (rTemp1,iTemp1) | = | (rTemp1 + iTemp3, iTemp1 - rTemp3) |
    //          | 1  0 -1  0 |   | (rTemp2,iTemp2) |   | (rTemp0 - rTemp2, iTemp0 - iTemp2) |
    //          | 0  1  0  j |   | (rTemp3,iTemp3) |   | (rTemp1 - iTemp3, iTemp1 + rTemp3) |
    //----------------------------------------------------------------------------------
    __forceinline void ButterflyDIT4_4 (_Inout_ XMVECTOR& r0,
                                        _Inout_ XMVECTOR& r1,
                                        _Inout_ XMVECTOR& r2,
                                        _Inout_ XMVECTOR& r3,
                                        _Inout_ XMVECTOR& i0,
                                        _Inout_ XMVECTOR& i1,
                                        _Inout_ XMVECTOR& i2,
                                        _Inout_ XMVECTOR& i3,
                                        _In_reads_(uStride*4) const XMVECTOR* __restrict pUnityTableReal,
                                        _In_reads_(uStride*4) const XMVECTOR* __restrict pUnityTableImaginary,
                                        _In_ size_t uStride,
                                        _In_ const bool fLast)
    {
        using namespace DirectX;

        assert(pUnityTableReal);
        assert(pUnityTableImaginary);
        assert((uintptr_t)pUnityTableReal % 16 == 0);
        assert((uintptr_t)pUnityTableImaginary % 16 == 0);
        assert(ISPOWEROF2(uStride));

        // calculating Temp
        XMVECTOR rTemp0 = XMVectorAdd(r0, r2);
        XMVECTOR iTemp0 = XMVectorAdd(i0, i2);

        XMVECTOR rTemp2 = XMVectorAdd(r1, r3);
        XMVECTOR iTemp2 = XMVectorAdd(i1, i3);

        XMVECTOR rTemp1 = XMVectorSubtract(r0, r2);
        XMVECTOR iTemp1 = XMVectorSubtract(i0, i2);

        XMVECTOR rTemp3 = XMVectorSubtract(r1, r3);
        XMVECTOR iTemp3 = XMVectorSubtract(i1, i3);

        XMVECTOR rTemp4 = XMVectorAdd(rTemp0, rTemp2); 
        XMVECTOR iTemp4 = XMVectorAdd(iTemp0, iTemp2);

        XMVECTOR rTemp5 = XMVectorAdd(rTemp1, iTemp3); 
        XMVECTOR iTemp5 = XMVectorSubtract(iTemp1, rTemp3);

        XMVECTOR rTemp6 = XMVectorSubtract(rTemp0, rTemp2);
        XMVECTOR iTemp6 = XMVectorSubtract(iTemp0, iTemp2);

        XMVECTOR rTemp7 = XMVectorSubtract(rTemp1, iTemp3);
        XMVECTOR iTemp7 = XMVectorAdd(iTemp1, rTemp3);

        // calculating Result
        // vmulComplex(rTemp0, iTemp0, rTemp0, iTemp0, pUnityTableReal[0], pUnityTableImaginary[0]); // first one is always trivial
        vmulComplex(rTemp5, iTemp5, pUnityTableReal[uStride], pUnityTableImaginary[uStride]);
        vmulComplex(rTemp6, iTemp6, pUnityTableReal[uStride*2], pUnityTableImaginary[uStride*2]);
        vmulComplex(rTemp7, iTemp7, pUnityTableReal[uStride*3], pUnityTableImaginary[uStride*3]);
        
        if (fLast)
        {
            ButterflyDIT4_1(rTemp4, iTemp4);
            ButterflyDIT4_1(rTemp5, iTemp5);
            ButterflyDIT4_1(rTemp6, iTemp6);
            ButterflyDIT4_1(rTemp7, iTemp7);
        }

        r0 = rTemp4;    i0 = iTemp4;
        r1 = rTemp5;    i1 = iTemp5;
        r2 = rTemp6;    i2 = iTemp6;
        r3 = rTemp7;    i3 = iTemp7;
    }

    //==================================================================================
    // F-U-N-C-T-I-O-N-S
    //==================================================================================

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  4-sample FFT.
    //
    // PARAMETERS:
    //  pReal      - [inout] real components, must have at least uCount elements
    //  pImaginary - [inout] imaginary components, must have at least uCount elements
    //  uCount     - [in]    number of FFT iterations
    //----------------------------------------------------------------------------------
    __forceinline void FFT4(_Inout_updates_(uCount) XMVECTOR* __restrict pReal,
                            _Inout_updates_(uCount) XMVECTOR* __restrict pImaginary,
                            _In_ const size_t uCount=1)
    {
        assert(pReal);
        assert(pImaginary);
        assert((uintptr_t)pReal % 16 == 0);
        assert((uintptr_t)pImaginary % 16 == 0);
        assert(ISPOWEROF2(uCount));

        for (size_t uIndex=0; uIndex < uCount; ++uIndex)
        {
            ButterflyDIT4_1(pReal[uIndex], pImaginary[uIndex]);
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  8-sample FFT.
    //
    // PARAMETERS:
    //  pReal      - [inout] real components, must have at least uCount*2 elements
    //  pImaginary - [inout] imaginary components, must have at least uCount*2 elements
    //  uCount     - [in]    number of FFT iterations
    //----------------------------------------------------------------------------------
    __forceinline void FFT8 (_Inout_updates_(uCount*2) XMVECTOR* __restrict pReal,
                             _Inout_updates_(uCount*2) XMVECTOR* __restrict pImaginary,
                             _In_ const size_t uCount=1)
    {
        using namespace DirectX;

        assert(pReal);
        assert(pImaginary);
        assert((uintptr_t)pReal % 16 == 0);
        assert((uintptr_t)pImaginary % 16 == 0);
        assert(ISPOWEROF2(uCount));

        static const XMVECTORF32 wr1 = {  1.0f,  0.70710677f,  0.0f, -0.70710677f };
        static const XMVECTORF32 wi1 = {  0.0f, -0.70710677f, -1.0f, -0.70710677f };
        static const XMVECTORF32 wr2 = { -1.0f, -0.70710677f,  0.0f,  0.70710677f };
        static const XMVECTORF32 wi2 = {  0.0f,  0.70710677f,  1.0f,  0.70710677f };

        for (size_t uIndex=0; uIndex < uCount; ++uIndex)
        {
            XMVECTOR* __restrict pR = pReal      + uIndex*2;
            XMVECTOR* __restrict pI = pImaginary + uIndex*2;

            XMVECTOR oddsR  = XMVectorPermute<1,3,5,7>(pR[0], pR[1]);
            XMVECTOR evensR = XMVectorPermute<0,2,4,6>(pR[0], pR[1]);
            XMVECTOR oddsI  = XMVectorPermute<1,3,5,7>(pI[0], pI[1]);
            XMVECTOR evensI = XMVectorPermute<0,2,4,6>(pI[0], pI[1]);
            ButterflyDIT4_1(oddsR, oddsI);
            ButterflyDIT4_1(evensR, evensI);

            XMVECTOR r, i;
            vmulComplex(r, i, oddsR, oddsI, wr1, wi1);
            pR[0] = XMVectorAdd(evensR, r);
            pI[0] = XMVectorAdd(evensI, i);

            vmulComplex(r, i, oddsR, oddsI, wr2, wi2);
            pR[1] = XMVectorAdd(evensR, r);
            pI[1] = XMVectorAdd(evensI, i);
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  16-sample FFT.
    //
    // PARAMETERS:
    //  pReal      - [inout] real components, must have at least uCount*4 elements
    //  pImaginary - [inout] imaginary components, must have at least uCount*4 elements
    //  uCount     - [in]    number of FFT iterations
    //----------------------------------------------------------------------------------
    __forceinline void FFT16 (_Inout_updates_(uCount*4) XMVECTOR* __restrict pReal,
                              _Inout_updates_(uCount*4) XMVECTOR* __restrict pImaginary,
                              _In_ const size_t uCount=1)
    {
        using namespace DirectX;

        assert(pReal);
        assert(pImaginary);
        assert((uintptr_t)pReal % 16 == 0);
        assert((uintptr_t)pImaginary % 16 == 0);
        assert(ISPOWEROF2(uCount));

        static const XMVECTORF32 aUnityTableReal[4]      = { { 1.0f, 1.0f, 1.0f, 1.0f },
                                                             { 1.0f, 0.92387950f, 0.70710677f, 0.38268343f },
                                                             { 1.0f, 0.70710677f, -4.3711388e-008f, -0.70710677f },
                                                             { 1.0f, 0.38268343f, -0.70710677f, -0.92387950f } };
        static const XMVECTORF32 aUnityTableImaginary[4] = { { -0.0f, -0.0f, -0.0f, -0.0f },
                                                             { -0.0f, -0.38268343f, -0.70710677f, -0.92387950f },
                                                             { -0.0f, -0.70710677f, -1.0f, -0.70710677f },
                                                             { -0.0f, -0.92387950f, -0.70710677f, 0.38268343f } };

        for (size_t uIndex=0; uIndex < uCount; ++uIndex)
        {
            ButterflyDIT4_4(pReal[uIndex*4],
                            pReal[uIndex*4 + 1],
                            pReal[uIndex*4 + 2],
                            pReal[uIndex*4 + 3],
                            pImaginary[uIndex*4],
                            pImaginary[uIndex*4 + 1],
                            pImaginary[uIndex*4 + 2],
                            pImaginary[uIndex*4 + 3],
                            reinterpret_cast<const XMVECTOR*>(aUnityTableReal),
                            reinterpret_cast<const XMVECTOR*>(aUnityTableImaginary),
                            1, true);
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  2^N-sample FFT.
    //
    // REMARKS:
    //  For FFTs length 16 and below, call FFT16(), FFT8(), or FFT4().
    //
    // PARAMETERS:
    //  pReal       - [inout] real components, must have at least (uLength*uCount)/4 elements
    //  pImaginary  - [inout] imaginary components, must have at least (uLength*uCount)/4 elements
    //  pUnityTable - [in]    unity table, must have at least uLength*uCount elements, see FFTInitializeUnityTable()
    //  uLength     - [in]    FFT length in samples, must be a power of 2 > 16
    //  uCount      - [in]    number of FFT iterations
    //----------------------------------------------------------------------------------
    inline void FFT (_Inout_updates_((uLength*uCount)/4) XMVECTOR* __restrict pReal,
                     _Inout_updates_((uLength*uCount)/4) XMVECTOR* __restrict pImaginary,
                     _In_reads_(uLength*uCount) const XMVECTOR* __restrict pUnityTable,
                     _In_ const size_t uLength,
                     _In_ const size_t uCount=1)
    {
        assert(pReal);
        assert(pImaginary);
        assert(pUnityTable);
        assert((uintptr_t)pReal % 16 == 0);
        assert((uintptr_t)pImaginary % 16 == 0);
        assert((uintptr_t)pUnityTable % 16 == 0);
        assert(uLength > 16);
        _Analysis_assume_(uLength > 16);
        assert(ISPOWEROF2(uLength));
        assert(ISPOWEROF2(uCount));

        const XMVECTOR* __restrict pUnityTableReal      = pUnityTable;
        const XMVECTOR* __restrict pUnityTableImaginary = pUnityTable + (uLength>>2);
        const size_t uTotal              = uCount * uLength;
        const size_t uTotal_vectors      = uTotal >> 2;
        const size_t uStage_vectors      = uLength >> 2;
        const size_t uStage_vectors_mask = uStage_vectors - 1;
        const size_t uStride        = uLength >> 4; // stride between butterfly elements
        const size_t uStrideMask    = uStride - 1;
        const size_t uStride2       = uStride * 2;
        const size_t uStride3       = uStride * 3;
        const size_t uStrideInvMask = ~uStrideMask;

        for (size_t uIndex=0; uIndex < (uTotal_vectors>>2); ++uIndex)
        {
            const size_t n = ((uIndex & uStrideInvMask) << 2) + (uIndex & uStrideMask);
            ButterflyDIT4_4(pReal[n],
                            pReal[n + uStride],
                            pReal[n + uStride2],
                            pReal[n + uStride3],
                            pImaginary[n ],
                            pImaginary[n + uStride],
                            pImaginary[n + uStride2],
                            pImaginary[n + uStride3],
                            pUnityTableReal      + (n & uStage_vectors_mask),
                            pUnityTableImaginary + (n & uStage_vectors_mask),
                            uStride, false);
        }

        if (uLength > 16*4)
        {
            FFT(pReal, pImaginary, pUnityTable+(uLength>>1), uLength>>2, uCount*4);
        }
        else if (uLength == 16*4)
        {
            FFT16(pReal, pImaginary, uCount*4);
        }
        else if (uLength == 8*4)
        {
            FFT8(pReal, pImaginary, uCount*4);
        }
        else if (uLength == 4*4)
        {
            FFT4(pReal, pImaginary, uCount*4);
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  Initializes unity roots lookup table used by FFT functions.
    //  Once initialized, the table need not be initialized again unless a
    //  different FFT length is desired.
    //
    // REMARKS:
    //  The unity tables of FFT length 16 and below are hard coded into the
    //  respective FFT functions and so need not be initialized.
    //
    // PARAMETERS:
    //  pUnityTable - [out] unity table, receives unity roots lookup table, must have at least uLength elements
    //  uLength     - [in]  FFT length in frames, must be a power of 2 > 16
    //----------------------------------------------------------------------------------
    inline void FFTInitializeUnityTable (_Out_writes_(uLength) XMVECTOR* __restrict pUnityTable, _In_ size_t uLength)
    {
        assert(pUnityTable);
        assert(uLength > 16);
        _Analysis_assume_(uLength > 16);
        assert(ISPOWEROF2(uLength));

        float* __restrict pfUnityTable = reinterpret_cast<float* __restrict>(pUnityTable);

        // initialize unity table for recursive FFT lengths: uLength, uLength/4, uLength/16... > 16
        do
        {
            float flStep = 6.283185307f / uLength; // 2PI / FFT length
            uLength >>= 2;

            // pUnityTable[0 to uLength*4-1] contains real components for current FFT length
            // pUnityTable[uLength*4 to uLength*8-1] contains imaginary components for current FFT length
            for (size_t i=0; i<4; ++i)
            {
                for (size_t j=0; j<uLength; ++j)
                {
                    size_t uIndex = (i*uLength) + j;
                    pfUnityTable[uIndex]             = cosf(float(i)*float(j)*flStep);  // real component
#pragma warning(suppress: 6386)
                    pfUnityTable[uIndex + uLength*4] = -sinf(float(i)*float(j)*flStep); // imaginary component
                }
            }
            pfUnityTable += uLength*8;
        }
        while (uLength > 16);
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  The FFT functions generate output in bit reversed order.
    //  Use this function to re-arrange them into order of increasing frequency.
    //
    // REMARKS:
    //
    // PARAMETERS:
    //  pOutput     - [out] output buffer, receives samples in order of increasing frequency, cannot overlap pInput, must have at least (1<<uLog2Length)/4 elements
    //  pInput      - [in]  input buffer, samples in bit reversed order as generated by FFT functions, cannot overlap pOutput, must have at least (1<<uLog2Length)/4 elements
    //  uLog2Length - [in]  LOG (base 2) of FFT length in samples, must be >= 2
    //----------------------------------------------------------------------------------
    inline void FFTUnswizzle (_Out_writes_((1<<uLog2Length)/4) XMVECTOR* __restrict pOutput,
                              _In_reads_((1<<uLog2Length)/4) const XMVECTOR* __restrict pInput,
                              _In_ const size_t uLog2Length)
    {
        assert(pOutput);
        assert(pInput);
        assert(uLog2Length >= 2);
        _Analysis_assume_(uLog2Length >= 2);

        float*       __restrict pfOutput = (float* __restrict)pOutput;
        const float* __restrict pfInput  = (const float* __restrict)pInput;
        const size_t uLength = size_t(1) << uLog2Length;

        if ((uLog2Length & 0x1) == 0)
        {
            // even powers of two
            for (size_t uIndex=0; uIndex < uLength; ++uIndex)
            {
                size_t n = uIndex;
                n = ( (n & 0xcccccccc) >> 2 )  | ( (n & 0x33333333) << 2 );
                n = ( (n & 0xf0f0f0f0) >> 4 )  | ( (n & 0x0f0f0f0f) << 4 );
                n = ( (n & 0xff00ff00) >> 8 )  | ( (n & 0x00ff00ff) << 8 );
                n = ( (n & 0xffff0000) >> 16 ) | ( (n & 0x0000ffff) << 16 );
                n >>= (32 - uLog2Length);
                pfOutput[n] = pfInput[uIndex];
            }
        }
        else
        {
            // odd powers of two
            for (size_t uIndex=0; uIndex < uLength; ++uIndex)
            {
                size_t n = (uIndex>>3);
                n = ( (n & 0xcccccccc) >> 2 )  | ( (n & 0x33333333) << 2 );
                n = ( (n & 0xf0f0f0f0) >> 4 )  | ( (n & 0x0f0f0f0f) << 4 );
                n = ( (n & 0xff00ff00) >> 8 )  | ( (n & 0x00ff00ff) << 8 );
                n = ( (n & 0xffff0000) >> 16 ) | ( (n & 0x0000ffff) << 16 );
                n >>= (32 - (uLog2Length-3));
                n |= ((uIndex & 0x7) << (uLog2Length - 3));
                pfOutput[n] = pfInput[uIndex];
            }
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  Convert complex components to polar form.
    //
    // PARAMETERS:
    //  pOutput         - [out] output buffer, receives samples in polar form, must have at least uLength/4 elements
    //  pInputReal      - [in]  input buffer (real components), must have at least uLength/4 elements
    //  pInputImaginary - [in]  input buffer (imaginary components), must have at least uLength/4 elements
    //  uLength         - [in]  FFT length in samples, must be a power of 2 >= 4
    //----------------------------------------------------------------------------------
#pragma warning(suppress: 6101)
    inline void FFTPolar (_Out_writes_(uLength/4) XMVECTOR* __restrict pOutput,
                          _In_reads_(uLength/4) const XMVECTOR* __restrict pInputReal,
                          _In_reads_(uLength/4) const XMVECTOR* __restrict pInputImaginary,
                          _In_ const size_t uLength)
    {
        using namespace DirectX;

        assert(pOutput);
        assert(pInputReal);
        assert(pInputImaginary);
        assert(uLength >= 4);
        _Analysis_assume_(uLength >= 4);
        assert(ISPOWEROF2(uLength));

        float flOneOverLength = 1.0f / uLength;

        // result = sqrtf((real/uLength)^2 + (imaginary/uLength)^2) * 2
        XMVECTOR vOneOverLength = XMVectorReplicate( flOneOverLength );

        for (size_t uIndex=0; uIndex < (uLength>>2); ++uIndex)
        {
            XMVECTOR vReal      = XMVectorMultiply(pInputReal[uIndex], vOneOverLength);
            XMVECTOR vImaginary = XMVectorMultiply(pInputImaginary[uIndex], vOneOverLength);
            XMVECTOR vRR        = XMVectorMultiply(vReal, vReal);
            XMVECTOR vII        = XMVectorMultiply(vImaginary, vImaginary);
            XMVECTOR vRRplusII  = XMVectorAdd(vRR, vII);
            XMVECTOR vTotal     = XMVectorSqrt(vRRplusII);
            pOutput[uIndex]     = XMVectorAdd(vTotal, vTotal);
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  Deinterleaves audio samples
    //
    // REMARKS:
    //  For example, audio of the form [LRLRLR] becomes [LLLRRR].
    //
    // PARAMETERS:
    //  pOutput       - [out] output buffer, receives samples in deinterleaved form, cannot overlap pInput, must have at least (uChannelCount*uFrameCount)/4 elements
    //  pInput        - [in]  input buffer, cannot overlap pOutput, must have at least (uChannelCount*uFrameCount)/4 elements
    //  uChannelCount - [in]  number of channels, must be > 1
    //  uFrameCount   - [in]  number of frames of valid data, must be > 0
    //----------------------------------------------------------------------------------
    inline void Deinterleave (_Out_writes_((uChannelCount*uFrameCount)/4) XMVECTOR* __restrict pOutput,
                              _In_reads_((uChannelCount*uFrameCount)/4) const XMVECTOR* __restrict pInput,
                              _In_ const size_t uChannelCount,
                              _In_ const size_t uFrameCount)
    {
        assert(pOutput);
        assert(pInput);
        assert(uChannelCount > 1);
        assert(uFrameCount > 0);

        float*       __restrict pfOutput = reinterpret_cast<float* __restrict>(pOutput);
        const float* __restrict pfInput  = reinterpret_cast<const float* __restrict>(pInput);

        for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
        {
            for (size_t uFrame=0; uFrame < uFrameCount; ++uFrame)
            {
                pfOutput[uChannel * uFrameCount + uFrame] = pfInput[uFrame * uChannelCount + uChannel];
            }
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  Interleaves audio samples
    //
    // REMARKS:
    //  For example, audio of the form [LLLRRR] becomes [LRLRLR].
    //
    // PARAMETERS:
    //  pOutput       - [out] output buffer, receives samples in interleaved form, cannot overlap pInput, must have at least (uChannelCount*uFrameCount)/4 elements
    //  pInput        - [in]  input buffer, cannot overlap pOutput, must have at least (uChannelCount*uFrameCount)/4 elements
    //  uChannelCount - [in]  number of channels, must be > 1
    //  uFrameCount   - [in]  number of frames of valid data, must be > 0
    //----------------------------------------------------------------------------------
    inline void Interleave (_Out_writes_((uChannelCount*uFrameCount)/4) XMVECTOR* __restrict pOutput,
                            _In_reads_((uChannelCount*uFrameCount)/4) const XMVECTOR* __restrict pInput,
                            _In_ const size_t uChannelCount,
                            _In_ const size_t uFrameCount)
    {
        assert(pOutput);
        assert(pInput);
        assert(uChannelCount > 1);
        assert(uFrameCount > 0);

        float*       __restrict pfOutput = reinterpret_cast<float* __restrict>(pOutput);
        const float* __restrict pfInput  = reinterpret_cast<const float* __restrict>(pInput);

        for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
        {
            for (size_t uFrame=0; uFrame < uFrameCount; ++uFrame)
            {
                pfOutput[uFrame * uChannelCount + uChannel] = pfInput[uChannel * uFrameCount + uFrame];
            }
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  This function applies a 2^N-sample FFT and unswizzles the result such
    //  that the samples are in order of increasing frequency.
    //  Audio is first deinterleaved if multichannel.
    //
    // PARAMETERS:
    //  pReal         - [inout] real components, must have at least (1<<uLog2Length*uChannelCount)/4 elements
    //  pImaginary    - [out]   imaginary components, must have at least (1<<uLog2Length*uChannelCount)/4 elements
    //  pUnityTable   - [in]    unity table, must have at least (1<<uLog2Length) elements, see FFTInitializeUnityTable()
    //  uChannelCount - [in]    number of channels, must be within [1, 6]
    //  uLog2Length   - [in]    LOG (base 2) of FFT length in frames, must within [2, 9]
    //----------------------------------------------------------------------------------
    inline void FFTInterleaved (_Inout_updates_(((1<<uLog2Length)*uChannelCount)/4) XMVECTOR* __restrict pReal,
                                _Out_writes_(((1<<uLog2Length)*uChannelCount)/4) XMVECTOR* __restrict pImaginary,
                                _In_reads_(1<<uLog2Length) const XMVECTOR* __restrict pUnityTable,
                                _In_ const size_t uChannelCount,
                                _In_ const size_t uLog2Length)
    {
        assert(pReal);
        assert(pImaginary);
        assert(pUnityTable);
        assert((uintptr_t)pReal % 16 == 0);
        assert((uintptr_t)pImaginary % 16 == 0);
        assert((uintptr_t)pUnityTable % 16 == 0);
        assert(uChannelCount > 0 && uChannelCount <= 6);
        assert(uLog2Length >= 2 && uLog2Length <= 9);

        XMVECTOR vRealTemp[768];
        XMVECTOR vImaginaryTemp[768];
        const size_t uLength = size_t(1) << uLog2Length;

        if (uChannelCount > 1)
        {
            Deinterleave(vRealTemp, pReal, uChannelCount, uLength);
        }
        else
        {
            memcpy_s(vRealTemp, sizeof(vRealTemp), pReal, (uLength>>2)*sizeof(XMVECTOR));
        }

        memset( vImaginaryTemp, 0, (uChannelCount*(uLength>>2)) * sizeof(XMVECTOR) );

        if (uLength > 16)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)], pUnityTable, uLength);
            }
        }
        else if (uLength == 16)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT16(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)]);
            }
        }
        else if (uLength == 8)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT8(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)]);
            }
        }
        else if (uLength == 4)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT4(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)]);
            }
        }

        for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
        {
            FFTUnswizzle(&pReal[uChannel*(uLength>>2)], &vRealTemp[uChannel*(uLength>>2)], uLog2Length);
            FFTUnswizzle(&pImaginary[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)], uLog2Length);
        }
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  This function applies a 2^N-sample inverse FFT.
    //  Audio is interleaved if multichannel.
    //
    // PARAMETERS:
    //  pReal         - [inout] real components, must have at least (1<<uLog2Length*uChannelCount)/4 elements
    //  pImaginary    - [in]    imaginary components, must have at least (1<<uLog2Length*uChannelCount)/4 elements
    //  pUnityTable   - [in]    unity table, must have at least (1<<uLog2Length) elements, see FFTInitializeUnityTable()
    //  uChannelCount - [in]    number of channels, must be > 0
    //  uLog2Length   - [in]    LOG (base 2) of FFT length in frames, must within [2, 9]
    //----------------------------------------------------------------------------------
    inline void IFFTDeinterleaved (_Inout_updates_(((1<<uLog2Length)*uChannelCount)/4) XMVECTOR* __restrict pReal,
                                   _In_reads_(((1<<uLog2Length)*uChannelCount)/4) const XMVECTOR* __restrict pImaginary,
                                   _In_reads_(1<<uLog2Length) const XMVECTOR* __restrict pUnityTable,
                                   _In_ const size_t uChannelCount,
                                   _In_ const size_t uLog2Length)
    {
        using namespace DirectX;

        assert(pReal);
        assert(pImaginary);
        assert(pUnityTable);
        assert((uintptr_t)pReal % 16 == 0);
        assert((uintptr_t)pImaginary % 16 == 0);
        assert((uintptr_t)pUnityTable % 16 == 0);
        assert(uChannelCount > 0 && uChannelCount <= 6);
        _Analysis_assume_(uChannelCount > 0 && uChannelCount <= 6);
        assert(uLog2Length >= 2 && uLog2Length <= 9);
        _Analysis_assume_(uLog2Length >= 2 && uLog2Length <= 9);

        XMVECTOR vRealTemp[768] = { 0 };
        XMVECTOR vImaginaryTemp[768] = { 0 };

        const size_t uLength = size_t(1) << uLog2Length;

        const XMVECTOR vRnp = XMVectorReplicate(1.0f/uLength);
        const XMVECTOR vRnm = XMVectorReplicate(-1.0f/uLength);
        for (size_t u=0; u < uChannelCount*(uLength>>2); u++)
        {
            vRealTemp[u]      = XMVectorMultiply(pReal[u], vRnp);
            vImaginaryTemp[u] = XMVectorMultiply(pImaginary[u], vRnm);
        }

        if (uLength > 16)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)], pUnityTable, uLength);
            }
        }
        else if (uLength == 16)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT16(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)]);
            }
        }
        else if (uLength == 8)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT8(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)]);
            }
        }
        else if (uLength == 4)
        {
            for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
            {
                FFT4(&vRealTemp[uChannel*(uLength>>2)], &vImaginaryTemp[uChannel*(uLength>>2)]);
            }
        }

        for (size_t uChannel=0; uChannel < uChannelCount; ++uChannel)
        {
            FFTUnswizzle(&vImaginaryTemp[uChannel*(uLength>>2)], &vRealTemp[uChannel*(uLength>>2)], uLog2Length);
        }

        if (uChannelCount > 1)
        {
            Interleave(pReal, vImaginaryTemp, uChannelCount, uLength);
        }
        else
        {
            memcpy_s(pReal, uLength*uChannelCount*sizeof(float), vImaginaryTemp, (uLength>>2)*sizeof(XMVECTOR));
        }
    }

}; // namespace XDSP

#pragma warning(pop)
