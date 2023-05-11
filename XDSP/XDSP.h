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
// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
//
// http://go.microsoft.com/fwlink/?LinkID=615557
//--------------------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <DirectXMath.h>

#include <cstdint>
#include <cstring>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 6001 6262)
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

namespace XDSP
{
    using XMVECTOR = DirectX::XMVECTOR;
    using FXMVECTOR = DirectX::FXMVECTOR;
    using GXMVECTOR = DirectX::GXMVECTOR;
    using CXMVECTOR = DirectX::CXMVECTOR;
    using XMFLOAT4A = DirectX::XMFLOAT4A;

    inline bool ISPOWEROF2(size_t n) { return (((n)&((n)-1)) == 0 && (n) != 0); }

    // Parallel multiplication of four complex numbers, assuming real and imaginary values are stored in separate vectors.
    inline void XM_CALLCONV vmulComplex(
        _Out_ XMVECTOR& rResult, _Out_ XMVECTOR& iResult,
        _In_ FXMVECTOR r1, _In_ FXMVECTOR i1, _In_ FXMVECTOR r2, _In_ GXMVECTOR i2) noexcept
    {
        using namespace DirectX;
        // (r1, i1) * (r2, i2) = (r1r2 - i1i2, r1i2 + r2i1)
        const XMVECTOR vr1r2 = XMVectorMultiply(r1, r2);
        const XMVECTOR vr1i2 = XMVectorMultiply(r1, i2);
        rResult = XMVectorNegativeMultiplySubtract(i1, i2, vr1r2); // real: (r1*r2 - i1*i2)
        iResult = XMVectorMultiplyAdd(r2, i1, vr1i2); // imaginary: (r1*i2 + r2*i1)
    }

    inline void XM_CALLCONV vmulComplex(
        _Inout_ XMVECTOR& r1, _Inout_ XMVECTOR& i1, _In_ FXMVECTOR r2, _In_ FXMVECTOR i2) noexcept
    {
        using namespace DirectX;
        // (r1, i1) * (r2, i2) = (r1r2 - i1i2, r1i2 + r2i1)
        const XMVECTOR vr1r2 = XMVectorMultiply(r1, r2);
        const XMVECTOR vr1i2 = XMVectorMultiply(r1, i2);
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
    inline void ButterflyDIT4_1 (_Inout_ XMVECTOR& r1, _Inout_ XMVECTOR& i1) noexcept
    {
        using namespace DirectX;

        // sign constants for radix-4 butterflies
        static const XMVECTORF32 vDFT4SignBits1 = { { { 1.0f, -1.0f, 1.0f, -1.0f } } };
        static const XMVECTORF32 vDFT4SignBits2 = { { { 1.0f, 1.0f, -1.0f, -1.0f } } };
        static const XMVECTORF32 vDFT4SignBits3 = { { { 1.0f, -1.0f, -1.0f, 1.0f } } };

        // calculating Temp
        // [r1X| r1X|r1Y| r1Y] + [r1Z|-r1Z|r1W|-r1W]
        // [i1X| i1X|i1Y| i1Y] + [i1Z|-i1Z|i1W|-i1W]
        const XMVECTOR r1L = XMVectorSwizzle<0, 0, 1, 1>(r1);
        const XMVECTOR r1H = XMVectorSwizzle<2, 2, 3, 3>(r1);

        const XMVECTOR i1L = XMVectorSwizzle<0, 0, 1, 1>(i1);
        const XMVECTOR i1H = XMVectorSwizzle<2, 2, 3, 3>(i1);

        const XMVECTOR rTemp = XMVectorMultiplyAdd(r1H, vDFT4SignBits1, r1L);
        const XMVECTOR iTemp = XMVectorMultiplyAdd(i1H, vDFT4SignBits1, i1L);

        // calculating Result
        const XMVECTOR rZrWiZiW = XMVectorPermute<2, 3, 6, 7>(rTemp, iTemp);   // [rTempZ|rTempW|iTempZ|iTempW]
        const XMVECTOR rZiWrZiW = XMVectorSwizzle<0, 3, 0, 3>(rZrWiZiW);       // [rTempZ|iTempW|rTempZ|iTempW]
        const XMVECTOR iZrWiZrW = XMVectorSwizzle<2, 1, 2, 1>(rZrWiZiW);       // [rTempZ|iTempW|rTempZ|iTempW]

        // [rTempX| rTempY| rTempX| rTempY] + [rTempZ| iTempW|-rTempZ|-iTempW]
        // [iTempX| iTempY| iTempX| iTempY] + // [iTempZ|-rTempW|-iTempZ| rTempW]
        const XMVECTOR rTempL = XMVectorSwizzle<0, 1, 0, 1>(rTemp);
        const XMVECTOR iTempL = XMVectorSwizzle<0, 1, 0, 1>(iTemp);

        r1 = XMVectorMultiplyAdd(rZiWrZiW, vDFT4SignBits2, rTempL);
        i1 = XMVectorMultiplyAdd(iZrWiZrW, vDFT4SignBits3, iTempL);
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
    inline void ButterflyDIT4_4(
        _Inout_ XMVECTOR& r0,
        _Inout_ XMVECTOR& r1,
        _Inout_ XMVECTOR& r2,
        _Inout_ XMVECTOR& r3,
        _Inout_ XMVECTOR& i0,
        _Inout_ XMVECTOR& i1,
        _Inout_ XMVECTOR& i2,
        _Inout_ XMVECTOR& i3,
        _In_reads_(uStride * 4) const XMVECTOR* __restrict pUnityTableReal,
        _In_reads_(uStride * 4) const XMVECTOR* __restrict pUnityTableImaginary,
        _In_ size_t uStride,
        _In_ const bool fLast) noexcept
    {
        using namespace DirectX;

        assert(pUnityTableReal);
        assert(pUnityTableImaginary);
        assert(reinterpret_cast<uintptr_t>(pUnityTableReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pUnityTableImaginary) % 16 == 0);
        assert(ISPOWEROF2(uStride));

        // calculating Temp
        const XMVECTOR rTemp0 = XMVectorAdd(r0, r2);
        const XMVECTOR iTemp0 = XMVectorAdd(i0, i2);

        const XMVECTOR rTemp2 = XMVectorAdd(r1, r3);
        const XMVECTOR iTemp2 = XMVectorAdd(i1, i3);

        const XMVECTOR rTemp1 = XMVectorSubtract(r0, r2);
        const XMVECTOR iTemp1 = XMVectorSubtract(i0, i2);

        const XMVECTOR rTemp3 = XMVectorSubtract(r1, r3);
        const XMVECTOR iTemp3 = XMVectorSubtract(i1, i3);

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
        vmulComplex(rTemp6, iTemp6, pUnityTableReal[uStride * 2], pUnityTableImaginary[uStride * 2]);
        vmulComplex(rTemp7, iTemp7, pUnityTableReal[uStride * 3], pUnityTableImaginary[uStride * 3]);

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
    inline void FFT4(
        _Inout_updates_(uCount) XMVECTOR* __restrict pReal,
        _Inout_updates_(uCount) XMVECTOR* __restrict pImaginary,
        const size_t uCount = 1) noexcept
    {
        assert(pReal);
        assert(pImaginary);
        assert(reinterpret_cast<uintptr_t>(pReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pImaginary) % 16 == 0);
        assert(ISPOWEROF2(uCount));

        for (size_t uIndex = 0; uIndex < uCount; ++uIndex)
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
    inline void FFT8(
        _Inout_updates_(uCount * 2) XMVECTOR* __restrict pReal,
        _Inout_updates_(uCount * 2) XMVECTOR* __restrict pImaginary,
        _In_ const size_t uCount = 1) noexcept
    {
        using namespace DirectX;

        assert(pReal);
        assert(pImaginary);
        assert(reinterpret_cast<uintptr_t>(pReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pImaginary) % 16 == 0);
        assert(ISPOWEROF2(uCount));

        static const XMVECTORF32 wr1 = { { { 1.0f, 0.70710677f, 0.0f, -0.70710677f } } };
        static const XMVECTORF32 wi1 = { { { 0.0f, -0.70710677f, -1.0f, -0.70710677f } } };
        static const XMVECTORF32 wr2 = { { { -1.0f, -0.70710677f, 0.0f, 0.70710677f } } };
        static const XMVECTORF32 wi2 = { { { 0.0f, 0.70710677f, 1.0f, 0.70710677f } } };

        for (size_t uIndex = 0; uIndex < uCount; ++uIndex)
        {
            XMVECTOR* __restrict pR = pReal + uIndex * 2;
            XMVECTOR* __restrict pI = pImaginary + uIndex * 2;

            XMVECTOR oddsR = XMVectorPermute<1, 3, 5, 7>(pR[0], pR[1]);
            XMVECTOR evensR = XMVectorPermute<0, 2, 4, 6>(pR[0], pR[1]);
            XMVECTOR oddsI = XMVectorPermute<1, 3, 5, 7>(pI[0], pI[1]);
            XMVECTOR evensI = XMVectorPermute<0, 2, 4, 6>(pI[0], pI[1]);
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
    inline void FFT16(
        _Inout_updates_(uCount * 4) XMVECTOR* __restrict pReal,
        _Inout_updates_(uCount * 4) XMVECTOR* __restrict pImaginary,
        _In_ const size_t uCount = 1) noexcept
    {
        using namespace DirectX;

        assert(pReal);
        assert(pImaginary);
        assert(reinterpret_cast<uintptr_t>(pReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pImaginary) % 16 == 0);
        assert(ISPOWEROF2(uCount));

        static const XMVECTORF32 aUnityTableReal[4] = {
            { { { 1.0f, 1.0f, 1.0f, 1.0f } } },
            { { { 1.0f, 0.92387950f, 0.70710677f, 0.38268343f } } },
            { { { 1.0f, 0.70710677f, -4.3711388e-008f, -0.70710677f } } },
            { { { 1.0f, 0.38268343f, -0.70710677f, -0.92387950f } } }
        };
        static const XMVECTORF32 aUnityTableImaginary[4] =
        {
            { { { -0.0f, -0.0f, -0.0f, -0.0f } } },
            { { { -0.0f, -0.38268343f, -0.70710677f, -0.92387950f } } },
            { { { -0.0f, -0.70710677f, -1.0f, -0.70710677f } } },
            { { { -0.0f, -0.92387950f, -0.70710677f, 0.38268343f } } }
        };

        for (size_t uIndex = 0; uIndex < uCount; ++uIndex)
        {
            ButterflyDIT4_4(pReal[uIndex * 4],
                pReal[uIndex * 4 + 1],
                pReal[uIndex * 4 + 2],
                pReal[uIndex * 4 + 3],
                pImaginary[uIndex * 4],
                pImaginary[uIndex * 4 + 1],
                pImaginary[uIndex * 4 + 2],
                pImaginary[uIndex * 4 + 3],
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
    inline void FFT (
        _Inout_updates_((uLength * uCount) / 4) XMVECTOR* __restrict pReal,
        _Inout_updates_((uLength * uCount) / 4) XMVECTOR* __restrict pImaginary,
        _In_reads_(uLength * uCount) const XMVECTOR* __restrict pUnityTable,
        _In_ const size_t uLength,
        _In_ const size_t uCount = 1) noexcept
    {
        assert(pReal);
        assert(pImaginary);
        assert(pUnityTable);
        assert(reinterpret_cast<uintptr_t>(pReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pImaginary) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pUnityTable) % 16 == 0);
        assert(uLength > 16);
        _Analysis_assume_(uLength > 16);
        assert(ISPOWEROF2(uLength));
        assert(ISPOWEROF2(uCount));

        const XMVECTOR* __restrict pUnityTableReal = pUnityTable;
        const XMVECTOR* __restrict pUnityTableImaginary = pUnityTable + (uLength >> 2);
        const size_t uTotal              = uCount * uLength;
        const size_t uTotal_vectors      = uTotal >> 2;
        const size_t uStage_vectors      = uLength >> 2;
        const size_t uStage_vectors_mask = uStage_vectors - 1;
        const size_t uStride        = uLength >> 4; // stride between butterfly elements
        const size_t uStrideMask    = uStride - 1;
        const size_t uStride2       = uStride * 2;
        const size_t uStride3       = uStride * 3;
        const size_t uStrideInvMask = ~uStrideMask;

        for (size_t uIndex=0; uIndex < (uTotal_vectors >> 2); ++uIndex)
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
                            pUnityTableReal + (n & uStage_vectors_mask),
                            pUnityTableImaginary + (n & uStage_vectors_mask),
                            uStride, false);
        }

        if (uLength > 16 * 4)
        {
            FFT(pReal, pImaginary, pUnityTable + (uLength >> 1), uLength >> 2, uCount * 4);
        }
        else if (uLength == 16 * 4)
        {
            FFT16(pReal, pImaginary, uCount * 4);
        }
        else if (uLength == 8 * 4)
        {
            FFT8(pReal, pImaginary, uCount * 4);
        }
        else if (uLength == 4 * 4)
        {
            FFT4(pReal, pImaginary, uCount * 4);
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
    inline void FFTInitializeUnityTable (_Out_writes_(uLength) XMVECTOR* __restrict pUnityTable, _In_ size_t uLength) noexcept
    {
        using namespace DirectX;

        assert(pUnityTable);
        assert(uLength > 16);
        _Analysis_assume_(uLength > 16);
        assert(ISPOWEROF2(uLength));

        // initialize unity table for recursive FFT lengths: uLength, uLength/4, uLength/16... > 16
        // pUnityTable[0 to uLength*4-1] contains real components for current FFT length
        // pUnityTable[uLength*4 to uLength*8-1] contains imaginary components for current FFT length
        static const XMVECTORF32 vXM0123 = { { { 0.0f, 1.0f, 2.0f, 3.0f } } };
        uLength >>= 2;
        XMVECTOR vlStep = XMVectorReplicate(XM_PIDIV2 / float(uLength));
        do
        {
            uLength >>= 2;
            XMVECTOR vJP = vXM0123;
            for (size_t j = 0; j < uLength; ++j)
            {
                XMVECTOR vSin, vCos;
                XMVECTOR viJP, vlS;

                pUnityTable[j] = g_XMOne;
                pUnityTable[j + uLength * 4] = XMVectorZero();

                vlS = XMVectorMultiply(vJP, vlStep);
                XMVectorSinCos(&vSin, &vCos, vlS);
                pUnityTable[j + uLength] = vCos;
                pUnityTable[j + uLength * 5] = XMVectorMultiply(vSin, g_XMNegativeOne);

                viJP = XMVectorAdd(vJP, vJP);
                vlS = XMVectorMultiply(viJP, vlStep);
                XMVectorSinCos(&vSin, &vCos, vlS);
                pUnityTable[j + uLength * 2] = vCos;
                pUnityTable[j + uLength * 6] = XMVectorMultiply(vSin, g_XMNegativeOne);

                viJP = XMVectorAdd(viJP, vJP);
                vlS = XMVectorMultiply(viJP, vlStep);
                XMVectorSinCos(&vSin, &vCos, vlS);
                pUnityTable[j + uLength * 3] = vCos;
                pUnityTable[j + uLength * 7] = XMVectorMultiply(vSin, g_XMNegativeOne);

                vJP = XMVectorAdd(vJP, g_XMFour);
            }
            vlStep = XMVectorMultiply(vlStep, g_XMFour);
            pUnityTable += uLength * 8;
        } while (uLength > 4);
    }

    //----------------------------------------------------------------------------------
    // DESCRIPTION:
    //  The FFT functions generate output in bit reversed order.
    //  Use this function to re-arrange them into order of increasing frequency.
    //
    // REMARKS:
    //  Exponential values and bits correspond, so the reversed upper index can be omitted depending on the number of exponents.
    //
    // PARAMETERS:
    //  pOutput     - [out] output buffer, receives samples in order of increasing frequency, cannot overlap pInput, must have at least (1<<uLog2Length)/4 elements
    //  pInput      - [in]  input buffer, samples in bit reversed order as generated by FFT functions, cannot overlap pOutput, must have at least (1<<uLog2Length)/4 elements
    //  uLog2Length - [in]  LOG (base 2) of FFT length in samples, must be >= 2
    //----------------------------------------------------------------------------------
    inline void FFTUnswizzle (
        _Out_writes_((1 << uLog2Length) / 4) XMVECTOR* __restrict pOutput,
        _In_reads_((1 << uLog2Length) / 4) const XMVECTOR* __restrict pInput,
        _In_ const size_t uLog2Length) noexcept
    {
        assert(pOutput);
        assert(pInput);
        assert(uLog2Length >= 2);
        _Analysis_assume_(uLog2Length >= 2);

        float* __restrict pfOutput = reinterpret_cast<float*>(pOutput);
        const size_t uLength = size_t(1) << (uLog2Length - 2);

        static const unsigned char cSwizzleTable[256] = {
            0x00, 0x40, 0x80, 0xC0, 0x10, 0x50, 0x90, 0xD0, 0x20, 0x60, 0xA0, 0xE0, 0x30, 0x70, 0xB0, 0xF0,
            0x04, 0x44, 0x84, 0xC4, 0x14, 0x54, 0x94, 0xD4, 0x24, 0x64, 0xA4, 0xE4, 0x34, 0x74, 0xB4, 0xF4,
            0x08, 0x48, 0x88, 0xC8, 0x18, 0x58, 0x98, 0xD8, 0x28, 0x68, 0xA8, 0xE8, 0x38, 0x78, 0xB8, 0xF8,
            0x0C, 0x4C, 0x8C, 0xCC, 0x1C, 0x5C, 0x9C, 0xDC, 0x2C, 0x6C, 0xAC, 0xEC, 0x3C, 0x7C, 0xBC, 0xFC,
            0x01, 0x41, 0x81, 0xC1, 0x11, 0x51, 0x91, 0xD1, 0x21, 0x61, 0xA1, 0xE1, 0x31, 0x71, 0xB1, 0xF1,
            0x05, 0x45, 0x85, 0xC5, 0x15, 0x55, 0x95, 0xD5, 0x25, 0x65, 0xA5, 0xE5, 0x35, 0x75, 0xB5, 0xF5,
            0x09, 0x49, 0x89, 0xC9, 0x19, 0x59, 0x99, 0xD9, 0x29, 0x69, 0xA9, 0xE9, 0x39, 0x79, 0xB9, 0xF9,
            0x0D, 0x4D, 0x8D, 0xCD, 0x1D, 0x5D, 0x9D, 0xDD, 0x2D, 0x6D, 0xAD, 0xED, 0x3D, 0x7D, 0xBD, 0xFD,
            0x02, 0x42, 0x82, 0xC2, 0x12, 0x52, 0x92, 0xD2, 0x22, 0x62, 0xA2, 0xE2, 0x32, 0x72, 0xB2, 0xF2,
            0x06, 0x46, 0x86, 0xC6, 0x16, 0x56, 0x96, 0xD6, 0x26, 0x66, 0xA6, 0xE6, 0x36, 0x76, 0xB6, 0xF6,
            0x0A, 0x4A, 0x8A, 0xCA, 0x1A, 0x5A, 0x9A, 0xDA, 0x2A, 0x6A, 0xAA, 0xEA, 0x3A, 0x7A, 0xBA, 0xFA,
            0x0E, 0x4E, 0x8E, 0xCE, 0x1E, 0x5E, 0x9E, 0xDE, 0x2E, 0x6E, 0xAE, 0xEE, 0x3E, 0x7E, 0xBE, 0xFE,
            0x03, 0x43, 0x83, 0xC3, 0x13, 0x53, 0x93, 0xD3, 0x23, 0x63, 0xA3, 0xE3, 0x33, 0x73, 0xB3, 0xF3,
            0x07, 0x47, 0x87, 0xC7, 0x17, 0x57, 0x97, 0xD7, 0x27, 0x67, 0xA7, 0xE7, 0x37, 0x77, 0xB7, 0xF7,
            0x0B, 0x4B, 0x8B, 0xCB, 0x1B, 0x5B, 0x9B, 0xDB, 0x2B, 0x6B, 0xAB, 0xEB, 0x3B, 0x7B, 0xBB, 0xFB,
            0x0F, 0x4F, 0x8F, 0xCF, 0x1F, 0x5F, 0x9F, 0xDF, 0x2F, 0x6F, 0xAF, 0xEF, 0x3F, 0x7F, 0xBF, 0xFF
        };
        if ((uLog2Length & 1) == 0)
        {
            // even powers of two
            const size_t uRev32 = 32 - uLog2Length;
            for (size_t uIndex = 0; uIndex < uLength; ++uIndex)
            {
                XMFLOAT4A f4a;
                XMStoreFloat4A(&f4a, pInput[uIndex]);
                const size_t n = uIndex * 4;
                const size_t uAddr = (static_cast<size_t>(cSwizzleTable[n & 0xff]) << 24) |
                    (static_cast<size_t>(cSwizzleTable[(n >> 8) & 0xff]) << 16) |
                    (static_cast<size_t>(cSwizzleTable[(n >> 16) & 0xff]) << 8) |
                    (static_cast<size_t>(cSwizzleTable[(n >> 24)]));
                pfOutput[uAddr >> uRev32] = f4a.x;
                pfOutput[(0x40000000 | uAddr) >> uRev32] = f4a.y;
                pfOutput[(0x80000000 | uAddr) >> uRev32] = f4a.z;
                pfOutput[(0xC0000000 | uAddr) >> uRev32] = f4a.w;
            }
        }
        else
        {
            // odd powers of two
            const size_t uRev7 = size_t(1) << (uLog2Length - 3);
            const size_t uRev32 = 32 - (uLog2Length - 3);
            for (size_t uIndex = 0; uIndex < uLength; ++uIndex)
            {
                XMFLOAT4A f4a;
                XMStoreFloat4A(&f4a, pInput[uIndex]);
                const size_t n = (uIndex >> 1);
                size_t uAddr = (((static_cast<size_t>(cSwizzleTable[n & 0xff]) << 24) |
                    (static_cast<size_t>(cSwizzleTable[(n >> 8) & 0xff]) << 16) |
                    (static_cast<size_t>(cSwizzleTable[(n >> 16) & 0xff]) << 8) |
                    (static_cast<size_t>(cSwizzleTable[(n >> 24)]))) >> uRev32) |
                    ((uIndex & 1) * uRev7 * 4);
                pfOutput[uAddr] = f4a.x;
                uAddr += uRev7;
                pfOutput[uAddr] = f4a.y;
                uAddr += uRev7;
                pfOutput[uAddr] = f4a.z;
                uAddr += uRev7;
                pfOutput[uAddr] = f4a.w;
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
#ifdef _MSC_VER
#pragma warning(suppress: 6101)
#endif
    inline void FFTPolar(
        _Out_writes_(uLength / 4) XMVECTOR* __restrict pOutput,
        _In_reads_(uLength / 4) const XMVECTOR* __restrict pInputReal,
        _In_reads_(uLength / 4) const XMVECTOR* __restrict pInputImaginary,
        _In_ const size_t uLength) noexcept
    {
        using namespace DirectX;

        assert(pOutput);
        assert(pInputReal);
        assert(pInputImaginary);
        assert(uLength >= 4);
        _Analysis_assume_(uLength >= 4);
        assert(ISPOWEROF2(uLength));

        const float flOneOverLength = 1.0f / float(uLength);

        // result = sqrtf((real/uLength)^2 + (imaginary/uLength)^2) * 2
        const XMVECTOR vOneOverLength = XMVectorReplicate(flOneOverLength);

        for (size_t uIndex = 0; uIndex < (uLength >> 2); ++uIndex)
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
    inline void Deinterleave (
        _Out_writes_((uChannelCount * uFrameCount) / 4) XMVECTOR* __restrict pOutput,
        _In_reads_((uChannelCount * uFrameCount) / 4) const XMVECTOR* __restrict pInput,
        _In_ const size_t uChannelCount,
        _In_ const size_t uFrameCount) noexcept
    {
        assert(pOutput);
        assert(pInput);
        assert(uChannelCount > 1);
        assert(uFrameCount > 0);

        float* __restrict pfOutput = reinterpret_cast<float* __restrict>(pOutput);
        const float* __restrict pfInput  = reinterpret_cast<const float* __restrict>(pInput);

        for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
        {
            for (size_t uFrame = 0; uFrame < uFrameCount; ++uFrame)
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
    inline void Interleave(
        _Out_writes_((uChannelCount * uFrameCount) / 4) XMVECTOR* __restrict pOutput,
        _In_reads_((uChannelCount * uFrameCount) / 4) const XMVECTOR* __restrict pInput,
        _In_ const size_t uChannelCount,
        _In_ const size_t uFrameCount) noexcept
    {
        assert(pOutput);
        assert(pInput);
        assert(uChannelCount > 1);
        assert(uFrameCount > 0);

        float* __restrict pfOutput = reinterpret_cast<float* __restrict>(pOutput);
        const float* __restrict pfInput  = reinterpret_cast<const float* __restrict>(pInput);

        for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
        {
            for (size_t uFrame = 0; uFrame < uFrameCount; ++uFrame)
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
    inline void FFTInterleaved(
        _Inout_updates_(((1 << uLog2Length) * uChannelCount) / 4) XMVECTOR* __restrict pReal,
        _Out_writes_(((1 << uLog2Length) * uChannelCount) / 4) XMVECTOR* __restrict pImaginary,
        _In_reads_(1 << uLog2Length) const XMVECTOR* __restrict pUnityTable,
        _In_ const size_t uChannelCount,
        _In_ const size_t uLog2Length) noexcept
    {
        assert(pReal);
        assert(pImaginary);
        assert(pUnityTable);
        assert(reinterpret_cast<uintptr_t>(pReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pImaginary) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pUnityTable) % 16 == 0);
        assert(uChannelCount > 0 && uChannelCount <= 6);
        assert(uLog2Length >= 2 && uLog2Length <= 9);

        XM_ALIGNED_DATA(16) XMVECTOR vRealTemp[768];
        XM_ALIGNED_DATA(16) XMVECTOR vImaginaryTemp[768];
        const size_t uLength = size_t(1) << uLog2Length;

        if (uChannelCount > 1)
        {
            Deinterleave(vRealTemp, pReal, uChannelCount, uLength);
        }
        else
        {
            memcpy_s(vRealTemp, sizeof(vRealTemp), pReal, (uLength >> 2) * sizeof(XMVECTOR));
        }

        memset(vImaginaryTemp, 0, (uChannelCount * (uLength >> 2)) * sizeof(XMVECTOR));

        if (uLength > 16)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)], pUnityTable, uLength);
            }
        }
        else if (uLength == 16)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT16(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)]);
            }
        }
        else if (uLength == 8)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT8(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)]);
            }
        }
        else if (uLength == 4)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT4(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)]);
            }
        }

        for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
        {
            FFTUnswizzle(&pReal[uChannel * (uLength >> 2)], &vRealTemp[uChannel * (uLength >> 2)], uLog2Length);
            FFTUnswizzle(&pImaginary[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)], uLog2Length);
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
    inline void IFFTDeinterleaved(
        _Inout_updates_(((1 << uLog2Length) * uChannelCount) / 4) XMVECTOR* __restrict pReal,
        _In_reads_(((1 << uLog2Length) * uChannelCount) / 4) const XMVECTOR* __restrict pImaginary,
        _In_reads_(1 << uLog2Length) const XMVECTOR* __restrict pUnityTable,
        _In_ const size_t uChannelCount,
        _In_ const size_t uLog2Length) noexcept
    {
        using namespace DirectX;

        assert(pReal);
        assert(pImaginary);
        assert(pUnityTable);
        assert(reinterpret_cast<uintptr_t>(pReal) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pImaginary) % 16 == 0);
        assert(reinterpret_cast<uintptr_t>(pUnityTable) % 16 == 0);
        assert(uChannelCount > 0 && uChannelCount <= 6);
        _Analysis_assume_(uChannelCount > 0 && uChannelCount <= 6);
        assert(uLog2Length >= 2 && uLog2Length <= 9);
        _Analysis_assume_(uLog2Length >= 2 && uLog2Length <= 9);

        XM_ALIGNED_DATA(16) XMVECTOR vRealTemp[768] = {};
        XM_ALIGNED_DATA(16) XMVECTOR vImaginaryTemp[768] = {};

        const size_t uLength = size_t(1) << uLog2Length;

        const XMVECTOR vRnp = XMVectorReplicate(1.0f / float(uLength));
        const XMVECTOR vRnm = XMVectorReplicate(-1.0f / float(uLength));
        for (size_t u = 0; u < uChannelCount * (uLength >> 2); u++)
        {
            vRealTemp[u] = XMVectorMultiply(pReal[u], vRnp);
            vImaginaryTemp[u] = XMVectorMultiply(pImaginary[u], vRnm);
        }

        if (uLength > 16)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)], pUnityTable, uLength);
            }
        }
        else if (uLength == 16)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT16(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)]);
            }
        }
        else if (uLength == 8)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT8(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)]);
            }
        }
        else if (uLength == 4)
        {
            for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
            {
                FFT4(&vRealTemp[uChannel * (uLength >> 2)], &vImaginaryTemp[uChannel * (uLength >> 2)]);
            }
        }

        for (size_t uChannel = 0; uChannel < uChannelCount; ++uChannel)
        {
            FFTUnswizzle(&vImaginaryTemp[uChannel * (uLength >> 2)], &vRealTemp[uChannel * (uLength >> 2)], uLog2Length);
        }

        if (uChannelCount > 1)
        {
            Interleave(pReal, vImaginaryTemp, uChannelCount, uLength);
        }
        else
        {
            memcpy_s(pReal, uLength * uChannelCount * sizeof(float), vImaginaryTemp, (uLength >> 2) * sizeof(XMVECTOR));
        }
    }

} // namespace XDSP

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif
