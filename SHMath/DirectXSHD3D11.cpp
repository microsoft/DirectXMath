//-------------------------------------------------------------------------------------
// DirectXSHD3D11.cpp -- C++ Spherical Harmonics Math Library
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//  
// Copyright (c) Microsoft Corporation. All rights reserved.
//
// http://go.microsoft.com/fwlink/p/?LinkId=262885
//-------------------------------------------------------------------------------------

#include "DirectXSH.h"

#include <d3d11.h>

#include <DirectXPackedVector.h>

#include <assert.h>
#include <memory>
#include <malloc.h>

namespace
{
struct aligned_deleter { void operator()(void* p) { _aligned_free(p); } };

typedef std::unique_ptr<DirectX::XMVECTOR, aligned_deleter> ScopedAlignedArrayXMVECTOR;

template<class T> class ScopedObject
{
public:
    explicit ScopedObject( T *p = 0 ) : _pointer(p) {}
    ~ScopedObject()
    {
        if ( _pointer )
        {
            _pointer->Release();
            _pointer = nullptr;
        }
    }

    bool IsNull() const { return (!_pointer); }

    T& operator*() { return *_pointer; }
    T* operator->() { return _pointer; }
    T** operator&() { return &_pointer; }

    void Reset(T *p = 0) { if ( _pointer ) { _pointer->Release(); } _pointer = p; }

    T* Get() const { return _pointer; }

private:
    ScopedObject(const ScopedObject&);
    ScopedObject& operator=(const ScopedObject&);
        
    T* _pointer;
};

//-------------------------------------------------------------------------------------
// This code is lifted from DirectXTex http://directxtex.codeplex.com/
// If you need additional DXGI format support, see DirectXTexConvert.cpp
//-------------------------------------------------------------------------------------
#define LOAD_SCANLINE( type, func )\
        if ( size >= sizeof(type) )\
        {\
            const type * __restrict sPtr = reinterpret_cast<const type*>(pSource);\
            for( size_t icount = 0; icount < ( size - sizeof(type) + 1 ); icount += sizeof(type) )\
            {\
                if ( dPtr >= ePtr ) break;\
                *(dPtr++) = func( sPtr++ );\
            }\
            return true;\
        }\
        return false;

#define LOAD_SCANLINE3( type, func, defvec )\
        if ( size >= sizeof(type) )\
        {\
            const type * __restrict sPtr = reinterpret_cast<const type*>(pSource);\
            for( size_t icount = 0; icount < ( size - sizeof(type) + 1 ); icount += sizeof(type) )\
            {\
                XMVECTOR v = func( sPtr++ );\
                if ( dPtr >= ePtr ) break;\
                *(dPtr++) = XMVectorSelect( defvec, v, g_XMSelect1110 );\
            }\
            return true;\
        }\
        return false;

#define LOAD_SCANLINE2( type, func, defvec )\
        if ( size >= sizeof(type) )\
        {\
            const type * __restrict sPtr = reinterpret_cast<const type*>(pSource);\
            for( size_t icount = 0; icount < ( size - sizeof(type) + 1 ); icount += sizeof(type) )\
            {\
                XMVECTOR v = func( sPtr++ );\
                if ( dPtr >= ePtr ) break;\
                *(dPtr++) = XMVectorSelect( defvec, v, g_XMSelect1100 );\
            }\
            return true;\
        }\
        return false;

#pragma warning(push)
#pragma warning(disable : 6101)
_Success_(return)
static bool _LoadScanline( _Out_writes_(count) DirectX::XMVECTOR* pDestination, _In_ size_t count,
                           _In_reads_bytes_(size) LPCVOID pSource, _In_ size_t size, _In_ DXGI_FORMAT format )
{
    assert( pDestination && count > 0 && (((uintptr_t)pDestination & 0xF) == 0) );
    assert( pSource && size > 0 );

    using namespace DirectX;
    using namespace DirectX::PackedVector;

    XMVECTOR* __restrict dPtr = pDestination;
    if ( !dPtr )
        return false;

    const XMVECTOR* ePtr = pDestination + count;

    switch( format )
    {
    case DXGI_FORMAT_R32G32B32A32_FLOAT:
        {
            size_t msize = (size > (sizeof(XMVECTOR)*count)) ? (sizeof(XMVECTOR)*count) : size;
            memcpy_s( dPtr, sizeof(XMVECTOR)*count, pSource, msize );
        }
        return true;

    case DXGI_FORMAT_R32G32B32_FLOAT:
        LOAD_SCANLINE3( XMFLOAT3, XMLoadFloat3, g_XMIdentityR3 )
            
    case DXGI_FORMAT_R16G16B16A16_FLOAT:
        LOAD_SCANLINE( XMHALF4, XMLoadHalf4 )

    case DXGI_FORMAT_R32G32_FLOAT:
        LOAD_SCANLINE2( XMFLOAT2, XMLoadFloat2, g_XMIdentityR3 )

    case DXGI_FORMAT_R11G11B10_FLOAT:
        LOAD_SCANLINE3( XMFLOAT3PK, XMLoadFloat3PK, g_XMIdentityR3 );

    case DXGI_FORMAT_R16G16_FLOAT:
        LOAD_SCANLINE2( XMHALF2, XMLoadHalf2, g_XMIdentityR3 )

    case DXGI_FORMAT_R32_FLOAT:
        if ( size >= sizeof(float) )
        {
            const float* __restrict sPtr = reinterpret_cast<const float*>(pSource);
            for( size_t icount = 0; icount < size; icount += sizeof(float) )
            {
                XMVECTOR v = XMLoadFloat( sPtr++ );
                if ( dPtr >= ePtr ) break;
                *(dPtr++) = XMVectorSelect( g_XMIdentityR3, v, g_XMSelect1000 );
            }
            return true;
        }
        return false;

    case DXGI_FORMAT_R16_FLOAT:
        if ( size >= sizeof(HALF) )
        {
            const HALF * __restrict sPtr = reinterpret_cast<const HALF*>(pSource);
            for( size_t icount = 0; icount < size; icount += sizeof(HALF) )
            {
                if ( dPtr >= ePtr ) break;
                *(dPtr++) = XMVectorSet( XMConvertHalfToFloat(*sPtr++), 0.f, 0.f, 1.f );
            }
            return true;
        }
        return false;

    default:
        return false;
    }
}
#pragma warning(pop)

}; // namespace anonymous

namespace DirectX
{

//-------------------------------------------------------------------------------------
// Projects a function represented in a cube map into spherical harmonics.
//
// http://msdn.microsoft.com/en-us/library/windows/desktop/ff476300.aspx
//-------------------------------------------------------------------------------------
HRESULT SHProjectCubeMap( _In_ ID3D11DeviceContext *context,
                          _In_ size_t order,
                          _In_ ID3D11Texture2D *cubeMap,
                          _Out_writes_opt_(order*order) float *resultR,
                          _Out_writes_opt_(order*order) float *resultG,
                          _Out_writes_opt_(order*order) float* resultB )
{
    if ( !context || !cubeMap )
        return E_INVALIDARG;

    if ( order < XM_SH_MINORDER || order > XM_SH_MAXORDER )
        return E_INVALIDARG;

    D3D11_TEXTURE2D_DESC desc;
    cubeMap->GetDesc( &desc );

    if ( (desc.ArraySize != 6)
         || (desc.Width != desc.Height)
         || (desc.SampleDesc.Count > 1) )
         return E_FAIL;

    switch( desc.Format )
    {
    case DXGI_FORMAT_R32G32B32A32_FLOAT:
    case DXGI_FORMAT_R32G32B32_FLOAT:
    case DXGI_FORMAT_R16G16B16A16_FLOAT:
    case DXGI_FORMAT_R32G32_FLOAT:
    case DXGI_FORMAT_R11G11B10_FLOAT:
    case DXGI_FORMAT_R16G16_FLOAT:
    case DXGI_FORMAT_R32_FLOAT:
    case DXGI_FORMAT_R16_FLOAT:
        // See _LoadScanline to support more pixel formats
        break;

    default:
        return E_FAIL;
    }

    //--- Create a staging resource copy (if needed) to be able to read data
    ID3D11Texture2D* texture = nullptr;

    ScopedObject<ID3D11Texture2D> staging;
    if ( !(desc.CPUAccessFlags & D3D11_CPU_ACCESS_READ) )
    {
        D3D11_TEXTURE2D_DESC sdesc = desc;
        sdesc.BindFlags = 0;
        sdesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        sdesc.Usage = D3D11_USAGE_STAGING;

        ScopedObject<ID3D11Device> device;
        context->GetDevice( &device );
        assert( !device.IsNull() );

        HRESULT hr = device->CreateTexture2D( &sdesc, nullptr, &staging );
        if ( FAILED(hr) )
            return hr;

        context->CopyResource( staging.Get(), cubeMap );
            
        texture = staging.Get();
    }
    else
        texture = cubeMap;

    assert( texture != 0 );

    //--- Setup for SH projection
    ScopedAlignedArrayXMVECTOR scanline( reinterpret_cast<XMVECTOR*>( _aligned_malloc( sizeof(XMVECTOR)*desc.Width, 16 ) ) );
    if ( !scanline )
        return E_OUTOFMEMORY;

    assert( desc.Width > 0 );
    float fSize = static_cast<float>( desc.Width );
    float fPicSize = 1.0f / fSize;

    // index from [0,W-1], f(0) maps to -1 + 1/W, f(W-1) maps to 1 - 1/w
    // linear function x*S +B, 1st constraint means B is (-1+1/W), plug into
    // second and solve for S: S = 2*(1-1/W)/(W-1). The old code that did 
    // this was incorrect - but only for computing the differential solid
    // angle, where the final value was 1.0 instead of 1-1/w...

    float fB = -1.0f + 1.0f/fSize;
    float fS = ( desc.Width > 1 ) ? (2.0f*(1.0f-1.0f/fSize)/(fSize-1.0f)) : 0.f;

    // clear out accumulation variables
    float fWt = 0.0f;

    if ( resultR )
        memset( resultR, 0, sizeof(float)*order*order );
    if ( resultG )
        memset( resultG, 0, sizeof(float)*order*order );
    if ( resultB )
        memset( resultB, 0, sizeof(float)*order*order );

    float shBuff[XM_SH_MAXORDER*XM_SH_MAXORDER];
    float shBuffB[XM_SH_MAXORDER*XM_SH_MAXORDER];

    //--- Process each face of the cubemap
    for (UINT face=0; face < 6; ++face )
    {
        UINT dindex = D3D11CalcSubresource( 0, face, desc.MipLevels );

        D3D11_MAPPED_SUBRESOURCE mapped;
        HRESULT hr = context->Map( texture, dindex, D3D11_MAP_READ, 0, &mapped );
        if ( FAILED(hr) )
            return hr;

        const uint8_t *pSrc = reinterpret_cast<const uint8_t*>(mapped.pData);
        for( UINT y=0; y < desc.Height; ++y )
        {
            XMVECTOR* ptr = scanline.get();
            if ( !_LoadScanline( ptr, desc.Width, pSrc, mapped.RowPitch, desc.Format ) )
            {
                context->Unmap( texture, dindex );
                return E_FAIL;
            }

            const float fV = y*fS + fB;

            XMVECTOR* pixel = ptr;
            for( UINT x=0; x < desc.Width; ++x, ++pixel )
            {
                const float fU = x*fS + fB;

                float ix, iy, iz;
                switch( face )
                {
                case 0: // Positive X
                    iz = 1.0f - (2.0f * (float)x + 1.0f) * fPicSize;
                    iy = 1.0f - (2.0f * (float)y + 1.0f) * fPicSize;
                    ix = 1.0f;
                    break;

                case 1: // Negative X
                    iz = -1.0f + (2.0f * (float)x + 1.0f) * fPicSize;
                    iy =  1.0f - (2.0f * (float)y + 1.0f) * fPicSize;
                    ix = -1;
                    break;

                case 2: // Positive Y
                    iz = -1.0f + (2.0f * (float)y + 1.0f) * fPicSize;
                    iy =  1.0f;
                    ix = -1.0f + (2.0f * (float)x + 1.0f) * fPicSize;
                    break;

                case 3: // Negative Y
                    iz =  1.0f - (2.0f * (float)y + 1.0f) * fPicSize;
                    iy = -1.0f;
                    ix = -1.0f + (2.0f * (float)x + 1.0f) * fPicSize;  
                    break;

                case 4: // Positive Z
                    iz =  1.0f;
                    iy =  1.0f - (2.0f * (float)y + 1.0f) * fPicSize;
                    ix = -1.0f + (2.0f * (float)x + 1.0f) * fPicSize;  
                    break;

                case 5: // Negative Z
                    iz = -1.0f;
                    iy =  1.0f - (2.0f * (float)y + 1.0f) * fPicSize;
                    ix =  1.0f - (2.0f * (float)x + 1.0f) * fPicSize;
                    break;

                default:
                    ix = iy = iz = 0.f;
                    assert(false);
                    break;
                }

                XMVECTOR dir = XMVectorSet( ix, iy, iz, 0 );
                dir = XMVector3Normalize( dir );

                const float fDiffSolid = 4.0f/((1.0f + fU*fU + fV*fV)*sqrtf(1.0f + fU*fU+fV*fV));
                fWt += fDiffSolid;

                XMSHEvalDirection(shBuff,order,dir);

                XMFLOAT3A clr;
                XMStoreFloat3A( &clr, *pixel );

                if ( resultR ) XMSHAdd(resultR,order,resultR, XMSHScale(shBuffB,order,shBuff,clr.x*fDiffSolid) );
                if ( resultG ) XMSHAdd(resultG,order,resultG, XMSHScale(shBuffB,order,shBuff,clr.y*fDiffSolid) );
                if ( resultB ) XMSHAdd(resultB,order,resultB, XMSHScale(shBuffB,order,shBuff,clr.z*fDiffSolid) );
           }

            pSrc += mapped.RowPitch;
        }

        context->Unmap( texture, dindex );
    }

    const float fNormProj = (4.0f*XM_PI)/fWt;

    if ( resultR ) XMSHScale(resultR,order,resultR,fNormProj);
    if ( resultG ) XMSHScale(resultG,order,resultG,fNormProj);
    if ( resultB ) XMSHScale(resultB,order,resultB,fNormProj);

    return S_OK;
}

}; // namespace DirectX
