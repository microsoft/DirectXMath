//-------------------------------------------------------------------------------------
// DirectXSH.h -- C++ Spherical Harmonics Math Library
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

#ifdef _MSC_VER
#pragma once
#endif

#define DIRECTX_SHMATH_VERSION 102

#include <DirectXMath.h>

#include <winerror.h>

struct ID3D11DeviceContext;
struct ID3D11Texture2D;

namespace DirectX
{
#if (DIRECTXMATH_VERSION < 305) && !defined(XM_CALLCONV)
#define XM_CALLCONV __fastcall
typedef const DirectX::XMVECTOR& HXMVECTOR;
typedef const DirectX::XMMATRIX& FXMMATRIX;
#endif

const size_t XM_SH_MINORDER = 2;
const size_t XM_SH_MAXORDER = 6;

float* XM_CALLCONV XMSHEvalDirection( _Out_writes_(order*order) float *result, _In_ size_t order, _In_ FXMVECTOR dir );

float* XM_CALLCONV XMSHRotate( _Out_writes_(order*order) float *result, _In_ size_t order, _In_ FXMMATRIX rotMatrix, _In_reads_(order*order) const float *input );

float* XMSHRotateZ( _Out_writes_(order*order) float *result, _In_ size_t order, _In_ float angle, _In_reads_(order*order) const float *input );

float* XMSHAdd( _Out_writes_(order*order) float *result, _In_ size_t order, _In_reads_(order*order) const float *inputA, _In_reads_(order*order) const float *inputB );

float* XMSHScale( _Out_writes_(order*order) float *result, _In_ size_t order, _In_reads_(order*order) const float *input, _In_ float scale );

float XMSHDot( _In_ size_t order, _In_reads_(order*order) const float *inputA, _In_reads_(order*order) const float *inputB );

float* XMSHMultiply( _Out_writes_(order*order) float *result, _In_ size_t order, _In_reads_(order*order) const float *inputF, _In_reads_(order*order) const float *inputG );

float* XMSHMultiply2( _Out_writes_(4) float *result, _In_reads_(4) const float *inputF, _In_reads_(4) const float *inputG );

float* XMSHMultiply3( _Out_writes_(9) float *result, _In_reads_(9) const float *inputF, _In_reads_(9) const float *inputG );

float* XMSHMultiply4( _Out_writes_(16) float *result, _In_reads_(16) const float *inputF, _In_reads_(16) const float *inputG );

float* XMSHMultiply5( _Out_writes_(25) float *result, _In_reads_(25) const float *inputF, _In_reads_(25) const float *inputG );

float* XMSHMultiply6( _Out_writes_(36) float *result, _In_reads_(36) const float *inputF, _In_reads_(36) const float *inputG );

bool XM_CALLCONV XMSHEvalDirectionalLight( _In_ size_t order, _In_ FXMVECTOR dir, _In_ FXMVECTOR color,
                                           _Out_writes_(order*order) float *resultR, _Out_writes_opt_(order*order) float *resultG, _Out_writes_opt_(order*order) float *resultB );

bool XM_CALLCONV XMSHEvalSphericalLight( _In_ size_t order, _In_ FXMVECTOR pos, _In_ float radius, _In_ FXMVECTOR color,
                                         _Out_writes_(order*order) float *resultR, _Out_writes_opt_(order*order) float *resultG, _Out_writes_opt_(order*order) float *resultB );

bool XM_CALLCONV XMSHEvalConeLight( _In_ size_t order, _In_ FXMVECTOR dir, _In_ float radius, _In_ FXMVECTOR color,
                                    _Out_writes_(order*order) float *resultR, _Out_writes_opt_(order*order) float *resultG, _Out_writes_opt_(order*order) float *resultB );

bool XM_CALLCONV XMSHEvalHemisphereLight( _In_ size_t order, _In_ FXMVECTOR dir, _In_ FXMVECTOR topColor, _In_ FXMVECTOR bottomColor,
                                          _Out_writes_(order*order) float *resultR, _Out_writes_opt_(order*order) float *resultG, _Out_writes_opt_(order*order) float *resultB );

HRESULT SHProjectCubeMap( _In_ ID3D11DeviceContext *context, _In_ size_t order, _In_ ID3D11Texture2D *cubeMap,
                          _Out_writes_opt_(order*order) float *resultR, _Out_writes_opt_(order*order) float *resultG, _Out_writes_opt_(order*order) float *resultB );

}; // namespace DirectX
