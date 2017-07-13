//-------------------------------------------------------------------------------------
// Stereo3DMatrixHelper.h -- SIMD C++ Math helper for Stereo 3D matricies
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//  
// Copyright (c) Microsoft Corporation. All rights reserved.
//-------------------------------------------------------------------------------------

#include "Stereo3DMatrixHelper.h"

using namespace DirectX;

//------------------------------------------------------------------------------

void StereoCreateDefaultParameters
(
    STEREO_PARAMETERS* pStereoParameters
)
{
    assert( pStereoParameters != nullptr );

    // Default assumption is 1920x1200 resolution, a 22" LCD monitor, and a 2' viewing distance
    pStereoParameters->fViewerDistanceInches = 24.0f;
    pStereoParameters->fPixelResolutionWidth = 1920.0f;
    pStereoParameters->fPixelResolutionHeight = 1200.0f;
    pStereoParameters->fDisplaySizeInches = 22.0f;

    pStereoParameters->fStereoSeparationFactor = 1.0f;
    pStereoParameters->fStereoExaggerationFactor = 1.0f;
}

//------------------------------------------------------------------------------

static inline bool StereoProjectionHelper
(
    const STEREO_PARAMETERS* pStereoParameters,
    _Out_ float* fVirtualProjection,
    _Out_ float* zNearWidth,
    _Out_ float* zNearHeight,
    float FovAngleY,
    float AspectHByW,
    float NearZ
)
{
// note that most people have difficulty fusing images into 3D
// if the separation equals even just the human average. by 
// reducing the separation (interocular distance) by 1/2, we
// guarantee a larger subset of people will see full 3D

// the conservative setting should always be used. the only problem
// with the conservative setting is that the 3D effect will be less 
// impressive on smaller screens (which makes sense, since your eye
// cannot be tricked as easily based on the smaller fov). to simulate
// the effect of a larger screen, use the liberal settings (debug only)

// Conservative Settings: * max acuity angle: 0.8f degrees * interoc distance: 1.25 inches

// Liberal Settings: * max acuity angle: 1.6f degrees * interoc distance: 2.5f inches

// maximum visual accuity angle allowed is 3.2 degrees for 
// a physical scene, and 1.6 degrees for a virtual one. 
// thus we cannot allow an object to appear any closer to
// the viewer than 1.6 degrees (divided by two for most 
// half-angle calculations)

    static const float fMaxStereoDistance = 780; // inches (should be between 10 and 20m)
    static const float fMaxVisualAcuityAngle = 1.6f * ( XM_PI / 180.0f );  // radians
    static const float fInterocularDistance = 1.25f; // inches
    
    bool ComfortableResult = true;
    float fDisplayHeight, fDisplayWidth, fHalfInterocular, fHalfPixelWidth, fHalfMaximumAcuityAngle, fHalfWidth;
    float fMaxSeparationAcuityAngle, fMaxSeparationDistance, fRefinedMaxStereoDistance, fFovHalfAngle;
    float fRefinedMaxSeparationAcuityAngle, fPhysicalZNearDistance, fScalingFactor, fNearZSeparation, fNearZSeparation2;

    fDisplayHeight = pStereoParameters->fDisplaySizeInches / sqrtf( AspectHByW * AspectHByW + 1.0f );
    fDisplayWidth = fDisplayHeight * AspectHByW;
    fHalfInterocular = 0.5f * fInterocularDistance * pStereoParameters->fStereoExaggerationFactor;
    fHalfPixelWidth  = fDisplayWidth / pStereoParameters->fPixelResolutionWidth * 0.5f;
    fHalfMaximumAcuityAngle = fMaxVisualAcuityAngle * 0.5f * pStereoParameters->fStereoExaggerationFactor;
    fHalfWidth = fDisplayWidth * 0.5f;

    fMaxSeparationAcuityAngle = atanf( fHalfInterocular / fMaxStereoDistance );
    fMaxSeparationDistance = fHalfPixelWidth / tanf ( fMaxSeparationAcuityAngle );
    fRefinedMaxStereoDistance = fMaxStereoDistance - fMaxSeparationDistance;
    fFovHalfAngle = FovAngleY / 2.0f;

    if ( fRefinedMaxStereoDistance < 0.0f || fMaxSeparationDistance > 0.1f * fMaxStereoDistance )
    {
        // Pixel resolution is too low to offer a comfortable stereo experience
        ComfortableResult = false;
    }

    fRefinedMaxSeparationAcuityAngle = atanf( fHalfInterocular / ( fRefinedMaxStereoDistance ) );
    fPhysicalZNearDistance = fHalfInterocular / tanf( fHalfMaximumAcuityAngle );
    fScalingFactor = fHalfMaximumAcuityAngle / atanf( fHalfInterocular / pStereoParameters->fViewerDistanceInches );

    fNearZSeparation = tanf( fRefinedMaxSeparationAcuityAngle ) * ( fRefinedMaxStereoDistance - fPhysicalZNearDistance );
    fNearZSeparation2 = fHalfInterocular * ( fRefinedMaxStereoDistance - fPhysicalZNearDistance ) / fRefinedMaxStereoDistance;

    (*zNearHeight) = cosf( fFovHalfAngle ) / sinf( fFovHalfAngle );
    (*zNearWidth)  = (*zNearHeight) / AspectHByW;
    (*fVirtualProjection) = ( fNearZSeparation * NearZ * (*zNearWidth * 4.0f) ) / ( 2.0f * NearZ );

    return ComfortableResult;
}

//------------------------------------------------------------------------------

XMMATRIX StereoProjectionFovLH
(
    const STEREO_PARAMETERS* pStereoParameters,
    STEREO_CHANNEL Channel,
    float FovAngleY,
    float AspectHByW,
    float NearZ,
    float FarZ,
    STEREO_MODE StereoMode
)
{
    float fVirtualProjection = 0.0f;
    float zNearWidth = 0.0f;
    float zNearHeight = 0.0f;
    float fInvertedAngle;
    XMMATRIX patchedProjection, proj;
    STEREO_PARAMETERS DefaultParameters;

    assert( Channel == STEREO_CHANNEL_LEFT || Channel == STEREO_CHANNEL_RIGHT );
    assert( StereoMode == STEREO_MODE_NORMAL || StereoMode == STEREO_MODE_INVERTED );
    assert(!XMScalarNearEqual(FovAngleY, 0.0f, 0.00001f * 2.0f));
    assert(!XMScalarNearEqual(AspectHByW, 0.0f, 0.00001f));
    assert(!XMScalarNearEqual(FarZ, NearZ, 0.00001f));

    proj = XMMatrixIdentity();

    if( pStereoParameters == nullptr )
    {
        StereoCreateDefaultParameters( &DefaultParameters );
        pStereoParameters = &DefaultParameters;
    }

    assert( pStereoParameters->fStereoSeparationFactor >= 0.0f && pStereoParameters->fStereoSeparationFactor <= 1.0f );
    assert( pStereoParameters->fStereoExaggerationFactor >= 1.0f && pStereoParameters->fStereoExaggerationFactor <= 2.0f );

    StereoProjectionHelper( pStereoParameters, &fVirtualProjection, &zNearWidth, &zNearHeight, FovAngleY, AspectHByW, NearZ );

    fVirtualProjection *= pStereoParameters->fStereoSeparationFactor; // incorporate developer defined bias

    //
    // By applying a translation, we are forcing our cameras to be parallel 
    //

    fInvertedAngle = atanf( fVirtualProjection / ( 2.0f * NearZ ) );

    proj = XMMatrixPerspectiveFovLH( FovAngleY, AspectHByW, NearZ, FarZ );

    if ( Channel == STEREO_CHANNEL_LEFT )
    {
        if ( StereoMode > STEREO_MODE_NORMAL )
        {
            XMMATRIX rots, trans;
            rots = XMMatrixRotationY( fInvertedAngle );
            trans = XMMatrixTranslation( -fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( XMMatrixMultiply( rots, trans ), proj );
        }
        else
        {
            XMMATRIX trans;
            trans = XMMatrixTranslation( -fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( trans, proj );
        }
    }
    else
    {
        if ( StereoMode > STEREO_MODE_NORMAL )
        {            
            XMMATRIX rots, trans;
            rots = XMMatrixRotationY( -fInvertedAngle );
            trans = XMMatrixTranslation( fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( XMMatrixMultiply( rots, trans), proj );
        }
        else
        {
            XMMATRIX trans;
            trans = XMMatrixTranslation( fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( trans, proj );
        }
    }

    return patchedProjection;
}

//------------------------------------------------------------------------------

XMMATRIX StereoProjectionFovRH
(
    const STEREO_PARAMETERS* pStereoParameters,
    STEREO_CHANNEL Channel,
    float FovAngleY,
    float AspectHByW,
    float NearZ,
    float FarZ,
    STEREO_MODE StereoMode
)
{
    float fVirtualProjection = 0.0f;
    float zNearWidth = 0.0f;
    float zNearHeight = 0.0f;
    float fInvertedAngle;
    XMMATRIX patchedProjection, proj;
    STEREO_PARAMETERS DefaultParameters;

    assert( Channel == STEREO_CHANNEL_LEFT || Channel == STEREO_CHANNEL_RIGHT );
    assert( StereoMode == STEREO_MODE_NORMAL || StereoMode == STEREO_MODE_INVERTED );
    assert(!XMScalarNearEqual(FovAngleY, 0.0f, 0.00001f * 2.0f));
    assert(!XMScalarNearEqual(AspectHByW, 0.0f, 0.00001f));
    assert(!XMScalarNearEqual(FarZ, NearZ, 0.00001f));

    proj = XMMatrixIdentity();

    if( pStereoParameters == nullptr )
    {
        StereoCreateDefaultParameters( &DefaultParameters );
        pStereoParameters = &DefaultParameters;
    }

    assert( pStereoParameters->fStereoSeparationFactor >= 0.0f && pStereoParameters->fStereoSeparationFactor <= 1.0f );
    assert( pStereoParameters->fStereoExaggerationFactor >= 1.0f && pStereoParameters->fStereoExaggerationFactor <= 2.0f );

    StereoProjectionHelper( pStereoParameters, &fVirtualProjection, &zNearWidth, &zNearHeight, FovAngleY, AspectHByW, NearZ );

    fVirtualProjection *= pStereoParameters->fStereoSeparationFactor; // incorporate developer defined bias

    //
    // By applying a translation, we are forcing our cameras to be parallel 
    //

    fInvertedAngle = atanf( fVirtualProjection / ( 2.0f * NearZ ) );

    proj = XMMatrixPerspectiveFovRH( FovAngleY, AspectHByW, NearZ, FarZ );

    //
    // By applying a translation, we are forcing our cameras to be parallel 
    //

    if ( Channel == STEREO_CHANNEL_LEFT )
    {
        if ( StereoMode > STEREO_MODE_NORMAL )
        {
            XMMATRIX rots, trans;
            rots = XMMatrixRotationY( fInvertedAngle );
            trans = XMMatrixTranslation( -fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( XMMatrixMultiply( rots, trans ), proj );
        }
        else
        {
            XMMATRIX trans;
            trans = XMMatrixTranslation( -fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( trans, proj );
        }
    }
    else
    {
        if ( StereoMode > STEREO_MODE_NORMAL )
        {   
            XMMATRIX rots, trans;
            rots = XMMatrixRotationY( -fInvertedAngle );
            trans = XMMatrixTranslation( fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( XMMatrixMultiply( rots, trans), proj );
        }
        else
        {
            XMMATRIX trans;
            trans = XMMatrixTranslation( fVirtualProjection, 0, 0);
            patchedProjection = XMMatrixMultiply( trans, proj );
        }
    }

    return patchedProjection;
}
