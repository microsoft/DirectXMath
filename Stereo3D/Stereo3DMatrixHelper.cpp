//// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
//// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
//// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
//// PARTICULAR PURPOSE.
////
//// Copyright (c) Microsoft Corporation. All rights reserved

#include "pch.h"
#include "Stereo3DMatrixHelper.h"

using namespace DirectX;

StereoParameters CreateDefaultStereoParameters(
    float viewportWidthInches,
    float viewportHeightInches,
    float worldScaleInInches,
    float stereoExaggeration
    )
{
    // The default stereo parameters produced by this method are based on two assumptions:
    // 1. The viewer's eyes are 24 inches from the display, and
    // 2. The viewer's eyes are separated by 1.25 inches (interocular distance.)
    const float DEFAULT_VIEWER_DISTANCE_IN_INCHES = 24.0f;
    const float DEFAULT_INTEROCULAR_DISTANCE_IN_INCHES = 1.25f;

    StereoParameters parameters;
    parameters.viewportWidth = viewportWidthInches / worldScaleInInches;
    parameters.viewportHeight = viewportHeightInches / worldScaleInInches;
    parameters.viewerDistance = DEFAULT_VIEWER_DISTANCE_IN_INCHES / worldScaleInInches;
    parameters.interocularDistance = DEFAULT_INTEROCULAR_DISTANCE_IN_INCHES / worldScaleInInches * stereoExaggeration;

    return parameters;
}

DirectX::XMMATRIX StereoProjectionFieldOfViewRightHand(
    const StereoParameters& parameters,
    float nearZ,
    float farZ,
    bool rightChannel
    )
{
    float yScale = 2.f * parameters.viewerDistance / parameters.viewportHeight;
    float xScale = 2.f * parameters.viewerDistance / parameters.viewportWidth;

    float mFactor = - parameters.interocularDistance / parameters.viewportWidth;

    if (!rightChannel)
    {
        mFactor = -mFactor;
    }

    float m22 = farZ / (nearZ - farZ);

    // Construct a stereo perspective projection matrix based on assumptions
    // about the viewer and specified stereo parameters. Note that compared
    // to a mono perspective projection matrix, there are two differences:
    //  - a non-zero x:z component (m20)
    //  - a non-zero x:w component (m30)
    // The values of these two factors affect both the x-offset between the
    // left and right eyes, as well as the depth at which they converge. The
    // math used to arrive at these values will often need to change depending
    // on the content being presented in order to ensure a comfortable viewing
    // experience. For example, the factors for rendering massive exterior
    // landscapes will be different than those used for rendering building
    // interiors. Because of this, developers are encouraged to experiment
    // with different techniques for generating these values.
    return XMMATRIX(
        xScale, 0, 0, 0,
        0, yScale, 0, 0,
        mFactor, 0, m22, -1,
        parameters.viewerDistance * mFactor, 0, nearZ * m22, 0
        );
}
