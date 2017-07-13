//// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
//// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
//// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
//// PARTICULAR PURPOSE.
////
//// Copyright (c) Microsoft Corporation. All rights reserved

#pragma once

// Stereo parameters are in the same units as the world.
struct StereoParameters
{
    float viewportWidth;        // viewport width
    float viewportHeight;       // viewport height
    float viewerDistance;       // distance from viewer
    float interocularDistance;  // interocular distance
};

StereoParameters CreateDefaultStereoParameters(
    float viewportWidthInches,
    float viewportHeightInches,
    float worldScaleInInches,
    float stereoExaggeration
    );

DirectX::XMMATRIX StereoProjectionFieldOfViewRightHand(
    const StereoParameters& parameters,
    float nearZ,
    float farZ,
    bool rightChannel
    );
