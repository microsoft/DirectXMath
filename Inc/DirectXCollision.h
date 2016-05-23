//-------------------------------------------------------------------------------------
// DirectXCollision.h -- C++ Collision Math library
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

#include "DirectXMath.h"

namespace DirectX
{

enum ContainmentType
{
    DISJOINT = 0,
    INTERSECTS = 1,
    CONTAINS = 2,
};

enum PlaneIntersectionType
{
    FRONT = 0,
    INTERSECTING = 1,
    BACK = 2,
};

struct BoundingBox;

#pragma warning(push)
#pragma warning(disable:4324)

//-------------------------------------------------------------------------------------
// Bounding sphere
//-------------------------------------------------------------------------------------
__declspec(align(16)) struct BoundingSphere
{
    XMFLOAT3 Center;            // Center of the sphere.
    float Radius;               // Radius of the sphere.

    // Creators
    BoundingSphere() : Center(0,0,0), Radius( 1.f ) {}
    BoundingSphere( _In_ const XMFLOAT3& center, _In_ float radius )
        : Center(center), Radius(radius) { assert( radius >= 0.f ); };
    BoundingSphere( _In_ const BoundingSphere& sp )
        : Center(sp.Center), Radius(sp.Radius) {}

    // Methods
    BoundingSphere& operator=( _In_ const BoundingSphere& sp ) { Center = sp.Center; Radius = sp.Radius; return *this; }

    void Transform( _Out_ BoundingSphere& Out, _In_ CXMMATRIX M ) const;
    void Transform( _Out_ BoundingSphere& Out, _In_ float Scale, _In_ FXMVECTOR Rotation, _In_ FXMVECTOR Translation ) const;
        // Transform the sphere

    ContainmentType Contains( _In_ FXMVECTOR Point ) const;
    ContainmentType Contains( _In_ FXMVECTOR V0, _In_ FXMVECTOR V1, _In_ FXMVECTOR V2 ) const;
    ContainmentType Contains( _In_ const BoundingSphere& sh ) const;
    ContainmentType Contains( _In_ const BoundingBox& box ) const;

    bool Intersects( _In_ const BoundingSphere& sh ) const;
    bool Intersects( _In_ const BoundingBox& box ) const;

    bool Intersects( _In_ FXMVECTOR V0, _In_ FXMVECTOR V1, _In_ FXMVECTOR V2 ) const;
        // Triangle-sphere test

    PlaneIntersectionType Intersects( _In_ FXMVECTOR Plane ) const;
        // Plane-sphere test
    
    bool Intersects( _In_ FXMVECTOR Origin, _In_ FXMVECTOR Direction, _Out_ float& Dist ) const;
        // Ray-sphere test

    // Static methods
    static void CreateMerged( _Out_ BoundingSphere& Out, _In_ const BoundingSphere& S1, _In_ const BoundingSphere& S2 );

    static void CreateFromBoundingBox( _Out_ BoundingSphere& Out, _In_ const BoundingBox& box );

    static void CreateFromPoints( _Out_ BoundingSphere& Out, _In_ size_t Count,
                                  _In_reads_bytes_(sizeof(XMFLOAT3)+Stride*(Count-1)) const XMFLOAT3* pPoints, _In_ size_t Stride );
};

//-------------------------------------------------------------------------------------
// Axis-aligned bounding box
//-------------------------------------------------------------------------------------
__declspec(align(16)) struct BoundingBox
{
    static const size_t CORNER_COUNT = 8;

    XMFLOAT3 Center;            // Center of the box.
    XMFLOAT3 Extents;           // Distance from the center to each side.

    // Creators
    BoundingBox() : Center(0,0,0), Extents( 1.f, 1.f, 1.f ) {}
    BoundingBox( _In_ const XMFLOAT3& center, _In_ const XMFLOAT3& extents )
        : Center(center), Extents(extents) { assert(extents.x >= 0 && extents.y >= 0 && extents.z >= 0); }
    BoundingBox( _In_ const BoundingBox& box ) : Center(box.Center), Extents(box.Extents) {}
    
    // Methods
    BoundingBox& operator=( _In_ const BoundingBox& box) { Center = box.Center; Extents = box.Extents; return *this; }

    void Transform( _Out_ BoundingBox& Out, _In_ CXMMATRIX M ) const;
    void Transform( _Out_ BoundingBox& Out, _In_ float Scale, _In_ FXMVECTOR Rotation, _In_ FXMVECTOR Translation ) const;

    void GetCorners( _Out_writes_(8) XMFLOAT3* Corners ) const;
        // Gets the 8 corners of the box

    ContainmentType Contains( _In_ FXMVECTOR Point ) const;
    ContainmentType Contains( _In_ FXMVECTOR V0, _In_ FXMVECTOR V1, _In_ FXMVECTOR V2 ) const;
    ContainmentType Contains( _In_ const BoundingSphere& sh ) const;
    ContainmentType Contains( _In_ const BoundingBox& box ) const;
    
    bool Intersects( _In_ const BoundingSphere& sh ) const;
    bool Intersects( _In_ const BoundingBox& box ) const;

    bool Intersects( _In_ FXMVECTOR V0, _In_ FXMVECTOR V1, _In_ FXMVECTOR V2 ) const;
        // Triangle-Box test

    PlaneIntersectionType Intersects( _In_ FXMVECTOR Plane ) const;
        // Plane-box test

    bool Intersects( _In_ FXMVECTOR Origin, _In_ FXMVECTOR Direction, _Out_ float& Dist ) const;
        // Ray-Box test

    // Static methods
    static void CreateMerged( _Out_ BoundingBox& Out, _In_ const BoundingBox& b1, _In_ const BoundingBox& b2 );

    static void CreateFromSphere( _Out_ BoundingBox& Out, _In_ const BoundingSphere& sh );

    static void CreateFromPoints( _Out_ BoundingBox& Out, _In_ FXMVECTOR pt1, _In_ FXMVECTOR pt2 );
    static void CreateFromPoints( _Out_ BoundingBox& Out, _In_ size_t Count,
                                  _In_reads_bytes_(sizeof(XMFLOAT3)+Stride*(Count-1)) const XMFLOAT3* pPoints, _In_ size_t Stride );
};

#pragma warning(pop)

/****************************************************************************
 *
 * Implementation
 *
 ****************************************************************************/

#pragma warning(push)
#pragma warning(disable:4068)

#pragma prefast(push)
#pragma prefast(disable : 25000, "FXMVECTOR is 16 bytes")

#include "DirectXCollision.inl"

#pragma prefast(pop)
#pragma warning(pop)

}; // namespace DirectX

