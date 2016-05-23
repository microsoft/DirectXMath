//-------------------------------------------------------------------------------------
// XNACollision.h
//  
// An opimtized collision library based on XNAMath
//  
// Microsoft XNA Developer Connection
// Copyright (c) Microsoft Corporation. All rights reserved.
//-------------------------------------------------------------------------------------

#pragma once

#ifndef _XNA_COLLISION_H_
#define _XNA_COLLISION_H_

#include <xnamath.h>

namespace XNA
{

//-----------------------------------------------------------------------------
// Bounding volumes structures.
//
// The bounding volume structures are setup for near minimum size because there
// are likely to be many of them, and memory bandwidth and space will be at a
// premium relative to CPU cycles on Xbox 360.
//-----------------------------------------------------------------------------

#pragma warning(push)
#pragma warning(disable: 4324)

_DECLSPEC_ALIGN_16_ struct Sphere
{
    XMFLOAT3 Center;            // Center of the sphere.
    FLOAT Radius;               // Radius of the sphere.
};

_DECLSPEC_ALIGN_16_ struct AxisAlignedBox
{
    XMFLOAT3 Center;            // Center of the box.
    XMFLOAT3 Extents;           // Distance from the center to each side.
};

_DECLSPEC_ALIGN_16_ struct OrientedBox
{
    XMFLOAT3 Center;            // Center of the box.
    XMFLOAT3 Extents;           // Distance from the center to each side.
    XMFLOAT4 Orientation;       // Unit quaternion representing rotation (box -> world).
};

_DECLSPEC_ALIGN_16_ struct Frustum
{
    XMFLOAT3 Origin;            // Origin of the frustum (and projection).
    XMFLOAT4 Orientation;       // Unit quaternion representing rotation.

    FLOAT RightSlope;           // Positive X slope (X/Z).
    FLOAT LeftSlope;            // Negative X slope.
    FLOAT TopSlope;             // Positive Y slope (Y/Z).
    FLOAT BottomSlope;          // Negative Y slope.
    FLOAT Near, Far;            // Z of the near plane and far plane.
};

#pragma warning(pop)

//-----------------------------------------------------------------------------
// Bounding volume construction.
//-----------------------------------------------------------------------------
VOID ComputeBoundingSphereFromPoints( Sphere* pOut, UINT Count, const XMFLOAT3* pPoints, UINT Stride );
VOID ComputeBoundingAxisAlignedBoxFromPoints( AxisAlignedBox* pOut, UINT Count, const XMFLOAT3* pPoints, UINT Stride );
VOID ComputeBoundingOrientedBoxFromPoints( OrientedBox* pOut, UINT Count, const XMFLOAT3* pPoints, UINT Stride );
VOID ComputeFrustumFromProjection( Frustum* pOut, XMMATRIX* pProjection );
VOID ComputePlanesFromFrustum( const Frustum* pVolume, XMVECTOR* pPlane0, XMVECTOR* pPlane1, XMVECTOR* pPlane2,
                               XMVECTOR* pPlane3, XMVECTOR* pPlane4, XMVECTOR* pPlane5 );



//-----------------------------------------------------------------------------
// Bounding volume transforms.
//-----------------------------------------------------------------------------
VOID TransformSphere( Sphere* pOut, const Sphere* pIn, FLOAT Scale, FXMVECTOR Rotation, FXMVECTOR Translation );
VOID TransformAxisAlignedBox( AxisAlignedBox* pOut, const AxisAlignedBox* pIn, FLOAT Scale, FXMVECTOR Rotation,
                              FXMVECTOR Translation );
VOID TransformOrientedBox( OrientedBox* pOut, const OrientedBox* pIn, FLOAT Scale, FXMVECTOR Rotation,
                           FXMVECTOR Translation );
VOID TransformFrustum( Frustum* pOut, const Frustum* pIn, FLOAT Scale, FXMVECTOR Rotation, FXMVECTOR Translation );



//-----------------------------------------------------------------------------
// Intersection testing routines.
//-----------------------------------------------------------------------------
BOOL IntersectPointSphere( FXMVECTOR Point, const Sphere* pVolume );
BOOL IntersectPointAxisAlignedBox( FXMVECTOR Point, const AxisAlignedBox* pVolume );
BOOL IntersectPointOrientedBox( FXMVECTOR Point, const OrientedBox* pVolume );
BOOL IntersectPointFrustum( FXMVECTOR Point, const Frustum* pVolume );
BOOL IntersectRayTriangle( FXMVECTOR Origin, FXMVECTOR Direction, FXMVECTOR V0, CXMVECTOR V1, CXMVECTOR V2,
                           FLOAT* pDist );
BOOL IntersectRaySphere( FXMVECTOR Origin, FXMVECTOR Direction, const Sphere* pVolume, FLOAT* pDist );
BOOL IntersectRayAxisAlignedBox( FXMVECTOR Origin, FXMVECTOR Direction, const AxisAlignedBox* pVolume, FLOAT* pDist );
BOOL IntersectRayOrientedBox( FXMVECTOR Origin, FXMVECTOR Direction, const OrientedBox* pVolume, FLOAT* pDist );
BOOL IntersectTriangleTriangle( FXMVECTOR A0, FXMVECTOR A1, FXMVECTOR A2, CXMVECTOR B0, CXMVECTOR B1, CXMVECTOR B2 );
BOOL IntersectTriangleSphere( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, const Sphere* pVolume );
BOOL IntersectTriangleAxisAlignedBox( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, const AxisAlignedBox* pVolume );
BOOL IntersectTriangleOrientedBox( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, const OrientedBox* pVolume );
BOOL IntersectSphereSphere( const Sphere* pVolumeA, const Sphere* pVolumeB );
BOOL IntersectSphereAxisAlignedBox( const Sphere* pVolumeA, const AxisAlignedBox* pVolumeB );
BOOL IntersectSphereOrientedBox( const Sphere* pVolumeA, const OrientedBox* pVolumeB );
BOOL IntersectAxisAlignedBoxAxisAlignedBox( const AxisAlignedBox* pVolumeA, const AxisAlignedBox* pVolumeB );
BOOL IntersectAxisAlignedBoxOrientedBox( const AxisAlignedBox* pVolumeA, const OrientedBox* pVolumeB );
BOOL IntersectOrientedBoxOrientedBox( const OrientedBox* pVolumeA, const OrientedBox* pVolumeB );



//-----------------------------------------------------------------------------
// Frustum intersection testing routines.
// Return values: 0 = no intersection, 
//                1 = intersection, 
//                2 = A is completely inside B
//-----------------------------------------------------------------------------
INT IntersectTriangleFrustum( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, const Frustum* pVolume );
INT IntersectSphereFrustum( const Sphere* pVolumeA, const Frustum* pVolumeB );
INT IntersectAxisAlignedBoxFrustum( const AxisAlignedBox* pVolumeA, const Frustum* pVolumeB );
INT IntersectOrientedBoxFrustum( const OrientedBox* pVolumeA, const Frustum* pVolumeB );
INT IntersectFrustumFrustum( const Frustum* pVolumeA, const Frustum* pVolumeB );




//-----------------------------------------------------------------------------
// Test vs six planes (usually forming a frustum) intersection routines.
// The intended use for these routines is for fast culling to a view frustum.  
// When the volume being tested against a view frustum is small relative to the
// view frustum it is usually either inside all six planes of the frustum or 
// outside one of the planes of the frustum. If neither of these cases is true
// then it may or may not be intersecting the frustum. Outside a plane is 
// defined as being on the positive side of the plane (and inside negative).
// Return values: 0 = volume is outside one of the planes (no intersection),
//                1 = not completely inside or completely outside (intersecting),
//                2 = volume is inside all the planes (completely inside)
//-----------------------------------------------------------------------------
INT IntersectTriangle6Planes( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, CXMVECTOR Plane0, CXMVECTOR Plane1,
                              CXMVECTOR Plane2, CXMVECTOR Plane3, CXMVECTOR Plane4, CXMVECTOR Plane5 );
INT IntersectSphere6Planes( const Sphere* pVolume, FXMVECTOR Plane0, FXMVECTOR Plane1, FXMVECTOR Plane2,
                            CXMVECTOR Plane3, CXMVECTOR Plane4, CXMVECTOR Plane5 );
INT IntersectAxisAlignedBox6Planes( const AxisAlignedBox* pVolume, FXMVECTOR Plane0, FXMVECTOR Plane1,
                                    FXMVECTOR Plane2, CXMVECTOR Plane3, CXMVECTOR Plane4, CXMVECTOR Plane5 );
INT IntersectOrientedBox6Planes( const OrientedBox* pVolume, FXMVECTOR Plane0, FXMVECTOR Plane1, FXMVECTOR Plane2,
                                 CXMVECTOR Plane3, CXMVECTOR Plane4, CXMVECTOR Plane5 );
INT IntersectFrustum6Planes( const Frustum* pVolume, FXMVECTOR Plane0, FXMVECTOR Plane1, FXMVECTOR Plane2,
                             CXMVECTOR Plane3, CXMVECTOR Plane4, CXMVECTOR Plane5 );


//-----------------------------------------------------------------------------
// Volume vs plane intersection testing routines.
// Return values: 0 = volume is outside the plane (on the positive sideof the plane),
//                1 = volume intersects the plane,
//                2 = volume is inside the plane (on the negative side of the plane) 
//-----------------------------------------------------------------------------
INT IntersectTrianglePlane( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, CXMVECTOR Plane );
INT IntersectSpherePlane( const Sphere* pVolume, FXMVECTOR Plane );
INT IntersectAxisAlignedBoxPlane( const AxisAlignedBox* pVolume, FXMVECTOR Plane );
INT IntersectOrientedBoxPlane( const OrientedBox* pVolume, FXMVECTOR Plane );
INT IntersectFrustumPlane( const Frustum* pVolume, FXMVECTOR Plane );

}; // namespace

#endif