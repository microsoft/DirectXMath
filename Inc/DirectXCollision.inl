//-------------------------------------------------------------------------------------
// DirectXCollision.inl -- C++ Collision Math library
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

XMGLOBALCONST XMVECTORF32 g_BoxOffset[8] =
{
    { -1.0f, -1.0f,  1.0f, 0.0f },
    {  1.0f, -1.0f,  1.0f, 0.0f },
    {  1.0f,  1.0f,  1.0f, 0.0f },
    { -1.0f,  1.0f,  1.0f, 0.0f },
    { -1.0f, -1.0f, -1.0f, 0.0f },
    {  1.0f, -1.0f, -1.0f, 0.0f },
    {  1.0f,  1.0f, -1.0f, 0.0f },
    { -1.0f,  1.0f, -1.0f, 0.0f },
};

XMGLOBALCONST XMVECTORF32 g_RayEpsilon = { 1e-20f, 1e-20f, 1e-20f, 1e-20f };
XMGLOBALCONST XMVECTORF32 g_RayNegEpsilon = { -1e-20f, -1e-20f, -1e-20f, -1e-20f };

namespace Internal
{

//-----------------------------------------------------------------------------
// Return true if any of the elements of a 3 vector are equal to 0xffffffff.
// Slightly more efficient than using XMVector3EqualInt.
//-----------------------------------------------------------------------------
inline bool XMVector3AnyTrue( _In_ FXMVECTOR V )
{
    XMVECTOR C;

    // Duplicate the fourth element from the first element.
    C = XMVectorSwizzle( V, 0, 1, 2, 0 );

    return XMComparisonAnyTrue( XMVector4EqualIntR( C, XMVectorTrueInt() ) );
}


//-----------------------------------------------------------------------------
// Return true if all of the elements of a 3 vector are equal to 0xffffffff.
// Slightly more efficient than using XMVector3EqualInt.
//-----------------------------------------------------------------------------
inline bool XMVector3AllTrue( _In_ FXMVECTOR V )
{
    XMVECTOR C;

    // Duplicate the fourth element from the first element.
    C = XMVectorSwizzle( V, 0, 1, 2, 0 );

    return XMComparisonAllTrue( XMVector4EqualIntR( C, XMVectorTrueInt() ) );
}

#if defined(_PREFAST) || !defined(NDEBUG)

XMGLOBALCONST XMVECTOR g_UnitVectorEpsilon = { 1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f };
XMGLOBALCONST XMVECTOR g_UnitQuaternionEpsilon = { 1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f };
XMGLOBALCONST XMVECTOR g_UnitPlaneEpsilon = { 1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f };

//-----------------------------------------------------------------------------
// Return true if the vector is a unit vector (length == 1).
//-----------------------------------------------------------------------------
inline bool XMVector3IsUnit( _In_ FXMVECTOR V )
{
    XMVECTOR Difference = XMVector3Length( V ) - XMVectorSplatOne();
    return XMVector4Less( XMVectorAbs( Difference ), g_UnitVectorEpsilon );
}

//-----------------------------------------------------------------------------
// Return true if the quaterion is a unit quaternion.
//-----------------------------------------------------------------------------
inline bool XMQuaternionIsUnit( _In_ FXMVECTOR Q )
{
    XMVECTOR Difference = XMVector4Length( Q ) - XMVectorSplatOne();
    return XMVector4Less( XMVectorAbs( Difference ), g_UnitQuaternionEpsilon );
}

//-----------------------------------------------------------------------------
// Return true if the plane is a unit plane.
//-----------------------------------------------------------------------------
inline bool XMPlaneIsUnit( _In_ FXMVECTOR Plane )
{
    XMVECTOR Difference = XMVector3Length( Plane ) - XMVectorSplatOne();
    return XMVector4Less( XMVectorAbs( Difference ), g_UnitPlaneEpsilon );
}

#endif // __PREFAST__ || !NDEBUG

//-----------------------------------------------------------------------------
inline XMVECTOR XMPlaneTransform( _In_ FXMVECTOR Plane, _In_ FXMVECTOR Rotation, _In_ FXMVECTOR Translation )
{
    XMVECTOR vNormal = XMVector3Rotate( Plane, Rotation );
    XMVECTOR vD = XMVectorSplatW( Plane ) - XMVector3Dot( vNormal, Translation );

    return XMVectorInsert( vNormal, vD, 0, 0, 0, 0, 1 );
}

//-----------------------------------------------------------------------------
// Return the point on the line segement (S1, S2) nearest the point P.
//-----------------------------------------------------------------------------
inline XMVECTOR PointOnLineSegmentNearestPoint( _In_ FXMVECTOR S1, _In_ FXMVECTOR S2, _In_ FXMVECTOR P )
{
    XMVECTOR Dir = S2 - S1;
    XMVECTOR Projection = ( XMVector3Dot( P, Dir ) - XMVector3Dot( S1, Dir ) );
    XMVECTOR LengthSq = XMVector3Dot( Dir, Dir );

    XMVECTOR t = Projection * XMVectorReciprocal( LengthSq );
    XMVECTOR Point = S1 + t * Dir;

    // t < 0
    XMVECTOR SelectS1 = XMVectorLess( Projection, XMVectorZero() );
    Point = XMVectorSelect( Point, S1, SelectS1 );

    // t > 1
    XMVECTOR SelectS2 = XMVectorGreater( Projection, LengthSq );
    Point = XMVectorSelect( Point, S2, SelectS2 );

    return Point;
}


//-----------------------------------------------------------------------------
// Test if the point (P) on the plane of the triangle is inside the triangle 
// (V0, V1, V2).
//-----------------------------------------------------------------------------
inline XMVECTOR PointOnPlaneInsideTriangle( _In_ FXMVECTOR P, _In_ FXMVECTOR V0, _In_ FXMVECTOR V1, _In_ CXMVECTOR V2 )
{
    // Compute the triangle normal.
    XMVECTOR N = XMVector3Cross( V2 - V0, V1 - V0 );

    // Compute the cross products of the vector from the base of each edge to 
    // the point with each edge vector.
    XMVECTOR C0 = XMVector3Cross( P - V0, V1 - V0 );
    XMVECTOR C1 = XMVector3Cross( P - V1, V2 - V1 );
    XMVECTOR C2 = XMVector3Cross( P - V2, V0 - V2 );

    // If the cross product points in the same direction as the normal the the
    // point is inside the edge (it is zero if is on the edge).
    XMVECTOR Zero = XMVectorZero();
    XMVECTOR Inside0 = XMVectorGreaterOrEqual( XMVector3Dot( C0, N ), Zero );
    XMVECTOR Inside1 = XMVectorGreaterOrEqual( XMVector3Dot( C1, N ), Zero );
    XMVECTOR Inside2 = XMVectorGreaterOrEqual( XMVector3Dot( C2, N ), Zero );

    // If the point inside all of the edges it is inside.
    return XMVectorAndInt( XMVectorAndInt( Inside0, Inside1 ), Inside2 );
}

//-----------------------------------------------------------------------------
inline void FastIntersectSpherePlane( _In_ FXMVECTOR Center, _In_ FXMVECTOR Radius, _In_ FXMVECTOR Plane,
                                      _Out_ XMVECTOR& Outside, _Out_ XMVECTOR& Inside )
{
    XMVECTOR Dist = XMVector4Dot( Center, Plane );

    // Outside the plane?
    Outside = XMVectorGreater( Dist, Radius );

    // Fully inside the plane?
    Inside = XMVectorLess( Dist, -Radius );
}


//-----------------------------------------------------------------------------
inline void FastIntersectAxisAlignedBoxPlane( _In_ FXMVECTOR Center, _In_ FXMVECTOR Extents, _In_ FXMVECTOR Plane,
                                              _Out_ XMVECTOR& Outside, _Out_ XMVECTOR& Inside )
{
    // Compute the distance to the center of the box.
    XMVECTOR Dist = XMVector4Dot( Center, Plane );

    // Project the axes of the box onto the normal of the plane.  Half the
    // length of the projection (sometime called the "radius") is equal to
    // h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))
    // where h(i) are extents of the box, n is the plane normal, and b(i) are the 
    // axes of the box. In this case b(i) = [(1,0,0), (0,1,0), (0,0,1)].
    XMVECTOR Radius = XMVector3Dot( Extents, XMVectorAbs( Plane ) );

    // Outside the plane?
    Outside = XMVectorGreater( Dist, Radius );

    // Fully inside the plane?
    Inside = XMVectorLess( Dist, -Radius );
}

}; // namespace Internal


/****************************************************************************
 *
 * BoundingSphere
 *
 ****************************************************************************/

//-----------------------------------------------------------------------------
// Transform a sphere by an angle preserving transform.
//-----------------------------------------------------------------------------
inline void BoundingSphere::Transform( BoundingSphere& Out, CXMMATRIX M ) const
{
    // Load the center of the sphere.
    XMVECTOR vCenter = XMLoadFloat3( &Center );

    // Transform the center of the sphere.
    XMVECTOR C = XMVector3Transform( vCenter, M );
    
    XMVECTOR dX = XMVector3Dot( M.r[0], M.r[0] );
    XMVECTOR dY = XMVector3Dot( M.r[1], M.r[1] );
    XMVECTOR dZ = XMVector3Dot( M.r[2], M.r[2] );

    XMVECTOR d = XMVectorMax( dX, XMVectorMax( dY, dZ ) );

    // Store the center sphere.
    XMStoreFloat3( &Out.Center, C );

    // Scale the radius of the pshere.
    float Scale = sqrtf( XMVectorGetX(d) );
    Out.Radius = Radius * Scale;
}

inline void BoundingSphere::Transform( BoundingSphere& Out, float Scale, FXMVECTOR Rotation, FXMVECTOR Translation ) const
{
    // Load the center of the sphere.
    XMVECTOR vCenter = XMLoadFloat3( &Center );

    // Transform the center of the sphere.
    vCenter = XMVector3Rotate( vCenter * XMVectorReplicate( Scale ), Rotation ) + Translation;

    // Store the center sphere.
    XMStoreFloat3( &Out.Center, vCenter );

    // Scale the radius of the pshere.
    Out.Radius = Radius * Scale;
}


//-----------------------------------------------------------------------------
// Point in sphere test.
//-----------------------------------------------------------------------------
inline ContainmentType BoundingSphere::Contains( FXMVECTOR Point ) const
{
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vRadius = XMVectorReplicatePtr( &Radius );

    XMVECTOR DistanceSquared = XMVector3LengthSq( Point - vCenter );
    XMVECTOR RadiusSquared = XMVectorMultiply( vRadius, vRadius );

    return XMVector3LessOrEqual( DistanceSquared, RadiusSquared ) ? CONTAINS : DISJOINT;
}


//-----------------------------------------------------------------------------
// Triangle in sphere test
//-----------------------------------------------------------------------------
inline ContainmentType BoundingSphere::Contains( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2 ) const
{
    if ( !Intersects(V0,V1,V2) )
        return DISJOINT;

    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vRadius = XMVectorReplicatePtr( &Radius );
    XMVECTOR RadiusSquared = XMVectorMultiply( vRadius, vRadius );

    XMVECTOR DistanceSquared = XMVector3LengthSq( V0 - vCenter );
    XMVECTOR Inside = XMVectorLessOrEqual(DistanceSquared, RadiusSquared);

    DistanceSquared = XMVector3LengthSq( V1 - vCenter );
    Inside = XMVectorAndInt( Inside, XMVectorLessOrEqual(DistanceSquared, RadiusSquared) );

    DistanceSquared = XMVector3LengthSq( V2 - vCenter );
    Inside = XMVectorAndInt( Inside, XMVectorLessOrEqual(DistanceSquared, RadiusSquared) );

    return ( XMVector3EqualInt( Inside, XMVectorTrueInt() ) ) ? CONTAINS : INTERSECTS;
}


//-----------------------------------------------------------------------------
// Sphere in sphere test.
//-----------------------------------------------------------------------------
inline ContainmentType BoundingSphere::Contains( const BoundingSphere& sh ) const
{
    XMVECTOR Center1 = XMLoadFloat3( &Center );
    float r1 = Radius;

    XMVECTOR Center2 = XMLoadFloat3( &sh.Center );
    float r2 = sh.Radius;

    XMVECTOR V = XMVectorSubtract( Center2, Center1 );

    XMVECTOR Dist = XMVector3Length( V );

    float d = XMVectorGetX( Dist );

    return (r1 + r2 >= d) ? ((r1 - r2 >= d) ? CONTAINS : INTERSECTS) : DISJOINT;
}


//-----------------------------------------------------------------------------
// Axis-aligned box in sphere test
//-----------------------------------------------------------------------------
inline ContainmentType BoundingSphere::Contains( const BoundingBox& box ) const
{
    if ( !box.Intersects(*this) )
        return DISJOINT;

    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vRadius = XMVectorReplicatePtr( &Radius );
    XMVECTOR RadiusSq = vRadius * vRadius;

    XMVECTOR boxCenter = XMLoadFloat3( &box.Center );
    XMVECTOR boxExtents = XMLoadFloat3( &box.Extents );

    XMVECTOR InsideAll = XMVectorTrueInt();

    XMVECTOR offset = boxCenter - vCenter;

    for( size_t i = 0; i < BoundingBox::CORNER_COUNT; ++i )
    {
        XMVECTOR C = XMVectorMultiplyAdd( boxExtents, g_BoxOffset[i], offset );
        XMVECTOR d = XMVector3LengthSq( C );
        InsideAll = XMVectorAndInt( InsideAll, XMVectorLessOrEqual( d, RadiusSq ) );
    }

    return ( XMVector3EqualInt( InsideAll, XMVectorTrueInt() ) ) ? CONTAINS : INTERSECTS;
}


//-----------------------------------------------------------------------------
// Sphere vs. sphere test.
//-----------------------------------------------------------------------------
inline bool BoundingSphere::Intersects( const BoundingSphere& sh ) const
{
    // Load A.
    XMVECTOR vCenterA = XMLoadFloat3( &Center );
    XMVECTOR vRadiusA = XMVectorReplicatePtr( &Radius );

    // Load B.
    XMVECTOR vCenterB = XMLoadFloat3( &sh.Center );
    XMVECTOR vRadiusB = XMVectorReplicatePtr( &sh.Radius );

    // Distance squared between centers.    
    XMVECTOR Delta = vCenterB - vCenterA;
    XMVECTOR DistanceSquared = XMVector3LengthSq( Delta );

    // Sum of the radii squared.
    XMVECTOR RadiusSquared = XMVectorAdd( vRadiusA, vRadiusB );
    RadiusSquared = XMVectorMultiply( RadiusSquared, RadiusSquared );

    return XMVector3LessOrEqual( DistanceSquared, RadiusSquared );
}


//-----------------------------------------------------------------------------
// Box vs. sphere test.
//-----------------------------------------------------------------------------
inline bool BoundingSphere::Intersects( const BoundingBox& box ) const
{
    return box.Intersects( *this );
}


//-----------------------------------------------------------------------------
// Triangle vs sphere test
//-----------------------------------------------------------------------------
inline bool BoundingSphere::Intersects( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2 ) const
{
    // Load the sphere.    
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vRadius = XMVectorReplicatePtr( &Radius );

    // Compute the plane of the triangle (has to be normalized).
    XMVECTOR N = XMVector3Normalize( XMVector3Cross( V1 - V0, V2 - V0 ) );

    // Assert that the triangle is not degenerate.
    assert( !XMVector3Equal( N, XMVectorZero() ) );

    // Find the nearest feature on the triangle to the sphere.
    XMVECTOR Dist = XMVector3Dot( vCenter - V0, N );

    // If the center of the sphere is farther from the plane of the triangle than
    // the radius of the sphere, then there cannot be an intersection.
    XMVECTOR NoIntersection = XMVectorLess( Dist, -vRadius );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Dist, vRadius ) );

    // Project the center of the sphere onto the plane of the triangle.
    XMVECTOR Point = vCenter - ( N * Dist );

    // Is it inside all the edges? If so we intersect because the distance 
    // to the plane is less than the radius.
    XMVECTOR Intersection = DirectX::Internal::PointOnPlaneInsideTriangle( Point, V0, V1, V2 );

    // Find the nearest point on each edge.
    XMVECTOR RadiusSq = vRadius * vRadius;

    // Edge 0,1
    Point = DirectX::Internal::PointOnLineSegmentNearestPoint( V0, V1, vCenter );

    // If the distance to the center of the sphere to the point is less than 
    // the radius of the sphere then it must intersect.
    Intersection = XMVectorOrInt( Intersection, XMVectorLessOrEqual( XMVector3LengthSq( vCenter - Point ), RadiusSq ) );

    // Edge 1,2
    Point = DirectX::Internal::PointOnLineSegmentNearestPoint( V1, V2, vCenter );

    // If the distance to the center of the sphere to the point is less than 
    // the radius of the sphere then it must intersect.
    Intersection = XMVectorOrInt( Intersection, XMVectorLessOrEqual( XMVector3LengthSq( vCenter - Point ), RadiusSq ) );

    // Edge 2,0
    Point = DirectX::Internal::PointOnLineSegmentNearestPoint( V2, V0, vCenter );

    // If the distance to the center of the sphere to the point is less than 
    // the radius of the sphere then it must intersect.
    Intersection = XMVectorOrInt( Intersection, XMVectorLessOrEqual( XMVector3LengthSq( vCenter - Point ), RadiusSq ) );

    return XMVector4EqualInt( XMVectorAndCInt( Intersection, NoIntersection ), XMVectorTrueInt() );
}


//-----------------------------------------------------------------------------
// Sphere-plane intersection
//-----------------------------------------------------------------------------
inline PlaneIntersectionType BoundingSphere::Intersects( FXMVECTOR Plane ) const
{
    assert( DirectX::Internal::XMPlaneIsUnit( Plane ) );

    // Load the sphere.
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vRadius = XMVectorReplicatePtr( &Radius );

    // Set w of the center to one so we can dot4 with a plane.
    vCenter = XMVectorInsert( vCenter, XMVectorSplatOne(), 0, 0, 0, 0, 1 );

    XMVECTOR Outside, Inside;
    DirectX::Internal::FastIntersectSpherePlane( vCenter, vRadius, Plane, Outside, Inside );

    // If the sphere is outside any plane it is outside.
    if ( XMVector4EqualInt( Outside, XMVectorTrueInt() ) )
        return FRONT;

    // If the sphere is inside all planes it is inside.
    if ( XMVector4EqualInt( Inside, XMVectorTrueInt() ) )
        return BACK;

    // The sphere is not inside all planes or outside a plane it intersects.
    return INTERSECTING;
}


//-----------------------------------------------------------------------------
// Compute the intersection of a ray (Origin, Direction) with a sphere.
//-----------------------------------------------------------------------------
inline bool BoundingSphere::Intersects( FXMVECTOR Origin, FXMVECTOR Direction, float& Dist ) const
{
    assert( DirectX::Internal::XMVector3IsUnit( Direction ) );

    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vRadius = XMVectorReplicatePtr( &Radius );

    // l is the vector from the ray origin to the center of the sphere.
    XMVECTOR l = vCenter - Origin;

    // s is the projection of the l onto the ray direction.
    XMVECTOR s = XMVector3Dot( l, Direction );

    XMVECTOR l2 = XMVector3Dot( l, l );

    XMVECTOR r2 = vRadius * vRadius;

    // m2 is squared distance from the center of the sphere to the projection.
    XMVECTOR m2 = l2 - s * s;

    XMVECTOR NoIntersection;

    // If the ray origin is outside the sphere and the center of the sphere is 
    // behind the ray origin there is no intersection.
    NoIntersection = XMVectorAndInt( XMVectorLess( s, XMVectorZero() ), XMVectorGreater( l2, r2 ) );

    // If the squared distance from the center of the sphere to the projection
    // is greater than the radius squared the ray will miss the sphere.
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( m2, r2 ) );

    // The ray hits the sphere, compute the nearest intersection point.
    XMVECTOR q = XMVectorSqrt( r2 - m2 );
    XMVECTOR t1 = s - q;
    XMVECTOR t2 = s + q;

    XMVECTOR OriginInside = XMVectorLessOrEqual( l2, r2 );
    XMVECTOR t = XMVectorSelect( t1, t2, OriginInside );

    if( XMVector4NotEqualInt( NoIntersection, XMVectorTrueInt() ) )
    {
        // Store the x-component to *pDist.
        XMStoreFloat( &Dist, t );
        return true;
    }

    return false;
}


//-----------------------------------------------------------------------------
// Creates a bounding sphere that contains two other bounding spheres
//-----------------------------------------------------------------------------
inline void BoundingSphere::CreateMerged( BoundingSphere& Out, const BoundingSphere& S1, const BoundingSphere& S2 )
{
    XMVECTOR Center1 = XMLoadFloat3( &S1.Center );
    float r1 = S1.Radius;

    XMVECTOR Center2 = XMLoadFloat3( &S2.Center );
    float r2 = S2.Radius;

    XMVECTOR V = XMVectorSubtract( Center2, Center1 );

    XMVECTOR Dist = XMVector3Length( V );

    float d = XMVectorGetX(Dist);

    if ( r1 + r2 >= d )
    {
        if ( r1 - r2 >= d )
        {
            Out = S1;
            return;
        }
        else if ( r2 - r1 >= d )
        {
            Out = S2;
            return;
        }
    }

    XMVECTOR N = XMVectorDivide( V, Dist );

    float t1 = XMMin( -r1, d-r2 );
    float t2 = XMMax( r1, d+r2 );
    float t_5 = (t2 - t1) * 0.5f;
    
    XMVECTOR NCenter = XMVectorAdd( Center1, XMVectorMultiply( N, XMVectorReplicate( t_5 + t1 ) ) );

    XMStoreFloat3( &Out.Center, NCenter );
    Out.Radius = t_5;
}


//-----------------------------------------------------------------------------
// Create sphere enscribing bounding box
//-----------------------------------------------------------------------------
inline void BoundingSphere::CreateFromBoundingBox( BoundingSphere& Out, const BoundingBox& box )
{
    Out.Center = box.Center;
    XMVECTOR vExtents = XMLoadFloat3( &box.Extents );
    Out.Radius = XMVectorGetX( XMVector3Length( vExtents ) );
}


//-----------------------------------------------------------------------------
// Find the approximate smallest enclosing bounding sphere for a set of 
// points. Exact computation of the smallest enclosing bounding sphere is 
// possible but is slower and requires a more complex algorithm.
// The algorithm is based on  Jack Ritter, "An Efficient Bounding Sphere", 
// Graphics Gems.
//-----------------------------------------------------------------------------
inline void BoundingSphere::CreateFromPoints( BoundingSphere& Out, size_t Count, const XMFLOAT3* pPoints, size_t Stride )
{
    assert( Count > 0 );
    assert( pPoints );

    // Find the points with minimum and maximum x, y, and z
    XMVECTOR MinX, MaxX, MinY, MaxY, MinZ, MaxZ;

    MinX = MaxX = MinY = MaxY = MinZ = MaxZ = XMLoadFloat3( pPoints );

    for( size_t i = 1; i < Count; ++i )
    {
        XMVECTOR Point = XMLoadFloat3( reinterpret_cast<const XMFLOAT3*>( reinterpret_cast<const uint8_t*>(pPoints) + i * Stride ) );

        float px = XMVectorGetX( Point );
        float py = XMVectorGetY( Point );
        float pz = XMVectorGetZ( Point );

        if( px < XMVectorGetX( MinX ) )
            MinX = Point;

        if( px > XMVectorGetX( MaxX ) )
            MaxX = Point;

        if( py < XMVectorGetY( MinY ) )
            MinY = Point;

        if( py > XMVectorGetY( MaxY ) )
            MaxY = Point;

        if( pz < XMVectorGetZ( MinZ ) )
            MinZ = Point;

        if( pz > XMVectorGetZ( MaxZ ) )
            MaxZ = Point;
    }

    // Use the min/max pair that are farthest apart to form the initial sphere.
    XMVECTOR DeltaX = MaxX - MinX;
    XMVECTOR DistX = XMVector3Length( DeltaX );

    XMVECTOR DeltaY = MaxY - MinY;
    XMVECTOR DistY = XMVector3Length( DeltaY );

    XMVECTOR DeltaZ = MaxZ - MinZ;
    XMVECTOR DistZ = XMVector3Length( DeltaZ );

    XMVECTOR vCenter;
    XMVECTOR vRadius;

    if( XMVector3Greater( DistX, DistY ) )
    {
        if( XMVector3Greater( DistX, DistZ ) )
        {
            // Use min/max x.
            vCenter = XMVectorLerp(MaxX,MinX,0.5f);
            vRadius = DistX * 0.5f;
        }
        else
        {
            // Use min/max z.
            vCenter = XMVectorLerp(MaxZ,MinZ,0.5f);
            vRadius = DistZ * 0.5f;
        }
    }
    else // Y >= X
    {
        if( XMVector3Greater( DistY, DistZ ) )
        {
            // Use min/max y.
            vCenter = XMVectorLerp(MaxY,MinY,0.5f);
            vRadius = DistY * 0.5f;
        }
        else
        {
            // Use min/max z.
            vCenter = XMVectorLerp(MaxZ,MinZ,0.5f);
            vRadius = DistZ * 0.5f;
        }
    }

    // Add any points not inside the sphere.
    for( size_t i = 0; i < Count; ++i )
    {
        XMVECTOR Point = XMLoadFloat3( reinterpret_cast<const XMFLOAT3*>( reinterpret_cast<const uint8_t*>(pPoints) + i * Stride ) );

        XMVECTOR Delta = Point - vCenter;

        XMVECTOR Dist = XMVector3Length( Delta );

        if( XMVector3Greater( Dist, vRadius ) )
        {
            // Adjust sphere to include the new point.
            vRadius = ( vRadius + Dist ) * 0.5f;
            vCenter += ( XMVectorReplicate( 1.0f ) - XMVectorDivide(vRadius,Dist) ) * Delta;
        }
    }

    XMStoreFloat3( &Out.Center, vCenter );
    XMStoreFloat( &Out.Radius, vRadius );
}


/****************************************************************************
 *
 * BoundingBox
 *
 ****************************************************************************/

//-----------------------------------------------------------------------------
// Transform an axis aligned box by an angle preserving transform.
//-----------------------------------------------------------------------------
inline void BoundingBox::Transform( BoundingBox& Out, CXMMATRIX M ) const
{
    // Load center and extents.
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    // Compute and transform the corners and find new min/max bounds.
    XMVECTOR Corner = XMVectorMultiplyAdd( vExtents, g_BoxOffset[0], vCenter );
    Corner = XMVector3Transform( Corner, M );

    XMVECTOR Min, Max;
    Min = Max = Corner;

    for( size_t i = 1; i < CORNER_COUNT; ++i )
    {
        Corner = XMVectorMultiplyAdd( vExtents, g_BoxOffset[i], vCenter );
        Corner = XMVector3Transform( Corner, M );

        Min = XMVectorMin( Min, Corner );
        Max = XMVectorMax( Max, Corner );
    }

    // Store center and extents.
    XMStoreFloat3( &Out.Center, ( Min + Max ) * 0.5f );
    XMStoreFloat3( &Out.Extents, ( Max - Min ) * 0.5f );
}

inline void BoundingBox::Transform( BoundingBox& Out, float Scale, FXMVECTOR Rotation, FXMVECTOR Translation ) const
{
    assert( DirectX::Internal::XMQuaternionIsUnit( Rotation ) );

    // Load center and extents.
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    XMVECTOR VectorScale = XMVectorReplicate( Scale );

    // Compute and transform the corners and find new min/max bounds.
    XMVECTOR Corner = XMVectorMultiplyAdd( vExtents, g_BoxOffset[0], vCenter );
    Corner = XMVector3Rotate( Corner * VectorScale, Rotation ) + Translation;

    XMVECTOR Min, Max;
    Min = Max = Corner;

    for( size_t i = 1; i < CORNER_COUNT; ++i )
    {
        Corner = XMVectorMultiplyAdd( vExtents, g_BoxOffset[i], vCenter );
        Corner = XMVector3Rotate( Corner * VectorScale, Rotation ) + Translation;

        Min = XMVectorMin( Min, Corner );
        Max = XMVectorMax( Max, Corner );
    }

    // Store center and extents.
    XMStoreFloat3( &Out.Center, ( Min + Max ) * 0.5f );
    XMStoreFloat3( &Out.Extents, ( Max - Min ) * 0.5f );
}


//-----------------------------------------------------------------------------
// Get the corner points of the box
//-----------------------------------------------------------------------------
inline void BoundingBox::GetCorners( XMFLOAT3* Corners ) const
{
    assert( Corners != nullptr );
    assert( CORNER_COUNT == 8 );

    // Load the box
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    for( size_t i = 0; i < CORNER_COUNT; ++i )
    {
        XMVECTOR C = XMVectorMultiplyAdd( vExtents, g_BoxOffset[i], vCenter );
        XMStoreFloat3( &Corners[i], C );
    }
}


//-----------------------------------------------------------------------------
// Point in axis-aligned box test
//-----------------------------------------------------------------------------
inline ContainmentType BoundingBox::Contains( FXMVECTOR Point ) const
{
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    return XMVector3InBounds( Point - vCenter, vExtents ) ? CONTAINS : DISJOINT;
}


//-----------------------------------------------------------------------------
// Triangle in axis-aligned box test
//-----------------------------------------------------------------------------
inline ContainmentType BoundingBox::Contains( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2 ) const
{
    if ( !Intersects(V0,V1,V2) )
        return DISJOINT;

    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    XMVECTOR d = XMVector3LengthSq( V0 - vCenter );
    XMVECTOR Inside = XMVectorLessOrEqual( d, vExtents );

    d = XMVector3LengthSq( V1 - vCenter );
    Inside = XMVectorAndInt( Inside, XMVectorLessOrEqual( d, vExtents ) );

    d = XMVector3LengthSq( V2 - vCenter );
    Inside = XMVectorAndInt( Inside, XMVectorLessOrEqual( d, vExtents ) );

    return ( XMVector3EqualInt( Inside, XMVectorTrueInt() ) ) ? CONTAINS : INTERSECTS;
}


//-----------------------------------------------------------------------------
// Sphere in axis-aligned box test
//-----------------------------------------------------------------------------
inline ContainmentType BoundingBox::Contains( const BoundingSphere& sh ) const
{
    XMVECTOR SphereCenter = XMLoadFloat3( &sh.Center );
    XMVECTOR SphereRadius = XMVectorReplicatePtr( &sh.Radius );

    XMVECTOR BoxCenter = XMLoadFloat3( &Center );
    XMVECTOR BoxExtents = XMLoadFloat3( &Extents );

    XMVECTOR BoxMin = BoxCenter - BoxExtents;
    XMVECTOR BoxMax = BoxCenter + BoxExtents;

    // Find the distance to the nearest point on the box.
    // for each i in (x, y, z)
    // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
    // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

    XMVECTOR d = XMVectorZero();

    // Compute d for each dimension.
    XMVECTOR LessThanMin = XMVectorLess( SphereCenter, BoxMin );
    XMVECTOR GreaterThanMax = XMVectorGreater( SphereCenter, BoxMax );

    XMVECTOR MinDelta = SphereCenter - BoxMin;
    XMVECTOR MaxDelta = SphereCenter - BoxMax;

    // Choose value for each dimension based on the comparison.
    d = XMVectorSelect( d, MinDelta, LessThanMin );
    d = XMVectorSelect( d, MaxDelta, GreaterThanMax );

    // Use a dot-product to square them and sum them together.
    XMVECTOR d2 = XMVector3Dot( d, d );

    if ( XMVector3Greater( d2, XMVectorMultiply( SphereRadius, SphereRadius ) ) )
        return DISJOINT;

    XMVECTOR InsideAll = XMVectorLessOrEqual( BoxMin + SphereRadius, SphereCenter );
    InsideAll = XMVectorAndInt( InsideAll, XMVectorLessOrEqual( SphereCenter, BoxMax - SphereRadius ) );
    InsideAll = XMVectorAndInt( InsideAll, XMVectorGreater( BoxMax - BoxMin, SphereRadius ) );

    return ( XMVector3EqualInt( InsideAll, XMVectorTrueInt() ) ) ? CONTAINS : INTERSECTS;
}


//-----------------------------------------------------------------------------
// Axis-aligned box in axis-aligned box test
//-----------------------------------------------------------------------------
inline ContainmentType BoundingBox::Contains( const BoundingBox& box ) const
{
    XMVECTOR CenterA = XMLoadFloat3( &Center );
    XMVECTOR ExtentsA = XMLoadFloat3( &Extents );

    XMVECTOR CenterB = XMLoadFloat3( &box.Center );
    XMVECTOR ExtentsB = XMLoadFloat3( &box.Extents );

    XMVECTOR MinA = CenterA - ExtentsA;
    XMVECTOR MaxA = CenterA + ExtentsA;

    XMVECTOR MinB = CenterB - ExtentsB;
    XMVECTOR MaxB = CenterB + ExtentsB;

    // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then return false
    XMVECTOR Disjoint = XMVectorOrInt( XMVectorGreater( MinA, MaxB ), XMVectorGreater( MinB, MaxA ) );

    if ( DirectX::Internal::XMVector3AnyTrue( Disjoint ) )
        return DISJOINT;

    // for each i in (x, y, z) if a_min(i) <= b_min(i) and b_max(i) <= a_max(i) then A contains B
    XMVECTOR Inside = XMVectorAndInt( XMVectorLessOrEqual( MinA, MinB ), XMVectorLessOrEqual( MaxB, MaxA ) );

    return DirectX::Internal::XMVector3AllTrue( Inside ) ? CONTAINS : INTERSECTS;
}


//-----------------------------------------------------------------------------
// Sphere vs axis-aligned box test
//-----------------------------------------------------------------------------
inline bool BoundingBox::Intersects( const BoundingSphere& sh ) const
{
    XMVECTOR SphereCenter = XMLoadFloat3( &sh.Center );
    XMVECTOR SphereRadius = XMVectorReplicatePtr( &sh.Radius );

    XMVECTOR BoxCenter = XMLoadFloat3( &Center );
    XMVECTOR BoxExtents = XMLoadFloat3( &Extents );

    XMVECTOR BoxMin = BoxCenter - BoxExtents;
    XMVECTOR BoxMax = BoxCenter + BoxExtents;

    // Find the distance to the nearest point on the box.
    // for each i in (x, y, z)
    // if (SphereCenter(i) < BoxMin(i)) d2 += (SphereCenter(i) - BoxMin(i)) ^ 2
    // else if (SphereCenter(i) > BoxMax(i)) d2 += (SphereCenter(i) - BoxMax(i)) ^ 2

    XMVECTOR d = XMVectorZero();

    // Compute d for each dimension.
    XMVECTOR LessThanMin = XMVectorLess( SphereCenter, BoxMin );
    XMVECTOR GreaterThanMax = XMVectorGreater( SphereCenter, BoxMax );

    XMVECTOR MinDelta = SphereCenter - BoxMin;
    XMVECTOR MaxDelta = SphereCenter - BoxMax;

    // Choose value for each dimension based on the comparison.
    d = XMVectorSelect( d, MinDelta, LessThanMin );
    d = XMVectorSelect( d, MaxDelta, GreaterThanMax );

    // Use a dot-product to square them and sum them together.
    XMVECTOR d2 = XMVector3Dot( d, d );

    return XMVector3LessOrEqual( d2, XMVectorMultiply( SphereRadius, SphereRadius ) );
}


//-----------------------------------------------------------------------------
// Axis-aligned box vs. axis-aligned box test
//-----------------------------------------------------------------------------
inline bool BoundingBox::Intersects( const BoundingBox& box ) const
{
    XMVECTOR CenterA = XMLoadFloat3( &Center );
    XMVECTOR ExtentsA = XMLoadFloat3( &Extents );

    XMVECTOR CenterB = XMLoadFloat3( &box.Center );
    XMVECTOR ExtentsB = XMLoadFloat3( &box.Extents );

    XMVECTOR MinA = CenterA - ExtentsA;
    XMVECTOR MaxA = CenterA + ExtentsA;

    XMVECTOR MinB = CenterB - ExtentsB;
    XMVECTOR MaxB = CenterB + ExtentsB;

    // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then return false
    XMVECTOR Disjoint = XMVectorOrInt( XMVectorGreater( MinA, MaxB ), XMVectorGreater( MinB, MaxA ) );

    return !DirectX::Internal::XMVector3AnyTrue( Disjoint );
}


//-----------------------------------------------------------------------------
// Triangle vs. axis aligned box test
//-----------------------------------------------------------------------------
inline bool BoundingBox::Intersects( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2 ) const
{
    static const XMVECTORI32 Permute0W1Z0Y0X =
    {
        XM_PERMUTE_0W, XM_PERMUTE_1Z, XM_PERMUTE_0Y, XM_PERMUTE_0X
    };
    static const XMVECTORI32 Permute0Z0W1X0Y =
    {
        XM_PERMUTE_0Z, XM_PERMUTE_0W, XM_PERMUTE_1X, XM_PERMUTE_0Y
    };
    static const XMVECTORI32 Permute1Y0X0W0Z =
    {
        XM_PERMUTE_1Y, XM_PERMUTE_0X, XM_PERMUTE_0W, XM_PERMUTE_0Z
    };

    XMVECTOR Zero = XMVectorZero();

    // Load the box.
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    XMVECTOR BoxMin = vCenter - vExtents;
    XMVECTOR BoxMax = vCenter + vExtents;

    // Test the axes of the box (in effect test the AAB against the minimal AAB 
    // around the triangle).
    XMVECTOR TriMin = XMVectorMin( XMVectorMin( V0, V1 ), V2 );
    XMVECTOR TriMax = XMVectorMax( XMVectorMax( V0, V1 ), V2 );

    // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then disjoint
    XMVECTOR Disjoint = XMVectorOrInt( XMVectorGreater( TriMin, BoxMax ), XMVectorGreater( BoxMin, TriMax ) );
    if( DirectX::Internal::XMVector3AnyTrue( Disjoint ) )
        return false;

    // Test the plane of the triangle.
    XMVECTOR Normal = XMVector3Cross( V1 - V0, V2 - V0 );
    XMVECTOR Dist = XMVector3Dot( Normal, V0 );

    // Assert that the triangle is not degenerate.
    assert( !XMVector3Equal( Normal, Zero ) );

    // for each i in (x, y, z) if n(i) >= 0 then v_min(i)=b_min(i), v_max(i)=b_max(i)
    // else v_min(i)=b_max(i), v_max(i)=b_min(i)
    XMVECTOR NormalSelect = XMVectorGreater( Normal, Zero );
    XMVECTOR V_Min = XMVectorSelect( BoxMax, BoxMin, NormalSelect );
    XMVECTOR V_Max = XMVectorSelect( BoxMin, BoxMax, NormalSelect );

    // if n dot v_min + d > 0 || n dot v_max + d < 0 then disjoint
    XMVECTOR MinDist = XMVector3Dot( V_Min, Normal );
    XMVECTOR MaxDist = XMVector3Dot( V_Max, Normal );

    XMVECTOR NoIntersection = XMVectorGreater( MinDist, Dist );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( MaxDist, Dist ) );

    // Move the box center to zero to simplify the following tests.
    XMVECTOR TV0 = V0 - vCenter;
    XMVECTOR TV1 = V1 - vCenter;
    XMVECTOR TV2 = V2 - vCenter;

    // Test the edge/edge axes (3*3).
    XMVECTOR e0 = TV1 - TV0;
    XMVECTOR e1 = TV2 - TV1;
    XMVECTOR e2 = TV0 - TV2;

    // Make w zero.
    e0 = XMVectorInsert( e0, Zero, 0, 0, 0, 0, 1 );
    e1 = XMVectorInsert( e1, Zero, 0, 0, 0, 0, 1 );
    e2 = XMVectorInsert( e2, Zero, 0, 0, 0, 0, 1 );

    XMVECTOR Axis;
    XMVECTOR p0, p1, p2;
    XMVECTOR Min, Max;
    XMVECTOR Radius;

    // Axis == (1,0,0) x e0 = (0, -e0.z, e0.y)
    Axis = XMVectorPermute( e0, -e0, Permute0W1Z0Y0X );
    p0 = XMVector3Dot( TV0, Axis );
    // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
    p2 = XMVector3Dot( TV2, Axis );
    Min = XMVectorMin( p0, p2 );
    Max = XMVectorMax( p0, p2 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (1,0,0) x e1 = (0, -e1.z, e1.y)
    Axis = XMVectorPermute( e1, -e1, Permute0W1Z0Y0X );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (1,0,0) x e2 = (0, -e2.z, e2.y)
    Axis = XMVectorPermute( e2, -e2, Permute0W1Z0Y0X );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,1,0) x e0 = (e0.z, 0, -e0.x)
    Axis = XMVectorPermute( e0, -e0, Permute0Z0W1X0Y );
    p0 = XMVector3Dot( TV0, Axis );
    // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
    p2 = XMVector3Dot( TV2, Axis );
    Min = XMVectorMin( p0, p2 );
    Max = XMVectorMax( p0, p2 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,1,0) x e1 = (e1.z, 0, -e1.x)
    Axis = XMVectorPermute( e1, -e1, Permute0Z0W1X0Y );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e2 = (e2.z, 0, -e2.x)
    Axis = XMVectorPermute( e2, -e2, Permute0Z0W1X0Y );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e0 = (-e0.y, e0.x, 0)
    Axis = XMVectorPermute( e0, -e0, Permute1Y0X0W0Z );
    p0 = XMVector3Dot( TV0, Axis );
    // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
    p2 = XMVector3Dot( TV2, Axis );
    Min = XMVectorMin( p0, p2 );
    Max = XMVectorMax( p0, p2 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e1 = (-e1.y, e1.x, 0)
    Axis = XMVectorPermute( e1, -e1, Permute1Y0X0W0Z );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e2 = (-e2.y, e2.x, 0)
    Axis = XMVectorPermute( e2, -e2, Permute1Y0X0W0Z );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( vExtents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    return XMVector4NotEqualInt( NoIntersection, XMVectorTrueInt() );
}


//-----------------------------------------------------------------------------
inline PlaneIntersectionType BoundingBox::Intersects( FXMVECTOR Plane ) const
{
    assert( DirectX::Internal::XMPlaneIsUnit( Plane ) );

    // Load the box.
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    // Set w of the center to one so we can dot4 with a plane.
    vCenter = XMVectorInsert( vCenter, XMVectorSplatOne(), 0, 0, 0, 0, 1);

    XMVECTOR Outside, Inside;
    DirectX::Internal::FastIntersectAxisAlignedBoxPlane( vCenter, vExtents, Plane, Outside, Inside );

    // If the box is outside any plane it is outside.
    if ( XMVector4EqualInt( Outside, XMVectorTrueInt() ) )
        return FRONT;

    // If the box is inside all planes it is inside.
    if ( XMVector4EqualInt( Inside, XMVectorTrueInt() ) )
        return BACK;

    // The box is not inside all planes or outside a plane it intersects.
    return INTERSECTING;
}


//-----------------------------------------------------------------------------
// Compute the intersection of a ray (Origin, Direction) with an axis aligned 
// box using the slabs method.
//-----------------------------------------------------------------------------
inline bool BoundingBox::Intersects( FXMVECTOR Origin, FXMVECTOR Direction, float& Dist ) const
{
    assert( DirectX::Internal::XMVector3IsUnit( Direction ) );

    static const XMVECTOR FltMin =
    {
        -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX
    };
    static const XMVECTOR FltMax =
    {
        FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX
    };

    // Load the box.
    XMVECTOR vCenter = XMLoadFloat3( &Center );
    XMVECTOR vExtents = XMLoadFloat3( &Extents );

    // Adjust ray origin to be relative to center of the box.
    XMVECTOR TOrigin = vCenter - Origin;

    // Compute the dot product againt each axis of the box.
    // Since the axii are (1,0,0), (0,1,0), (0,0,1) no computation is necessary.
    XMVECTOR AxisDotOrigin = TOrigin;
    XMVECTOR AxisDotDirection = Direction;

    // if (fabs(AxisDotDirection) <= Epsilon) the ray is nearly parallel to the slab.
    XMVECTOR IsParallel = XMVectorLessOrEqual( XMVectorAbs( AxisDotDirection ), g_RayEpsilon );

    // Test against all three axii simultaneously.
    XMVECTOR InverseAxisDotDirection = XMVectorReciprocal( AxisDotDirection );
    XMVECTOR t1 = ( AxisDotOrigin - vExtents ) * InverseAxisDotDirection;
    XMVECTOR t2 = ( AxisDotOrigin + vExtents ) * InverseAxisDotDirection;

    // Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
    // use the results from any directions parallel to the slab.
    XMVECTOR t_min = XMVectorSelect( XMVectorMin( t1, t2 ), FltMin, IsParallel );
    XMVECTOR t_max = XMVectorSelect( XMVectorMax( t1, t2 ), FltMax, IsParallel );

    // t_min.x = maximum( t_min.x, t_min.y, t_min.z );
    // t_max.x = minimum( t_max.x, t_max.y, t_max.z );
    t_min = XMVectorMax( t_min, XMVectorSplatY( t_min ) );  // x = max(x,y)
    t_min = XMVectorMax( t_min, XMVectorSplatZ( t_min ) );  // x = max(max(x,y),z)
    t_max = XMVectorMin( t_max, XMVectorSplatY( t_max ) );  // x = min(x,y)
    t_max = XMVectorMin( t_max, XMVectorSplatZ( t_max ) );  // x = min(min(x,y),z)

    // if ( t_min > t_max ) return false;
    XMVECTOR NoIntersection = XMVectorGreater( XMVectorSplatX( t_min ), XMVectorSplatX( t_max ) );

    // if ( t_max < 0.0f ) return false;
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( XMVectorSplatX( t_max ), XMVectorZero() ) );

    // if (IsParallel && (-Extents > AxisDotOrigin || Extents < AxisDotOrigin)) return false;
    XMVECTOR ParallelOverlap = XMVectorInBounds( AxisDotOrigin, vExtents );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorAndCInt( IsParallel, ParallelOverlap ) );

    if( !DirectX::Internal::XMVector3AnyTrue( NoIntersection ) )
    {
        // Store the x-component to *pDist
        XMStoreFloat( &Dist, t_min );
        return true;
    }

    return false;
}


//-----------------------------------------------------------------------------
// Create axis-aligned box that contains two other bounding boxes
//-----------------------------------------------------------------------------
inline void BoundingBox::CreateMerged( BoundingBox& Out, const BoundingBox& b1, const BoundingBox& b2 )
{
    XMVECTOR b1Center = XMLoadFloat3( &b1.Center );
    XMVECTOR b1Extents = XMLoadFloat3( &b1.Extents );

    XMVECTOR b2Center = XMLoadFloat3( &b2.Center );
    XMVECTOR b2Extents = XMLoadFloat3( &b2.Extents );

    XMVECTOR Min = XMVectorSubtract( b1Center, b1Extents );
    Min = XMVectorMin( Min, XMVectorSubtract( b2Center, b2Extents ) );

    XMVECTOR Max = XMVectorAdd( b1Center, b1Extents );
    Max = XMVectorMax( Max, XMVectorAdd( b2Center, b2Extents ) );

    assert( XMVector3LessOrEqual( Min, Max ) );

    XMStoreFloat3( &Out.Center, ( Min + Max ) * 0.5f );
    XMStoreFloat3( &Out.Extents, ( Max - Min ) * 0.5f );
}


//-----------------------------------------------------------------------------
// Create axis-aligned box that contains a bounding sphere
//-----------------------------------------------------------------------------
inline void BoundingBox::CreateFromSphere( BoundingBox& Out, const BoundingSphere& sh )
{
    XMVECTOR spCenter = XMLoadFloat3( &sh.Center );
    XMVECTOR shRadius = XMVectorReplicatePtr( &sh.Radius );

    XMVECTOR Min = XMVectorSubtract( spCenter, shRadius );
    XMVECTOR Max = XMVectorAdd( spCenter, shRadius );

    assert( XMVector3LessOrEqual( Min, Max ) );

    XMStoreFloat3( &Out.Center, ( Min + Max ) * 0.5f );
    XMStoreFloat3( &Out.Extents, ( Max - Min ) * 0.5f );
}


//-----------------------------------------------------------------------------
// Create axis-aligned box from min/max points
//-----------------------------------------------------------------------------
inline void BoundingBox::CreateFromPoints( BoundingBox& Out, FXMVECTOR pt1, FXMVECTOR pt2 )
{
    XMVECTOR Min = XMVectorMin( pt1, pt2 );
    XMVECTOR Max = XMVectorMax( pt1, pt2 );

    // Store center and extents.
    XMStoreFloat3( &Out.Center, ( Min + Max ) * 0.5f );
    XMStoreFloat3( &Out.Extents, ( Max - Min ) * 0.5f );
}


//-----------------------------------------------------------------------------
// Find the minimum axis aligned bounding box containing a set of points.
//-----------------------------------------------------------------------------
inline void BoundingBox::CreateFromPoints( BoundingBox& Out, size_t Count, const XMFLOAT3* pPoints, size_t Stride )
{
    assert( Count > 0 );
    assert( pPoints );

    // Find the minimum and maximum x, y, and z
    XMVECTOR vMin, vMax;

    vMin = vMax = XMLoadFloat3( pPoints );

    for( size_t i = 1; i < Count; ++i )
    {
        XMVECTOR Point = XMLoadFloat3( ( const XMFLOAT3* )( ( const uint8_t* )pPoints + i * Stride ) );

        vMin = XMVectorMin( vMin, Point );
        vMax = XMVectorMax( vMax, Point );
    }

    // Store center and extents.
    XMStoreFloat3( &Out.Center, ( vMin + vMax ) * 0.5f );
    XMStoreFloat3( &Out.Extents, ( vMax - vMin ) * 0.5f );
}

