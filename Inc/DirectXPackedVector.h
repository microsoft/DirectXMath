//-------------------------------------------------------------------------------------
// DirectXPackedVector.h -- SIMD C++ Math library
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
    
namespace PackedVector
{

#ifdef _XM_BIGENDIAN_
#pragma bitfield_order(push)
#pragma bitfield_order(lsb_to_msb)
#endif

#pragma warning(push)
#pragma warning(disable:4201 4365 4324)

//------------------------------------------------------------------------------
// ARGB Color; 8-8-8-8 bit unsigned normalized integer components packed into
// a 32 bit integer.  The normalized color is packed into 32 bits using 8 bit
// unsigned, normalized integers for the alpha, red, green, and blue components.
// The alpha component is stored in the most significant bits and the blue
// component in the least significant bits (A8R8G8B8):
// [32] aaaaaaaa rrrrrrrr gggggggg bbbbbbbb [0]
struct XMCOLOR
{
    union
    {
        struct
        {
            uint8_t b;  // Blue:    0/255 to 255/255
            uint8_t g;  // Green:   0/255 to 255/255
            uint8_t r;  // Red:     0/255 to 255/255
            uint8_t a;  // Alpha:   0/255 to 255/255
        };
        uint32_t c;
    };

    XMCOLOR() {}
    XMCOLOR(uint32_t Color) : c(Color) {}
    XMCOLOR(float _r, float _g, float _b, float _a);
    explicit XMCOLOR(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return c; }

    XMCOLOR& operator= (const XMCOLOR& Color) { c = Color.c; return *this; }
    XMCOLOR& operator= (const uint32_t Color) { c = Color; return *this; }
};

//------------------------------------------------------------------------------
// 16 bit floating point number consisting of a sign bit, a 5 bit biased 
// exponent, and a 10 bit mantissa
typedef uint16_t HALF;

//------------------------------------------------------------------------------
// 2D Vector; 16 bit floating point components
struct XMHALF2
{
    union
    {
        struct
        {
            HALF x;
            HALF y;
        };
        uint32_t v;
    };

    XMHALF2() {}
    explicit XMHALF2(uint32_t Packed) : v(Packed) {}
    XMHALF2(HALF _x, HALF _y) : x(_x), y(_y) {}
    explicit XMHALF2(_In_reads_(2) const HALF *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMHALF2(float _x, float _y);
    explicit XMHALF2(_In_reads_(2) const float *pArray);

    XMHALF2& operator= (const XMHALF2& Half2) { x = Half2.x; y = Half2.y; return *this; }
    XMHALF2& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 2D Vector; 16 bit signed normalized integer components
struct XMSHORTN2
{
    union
    {
        struct
        {
            int16_t x;
            int16_t y;
        };
        uint32_t v;
    };

    XMSHORTN2() {}
    explicit XMSHORTN2(uint32_t Packed) : v(Packed) {}
    XMSHORTN2(int16_t _x, int16_t _y) : x(_x), y(_y) {}
    explicit XMSHORTN2(_In_reads_(2) const int16_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMSHORTN2(float _x, float _y);
    explicit XMSHORTN2(_In_reads_(2) const float *pArray);

    XMSHORTN2& operator= (const XMSHORTN2& ShortN2) { x = ShortN2.x; y = ShortN2.y; return *this; }
    XMSHORTN2& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 2D Vector; 16 bit signed integer components
struct XMSHORT2
{
    union
    {
        struct
        {
            int16_t x;
            int16_t y;
        };
        uint32_t v;
    };

    XMSHORT2() {}
    explicit XMSHORT2(uint32_t Packed) : v(Packed) {}
    XMSHORT2(int16_t _x, int16_t _y) : x(_x), y(_y) {}
    explicit XMSHORT2(_In_reads_(2) const int16_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMSHORT2(float _x, float _y);
    explicit XMSHORT2(_In_reads_(2) const float *pArray);

    XMSHORT2& operator= (const XMSHORT2& Short2) { x = Short2.x; y = Short2.y; return *this; }
    XMSHORT2& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 2D Vector; 16 bit unsigned normalized integer components
struct XMUSHORTN2
{
    union
    {
        struct
        {
            uint16_t x;
            uint16_t y;
        };
        uint32_t v;
    };

    XMUSHORTN2() {}
    explicit XMUSHORTN2(uint32_t Packed) : v(Packed) {}
    XMUSHORTN2(uint16_t _x, uint16_t _y) : x(_x), y(_y) {}
    explicit XMUSHORTN2(_In_reads_(2) const uint16_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMUSHORTN2(float _x, float _y);
    explicit XMUSHORTN2(_In_reads_(2) const float *pArray);

    XMUSHORTN2& operator= (const XMUSHORTN2& UShortN2) { x = UShortN2.x; y = UShortN2.y; return *this; }
    XMUSHORTN2& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 2D Vector; 16 bit unsigned integer components
struct XMUSHORT2
{
    union
    {
        struct
        {
            uint16_t x;
            uint16_t y;
        };
        uint32_t v;
    };

    XMUSHORT2() {}
    explicit XMUSHORT2(uint32_t Packed) : v(Packed) {}
    XMUSHORT2(uint16_t _x, uint16_t _y) : x(_x), y(_y) {}
    explicit XMUSHORT2(_In_reads_(2) const uint16_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMUSHORT2(float _x, float _y);
    explicit XMUSHORT2(_In_reads_(2) const float *pArray);

    XMUSHORT2& operator= (const XMUSHORT2& UShort2) { x = UShort2.x; y = UShort2.y; return *this; }
    XMUSHORT2& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 2D Vector; 8 bit signed normalized integer components
struct XMBYTEN2
{
    union
    {
        struct
        {
            int8_t x;
            int8_t y;
        };
        uint16_t v;
    };

    XMBYTEN2() {}
    explicit XMBYTEN2(uint16_t Packed) : v(Packed) {}
    XMBYTEN2(int8_t _x, int8_t _y) : x(_x), y(_y) {}
    explicit XMBYTEN2(_In_reads_(2) const int8_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMBYTEN2(float _x, float _y);
    explicit XMBYTEN2(_In_reads_(2) const float *pArray);

    XMBYTEN2& operator= (const XMBYTEN2& ByteN2) { x = ByteN2.x; y = ByteN2.y; return *this; }
    XMBYTEN2& operator= (uint16_t Packed) { v = Packed; return *this; }
};

// 2D Vector; 8 bit signed integer components
struct XMBYTE2
{
    union
    {
        struct
        {
            int8_t x;
            int8_t y;
        };
        uint16_t v;
    };

    XMBYTE2() {}
    explicit XMBYTE2(uint16_t Packed) : v(Packed) {}
    XMBYTE2(int8_t _x, int8_t _y) : x(_x), y(_y) {}
    explicit XMBYTE2(_In_reads_(2) const int8_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMBYTE2(float _x, float _y);
    explicit XMBYTE2(_In_reads_(2) const float *pArray);

    XMBYTE2& operator= (const XMBYTE2& Byte2) { x = Byte2.x; y = Byte2.y; return *this; }
    XMBYTE2& operator= (uint16_t Packed) { v = Packed; return *this; }
};

// 2D Vector; 8 bit unsigned normalized integer components
struct XMUBYTEN2
{
    union
    {
        struct
        {
            uint8_t x;
            uint8_t y;
        };
        uint16_t v;
    };

    XMUBYTEN2() {}
    explicit XMUBYTEN2(uint16_t Packed) : v(Packed) {}
    XMUBYTEN2(uint8_t _x, uint8_t _y) : x(_x), y(_y) {}
    explicit XMUBYTEN2(_In_reads_(2) const uint8_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMUBYTEN2(float _x, float _y);
    explicit XMUBYTEN2(_In_reads_(2) const float *pArray);

    XMUBYTEN2& operator= (const XMUBYTEN2& UByteN2) { x = UByteN2.x; y = UByteN2.y; return *this; }
    XMUBYTEN2& operator= (uint16_t Packed) { v = Packed; return *this; }
};

// 2D Vector; 8 bit unsigned integer components
struct XMUBYTE2
{
    union
    {
        struct
        {
            uint8_t x;
            uint8_t y;
        };
        uint16_t v;
    };

    XMUBYTE2() {}
    explicit XMUBYTE2(uint16_t Packed) : v(Packed) {}
    XMUBYTE2(uint8_t _x, uint8_t _y) : x(_x), y(_y) {}
    explicit XMUBYTE2(_In_reads_(2) const uint8_t *pArray) : x(pArray[0]), y(pArray[1]) {}
    XMUBYTE2(float _x, float _y);
    explicit XMUBYTE2(_In_reads_(2) const float *pArray);

    XMUBYTE2& operator= (const XMUBYTE2& UByte2) { x = UByte2.x; y = UByte2.y; return *this; }
    XMUBYTE2& operator= (uint16_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 3D vector: 5/6/5 unsigned integer components
struct XMU565
{
    union
    {
        struct
        {
            uint16_t x  : 5;    // 0 to 31
            uint16_t y  : 6;    // 0 to 63
            uint16_t z  : 5;    // 0 to 31
        };
        uint16_t v;
    };

    XMU565() {}
    explicit XMU565(uint16_t Packed) : v(Packed) {}
    XMU565(uint8_t _x, uint8_t _y, uint8_t _z) : x(_x), y(_y), z(_z) {}
    explicit XMU565(_In_reads_(3) const int8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}
    XMU565(float _x, float _y, float _z);
    explicit XMU565(_In_reads_(3) const float *pArray);

    operator uint16_t () const { return v; }

    XMU565& operator= (const XMU565& U565) { v = U565.v; return *this; }
    XMU565& operator= (uint16_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 3D vector: 11/11/10 floating-point components
// The 3D vector is packed into 32 bits as follows: a 5-bit biased exponent
// and 6-bit mantissa for x component, a 5-bit biased exponent and
// 6-bit mantissa for y component, a 5-bit biased exponent and a 5-bit
// mantissa for z. The z component is stored in the most significant bits
// and the x component in the least significant bits. No sign bits so
// all partial-precision numbers are positive.
// (Z10Y11X11): [32] ZZZZZzzz zzzYYYYY yyyyyyXX XXXxxxxx [0]
struct XMFLOAT3PK
{
    union
    {
        struct
        {
            uint32_t xm : 6; // x-mantissa
            uint32_t xe : 5; // x-exponent
            uint32_t ym : 6; // y-mantissa
            uint32_t ye : 5; // y-exponent
            uint32_t zm : 5; // z-mantissa
            uint32_t ze : 5; // z-exponent
        };
        uint32_t v;
    };

    XMFLOAT3PK() {}
    explicit XMFLOAT3PK(uint32_t Packed) : v(Packed) {}
    XMFLOAT3PK(float _x, float _y, float _z);
    explicit XMFLOAT3PK(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMFLOAT3PK& operator= (const XMFLOAT3PK& float3pk) { v = float3pk.v; return *this; }
    XMFLOAT3PK& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 3D vector: 9/9/9 floating-point components with shared 5-bit exponent
// The 3D vector is packed into 32 bits as follows: a 5-bit biased exponent
// with 9-bit mantissa for the x, y, and z component. The shared exponent
// is stored in the most significant bits and the x component mantissa is in
// the least significant bits. No sign bits so all partial-precision numbers
// are positive.
// (E5Z9Y9X9): [32] EEEEEzzz zzzzzzyy yyyyyyyx xxxxxxxx [0]
struct XMFLOAT3SE
{
    union
    {
        struct
        {
            uint32_t xm : 9; // x-mantissa
            uint32_t ym : 9; // y-mantissa
            uint32_t zm : 9; // z-mantissa
            uint32_t e  : 5; // shared exponent
        };
        uint32_t v;
    };

    XMFLOAT3SE() {}
    explicit XMFLOAT3SE(uint32_t Packed) : v(Packed) {}
    XMFLOAT3SE(float _x, float _y, float _z);
    explicit XMFLOAT3SE(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMFLOAT3SE& operator= (const XMFLOAT3SE& float3se) { v = float3se.v; return *this; }
    XMFLOAT3SE& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D Vector; 16 bit floating point components
struct XMHALF4
{
    union
    {
        struct
        {
            HALF x;
            HALF y;
            HALF z;
            HALF w;
        };
        uint64_t v;
    };

    XMHALF4() {}
    explicit XMHALF4(uint64_t Packed) : v(Packed) {}
    XMHALF4(HALF _x, HALF _y, HALF _z, HALF _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMHALF4(_In_reads_(4) const HALF *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMHALF4(float _x, float _y, float _z, float _w);
    explicit XMHALF4(_In_reads_(4) const float *pArray);

    XMHALF4& operator= (const XMHALF4& Half4) { x = Half4.x; y = Half4.y; z = Half4.z; w = Half4.w; return *this; }
    XMHALF4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D Vector; 16 bit signed normalized integer components
struct XMSHORTN4
{
    union
    {
        struct
        {
            int16_t x;
            int16_t y;
            int16_t z;
            int16_t w;
        };
        uint64_t v;
    };

    XMSHORTN4() {}
    explicit XMSHORTN4(uint64_t Packed) : v(Packed) {}
    XMSHORTN4(int16_t _x, int16_t _y, int16_t _z, int16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMSHORTN4(_In_reads_(4) const int16_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMSHORTN4(float _x, float _y, float _z, float _w);
    explicit XMSHORTN4(_In_reads_(4) const float *pArray);

    XMSHORTN4& operator= (const XMSHORTN4& ShortN4) { x = ShortN4.x; y = ShortN4.y; z = ShortN4.z; w = ShortN4.w; return *this; }
    XMSHORTN4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 16 bit signed integer components
struct XMSHORT4
{
    union
    {
        struct
        {
            int16_t x;
            int16_t y;
            int16_t z;
            int16_t w;
        };
        uint64_t v;
    };

    XMSHORT4() {}
    explicit XMSHORT4(uint64_t Packed) : v(Packed) {}
    XMSHORT4(int16_t _x, int16_t _y, int16_t _z, int16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMSHORT4(_In_reads_(4) const int16_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMSHORT4(float _x, float _y, float _z, float _w);
    explicit XMSHORT4(_In_reads_(4) const float *pArray);

    XMSHORT4& operator= (const XMSHORT4& Short4) { x = Short4.x; y = Short4.y; z = Short4.z; w = Short4.w; return *this; }
    XMSHORT4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 16 bit unsigned normalized integer components
struct XMUSHORTN4
{
    union
    {
        struct
        {
            uint16_t x;
            uint16_t y;
            uint16_t z;
            uint16_t w;
        };
        uint64_t v;
    };

    XMUSHORTN4() {}
    explicit XMUSHORTN4(uint64_t Packed) : v(Packed) {}
    XMUSHORTN4(uint16_t _x, uint16_t _y, uint16_t _z, uint16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMUSHORTN4(_In_reads_(4) const uint16_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMUSHORTN4(float _x, float _y, float _z, float _w);
    explicit XMUSHORTN4(_In_reads_(4) const float *pArray);

    XMUSHORTN4& operator= (const XMUSHORTN4& UShortN4) { x = UShortN4.x; y = UShortN4.y; z = UShortN4.z; w = UShortN4.w; return *this; }
    XMUSHORTN4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 16 bit unsigned integer components
struct XMUSHORT4
{
    union
    {
        struct
        {
            uint16_t x;
            uint16_t y;
            uint16_t z;
            uint16_t w;
        };
        uint64_t v;
    };

    XMUSHORT4() {}
    explicit XMUSHORT4(uint64_t Packed) : v(Packed) {}
    XMUSHORT4(uint16_t _x, uint16_t _y, uint16_t _z, uint16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMUSHORT4(_In_reads_(4) const uint16_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMUSHORT4(float _x, float _y, float _z, float _w);
    explicit XMUSHORT4(_In_reads_(4) const float *pArray);

    XMUSHORT4& operator= (const XMUSHORT4& UShort4) { x = UShort4.x; y = UShort4.y; z = UShort4.z; w = UShort4.w; return *this; }
    XMUSHORT4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D Vector; 10-10-10-2 bit normalized components packed into a 32 bit integer
// The normalized 4D Vector is packed into 32 bits as follows: a 2 bit unsigned, 
// normalized integer for the w component and 10 bit signed, normalized 
// integers for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
struct XMXDECN4
{
    union
    {
        struct
        {
            int32_t x   : 10;    // -511/511 to 511/511
            int32_t y   : 10;    // -511/511 to 511/511
            int32_t z   : 10;    // -511/511 to 511/511
            uint32_t w  : 2;     //      0/3 to     3/3
        };
        uint32_t v;
    };

    XMXDECN4() {}
    explicit XMXDECN4(uint32_t Packed) : v(Packed) {}
    XMXDECN4(float _x, float _y, float _z, float _w);
    explicit XMXDECN4(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return v; }

    XMXDECN4& operator= (const XMXDECN4& XDecN4) { v = XDecN4.v; return *this; }
    XMXDECN4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 10-10-10-2 bit components packed into a 32 bit integer
// The normalized 4D Vector is packed into 32 bits as follows: a 2 bit unsigned
// integer for the w component and 10 bit signed integers for the 
// z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
struct XMXDEC4
{
    union
    {
        struct
        {
            int32_t x   : 10;    // -511 to 511
            int32_t y   : 10;    // -511 to 511
            int32_t z   : 10;    // -511 to 511
            uint32_t w  : 2;     // 0 to 3
        };
        uint32_t v;
    };

    XMXDEC4() {}
    explicit XMXDEC4(uint32_t Packed) : v(Packed) {}
    XMXDEC4(float _x, float _y, float _z, float _w);
    explicit XMXDEC4(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return v; }

    XMXDEC4& operator= (const XMXDEC4& XDec4) { v = XDec4.v; return *this; }
    XMXDEC4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 10-10-10-2 bit normalized components packed into a 32 bit integer
// The normalized 4D Vector is packed into 32 bits as follows: a 2 bit signed, 
// normalized integer for the w component and 10 bit signed, normalized 
// integers for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
struct XMDECN4
{
    union
    {
        struct
        {
            int32_t x   : 10;    // -511/511 to 511/511
            int32_t y   : 10;    // -511/511 to 511/511
            int32_t z   : 10;    // -511/511 to 511/511
            int32_t w   : 2;     //     -1/1 to     1/1
        };
        uint32_t v;
    };

    XMDECN4() {}
    explicit XMDECN4(uint32_t Packed) : v(Packed) {}
    XMDECN4(float _x, float _y, float _z, float _w);
    explicit XMDECN4(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return v; }

    XMDECN4& operator= (const XMDECN4& DecN4) { v = DecN4.v; return *this; }
    XMDECN4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 10-10-10-2 bit components packed into a 32 bit integer
// The 4D Vector is packed into 32 bits as follows: a 2 bit signed, 
// integer for the w component and 10 bit signed integers for the 
// z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
struct XMDEC4
{
    union
    {
        struct
        {
            int32_t  x  : 10;    // -511 to 511
            int32_t  y  : 10;    // -511 to 511
            int32_t  z  : 10;    // -511 to 511
            int32_t  w  : 2;     //   -1 to   1
        };
        uint32_t v;
    };

    XMDEC4() {}
    explicit XMDEC4(uint32_t Packed) : v(Packed) {}
    XMDEC4(float _x, float _y, float _z, float _w);
    explicit XMDEC4(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return v; }

    XMDEC4& operator= (const XMDEC4& Dec4) { v = Dec4.v; return *this; }
    XMDEC4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 10-10-10-2 bit normalized components packed into a 32 bit integer
// The normalized 4D Vector is packed into 32 bits as follows: a 2 bit unsigned, 
// normalized integer for the w component and 10 bit unsigned, normalized 
// integers for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
struct XMUDECN4
{
    union
    {
        struct
        {
            uint32_t x  : 10;    // 0/1023 to 1023/1023
            uint32_t y  : 10;    // 0/1023 to 1023/1023
            uint32_t z  : 10;    // 0/1023 to 1023/1023
            uint32_t w  : 2;     //    0/3 to       3/3
        };
        uint32_t v;
    };

    XMUDECN4() {}
    explicit XMUDECN4(uint32_t Packed) : v(Packed) {}
    XMUDECN4(float _x, float _y, float _z, float _w);
    explicit XMUDECN4(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return v; }

    XMUDECN4& operator= (const XMUDECN4& UDecN4) { v = UDecN4.v; return *this; }
    XMUDECN4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 10-10-10-2 bit components packed into a 32 bit integer
// The 4D Vector is packed into 32 bits as follows: a 2 bit unsigned, 
// integer for the w component and 10 bit unsigned integers 
// for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
struct XMUDEC4
{
    union
    {
        struct
        {
            uint32_t x  : 10;    // 0 to 1023
            uint32_t y  : 10;    // 0 to 1023
            uint32_t z  : 10;    // 0 to 1023
            uint32_t w  : 2;     // 0 to    3
        };
        uint32_t v;
    };

    XMUDEC4() {}
    explicit XMUDEC4(uint32_t Packed) : v(Packed) {}
    XMUDEC4(float _x, float _y, float _z, float _w);
    explicit XMUDEC4(_In_reads_(4) const float *pArray);

    operator uint32_t () const { return v; }

    XMUDEC4& operator= (const XMUDEC4& UDec4) { v = UDec4.v; return *this; }
    XMUDEC4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D Vector; 8 bit signed normalized integer components
struct XMBYTEN4
{
    union
    {
        struct
        {
            int8_t x;
            int8_t y;
            int8_t z;
            int8_t w;
        };
        uint32_t v;
    };

    XMBYTEN4() {}
    XMBYTEN4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMBYTEN4(uint32_t Packed) : v(Packed) {}
    explicit XMBYTEN4(_In_reads_(4) const int8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMBYTEN4(float _x, float _y, float _z, float _w);
    explicit XMBYTEN4(_In_reads_(4) const float *pArray);

    XMBYTEN4& operator= (const XMBYTEN4& ByteN4) { x = ByteN4.x; y = ByteN4.y; z = ByteN4.z; w = ByteN4.w; return *this; }
    XMBYTEN4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 8 bit signed integer components
struct XMBYTE4
{
    union
    {
        struct
        {
            int8_t x;
            int8_t y;
            int8_t z;
            int8_t w;
        };
        uint32_t v;
    };

    XMBYTE4() {}
    XMBYTE4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMBYTE4(uint32_t Packed) : v(Packed) {}
    explicit XMBYTE4(_In_reads_(4) const int8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMBYTE4(float _x, float _y, float _z, float _w);
    explicit XMBYTE4(_In_reads_(4) const float *pArray);

    XMBYTE4& operator= (const XMBYTE4& Byte4) { x = Byte4.x; y = Byte4.y; z = Byte4.z; w = Byte4.w; return *this; }
    XMBYTE4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 8 bit unsigned normalized integer components
struct XMUBYTEN4
{
    union
    {
        struct
        {
            uint8_t x;
            uint8_t y;
            uint8_t z;
            uint8_t w;
        };
        uint32_t v;
    };

    XMUBYTEN4() {}
    XMUBYTEN4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMUBYTEN4(uint32_t Packed) : v(Packed) {}
    explicit XMUBYTEN4(_In_reads_(4) const uint8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMUBYTEN4(float _x, float _y, float _z, float _w);
    explicit XMUBYTEN4(_In_reads_(4) const float *pArray);

    XMUBYTEN4& operator= (const XMUBYTEN4& UByteN4) { x = UByteN4.x; y = UByteN4.y; z = UByteN4.z; w = UByteN4.w; return *this; }
    XMUBYTEN4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 8 bit unsigned integer components
struct XMUBYTE4
{
    union
    {
        struct
        {
            uint8_t x;
            uint8_t y;
            uint8_t z;
            uint8_t w;
        };
        uint32_t v;
    };

    XMUBYTE4() {}
    XMUBYTE4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMUBYTE4(uint32_t Packed) : v(Packed) {}
    explicit XMUBYTE4(_In_reads_(4) const uint8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMUBYTE4(float _x, float _y, float _z, float _w);
    explicit XMUBYTE4(_In_reads_(4) const float *pArray);

    XMUBYTE4& operator= (const XMUBYTE4& UByte4) { x = UByte4.x; y = UByte4.y; z = UByte4.z; w = UByte4.w; return *this; }
    XMUBYTE4& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D vector; 4 bit unsigned integer components
struct XMUNIBBLE4
{
    union
    {
        struct
        {
            uint16_t x  : 4;    // 0 to 15
            uint16_t y  : 4;    // 0 to 15
            uint16_t z  : 4;    // 0 to 15
            uint16_t w  : 4;    // 0 to 15
        };
        uint16_t v;
    };

    XMUNIBBLE4() {}
    explicit XMUNIBBLE4(uint16_t Packed) : v(Packed) {}
    XMUNIBBLE4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMUNIBBLE4(_In_reads_(4) const int8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    XMUNIBBLE4(float _x, float _y, float _z, float _w);
    explicit XMUNIBBLE4(_In_reads_(4) const float *pArray);

    operator uint16_t () const { return v; }

    XMUNIBBLE4& operator= (const XMUNIBBLE4& UNibble4) { v = UNibble4.v; return *this; }
    XMUNIBBLE4& operator= (uint16_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D vector: 5/5/5/1 unsigned integer components
struct XMU555
{
    union
    {
        struct
        {
            uint16_t x  : 5;    // 0 to 31
            uint16_t y  : 5;    // 0 to 31
            uint16_t z  : 5;    // 0 to 31
            uint16_t w  : 1;    // 0 or 1
        };
        uint16_t v;
    };

    XMU555() {}
    explicit XMU555(uint16_t Packed) : v(Packed) {}
    XMU555(int8_t _x, int8_t _y, int8_t _z, bool _w) : x(_x), y(_y), z(_z), w(_w ? 0x1 : 0) {}
    XMU555(_In_reads_(3) const int8_t *pArray, _In_ bool _w) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(_w ? 0x1 : 0) {}
    XMU555(float _x, float _y, float _z, bool _w);
    XMU555(_In_reads_(3) const float *pArray, _In_ bool _w);

    operator uint16_t () const { return v; }

    XMU555& operator= (const XMU555& U555) { v = U555.v; return *this; }
    XMU555& operator= (uint16_t Packed) { v = Packed; return *this; }
};

///begin_xbox360
////////////////////////////////////////////////////////////////////////////////
// PackedVector::Xbox
////////////////////////////////////////////////////////////////////////////////

namespace Xbox
{

//------------------------------------------------------------------------------
// 3D Vector; 11-11-10 bit normalized components packed into a 32 bit integer
// The normalized 3D Vector is packed into 32 bits as follows: a 10 bit signed, 
// normalized integer for the z component and 11 bit signed, normalized 
// integers for the x and y components.  The z component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z10Y11X11): [32] zzzzzzzz zzyyyyyy yyyyyxxx xxxxxxxx [0]
struct XMHENDN3
{
    union
    {
        struct
        {
            int32_t x   : 11;    // -1023/1023 to 1023/1023
            int32_t y   : 11;    // -1023/1023 to 1023/1023
            int32_t z   : 10;    // -511/511 to 511/511
        };
        uint32_t v;
    };

    XMHENDN3() {}
    explicit XMHENDN3(uint32_t Packed) : v(Packed) {}
    XMHENDN3(float _x, float _y, float _z);
    explicit XMHENDN3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMHENDN3& operator= (const XMHENDN3& HenDN3) { v = HenDN3.v; return *this; }
    XMHENDN3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 11-11-10 bit components packed into a 32 bit integer
// The 3D Vector is packed into 32 bits as follows: a 10 bit signed, 
// integer for the z component and 11 bit signed integers for the 
// x and y components.  The z component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z10Y11X11): [32] zzzzzzzz zzyyyyyy yyyyyxxx xxxxxxxx [0]
struct XMHEND3
{
    union
    {
        struct
        {
            int32_t x   : 11;    // -1023 to 1023
            int32_t y   : 11;    // -1023 to 1023
            int32_t z   : 10;    // -511 to 511
        };
        uint32_t v;
    };

    XMHEND3() {}
    explicit XMHEND3(uint32_t Packed) : v(Packed) {}
    XMHEND3(float _x, float _y, float _z);
    explicit XMHEND3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMHEND3& operator= (const XMHEND3& HenD3) { v = HenD3.v; return *this; }
    XMHEND3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 11-11-10 bit normalized components packed into a 32 bit integer
// The normalized 3D Vector is packed into 32 bits as follows: a 10 bit unsigned, 
// normalized integer for the z component and 11 bit unsigned, normalized 
// integers for the x and y components.  The z component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z10Y11X11): [32] zzzzzzzz zzyyyyyy yyyyyxxx xxxxxxxx [0]
struct XMUHENDN3
{
    union
    {
        struct
        {
            uint32_t x  : 11;    // 0/2047 to 2047/2047
            uint32_t y  : 11;    // 0/2047 to 2047/2047
            uint32_t z  : 10;    // 0/1023 to 1023/1023
        };
        uint32_t v;
    };

    XMUHENDN3() {}
    explicit XMUHENDN3(uint32_t Packed) : v(Packed) {}
    XMUHENDN3(float _x, float _y, float _z);
    explicit XMUHENDN3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMUHENDN3& operator= (const XMUHENDN3& UHenDN3) { v = UHenDN3.v; return *this; }
    XMUHENDN3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 11-11-10 bit components packed into a 32 bit integer
// The 3D Vector is packed into 32 bits as follows: a 10 bit unsigned
// integer for the z component and 11 bit unsigned integers 
// for the x and y components.  The z component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z10Y11X11): [32] zzzzzzzz zzyyyyyy yyyyyxxx xxxxxxxx [0]
struct XMUHEND3
{
    union
    {
        struct
        {
            uint32_t x  : 11;    // 0 to 2047
            uint32_t y  : 11;    // 0 to 2047
            uint32_t z  : 10;    // 0 to 1023
        };
        uint32_t v;
    };

    XMUHEND3() {}
    explicit XMUHEND3(uint32_t Packed) : v(Packed) {}
    XMUHEND3(float _x, float _y, float _z);
    explicit XMUHEND3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMUHEND3& operator= (const XMUHEND3& UHenD3) { v = UHenD3.v; return *this; }
    XMUHEND3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 10-11-11 bit normalized components packed into a 32 bit integer
// The normalized 3D Vector is packed into 32 bits as follows: a 10 bit signed, 
// normalized integer for the x component and 11 bit signed, normalized 
// integers for the y and z components.  The z component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z11Y11X10): [32] zzzzzzzz zzzyyyyy yyyyyyxx xxxxxxxx [0]
struct XMDHENN3
{
    union
    {
        struct
        {
            int32_t x   : 10;    // -511/511 to 511/511
            int32_t y   : 11;    // -1023/1023 to 1023/1023
            int32_t z   : 11;    // -1023/1023 to 1023/1023
        };
        uint32_t v;
    };

    XMDHENN3() {};
    explicit XMDHENN3(uint32_t Packed) : v(Packed) {};
    XMDHENN3(float _x, float _y, float _z);
    explicit XMDHENN3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMDHENN3& operator= (const XMDHENN3& DHenN3) { v = DHenN3.v; return *this; }
    XMDHENN3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 10-11-11 bit components packed into a 32 bit integer
// The 3D Vector is packed into 32 bits as follows: a 10 bit signed, 
// integer for the x component and 11 bit signed integers for the 
// y and z components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z11Y11X10): [32] zzzzzzzz zzzyyyyy yyyyyyxx xxxxxxxx [0]
struct XMDHEN3
{
    union
    {
        struct
        {
            int32_t x   : 10;    // -511 to 511
            int32_t y   : 11;    // -1023 to 1023
            int32_t z   : 11;    // -1023 to 1023
        };
        uint32_t v;
    };

    XMDHEN3() {}
    explicit XMDHEN3(uint32_t Packed) : v(Packed) {}
    XMDHEN3(float _x, float _y, float _z);
    explicit XMDHEN3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMDHEN3& operator= (const XMDHEN3& DHen3) { v = DHen3.v; return *this; }
    XMDHEN3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 10-11-11 bit normalized components packed into a 32 bit integer
// The normalized 3D Vector is packed into 32 bits as follows: a 10 bit unsigned, 
// normalized integer for the x component and 11 bit unsigned, normalized 
// integers for the y and z components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z11Y11X10): [32] zzzzzzzz zzzyyyyy yyyyyyxx xxxxxxxx [0]
struct XMUDHENN3
{
    union
    {
        struct
        {
            uint32_t x  : 10;    // 0/1023 to 1023/1023
            uint32_t y  : 11;    // 0/2047 to 2047/2047
            uint32_t z  : 11;    // 0/2047 to 2047/2047
        };
        uint32_t v;
    };

    XMUDHENN3() {}
    explicit XMUDHENN3(uint32_t Packed) : v(Packed) {}
    XMUDHENN3(float _x, float _y, float _z);
    explicit XMUDHENN3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMUDHENN3& operator= (const XMUDHENN3& UDHenN3) { v = UDHenN3.v; return *this; }
    XMUDHENN3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

// 3D Vector; 10-11-11 bit components packed into a 32 bit integer
// The 3D Vector is packed into 32 bits as follows: a 10 bit unsigned, 
// integer for the x component and 11 bit unsigned integers 
// for the y and z components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (Z11Y11X10): [32] zzzzzzzz zzzyyyyy yyyyyyxx xxxxxxxx [0]
struct XMUDHEN3
{
    union
    {
        struct
        {
            uint32_t x  : 10;    // 0 to 1023
            uint32_t y  : 11;    // 0 to 2047
            uint32_t z  : 11;    // 0 to 2047
        };
        uint32_t v;
    };

    XMUDHEN3() {}
    explicit XMUDHEN3(uint32_t Packed) : v(Packed) {}
    XMUDHEN3(float _x, float _y, float _z);
    explicit XMUDHEN3(_In_reads_(3) const float *pArray);

    operator uint32_t () const { return v; }

    XMUDHEN3& operator= (const XMUDHEN3& UDHen3) { v = UDHen3.v; return *this; }
    XMUDHEN3& operator= (uint32_t Packed) { v = Packed; return *this; }
};

//------------------------------------------------------------------------------
// 4D Vector; 20-20-20-4 bit normalized components packed into a 64 bit integer
// The normalized 4D Vector is packed into 64 bits as follows: a 4 bit unsigned, 
// normalized integer for the w component and 20 bit signed, normalized 
// integers for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W4Z20Y20X20): [64] wwwwzzzz zzzzzzzz zzzzzzzz yyyyyyyy yyyyyyyy yyyyxxxx xxxxxxxx xxxxxxxx [0]
struct XMXICON4
{
    union
    {
        struct
        {
            int64_t x           : 20;    // -524287/524287 to 524287/524287
            int64_t y           : 20;    // -524287/524287 to 524287/524287
            int64_t z           : 20;    // -524287/524287 to 524287/524287
            uint64_t w          : 4;     //           0/15 to         15/15
        };
        uint64_t v;
    };

    XMXICON4() {}
    explicit XMXICON4(uint64_t Packed) : v(Packed) {}
    XMXICON4(float _x, float _y, float _z, float _w);
    explicit XMXICON4(_In_reads_(4) const float *pArray);

    operator uint64_t () const { return v; }

    XMXICON4& operator= (const XMXICON4& XIcoN4) { v = XIcoN4.v; return *this; }
    XMXICON4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 20-20-20-4 bit components packed into a 64 bit integer
// The 4D Vector is packed into 64 bits as follows: a 4 bit unsigned
// integer for the w component and 20 bit signed integers for the 
// z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W4Z20Y20X20): [64] wwwwzzzz zzzzzzzz zzzzzzzz yyyyyyyy yyyyyyyy yyyyxxxx xxxxxxxx xxxxxxxx [0]
struct XMXICO4
{
    union
    {
        struct
        {
            int64_t x           : 20;    // -524287 to 524287
            int64_t y           : 20;    // -524287 to 524287
            int64_t z           : 20;    // -524287 to 524287
            uint64_t w          : 4;     //       0 to     15
        };
        uint64_t v;
    };

    XMXICO4() {}
    explicit XMXICO4(uint64_t Packed) : v(Packed) {}
    XMXICO4(float _x, float _y, float _z, float _w);
    explicit XMXICO4(_In_reads_(4) const float *pArray);

    operator uint64_t () const { return v; }

    XMXICO4& operator= (const XMXICO4& XIco4) { v = XIco4.v; return *this; }
    XMXICO4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 20-20-20-4 bit normalized components packed into a 64 bit integer
// The normalized 4D Vector is packed into 64 bits as follows: a 4 bit signed, 
// normalized integer for the w component and 20 bit signed, normalized 
// integers for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W4Z20Y20X20): [64] wwwwzzzz zzzzzzzz zzzzzzzz yyyyyyyy yyyyyyyy yyyyxxxx xxxxxxxx xxxxxxxx [0]
struct XMICON4
{
    union
    {
        struct
        {
            int64_t x   : 20;    // -524287/524287 to 524287/524287
            int64_t y   : 20;    // -524287/524287 to 524287/524287
            int64_t z   : 20;    // -524287/524287 to 524287/524287
            int64_t w   : 4;     //           -7/7 to           7/7
        };
        uint64_t v;
    };

    XMICON4() {}
    explicit XMICON4(uint64_t Packed) : v(Packed) {}
    XMICON4(float _x, float _y, float _z, float _w);
    explicit XMICON4(_In_reads_(4) const float *pArray);

    operator uint64_t () const { return v; }

    XMICON4& operator= (const XMICON4& IcoN4) { v = IcoN4.v; return *this; }
    XMICON4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 20-20-20-4 bit components packed into a 64 bit integer
// The 4D Vector is packed into 64 bits as follows: a 4 bit signed, 
// integer for the w component and 20 bit signed integers for the 
// z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W4Z20Y20X20): [64] wwwwzzzz zzzzzzzz zzzzzzzz yyyyyyyy yyyyyyyy yyyyxxxx xxxxxxxx xxxxxxxx [0]
struct XMICO4
{
    union
    {
        struct
        {
            int64_t x   : 20;    // -524287 to 524287
            int64_t y   : 20;    // -524287 to 524287
            int64_t z   : 20;    // -524287 to 524287
            int64_t w   : 4;     //      -7 to      7
        };
        uint64_t v;
    };

    XMICO4() {}
    explicit XMICO4(uint64_t Packed) : v(Packed) {}
    XMICO4(float _x, float _y, float _z, float _w);
    explicit XMICO4(_In_reads_(4) const float *pArray);

    operator uint64_t () const { return v; }

    XMICO4& operator= (const XMICO4& Ico4) { v = Ico4.v; return *this; }
    XMICO4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 20-20-20-4 bit normalized components packed into a 64 bit integer
// The normalized 4D Vector is packed into 64 bits as follows: a 4 bit unsigned, 
// normalized integer for the w component and 20 bit unsigned, normalized 
// integers for the z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W4Z20Y20X20): [64] wwwwzzzz zzzzzzzz zzzzzzzz yyyyyyyy yyyyyyyy yyyyxxxx xxxxxxxx xxxxxxxx [0]
struct XMUICON4
{
    union
    {
        struct
        {
            uint64_t x  : 20;    // 0/1048575 to 1048575/1048575
            uint64_t y  : 20;    // 0/1048575 to 1048575/1048575
            uint64_t z  : 20;    // 0/1048575 to 1048575/1048575
            uint64_t w  : 4;     //      0/15 to           15/15
        };
        uint64_t v;
    };

    XMUICON4() {}
    explicit XMUICON4(uint64_t Packed) : v(Packed) {}
    XMUICON4(float _x, float _y, float _z, float _w);
    explicit XMUICON4(_In_reads_(4) const float *pArray);

    operator uint64_t () const { return v; }

    XMUICON4& operator= (const XMUICON4& UIcoN4) { v = UIcoN4.v; return *this; }
    XMUICON4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

// 4D Vector; 20-20-20-4 bit components packed into a 64 bit integer
// The 4D Vector is packed into 64 bits as follows: a 4 bit unsigned 
// integer for the w component and 20 bit unsigned integers for the 
// z, y, and x components.  The w component is stored in the 
// most significant bits and the x component in the least significant bits
// (W4Z20Y20X20): [64] wwwwzzzz zzzzzzzz zzzzzzzz yyyyyyyy yyyyyyyy yyyyxxxx xxxxxxxx xxxxxxxx [0]
struct XMUICO4
{
    union
    {
        struct
        {
            uint64_t x  : 20;    // 0 to 1048575
            uint64_t y  : 20;    // 0 to 1048575
            uint64_t z  : 20;    // 0 to 1048575
            uint64_t w  : 4;     // 0 to      15
        };
        uint64_t v;
    };

    XMUICO4() {}
    explicit XMUICO4(uint64_t Packed) : v(Packed) {}
    XMUICO4(float _x, float _y, float _z, float _w);
    explicit XMUICO4(_In_reads_(4) const float *pArray);

    operator uint64_t () const { return v; }

    XMUICO4& operator= (const XMUICO4& UIco4) { v = UIco4.v; return *this; }
    XMUICO4& operator= (uint64_t Packed) { v = Packed; return *this; }
};

}; // namespace Xbox
///end_xbox360

#pragma warning(pop)

#ifdef _XM_BIGENDIAN_
#pragma bitfield_order(pop)
#endif


/****************************************************************************
 *
 * Data conversion operations
 *
 ****************************************************************************/

float           XMConvertHalfToFloat(HALF Value);
float*          XMConvertHalfToFloatStream(_Out_writes_bytes_(sizeof(float)+OutputStride*(HalfCount-1)) float* pOutputStream,
                                           _In_ size_t OutputStride,
                                           _In_reads_bytes_(sizeof(HALF)+InputStride*(HalfCount-1)) const HALF* pInputStream,
                                           _In_ size_t InputStride, _In_ size_t HalfCount);
HALF            XMConvertFloatToHalf(float Value);
HALF*           XMConvertFloatToHalfStream(_Out_writes_bytes_(sizeof(HALF)+OutputStride*(FloatCount-1)) HALF* pOutputStream,
                                           _In_ size_t OutputStride,
                                           _In_reads_bytes_(sizeof(float)+InputStride*(FloatCount-1)) const float* pInputStream,
                                           _In_ size_t InputStride, _In_ size_t FloatCount);

/****************************************************************************
 *
 * Load operations
 *
 ****************************************************************************/

XMVECTOR        XMLoadColor(_In_ const XMCOLOR* pSource);

XMVECTOR        XMLoadHalf2(_In_ const XMHALF2* pSource);
XMVECTOR        XMLoadShortN2(_In_ const XMSHORTN2* pSource);
XMVECTOR        XMLoadShort2(_In_ const XMSHORT2* pSource);
XMVECTOR        XMLoadUShortN2(_In_ const XMUSHORTN2* pSource);
XMVECTOR        XMLoadUShort2(_In_ const XMUSHORT2* pSource);
XMVECTOR        XMLoadByteN2(_In_ const XMBYTEN2* pSource);
XMVECTOR        XMLoadByte2(_In_ const XMBYTE2* pSource);
XMVECTOR        XMLoadUByteN2(_In_ const XMUBYTEN2* pSource);
XMVECTOR        XMLoadUByte2(_In_ const XMUBYTE2* pSource);

XMVECTOR        XMLoadU565(_In_ const XMU565* pSource);
XMVECTOR        XMLoadFloat3PK(_In_ const XMFLOAT3PK* pSource);
XMVECTOR        XMLoadFloat3SE(_In_ const XMFLOAT3SE* pSource);

XMVECTOR        XMLoadHalf4(_In_ const XMHALF4* pSource);
XMVECTOR        XMLoadShortN4(_In_ const XMSHORTN4* pSource);
XMVECTOR        XMLoadShort4(_In_ const XMSHORT4* pSource);
XMVECTOR        XMLoadUShortN4(_In_ const XMUSHORTN4* pSource);
XMVECTOR        XMLoadUShort4(_In_ const XMUSHORT4* pSource);
XMVECTOR        XMLoadXDecN4(_In_ const XMXDECN4* pSource);
XMVECTOR        XMLoadXDec4(_In_ const XMXDEC4* pSource);
XMVECTOR        XMLoadDecN4(_In_ const XMDECN4* pSource);
XMVECTOR        XMLoadDec4(_In_ const XMDEC4* pSource);
XMVECTOR        XMLoadUDecN4(_In_ const XMUDECN4* pSource);
XMVECTOR        XMLoadUDec4(_In_ const XMUDEC4* pSource);
XMVECTOR        XMLoadByteN4(_In_ const XMBYTEN4* pSource);
XMVECTOR        XMLoadByte4(_In_ const XMBYTE4* pSource);
XMVECTOR        XMLoadUByteN4(_In_ const XMUBYTEN4* pSource);
XMVECTOR        XMLoadUByte4(_In_ const XMUBYTE4* pSource);
XMVECTOR        XMLoadUNibble4(_In_ const XMUNIBBLE4* pSource);
XMVECTOR        XMLoadU555(_In_ const XMU555* pSource);

///begin_xbox360
namespace Xbox
{

XMVECTOR        XMLoadHenDN3(_In_ const XMHENDN3* pSource);
XMVECTOR        XMLoadHenD3(_In_ const XMHEND3* pSource);
XMVECTOR        XMLoadUHenDN3(_In_ const XMUHENDN3* pSource);
XMVECTOR        XMLoadUHenD3(_In_ const XMUHEND3* pSource);
XMVECTOR        XMLoadDHenN3(_In_ const XMDHENN3* pSource);
XMVECTOR        XMLoadDHen3(_In_ const XMDHEN3* pSource);
XMVECTOR        XMLoadUDHenN3(_In_ const XMUDHENN3* pSource);
XMVECTOR        XMLoadUDHen3(_In_ const XMUDHEN3* pSource);

XMVECTOR        XMLoadXIcoN4(_In_ const XMXICON4* pSource);
XMVECTOR        XMLoadXIco4(_In_ const XMXICO4* pSource);
XMVECTOR        XMLoadIcoN4(_In_ const XMICON4* pSource);
XMVECTOR        XMLoadIco4(_In_ const XMICO4* pSource);
XMVECTOR        XMLoadUIcoN4(_In_ const XMUICON4* pSource);
XMVECTOR        XMLoadUIco4(_In_ const XMUICO4* pSource);

}; // namespace Xbox
///end_xbox360

/****************************************************************************
 *
 * Store operations
 *
 ****************************************************************************/

void            XMStoreColor(_Out_ XMCOLOR* pDestination, _In_ FXMVECTOR V);

void            XMStoreHalf2(_Out_ XMHALF2* pDestination, _In_ FXMVECTOR V);
void            XMStoreShortN2(_Out_ XMSHORTN2* pDestination, _In_ FXMVECTOR V);
void            XMStoreShort2(_Out_ XMSHORT2* pDestination, _In_ FXMVECTOR V);
void            XMStoreUShortN2(_Out_ XMUSHORTN2* pDestination, _In_ FXMVECTOR V);
void            XMStoreUShort2(_Out_ XMUSHORT2* pDestination, _In_ FXMVECTOR V);
void            XMStoreByteN2(_Out_ XMBYTEN2* pDestination, _In_ FXMVECTOR V);
void            XMStoreByte2(_Out_ XMBYTE2* pDestination, _In_ FXMVECTOR V);
void            XMStoreUByteN2(_Out_ XMUBYTEN2* pDestination, _In_ FXMVECTOR V);
void            XMStoreUByte2(_Out_ XMUBYTE2* pDestination, _In_ FXMVECTOR V);

void            XMStoreU565(_Out_ XMU565* pDestination, _In_ FXMVECTOR V);
void            XMStoreFloat3PK(_Out_ XMFLOAT3PK* pDestination, _In_ FXMVECTOR V);
void            XMStoreFloat3SE(_Out_ XMFLOAT3SE* pDestination, _In_ FXMVECTOR V);

void            XMStoreHalf4(_Out_ XMHALF4* pDestination, _In_ FXMVECTOR V);
void            XMStoreShortN4(_Out_ XMSHORTN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreShort4(_Out_ XMSHORT4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUShortN4(_Out_ XMUSHORTN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUShort4(_Out_ XMUSHORT4* pDestination, _In_ FXMVECTOR V);
void            XMStoreXDecN4(_Out_ XMXDECN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreXDec4(_Out_ XMXDEC4* pDestination, _In_ FXMVECTOR V);
void            XMStoreDecN4(_Out_ XMDECN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreDec4(_Out_ XMDEC4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUDecN4(_Out_ XMUDECN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUDec4(_Out_ XMUDEC4* pDestination, _In_ FXMVECTOR V);
void            XMStoreByteN4(_Out_ XMBYTEN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreByte4(_Out_ XMBYTE4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUByteN4(_Out_ XMUBYTEN4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUByte4(_Out_ XMUBYTE4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUNibble4(_Out_ XMUNIBBLE4* pDestination, _In_ FXMVECTOR V);
void            XMStoreU555(_Out_ XMU555* pDestination, _In_ FXMVECTOR V);

///begin_xbox360
namespace Xbox
{

void            XMStoreHenDN3(_Out_ XMHENDN3* pDestination, _In_ FXMVECTOR V);
void            XMStoreHenD3(_Out_ XMHEND3* pDestination, _In_ FXMVECTOR V);
void            XMStoreUHenDN3(_Out_ XMUHENDN3* pDestination, _In_ FXMVECTOR V);
void            XMStoreUHenD3(_Out_ XMUHEND3* pDestination, _In_ FXMVECTOR V);
void            XMStoreDHenN3(_Out_ XMDHENN3* pDestination, _In_ FXMVECTOR V);
void            XMStoreDHen3(_Out_ XMDHEN3* pDestination, _In_ FXMVECTOR V);
void            XMStoreUDHenN3(_Out_ XMUDHENN3* pDestination, _In_ FXMVECTOR V);
void            XMStoreUDHen3(_Out_ XMUDHEN3* pDestination, _In_ FXMVECTOR V);

void            XMStoreXIcoN4(_Out_ XMXICON4* pDestination, _In_ FXMVECTOR V);
void            XMStoreXIco4(_Out_ XMXICO4* pDestination, _In_ FXMVECTOR V);
void            XMStoreIcoN4(_Out_ XMICON4* pDestination, _In_ FXMVECTOR V);
void            XMStoreIco4(_Out_ XMICO4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUIcoN4(_Out_ XMUICON4* pDestination, _In_ FXMVECTOR V);
void            XMStoreUIco4(_Out_ XMUICO4* pDestination, _In_ FXMVECTOR V);

}; // namespace Xbox
///end_xbox360

/****************************************************************************
 *
 * Implementation
 *
 ****************************************************************************/

#pragma warning(push)
#pragma warning(disable:4068 4214 4204 4365 4616 6001)

#pragma prefast(push)
#pragma prefast(disable : 25000, "FXMVECTOR is 16 bytes")

#include "DirectXPackedVector.inl"

#pragma prefast(pop)
#pragma warning(pop)

}; // namespace PackedVector

}; // namespace DirectX


