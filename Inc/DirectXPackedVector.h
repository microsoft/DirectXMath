//-------------------------------------------------------------------------------------
// DirectXPackedVector.h -- SIMD C++ Math library
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//  
// Copyright (c) Microsoft Corporation. All rights reserved.
//
// http://go.microsoft.com/fwlink/?LinkID=615560
//-------------------------------------------------------------------------------------

#pragma once

#include "DirectXMath.h"

namespace DirectX
{
    
namespace PackedVector
{

#pragma warning(push)
#pragma warning(disable:4201 4365 4324 4996)
// C4201: nonstandard extension used
// C4365: Off by default noise
// C4324: alignment padding warnings
// C4996: deprecation warnings

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

    XMCOLOR() =default;
    XM_CONSTEXPR XMCOLOR(uint32_t Color) : c(Color) {}
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

    XMHALF2() =default;
    explicit XM_CONSTEXPR XMHALF2(uint32_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMHALF2(HALF _x, HALF _y) : x(_x), y(_y) {}
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

    XMSHORTN2() =default;
    explicit XM_CONSTEXPR XMSHORTN2(uint32_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMSHORTN2(int16_t _x, int16_t _y) : x(_x), y(_y) {}
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

    XMSHORT2() =default;
    explicit XM_CONSTEXPR XMSHORT2(uint32_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMSHORT2(int16_t _x, int16_t _y) : x(_x), y(_y) {}
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

    XMUSHORTN2() =default;
    explicit XM_CONSTEXPR XMUSHORTN2(uint32_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUSHORTN2(uint16_t _x, uint16_t _y) : x(_x), y(_y) {}
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

    XMUSHORT2() =default;
    explicit XM_CONSTEXPR XMUSHORT2(uint32_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUSHORT2(uint16_t _x, uint16_t _y) : x(_x), y(_y) {}
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

    XMBYTEN2() =default;
    explicit XM_CONSTEXPR XMBYTEN2(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMBYTEN2(int8_t _x, int8_t _y) : x(_x), y(_y) {}
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

    XMBYTE2() =default;
    explicit XM_CONSTEXPR XMBYTE2(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMBYTE2(int8_t _x, int8_t _y) : x(_x), y(_y) {}
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

    XMUBYTEN2() =default;
    explicit XM_CONSTEXPR XMUBYTEN2(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUBYTEN2(uint8_t _x, uint8_t _y) : x(_x), y(_y) {}
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

    XMUBYTE2() =default;
    explicit XM_CONSTEXPR XMUBYTE2(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUBYTE2(uint8_t _x, uint8_t _y) : x(_x), y(_y) {}
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

    XMU565() =default;
    explicit XM_CONSTEXPR XMU565(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMU565(uint8_t _x, uint8_t _y, uint8_t _z) : x(_x), y(_y), z(_z) {}
    explicit XMU565(_In_reads_(3) const uint8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}
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

    XMFLOAT3PK() =default;
    explicit XM_CONSTEXPR XMFLOAT3PK(uint32_t Packed) : v(Packed) {}
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

    XMFLOAT3SE() =default;
    explicit XM_CONSTEXPR XMFLOAT3SE(uint32_t Packed) : v(Packed) {}
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

    XMHALF4() =default;
    explicit XM_CONSTEXPR XMHALF4(uint64_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMHALF4(HALF _x, HALF _y, HALF _z, HALF _w) : x(_x), y(_y), z(_z), w(_w) {}
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

    XMSHORTN4() =default;
    explicit XM_CONSTEXPR XMSHORTN4(uint64_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMSHORTN4(int16_t _x, int16_t _y, int16_t _z, int16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
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

    XMSHORT4() =default;
    explicit XM_CONSTEXPR XMSHORT4(uint64_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMSHORT4(int16_t _x, int16_t _y, int16_t _z, int16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
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

    XMUSHORTN4() =default;
    explicit XM_CONSTEXPR XMUSHORTN4(uint64_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUSHORTN4(uint16_t _x, uint16_t _y, uint16_t _z, uint16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
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

    XMUSHORT4() =default;
    explicit XM_CONSTEXPR XMUSHORT4(uint64_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUSHORT4(uint16_t _x, uint16_t _y, uint16_t _z, uint16_t _w) : x(_x), y(_y), z(_z), w(_w) {}
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

    XMXDECN4() =default;
    explicit XM_CONSTEXPR XMXDECN4(uint32_t Packed) : v(Packed) {}
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
struct XM_DEPRECATED XMXDEC4
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

    XMXDEC4() =default;
    explicit XM_CONSTEXPR XMXDEC4(uint32_t Packed) : v(Packed) {}
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
struct XM_DEPRECATED XMDECN4
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

    XMDECN4() =default;
    explicit XM_CONSTEXPR XMDECN4(uint32_t Packed) : v(Packed) {}
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
struct XM_DEPRECATED XMDEC4
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

    XMDEC4() =default;
    explicit XM_CONSTEXPR XMDEC4(uint32_t Packed) : v(Packed) {}
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

    XMUDECN4() =default;
    explicit XM_CONSTEXPR XMUDECN4(uint32_t Packed) : v(Packed) {}
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

    XMUDEC4() =default;
    explicit XM_CONSTEXPR XMUDEC4(uint32_t Packed) : v(Packed) {}
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

    XMBYTEN4() =default;
    XM_CONSTEXPR XMBYTEN4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XM_CONSTEXPR XMBYTEN4(uint32_t Packed) : v(Packed) {}
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

    XMBYTE4() =default;
    XM_CONSTEXPR XMBYTE4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XM_CONSTEXPR XMBYTE4(uint32_t Packed) : v(Packed) {}
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

    XMUBYTEN4() =default;
    XM_CONSTEXPR XMUBYTEN4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XM_CONSTEXPR XMUBYTEN4(uint32_t Packed) : v(Packed) {}
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

    XMUBYTE4() =default;
    XM_CONSTEXPR XMUBYTE4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XM_CONSTEXPR XMUBYTE4(uint32_t Packed) : v(Packed) {}
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

    XMUNIBBLE4() =default;
    explicit XM_CONSTEXPR XMUNIBBLE4(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMUNIBBLE4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) : x(_x), y(_y), z(_z), w(_w) {}
    explicit XMUNIBBLE4(_In_reads_(4) const uint8_t *pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
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

    XMU555() =default;
    explicit XM_CONSTEXPR XMU555(uint16_t Packed) : v(Packed) {}
    XM_CONSTEXPR XMU555(uint8_t _x, uint8_t _y, uint8_t _z, bool _w) : x(_x), y(_y), z(_z), w(_w ? 0x1 : 0) {}
    XMU555(_In_reads_(3) const uint8_t *pArray, _In_ bool _w) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(_w ? 0x1 : 0) {}
    XMU555(float _x, float _y, float _z, bool _w);
    XMU555(_In_reads_(3) const float *pArray, _In_ bool _w);

    operator uint16_t () const { return v; }

    XMU555& operator= (const XMU555& U555) { v = U555.v; return *this; }
    XMU555& operator= (uint16_t Packed) { v = Packed; return *this; }
};

#pragma warning(pop)


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

XMVECTOR    XM_CALLCONV     XMLoadColor(_In_ const XMCOLOR* pSource);

XMVECTOR    XM_CALLCONV     XMLoadHalf2(_In_ const XMHALF2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadShortN2(_In_ const XMSHORTN2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadShort2(_In_ const XMSHORT2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUShortN2(_In_ const XMUSHORTN2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUShort2(_In_ const XMUSHORT2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadByteN2(_In_ const XMBYTEN2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadByte2(_In_ const XMBYTE2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUByteN2(_In_ const XMUBYTEN2* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUByte2(_In_ const XMUBYTE2* pSource);

XMVECTOR    XM_CALLCONV     XMLoadU565(_In_ const XMU565* pSource);
XMVECTOR    XM_CALLCONV     XMLoadFloat3PK(_In_ const XMFLOAT3PK* pSource);
XMVECTOR    XM_CALLCONV     XMLoadFloat3SE(_In_ const XMFLOAT3SE* pSource);

XMVECTOR    XM_CALLCONV     XMLoadHalf4(_In_ const XMHALF4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadShortN4(_In_ const XMSHORTN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadShort4(_In_ const XMSHORT4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUShortN4(_In_ const XMUSHORTN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUShort4(_In_ const XMUSHORT4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadXDecN4(_In_ const XMXDECN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUDecN4(_In_ const XMUDECN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUDecN4_XR(_In_ const XMUDECN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUDec4(_In_ const XMUDEC4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadByteN4(_In_ const XMBYTEN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadByte4(_In_ const XMBYTE4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUByteN4(_In_ const XMUBYTEN4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUByte4(_In_ const XMUBYTE4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadUNibble4(_In_ const XMUNIBBLE4* pSource);
XMVECTOR    XM_CALLCONV     XMLoadU555(_In_ const XMU555* pSource);

#pragma warning(push)
#pragma warning(disable : 4996)
// C4996: ignore deprecation warning

XMVECTOR    XM_DEPRECATED XM_CALLCONV XMLoadDecN4(_In_ const XMDECN4* pSource);
XMVECTOR    XM_DEPRECATED XM_CALLCONV XMLoadDec4(_In_ const XMDEC4* pSource);
XMVECTOR    XM_DEPRECATED XM_CALLCONV XMLoadXDec4(_In_ const XMXDEC4* pSource);
#pragma warning(pop)

/****************************************************************************
 *
 * Store operations
 *
 ****************************************************************************/

void    XM_CALLCONV     XMStoreColor(_Out_ XMCOLOR* pDestination, _In_ FXMVECTOR V);

void    XM_CALLCONV     XMStoreHalf2(_Out_ XMHALF2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreShortN2(_Out_ XMSHORTN2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreShort2(_Out_ XMSHORT2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUShortN2(_Out_ XMUSHORTN2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUShort2(_Out_ XMUSHORT2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreByteN2(_Out_ XMBYTEN2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreByte2(_Out_ XMBYTE2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUByteN2(_Out_ XMUBYTEN2* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUByte2(_Out_ XMUBYTE2* pDestination, _In_ FXMVECTOR V);

void    XM_CALLCONV     XMStoreU565(_Out_ XMU565* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreFloat3PK(_Out_ XMFLOAT3PK* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreFloat3SE(_Out_ XMFLOAT3SE* pDestination, _In_ FXMVECTOR V);

void    XM_CALLCONV     XMStoreHalf4(_Out_ XMHALF4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreShortN4(_Out_ XMSHORTN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreShort4(_Out_ XMSHORT4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUShortN4(_Out_ XMUSHORTN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUShort4(_Out_ XMUSHORT4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreXDecN4(_Out_ XMXDECN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUDecN4(_Out_ XMUDECN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUDecN4_XR(_Out_ XMUDECN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUDec4(_Out_ XMUDEC4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreByteN4(_Out_ XMBYTEN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreByte4(_Out_ XMBYTE4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUByteN4(_Out_ XMUBYTEN4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUByte4(_Out_ XMUBYTE4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreUNibble4(_Out_ XMUNIBBLE4* pDestination, _In_ FXMVECTOR V);
void    XM_CALLCONV     XMStoreU555(_Out_ XMU555* pDestination, _In_ FXMVECTOR V);

#pragma warning(push)
#pragma warning(disable : 4996)
// C4996: ignore deprecation warning

void    XM_DEPRECATED XM_CALLCONV XMStoreDecN4(_Out_ XMDECN4* pDestination, _In_ FXMVECTOR V);
void    XM_DEPRECATED XM_CALLCONV XMStoreDec4(_Out_ XMDEC4* pDestination, _In_ FXMVECTOR V);
void    XM_DEPRECATED XM_CALLCONV XMStoreXDec4(_Out_ XMXDEC4* pDestination, _In_ FXMVECTOR V);
#pragma warning(pop)

/****************************************************************************
 *
 * Implementation
 *
 ****************************************************************************/

#pragma warning(push)
#pragma warning(disable:4068 4214 4204 4365 4616 6001 6101)
// C4068/4616: ignore unknown pragmas
// C4214/4204: nonstandard extension used
// C4365: Off by default noise
// C6001/6101: False positives

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 25000, "FXMVECTOR is 16 bytes")
#endif

#include "DirectXPackedVector.inl"

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

#pragma warning(pop)

} // namespace PackedVector

} // namespace DirectX

