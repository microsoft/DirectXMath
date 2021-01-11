![DirectX Logo](https://github.com/Microsoft/DirectXMath/wiki/X_jpg.jpg)

# DirectXMath

https://github.com/Microsoft/DirectXMath

Copyright (c) Microsoft Corporation. All rights reserved.

**January 2021**

This package contains the DirectXMath library, an all inline SIMD C++ linear algebra library for use in games and graphics apps.

This code is designed to build with Visual Studio 2017, Visual Studio 2019, or clang/LLVM for Windows. It is recommended that you make use of the latest updates (VS 2017 15.9, or VS 2019 16.4 or later).

These components are designed to work without requiring any content from the legacy DirectX SDK. For details, see [Where is the DirectX SDK?](https://aka.ms/dxsdk).

## Directory Layout

* ``Inc\``

  + DirectXMath Files (in the DirectX C++ namespace)

    * DirectXMath.h - Core library
    * DirectXPackedVector.h - Load/Store functions and types for working with various compressed GPU formats
    * DirectXColors.h - .NET-style Color defines in sRGB color space
    * DirectXCollision.h - Bounding volume collision library

* ``Extentions\``

  + Advanced instruction set variants for guarded codepaths

    * DirectXMathSSE3.h - SSE3
    * DirectXMathBE.h - Supplemental SSE3 (SSSE3)
    * DirectXMathSSE4.h - SSE4.1
    * DirectXMathAVX.h - Advanced Vector Extensions (AVX)
    * DirectXMathAVX2.h - Advanced Vector Extensions 2 (AVX2)
    * DirectXMathF16C.h - Half-precision conversions (F16C)
    * DirectXMathFMA3.h - Fused multiply-accumulate (FMA3)
    * DirectXMathFMA4.h - Fused multiply-accumulate (FMA4)

* ``SHMath\``

  + Spherical Harmonics math functions

    * DirectXSH.h - Header for SHMath functions
    * DirectXSH.cpp, DirectXSHD3D11.cpp, DirectXSHD3D12.cpp - Implementation

* ``XDSP\``

  + XDSP.h - Digital Signal Processing helper functions

## Documentation

Documentation is available on the [Microsoft Docs](https://docs.microsoft.com/en-us/windows/desktop/dxmath/directxmath-portal). Additional information can be found on the [project wiki](https://github.com/microsoft/DirectXMath/wiki).

## Compiler support

Officially the library is supported with Microsoft Visual C++ and clang/LLVM. It should also compile with the Intel C++ compiler, GCC, and MinGW compilers.

To build for non-Windows platforms, you need to provide a ``sal.h`` header in your include path. You can obtain an open source version from [GitHub](https://github.com/dotnet/corert/blob/master/src/Native/inc/unix/sal.h).

## Notices

All content and source code for this package are subject to the terms of the [MIT License](http://opensource.org/licenses/MIT).

For the latest version of DirectXMath, bug reports, etc. please visit the project site on [GitHub](https://github.com/microsoft/DirectXMath).

## Contributing

This project welcomes contributions and suggestions. Most contributions require you to agree to a Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/). For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## Trademarks

This project may contain trademarks or logos for projects, products, or services. Authorized use of Microsoft trademarks or logos is subject to and must follow [Microsoft's Trademark & Brand Guidelines](https://www.microsoft.com/en-us/legal/intellectualproperty/trademarks/usage/general). Use of Microsoft trademarks or logos in modified versions of this project must not cause confusion or imply Microsoft sponsorship. Any use of third-party trademarks or logos are subject to those third-party's policies.
