# GitHub Copilot Instructions

These instructions define how GitHub Copilot should assist with this project. The goal is to ensure consistent, high-quality code generation aligned with our conventions, stack, and best practices.

## Context

- **Project Type**: Math Library / DirectX / Direct3D
- **Project Name**: DirectXMath SIMD C++ linear algebra library
- **Language**: C++
- **Framework / Libraries**: STL / CMake / CTest

## Getting Started

- See the Getting Started guide on [Microsoft Learn](https://learn.microsoft.com/windows/win32/dxmath/pg-xnamath-getting-started).
- The recommended way to integrate *DirectXMath* into your project is by using the *vcpkg* Package Manager.
- You can make use of the nuget.org package **directxmath**.
- You can also use the library source code directly in your project or as a git submodule.
- The Windows SDK includes the DirectXMath library, although that version is not as up-to-date as other integration methods.

## General Guidelines

- **Code Style**: The project uses an .editorconfig file to enforce coding standards. Follow the rules defined in `.editorconfig` for indentation, line endings, and other formatting. Additional information can be found on the wiki at [Implementation](https://github.com/microsoft/DirectXMath/wiki/Implementation). The code requires C++11/C++14 features.
- **Documentation**: The project provides documentation on [Microsoft Learn](https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal) with additional wiki pages available on [GitHub](https://github.com/microsoft/DirectXMath/wiki/).
- **Error Handling**: The majority of functions have no error conditions and do not throw C++ exceptions which is why they are marked `noexcept`. A few functions have `bool` results to indicate success or failure.
- **Testing**: Unit tests for this project are implemented in this repository [Test Suite](https://github.com/walbourn/directxmathtest/) and can be run using CTest per the instructions at [Test Documentation](https://github.com/walbourn/directxmathtest/wiki).
- **Security**: This project uses secure coding practices from the Microsoft Secure Coding Guidelines, and is subject to the `SECURITY.md` file in the root of the repository.
- **Dependencies**: The project has minimal dependencies, primarily relying on compiler intrinsics. It is designed to be self-contained and portable across different platforms and toolsets.
- **Continuous Integration**: This project implements GitHub Actions for continuous integration, ensuring that all code changes are tested and validated before merging. This includes building the project for a number of configurations and toolsets, running unit tests, and static code analysis including GitHub super-linter, CodeQL, and MSVC Code Analysis.
- **Code of Conduct**: The project adheres to the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/). All contributors are expected to follow this code of conduct in all interactions related to the project.

## File Structure

```txt
.azuredevops/ # Azure DevOps pipeline configuration and policy files.
.github/      # GitHub Actions workflow files and linter configuration files.
.nuget/       # NuGet package configuration files.
build/        # Miscellaneous build files and scripts.
Inc/          # DirectXMath public and implementation files. The library is header-only, so all files are in this directory.
Extensions/   # Extensions to the DirectXMath library with standalone SSE-level specific functions for runtime selection of SIMD instruction set.
MatrixStack/  # D3DX9-like matrix stack implementation for DirectXMath.
SHMath/       # Spherical harmonic functions using DirectXMath.
Stereo3D/     # Stereo 3D projection matrix functions using DirectXMath created for HoloLens.
XDSP/         # Digital Signal Processing (DSP) functions using DirectXMath.
Tests/        # Tests are designed to be cloned from a separate repository at this location.
```

> The `Extensions` are not needed if building the library using `/arch:AVX` or `/arch:AVX2` which causes the DirectXMath library to build utilizing the additional SIMD instructions.

## References

- [Source git repository on GitHub](https://github.com/microsoft/DirectXMath.git)
- [DirectXMath wiki git repository on GitHub](https://github.com/microsoft/DirectXMath.wiki.git)
- [DirectXMath test suite git repository on GitHub](https://github.com/walbourn/directxmathtest.wiki.git).
- [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
- [Microsoft Secure Coding Guidelines](https://learn.microsoft.com/en-us/security/develop/secure-coding-guidelines)
- [CMake Documentation](https://cmake.org/documentation/)
- [VCPKG Documentation](https://learn.microsoft.com/vcpkg/)
- [DirectXMath Documentation](https://learn.microsoft.com/windows/win32/dxmath/directxmath-portal)
- [Games for Windows and the DirectX SDK blog - March 2012](https://walbourn.github.io/introducing-directxmath/)

## No speculation

When creating documentation:

### Document Only What Exists

- Only document features, patterns, and decisions that are explicitly present in the source code.
- Only include configurations and requirements that are clearly specified.
- Do not make assumptions about implementation details.

### Handle Missing Information

- Ask the user questions to gather missing information.
- Document gaps in current implementation or specifications.
- List open questions that need to be addressed.

### Source Material

- Always cite the specific source file and line numbers for documented features.
- Link directly to relevant source code when possible.
- Indicate when information comes from requirements vs. implementation.

### Verification Process

- Review each documented item against source code whenever related to the task.
- Remove any speculative content.
- Ensure all documentation is verifiable against the current state of the codebase.

## Cross-platform Support Notes

- The code supports building for Windows and Linux.
- Portability and conformance of the code is validated by building with Visual C++, clang/LLVM for Windows, MinGW, and GCC for Linux.

## Code Review Instructions

When reviewing code, focus on the following aspects:

- Adherence to coding standards defined in `.editorconfig` and on the [wiki](https://github.com/microsoft/DirectXTK/wiki/Implementation).
- Make coding recommendations based on the *C++ Core Guidelines*.
- Proper use of RAII and smart pointers.
- Correct error handling practices and C++ Exception safety.
- Clarity and maintainability of the code.
- Adequate comments where necessary.
- Public interfaces located in `Inc\*.h` should be clearly defined and documented on the Microsoft Docs pages.
- Optional functions are available in headers in the `SHMath`, `Stereo3D`, `MatrixStack`, and `XDSP` folders.
- Compliance with the project's architecture and design patterns.
- Ensure that all public functions and classes are covered by unit tests located on [GitHub](https://github.com/walbourn/directxmathtest.git) where applicable. Report any gaps in test coverage.
- Check for performance implications, especially in geometry processing algorithms.
- Provide brutally honest feedback on code quality, design, and potential improvements as needed.
