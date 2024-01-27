# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

cmake_minimum_required (VERSION 3.20)

set(DIRECTXMATH_VERSION 3.1.8)

project(DirectXMath
  VERSION ${DIRECTXMATH_VERSION}
  DESCRIPTION "DirectXMath SIMD C++ math library"
  HOMEPAGE_URL "https://go.microsoft.com/fwlink/?LinkID=615560"
  LANGUAGES CXX)

include(GNUInstallDirs)

#--- Library
set(LIBRARY_HEADERS
    Inc/DirectXCollision.h
    Inc/DirectXCollision.inl
    Inc/DirectXColors.h
    Inc/DirectXMath.h
    Inc/DirectXMathConvert.inl
    Inc/DirectXMathMatrix.inl
    Inc/DirectXMathMisc.inl
    Inc/DirectXMathVector.inl
    Inc/DirectXPackedVector.h
    Inc/DirectXPackedVector.inl
	${CMAKE_CURRENT_BINARY_DIR}/DirectXMathConfig.h)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Inc>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/directxmath>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

#--- configure DirectXMathConfig.h

#---------- 2 parameter __cpuid from intrin.h
file(WRITE ${CMAKE_BINARY_DIR}/CMakeTmp/testcpuid.cpp
"#include <intrin.h>
int main () {
	int CPUInfo[4] = { -1 };
	__cpuid(CPUInfo, 0);
	return(0);
}")

try_compile(USE_INTRIN_H_CPUID
   ${CMAKE_BINARY_DIR}/CMakeTmp
   ${CMAKE_BINARY_DIR}/CMakeTmp/testcpuid.cpp)

if(USE_INTRIN_H_CPUID)
   set(USE_CPUID_H_CPUID FALSE)
   message(STATUS "Configure __cpuid:  intrin.h - 2 parameters")
endif()

if(NOT(USE_INTRIN_H_CPUID))
#--------- 5 parameter __cpuid.h
  file(WRITE ${CMAKE_BINARY_DIR}/CMakeTmp/testcpuid.cpp
"#include <cpuid.h>
int main () {
	int CPUInfo[4] = { -1 };
	__cpuid(0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
	return(0);
}")

  try_compile(USE_CPUID_H_CPUID
    ${CMAKE_BINARY_DIR}/CMakeTmp
    ${CMAKE_BINARY_DIR}/CMakeTmp/testcpuid.cpp)

  if(USE_CPUID_H_CPUID)
     message(STATUS "Configure __cpuid:  cpuid.h - 5 parameters")
  endif()
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/build/DirectXMathConfig.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/DirectXMathConfig.h" @ONLY)

#--- Package
include(CMakePackageConfigHelpers)

string(TOLOWER ${PROJECT_NAME} PACKAGE_NAME)

write_basic_package_version_file(
  ${PACKAGE_NAME}-config-version.cmake
  VERSION ${DIRECTXMATH_VERSION}
  COMPATIBILITY AnyNewerVersion
  ARCH_INDEPENDENT)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/build/${PROJECT_NAME}-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PACKAGE_NAME})

install(EXPORT ${PROJECT_NAME}-targets
  FILE ${PROJECT_NAME}-targets.cmake
  NAMESPACE Microsoft::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PACKAGE_NAME})

install(FILES ${LIBRARY_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/directxmath)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PACKAGE_NAME})

# Create pkg-config file
include(build/JoinPaths.cmake)
# from: https://github.com/jtojnar/cmake-snips#concatenating-paths-when-building-pkg-config-files
join_paths(DIRECTXMATH_INCLUDEDIR_FOR_PKG_CONFIG "\${prefix}" "${CMAKE_INSTALL_INCLUDEDIR}")
join_paths(DIRECTXMATH_LIBDIR_FOR_PKG_CONFIG "\${prefix}"     "${CMAKE_INSTALL_LIBDIR}")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/build/DirectXMath.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/DirectXMath.pc" @ONLY)

# Install the pkg-config file
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/DirectXMath.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

#--- Test suite
if (DEFINED VCPKG_TARGET_ARCHITECTURE)
    set(DXMATH_ARCHITECTURE ${VCPKG_TARGET_ARCHITECTURE})
elseif(CMAKE_GENERATOR_PLATFORM MATCHES "^[Ww][Ii][Nn]32$")
    set(DXMATH_ARCHITECTURE x86)
elseif(CMAKE_GENERATOR_PLATFORM MATCHES "^[Xx]64$")
    set(DXMATH_ARCHITECTURE x64)
elseif(CMAKE_GENERATOR_PLATFORM MATCHES "^[Aa][Rr][Mm]$")
    set(DXMATH_ARCHITECTURE arm)
elseif(CMAKE_GENERATOR_PLATFORM MATCHES "^[Aa][Rr][Mm]64$")
    set(DXMATH_ARCHITECTURE arm64)
elseif(NOT DXMATH_ARCHITECTURE)
    set(DXMATH_ARCHITECTURE "x64")
endif()
set(DXMATH_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}")

#--- Test suite
include(CTest)
if(BUILD_TESTING AND WIN32 AND (NOT WINDOWS_STORE) AND (EXISTS "${CMAKE_CURRENT_LIST_DIR}/Tests/CMakeLists.txt"))
  enable_testing()
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Tests)
endif()
