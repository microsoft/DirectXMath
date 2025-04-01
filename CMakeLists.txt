# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

cmake_minimum_required (VERSION 3.20)

set(DIRECTXMATH_VERSION 3.20)

if(WINDOWS_STORE OR (DEFINED XBOX_CONSOLE_TARGET))
  set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
endif()

project(DirectXMath
  VERSION ${DIRECTXMATH_VERSION}
  DESCRIPTION "DirectXMath SIMD C++ math library"
  HOMEPAGE_URL "https://go.microsoft.com/fwlink/?LinkID=615560"
  LANGUAGES CXX)

option(BUILD_XDSP "Build XDSP math" OFF)

option(BUILD_SHMATH "Build Spherical Harmonics math" OFF)

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
    Inc/DirectXPackedVector.inl)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Inc>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

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
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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

#--- Optional extension libraries
if(BUILD_XDSP)
  message(STATUS "Including XDSP Digital Signal Processing (DSP)")
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/XDSP)
endif()

if(BUILD_SHMATH)
  message(STATUS "Including C++ Spherical Harmonics Math")
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/SHMath)
endif()

#--- Test suite
if((NOT WINDOWS_STORE) AND (NOT (DEFINED XBOX_CONSOLE_TARGET)))
  include(CTest)
  if(BUILD_TESTING AND (EXISTS "${CMAKE_CURRENT_LIST_DIR}/Tests/CMakeLists.txt"))
    enable_testing()
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Tests)
  endif()
endif()
