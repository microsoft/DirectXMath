# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

cmake_minimum_required (VERSION 3.21)

if(DEFINED DIRECTXMATH_VERSION)
   set(SHMATH_VERSION ${DIRECTXMATH_VERSION})
else()
   set(SHMATH_VERSION 1.0.6)
endif()

project(DirectXSH
  VERSION ${SHMATH_VERSION}
  DESCRIPTION "C++ Spherical Harmonics Math Library"
  HOMEPAGE_URL "https://go.microsoft.com/fwlink/?LinkID=615560"
  LANGUAGES CXX)

if(PROJECT_IS_TOP_LEVEL)
  message(FATAL_ERROR "SHMath should be built by the main CMakeLists using BUILD_SHMATH")
endif()

option(BUILD_DX11 "Build with DirectX11 support" OFF)

option(BUILD_DX12 "Build with DirectX12 support" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_CURRENT_LIST_DIR}/../build/CompilerAndLinker.cmake)


#--- Library
set(LIBRARY_HEADERS DirectXSH.h)
set(LIBRARY_SOURCES DirectXSH.cpp)

if(WIN32 AND BUILD_DX11)
  set(LIBRARY_SOURCES ${LIBRARY_SOURCES} DirectXSHD3D11.cpp)
endif()

if(WIN32 AND BUILD_DX12)
  set(LIBRARY_SOURCES ${LIBRARY_SOURCES} DirectXSHD3D12.cpp)
endif()

add_library(${PROJECT_NAME} STATIC)

target_sources(${PROJECT_NAME} PRIVATE ${LIBRARY_HEADERS} ${LIBRARY_SOURCES})

if(COMPILER_DEFINES)
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${COMPILER_DEFINES})
endif()
if(COMPILER_SWITCHES)
    target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILER_SWITCHES})
endif()
if(LINKER_SWITCHES)
    target_link_options(${PROJECT_NAME} PRIVATE ${LINKER_SWITCHES})
endif()

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PRIVATE DirectXMath)

if(MINGW AND BUILD_DX12)
    find_package(directx-headers CONFIG REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC Microsoft::DirectX-Headers)
    target_compile_definitions(${PROJECT_NAME} PUBLIC USING_DIRECTX_HEADERS)
else()
    find_package(directx-headers CONFIG QUIET)
endif()

if(directx-headers_FOUND)
    message(STATUS "Using DirectX-Headers package")
    target_link_libraries(${PROJECT_NAME} PRIVATE Microsoft::DirectX-Headers)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USING_DIRECTX_HEADERS)
endif()

#--- Package
include(CMakePackageConfigHelpers)

string(TOLOWER ${PROJECT_NAME} PACKAGE_NAME)
cmake_path(GET CMAKE_CURRENT_LIST_DIR PARENT_PATH DIRECTXMATH_PATH)

write_basic_package_version_file(
  ${PACKAGE_NAME}-config-version.cmake
  VERSION ${SHMATH_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

configure_package_config_file(${DIRECTXMATH_PATH}/build/DirectXMath-config.cmake.in
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


if((NOT WIN32) AND (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/sal/sal.h"))
    file(DOWNLOAD
        https://raw.githubusercontent.com/dotnet/runtime/v9.0.2/src/coreclr/pal/inc/rt/sal.h
        "${CMAKE_CURRENT_BINARY_DIR}/sal/sal.h"
        EXPECTED_HASH SHA512=8085f67bfa4ce01ae89461cadf72454a9552fde3f08b2dcc3de36b9830e29ce7a6192800f8a5cb2a66af9637be0017e85719826a4cfdade508ae97f319e0ee8e
        STATUS DOWNLOAD_STATUS
    )
    list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
    list(GET DOWNLOAD_STATUS 1 DOWNLOAD_ERROR)
    if(NOT (STATUS_CODE EQUAL 0))
        message(FATAL_ERROR "Failed to download sal.h: ${DOWNLOAD_ERROR}")
    endif()

    target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/sal")
endif()
