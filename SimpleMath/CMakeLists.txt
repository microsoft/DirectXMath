# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

cmake_minimum_required (VERSION 3.21)

if(DEFINED DIRECTXMATH_VERSION)
   set(SIMPLEMATH_VERSION ${DIRECTXMATH_VERSION})
else()
   set(SIMPLEMATH_VERSION 1.0.0)
endif()

project(SimpleMath
  VERSION ${SIMPLEMATH_VERSION}
  DESCRIPTION "Simplified C++ Math wrapper for DirectXMath"
  HOMEPAGE_URL "https://go.microsoft.com/fwlink/?LinkID=248929"
  LANGUAGES CXX)

if(PROJECT_IS_TOP_LEVEL)
  message(FATAL_ERROR "SimpleMath should be built by the main CMakeLists using BUILD_SIMPLEMATH")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_CURRENT_LIST_DIR}/../build/CompilerAndLinker.cmake)

#--- Library
set(LIBRARY_HEADERS SimpleMath.h SimpleMath.inl)
set(LIBRARY_SOURCES SimpleMath.cpp)

option(BUILD_SHARED_LIBS "Build SimpleMath as a shared library" OFF)

add_library(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE ${LIBRARY_HEADERS} ${LIBRARY_SOURCES})

if(WIN32 AND BUILD_SHARED_LIBS)
  message(STATUS "Build SimpleMath library as a DLL")

  configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/build/SimpleMath.rc.in"
      "${CMAKE_CURRENT_BINARY_DIR}/SimpleMath.rc" @ONLY)

  target_sources(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/SimpleMath.rc")

  target_compile_definitions(${PROJECT_NAME} PRIVATE DIRECTX_TOOLKIT_EXPORT)
  target_compile_definitions(${PROJECT_NAME} INTERFACE DIRECTX_TOOLKIT_IMPORT)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PRIVATE DirectXMath)

#--- Package
include(CMakePackageConfigHelpers)

string(TOLOWER ${PROJECT_NAME} PACKAGE_NAME)
cmake_path(GET CMAKE_CURRENT_LIST_DIR PARENT_PATH DIRECTXMATH_PATH)

write_basic_package_version_file(
  ${PACKAGE_NAME}-config-version.cmake
  VERSION ${SIMPLEMATH_VERSION}
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

