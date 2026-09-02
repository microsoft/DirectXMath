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

#--- Library
set(LIBRARY_HEADERS SimpleMath.h SimpleMath.inl)
set(LIBRARY_SOURCES SimpleMath.cpp)

add_library(${PROJECT_NAME} STATIC)

target_sources(${PROJECT_NAME} PRIVATE ${LIBRARY_HEADERS} ${LIBRARY_SOURCES})

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

#--- Compiler switches
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /Wall /EHsc /GR /Gd /GS /Zc:forScope /Zc:inline /Zc:wchar_t "$<$<NOT:$<CONFIG:DEBUG>>:/guard:cf>")

    if((MSVC_VERSION GREATER_EQUAL 1928)
       AND (CMAKE_SIZEOF_VOID_P EQUAL 8)
       AND ((NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")) OR (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13.0)))
        target_compile_options(${PROJECT_NAME} PRIVATE "$<$<NOT:$<CONFIG:DEBUG>>:/guard:ehcont>")
    endif()
else()
    target_compile_definitions(${PROJECT_NAME} PRIVATE $<IF:$<CONFIG:DEBUG>,_DEBUG,NDEBUG>)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16.0)
        target_compile_options(${PROJECT_NAME} PRIVATE /ZH:SHA_256 "-Wno-unsafe-buffer-usage")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wno-ignored-attributes)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(${PROJECT_NAME} PRIVATE /Zc:__cplusplus /Zc:inline /fp:fast /Qdiag-disable:161)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${PROJECT_NAME} PRIVATE
         /sdl /fp:fast
         /wd4061 /wd4365 /wd4514 /wd4571 /wd4668 /wd4710 /wd4820 /wd5039 /wd5045)

    if(CMAKE_INTERPROCEDURAL_OPTIMIZATION)
        target_compile_options(${PROJECT_NAME} PRIVATE $<$<NOT:$<CONFIG:DEBUG>>:/Gy /Gw>)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.10)
        target_compile_options(${PROJECT_NAME} PRIVATE /permissive-)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.14)
        target_compile_options(${PROJECT_NAME} PRIVATE /Zc:__cplusplus)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.15)
        target_compile_options(${PROJECT_NAME} PRIVATE /JMC-)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.24)
        target_compile_options(${PROJECT_NAME} PRIVATE /ZH:SHA_256)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.26)
        target_compile_options(${PROJECT_NAME} PRIVATE /Zc:preprocessor /wd5105)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.28)
        target_compile_options(${PROJECT_NAME} PRIVATE /Zc:lambda)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.29)
        target_compile_options(${PROJECT_NAME} PRIVATE /external:W4)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.34)
        target_compile_options(${PROJECT_NAME} PRIVATE /wd5262 /wd5264)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.50)
      target_compile_options(${PROJECT_NAME} PRIVATE /wd4865)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.35)
        if(CMAKE_INTERPROCEDURAL_OPTIMIZATION)
          target_compile_options(${PROJECT_NAME} PRIVATE $<$<NOT:$<CONFIG:DEBUG>>:/Zc:checkGwOdr>)
        endif()

        target_compile_options(${PROJECT_NAME} PRIVATE $<$<VERSION_GREATER_EQUAL:${CMAKE_VS_WINDOWS_TARGET_PLATFORM_VERSION},10.0.22000>:/Zc:templateScope>)
    endif()

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.50)
      target_compile_options(${PROJECT_NAME} PRIVATE /Zc:u8EscapeEncoding)
    endif()
endif()


