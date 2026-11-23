cmake_minimum_required(VERSION 3.21)

project(cnpy
    VERSION 2.0.0
    DESCRIPTION "Read and write NumPy .npy and .npz files from C++20"
    HOMEPAGE_URL "https://github.com/lelapinagile/cnpy"
    LANGUAGES CXX)

option(CNPY_BUILD_TESTS "Build the unit and NumPy interop tests" ${PROJECT_IS_TOP_LEVEL})
option(CNPY_BUILD_EXAMPLES "Build the examples" ${PROJECT_IS_TOP_LEVEL})
option(CNPY_BUILD_BENCHMARKS "Build the benchmark" OFF)
option(CNPY_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(CNPY_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(CNPY_INSTALL "Generate install rules" ${PROJECT_IS_TOP_LEVEL})

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin) # DLLs next to the executables on Windows
endif()

if(CNPY_ENABLE_SANITIZERS AND NOT MSVC)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    add_link_options(-fsanitize=address,undefined)
endif()

find_package(ZLIB REQUIRED)

# ---------------------------------------------------------------------------------------- library

add_library(cnpy
    src/array.cpp
    src/io.cpp
    src/npy.cpp
    src/npz.cpp)
add_library(cnpy::cnpy ALIAS cnpy)

target_compile_features(cnpy PUBLIC cxx_std_20)
set_target_properties(cnpy PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

include(GenerateExportHeader)
generate_export_header(cnpy
    BASE_NAME CNPY
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/include/cnpy/export.hpp)
get_target_property(_cnpy_type cnpy TYPE)
if(_cnpy_type STREQUAL "STATIC_LIBRARY")
    target_compile_definitions(cnpy PUBLIC CNPY_STATIC_DEFINE)
endif()

include(GNUInstallDirs)
target_include_directories(cnpy PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(cnpy PRIVATE ZLIB::ZLIB)

if(UNIX)
    target_compile_definitions(cnpy PRIVATE _FILE_OFFSET_BITS=64)
endif()

# Warnings for everything cnpy builds itself.
add_library(cnpy_warnings INTERFACE)
if(MSVC)
    target_compile_options(cnpy_warnings INTERFACE /W4 /permissive- /utf-8 /Zc:__cplusplus
        $<$<BOOL:${CNPY_WARNINGS_AS_ERRORS}>:/WX>)
else()
    target_compile_options(cnpy_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wnon-virtual-dtor -Wold-style-cast -Wcast-align
        $<$<BOOL:${CNPY_WARNINGS_AS_ERRORS}>:-Werror>)
endif()
target_link_libraries(cnpy PRIVATE $<BUILD_INTERFACE:cnpy_warnings>)

# ---------------------------------------------------------------------------------------- extras

if(CNPY_BUILD_EXAMPLES)
    add_executable(cnpy_example examples/example.cpp)
    target_link_libraries(cnpy_example PRIVATE cnpy::cnpy cnpy_warnings)
endif()

if(CNPY_BUILD_BENCHMARKS)
    add_executable(cnpy_bench bench/bench.cpp)
    target_link_libraries(cnpy_bench PRIVATE cnpy::cnpy)
endif()

if(CNPY_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------------------- install

if(CNPY_INSTALL)
    include(CMakePackageConfigHelpers)

    install(TARGETS cnpy EXPORT cnpyTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(DIRECTORY include/cnpy DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/cnpy/export.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cnpy)

    set(_cnpy_cmake_dir ${CMAKE_INSTALL_LIBDIR}/cmake/cnpy)
    install(EXPORT cnpyTargets NAMESPACE cnpy:: DESTINATION ${_cnpy_cmake_dir})
    configure_package_config_file(cmake/cnpyConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/cnpyConfig.cmake
        INSTALL_DESTINATION ${_cnpy_cmake_dir})
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/cnpyConfigVersion.cmake
        COMPATIBILITY SameMajorVersion)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cnpyConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cnpyConfigVersion.cmake
        DESTINATION ${_cnpy_cmake_dir})
endif()
