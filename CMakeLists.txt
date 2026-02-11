cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif()

project(CNPY)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CNPY_WITH_ZLIB "Enable .npz support via zlib" ON)

add_library(cnpy INTERFACE)
target_include_directories(cnpy
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

if(CNPY_WITH_ZLIB)
    find_package(ZLIB)
    if(ZLIB_FOUND)
        target_link_libraries(cnpy INTERFACE ZLIB::ZLIB)
        target_compile_definitions(cnpy INTERFACE CNPY_HAS_ZLIB=1)
    else()
        message(FATAL_ERROR
            "CNPY_WITH_ZLIB is ON, but zlib was not found. Install the zlib development package (e.g., zlib1g-dev), or configure with -DCNPY_WITH_ZLIB=OFF for .npy-only support."
        )
    endif()
else()
    target_compile_definitions(cnpy INTERFACE CNPY_HAS_ZLIB=0)
endif()

install(TARGETS cnpy EXPORT cnpyTargets)
install(DIRECTORY include/ DESTINATION include)
install(PROGRAMS tools/mat2npz tools/npy2mat tools/npz2mat DESTINATION bin)

add_executable(example1 examples/example1.cpp)
target_link_libraries(example1 PRIVATE cnpy)

enable_testing()

file(GLOB CNPY_TEST_SOURCES CONFIGURE_DEPENDS tests/test_*.cpp)
foreach(test_source ${CNPY_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} PRIVATE cnpy)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
