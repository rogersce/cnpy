cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif()

project(CNPY)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ZLIB REQUIRED)

add_library(cnpy INTERFACE)
target_include_directories(cnpy
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(cnpy INTERFACE ZLIB::ZLIB)

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
