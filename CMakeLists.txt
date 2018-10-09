CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

project(CNPY)

# Update submodules as needed
find_package(Git QUIET)
if(GIT_FOUND)
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

option(ENABLE_STATIC "Build static (.a) library" ON)

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
FIND_PACKAGE(LibZip REQUIRED)

add_library(cnpy SHARED "cnpy.cpp")
target_include_directories(cnpy PUBLIC ${LIBZIP_INCLUDE_DIR_ZIP})
target_include_directories(cnpy PUBLIC "extern/libzippp/")
target_link_libraries(cnpy ${LIBZIP_LIBRARIES})
target_compile_features(cnpy PRIVATE cxx_std_14)

install(TARGETS "cnpy" LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if(ENABLE_STATIC)
    add_library(cnpy-static STATIC "cnpy.cpp")
    set_target_properties(cnpy-static PROPERTIES OUTPUT_NAME "cnpy")
    install(TARGETS "cnpy-static" ARCHIVE DESTINATION lib)
endif(ENABLE_STATIC)

install(FILES "cnpy.h" DESTINATION include)
install(FILES "mat2npz" "npy2mat" "npz2mat" DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

add_executable(example1 example1.cpp)
target_link_libraries(example1 cnpy)
