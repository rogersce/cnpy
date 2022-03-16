CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(CNPY)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

option(ENABLE_STATIC "Build static (.a) library" ON)

find_package(ZLIB REQUIRED)

include_directories(${ZLIB_INCLUDE_DIRS})

add_library(cnpy SHARED "cnpy.cpp")
# create alias with exported name for including via find_package or FetchContent
add_library(Cnpy::Cnpy ALIAS cnpy)
target_include_directories(cnpy
    PUBLIC
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(cnpy ${ZLIB_LIBRARIES})

include(GNUInstallDirs)

install(TARGETS cnpy
    EXPORT cnpy-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

#This is required so that the exported target has the name JSONUtils and not jsonutils
set_target_properties(cnpy PROPERTIES EXPORT_NAME Cnpy)

if(ENABLE_STATIC)
    add_library(cnpy-static STATIC "cnpy.cpp")
    set_target_properties(cnpy-static PROPERTIES OUTPUT_NAME "cnpy")
    install(TARGETS "cnpy-static"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
endif(ENABLE_STATIC)

install(FILES "cnpy.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "mat2npz" "npy2mat" "npz2mat"
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

add_executable(example1 example1.cpp)
target_link_libraries(example1 cnpy)

install(EXPORT cnpy-targets
  FILE
    CnpyTargets.cmake
  NAMESPACE
    Cnpy::
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/Cnpy
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/CnpyConfigVersion.cmake
    VERSION 0.1.0
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/CnpyConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CnpyConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Cnpy
)

#Install the config, configversion and custom find modules
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CnpyConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/CnpyConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Cnpy
)

