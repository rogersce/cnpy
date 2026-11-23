find_package(doctest REQUIRED)

add_executable(cnpy_tests
    main.cpp
    test_dtype.cpp
    test_npy.cpp
    test_mmap.cpp
    test_npz.cpp)
target_link_libraries(cnpy_tests PRIVATE cnpy::cnpy doctest::doctest cnpy_warnings)
add_test(NAME unit COMMAND cnpy_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# NumPy interop: Python writes files that C++ checks, C++ writes files that Python checks
# (including byte-for-byte equality with numpy.save).
add_executable(cnpy_interop interop/interop.cpp)
target_link_libraries(cnpy_interop PRIVATE cnpy::cnpy cnpy_warnings)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE _cnpy_numpy_missing OUTPUT_QUIET ERROR_QUIET)
endif()
if(Python3_Interpreter_FOUND AND NOT _cnpy_numpy_missing)
    add_test(NAME numpy_interop
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/interop/interop.py
                --exe $<TARGET_FILE:cnpy_interop>
                --workdir ${CMAKE_CURRENT_BINARY_DIR}/interop_data)
else()
    message(STATUS "cnpy: Python 3 with NumPy not found, skipping the NumPy interop test")
endif()
