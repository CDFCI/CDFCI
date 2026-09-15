cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

project(
    CDFCI
    VERSION 1.0.2
    DESCRIPTION "Coordinate Descent Full Configuration Interaction"
    LANGUAGES CXX
)

include(CTest)
include(CheckIncludeFileCXX)
include(FetchContent)
include(GNUInstallDirs)

option(CDFCI_ENABLE_OPENMP "Build OpenMP-enabled executables when OpenMP is available" ON)
option(CDFCI_BUILD_PYTHON "Build the optional Python extension module" OFF)
set(CDFCI_ARCH_FLAGS "" CACHE STRING "Additional architecture-specific compiler flags")

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Prefer an existing Eigen installation.  A pinned release is downloaded when
# Eigen is not installed, so a submodule checkout is not required.
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/eigen/Eigen/Core")
        message(STATUS "Using the vendored Eigen checkout")
        set(eigen_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/eigen")
    else()
        message(STATUS "Eigen3 was not found; fetching Eigen 3.4.0")
        FetchContent_Declare(
            eigen
            URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
            URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
            # Eigen is header-only here.  Skipping its top-level CMake project
            # avoids importing Eigen's tests and install rules into this project.
            SOURCE_SUBDIR _cdfci_header_only
        )
        FetchContent_MakeAvailable(eigen)
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
    set_target_properties(
        Eigen3::Eigen
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}"
    )
endif()

add_library(cdfci_headers INTERFACE)
add_library(CDFCI::headers ALIAS cdfci_headers)
target_compile_features(cdfci_headers INTERFACE cxx_std_17)
target_include_directories(
    cdfci_headers
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(cdfci_headers INTERFACE Eigen3::Eigen)

# config.h selects __float128 for GCC when quadmath.h is available.  GCC ships
# the matching runtime, which must be linked by every executable/module.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    check_include_file_cxx(quadmath.h CDFCI_HAVE_QUADMATH_HEADER)
    if(CDFCI_HAVE_QUADMATH_HEADER)
        target_link_libraries(cdfci_headers INTERFACE quadmath)
    endif()
endif()

if(CDFCI_ARCH_FLAGS)
    separate_arguments(CDFCI_ARCH_FLAGS_LIST NATIVE_COMMAND "${CDFCI_ARCH_FLAGS}")
    target_compile_options(cdfci_headers INTERFACE ${CDFCI_ARCH_FLAGS_LIST})
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(
        cdfci_headers
        INTERFACE
            -Wno-deprecated-builtins
            -Wno-tautological-constant-compare
    )
endif()

set(CDFCI_HAVE_OPENMP FALSE)
if(CDFCI_ENABLE_OPENMP)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        set(CDFCI_HAVE_OPENMP TRUE)
    else()
        message(WARNING "OpenMP was not found; OpenMP targets will not be built")
    endif()
endif()

function(cdfci_add_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE CDFCI::headers)
endfunction()

function(cdfci_add_openmp_executable target source)
    if(CDFCI_HAVE_OPENMP)
        cdfci_add_executable(${target} ${source})
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endfunction()

cdfci_add_executable(cdfci src/main.cpp)
cdfci_add_executable(cdfci_tools src/tools.cpp)
cdfci_add_executable(optorbfci src/optorbfci.cpp)
cdfci_add_executable(xcdfci src/xcdfci.cpp)
cdfci_add_openmp_executable(cdfci_omp src/main.cpp)
cdfci_add_openmp_executable(xcdfci_omp src/xcdfci.cpp)

set(CDFCI_INSTALL_TARGETS cdfci cdfci_tools optorbfci xcdfci)
if(CDFCI_HAVE_OPENMP)
    list(APPEND CDFCI_INSTALL_TARGETS cdfci_omp xcdfci_omp)
endif()
install(TARGETS ${CDFCI_INSTALL_TARGETS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(CDFCI_BUILD_PYTHON)
    find_package(Python 3 REQUIRED COMPONENTS Interpreter Development.Module)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 2.11 QUIET CONFIG)
    if(NOT pybind11_FOUND)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/pybind11/CMakeLists.txt")
            message(STATUS "Using the vendored pybind11 checkout")
            add_subdirectory(src/lib/pybind11 EXCLUDE_FROM_ALL)
        else()
            message(STATUS "pybind11 was not found; fetching pybind11 2.13.6")
            FetchContent_Declare(
                pybind11
                URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.tar.gz
            )
            FetchContent_MakeAvailable(pybind11)
        endif()
    endif()

    pybind11_add_module(cdfci_python MODULE src/pycdfci.cpp)
    set_target_properties(
        cdfci_python
        PROPERTIES
            OUTPUT_NAME _cdfci
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    )
    target_link_libraries(cdfci_python PRIVATE CDFCI::headers)
    if(CDFCI_HAVE_OPENMP)
        target_link_libraries(cdfci_python PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

if(BUILD_TESTING)
    set(
        CDFCI_REGRESSION_TEST_CASES
        c2_ccpvdz_psi4
        cr2_ahlrichs_psi4
        h2o_ccpvdz_psi4
        h2o_ccpvdz_pyscf
        h2o_sto3g_psi4
        hubbard
        n2_ccpvdz_psi4
        n2_ccpvdz_psi4_eps1e-2
        n2_ccpvdz_psi4_triplet
    )

    function(cdfci_add_regression_tests target variant)
        foreach(test_case IN LISTS CDFCI_REGRESSION_TEST_CASES)
            set(test_name "cdfci.${variant}.${test_case}")
            add_test(
                NAME "${test_name}"
                COMMAND "$<TARGET_FILE:${target}>" -D "--test-case=${test_case}"
            )
            set_tests_properties(
                "${test_name}"
                PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test"
            )
            if(ARGC GREATER 2)
                set_tests_properties("${test_name}" PROPERTIES ENVIRONMENT "${ARGV2}")
            endif()
        endforeach()
    endfunction()

    cdfci_add_executable(cdfci_test test/test_system.cpp)
    cdfci_add_regression_tests(cdfci_test serial)

    if(CDFCI_HAVE_OPENMP)
        cdfci_add_openmp_executable(cdfci_test_omp test/test_system.cpp)
        cdfci_add_regression_tests(cdfci_test_omp openmp "OMP_NUM_THREADS=2")
    endif()
endif()

if(UNIX AND CDFCI_HAVE_OPENMP)
    find_program(BASH_EXECUTABLE bash REQUIRED)
    add_custom_target(
        examples
        COMMAND
            "${CMAKE_COMMAND}" -E env
            "CDFCI_BIN_DIR=$<TARGET_FILE_DIR:cdfci>"
            "${BASH_EXECUTABLE}" ./run_all.sh
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples"
        DEPENDS cdfci cdfci_omp cdfci_tools optorbfci xcdfci
        USES_TERMINAL
        COMMENT "Running CDFCI examples"
    )
endif()
