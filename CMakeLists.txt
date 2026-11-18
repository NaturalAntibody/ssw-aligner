cmake_minimum_required(VERSION 3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(Cython MODULE REQUIRED VERSION 3.0)
include(UseCython)

cython_transpile(ssw_aligner/_ssw_wrapper.pyx LANGUAGE C OUTPUT_VARIABLE ssw_wrapper_c)

python_add_library(_ssw_wrapper MODULE WITH_SOABI
    ${ssw_wrapper_c}
    ssw_aligner/_lib/ssw.c
)
# simde-sse2.h lives in the repo root; it maps SSE2 intrinsics to NEON etc. on non-x86 targets.
# The wrapper includes "_lib/ssw.h" relative to the package directory.
target_include_directories(_ssw_wrapper PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ssw_aligner
)
target_link_libraries(_ssw_wrapper PRIVATE Python::NumPy)

# SIMD hints, as in the scikit-bio build this extension was extracted from.
if(MSVC)
    target_compile_options(_ssw_wrapper PRIVATE /openmp:experimental)
elseif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_compile_options(_ssw_wrapper PRIVATE -fopenmp-simd)
    target_compile_definitions(_ssw_wrapper PRIVATE SIMDE_ENABLE_OPENMP)
elseif(CMAKE_C_COMPILER_ID MATCHES "^Intel")
    target_compile_options(_ssw_wrapper PRIVATE -qopenmp-simd)
    target_compile_definitions(_ssw_wrapper PRIVATE SIMDE_ENABLE_OPENMP)
endif()
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "i686")
    target_compile_options(_ssw_wrapper PRIVATE -msse2)
endif()

install(TARGETS _ssw_wrapper DESTINATION ssw_aligner)
