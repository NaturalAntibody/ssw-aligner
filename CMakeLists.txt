cmake_minimum_required(VERSION 3.14)
project(ssw_aligner CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Options
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(HAVE_MPI "Enable MPI support" OFF)
option(HAVE_AVX2 "Enable AVX2 SIMD" OFF)
option(HAVE_SSE4_1 "Enable SSE4.1 SIMD" OFF)
option(HAVE_SSE2 "Enable SSE2 SIMD" ON)

# Compiler flags
add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare)
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR NOT CMAKE_BUILD_TYPE)
    add_compile_options(-O3 -DNDEBUG)
endif()

# SIMD flags
if(HAVE_AVX2)
    add_compile_options(-mavx2 -DAVX2)
    add_compile_definitions(AVX2=1)
elseif(HAVE_SSE4_1)
    add_compile_options(-msse4.1)
    add_compile_definitions(SSE=1)
elseif(HAVE_SSE2)
    # SSE2 is default on x86_64 - no flag needed
    add_compile_definitions(SSE=1)
endif()

# Use fmt as header-only
add_compile_definitions(FMT_HEADER_ONLY=1)

# ---------- Include directories ----------
set(INC_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/api
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alignment
    ${CMAKE_CURRENT_SOURCE_DIR}/src/commons
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/alp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/simd
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/simde
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/fmt
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/block-aligner/c
    ${CMAKE_CURRENT_SOURCE_DIR}/generated
)

# ---------- Source files ----------
# Core alignment
set(ALIGNMENT_SRCS
    src/alignment/StripedSmithWaterman.cpp
)

# Commons
set(COMMONS_SRCS
    src/commons/BaseMatrix.cpp
    src/commons/Debug.cpp
    src/commons/LambdaCalculation.cpp
    src/commons/LibraryReader.cpp
    src/commons/MultiParam.cpp
    src/commons/Parameters.cpp
    src/commons/ProfileStates.cpp
    src/commons/Sequence.cpp
    src/commons/SubstitutionMatrix.cpp
    src/commons/Util.cpp
)

# ALP library (Gumbel statistics)
set(ALP_SRCS
    lib/alp/njn_dynprogprob.cpp
    lib/alp/njn_dynprogproblim.cpp
    lib/alp/njn_dynprogprobproto.cpp
    lib/alp/njn_ioutil.cpp
    lib/alp/njn_localmaxstat.cpp
    lib/alp/njn_localmaxstatmatrix.cpp
    lib/alp/njn_localmaxstatutil.cpp
    lib/alp/njn_random.cpp
    lib/alp/sls_alignment_evaluer.cpp
    lib/alp/sls_alp.cpp
    lib/alp/sls_alp_data.cpp
    lib/alp/sls_alp_regression.cpp
    lib/alp/sls_alp_sim.cpp
    lib/alp/sls_basic.cpp
    lib/alp/sls_pvalues.cpp
)

# C API adapter
set(API_SRCS
    src/api/ssw_api.cpp
)

# Block-aligner stub (C)
set(BLOCK_ALIGNER_SRCS
    lib/block-aligner/c/block_aligner_stub.c
)

# All sources
set(ALL_SRCS
    ${ALIGNMENT_SRCS}
    ${COMMONS_SRCS}
    ${API_SRCS}
    ${ALP_SRCS}
    ${BLOCK_ALIGNER_SRCS}
)

# Exclude MPI file if not using MPI
if(NOT HAVE_MPI)
    list(REMOVE_ITEM ALL_SRCS src/commons/MMseqsMPI.cpp)
endif()

# ---------- Library target ----------
add_library(ssw_aligner ${ALL_SRCS})

target_include_directories(ssw_aligner PUBLIC ${INC_DIRS})

# Link system libraries
target_link_libraries(ssw_aligner PUBLIC pthread m)

if(HAVE_MPI)
    find_package(MPI REQUIRED)
    target_compile_definitions(ssw_aligner PUBLIC HAVE_MPI=1)
    target_link_libraries(ssw_aligner PUBLIC MPI::MPI_CXX)
endif()

# Provide binary_name and version symbols required by Parameters.cpp
target_sources(ssw_aligner PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
    "const char* binary_name = \"ssw_aligner\";\nconst char* version = \"1.0.0\";\n")

# ---------- Install ----------
install(TARGETS ssw_aligner
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES
    src/api/ssw_api.h
    src/alignment/StripedSmithWaterman.h
    src/alignment/EvalueComputation.h
    src/alignment/DistanceCalculator.h
    DESTINATION include/ssw_aligner/alignment
)

install(FILES
    src/commons/Sequence.h
    src/commons/Parameters.h
    src/commons/BaseMatrix.h
    src/commons/SubstitutionMatrix.h
    src/commons/Debug.h
    src/commons/Util.h
    src/commons/MathUtil.h
    src/commons/ScoreMatrix.h
    DESTINATION include/ssw_aligner/commons
)
