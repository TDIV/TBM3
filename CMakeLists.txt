cmake_minimum_required(VERSION 2.8)

#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
#SET(MACOSX_DEPLOYMENT_TARGET ${DARWIN_MAJOR_SDK_VERSION}.${DARWIN_MINOR_SDK_VERSION})
#SET(CMAKE_CXX_FLAGS "--sysroot ${CMAKE_OSX_SYSROOT} ${CMAKE_CXX_FLAGS}")
MESSAGE("Setting MACOSX_DEPLOYMENT_TARGET to '${MACOSX_DEPLOYMENT_TARGET}'.")

project(TBM3)


add_definitions(-DADD_ -DHAVE_CUBLAS)

# --- ACCELERATE framework (for lapack interface) -----
if(APPLE)
   find_library(ACCELERATE_LIBRARY accelerate)
   mark_as_advanced (ACCELERATE_LIBRARY )
   set(EXTRA_LIBS ${ACCELERATE_LIBRARY} )
endif (APPLE)



if (NOT DEFINED     BOOST_ROOT        AND
    NOT DEFINED ENV{BOOST_ROOT}       AND
    NOT DEFINED     BOOST_INCLUDEDIR  AND
    NOT DEFINED ENV{BOOST_INCLUDEDIR} AND
    NOT DEFINED     BOOST_LIBRARYDIR  AND
    NOT DEFINED ENV{BOOST_LIBRARYDIR})
    if (APPLE)
        set(BOOST_ROOT "../../../boost/boost_1_54_0/mac")
    elseif (WIN32)
        set(BOOST_INCLUDEDIR "C:/local/boost_1_55_0")
        set(BOOST_LIBRARYDIR "C:/local/boost_1_55_0/lib32-msvc-10.0")
    endif()
endif()
if (APPLE OR WIN32)
    set(Boost_USE_STATIC_LIBS TRUE)
endif()
find_package(Boost 1.32 REQUIRED COMPONENTS program_options)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# --- MAGMA -----------------------------------------
set(LIBS ${LIBS} magma )
link_directories(/Users/anrris/opt/magma/lib)
include_directories(SYSTEM /Users/anrris/opt/magma/include)

# --- CUDA -----------------------------------------
find_package(CUDA 6.5)
if (CUDA_FOUND)
  message(STATUS "Found CUDA " ${CUDA_VERSION})
  include_directories(SYSTEM ${CUDA_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY})
  add_definitions(-DWITH_CUDA)
  if (APPLE)
    # With recent versions of Xcode, must explicitly use clang
    set(CUDA_HOST_COMPILER /usr/bin/clang)
  endif()
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O2;
    -gencode=arch=compute_20,code=sm_20;
    -gencode=arch=compute_30,code=sm_30;
    -gencode=arch=compute_35,code=sm_35;
    -gencode=arch=compute_50,code=sm_50)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
endif()


# --- Execute -----------------------------------------
include_directories(SYSTEM ~/opt/TBM3/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

add_executable(tbm-run src/main-tbm-run.cpp)
add_executable(tbm-wannier src/main-tbm-wannier.cpp)

target_link_libraries(tbm-run ${LIBS} ${EXTRA_LIBS} stdc++)
target_link_libraries(tbm-wannier ${LIBS} ${EXTRA_LIBS} stdc++)


# -- REPORT -------------------------------------------
message("")
message("Configuration results")
message("---------------------")
message("C++ compiler        : ${CMAKE_CXX_COMPILER}")
message("C++ compiler flags  : ${CMAKE_CXX_FLAGS}")

get_property(INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("LIBS:     ${LIBS}")
message("INCLUDES: ${INC_DIRS}")
message("WITH_MPI: ${WITH_MPI}")
message("")
