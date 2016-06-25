cmake_minimum_required(VERSION 2.8)

project(GraMat)


add_definitions(-DADD_ -DHAVE_CUBLAS)


# --- ACCELERATE framework -------------------------
if(APPLE)
   find_library(ACCELERATE_LIBRARY accelerate)
   mark_as_advanced (ACCELERATE_LIBRARY )
   set(EXTRA_LIBS ${ACCELERATE_LIBRARY})
endif (APPLE)



# --- MAGMA -----------------------------------------
set(LIBS ${LIBS} magma)
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
include_directories(SYSTEM /Users/anrris/Dropbox/code-gitsource/Gramat/include)

add_executable(run main.cpp)
target_link_libraries(run ${LIBS} ${EXTRA_LIBS})







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
