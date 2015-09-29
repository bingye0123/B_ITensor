cmake_minimum_required(VERSION 2.8)

project(it_dmrg_sandbox)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(BOOST_DIR "~/opt/boost_1_54_0")

add_definitions(-DPLATFORM_mkl)
set(BLAS_LAPACK_INCLUDE_DIR "/opt/intel/mkl/include")
set(BLAS_LAPACK_LIB_DIR "/opt/intel/lib/intel64;/opt/intel/mkl/lib/intel64")

set(IT_DMRG_INCLUDE_DIR "/home/ranying/projects/IT_DMRG/utilities;/home/ranying/projects/IT_DMRG/matrix;/home/ranying/projects/IT_DMRG/itensor")

set(IT_DMRG_LIB_DIR "/home/ranying/projects/IT_DMRG/build/utilities;/home/ranying/projects/IT_DMRG/build/matrix;/home/ranying/projects/IT_DMRG/build/itensor")

include_directories(${IT_DMRG_INCLUDE_DIR} ${BOOST_DIR} ${BLAS_LAPACK_INCLUDE_DIR})

link_directories(${IT_DMRG_LIB_DIR} ${BLAS_LAPACK_LIB_DIR})

add_executable(dmrg dmrg.cc)
target_link_libraries(dmrg itensor matrix utilities mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 gfortran pthread)

add_executable(dmrgj1j2 dmrgj1j2.cc)
target_link_libraries(dmrgj1j2 itensor matrix utilities mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 gfortran pthread)

add_executable(dmrg_tutorial dmrg_tutorial.cc)
target_link_libraries(dmrg_tutorial itensor matrix utilities mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 gfortran pthread)

add_executable(iqdmrg iqdmrg.cc)
target_link_libraries(iqdmrg itensor matrix utilities mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 gfortran pthread)


