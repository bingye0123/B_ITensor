
cmake_minimum_required(VERSION 2.8)

project(it_dmrg)

set(CMAKE_VERBOSE_MAKEFILE ON)

# set(BOOST_DIR "~/opt/boost_1_54_0")
# set(BOOST_DIR "~/code/boost")
set(BOOST_DIR "~/bingye/code/boost_1_55_0")

set(CMAKE_CXX_FLAGS "-O2 -DNDEBUG -m64")

set(BUILD_SHARED_LIBS OFF)

set(CMAKE_VERBOSE_MAKEFILE ON)

# set(BOOST_DIR "~/opt/boost_1_54_0")
set(BOOST_DIR "~/bingye/code/boost_1_55_0")

add_definitions(-DPLATFORM_mkl)
# set(BLAS_LAPACK_INCLUDE_DIR "/home/rany/opt/intel/composer_xe_2013_sp1.0.080/composer_xe_2013_sp1.0.080/mkl/include")
# set(BLAS_LAPACK_LIB_DIR "/home/rany/opt/intel/composer_xe_2013_sp1.0.080/composer_xe_2013_sp1.0.080/compiler/lib/intel64;/home/rany/opt/intel/composer_xe_2013_sp1.0.080/composer_xe_2013_sp1.0.080/mkl/lib/intel64")
# set(BLAS_LAPACK_INCLUDE_DIR "/usr/public/intel/composer_xe_2013.1.117/mkl/include")
# set(BLAS_LAPACK_LIB_DIR "/usr/public/intel/composer_xe_2013.1.117/compiler/lib/intel64;/usr/public/intel/composer_xe_2013.1.117/mkl/lib/intel64")
set(BLAS_LAPACK_INCLUDE_DIR "/home/yeba/intel/composer_xe_2013_sp1.0.080/mkl/include")
set(BLAS_LAPACK_LIB_DIR
"/home/yeba/intel/composer_xe_2013_sp1.0.080/compiler/lib/intel64;/home/yeba/intel/composer_xe_2013_sp1.0.080/mkl/lib/intel64")

# Require MPI for this project:
find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})


# set(TRNG_INCLUDE_DIR "~/opt/trng/include")
# set(TRNG_LIB_DIR "~/opt/trng/lib")
set(TRNG_INCLUDE_DIR "~/bingye/code/trng-4.14/include")
set(TRNG_LIB_DIR "~/bingye/code/trng-4.14/lib")

# set(ITPP_INCLUDE_DIR "~/opt/itpp/include")
# set(ITPP_LIB_DIR "~/opt/itpp/lib")
# set(ITPP_INCLUDE_DIR "/usr/public/itpp/4.3.1/include")
# set(ITPP_LIB_DIR "/usr/public/itpp/4.3.1/lib")
# set(ITPP_INCLUDE_DIR "~/code/itpp4.3.1intel/include")
# set(ITPP_LIB_DIR "~/code/itpp4.3.1intel/lib")
set(ITPP_INCLUDE_DIR "~/bingye/code/itpp4.2intel/include")
set(ITPP_LIB_DIR "~/bingye/code/itpp4.2intel/lib")


#--------------------------------------------------
set(ITENSOR_INCLUDE_DIR "/home/yeba/bingye/code/itensor/include")
set(ITENSOR_LIB_DIR "/home/yeba/bingye/code/itensor/lib")
#--------------------------------------------------


include_directories(sym_MonteCarlo utilities matrix itensor ${BOOST_DIR} ${BLAS_LAPACK_INCLUDE_DIR} ${ITPP_INCLUDE_DIR} ${TRNG_INCLUDE_DIR} ${ITENSOR_INCLUDE_DIR})
#include_directories(include)

add_subdirectory(utilities)
add_subdirectory(matrix)
add_subdirectory(itensor)
add_subdirectory(sym_MonteCarlo)

link_directories(sym_MonteCarlo utilities matrix itensor ${BLAS_LAPACK_LIB_DIR} ${ITPP_LIB_DIR} ${TRNG_LIB_DIR} ${ITENSOR_LIB_DIR})

add_executable(it_dmrg it_dmrg.cc)
target_link_libraries(it_dmrg symMonteCarlo -litpp itensor matrix utilities -Wl,--start-group mkl_intel_lp64.a mkl_intel_thread.a mkl_core.a --end-group iomp5 gfortran pthread trng4.a ${MPI_LIBRARIES})


