cmake_minimum_required (VERSION 3.10.0)


SET(ALSVINN_USE_HUNTER OFF CACHE BOOL "Use Hunter for automatic compiling of additional packages")
SET(ALSVINN_USE_HUNTER_FOR_BOOST ON CACHE BOOL "Use hunter to build boost")
SET(HUNTER_ROOT ${CMAKE_BINARY_DIR}/hunter_root)

IF(${ALSVINN_USE_HUNTER})
	include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/cpp-pm/hunter/archive/v0.23.297.tar.gz"
    SHA1 "3319fe6a3b08090df7df98dee75134d68e2ef5a3"
    LOCAL
)

ELSE()
	MACRO(hunter_add_package args)
		# do nothing
	ENDMACRO()
ENDIF()

project (alsvinn VERSION 0.5.4 LANGUAGES C CXX )
SET(ALSVINN_USE_FLOAT OFF CACHE BOOL "Do all computations with 32 bits
floating point numbers, instead of 64 bit double precision")
SET(ALSVINN_USE_CUDA ON CACHE BOOL "Build with CUDA")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)




if(${ALSVINN_USE_CUDA})
        enable_language(CUDA)
	find_package(CUDA REQUIRED)
endif()
INCLUDE("cmake/version.cmake")      


MACRO(set_source_groups  folder_name directory source_files)
	foreach(f ${source_files})
		# Get the path of the file relative to ${DIRECTORY},
		# then alter it (not compulsory)

		file(RELATIVE_PATH SRCGR ${directory} ${f})
		set(SRCGR "${folder_name}/${SRCGR}")

		# Extract the folder, ie remove the filename part
		string(REGEX REPLACE "(.*)(/[^/]*)$" "\\1" SRCGR ${SRCGR})

		# Source_group expects \\ (double antislash), not / (slash)
		string(REPLACE / \\ SRCGR ${SRCGR})
		source_group("${SRCGR}" FILES ${f})
	endforeach()
ENDMACRO(set_source_groups)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)
SET(ALSVINN_BUILD_TESTS ON CACHE BOOL "Build unittests (requires Gtest). HIGHLY recommended!")
SET(ALSVINN_USE_MPI ON CACHE BOOL "Build with MPI")
SET(ALSVINN_CXX_FLAGS "" CACHE STRING "CXX flags to use")
SET(ALSVINN_USE_TIMERS ON CACHE BOOL "Enable profiling timings (should not hurt performance)")
SET(ALSVINN_USE_CUDA_TIMERS OFF CACHE BOOL "Enable CUDA profiling timings (hurts performance)")
set(ALSVINN_BUILD_FBM ON CACHE BOOL "Build the fast Fractional Brownian Motion addon for python initial data")
set(ALSVINN_BUILD_QMC ON CACHE BOOL "Build the QMC module for Quasi Monte Carlo")
SET(ALSVINN_IN_DOCKER OFF CACHE BOOL "Only ON if Alsvinn is built within a docker container")
SET(ALSVINN_DOCKER_CONTAINER "Not in docker" CACHE STRING "The docker container used to build docker")
SET(ALSVINN_PRINT_MEMORY_ALLOCATIONS OFF CACHE BOOL "Log every memory allocation")


SET(ALSVINN_BUILD_DOXYGEN ON CACHE BOOL "Build doxygen documentation")
SET(ALSVINN_PROFILE OFF CACHE BOOL "Add profiling information")
SET(ALSVINN_MPI_GPU_DIRECT OFF CACHE BOOL "Use NVidia's GPU direct to do MPI transfers (this is faster, but not only available)")

SET(ALSVINN_PYTHON_VERSION "3.7" CACHE STRING "Python version to use (either 2.7 or 3.6, 3.7, etc)")

STRING(REPLACE "." "" ALSVINN_PYTHON_VERSION_PURE ${ALSVINN_PYTHON_VERSION})


MARK_AS_ADVANCED(ALSVINN_PYTHON_VERSION_PURE)

SET(ALSVINN_BOOST_PYTHON "python${ALSVINN_PYTHON_VERSION_PURE}")
MARK_AS_ADVANCED(ALSVINN_BOOST_PYTHON)

SET(ALSVINN_BOOST_NUMPY "numpy${ALSVINN_PYTHON_VERSION_PURE}")
MARK_AS_ADVANCED(ALSVINN_BOOST_NUMPY)


#ADD_DEFINITIONS(-DMPICH_RDMA_ENABLED_CUDA=1)
FIND_PACKAGE(OpenMP REQUIRED)

add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>)



if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-undefined-var-template>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-gnu-zero-variadic-macro-arguments>)
endif()

if (NOT WIN32)
	add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wall>)
	add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-comment>)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>)
endif()

# Python dependency could be made optional eventually, but for now it's
# the only way of specify initial data
find_package(Python ${ALSVINN_PYTHON_VERSION} REQUIRED COMPONENTS Development)

# We want to specify python versions, which the cmake filees from boost do not allow us to

#set(Boost_NO_BOOST_CMAKE ON)

if(ALSVINN_USE_HUNTER_FOR_BOOST AND ALSVINN_USE_HUNTER)
  hunter_add_package(Boost COMPONENTS regex thread system filesystem chrono date_time program_options log python)
   find_package(Boost CONFIG REQUIRED filesystem thread system chrono date_time program_options log log_setup ${ALSVINN_BOOST_PYTHON} ${ALSVINN_BOOST_NUMPY})

else()
    message(STATUS "Not using hunter to build boost")
    find_package(Boost REQUIRED COMPONENTS filesystem thread system chrono date_time program_options log log_setup ${ALSVINN_BOOST_PYTHON} ${ALSVINN_BOOST_NUMPY})
endif()
    find_package(Boost REQUIRED COMPONENTS filesystem thread system chrono date_time program_options log log_setup ${ALSVINN_BOOST_PYTHON} ${ALSVINN_BOOST_NUMPY} filesystem)
    add_compile_options("-DBOOST_NO_CXX11_SCOPED_ENUMS")
    add_compile_options("-DBOOST_NO_SCOPED_ENUMS")
    add_compile_options("-DBOOST_ALL_DYN_LINK")


include(FetchContent)
find_package(MPI REQUIRED COMPONENTS C)
include_directories(${MPI_C_INCLUDE_PATH})
set(HDF5_PREFER_PARALLEL ON)
set(BUILD_UTILITIES OFF)
set(ENABLE_EXAMPLES OFF)
set(ENABLE_TESTS OFF)
set(ENABLE_DAP OFF)
set(ENABLE_HDF5 ON)
set(HDF5_PARALLEL ON)
set(ENABLE_NETCDF4 ON)
set(ENABLE_PNETCDF ON)
set(ENABLE_PARALLEL4 ON)
set(ENABLE_PARALLEL ON)
set(USE_PARALLEL ON)
FetchContent_Declare(
  NetCDF
  GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
  GIT_TAG main
  OVERRIDE_FIND_PACKAGE
)
find_package(NetCDF REQUIRED)

find_package (Threads)
include("cmake/hdf5.cmake")
#include("cmake/netcdf.cmake")
if (ALSVINN_USE_CUDA)
    include("cmake/cub.cmake")
endif()

if(ALSVINN_USE_CUDA)
	SET(CUDA_TURING_FLAGS "")
	IF (${CMAKE_CUDA_COMPILER_VERSION} VERSION_GREATER_EQUAL 10.0.0)
	   SET(CUDA_TURING_FLAGS "-gencode arch=compute_75,code=sm_75 -gencode arch=compute_70,code=sm_70")
	ENDIF()

        set(CMAKE_CUDA_FLAGS "--use_fast_math -Xptxas -warn-spills -gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=sm_60 ${CUDA_TURING_FLAGS} -gencode arch=compute_61,code=sm_61" )
	set(CMAKE_CUDA_FLAGS_DEBUG "-G -g -O0 -lineinfo -src-in-ptx -Xptxas --device-debug,-O0")
	set(ALSVINN_HAVE_CUDA On)
endif()


if(ALSVINN_USE_MPI)
    include("cmake/mpi.cmake")
endif()

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/alsvinn)

add_subdirectory("addons")
add_subdirectory("alsutils")

add_subdirectory("alsfvm")
add_subdirectory("alsvinncli")

IF(ALSVINN_USE_MPI)
	add_subdirectory("alsuq")
	add_subdirectory("alsuqcli")
ENDIF()
if(ALSVINN_BUILD_TESTS)
  #hunter_add_package(GTest)
	find_package(GTest REQUIRED)
	add_subdirectory("test")


endif()

if(ALSVINN_BUILD_DOXYGEN)
   add_subdirectory("documentation")
endif()
add_subdirectory("python")


install(EXPORT alsvinn-targets
  FILE
    AlsvinnTargets.cmake
  NAMESPACE
    alsvinn::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)


include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/AlsvinnConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/AlsvinnConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/AlsvinnConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)


install(FILES
    ${CMAKE_CURRENT_LIST_DIR}/cmake/netcdf.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/hdf5.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/mpi.cmake
    ${CMAKE_CURRENT_LIST_DIR}/cmake/cub.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AlsvinnConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AlsvinnConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)


export(EXPORT alsvinn-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/AlsvinnTargets.cmake NAMESPACE alsvinn::)
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/netcdf.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/netcdf.cmake
    COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/hdf5.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hdf5.cmake
    COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/cub.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cub.cmake
    COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/mpi.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/mpi.cmake
    COPYONLY)


export(PACKAGE alsvinn)
