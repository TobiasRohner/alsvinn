cmake_minimum_required (VERSION 3.1.0)



FILE(GLOB_RECURSE SRC src/*.cpp)
FILE(GLOB_RECURSE HEADERS include/*.hpp)

set_source_groups( "Source Files" "${CMAKE_CURRENT_SOURCE_DIR}/src" "${SRC}")
set_source_groups( "Header Files" "${CMAKE_CURRENT_SOURCE_DIR}/include/alsuq" "${HEADERS}")


INCLUDE_DIRECTORIES("include" ${CUDA_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

IF (NOT ALSVINN_USE_MPI)
    FILE(GLOB_RECURSE MPI_SRC src/*MPI*.cpp src/*mpi*.cpp)
    LIST(REMOVE_ITEM SRC ${MPI_SRC})
ENDIF()

FILE(GLOB_RECURSE HEADERS include/*.hpp)

IF (NOT ALSVINN_USE_MPI)
    FILE(GLOB_RECURSE MPI_HEADERS include/*MPI*.hpp include/*mpi*.hpp)
    LIST(REMOVE_ITEM HEADERS ${MPI_HEADERS})
ENDIF()

IF(ALSVINN_HAVE_CUDA)
	FILE(GLOB_RECURSE CUDA_SRC src/*.cu)
	set_source_groups( "CUDA Files" "${CMAKE_CURRENT_SOURCE_DIR}/src/" "${CUDA_SRC}")
  ADD_LIBRARY(alsuq_cuda SHARED ${CUDA_SRC})
  target_link_libraries(alsuq_cuda PUBLIC alsutils_include alsfvm_include MPI::MPI_C_include Boost::log)
ELSE()
	ADD_LIBRARY(alsuq_cuda INTERFACE)
ENDIF()

add_library(alsvinn::alsuq_cuda ALIAS alsuq_cuda)
ADD_LIBRARY( alsuq SHARED ${SRC} ${HEADERS})
add_library(alsvinn::alsuq ALIAS alsuq)

target_include_directories(alsuq PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

TARGET_LINK_LIBRARIES(alsuq
  PUBLIC
  alsfvm
  alsutils
  alsuq_cuda
  QMC::qmc_base
  Boost::filesystem
  Boost::system
  Boost::log
  Boost::date_time	
  MPI::MPI_C
  hdf5::HDF5_C
  Python::Python
  netCDF::netcdf
  ${CUDA_LIBRARIES})


install(TARGETS alsuq_cuda
    EXPORT alsvinn-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT alsuq)

install(TARGETS alsuq
    EXPORT alsvinn-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT alsuq)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT alsuq)
