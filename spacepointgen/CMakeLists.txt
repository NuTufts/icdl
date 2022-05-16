cmake_minimum_required(VERSION 3.6)

project(spacepointgen CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

set( larlite_DIR $ENV{LARLITE_LIBDIR}/cmake/larlite )
set( larcv_DIR $ENV{LARCV_LIBDIR}/cmake/larcv )
find_package( ROOT REQUIRED )
find_package( larlite REQUIRED )
find_package( larcv REQUIRED )

#add_executable
add_executable(gen_wireoverlap_tensors gen_wireoverlap_tensors.cxx)
add_executable(gen_spacepoints_dev gen_spacepoints_dev.cxx)

# To set the C++ standard to be used, only one of the following two lines should
# be un-commented.
#set_property(TARGET demo PROPERTY CXX_STANDARD 14)
set_property(TARGET gen_wireoverlap_tensors PROPERTY CXX_STANDARD 17)
set_property(TARGET gen_spacepoints_dev PROPERTY CXX_STANDARD 17)

include_directories( ${LARLITE_INCLUDE_DIR} )
include_directories( ${LARCV_INCLUDE_DIR} )
include_directories( ${ROOT_INCLUDE_DIRS} )

set( LINK_LIBS larlite_LArUtil larlite_Base LArCVCoreDataFormat LArCVCoreBase ${ROOT_LIBRARIES} )

target_link_libraries(gen_wireoverlap_tensors ${LINK_LIBS})
target_link_libraries(gen_spacepoints_dev ${LINK_LIBS})

