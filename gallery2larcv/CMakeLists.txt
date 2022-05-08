cmake_minimum_required(VERSION 3.6)
project(gallery2larcv CXX)


# dependency hell
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} )
find_package(PkgConfig)
pkg_check_modules(SQLite3 REQUIRED sqlite3)
set( SQLite3_LIBRARY $ENV{SQLITE_LIB}/libsqlite3_ups.so )
find_package( Boost REQUIRED )
find_package( ROOT REQUIRED )
set( Range_DIR $CMAKE_CURRENT_SOURCE_DIR )
find_package( Range )
find_package( icarusalg REQUIRED QUIET )
set( larlite_DIR $ENV{LARLITE_LIBDIR}/cmake/larlite )
find_package( larlite REQUIRED )

add_executable(gallery2larcv gallery2larcv.cxx)

# To set the C++ standard to be used, only one of the following two lines should
# be un-commented.
#set_property(TARGET demo PROPERTY CXX_STANDARD 14)
set_property(TARGET gallery2larcv PROPERTY CXX_STANDARD 17)

include_directories( ${fhiclcpp_INCLUDE_DIRS} )
include_directories( ${messagefacility_INCLUDE_DIRS} )
include_directories( ${lardataobj_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${icarusalg_INCLUDE_DIRS} )
include_directories( ${LARLITE_INCLUDE_DIR} )

message( status "${larlite_INCLUDE_DIRS}" )

set( LINK_LIBS fhiclcpp::fhiclcpp larcorealg::Geometry icarusalg::Geometry larlite_Geometry ${ROOT_LIBRARIES} )

#target_link_libraries(demo ${LINK_LIBS} )
target_link_libraries(gallery2larcv ${LINK_LIBS})

