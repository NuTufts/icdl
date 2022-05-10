cmake_minimum_required(VERSION 3.6)
project(gallerydemo CXX)

# include_directories($ENV{BOOST_INC})
# include_directories($ENV{CANVAS_INC})
# include_directories($ENV{CETLIB_EXCEPT_INC})
# include_directories($ENV{CETLIB_INC})
# include_directories($ENV{CLHEP_INC})
# include_directories($ENV{FHICLCPP_INC})
# include_directories($ENV{GALLERY_INC})
# include_directories($ENV{CANVAS_ROOT_IO_INC})
# include_directories($ENV{LARCOREOBJ_INC})
# include_directories($ENV{LARCOREALG_INC})
# include_directories($ENV{LARDATAOBJ_INC})
# include_directories($ENV{LARDATAALG_INC})
# include_directories($ENV{NUSIMDATA_INC})
# include_directories($ENV{ROOT_INC})
# include_directories($ENV{MESSAGEFACILITY_INC})
# include_directories($ENV{ICARUSCODE_INC})
# link_directories($ENV{ROOTSYS}/lib $ENV{CANVAS_LIB} $ENV{CETLIB_EXCEPT_LIB} $ENV{CETLIB_LIB} $ENV{GALLERY_LIB} $ENV{NUSIMDATA_LIB} $ENV{LARCOREOBJ_LIB} $ENV{LARDATAOBJ_LIB})

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


#add_executable(demo demo.cc analyze.cc)
add_executable(servicesdemo servicesdemo.cc)
add_executable(gen_geodata_rootfile gen_geodata_rootfile.cc)

# To set the C++ standard to be used, only one of the following two lines should
# be un-commented.
#set_property(TARGET demo PROPERTY CXX_STANDARD 14)
set_property(TARGET servicesdemo PROPERTY CXX_STANDARD 17)
set_property(TARGET gen_geodata_rootfile PROPERTY CXX_STANDARD 17)

include_directories( ${fhiclcpp_INCLUDE_DIRS} )
include_directories( ${messagefacility_INCLUDE_DIRS} )
include_directories( ${lardataobj_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${larcorealg_INCLUDE_DIRS} )
include_directories( ${icarusalg_INCLUDE_DIRS} )
include_directories( ${LARLITE_INCLUDE_DIR} )

message( status "${larlite_INCLUDE_DIRS}" )

set( LINK_LIBS fhiclcpp::fhiclcpp larcorealg::Geometry icarusalg::Geometry larlite_LArUtil larlite_Base ${ROOT_LIBRARIES} )

#target_link_libraries(demo ${LINK_LIBS} )
target_link_libraries(servicesdemo ${LINK_LIBS} )
target_link_libraries(gen_geodata_rootfile ${LINK_LIBS})

