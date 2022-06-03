cmake_minimum_required(VERSION 3.10)
project(litemaker VERSION 1.0.0)

#---Install location by default is in the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/installed" CACHE PATH "${CMAKE_BUILD_DIR}" FORCE )
endif()

#---Module directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "$ENV{IFDHC_LIB}/ifdhc/cmake")
message( "${CMAKE_MODULE_PATH}" )

## Dependencies

#---CERN ROOT
find_package(ROOT CONFIG REQUIRED)
# we set the compiler version based on the one used for ROOT
execute_process(
    COMMAND bash "-c" "root-config --features | awk '{print $1}'"
    OUTPUT_VARIABLE ROOT_CXX_STANDARD )
if (${ROOT_CXX_STANDARD} MATCHES "cxx17")
  message( STATUS "using cxx17")
  set(CMAKE_CXX_STANDARD 17)
elseif (${ROOT_CXX_STANDARD} MATCHES "cxx14")
  message( STATUS "ssing cxx14")
  set(CMAKE_CXX_STANDARD 14)
else ()
  message( STATUS "using cxx11")  
  set(CMAKE_CXX_STANDARD 11)
endif()

#---PYTHON

if (USE_PYTHON3)
  find_package( PythonInterp 3 REQUIRED )
  find_package( PythonLibs 3 REQUIRED )
  # warning gross hack
  execute_process( COMMAND python3 -c "import numpy; print(numpy.get_include())" OUTPUT_VARIABLE NUMPY_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE )
  #message(${NUMPY_INCLUDES})
  set(LITEMAKER_USE_PYTHON 1)
  set(LITEMAKER_PYTHON_VERSION 3)
elseif(USE_PYTHON2)
  find_package( PythonInterp 2.7 REQUIRED )
  find_package( PythonLibs 2.7 REQUIRED )
  find_package( Numpy REQUIRED )
  set(LITEMAKER_USE_PYTHON 1)
  set(LITEMAKER_PYTHON_VERSION 2)
else()
  set(LITEMAKER_USE_PYTHON 0)
endif()
message( STATUS "Likemaker python version: ${LITEMAKER_PYTHON_VERSION}")

# -- LARLITE --
set( larlite_DIR $ENV{LARLITE_LIBDIR}/cmake/larlite )
find_package( larlite REQUIRED )

# -- ART/GALLERY/CANVAS --
find_package(PkgConfig)
pkg_check_modules(SQLite3 REQUIRED sqlite3)
set( SQLite3_LIBRARY $ENV{SQLITE_LIB}/libsqlite3_ups.so )
find_package( Boost REQUIRED )
find_package( Range )
find_package( cetlib REQUIRED QUIET )
find_package( cetlib_except REQUIRED QUIET )
find_package( messagefacility REQUIRED QUIET )
find_package( art REQUIRED QUIET )
find_package( hep_concurrency REQUIRED QUIET )
find_package( canvas REQUIRED QUIET )
find_package( gallery REQUIRED QUIET )
#find_package( icarusalg REQUIRED QUIET )
find_package( lardataobj REQUIRED QUIET )
find_package( lardata REQUIRED QUIET )
find_package( lardataalg REQUIRED QUIET )
find_package( larcoreobj REQUIRED QUIET )
find_package( larcorealg REQUIRED QUIET )
find_package( larcore REQUIRED QUIET )
find_package( ifdhc REQUIRED QUIET )
find_package( nusimdata REQUIRED QUIET )
find_package( sbnobj REQUIRED QUIET )
#find_package( larsim REQUIRED QUIET )


## Setup

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Define Variables for subdirectory cmakelists
set(LITEMAKER_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

## Add Source

#---Add base of source
add_subdirectory(litemaker)

#---Configuration export
include(CMakePackageConfigHelpers)
set(LITEMAKER_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(LITEMAKER_LIBRARY_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(TARGETS_EXPORT_NAME "litemakerTargets")
configure_package_config_file("cmake/Config.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake/litemakerConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/litemaker"
  PATH_VARS LITEMAKER_INCLUDE_DIR LITEMAKER_LIBRARY_DIR LITEMAKER_BASE_DIR LITEMAKER_USE_PYTHON LITEMAKER_PYTHON_VERSION)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/cmake/litemakerConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)


## Export Targets

#---Create info to include Litemaker in other cmake files
install(EXPORT litemaker
        FILE litemakerTargets.cmake
        DESTINATION "${CMAKE_BINARY_DIR}/installed/lib/cmake/litemaker")

#---Create info to include Litemaker in other cmake files
install(FILES "${CMAKE_BINARY_DIR}/cmake/litemakerConfig.cmake"
  "${CMAKE_BINARY_DIR}/cmake/litemakerConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/litemaker")
