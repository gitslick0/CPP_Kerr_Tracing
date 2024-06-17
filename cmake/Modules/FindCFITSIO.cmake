# FindCFITSIO.cmake
# Locate CFITSIO library
# This module defines
#  CFITSIO_FOUND, if false, do not try to use CFITSIO
#  CFITSIO_INCLUDE_DIR, where to find fitsio.h
#  CFITSIO_LIBRARIES, the libraries to link against to use CFITSIO

find_path(CFITSIO_INCLUDE_DIR fitsio.h
  /usr/include
  /usr/local/include
)

find_library(CFITSIO_LIBRARIES NAMES cfitsio
  /usr/lib
  /usr/local/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFITSIO DEFAULT_MSG
                                  CFITSIO_INCLUDE_DIR CFITSIO_LIBRARIES)

mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARIES)
