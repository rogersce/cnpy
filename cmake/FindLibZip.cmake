find_package(PkgConfig)
pkg_check_modules(PC_LIBZIP QUIET libzip)
set(LIBZIP_DEFINITIONS ${PC_LIBZIP_CFLAGS_OTHER})

find_path(LIBZIP_INCLUDE_DIR zip.h
          HINTS ${PC_LIBZIP_INCLUDEDIR} ${PC_LIBZIPP_INCLUDE_DIRS}
          PATH_SUFFIXES zip )

find_library(LIBZIP_LIBRARY NAMES zip libzip
             HINTS ${PC_LIBZIP_LIBDIR} ${PC_LIBZIP_LIBRARY_DIRS} )

set(LIBZIP_LIBRARIES ${LIBZIP_LIBRARY} )
set(LIBZIP_INCLUDE_DIRS ${LIBZIP_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libzip  DEFAULT_MSG LIBZIP_LIBRARY LIBZIP_INCLUDE_DIR)

mark_as_advanced(LIBZIP_INCLUDE_DIR LIBZIP_LIBRARY )
