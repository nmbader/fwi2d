# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)
find_path (FFTWF_INCLUDES fftw3.f03)

find_library (FFTW_LIB NAMES fftw3)
find_library (FFTWF_LIB NAMES fftw3f)
if(FFTWF_INCLUDES)
if(FFTW_LIB)
set(FFTW_LIBRARIES ${FFTW_LIB})
endif(FFTW_LIB)
if(FFTWF_LIB)
set(FFTW_LIBRARIES ${FFTWF_LIB} ${FFTW_LIB})
endif(FFTWF_LIB)
endif(FFTWF_INCLUDES)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES )
