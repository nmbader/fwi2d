#
# - SymLink
#
# defines the macro install_symlink due to Rian Quinn
#
set (SYMLINK_FOUND FALSE)
if (CMAKE_HOST_UNIX)
set (SYMLINK_FOUND TRUE)
macro(install_symlink filepath sympath)
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${filepath} ${sympath})")
    install(CODE "message(\"-- Created symlink: ${sympath} -> ${filepath}\")")
endmacro(install_symlink)
    MESSAGE(STATUS "Symlinks are available")
endif (CMAKE_HOST_UNIX)
