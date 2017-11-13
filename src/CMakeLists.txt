# Define source files
set( KSVD_SRC pdgeqsvd.c )


list( APPEND CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-new-dtags")
# Transform the list into a space separate string
string(REPLACE ";" " " CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")

# Build KSVD
add_library( ${KSVD_LIBNAME} ${KSVD_SRC} )
# Libraries to link
target_link_libraries( ${KSVD_LIBNAME} ${MPI_C_LIBRARIES} ${SCALAPACK_LIBRARIES} ${ELPA_LIBRARIES} ${EXTRA_LIBS} )

# Add ELPA depencency, to build ELPA first
if( NOT ELPA_FOUND )
    add_dependencies(${KSVD_LIBNAME} ELPA)
endif()

# Installation
install( TARGETS ${KSVD_LIBNAME} DESTINATION lib )
