# Defines the source files of the tests. Each file generates 1 test
file(GLOB tests_files "*.c")

include_directories( ${CMAKE_SOURCE_DIR}/include )
link_directories( ${CMAKE_BINARY_DIR}/lib )
link_directories( ${POLAR_DIR}/lib )

# Uses RUNPATH instead of RPATH
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-new-dtags")
list( APPEND CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")
# Transform the list into a space separate string
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

foreach( test_src ${tests_files} )
    get_filename_component( test_exe ${test_src}  NAME_WE )
    if( "${test_exe}" STREQUAL "testing" )
        # skip testing.c
        continue()
    endif()

    add_executable( ${test_exe} testing.c ${test_src}  )
    add_dependencies( ${test_exe} ${KSVD_LIBNAME} )
    if( _build_netlib-scalapack )
        add_dependencies( ${test_exe} NETLIB-SCALAPACK )
    endif()

    target_link_libraries( ${test_exe} ${KSVD_LIBNAME} polar
                                       ${MPI_C_LIBRARIES}
                                       ${ELPA_LIBRARIES}
                                       ${SLTMG_LIBRARIES}
                                       ${SCALAPACK_LIBRARIES}
                                       ${LAPACKE_LIBRARIES}
                                       ${LAPACK_LIBRARIES}
                                       ${BLAS_LIBRARIES}
                                       ${EXTRA_LIBS} )
    install( TARGETS ${test_exe} DESTINATION lib/testing )
endforeach()
set( _p 4 )
set( _q 4 )
math( EXPR _pq "${_p} * ${_q}" )

# Add basic tests. To run them type: 'make test' in build directory
add_test( NAME "testing_pdgesvd:QDWH-MR" COMMAND mpirun -np ${_pq} ./testing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --ksvdmr --polarqdwh )
add_test( NAME "testing_pdgesvd:QDWH-DC" COMMAND mpirun -np ${_pq} ./testing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --ksvddc --polarqdwh )
add_test( NAME "testing_pdgesvd:QDWH-EL" COMMAND mpirun -np ${_pq} ./testing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --ksvdel )
add_test( NAME "testing_pdgesvd:ZOLOPD-MR" COMMAND mpirun -np ${_pq} ./testing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --ksvdmr --polarzolopd )
add_test( NAME "testing_pdgesvd:ZOLOPD-DC" COMMAND mpirun -np ${_pq} ./testing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --ksvddc --polarzolopd )
add_test( NAME "testing_pdgesvd:ZOLOPD-EL" COMMAND mpirun -np ${_pq} ./testing_pdgesvd --nprow ${_p} --npcol ${_q} --nb 64 --n_range 512:512:512 --niter 3 --check --ksvdel --polarzolopd )
