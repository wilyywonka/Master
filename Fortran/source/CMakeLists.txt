include_directories(${CMAKE_BINARY_DIR}/source/)

set(hdf5_i /usr/lib/x86_64-linux-gnu/hdf5/serial)
set(hdf5 /usr/lib/x86_64-linux-gnu/hdf5/serial)
include_directories(${hdf5_i})
link_directories (${hdf5})

set(hdf5_fortran_i /usr/lib/x86_64-linux-gnu/hdf5/serial)
set(hdf5_fortran /usr/lib/x86_64-linux-gnu/hdf5/serial)
include_directories(${hdf5_fortran_i})
link_directories (${hdf5_fortran})

set(hdf5hl_fortran_i /usr/lib/x86_64-linux-gnu/hdf5/serial)
set(hdf5hl_fortran /usr/lib/x86_64-linux-gnu/hdf5/serial)
include_directories(${hdf5hl_fortran_i})
link_directories (${hdf5hl_fortran})


set(h5fortran_i /home/wilhelm/Documents/h5fortran/Current/include)
set(h5fortran /home/wilhelm/Documents/h5fortran/Current/lib)
include_directories(${h5fortran_i})
link_directories (${h5fortran})

add_executable(my_program
	h5tests.f90
	testmodfile.f90
	subsubroutines.f90
	)

target_link_libraries(my_program h5fortran hdf5 hdf5_fortran hdf5hl_fortran)