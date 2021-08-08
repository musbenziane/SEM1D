cmake_minimum_required(VERSION 3.19)
project(SEM1D Fortran)



if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "-fopenmp")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()

enable_language(Fortran)
add_executable(SEM1D main.f90 gll.f90 lagrangep.f90 legendrep.f90 lagrangeprime.f90 connectivity_matrix.f90 readmodelfiles1D.f90 shapefunc.f90 mapmodel.f90 ricker.f90
    compute_forcemat.f90 read_sim_parameters.f90)
#add_executable(SEM1D test.f90)

