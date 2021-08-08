! Created by mus on 06/08/2021.

program test
    implicit none
    integer :: i, N, ne, nt, isnap, isrc
    real ::  h, f0, dt


    open (2, file = 'parameters.in', status = 'old')

    read(2,*) N
    read(2,*) ne
    read(2,*) h
    read(2,*) f0
    read(2,*) dt
    read(2,*) nt
    read(2,*) isrc
    read(2,*) isnap



    close(2)


    write(*,*) N
    write(*,*) ne
    write(*,*) h
    write(*,*) f0
    write(*,*) dt
    write(*,*) nt
    write(*,*) isrc
    write(*,*) isnap


end program test