! Created by mus on 29/07/2021.

program create1Dmodel_files
    implicit none
    integer :: modtype, ns, nend, nel
    real (kind=8) :: vs, rho,rat
    real (kind=8), dimension(:), allocatable :: v1D, rho1D

    write(*,*) "Simple model -> 0 : vs = 2000 m/s, rho = 1800 kg/m³"
    write(*,*) "Homogeneous -> 1 | Heterogenous  [1 diffrent layer ] -> 2"
    read (*,*) modtype
    write(*,*) "Give the number of elements "
    read (*,*) nel

    allocate(v1D(nel))
    allocate(rho1D(nel))


    select case (modtype)

        case (0)
        v1D = 2000
        rho1D = 1800
        open(9,file="testing_vs",access='direct',recl=nel*8)
        write(9,rec=1) v1D
        close(9)

        open(8,file="testing_rho",access='direct',recl=nel*8)
        write(8,rec=1) rho1D
        close(8)

        deallocate(v1D)
        deallocate(rho1D)

        case (1)
        write(*,*) "Give velocity value & density value"
        read (*,*) vs, rho
        v1D = vs
        rho1D = rho

        open(10,file="testing_vs",access='direct',recl=nel*8)
        write(10,rec=1) v1D
        close(10)

        open(11,file="testing_rho",access='direct',recl=nel*8)
        write(11,rec=1) rho1D
        close(11)

        deallocate(v1D)
        deallocate(rho1D)

        case (2)
        write(*,*) "Give velocity value & density value of the background"
        read (*,*) vs, rho
        write(*,*) "Give ratio of the low/high velocity zone not in %"
        read (*,*) rat
        write(*,*) "give position in element number of the low/high velocity zone [nstart, nend]"
        read (*,*) ns, nend

        v1D = vs
        v1D(ns:nend) = vs * rat
        rho1D = rho

        open(12,file="testing_vs",access='direct',recl=nel*8)
        write(12,rec=1) v1D
        close(12)

        open(13,file="testing_rho",access='direct',recl=nel*8)
        write(13,rec=1) rho1D
        close(13)

        deallocate(v1D)
        deallocate(rho1D)

        case default
        write (*,*) "There has been an issue, start over"

    end select

    write(*,*) "Model files names are: VP1D.bin & RHOID.bin"

end program create1Dmodel_files