! Created by mus on 30/07/2021.

subroutine readmodelfiles1D(v1D, rho1D, ne)
    implicit none
    integer, intent(in) :: ne
    real (kind=4), dimension(ne), intent(out) :: v1D, rho1D

    open(10,file="V1D.bin",access='direct',recl=ne*4)
    read(10,rec=1) v1D
    close(10)


    open(11,file="RHO1D.bin",access='direct',recl=ne*4)
    read(11,rec=1) rho1D
    close(11)



end subroutine readmodelfiles1D