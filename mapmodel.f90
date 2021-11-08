! Created by mus on 30/07/2021.

subroutine mapmodel(N,ne,rho1D,v1D,rho1Dgll,v1Dgll)
    implicit none

    interface connectivity_matrix
        function connectivity_matrix(N,ne)
            integer, intent(in)                                :: N, ne
            integer connectivity_matrix(N+1,ne)
        end function
    end interface

    integer, intent(in)               :: N, ne
    real (kind=8), intent(in)         :: rho1D(ne), v1D(ne)
    real (kind=8), intent(out)        :: rho1Dgll(N*ne+1), v1Dgll(N*ne+1)
    integer, dimension(N+1,ne)        :: Cij
    integer                           :: i, j

    Cij = connectivity_matrix(N,ne)

    do i=1,ne
        do j=1,N+1
            rho1Dgll(Cij(j,i)) = rho1D(i)
            v1Dgll(Cij(j,i))   = v1D(i)
        end do
    end do

end subroutine mapmodel