! Created by mus on 29/07/2021.
! Connectivty Matrix for global assembly
! See Bernhard, Schuberth. (2003) Thesis

subroutine connectivity_matrix(N,ne,Cij)
    implicit none
    integer, intent(in) :: N, ne
    integer, dimension(N+1,ne), intent(out) :: Cij
    integer :: i, j, c

    c = 1
    do j=1,ne
        do i=1,N+1
            Cij(i,j) = c
            c = c + 1
        end do
        c = c -1
    end do



end subroutine connectivity_matrix