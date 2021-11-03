! Created by mus on 29/07/2021.
! Connectivty Matrix for global assembly
! See Bernhard, Schuberth. (2003) Thesis

function connectivity_matrix(N,ne)
    integer, intent(in)                             :: N, ne
    integer                                         :: i, j, c
    integer                                            connectivity_matrix(N+1,ne)

    c = 1
    do j=1,ne
        do i=1,N+1
            connectivity_matrix(i,j) = c
            c = c + 1
        end do
        c = c -1
    end do

end function