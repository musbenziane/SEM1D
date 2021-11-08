! Created by mus on 29/07/2021.

subroutine shapefunc(N,h,ne,Cij,xg)
    implicit none
    real (kind=8), intent(in)    :: h
    integer, intent(in)          :: N, ne, Cij(N+1,ne)
    real (kind=8), intent(out)   :: xg(ne*N+1)
    integer                      :: i, j
    real (kind=8)                :: xi(N+1), wi(N+1)

    call zwgljd(xi,wi,N+1,0.,0.)

    do i=1,ne
        do j=1,N+1
            xg(Cij(j,i)) = h*(i-1) +  h * ((xi(j) + 1) / 2)
        end do
    end do

end subroutine shapefunc