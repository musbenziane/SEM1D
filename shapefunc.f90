! Created by mus on 29/07/2021.

subroutine shapefunc(N,h,ne,xg)
    implicit none
    real (kind=8), intent(in) :: h
    integer, intent(in) :: N, ne
    real (kind=8), intent(out) :: xg(N*ne+1)
    integer :: i, j, c
    real (kind=8) :: xi(N+1), wi(N+1)

    call gll(N,xi,wi)

    c = 1
    do i=1,ne
        do j=1,N
            c = c + 1
            xg(c) = h*(i-1) +  h * ((xi(j+1) + 1) / 2)
        end do
    end do

end subroutine shapefunc