! Created by mus on 29/07/2021.
! Recursive calculation of Legendre polynomials

function legendrep(N,x) result(fac)
    !implicit none
    real (kind=8) :: fac
    real (kind=8), intent(in) :: x
    integer, intent(in) :: N
    real (kind=8), dimension(N+1) :: poly
    integer :: i

    select case (N)
        case (0)
        poly(1) = 1
        fac = poly(1)

        case (1)
        poly(2) = x
        fac = poly(2)

        case default
        poly(1) = 1
        poly(2) = x
    end select

    do i=3,N+1
        poly(i) = ((2*(i-1)-1)*x*poly(i-1) - (i-2)*poly(i-2))/(i-1)
    end do
    fac = poly(N+1)
end function legendrep