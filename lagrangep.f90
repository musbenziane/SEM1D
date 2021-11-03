! Created by mus on 28/07/2021.
! l :: dimension of x input
! N :: Polynomial order
! i :: i index in the lagrange polynomial ell
! fac :: resulting polynomial array
! xi :: gll points
! wi :: gll quadrature weights
subroutine lagrangep(N,i,x,fac,l)
    implicit none
    integer, intent(in) :: N, i,l
    real (kind=8), dimension(l),intent(in) :: x
    real (kind=8), dimension(l),intent(out) :: fac
    integer :: j
    real (kind=8) , dimension(N+1) :: xi,wi

    call zwgljd(xi,wi,N+1,0.,0.)
    fac = 1
    do j=0,N
        if (i /= j) then
            fac = fac * ((x - xi(j+1))/(xi(i+1)-xi(j+1)))
        end if
    end do
end subroutine lagrangep