! Created by mus on 29/07/2021.
! Derivation of Lagrange polynomials at the GLL points based on results from Funaro 1993 || Check Heiner Igel pg 196
! N: Polynomial order
! primed: derivative of lagrange polynomial matrix, each "k" row is the derivative of "lk(xi)" Lagrange polynomial. Each row is a the value of the derivative at a GLL point
! dummy arguments: dij -> derivation matrix | xi, wi -> GLL points and weights ...
! ... Lxi, Lxj -> Legendre polynomial value at xi and xj | lk -> k-th lagrange polynomial | sum -> summation return for the last derivation sum

subroutine lagrangeprime(N,primed)
    implicit none
    integer, intent(in)                            :: N
    real (kind=8), dimension(N+1,N+1), intent(out) :: primed
    real (kind=8), dimension(N+1,N+1)              :: dij
    real (kind=8), dimension(N+1)                  :: xi, wi
    real (kind=8)                                  :: Lxi, Lxj, sum, lk
    integer                                        :: i,j,k
    real (kind=8), external                        :: legendrep


    call zwgljd(xi,wi,N+1,0.,0.)


    do i=0,N
        do j=0,N
            if (i/=j) then
                Lxi = legendrep(N,xi(i+1))
                Lxj = legendrep(N,xi(j+1))
                dij(i+1,j+1) = (Lxi / Lxj) / (xi(i+1) - xi(j+1))
            elseif (i==j) then
                if (i==0 .and. j==0) then
                    dij(i+1,j+1) = -0.25 * N * (N + 1)
                elseif (i>=1 .and. i<=N-1) then
                    dij(i+1,j+1) = 0
                elseif (i==N) then
                    dij(i+1,j+1) = 0.25 * N *(N + 1)
                endif
            end if
        end do
    end do


    do k=0,N
        do i=0,N
            sum = 0
            do j=0,N
                call lagrangep(N,k,xi(j+1),lk,1)
                sum = sum + dij(i+1,j+1) * lk
            end do
            primed(k+1,i+1) = sum
        end do
    end do
end subroutine lagrangeprime