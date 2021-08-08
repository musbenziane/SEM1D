! Created by mus on 30/07/2021.

subroutine compute_forcemat(N, ne, h, nt, f0, dt, t, isrc, rho1Dgll, v1Dgll, xgll,Cij, u, F)
    implicit none
    integer, intent(in)                                :: N, ne, nt, t
    real (kind=4), intent(in)                          :: h, f0, dt, xgll
    real (kind=4), dimension(N*ne+1), intent(in)       :: rho1Dgll, v1Dgll
    real (kind=4), dimension(N*ne+1), intent(out)      :: F,u
    integer, dimension(N+1,ne)                         :: Cij
    real (kind=4), dimension(N+1,N+1)                  :: lprime
    real (kind=4), dimension(nt)                       :: src
    real (kind=4), dimension(N+1)                      :: xi, wi, fe,  b_force
    real (kind=4), dimension(N*ne+1)                   :: mu, mu1Dgll
    real (kind=4)                                      :: temp, Jc, Jci, sigma
    integer                                            :: isrc, i, j, k, c, el, ngll


    call lagrangeprime(N, lprime)
    call ricker(nt,f0,dt,src)
    call gll(N,xi,wi)

    !print*,src

    c =0
    ngll = N * ne + 1
    Jc = h/2
    Jci = 1/Jc
    mu1Dgll = 0
    mu1Dgll = rho1Dgll * v1Dgll**2
    c = 1
    temp = 0
    u = 0
    do el=1,ne
        fe = 0
        b_force = 0
        do i=1,N+1
            temp = 0
            do j=1,N+1
                temp = temp + u(Cij(j,el)) * lprime(i,j)
            end do
            sigma = mu(Cij(i,el)) * temp * Jci
            b_force(i) = sigma * Jc * Jci
        end do

        do j=1,N+1
            temp = 0
            do k=1,N+1
                temp = temp + b_force(k) * lprime(j,k) * wi(k)
            end do
            fe(j) = -temp
            if (Cij(j,el)==isrc) then
                fe(j) = fe(j) + src(t)
                !print*,src(t)
            end if
        end do
        do i=1,N
            F(Cij(i,el)) = F(Cij(i,el)) + fe(i)
        end do
    end do

end subroutine compute_forcemat