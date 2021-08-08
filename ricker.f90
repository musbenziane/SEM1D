! Ricker wavelet
subroutine ricker(nt,f0,dt,source)
    implicit none
    integer :: nt
    real (kind=4):: f0
    real (kind=4),dimension(nt)::source

    real (kind=4)     	:: dt,t0,t,a,pt
    integer 	:: i
    pt = 1/f0
    !n = int(2 * pt / dt)
    t0 = pt / dt
    a = 4 / pt

    do i=1,nt
        t = ((i + 1) - t0) * dt
        source(i) = -2. * a * t * exp(-(a * t)**2)
    end do

END subroutine ricker