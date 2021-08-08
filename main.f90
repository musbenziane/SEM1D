! Created by mus on 28/07/2021.
! Spectral Elements 1D Solver with a regular mesh

program SEM1D
    USE OMP_LIB

    implicit none
    real (kind=4)                                              :: Jc, Jci, h, f0, dt, sum, CFL, mindist, lambdamin
    real (kind=4), dimension(:), allocatable                   :: xi, wi, v1D, rho1D, Me, M, xgll, rho1Dgll, v1Dgll, fe
    real (kind=4), dimension(:), allocatable                   :: mu1Dgll, u, uold, unew, F, src
    real (kind=4), dimension(:,:), allocatable                 :: lprime, Minv, Kg, Ke,Uout
    integer, dimension(:,:), allocatable                       :: Cij,kij
    integer                                                    :: N, ne, ngll, i, j, k, nt, isrc, t, el, isnap
    character(len=40)                                          :: filename, filecheck, outname

    outname = "snapshots.bin"

    write(*,*) "##########################################"
    write(*,*) "######## Reading parameters file #########"
    write(*,*) "##########################################"

    filename="parameters.in"

    print*,"Is the parameters input file (parameters.in) [Yes/no]"
    read(*,*) filecheck

    if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
            filecheck=="Y") then
        write(*,*) "Reading simulation parameters..."

    elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") then
        write(*,*) "Enter simulation parameters text file name with extension"
        write(*,*) "40 characters max"
        read(*,*) filename

    else
        write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
        write(*,*) "The program have been terminated, please star over"
        stop
    end if

    open (2, file=filename, status = 'old')
    read(2,*) N
    read(2,*) ne
    read(2,*) h
    read(2,*) f0
    read(2,*) dt
    read(2,*) nt
    read(2,*) isrc
    read(2,*) isnap
    close(2)

    print*,"Polynomial order         -> ",N
    print*,"Number of elements       -> ",ne
    print*,"Element size             -> ",h
    print*,"Wavelet's peak frequency -> ",f0
    print*,"Time step                -> ",dt
    print*,"Number of time steps     -> ",nt
    print*,"Source location          -> ",isrc
    print*,"Snapshot interval        -> ",isnap

    ngll = N * ne + 1                  ! Total GLL points
    Jc = h / 2                         ! Jacobian for structured 1D mesh
    Jci = 1 / Jc                       ! Jacobian inverse


    allocate(xi(N+1))                   ! GLL points
    allocate(wi(N+1))                   ! GLL Quadrature weights
    allocate(v1D(ne))                   ! 1D velocity model in elements
    allocate(rho1D(ne))                 ! Density velocity model in elements
    allocate(rho1Dgll(ngll))            ! 1D density model mapped
    allocate(v1Dgll(ngll))              ! 1D velocity mapped
    allocate(mu1Dgll(ngll))             ! Shear modulus mapped
    allocate(xgll(ngll))                ! Array for global mapping
    allocate(Cij(N+1,ne))               ! Connectivity matrix
    allocate(M(ngll))                   ! Global mass matrix in vector form
    allocate(Minv(ngll,ngll))           ! Inverse of the mass matrix
    allocate(Me(N+1))                   ! Elemental mass matrix
    allocate(Kg(ngll,ngll))             ! Global stifness matrix
    allocate(Ke(N+1,N+1))               ! Elemental stifness matrix
    allocate(lprime(N+1,N+1))           ! Dervatives of Lagrange polynomials
    allocate(u(ngll))                   ! Displacement vector at time t
    allocate(uold(ngll))                ! displacement vector at time t + dt
    allocate(unew(ngll))                ! displacement vecotr at time t - dt
    allocate(src(nt))                   ! Source time function
    allocate(F(ngll))                   ! External force
    allocate(Uout(nt,ngll))             ! Snapshots



    call gll(N,xi,wi)                              ! Getting GLL points and weights
    call readmodelfiles1D(v1D, rho1D, ne)          ! Reading model files
    call connectivity_matrix(N,ne,Cij)             ! Getting connectivity matrix
    call shapefunc(N,h,ne, xgll)                   ! Global domain mapping
    call mapmodel(N,ne, rho1D,v1D,rho1Dgll,v1Dgll) ! Mapping models
    call lagrangeprime(N,lprime)                   ! Lagrange polynomials derivatives
    call ricker(nt,f0,dt,src)                      ! Source time function



    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"
    write(*,*)"##########################################"

    mindist = xgll(2) - xgll(1)
    CFL = (dt/mindist) * maxval(v1D(:))
    if (CFL > .8) then
        print"(a14,f6.3)","CFL value is ",CFL
        print*,"Decrease time step, the program has been terminated"
        stop

    else
        print"(a14,f6.3)","CFL value is ",CFL
        print*,"Simulation is stable"
    end if



    write(*,*)"##########################################"
    write(*,*)"########## Space Sampling check ##########"
    write(*,*)"##########################################"

    lambdamin = minval(v1D)/(f0*2.5)

    print"(a32,f3.1)", " Elements per minimum wavelength ->", lambdamin/h
    if ((lambdamin/h)<1) then
        print*,"Element size is too large"
        print*,"Numerical dispersion might be present"
        print*,"Do you wish to continue anyways? Yes/no"
        read*, filecheck

        if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
                filecheck=="Y") then
            print*,"Proceeding..."

        elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") then
            write(*,*) "The program had been terminated"
            write(*,*) "Reduce element size or frequency as such"
            write(*,*) "you have at least 1 element per min wavelength"
            write(*,*) "This holds for N=5 to N=10 as per"
            write(*,*) "Komatitsch and Tromp 1999"
            stop
        else
            write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
            write(*,*) "The program have been terminated, please star over"
            stop
        end if
    else
        print*, "Spatial sampling as OK!"
    end if


    !##########################################
    !####### Construct the mass matrix ########
    !##########################################

    do i=1,ne
        do j=1,N+1
            Me(j) = rho1Dgll(Cij(j,i)) * wi(j) * Jc   ! Elemental mass matrix construction
            M(Cij(j,i)) =  M(Cij(j,i)) + Me(j)        ! Filling the global mass matrix
        end do
    end do

    Minv(:,:) = 0
    ! Invert mass matrix
    do i=1,ngll
        Minv(i,i) = 1 / M(i)
    end do


    !###############################################
    !####### Construct the Stiffness matrix ########
    !###############################################

    mu1Dgll(:) = rho1Dgll(:) * v1Dgll(:)**2.          ! Shear modulus
    Kg(:,:) = 0
    do el=1,ne
        do i=1,N+1                                    ! Elemental stifness matrix
            do j=1,N+1
                sum = 0.
                do k=1,N+1
                    sum = sum + mu1Dgll(Cij(k,el)) * wi(k) * lprime(j,k) * lprime(i,k) * Jc * Jci**2
                end do
                Ke(i,j) = sum
            end do
        end do

        do i=1,N+1                                      ! Global assembly
            do j=1,N+1
                Kg(Cij(i,el),Cij(j,el)) = Kg(Cij(i,el),Cij(j,el)) + Ke(i,j)
            end do
        end do
    end do


    write(*,*) "!##########################################"
    write(*,*) "############ BEGIN TIME LOOP ##############"
    write(*,*) "###########################################"

    uold(:) = 0.
    u(:)    = 0.
    unew(:) = 0.
    F(:)    = 0.

    do t=1,nt
        ! Injecting source
        F(isrc) = src(t)

        unew = (dt**2.) * matmul(Minv,(F - matmul(Kg,u)))  + 2. * u - uold
        uold = u
        u = unew


        if ( mod(t,isnap) == 0) then
            Uout(t,:) = u
            if (mod(t,50) == 0) then
                print*,"At time sample ->",t, "/",nt
            end if
        end if
    end do


    open(15,file=outname,access="direct",recl=nt*ngll*4)
    write(15,rec=1) Uout
    close(15)

    deallocate(xi,wi,v1D,v1Dgll,rho1D,rho1Dgll,mu1Dgll,xgll,Cij,M,Me,Minv)
    deallocate(u,uold,unew,lprime,Kg,Ke)
    deallocate(Uout)
end program

